module mod_image

  use iso_fortran_env
  use RCImageIO
  implicit none
  private

  public :: symmetrize_halfspace_rot_2d, symmetrize_halfspace_rot_3d
  public :: render_array
  public :: corner_to_center_2d, corner_to_center_3d
  public :: center_to_corner_2d, center_to_corner_3d
  public :: displacement_to_affine, scale_to_affine
  public :: affine_compose, affine_transform, affine_invert,affine_interpolate
  public :: lower_voigt_to_affine, upper_voigt_to_affine
!!$  public :: affine_transform_rot, affine_transform_shear
!!$  public :: affine_transform_scale, affine_transform_translate
!!$  public :: affine_transform, affine_compose

contains

  ! ----------------------------------------------------------------------
  ! Render an image from a 2D array with a given color scheme and
  ! threshold values.  There is an option to use logarithmic scaling.
  ! Current palettes are 'jt' and 'gs' for jet and grayscale.
  ! IMPORTANT: I assume that the 0,0 coordinate (or 1,1 index in Fortran)
  !            should be in the lower left-hand corner of the image.
  !            However, in images, pixel 0,0 is in the upper left-hand
  !            corner.  Therefore, I flip the vertical array index to
  !            get the pixel.  E.g., array index (81,2) goes to pixel
  !            index (80,48) for array/image dimensions of 100x50.
  ! ----------------------------------------------------------------------

  subroutine render_array(a, xmin, xmax, logflag, palette, fname)

    use RCImageIO

    real(real64), dimension(:,:), intent(inout) :: a
    real(real64), intent(inout) :: xmin,xmax
    logical, intent(in) :: logflag
    character(len=2), intent(in) :: palette
    character(len=*), intent(in) :: fname
    type(rgbimage) :: img
    type(rgb) :: color
    integer :: i,j,nx,ny
    real(real64), parameter :: tol = 0.000001

    nx = size(a,1)
    ny = size(a,2)

    if (xmin > xmax) then
       write(*, '(a)') 'xmin is greater than xmax in render_array'
       stop
    end if
    if (logflag) then
       do j = 1, ny
          do i = 1, nx
             if (a(i,j) < 0.0) then
                write(*, '(a)') 'cannot have positive value with logflag'
                stop
             end if
             if (a(i,j) < tol) then
                a(i,j) = tol
             end if
             a(i,j) = log(a(i,j))
          end do
       end do
       xmin = real(log(real(xmin)),real64)
       xmax = real(log(real(xmax)),real64)
    end if
    if ((palette /= 'jt') .and. (palette /= 'gs')) then
       write(*, '(a)') 'unknown palette'
       stop
    end if

    call init_img(img)
    call alloc_img(img,nx,ny)
    do j = 1, ny
       do i = 1, nx
          if (palette == 'jt') then
             color = val2rgb_jt(a(i,j), xmin, xmax)
          else
             color = val2rgb_gs(a(i,j), xmin, xmax)
          end if
          call put_pixel(img, i-1, ny-j, color)
       end do
    end do

    open(8, file=trim(adjustl(fname)), action='write')
    call output_ppm(8, img)
    call free_img(img)
    close(8)

  end subroutine render_array

  ! ----------------------------------------------------------------------
  ! Take the 'half-space' output (from a DFT from FFTW3, for example),
  ! and use rotational symmetry to fill in the other 'half-space'.
  ! Output array must be already allocated and correct dimensions.
  ! 2D version.
  ! Assume that the half space is in x, or the 1st index.
  ! If DC point is at 1,1 in input array, it is in 1,1 in output array.
  ! ----------------------------------------------------------------------

  subroutine symmetrize_halfspace_rot_2d(ai, ao)

    real(real64), dimension(:,:), intent(in) :: ai
    real(real64), dimension(:,:), allocatable, intent(inout) :: ao
    integer :: i,j,nxi,nyi

    nxi = size(ai,1)
    nyi = size(ai,2)
    allocate(ao(2*nxi-1,nyi))

    ao(1:nxi, 1:nyi/2+1) = ai(:,1:nyi/2+1)
    forall (i=2:nxi, j=2:nyi)
       ao(2*nxi-i+1,nyi+2-j) = ai(i,j)
    end forall
    ao(nxi+1:2*nxi-1,1) = ai(nxi:2:-1,1)

  end subroutine symmetrize_halfspace_rot_2d

  ! ----------------------------------------------------------------------
  ! Take the 'half-space' output (from a DFT from FFTW3, for example),
  ! and use rotational symmetry to fill in the other 'half-space'.
  ! Output array must be already allocated and correct dimensions.
  ! 3D version.
  ! Assume that the half space is in x, or the 1st index.
  ! If DC point is at 1,1,1 in input array, it is in 1,1,1 in output array.
  ! ----------------------------------------------------------------------

  subroutine symmetrize_halfspace_rot_3d(ai, ao)

    real(real64), dimension(:,:,:), intent(in) :: ai
    real(real64), dimension(:,:,:), allocatable, intent(inout) :: ao
    integer :: i,j,k,nxi,nyi,nzi

    nxi = size(ai,1)
    nyi = size(ai,2)
    nzi = size(ai,3)
    allocate(ao(2*nxi-1,nyi,nzi))

    ao(1:nxi, 1:nyi, 1:nzi) = ai(:,:,:)
    forall (i=2:nxi, j=2:nyi, k=2:nzi)
       ao(2*nxi-i+1,nyi+2-j,nzi+2-k) = ai(i,j,k)
    end forall
    ao(nxi+1:2*nxi-1,1,2:nzi) = ai(nxi:2:-1,1,nzi:2:-1)
    ao(nxi+1:2*nxi-1,2:nyi,1) = ai(nxi:2:-1,nyi:2:-1,1)
    ao(nxi+1:2*nxi-1,1,1) = ai(nxi:2:-1,1,1)

  end subroutine symmetrize_halfspace_rot_3d

  ! ----------------------------------------------------------------------
  ! Functions that turn real values into pixel colors (16-bit)
  ! ----------------------------------------------------------------------

  function val2rgb_gs(x, xmin, xmax) result (color)

    use RCImageIO

    real(real64), intent(in) :: x,xmin,xmax
    type(rgb) :: color
    integer :: c

    c = int((x-xmin)/(xmax-xmin)*256)
    if (c < 0) c = 0
    if (c > 255) c = 255
    color%red = c
    color%green = c
    color%blue = c

  end function val2rgb_gs

  function val2rgb_jt(x, xmin, xmax) result (color)

    real(real64), intent(in) :: x,xmin,xmax
    type(rgb) :: color
    integer :: c

    c = int((x-xmin)/(xmax-xmin)*64)
    select case (c)
    case (:0)
       color%red = 0
       color%green = 0
       color%blue = 143
    case (1)
       color%red = 0
       color%green = 0
       color%blue = 159
    case (2)
       color%red = 0
       color%green = 0
       color%blue = 175
    case (3)
       color%red = 0
       color%green = 0
       color%blue = 191
    case (4)
       color%red = 0
       color%green = 0
       color%blue = 207
    case (5)
       color%red = 0
       color%green = 0
       color%blue = 223
    case (6)
       color%red = 0
       color%green = 0
       color%blue = 239
    case (7)
       color%red = 0
       color%green = 0
       color%blue = 255
    case (8)
       color%red = 0
       color%green = 15
       color%blue = 255
    case (9)
       color%red = 0
       color%green = 31
       color%blue = 255
    case (10)
       color%red = 0
       color%green = 47
       color%blue = 255
    case (11)
       color%red = 0
       color%green = 63
       color%blue = 255
    case (12)
       color%red = 0
       color%green = 79
       color%blue = 255
    case (13)
       color%red = 0
       color%green = 95
       color%blue = 255
    case (14)
       color%red = 0
       color%green = 111
       color%blue = 255
    case (15)
       color%red = 0
       color%green = 127
       color%blue = 255
    case (16)
       color%red = 0
       color%green = 143
       color%blue = 255
    case (17)
       color%red = 0
       color%green = 159
       color%blue = 255
    case (18)
       color%red = 0
       color%green = 175
       color%blue = 255
    case (19)
       color%red = 0
       color%green = 191
       color%blue = 255
    case (20)
       color%red = 0
       color%green = 207
       color%blue = 255
    case (21)
       color%red = 0
       color%green = 223
       color%blue = 255
    case (22)
       color%red = 0
       color%green = 239
       color%blue = 255
    case (23)
       color%red = 0
       color%green = 255
       color%blue = 255
    case (24)
       color%red = 15
       color%green = 255
       color%blue = 239
    case (25)
       color%red = 31
       color%green = 255
       color%blue = 223
    case (26)
       color%red = 47
       color%green = 255
       color%blue = 207
    case (27)
       color%red = 63
       color%green = 255
       color%blue = 191
    case (28)
       color%red = 79
       color%green = 255
       color%blue = 175
    case (29)
       color%red = 95
       color%green = 255
       color%blue = 159
    case (30)
       color%red = 111
       color%green = 255
       color%blue = 143
    case (31)
       color%red = 127
       color%green = 255
       color%blue = 127
    case (32)
       color%red = 143
       color%green = 255
       color%blue = 111
    case (33)
       color%red = 159
       color%green = 255
       color%blue = 95
    case (34)
       color%red = 175
       color%green = 255
       color%blue = 79
    case (35)
       color%red = 191
       color%green = 255
       color%blue = 63
    case (36)
       color%red = 207
       color%green = 255
       color%blue = 47
    case (37)
       color%red = 223
       color%green = 255
       color%blue = 31
    case (38)
       color%red = 239
       color%green = 255
       color%blue = 15
    case (39)
       color%red = 255
       color%green = 255
       color%blue = 0
    case (40)
       color%red = 255
       color%green = 239
       color%blue = 0
    case (41)
       color%red = 255
       color%green = 223
       color%blue = 0
    case (42)
       color%red = 255
       color%green = 207
       color%blue = 0
    case (43)
       color%red = 255
       color%green = 191
       color%blue = 0
    case (44)
       color%red = 255
       color%green = 175
       color%blue = 0
    case (45)
       color%red = 255
       color%green = 159
       color%blue = 0
    case (46)
       color%red = 255
       color%green = 143
       color%blue = 0
    case (47)
       color%red = 255
       color%green = 127
       color%blue = 0
    case (48)
       color%red = 255
       color%green = 111
       color%blue = 0
    case (49)
       color%red = 255
       color%green = 95
       color%blue = 0
    case (50)
       color%red = 255
       color%green = 79
       color%blue = 0
    case (51)
       color%red = 255
       color%green = 63
       color%blue = 0
    case (52)
       color%red = 255
       color%green = 47
       color%blue = 0
    case (53)
       color%red = 255
       color%green = 31
       color%blue = 0
    case (54)
       color%red = 255
       color%green = 15
       color%blue = 0
    case (55)
       color%red = 255
       color%green = 0
       color%blue = 0
    case (56)
       color%red = 239
       color%green = 0
       color%blue = 0
    case (57)
       color%red = 223
       color%green = 0
       color%blue = 0
    case (58)
       color%red = 207
       color%green = 0
       color%blue = 0
    case (59)
       color%red = 191
       color%green = 0
       color%blue = 0
    case (60)
       color%red = 175
       color%green = 0
       color%blue = 0
    case (61)
       color%red = 159
       color%green = 0
       color%blue = 0
    case (62)
       color%red = 143
       color%green = 0
       color%blue = 0
    case (63:)
       color%red = 127
       color%green = 0
       color%blue = 0
    end select

  end function val2rgb_jt

  ! ----------------------------------------------------------------------
  ! Shift elements circularly so that element 1,1 is at element
  ! nx/2+1,ny/2+1.  This will be exactly in the center for
  ! dimensions with an odd number of elements.
  ! ----------------------------------------------------------------------

  subroutine corner_to_center_2d(a)

    real(real64), dimension(:,:), intent(inout) :: a
    integer :: is,js

    is = size(a,1)/2! + 1
    js = size(a,2)/2! + 1
    a = cshift(a, -is, 1)
    a = cshift(a, -js, 2)

  end subroutine corner_to_center_2d

  ! ----------------------------------------------------------------------
  ! Shift elements circularly so that element 1,1,1 is at element
  ! nx/2+1,ny/2+1,nz/2+1.  This will be exactly in the center for
  ! dimensions with an odd number of elements.
  ! ----------------------------------------------------------------------

  subroutine corner_to_center_3d(a)

    real(real64), dimension(:,:,:), intent(inout) :: a
    integer :: is,js,ks

    is = size(a,1)/2! + 1
    js = size(a,2)/2! + 1
    ks = size(a,3)/2! + 1
    a = cshift(a, -is, 1)
    a = cshift(a, -js, 2)
    a = cshift(a, -ks, 3)

  end subroutine corner_to_center_3d

  ! ----------------------------------------------------------------------
  ! Shift elements circularly so that element nx/2+1,ny/2+1 is at element
  ! 1,1.  This is the opposite of corner_to_center_2d.
  ! If the array has even indices in all directions, corner_to_center
  ! does the same thing as center_to_corner, but not in the most general
  ! case.
  ! ----------------------------------------------------------------------

  subroutine center_to_corner_2d(a)

    real(real64), dimension(:,:), intent(inout) :: a
    integer :: is,js
!
    is = size(a,1)/2! + 1
    js = size(a,2)/2! + 1
    a = cshift(a, is, 1)
    a = cshift(a, js, 2)

  end subroutine center_to_corner_2d

  ! ----------------------------------------------------------------------
  ! Shift elements circularly so that element nx/2+1,ny/2+1,nz/2+1 is at
  ! element 1,1,1.  This is the opposite of corner_to_center_3d.
  ! If the array has even indices in all directions, corner_to_center
  ! does the same thing as center_to_corner, but not in the most general
  ! case.
  ! ----------------------------------------------------------------------

  subroutine center_to_corner_3d(a)

    real(real64), dimension(:,:,:), intent(inout) :: a
    integer :: is,js,ks

    is = size(a,1)/2! + 1
    js = size(a,2)/2! + 1
    ks = size(a,3)/2! + 1
    a = cshift(a, is, 1)
    a = cshift(a, js, 2)
    a = cshift(a, ks, 3)

  end subroutine center_to_corner_3d

  ! ----------------------------------------------------------------------
  ! General affine transformation.  Inputs are a 3D array and a vector
  ! representation of a 4x4 matrix with thre zero entries used to define
  ! the affine transformation.
  ! Can limit the range of interpolation with optional bnds.
  ! E.g., if bnds = (/ 1, 10, 2, 2, 2, 1000 /) then interpolation is
  ! done only to get the index box of size 10x1x999 located between
  ! the index values described.
  ! This can be used to just perform interpolation for a plane of data
  ! rather than waste time performing interpolation for the whole box.
  ! ----------------------------------------------------------------------
  ! I change it to a point-wise transform (instead of matrix-wise)
  ! b/c the memory issue

  subroutine affine_transform(a, m, b, i, j, k)

    real(real64), dimension(:,:,:), intent(inout) :: a
    real(real64), dimension(12), intent(in) :: m
    real(real64), intent(inout) :: b
    !integer, dimension(6), intent(in), optional :: bnds
    integer, intent(in) :: i,j,k
    real(real64), dimension(12) :: minv
    real(real64) :: xr,yr,zr
    integer :: imin,imax,jmin,jmax,kmin,kmax

    !imin = 1
    !jmin = 1
    !kmin = 1
    !imax = size(b,1)
    !jmax = size(b,2)
    !kmax = size(b,3)
    !if (present(bnds)) then
    !   do n = 1, 6
    !      if (bnds(n) < 1) then
    !         write(*,'(a)') 'bound in bnds must be >= 1'
    !         stop
    !      end if
    !   end do
    !   do n = 1, 3
    !      if (bnds(2*n) < bnds(2*n-1)) then
    !         write(*,'(a)') 'bound in bnds must be low index, then high index'
    !         stop
    !      end if
    !      if (bnds(2*n) > size(b,n)) then
    !         write(*,'(a)') 'upper index bound must be < output array size'
    !         stop
    !      end if
    !   end do
    !   imin = bnds(1)
    !   imax = bnds(2)
    !   jmin = bnds(3)
    !   jmax = bnds(4)
    !   kmin = bnds(5)
    !   kmax = bnds(6)
    !end if

    minv = affine_invert(m)

    !do concurrent (i=imin:imax, j=jmin:jmax, k=kmin:kmax)
       xr = minv(1)*(i-1) + minv( 2)*(j-1) + minv( 3)*(k-1) + minv( 4) + 1.0
       yr = minv(5)*(i-1) + minv( 6)*(j-1) + minv( 7)*(k-1) + minv( 8) + 1.0
       zr = minv(9)*(i-1) + minv(10)*(j-1) + minv(11)*(k-1) + minv(12) + 1.0
       b = affine_interpolate(a, xr, yr, zr)
    !end do

  end subroutine affine_transform

  ! ----------------------------------------------------------------------
  ! Invert an affine transformation matrix, using vector representation.
  ! This is used for inverse mapping in affine_transform.
  ! Will exit with error if non-invertible: ae-bd=0 for
  !
  !     | a b c d |
  !     | e f g h |
  ! m = | i j k l | = (/ a, b, c, d, e, f, g, h, i j, k, l /)
  !     | 0 0 0 1 |
  !
  ! ----------------------------------------------------------------------

  function affine_invert(m) result(n)

    real(real64), dimension(12), intent(in) :: m
    real(real64), dimension(12) :: n
    real(real64) :: det

    det = m(1)*m(6)*m(11) - m(1)*m(7)*m(10) - m(2)*m(5)*m(11) &
        + m(2)*m(7)*m( 9) + m(3)*m(5)*m(10) - m(3)*m(6)*m(9)
    if (det == 0.0_real64) then
       write(*, '(a)') 'determinant zero in affine_transform_invert'
       stop
    end if

    n( 1) = m(6)*m(11) - m(7)*m(10)
    n( 2) = m(3)*m(10) - m(2)*m(11)
    n( 3) = m(2)*m( 7) - m(3)*m( 6)
    n( 5) = m(7)*m( 9) - m(5)*m(11)
    n( 6) = m(1)*m(11) - m(3)*m( 9)
    n( 7) = m(3)*m( 5) - m(1)*m( 7)
    n( 9) = m(5)*m(10) - m(6)*m( 9)
    n(10) = m(2)*m( 9) - m(1)*m(10)
    n(11) = m(1)*m( 6) - m(2)*m( 5)

    n( 4) =  m(4)*m(7)*m(10) - m(3)*m(8)*m(10) - m(4)*m(6)*m(11) &
           + m(2)*m(8)*m(11) + m(3)*m(6)*m(12) - m(2)*m(7)*m(12)
    n( 8) = -m(4)*m(7)*m( 9) + m(3)*m(8)*m( 9) + m(4)*m(5)*m(11) &
           - m(1)*m(8)*m(11) - m(3)*m(5)*m(12) + m(1)*m(7)*m(12)
    n(12) =  m(4)*m(6)*m( 9) - m(2)*m(8)*m( 9) - m(4)*m(5)*m(10) &
           + m(1)*m(8)*m(10) + m(2)*m(5)*m(12) - m(1)*m(6)*m(12)

    n(:) = n(:) / det

  end function affine_invert

  ! ----------------------------------------------------------------------
  ! Convert a lower triangular matrix in Voigt notation to an affine.
  !
  !              | a 0 0 0 |
  ! | a 0 0 |    | f b 0 0 |
  ! | f b 0 | -> | e d c 0 |
  ! | e d c |    | 0 0 0 1 |
  !
  ! ...in vectors: (/ a, b, ..., f /) -> (/ a, 0, 0, 0, f, b, ..., c, 0 /)
  ! 
  ! ----------------------------------------------------------------------

  function lower_voigt_to_affine(h) result (m)

    real(real64), dimension(6), intent(in) :: h
    real(real64), dimension(12) :: m

    m( 1) = h(1)
    m( 2) = 0.0_real64
    m( 3) = 0.0_real64
    m( 4) = 0.0_real64
    m( 5) = h(6)
    m( 6) = h(2)
    m( 7) = 0.0_real64
    m( 8) = 0.0_real64
    m( 9) = h(5)
    m(10) = h(4)
    m(11) = h(3)
    m(12) = 0.0_real64

  end function lower_voigt_to_affine

  ! ----------------------------------------------------------------------
  ! Convert an upper triangular matrix in Voigt notation to an affine.
  !
  !              | a f e 0 |
  ! | a f e |    | 0 b d 0 |
  ! | 0 b d | -> | 0 0 c 0 |
  ! | 0 0 c |    | 0 0 0 1 |
  !
  ! ...in vectors: (/ a, b, ..., f /) -> (/ a, 0, 0, 0, b, f, ..., e, 0 /)
  ! 
  ! ----------------------------------------------------------------------

  function upper_voigt_to_affine(h) result (m)

    real(real64), dimension(6), intent(in) :: h
    real(real64), dimension(12) :: m

    m( 1) = h(1)
    m( 2) = h(6)
    m( 3) = h(5)
    m( 4) = 0.0_real64
    m( 5) = 0.0_real64
    m( 6) = h(2)
    m( 7) = h(4)
    m( 8) = 0.0_real64
    m( 9) = 0.0_real64
    m(10) = 0.0_real64
    m(11) = h(3)
    m(12) = 0.0_real64

  end function upper_voigt_to_affine

  ! ----------------------------------------------------------------------
  ! Convert displacement in 3D to an affine.
  !
  !               | 1 0 0 x |
  !               | 0 1 0 y |
  ! (dx,dy,dz) -> | 0 0 1 z |
  !               | 0 0 0 1 |
  !
  ! ----------------------------------------------------------------------

  function displacement_to_affine(x, y, z) result (m)

    real(real64), intent(in) :: x,y,z
    real(real64), dimension(12) :: m

    m( 1) = 1.0_real64
    m( 2) = 0.0_real64
    m( 3) = 0.0_real64
    m( 4) = x
    m( 5) = 0.0_real64
    m( 6) = 1.0_real64
    m( 7) = 0.0_real64
    m( 8) = y
    m( 9) = 0.0_real64
    m(10) = 0.0_real64
    m(11) = 1.0_real64
    m(12) = z

  end function displacement_to_affine

  ! ----------------------------------------------------------------------
  ! Scale an image by factors of (sx,sy,xz).
  ! ----------------------------------------------------------------------

  function scale_to_affine(sx, sy, sz) result (m)

    real(real64), intent(in) :: sx,sy,sz
    real(real64), dimension(12) :: m

    m( 1) = sx
    m( 2) = 0.0_real64
    m( 3) = 0.0_real64
    m( 4) = 0.0_real64
    m( 5) = 0.0_real64
    m( 6) = sy
    m( 7) = 0.0_real64
    m( 8) = 0.0_real64
    m( 9) = 0.0_real64
    m(10) = 0.0_real64
    m(11) = sz
    m(12) = 0.0_real64

  end function scale_to_affine

  ! ----------------------------------------------------------------------
  ! Compose two affines
  ! ----------------------------------------------------------------------

  function affine_compose(m1, m2) result (m3)

    real(real64), dimension(12), intent(in) :: m1,m2
    real(real64), dimension(12) :: m3

    m3( 1) = m1(1)*m2(1) + m1( 2)*m2(5) + m1(3)*m2( 9)
    m3( 2) = m1(1)*m2(2) + m1( 2)*m2(6) + m1(3)*m2(10)
    m3( 3) = m1(1)*m2(3) + m1( 2)*m2(7) + m1(3)*m2(11)
    m3( 5) = m1(5)*m2(1) + m1( 6)*m2(5) + m1(7)*m2( 9)
    m3( 6) = m1(5)*m2(2) + m1( 6)*m2(6) + m1(7)*m2(10)
    m3( 7) = m1(5)*m2(3) + m1( 6)*m2(7) + m1(7)*m2(11)
    m3( 9) = m1(9)*m2(1) + m1(10)*m2(5) + m1(11)*m2( 9)
    m3(10) = m1(9)*m2(2) + m1(10)*m2(6) + m1(11)*m2(10)
    m3(11) = m1(9)*m2(3) + m1(10)*m2(7) + m1(11)*m2(11)

    m3( 4) = m1( 4) + m1(1)*m2(4) + m1( 2)*m2(8) + m1( 3)*m2(12)
    m3( 8) = m1( 8) + m1(5)*m2(4) + m1( 6)*m2(8) + m1( 7)*m2(12)
    m3(12) = m1(12) + m1(9)*m2(4) + m1(10)*m2(8) + m1(11)*m2(12)

  end function affine_compose

  ! ----------------------------------------------------------------------
  ! Use bilinear interpolation to find an intermediate value using 4
  ! neighboring pixels.  If a pixel is off of the grid, consider its
  ! value to be zero.
  ! ----------------------------------------------------------------------

  pure function affine_interpolate(a,xr,yr,zr) result (v)

    real(real64), dimension(:,:,:), intent(in) :: a
    real(real64), intent(in) :: xr,yr,zr
    real(real64) :: v
    real(real64) :: dx,dy,dz,dxm,dym,dzm
    real(real64), dimension(8) :: vals
    integer :: i,j,k,m,n,o

    m = size(a,1)
    n = size(a,2)
    o = size(a,3)
    i = int(xr)
    j = int(yr)
    k = int(zr)
    dx = xr - real(i,real64)
    dy = yr - real(j,real64)
    dz = zr - real(k,real64)
    dxm = 1_real64 - dx
    dym = 1_real64 - dy
    dzm = 1_real64 - dz

    ! vals(1) relates to a(  i,   j,   k)
    ! vals(2) relates to a(i+1,   j,   k)
    ! vals(3) relates to a(  i, j+1,   k)
    ! vals(4) relates to a(  i,   j, k+1)
    ! vals(5) relates to a(i+1,   j, k+1)
    ! vals(6) relates to a(  i, j+1, k+1)
    ! vals(7) relates to a(i+1, j+1,   k)
    ! vals(8) relates to a(i+1, j+1, k+1)

    vals(:) = 0.0
    if ((  i > 0) .and. (  i <= m) .and. &
        (  j > 0) .and. (  j <= n) .and. &
        (  k > 0) .and. (  k <= o)) then
       vals(1) = a(i,j,k)
    end if
    if ((i+1 > 0) .and. (i+1 <= m) .and. &
        (  j > 0) .and. (  j <= n) .and. &
        (  k > 0) .and. (  k <= o)) then
       vals(2) = a(i+1,j,k)
    end if
    if ((  i > 0) .and. (  i <= m) .and. &
        (j+1 > 0) .and. (j+1 <= n) .and. &
        (  k > 0) .and. (  k <= o)) then
       vals(3) = a(i,j+1,k)
    end if
    if ((  i > 0) .and. (  i <= m) .and. &
        (  j > 0) .and. (  j <= n) .and. &
        (k+1 > 0) .and. (k+1 <= o)) then
       vals(4) = a(i,j,k+1)
    end if
    if ((i+1 > 0) .and. (i+1 <= m) .and. &
        (  j > 0) .and. (  j <= n) .and. &
        (k+1 > 0) .and. (k+1 <= o)) then
       vals(5) = a(i+1,j,k+1)
    end if
    if ((  i > 0) .and. (  i <= m) .and. &
        (j+1 > 0) .and. (j+1 <= n) .and. &
        (k+1 > 0) .and. (k+1 <= o)) then
       vals(6) = a(i,j+1,k+1)
    end if
    if ((i+1 > 0) .and. (i+1 <= m) .and. &
        (j+1 > 0) .and. (j+1 <= n) .and. &
        (  k > 0) .and. (  k <= o)) then
       vals(7) = a(i+1,j+1,k)
    end if
    if ((i+1 > 0) .and. (i+1 <= m) .and. &
        (j+1 > 0) .and. (j+1 <= n) .and. &
        (k+1 > 0) .and. (k+1 <= o)) then
       vals(8) = a(i+1,j+1,k+1)
    end if

    v = vals(1) * dxm * dym * dzm &
      + vals(2) *  dx * dym * dzm &
      + vals(3) * dxm *  dy * dzm &
      + vals(4) * dxm * dym *  dz &
      + vals(5) *  dx * dym *  dz &
      + vals(6) * dxm *  dy *  dz &
      + vals(7) *  dx *  dy * dzm &
      + vals(8) *  dx *  dy *  dz

  end function affine_interpolate


!!$  ! ----------------------------------------------------------------------
!!$  ! Rotate an image by c radians.
!!$  ! ----------------------------------------------------------------------
!!$
!!$  function affine_transform_rot(c) result (m)
!!$
!!$    real(real64), intent(in) :: c
!!$    real(real64), dimension(6) :: m
!!$
!!$    m = (/ cos(c), -sin(c), 0.0_real64, sin(c), cos(c), 0.0_real64 /)
!!$
!!$  end function affine_transform_rot
!!$
!!$  ! ----------------------------------------------------------------------
!!$  ! Shear an image by an engineering strain of c.
!!$  ! ----------------------------------------------------------------------
!!$
!!$  function affine_transform_shear(c) result (m)
!!$
!!$    real(real64), intent(in) :: c
!!$    real(real64), dimension(6) :: m
!!$
!!$    m = (/ 1.0_real64, c, 0.0_real64, &
!!$         0.0_real64, 1.0_real64, 0.0_real64 /)
!!$
!!$  end function affine_transform_shear
!!$
!!$  ! ----------------------------------------------------------------------
!!$  ! Translate an array by (dx,dy), in pixel units.
!!$  ! ----------------------------------------------------------------------
!!$
!!$  function affine_transform_translate(dx, dy) result (m)
!!$
!!$    real(real64), intent(in) :: dx,dy
!!$    real(real64), dimension(6) :: m
!!$
!!$    m = (/ 1.0_real64, 0.0_real64, dx, 0.0_real64, 1.0_real64, dy /)
!!$
!!$  end function affine_transform_translate
!!$
!!$
!!$

end module mod_image
