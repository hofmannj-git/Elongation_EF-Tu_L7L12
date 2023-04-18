program plot_ssf

  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_c_binding
  use mod_particle, only : particle_t
  use mod_domain, only : domain_t
  use mod_group, only : group_t
  use mod_cells, only : cells_t
  use mod_statistics, only : running_stat_t
  use mod_stat_func, only : eval_ssf_3d,init_ffts,close_ffts
  use mod_data, only : source_data_t
  use mod_constants, only : pi
  use mod_image, only : render_array,symmetrize_halfspace_rot_3d,&
                        corner_to_center_3d,affine_transform,&
                        affine_invert,affine_interpolate
  implicit none

  integer :: i,j,k,n,m,nparticles,nargs,snap,nsnaps,nq,qbin
  type(particle_t) :: particle
  type(domain_t) :: domain
  type(group_t) :: group
  type(cells_t) :: cells
  type(source_data_t) :: source_data
  type(running_stat_t), dimension(:,:), allocatable :: rs_array_slice_yz
  type(running_stat_t), dimension(:,:), allocatable :: rs_array_slice_xz
  type(running_stat_t), dimension(:,:), allocatable :: rs_array_slice_xy
  type(running_stat_t), dimension(:), allocatable :: rs_array_rad, rs_ssf_rad
  type(running_stat_t) :: rs_ssf_rad_min
  character(len=300) :: fname,buf,directory
  real(real64) :: sqmin,sqmax,qx,v_shell,r,ssf_interp,ssf_rad_min,&
                  sqminlog,sqmaxlog
  real :: target_size,qmax,delq,qmin
  integer(int64) :: tstart,tdur,tinc,tstep,tout
  integer, parameter :: nyquist_multiplier = 1
  ! not using 3d array here to save memory
  ! real(real64), dimension(:,:,:), allocatable :: fin 
  real, dimension(:,:,:), allocatable :: ssf
  real(real64), dimension(:,:), allocatable :: sq_yz,sq_xz,sq_xy
  real(real64), dimension(:,:), allocatable :: sq_yz_center,sq_xz_center,&
                                               sq_xy_center
  real, dimension(:), allocatable :: sq_rad
  real(real64), dimension(:,:,:), allocatable :: qxx,qxy,qxz,sym
  real(real64), dimension(:), allocatable :: ssf_rad,qx_bin
  real(real64), dimension(6) :: h_inv
  real(real64), dimension(12) :: mcomp,minv
  real(real64) :: qval,dx,dy,dz,scale,xr,yr,zr,sq_rad_min

  ! read inputs

  nargs = command_argument_count()
  if (nargs .ne. 8) then
     write(*, '(a)') 'plot_ssf [directory] [tstart] [tdur] [tinc] &
                               &[qmax] [nq] [sqmin] [sqmax]'
     stop
  end if
  call get_command_argument(1, buf); directory = trim(adjustl(buf))
  call get_command_argument(2, buf); read(buf, *) tstart
  call get_command_argument(3, buf); read(buf, *) tdur
  call get_command_argument(4, buf); read(buf, *) tinc
  call get_command_argument(5, buf); read(buf, *) qmax
  call get_command_argument(6, buf); read(buf, *) nq
  call get_command_argument(7, buf); read(buf, *) sqmin
  call get_command_argument(8, buf); read(buf, *) sqmax

  if (modulo(tdur,tinc) /= 0) then
     write(*,'(a)') 'tdur is not a multiple of tinc'
     stop
  end if
  if (qmax <= 0.0) then
     write(*,'(a)') 'qmax must be > 0.0'
     stop
  end if
  if (nq < 1) then
     write(*,'(a)') 'nq must be > 0'
     stop
  end if
  if (modulo(tdur,tinc) /= 0) then
     write(*,'(a)') 'tdur is not a multiple of tinc'
     stop
  end if
  if (sqmax <= 0.0) then
     write(*,'(a)') 'sqmax must be > 0.0'
     stop
  end if
  if (sqmin < 0.0) then
     write(*,'(a)') 'sqmin must be >= 0.0'
     stop
  end if
  if (sqmin >= sqmax) then
     write(*,'(a)') 'sqmax must be > sqmin'
     stop
  end if

  ! Init the key classes.

  call source_data%init(trim(adjustl(directory)), particle, domain)
  if (.not. source_data%step_exists(tstart)) then
     write(*,'(a)') 'bad tstart value'
     stop
  end if
  if (.not. source_data%regular) then
     write(*,'(a)') 'hdf5 timesteps appear irregular'
     write(*,'(a)') 'not supporting nonzero tdur for irregular'
     stop
  end if
  nsnaps = int(tdur/tinc) + 1
  nparticles = particle%nparticles
  call group%init(domain)

  ! Get wavenumber spacing, and allocate output arrays.
  ! Number of samples in S(q) vs. q is nq.
  ! Number of pixels in output image is 2*nq-1.

  delq = qmax / (nq-1)
  allocate(qx_bin(nq))
  forall (i=1:nq) qx_bin(i) = (i-1) * delq
  allocate(rs_array_rad(nq-1))
  allocate(sq_rad(nq-1))
  allocate(rs_array_slice_yz(2*nq-1,2*nq-1))
  allocate(rs_array_slice_xz(2*nq-1,2*nq-1))
  allocate(rs_array_slice_xy(2*nq-1,2*nq-1))
  allocate(sq_yz(2*nq-1,2*nq-1))
  allocate(sq_xz(2*nq-1,2*nq-1))
  allocate(sq_xy(2*nq-1,2*nq-1))
  allocate(sq_yz_center(nq,nq))
  allocate(sq_xz_center(nq,nq))
  allocate(sq_xy_center(nq,nq))
  do i=1,2*nq-1
     do j=1,2*nq-1 
        call rs_array_slice_yz(i,j)%init
        call rs_array_slice_xz(i,j)%init
        call rs_array_slice_xy(i,j)%init
     end do
  end do
  allocate(rs_ssf_rad(nq-1)) ! This is used to get statistics for each step
  do i=1,nq-1 
     call rs_ssf_rad(i)%init
  end do
  call rs_ssf_rad_min%init ! This is also used to get statistics for each step

  ! Determine real-space sampling required to resolve maximum wavenumber.
  ! In 1D, need sample increment to be less than (1/2)*(2*pi/qmax).
  ! In 3D, seems okay to find the maximum triclinic sub-cell that fits
  ! inside of a sphere of diameter of (1/2)*(2*pi/qmax).
  ! The .true. option specifies to use the maximizing procedure noted above.

  ! In this code we will fix nyquist_multiplier to 1; To reduce the aliasing
  ! effect the user need to specify a proper qmax.
  ! For a reference, qmax=8.0, lx=125.419 and no-flip box, the fft nodes are 
  ! 452x452x452, which takes 452x452x452x8/1024/1024=704 MB (since here fft 
  ! is out-of-place, there might be a factor of 2. This is the biggest system
  ! I can do on my virtual machine b/c I only allocate 2 GB to RAM.)

  target_size = pi / qmax / nyquist_multiplier / sqrt(2.0) 
  call cells%init(group, target_size, .true.)

  ! Loop over snaps.

  do n = 1, nsnaps

     tstep = tstart + (n-1)*tinc
     snap = source_data%step_index(tstep)

     call source_data%read_snap(particle,domain,snap)
     call cells%reset_grid()
     qmin = 2.0*sqrt(2.0)*qmax/cells%ncellx
     write(*,'(a,3i6)') 'fft nodes:', cells%ncellx, cells%ncelly, cells%ncellz
     call cells%assign_particles()
     call eval_ssf_3d(cells, ssf)
     ! quick and dirty fix: set DC value to 0 so that the interpolation will not
     ! make results ugly
     ssf(1,1,1) = 0.0

     ! doing symmetrizing and circular shift outside render_ssf

     call symmetrize_halfspace_rot_3d(real(ssf,real64), sym)
     call corner_to_center_3d(sym)

     h_inv(:) = real(domain%h_inv(:),real64)

     ! Peform the affine transformation from the input to the image.

     call affine_transform_forward(sym, mcomp, h_inv, real(qmax,real64), nq)

     ! invert the transform matrix to calculate the real coords each pixel maps to
     ! easier to find interpolation points this way
 
     minv = affine_invert(mcomp)

     ! render image (not rendering the whole array to save memory)

     ! call render_ssf(ssf, fin, h_inv, tstep, &
     !      sqmin, sqmax, real(qmax,real64), nq)


     do i=1,2*nq-1
        do j=1,2*nq-1
           do k=1,2*nq-1
              m =     (i-nq+1) * (i-nq+1)
              m = m + (j-nq+1) * (j-nq+1)
              m = m + (k-nq+1) * (k-nq+1)
              qx = sqrt(real(m,real64))*real(qmax,real64)/real(nq-1,real64)
              ! qmax is the maximum wavenumber we want to see on the image.
              ! qmin is the wavenumber spacing from the raw data; any q below qmin is 
              ! not physical.
              if (qx < qmax .and. qx > qmin) then
                 call decide_qx_bin(qx,qx_bin,qbin)

                 ! Doing point-wise affine transformation b/c the memory issue

                 xr = minv(1)*(i-1) + minv( 2)*(j-1) + minv( 3)*(k-1) + minv( 4) + 1.0
                 yr = minv(5)*(i-1) + minv( 6)*(j-1) + minv( 7)*(k-1) + minv( 8) + 1.0
                 zr = minv(9)*(i-1) + minv(10)*(j-1) + minv(11)*(k-1) + minv(12) + 1.0
                 ssf_interp = affine_interpolate(sym, xr, yr, zr)
                 call rs_ssf_rad(qbin)%push(real(ssf_interp,real64))
              end if
           end do
        end do
     end do
     ! There are only 6 points at q=qmin, so for the qmin bin S(q) will be the average
     ! of these 6 points. 
     ssf_rad_min = 1.0/6.0*(sym(size(sym,1)/2+1+1,size(sym,2)/2+1,size(sym,3)/2+1) + &
                            sym(size(sym,1)/2+1-1,size(sym,2)/2+1,size(sym,3)/2+1) + &
                            sym(size(sym,1)/2+1,size(sym,2)/2+1+1,size(sym,3)/2+1) + &
                            sym(size(sym,1)/2+1,size(sym,2)/2+1-1,size(sym,3)/2+1) + &
                            sym(size(sym,1)/2+1,size(sym,2)/2+1,size(sym,3)/2+1+1) + &
                            sym(size(sym,1)/2+1,size(sym,2)/2+1,size(sym,3)/2+1-1))
     call rs_ssf_rad_min%push(real(ssf_rad_min,real64))

     ! Update pixel information on 3 slices (yz,xz,xy) at each time step

     do i=1,2*nq-1
        do j=1,2*nq-1
           call affine_transform(sym, mcomp, ssf_interp, nq, i, j)
           sq_yz(i,j) = ssf_interp

           call affine_transform(sym, mcomp, ssf_interp, i, nq, j)
           sq_xz(i,j) = ssf_interp

           call affine_transform(sym, mcomp, ssf_interp, i, j, nq)
           sq_xy(i,j) = ssf_interp
        end do
     end do

     ! Output images at each step if needed

     ! call image_ssf(sq_xy, 'yz', tstep, sqmin, sqmax, real(qmax,real64), nq)
     ! call image_ssf(sq_xy, 'xz', tstep, sqmin, sqmax, real(qmax,real64), nq)
     ! call image_ssf(sq_xy, 'xy', tstep, sqmin, sqmax, real(qmax,real64), nq)

     ! update running average for 3 slices

     do i=1,2*nq-1
        do j=1,2*nq-1
           call rs_array_slice_yz(i,j)%push(real(sq_yz(i,j), real64))
           call rs_array_slice_xz(i,j)%push(real(sq_xz(i,j), real64))
           call rs_array_slice_xy(i,j)%push(real(sq_xy(i,j), real64))
        end do
     end do

     ! update running average for S(q) (this is only needed for large qmax; 
     ! no need to do this for small qmax)

     do i = 1,nq-1
        call rs_array_rad(i)%push(real(rs_ssf_rad(i)%mean(), real64))
     end do

     ! Have to deallocate ssf and sym, since eval_ssf and symmetrize allocates them 

     deallocate(ssf,sym)

     ! Have to clear those two stats variables for the current step
     do i=1,nq-1 
        call rs_ssf_rad(i)%clear
     end do
     call rs_ssf_rad_min%clear 

  end do

  ! output running average

  ! images of 3 slices (yz,xz,xy)
  do i = 1, 2*nq-1
     do j = 1, 2*nq-1
        sq_yz(i,j) = rs_array_slice_yz(i,j)%mean()
        sq_xz(i,j) = rs_array_slice_xz(i,j)%mean()
        sq_xy(i,j) = rs_array_slice_xy(i,j)%mean()
     end do
  end do
  sq_yz_center=sq_yz(nq/2+1:nq/2+nq,nq/2+1:nq/2+nq)
  sq_xz_center=sq_xz(nq/2+1:nq/2+nq,nq/2+1:nq/2+nq)
  sq_xy_center=sq_xy(nq/2+1:nq/2+nq,nq/2+1:nq/2+nq)
  sqminlog = sqmin
  sqmaxlog = sqmax
  call image_ssf(sq_yz_center, 'yz-runningave', tstep, sqminlog, sqmaxlog, real(qmax,real64), nq)
  sqminlog = sqmin
  sqmaxlog = sqmax
  call image_ssf(sq_xz_center, 'xz-runningave', tstep, sqminlog, sqmaxlog, real(qmax,real64), nq)
  sqminlog = sqmin
  sqmaxlog = sqmax
  call image_ssf(sq_xy_center, 'xy-runningave', tstep, sqminlog, sqmaxlog, real(qmax,real64), nq)
  ! S(q) vs q
  sq_rad_min = real(rs_ssf_rad_min%mean())
  call decide_qx_bin(real(qmin,real64),qx_bin,qbin)
  qval = qmin
  ! Use upper bound for each bin as the bin location 
  ! to avoid overlap of the 1st bin (qmin)
  ! Previously used bin center as the bin location
  print*, qval, sq_rad_min
  do i = 1,nq-1
     if (qx_bin(i+1) > qmin) then
        sq_rad(i) = real(rs_array_rad(i)%mean())
        qval = qx_bin(i+1)
        print*, qval, sq_rad(i)
     end if 
  end do

contains

  ! combine all forward affine transform together 
  !   to save computation costs

  subroutine affine_transform_forward(sym,mcomp, h_inv, qmax, nq)
     
    use mod_image, only : displacement_to_affine, affine_compose, &
         lower_voigt_to_affine, scale_to_affine

    real(real64), dimension(:,:,:), intent(inout) :: sym 
    real(real64), dimension(6), intent(in) :: h_inv
    real(real64), dimension(12) :: m
    real(real64), dimension(12), intent(out) :: mcomp
    real(real64), intent(in) :: qmax
    integer, intent(in) :: nq

    dx = real(size(sym,1)/2,real64)
    dy = real(size(sym,2)/2,real64)
    dz = real(size(sym,3)/2,real64)
    mcomp = displacement_to_affine(-dx,-dy,-dz)

    ! h_inv is a voigt_upper matrix
    ! saying that it is a voigt_lower matrix is taking its transpose:
    ! H^{-T}

    m(:) = lower_voigt_to_affine(h_inv(:))
    mcomp(:) = affine_compose(m(:), mcomp(:))

    ! Above would be used to go from FT of Lambda coords
    ! to FT of real coords, but we also have to scale to
    ! get things in terms of pixels.

    scale = pi * (2*nq-1) / qmax
    m(:) = scale_to_affine(scale, scale, scale)
    mcomp(:) = affine_compose(m(:), mcomp(:))

    ! Shift the corner back to the center.
    ! The shift is in terms of pixels in the output image.

    dx = real(nq-1,real64)
    m(:) = displacement_to_affine(dx, dx, dx)
    mcomp(:) = affine_compose(m(:), mcomp(:))

  end subroutine affine_transform_forward

  ! ----------------------------------------------------------------------
  ! Render an image of the ssf
  ! 'flava' can be 'raw' for raw FFTW output, 'sym' for symmetrized
  ! version, or 'fin' for finalized version.  Only the finalized version
  ! has the correct final wavevectors.  'raw' and 'sym' are just used
  ! for debugging purposes.
  ! ----------------------------------------------------------------------

  subroutine render_ssf(sym, fin, i, j, k, h_inv, tstep, sqmin, sqmax, qmax, nq)

    use mod_image, only : render_array, symmetrize_halfspace_rot_3d, &
         corner_to_center_3d, displacement_to_affine, affine_compose, &
         lower_voigt_to_affine, affine_transform, scale_to_affine

    real(real64), dimension(:,:,:), intent(inout) :: sym 
    real(real64), dimension(6), intent(in) :: h_inv
    integer(int64), intent(in) :: tstep
    real(real64), intent(in) :: sqmin,sqmax
    real(real64), intent(in) :: qmax
    integer, intent(in) :: nq,i,j,k
    !real(real64), dimension(:,:,:), allocatable :: sym
    real(real64), intent(out) :: fin
    character(len=300) :: fname,buf
    real(real64), dimension(12) :: m,mcomp
    real(real64) :: dx,dy,dz,scale
    integer, dimension(6) :: bnds

    ! if ((plane /= 'yz') .and. (plane .ne. 'xz') .and. (plane /= 'xy')) then
    !    write(*,'(a)') 'need yz, xz, or xy for plane in render_ssf'
    !    stop
    ! end if

    !allocate(fin(2*nq-1,2*nq-1,2*nq-1))

    ! Symmetrize the ssf, allocating an array in the process.
    ! All this does is fills in the redundant (x) halfspace
    ! of S(q).

    !call symmetrize_halfspace_rot_3d(real(ssf,real64), sym)

    ! Temporarily bring corner, where DC is, to the center of the
    ! symmetric array.  This way we will get the S(q) rings in
    ! the center of the image.  Below performs cyclic shifting.

    !call corner_to_center_3d(sym)

    ! Compose linear transformations, to be performed in one shot
    ! to avoid losing pixels off of the edge of the image, to do
    ! the following:
    !
    !   1. Shift the DC point back to the edge.
    !   2. Map each lambda-wavevector to original wavevectors by
    !          q_x = H^{-T} q_\lambda
    !   3. Shift the DC point back to the center.
    !
    ! I'm doing this by composing these 'affines', finding their
    ! inverse, and finding which pixels to interpolate from in
    ! the raw data to generate the 2D image.

    dx = real(size(sym,1)/2,real64)
    dy = real(size(sym,2)/2,real64)
    dz = real(size(sym,3)/2,real64)
    mcomp = displacement_to_affine(-dx,-dy,-dz)

    ! h_inv is a voigt_upper matrix
    ! saying that it is a voigt_lower matrix is taking its transpose:
    ! H^{-T}

    m(:) = lower_voigt_to_affine(h_inv(:))
    mcomp(:) = affine_compose(m(:), mcomp(:))

    ! Above would be used to go from FT of Lambda coords
    ! to FT of real coords, but we also have to scale to
    ! get things in terms of pixels.

    scale = pi * (2*nq-1) / qmax
    m(:) = scale_to_affine(scale, scale, scale)
    mcomp(:) = affine_compose(m(:), mcomp(:))

    ! Shift the corner back to the center.
    ! The shift is in terms of pixels in the output image.

    dx = real(nq-1,real64)
    m(:) = displacement_to_affine(dx, dx, dx)
    mcomp(:) = affine_compose(m(:), mcomp(:))

    ! Peform the affine transformation from the input to the image.

    !bnds = (/1, 2*nq-1, 1, 2*nq-1, 1, 2*nq-1/) ! transform the whole array
    call affine_transform(sym, mcomp, fin, i, j, k)

  end subroutine render_ssf

  ! Suppose to transform pixel indices to qx coordinates
  ! Need further consideration

  subroutine pixel2x_ssf(sq_xyz, qxx, qxy, qxz, nq, qmax)

    use mod_image, only : displacement_to_affine, affine_compose, &
         scale_to_affine, affine_invert

    real(real64), dimension(:,:,:), intent(in) :: sq_xyz 
    real(real64), dimension(:,:,:), allocatable, intent(out) :: qxx,qxy,qxz
    real(real64), intent(in) :: qmax
    integer, intent(in) :: nq
    integer :: i,j,k
    real(real64), dimension(12) :: m,mcomp,minv
    real(real64) :: scale,dx

    ! Compose linear transformations, to be performed in one shot
    ! to avoid losing pixels off of the edge of the image, to do
    ! the following:
    !
    !   1. Shift the DC point back to the edge.
    !   2. Map each lambda-wavevector to original wavevectors by
    !          q_x = H^{-T} q_\lambda
    !   3. Shift the DC point back to the center.
    !
    ! I'm doing this by composing these 'affines', finding their
    ! inverse, and finding which pixels to interpolate from in
    ! the raw data to generate the 2D image.

    ! Above would be used to go from FT of Lambda coords
    ! to FT of real coords, but we also have to scale to
    ! get things in terms of pixels.

    scale = pi * (2*nq-1) / qmax
    mcomp(:) = scale_to_affine(scale, scale, scale)

    ! Shift the corner back to the center.
    ! The shift is in terms of pixels in the output image.

    m(:) = displacement_to_affine(dx, dx, dx)
    mcomp(:) = affine_compose(m(:), mcomp(:))

    minv = affine_invert(m)

    allocate(qxx(size(sq_xyz,1),size(sq_xyz,2),size(sq_xyz,3)))
    allocate(qxy(size(sq_xyz,1),size(sq_xyz,2),size(sq_xyz,3)))
    allocate(qxz(size(sq_xyz,1),size(sq_xyz,2),size(sq_xyz,3)))

    do concurrent (i=1:size(sq_xyz,1), j=1:size(sq_xyz,2), k=1:size(sq_xyz,3))
       qxx(i,j,k) = minv(1)*(i-1) + minv( 2)*(j-1) + minv( 3)*(k-1) + minv( 4)!)*2.0*pi/125.419
       qxy(i,j,k) = minv(5)*(i-1) + minv( 6)*(j-1) + minv( 7)*(k-1) + minv( 8)!)*2.0*pi/125.419
       qxz(i,j,k) = minv(9)*(i-1) + minv(10)*(j-1) + minv(11)*(k-1) + minv(12)!)*2.0*pi/125.419
    end do

  end subroutine pixel2x_ssf

  ! Output image

  subroutine image_ssf(fin, plane, tstep, sqmin, sqmax, qmax, nq)

    use mod_image, only : render_array

    real(real64), dimension(:,:), intent(inout) :: fin 
    character(len=*), intent(in) :: plane
    integer(int64), intent(in) :: tstep
    real(real64), intent(inout) :: sqmin,sqmax
    real(real64), intent(in) :: qmax
    integer, intent(in) :: nq
    character(len=300) :: fname,buf,buf1
    integer, dimension(6) :: bnds
   
    write(buf,'(i20)') tstep
    write(buf1,'(f10.1)') qmax
    write(fname,'(5a)') 'ssf-', plane, '-', trim(adjustl(buf)), '.ppm'
    call render_array(fin(:,:),sqmin,sqmax,.true.,'jt',fname)

  end subroutine image_ssf

  ! Determine which qx_bin it's in for all elements in 3D pixel array

  subroutine decide_qx_bin(qx, qx_bin, qbin)

    real(real64), intent(in) :: qx
    real(real64), dimension(:), intent(in) :: qx_bin
    integer, intent(out) :: qbin
    integer :: i,nbin

    nbin = size(qx_bin)-1
    do i = 1,nbin
       ! To avoid round-off error, here I used strictly '<' and '>'
       ! May miss a few points which are exactly on the bin boundaries
       if (qx < qx_bin(i+1) .and. qx > qx_bin(i) ) then
          qbin = i
          EXIT
       end if
    end do

  end subroutine decide_qx_bin

end program plot_ssf
