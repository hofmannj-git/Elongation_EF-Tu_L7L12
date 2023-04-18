module mod_statistics

  use, intrinsic :: iso_fortran_env
  use mod_group, only : group_t, MAX_GROUP
  implicit none
  private

  public :: group_stats

  type, public :: running_stat_t
     integer(int64), private :: n=0
     double precision, private :: m1=0,m2=0,m3=0,m4=0,a=0
     logical, private :: ema_flag = .false.
   contains
     procedure, public :: init => rs_init
     procedure, public :: push => rs_push
     procedure, public :: count => rs_count
     procedure, public :: clear => rs_clear
     procedure, public :: mean => rs_mean
     procedure, public :: variance => rs_variance
     procedure, public :: stdev => rs_stdev
     procedure, public :: sterr => rs_sterr
     procedure, public :: skewness => rs_skewness
     procedure, public :: kurtosis => rs_kurtosis
     procedure, public :: set_ema => rs_set_ema
     procedure, private :: push_ema => rs_push_ema
  end type running_stat_t

  type, public :: group_stat_t
     type(running_stat_t), dimension(MAX_GROUP) :: rs_vec
   contains
     procedure, public :: init => gs_init
     procedure, public :: push => gs_push
     procedure, public :: count => gs_count
     procedure, public :: clear => gs_clear
     procedure, public :: mean => gs_mean
     procedure, public :: variance => gs_variance
     procedure, public :: stdev => gs_stdev
     procedure, public :: sterr => gs_sterr
     procedure, public :: skewness => gs_skewness
     procedure, public :: kurtosis => gs_kurtosis
  end type group_stat_t

contains

  ! ----------------------------------------------------------------------
  ! set to ema mode with parameter a
  ! ----------------------------------------------------------------------

  subroutine rs_set_ema(this, a)

    class(running_stat_t), intent(inout) :: this
    double precision, intent(in) :: a

    this%ema_flag = .true.
    this%a = a

  end subroutine rs_set_ema

  ! ----------------------------------------------------------------------
  ! Update fields of exponential moving average running_stat.
  ! Non-supported fields stay zero.
  ! ----------------------------------------------------------------------

  subroutine rs_push_ema(this, x)

    class(running_stat_t), intent(inout) :: this
    real(real64), intent(in) :: x
    real(real64) :: delta,delta_alpha

    this%n = this%n + 1
    delta = x - this%m1
    delta_alpha = this%a*delta
    this%m1 = this%a*x + (1.0_real64 - this%a)*this%m1
    this%m2 = (1.0_real64 - this%a)*(this%m2 + delta*delta_alpha)

  end subroutine rs_push_ema

  ! ----------------------------------------------------------------------
  ! get stats for per-particle x array
  ! use only indices specified by 'group'
  ! return zeros if the group is empty
  ! below, vec(1) = mean
  !        vec(2) = stdev
  !        vec(3) = sterr
  ! ----------------------------------------------------------------------

  subroutine group_stats(group, igroup, x, vec)

    type(group_t), intent(in) :: group
    integer, intent(in) :: igroup
    real, dimension(:), intent(in) :: x
    real, dimension(3), intent(out) :: vec
    type(running_stat_t) :: rs
    integer :: nparticles,i,groupbit

    nparticles = group%particle%nparticles
    groupbit = group%bitmask(igroup)

    ! NOTE: skipping producing error if igroup does not exist

    call rs%init()
    do i = 1, nparticles
       if (iand(group%mask(i),groupbit) /= 0) call rs%push(real(x(i),real64))
    end do

    vec(1) = real(rs%mean())
    vec(2) = real(rs%stdev())
    vec(3) = real(rs%sterr())

  end subroutine group_stats

  ! ----------------------------------------------------------------------
  ! init the running stat type
  ! ----------------------------------------------------------------------

  subroutine rs_init(this)
    class(running_stat_t), intent(out) :: this
    call this%clear()
  end subroutine rs_init

  ! ----------------------------------------------------------------------
  ! insert new value to running stat type
  ! ----------------------------------------------------------------------

  subroutine rs_push(this, x)

     class(running_stat_t) :: this
     real(real64), intent(in) :: x
     real(real64) :: delta,delta_n,delta_n2,term1
     integer(int64) :: n1

     if (this%ema_flag) then
        call this%push_ema(x)
        return
     end if

     n1 = this%n
     this%n = this%n + 1
     delta = x - this%m1
     delta_n = delta / this%n
     delta_n2 = delta_n * delta_n
     term1 = delta * delta_n * n1
     this%m1 = this%m1 + delta_n
     this%m4 = this%m4 &
          + term1 * delta_n2 * ( this%n*this%n - 3*this%n + 3 ) &
          + 6 * delta_n2 * this%m2 &
          - 4 * delta_n * this%m3
     this%m3 = this%m3 &
          + term1 * delta_n * (this%n-2) &
          - 3 * delta_n * this%m2
     this%m2 = this%m2 + term1

   end subroutine rs_push

  ! ----------------------------------------------------------------------
  ! count number of entries added to running stat type
  ! ----------------------------------------------------------------------

   function rs_count(this)
     class(running_stat_t), intent(in) :: this
     integer(int64) :: rs_count
     rs_count = this%n
   end function rs_count

  ! ----------------------------------------------------------------------
  ! clear the running stat
  ! ----------------------------------------------------------------------

  subroutine rs_clear(this)

    class(running_stat_t), intent(inout) :: this

    this%n = 0
    this%m1 = 0.0
    this%m2 = 0.0
    this%m3 = 0.0
    this%m4 = 0.0

  end subroutine rs_clear

  ! ----------------------------------------------------------------------
  ! return running stats: mean, variance, standard deviation,
  ! standard error, skewness, kurtosis
  ! ----------------------------------------------------------------------

   function rs_mean(this) result (x)
     class(running_stat_t), intent(in) :: this
     real(real64) :: x
     x = this%m1
   end function rs_mean

   ! ----------------------------------------------------------------------

   function rs_variance(this) result (x)
     class(running_stat_t), intent(in) :: this
     real(real64) :: x
     if (this%n > 1) then
        x = this%m2 / (this%n-1)
     else
        x = 0.0_real64
     end if
   end function rs_variance

   ! ----------------------------------------------------------------------

   function rs_stdev(this) result (x)
     class(running_stat_t), intent(in) :: this
     real(real64) :: x
     x = sqrt(rs_variance(this))
   end function rs_stdev

   ! ----------------------------------------------------------------------

   function rs_sterr(this) result (x)
     class(running_stat_t), intent(in) :: this
     real(real64) :: x
     if (this%n > 0) then
        x = rs_stdev(this) / sqrt(real(this%n,real64))
     else
        x = 0.0_real64
     end if
   end function rs_sterr

   ! ----------------------------------------------------------------------

   function rs_skewness(this) result (x)
     class(running_stat_t), intent(in) :: this
     real(real64) :: x
     if (this%m2 > 0.0_real64) then
        x = sqrt(real(this%n,real64)) * this%m3 / this%m2**1.5_real64
     else
        x = 0.0_real64
     end if
   end function rs_skewness

   ! ----------------------------------------------------------------------

   function rs_kurtosis(this) result (x)
     class(running_stat_t), intent(in) :: this
     real(real64) :: x
     if (this%m2 > 0.0_real64) then
        x = this%n * this%m4 / (this%m2*this%m2) &
             - 3.0_real64
     else
        x = 0.0_real64
     end if
   end function rs_kurtosis

   ! ----------------------------------------------------------------------
   ! group stat functions
   ! ----------------------------------------------------------------------

   subroutine gs_init(this)
     class(group_stat_t), intent(inout) :: this
     call this%clear()
   end subroutine gs_init

   ! ----------------------------------------------------------------------

   subroutine gs_push(this, group, vals)

     class(group_stat_t), intent(inout) :: this
     type(group_t), intent(in) :: group
     real(real64), dimension(:), intent(in) :: vals
     integer :: i,j,groupbit,nparticles

     nparticles = group%particle%nparticles

    ! NOTE: skipping producing error if igroup does not exist

     do i = 1, group%ngroup
        groupbit = group%bitmask(i)
        do j = 1, nparticles
           if (iand(group%mask(j),groupbit) /= 0) then
              call this%rs_vec(i)%push(vals(j))
           end if
        end do
     end do

   end subroutine gs_push

   ! ----------------------------------------------------------------------

   function gs_count(this, igroup) result(n)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     integer(int64) :: n
     n = this%rs_vec(igroup)%count()
   end function gs_count

   ! ----------------------------------------------------------------------

   subroutine gs_clear(this)
     class(group_stat_t), intent(inout) :: this
     integer :: i
     do i = 1, size(this%rs_vec)
        call this%rs_vec(i)%clear()
     end do
   end subroutine gs_clear

   ! ----------------------------------------------------------------------

   function gs_mean(this, igroup) result(x)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     real(real64) :: x
     x = this%rs_vec(igroup)%mean()
   end function gs_mean

   ! ----------------------------------------------------------------------

   function gs_variance(this, igroup) result(x)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     real(real64) :: x
     x = this%rs_vec(igroup)%variance()
   end function gs_variance

   ! ----------------------------------------------------------------------

   function gs_stdev(this, igroup) result(x)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     real(real64) :: x
     x = this%rs_vec(igroup)%stdev()
   end function gs_stdev

   ! ----------------------------------------------------------------------

   function gs_sterr(this, igroup) result (x)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     real(real64) :: x
     x = this%rs_vec(igroup)%sterr()
   end function gs_sterr

   ! ----------------------------------------------------------------------

   function gs_skewness(this, igroup) result (x)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     real(real64) :: x
     x = this%rs_vec(igroup)%skewness()
   end function gs_skewness

   ! ----------------------------------------------------------------------

   function gs_kurtosis(this, igroup) result (x)
     class(group_stat_t), intent(in) :: this
     integer, intent(in) :: igroup
     real(real64) :: x
     x = this%rs_vec(igroup)%kurtosis()
   end function gs_kurtosis

   ! ----------------------------------------------------------------------

end module mod_statistics
