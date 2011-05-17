!------------------------------------------------------------------------------
! cswt_swav_mod -- CSWT library swav class
!
!! Provides functionality to support and manipulate a spherical wavelet.
!! The spherical wavelet data structure includes the wavelet sky plus its
!! parameters.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module cswt_swav_mod

  use cswt_error_mod
  use s2_types_mod
  use s2_sky_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------
  
  public :: &
    cswt_swav_init, &
    cswt_swav_free, &
    cswt_swav_dilate, &
    cswt_swav_rotate, &
    cswt_swav_azimuthal_bl, &
    cswt_swav_compute_alm, &
    cswt_swav_write_map_file, &
    cswt_swav_get_init, &
    cswt_swav_get_sky, &
    cswt_swav_get_sky_alm, &
    cswt_swav_get_map_status, &
    cswt_swav_get_alm_status, &
    cswt_swav_get_nside, &
    cswt_swav_get_pix_scheme, &
    cswt_swav_get_name, &
    cswt_swav_get_admiss, &
    cswt_swav_get_admiss_perc, &
    cswt_swav_get_dilation, &
    cswt_swav_get_alpha, &
    cswt_swav_get_beta, &
    cswt_swav_get_gamma


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface cswt_swav_init
    module procedure &
      cswt_swav_init_fun, &
      cswt_swav_init_sky, &
      cswt_swav_init_file, &
      cswt_swav_init_copy
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Default healpix pixelisation scheme.
  integer, public, parameter :: CSWT_SWAV_DEFAULT_PIX_SCHEME = S2_SKY_RING

  !! Default wavelet name if not specified.
  character(len=*), public, parameter :: &
     CSWT_SWAV_DEFAULT_NAME = 'Not specified'

  !! Dimension of dilation.
  integer, public, parameter :: CSWT_SWAV_DILATION_DIM = 2


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: cswt_swav
     private
     logical :: init = .false.
     type(s2_sky) :: sky
     character(len=S2_STRING_LEN) :: name = trim(CSWT_SWAV_DEFAULT_NAME)
     real(s2_sp) :: admiss = 0.0e0
     real(s2_sp) :: dilation(CSWT_SWAV_DILATION_DIM) = (/ 1.0e0, 1.0e0 /)
     real(s2_sp) :: alpha = 0.0e0
     real(s2_sp) :: beta = 0.0e0
     real(s2_sp) :: gamma = 0.0e0
  end type cswt_swav


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! cswt_swav_init_fun
    !
    !! Initialise a swav (spherical wavelet) from a template function defined
    !! either on the sphere or plane.
    !!
    !! Variables:
    !!   - fun: Function `pointer' to a function defined on either the sphere 
    !!     or the plane.
    !!   - fun_type_in: Integer specifing the function type (whether function 
    !!     defined directly on sphere or whether defined on plane and to be 
    !!     numericall projected onto the sphere.
    !!   - nside: Healpix sky resolution.
    !!   - [pix_scheme_in]: Pixelisation scheme to create sky in (if not
    !!     specified then CSWT_SWAV_DEFAULT_PIX_SCHEME is used).
    !!   - [lmax]: Alm maximum l.
    !!   - [mmax]: Alm maximum m.
    !!   - [param]: Parameter array specifying analytic parameters for the
    !!     function to be evaluated.
    !!   - [name]: String specifying name of wavelet.
    !!   - [dilation]: Dilation to scale template by.
    !!   - [alpha]: Alpha Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [beta]: Beta Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [gamma]: Gamma Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [norm_preserve]: Logical to specify whether or not a norm 
    !!     preserving dilation is to be performed.
    !!   - swav: Initialised swav.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function cswt_swav_init_fun(fun, nside, pix_scheme_in, lmax, mmax, &
      param, name, dilation, alpha, beta, gamma, fun_type, norm_preserve) &
      result(swav)

      interface 
         function fun(x, phi, param) result(val)
           use s2_types_mod
           real(s2_sp), intent(in) :: x, phi
           real(s2_sp), intent(in), optional :: param(:)
           real(s2_sp) :: val
         end function fun
      end interface
      integer, intent(in) :: nside
      integer, intent(in), optional :: pix_scheme_in, lmax, mmax
      real(s2_sp), intent(in), optional :: param(:)
      character(len=*), intent(in), optional :: name
      real(s2_sp), intent(in), optional :: dilation(CSWT_SWAV_DILATION_DIM)
      real(s2_sp), intent(in), optional :: alpha, beta, gamma
      integer, intent(in), optional :: fun_type
      logical, intent(in), optional :: norm_preserve
      type(cswt_swav) :: swav

      integer :: pix_scheme

      ! Check object not already initialised.
      if(swav%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_swav_init_fun')
        return
      end if

      ! Set healpix pixelisation scheme.
      if(present(pix_scheme_in)) then
         pix_scheme = pix_scheme_in
      else
         pix_scheme = CSWT_SWAV_DEFAULT_PIX_SCHEME
      end if

      ! Initialise sky.
      swav%sky = s2_sky_init(fun, nside, pix_scheme, lmax, mmax, param, &
        fun_type)

      ! Set object attributes.
      if(present(name)) swav%name = name
      ! Dilation values and Euler angles will be set by respective
      ! dilation and rotation routines.

      ! Set object as initialised.
      ! Must do this before call rotate and dilate routines since they 
      ! check that the swav object is already initialised.
      swav%init = .true.

      ! Perform dilation if required.
      if(present(dilation)) then
         call cswt_swav_dilate(swav, dilation, norm_preserve)  
      end if

      ! Perform rotation if required.
      if(present(alpha) .and. present(beta) .and. present(gamma)) then
         call cswt_swav_rotate(swav, alpha, beta, gamma)
      end if

      ! Calculate wavelet admissibility.
      swav%admiss = s2_sky_admiss(swav%sky)

    end function cswt_swav_init_fun


    !--------------------------------------------------------------------------
    ! cswt_swav_init_sky
    !
    !! Initialise a swav (spherical wavelet) from a predefined sky.
    !!
    !! Variables:
    !!   - sky: Function defined on the sky(/sphere) to initialise the
    !!     spherical wavelet.
    !!   - [name]: String specifying name of wavelet.
    !!   - [dilation]: Dilation to scale template by.
    !!   - [alpha]: Alpha Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [beta]: Beta Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [gamma]: Gamma Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [norm_preserve]: Logical to specify whether or not a norm 
    !!     preserving dilation is to be performed.
    !!   - swav: Initialised swav.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
     
    function cswt_swav_init_sky(sky, name, dilation, alpha, beta, &
      gamma, norm_preserve) result(swav)

      type(s2_sky), intent(in) :: sky
      character(len=*), intent(in), optional :: name
      real(s2_sp), intent(in), optional :: dilation(CSWT_SWAV_DILATION_DIM)
      real(s2_sp), intent(in), optional :: alpha, beta, gamma
      logical, intent(in), optional :: norm_preserve
      type(cswt_swav) :: swav

     ! Check object not already initialised.
      if(swav%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_swav_init_sky')
        return
      end if

      ! Initialise wavelet sky.
      swav%sky = s2_sky_init(sky)

      ! Set object attributes.
      if(present(name)) swav%name = name
      ! Dilation values and Euler angles will be set by respective
      ! dilation and rotation routines.

      ! Set object as initialised.
      ! Must do this before call rotate and dilate routines since they 
      ! check that the swav object is already initialised.
      swav%init = .true.

      ! Perform dilation if required.
      if(present(dilation)) then
         call cswt_swav_dilate(swav, dilation, norm_preserve)  
      end if

      ! Perform rotation if required.
      if(present(alpha) .and. present(beta) .and. present(gamma)) then
         call cswt_swav_rotate(swav, alpha, beta, gamma)
      end if

      ! Calculate wavelet admissibility.
      ! Extension 12/04/1005: Only compute admissibility if map is present.
      ! If wavelet initialised from file and defined as alms only then
      ! not worth computing map just to compute numerical admissibility.
      ! Besides, admissibility may be computed directly elsewhere if check required.
      if(s2_sky_get_map_status(swav%sky)) then
         swav%admiss = s2_sky_admiss(swav%sky)
      else
         swav%admiss = -9.0e0  ! Set as -9 to indicate not calculated.
                               ! Valid range is [-1,1] so obviously invalid.
      end if

    end function cswt_swav_init_sky


    !--------------------------------------------------------------------------
    ! cswt_swav_init_file
    !
    !! Initialise a swav (spherical wavelet) from a sky read from an input 
    !! file.
    !!
    !! Variables:
    !!   - filename: Name of input file containing the function defined on 
    !!     the sky that is used to initialise the spherical wavelet.
    !!   - [name]: String specifying name of wavelet.
    !!   - [dilation]: Dilation to scale template by.
    !!   - [alpha]: Alpha Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [beta]: Beta Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [gamma]: Gamma Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - [norm_preserve]: Logical to specify whether or not a norm 
    !!     preserving dilation is to be performed.
    !!   - swav: Initialised swav.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
     
    function cswt_swav_init_file(filename, name, dilation, alpha, beta, &
      gamma, norm_preserve) result(swav)

      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: name
      real(s2_sp), intent(in), optional :: dilation(CSWT_SWAV_DILATION_DIM)
      real(s2_sp), intent(in), optional :: alpha, beta, gamma
      logical, intent(in), optional :: norm_preserve
      type(cswt_swav) :: swav

      type(s2_sky) :: temp_sky

     ! Check object not already initialised.
      if(swav%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_swav_init_file')
        return
      end if

      temp_sky = s2_sky_init(filename, S2_SKY_FILE_TYPE_MAP)

      swav = cswt_swav_init_sky(temp_sky, name=name, dilation=dilation, &
        alpha=alpha, beta=beta, gamma=gamma, norm_preserve=norm_preserve)

      call s2_sky_free(temp_sky)

    end function cswt_swav_init_file


    !--------------------------------------------------------------------------
    ! cswt_swav_init_copy
    !
    !! Initialise a new swav (spherical wavelet) as a copy of another swav.
    !!
    !! Variables:
    !!   - orig: The original swav to copy.
    !!   - copy: The initialised swav copied from the original.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
     
    function cswt_swav_init_copy(orig) result(copy)

      type(cswt_swav), intent(in) :: orig
      type(cswt_swav) :: copy

      ! Check object initialised.
      if(.not. orig%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_init_copy')
      end if 

      ! Copy attributes.
      copy%sky = s2_sky_init(orig%sky)
      copy%name = orig%name
      copy%admiss = orig%admiss
      copy%dilation = orig%dilation
      copy%alpha = orig%alpha
      copy%beta = orig%beta
      copy%gamma = orig%gamma
      
      ! Set object as initialised.
      copy%init = .true.

    end function cswt_swav_init_copy


    !--------------------------------------------------------------------------
    ! cswt_swav_free
    !
    !! Free all data associated with an initialised swav and reset all other 
    !! attributes.
    !!
    !! Variables:
    !!   - swav: The swacv to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine cswt_swav_free(swav)

      type(cswt_swav), intent(inout) :: swav

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_free')
      end if 

      ! Free sky.
      call s2_sky_free(swav%sky)

      ! Reset all other variables.
      swav%name = trim(CSWT_SWAV_DEFAULT_NAME)
      swav%admiss = 0.0e0
      swav%dilation = 1.0e0
      swav%alpha = 0.0e0
      swav%beta = 0.0e0
      swav%gamma = 0.0e0

      swav%init = .false.

    end subroutine cswt_swav_free


    !--------------------------------------------------------------------------
    ! cswt_swav_dilate
    !
    !! Dilate the spherical wavelet and update the dilation parameter.  Only
    !! a mother wavelet with dilation of one may be dilated (otherwise a 
    !! different dilation would have to be performed to give the correct 
    !! final dilation).
    !!
    !! Notes:
    !!   - If dilation value is set to one then no dilation performed and 
    !!     routine simply exists.  This is useful for the case where only the
    !!     wavelet alms are defined and not the map.  Then routine may be 
    !!     called as normal but nothing will happen so easy to generalise
    !!     parent code.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet dilated.
    !!   - dilation: Dilation value.
    !!   - [norm_preserve_in]: Logical to specify whether or not a norm
    !!     preserving dilation is to be performed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine cswt_swav_dilate(swav, dilation, norm_preserve_in)
      
      type(cswt_swav), intent(inout) :: swav
      real(s2_sp), intent(in) :: dilation(CSWT_SWAV_DILATION_DIM)
      logical, intent(in), optional :: norm_preserve_in

      real(s2_sp), parameter :: TOL = 1e-5
      logical, parameter :: DEFAULT_NORM_PRESERVE = .true.

      logical :: norm_preserve

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_dilate')
      end if 

      ! Check wavelet not rotated.
      if( abs(swav%alpha)>TOL .or. abs(swav%beta)>TOL &
        .or. abs(swav%gamma)>TOL) then
         call cswt_error(CSWT_ERROR_SWAV_DIL_ERROR, 'cswt_swav_dilate', &
           comment_add='Wavelet already rotated before attempt to dilate')
      end if

      ! Check dilation originally one.
      if( sum(abs(swav%dilation - 1.0e0)) > TOL ) then
         call cswt_error(CSWT_ERROR_SWAV_DIL_ERROR, 'cswt_swav_dilate', &
           comment_add='Wavelet already dilated (can only dilate mother wavelet with original dilation of 1.0e0.')
      end if

      ! If new dilation is one then do nothing.
      if( sum(abs(dilation - 1.0e0)) <= TOL ) then
         return
      end if

      ! Set norm_preserve.
      if(present(norm_preserve_in)) then
         norm_preserve = norm_preserve_in
      else
         norm_preserve = DEFAULT_NORM_PRESERVE
      end if

      ! Dilate.
      call s2_sky_dilate(swav%sky, dilation(1), dilation(2), norm_preserve)

      swav%dilation = dilation

    end subroutine cswt_swav_dilate


    !--------------------------------------------------------------------------
    ! cswt_swav_rotate
    !
    !! Rotate a spherical wavelet.  Only an originally non-rotated wavelet may
    !! be rotated (otherwise a different rotation would have to be performed 
    !! to give the correct final rotation).
    !!
    !! Variables:
    !!   - swav: Spherical wavelet rotated.
    !!   - alpha: Alpha Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - beta: Beta Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !!   - gamma: Gamma Euler rotation angle applied to rotate the 
    !!     funciton defined on the sphere.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine cswt_swav_rotate(swav, alpha, beta, gamma)

      type(cswt_swav), intent(inout) :: swav
      real(s2_sp), intent(in) :: alpha, beta, gamma

      real(s2_sp), parameter :: TOL = 1e-5

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_rotate')
      end if 

      ! Check wavelet not already rotated.
      if( abs(swav%alpha)>TOL .or. abs(swav%beta)>TOL &
        .or. abs(swav%gamma)>TOL) then
         call cswt_error(CSWT_ERROR_SWAV_ROT_ERROR, 'cswt_swav_dilate', &
           comment_add='Wavelet already rotated (can only rotate mother wavelet with original 0.0e0 rotation.')
      end if

      call s2_sky_rotate(swav%sky, alpha, beta, gamma)

      swav%alpha = alpha
      swav%beta = beta
      swav%gamma = gamma

    end subroutine cswt_swav_rotate


    !--------------------------------------------------------------------------
    ! cswt_swav_azimuthal_bl
    !
    !! Find azimuthal band limit of wavelet.  Finds lowest m' value such that 
    !! cutoff_prop*100 percent of the cm power is contained in the alms 
    !! with m index below m'.
    !!
    !! Notes:
    !!   - Only successful if alms for the sky are already computed.  Cannot 
    !      automatically compute alms if not already computed since lmax and 
    !!     mmax may not be defined.
    !!
    !! Variables
    !!  - swav: Sky to find azimuthal band limit for.
    !!  - [cutoff_prop]: Cutoff proportion to use when finding azimuthal 
    !!    band limit.  If not specified default value is used.
    !!  - mmax_min: The azimuthal band limit found.
    !
    !! @author J. D. McEwen
    !! @version 0.1 May 2005
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_swav_azimuthal_bl(swav, cutoff_prop) result(mmax_min)

      type(cswt_swav), intent(inout) :: swav
      real(s2_sp), intent(in), optional :: cutoff_prop
      integer :: mmax_min

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_azimuthal_bl')
      end if 

      ! Find minimim band limit.
      mmax_min = s2_sky_azimuthal_bl(swav%sky, cutoff_prop)

    end function cswt_swav_azimuthal_bl


    !--------------------------------------------------------------------------
    ! cswt_swav_compute_alm
    !
    !! Compute the spherical wavelet sky alms.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet to compute alms of.
    !!   - lmax: Spherical harmonic lmax to compute alms up to.
    !!   - mmax: Spherical harmoinc mmax to compute alms up to.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine cswt_swav_compute_alm(swav, lmax, mmax, message)

      type(cswt_swav), intent(inout) :: swav
      integer, intent(in) :: lmax, mmax
      logical, intent(in), optional :: message

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_compute_alm')
      end if 

      call s2_sky_compute_alm(swav%sky, lmax, mmax, message)

    end subroutine cswt_swav_compute_alm


    !--------------------------------------------------------------------------
    ! cswt_swav_write_map_file
    !
    !! Write the spherical wavelet sky to an output file.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet containing sky to write to file.
    !!   - filename: Name of the output fits file.
    !!   - [comment]: Optional additional comment to be added to the fits 
    !!     file header.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine cswt_swav_write_map_file(swav, filename, comment)

      type(cswt_swav), intent(in) :: swav
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_write_map_file')
      end if 

      call s2_sky_write_map_file(swav%sky, filename, comment)

    end subroutine cswt_swav_write_map_file


    !--------------------------------------------------------------------------
    ! Get routines.
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_swav_get_init
    !
    !! Get init variable from the passed swav.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - init: Copy of object init variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_init(swav) result(init)

      type(cswt_swav), intent(in) :: swav
      logical :: init

      init = swav%init

    end function cswt_swav_get_init


    !--------------------------------------------------------------------------
    ! cswt_swav_get_sky
    !
    !! Get sky variable from the passed swav.
    !!
    !! Notes:
    !!   - Initialises a new sky as a copy of the swav sky.
    !!   - The returned sky is subsequently independed of the sky stored
    !!     herein and should be freed by the calling routine at some point.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - sky: Copy of sky extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_sky(swav) result(sky)

      type(cswt_swav), intent(in) :: swav
      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_sky')
      end if 

      ! Make a copy for the returned sky.
      sky = s2_sky_init(swav%sky)

    end function cswt_swav_get_sky

    
    !--------------------------------------------------------------------------
    ! cswt_swav_get_sky_alm
    !
    !! Get alm of the swav sky.
    !!
    !! Variables:
    !!   - swav: pherical wavelet object to get copy of variable of.
    !!   - alm(:,:): Object alm variable returned.  Space must
    !!     already be allocated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 July 2005
    !
    ! Revisions:
    !   July 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    subroutine cswt_swav_get_sky_alm(swav, alm)

      type(cswt_swav), intent(in) :: swav
      complex(s2_spc), intent(out) :: alm(:,:)

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_sky')
      end if 

      ! Get wavelet sky alm.
      call s2_sky_get_alm(swav%sky, alm)

    end subroutine cswt_swav_get_sky_alm

    
    !--------------------------------------------------------------------------
    ! cswt_swav_get_map_status
    !
    !! Get map_status of the sky of the passed swav.  Note that the map 
    !! status is stored within the wavelet sky object.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get map_status of.
    !!   - map_status: Map status of the wavelet sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_map_status(swav) result(map_status)

      type(cswt_swav), intent(in) :: swav
      logical :: map_status

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_map_status')
      end if 

      ! Get map status from the sky.
      map_status = s2_sky_get_map_status(swav%sky)

    end function cswt_swav_get_map_status


    !--------------------------------------------------------------------------
    ! cswt_swav_get_alm_status
    !
    !! Get alm_status of the sky of the passed swav.  Note that the alm 
    !! status is stored within the wavelet sky object.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get alm_status of.
    !!   - alm_status: Alm status of the wavelet sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_alm_status(swav) result(alm_status)

      type(cswt_swav), intent(in) :: swav
      logical :: alm_status

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_alm_status')
      end if 

      ! Get map status from the sky.
      alm_status = s2_sky_get_alm_status(swav%sky)

    end function cswt_swav_get_alm_status


    !--------------------------------------------------------------------------
    ! cswt_swav_get_nside
    !
    !! Get nside of the sky of the passed swav.  Note that the alm 
    !! status is stored within the wavelet sky object.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get nside of.
    !!   - nside: Nside of the wavelet sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_nside(swav) result(nside)

      type(cswt_swav), intent(in) :: swav
      integer :: nside

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_nside')
      end if 

      ! Get map status from the sky.
      nside = s2_sky_get_nside(swav%sky)

    end function cswt_swav_get_nside


    !--------------------------------------------------------------------------
    ! cswt_swav_get_pix_scheme
    !
    !! Get pix_scheme of the sky of the passed swav.  Note that the alm 
    !! status is stored within the wavelet sky object.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get pix_scheme of.
    !!   - pix_scheme: Pix_scheme of the wavelet sky.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_pix_scheme(swav) result(pix_scheme)

      type(cswt_swav), intent(in) :: swav
      integer :: pix_scheme

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_pix_scheme')
      end if 

      ! Get map status from the sky.
      pix_scheme = s2_sky_get_pix_scheme(swav%sky)

    end function cswt_swav_get_pix_scheme


    !--------------------------------------------------------------------------
    ! cswt_swav_get_name
    !
    !! Get name variable from the passed swav.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - name: Copy of object name variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_name(swav) result(name)

      type(cswt_swav), intent(in) :: swav
      character(len=S2_STRING_LEN) :: name

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_name')
      end if 

      name = swav%name

    end function cswt_swav_get_name


    !--------------------------------------------------------------------------
    ! cswt_swav_get_admiss
    !
    !! Get normalised numerical admissibility of the passed swav.
    !!
    !! Variables: 
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - admiss: Copy of object admiss varaible extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_admiss(swav) result(admiss)

      type(cswt_swav), intent(in) :: swav
      real(s2_sp) :: admiss

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_write_map_file')
      end if 

      admiss = swav%admiss

    end function cswt_swav_get_admiss


    !--------------------------------------------------------------------------
    ! cswt_swav_get_admiss_perc
    !
    !! Get admissibility percentage of the passed swav.
    !!
    !! Notes:
    !!   - Only the normalised admissibility value is stored as a spherical
    !!     wavelet attribute.  The admissibility percentage is calculated
    !!     from this when required.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet to calculate admissibility percentage of.
    !!   - admiss_perc: Admissibility percentage calculated.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_admiss_perc(swav) result(admiss_perc)

      type(cswt_swav), intent(in) :: swav
      real(s2_sp) :: admiss_perc
 
      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_write_map_file')
      end if 

      admiss_perc = (1 - abs(swav%admiss)) * 100e0

    end function cswt_swav_get_admiss_perc


    !--------------------------------------------------------------------------
    ! cswt_swav_get_dilation
    !
    !! Get dilation variable from the passed swav.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - dilation: Copy of object dilation variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_dilation(swav) result(dilation)

      type(cswt_swav), intent(in) :: swav
      real(s2_sp) :: dilation(CSWT_SWAV_DILATION_DIM)

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_dilation')
      end if 

      dilation = swav%dilation

    end function cswt_swav_get_dilation


    !--------------------------------------------------------------------------
    ! cswt_swav_get_alpha
    !
    !! Get alpha Euler angle from the passed swav.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - alpha: Copy of object alpha Euler angle variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_alpha(swav) result(alpha)

      type(cswt_swav), intent(in) :: swav
      real(s2_sp) :: alpha

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_alpha')
      end if 

      alpha = swav%alpha

    end function cswt_swav_get_alpha
    

    !--------------------------------------------------------------------------
    ! cswt_swav_get_beta
    !
    !! Get beta Euler angle from the passed swav.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - beta: Copy of object beta Euler angle variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_beta(swav) result(beta)

      type(cswt_swav), intent(in) :: swav
      real(s2_sp) :: beta

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_beta')
      end if 

      beta = swav%beta

    end function cswt_swav_get_beta


    !--------------------------------------------------------------------------
    ! cswt_swav_get_gamma
    !
    !! Get gamma Euler angle from the passed swav.
    !!
    !! Variables:
    !!   - swav: Spherical wavelet object to get copy of variable of.
    !!   - gamma: Copy of object gamma Euler angle variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_swav_get_gamma(swav) result(gamma)

      type(cswt_swav), intent(in) :: swav
      real(s2_sp) :: gamma

      ! Check object initialised.
      if(.not. swav%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_swav_get_gamma')
      end if 

      gamma = swav%gamma

    end function cswt_swav_get_gamma


end module cswt_swav_mod
