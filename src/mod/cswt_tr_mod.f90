!------------------------------------------------------------------------------
! cswt_tr_mod -- CSWT library tr class
!
!! Provides functionality to support and compute the directional spherical 
!! wavelet transform.  In addition the mother wavelet, the dilation values 
!! and all wavelet parameters are stored herein.  The resulting wavelet 
!! coefficients for each dilation are sampled on an ecp (equi-sampled) 
!! alpha-beta-gamma Euler angle grid.  Functionality is also provided to
!! convert the alpha-beta dimensions to Healpix sky representations and
!! vice versa.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

module cswt_tr_mod

  use s2_types_mod
  use s2_sky_mod
  use s2_vect_mod
  use cswt_error_mod
  use cswt_swav_mod

  implicit none

  private 


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    cswt_tr_init, &
    cswt_tr_free, &
    cswt_tr_analysis, &
    cswt_tr_norm, &
    cswt_tr_multiply, &
    cswt_tr_const_multiply, &
    cswt_tr_const_subtract, &
    cswt_tr_wcoeff_sigma, &
    cswt_tr_localmax_ab, &
    cswt_tr_wcoeff_thres_nsigma, &
    cswt_tr_mask_nonzero, &
    cswt_tr_mask_nonzero_weight, &
    cswt_tr_mask_invert, &
    cswt_tr_mask_apply, &
    cswt_tr_mask_thres_nsigma, &
    cswt_tr_mask_copy, &
    cswt_tr_mask_gen, &
    cswt_tr_smooth_gaussian, &
    cswt_tr_extract_wcoeff_sky, &
    cswt_tr_io_fits_write_wcoeff_sky, &
    cswt_tr_io_fits_write_wcoeff, &
    cswt_tr_io_txt_dilation_write, &
    cswt_tr_io_txt_dilation_read, &
    cswt_tr_get_init, &
    cswt_tr_get_swav_mother, &
    cswt_tr_get_wcoeff, &
    cswt_tr_get_dilation, &
    cswt_tr_get_lmax, &
    cswt_tr_get_mmax, &
    cswt_tr_get_n_alpha, &
    cswt_tr_get_n_beta, &
    cswt_tr_get_n_gamma, &
    cswt_tr_get_n_dilation, &
    cswt_tr_get_wcoeff_status, &
    cswt_tr_get_swav_mother_status


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface cswt_tr_init
     module procedure &
        cswt_tr_init_data, &
        cswt_tr_init_file_dil, &
        cswt_tr_init_file, &
        cswt_tr_init_read, &
        cswt_tr_init_copy, &
        cswt_tr_init_merg
  end interface

  interface cswt_tr_const_multiply
     module procedure &
        cswt_tr_const_multiply_array, &
        cswt_tr_const_multiply_val
  end interface
 
  interface cswt_tr_const_subtract
     module procedure &
        cswt_tr_const_subtract_array, &
        cswt_tr_const_subtract_val
  end interface

  interface cswt_tr_wcoeff_sigma
     module procedure &
        cswt_tr_wcoeff_sigma_each, &
        cswt_tr_wcoeff_sigma_all
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Radian dilation unit specifier.
  integer, public, parameter :: CSWT_TR_DIL_UNIT_RAD = 1

  !! Arcminute dilation unit specifier.
  integer, public, parameter :: CSWT_TR_DIL_UNIT_ARCMIN = 2

  ! Dilation unit strings (not public). 
  character(len=*), parameter :: CSWT_TR_DIL_UNIT_STR_RAD = 'radian'
  character(len=*), parameter :: CSWT_TR_DIL_UNIT_STR_ARCMIN = 'arcmin'

  ! Spherical wavelet analysis methods.
  !! Fast FFT analysis method specifier.
  integer, public, parameter :: CSWT_TR_METHOD_FAST_FFT_REAL = 1
  !! Fast FFT analysis method specifier.
  integer, public, parameter :: CSWT_TR_METHOD_FAST_FFT = 2
  !! Fast DFT analysis method specifier.
  integer, public, parameter :: CSWT_TR_METHOD_FAST_DFT = 3
  !! Fast isotropic analysis method specifier.
  integer, public, parameter :: CSWT_TR_METHOD_FAST_ISOTROPIC = 4
  !! Direct analysis method specifier.
  integer, public, parameter :: CSWT_TR_METHOD_DIRECT = 5

  ! Kernel types for morpholocial operations.
  !! Square kernel type specifier.
  integer, public, parameter :: CSWT_TR_KERNEL_TYPE_SQR = 1
  !! Circular kernel type specifier.
  integer, public, parameter :: CSWT_TR_KERNEL_TYPE_CIRC = 2

  ! Extended coefficient mask generation methods.
  !! Generate mask from wavelet coefficients specifier.
  integer, public, parameter :: CSWT_TR_MASK_GEN_CSWT = 1
  !! Generate mask by morphological operations specifier.
  integer, public, parameter :: CSWT_TR_MASK_GEN_MORPH = 2

  ! Nsigma threshold mask generation methods.
  !! Generate thresholded nsigma mask by considering absolte value.
  integer, public, parameter :: CSWT_TR_THRES_NSIGMA_MODE_ABS = 1
  !! Generate thresholded nsigma mask by leaving values above threshold.
  integer, public, parameter :: CSWT_TR_THRES_NSIGMA_MODE_ABOVE = 2
  !! Generate thresholded nsigma mask by leaving values below threshold.
  integer, public, parameter :: CSWT_TR_THRES_NSIGMA_MODE_BELOW = 3


  !---------------------------------------
  ! Data types
  !---------------------------------------

  type, public :: cswt_tr
     private
     logical :: init = .false.
     type(cswt_swav) :: swav_mother
     real(s2_sp), allocatable :: wcoeff(:,:,:,:)
     real(s2_sp), allocatable :: dilation(:,:)
     integer :: lmax = 0
     integer :: mmax = 0
     integer :: n_alpha = 0
     integer :: n_beta = 0
     integer :: n_gamma = 0
     integer :: n_dilation = 0
     logical :: wcoeff_status = .false.
     logical :: swav_mother_status = .false.  ! If read from file swav_mother
                                              ! will not be present.
  end type cswt_tr


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! cswt_tr_init_data
    !
    !! Initialise a tr structure from a mother wavelet, dilations and 
    !! resolution parameters.
    !!
    !! Notes:
    !!   - Space for wcoeff is allocated when tr is initialised (before 
    !!     actually calculate coefficients).
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !!   - n_gamma must be odd for the fast cswt.
    !!
    !! Variables:
    !!   - swav_mother: Mother wavelet.
    !!   - dilation: Dilation values.
    !!   - lmax: Alm lmax.
    !!   - mmax: Alm mmax.
    !!   - n_gamma: Number of gamma Euler orientation to consider in the 
    !!     directional transform.  
    !!   - tr: Returned tr structure initialised with a mother wavelet, 
    !!     dilations and resolution parameters.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function cswt_tr_init_data(swav_mother, dilation, lmax, mmax, n_gamma) &
      result(tr)

      type(cswt_swav), intent(in) :: swav_mother
      real(s2_sp), intent(in) :: dilation(:,:)
      integer, intent(in) :: lmax, mmax
      integer, intent(in) :: n_gamma
      type(cswt_tr) :: tr

      integer :: fail
      
      ! Check object not already initialised.
      if(tr%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_init_data')
        return
      end if

      ! Initialise object attributes.
      tr%swav_mother = cswt_swav_init(swav_mother)
      tr%swav_mother_status = .true.
      tr%n_dilation = size(dilation, 1)
      tr%lmax = lmax
      tr%mmax = mmax
      allocate(tr%dilation(1:tr%n_dilation, CSWT_SWAV_DILATION_DIM), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_init_data')
      end if
      tr%dilation = dilation
      tr%n_alpha = 2 * tr%lmax + 1
      tr%n_beta  = tr%lmax + 1
      tr%n_gamma = n_gamma

      ! Check n_gamma is odd.
      if(mod(tr%n_gamma,2) == 0) then
         call cswt_error(CSWT_ERROR_TR_NGAMMA_INVALID, 'cswt_tr_init_data', &
           comment_add='Number of orientations is even and must be odd')
      end if

      ! Allocate space for wavelet coefficients.
      allocate(tr%wcoeff(1:tr%n_dilation, &
        0:tr%n_alpha-1, 0:tr%n_beta-1, 0:tr%n_gamma-1), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_init_data')
      end if
      ! Initialise to zero.
      tr%wcoeff = 0.0e0

      tr%init = .true.

    end function cswt_tr_init_data


    !--------------------------------------------------------------------------
    ! cswt_tr_init_file_dil
    !
    !! Initialise a tr strucutre from a mother wavelet, dilations read from
    !! an input dilation text file and resolution parameters.
    !!
    !! Notes:
    !!   - Once dilations are read cswt_tr_init_data is called.
    !!   - Space for wcoeff is allocated when tr is initialised (before 
    !!     actually calculate coefficients).
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !!   - n_gamma must be odd for the fast cswt.
    !!
    !! Variables:
    !!   - swav_mother: Mother wavelet.
    !!   - filename_dil: Name of input dilation file.
    !!   - lmax: Alm lmax.
    !!   - mmax: Alm mmax.
    !!   - n_gamma: Number of gamma Euler orientation to consider in the 
    !!     directional transform.  
    !!   - tr: Returned tr structure initialised with a mother wavelet, 
    !!     dilations and resolution parameters.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_init_file_dil(swav_mother, filename_dil, lmax, mmax, &
      n_gamma) result(tr)

      type(cswt_swav), intent(in) :: swav_mother
      character(len=*), intent(in) :: filename_dil
      integer, intent(in) :: lmax, mmax
      integer, intent(in) :: n_gamma
      type(cswt_tr) :: tr

      real(s2_sp), allocatable :: dilation(:,:)

      ! Check object not already initialised.
      if(tr%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_init_file_dil')
        return
      end if

      ! Read dilations from file.
      call cswt_tr_io_txt_dilation_read(filename_dil, dilation)

      ! Construct tr data structure from data.
      tr = cswt_tr_init_data(swav_mother, dilation, lmax, mmax, n_gamma)

      ! Free temporary dilation storage array.
      deallocate(dilation)

    end function cswt_tr_init_file_dil


    !--------------------------------------------------------------------------
    ! cswt_tr_init_file
    !
    !! Initialise a tr strucutre from a mother wavelet read from an input 
    !! fits file, dilations read from an input dilation text file and 
    !! resolution parameters.
    !!
    !! Notes:
    !!   - Once the wavelet is read cswt_tr_init_dil is called.
    !!   - Space for wcoeff is allocated when tr is initialised (before 
    !!     actually calculate coefficients).
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !!   - n_gamma must be odd for the fast cswt.
    !!
    !! Variables:
    !!   - filename_swav: Name of input fits sky file containing the mother
    !!     wavelet.
    !!   - filename_dil: Name of input dilation file.
    !!   - lmax: Alm lmax.
    !!   - mmax: Alm mmax.
    !!   - n_gamma: Number of gamma Euler orientation to consider in the 
    !!     directional transform.  
    !!   - [wav_name]: Name of mother wavelet type.
    !!   - tr: Returned tr structure initialised with a mother wavelet, 
    !!     dilations and resolution parameters.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_init_file(filename_swav, filename_dil, lmax, mmax, &
      n_gamma, wav_name) result(tr)
 
      character(len=*), intent(in) :: filename_swav
      character(len=*), intent(in) :: filename_dil
      integer, intent(in) :: lmax, mmax
      integer, intent(in) :: n_gamma
      character(len=*), intent(in), optional :: wav_name
      type(cswt_tr) :: tr

      type(cswt_swav) :: temp_swav

      ! Check object not already initialised.
      if(tr%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_init_file')
        return
      end if

      ! Initialise temporary spherical wavelet from file.
      temp_swav = cswt_swav_init(filename_swav, name=wav_name)

      ! Construct tr data structure from swav and dilations contained in file.
      tr = cswt_tr_init_file_dil(temp_swav, filename_dil, lmax, mmax, n_gamma)

      ! Free temporary swav (copy made in tr data structure).
      call cswt_swav_free(temp_swav)

    end function cswt_tr_init_file


    !--------------------------------------------------------------------------
    ! cswt_tr_init_read
    !
    !! Wrapper to initalise a tr data structure from a file.  The tr 
    !! structure is read and initialised by the routine 
    !! cswt_tr_io_fits_read_wcoeff.
    !!
    !! Notes: 
    !!   - If a tr structure is initialised in this manner then no wavelet
    !!     is stored and further analysis calls may not be made.
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !!
    !! Variables:
    !!   - filename: Name of wavelet coefficient fits file containing the tr 
    !!     data to be read.
    !!   - tr: Returned tr structure initialised with the data contained in
    !!     the input fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_init_read(filename) result(tr)

      character(len=*), intent(in) :: filename
      type(cswt_tr) :: tr
    
      ! Check object not already initialised.
      if(tr%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_init_file')
        return
      end if

      call cswt_tr_io_fits_read_wcoeff(filename, tr)

      ! Set as initialised.
      tr%init = .true.

    end function cswt_tr_init_read


    !--------------------------------------------------------------------------
    ! cswt_tr_init_computed
    !
    !! Initialise a tr structure from precomputed wavelet coefficients
    !! (if, say, read from file).
    !!
    !! Notes:
    !!   - If a tr structure is initialised in this manner then no wavelet
    !!     is stored and further analysis calls may not be made.
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !!
    !! Variables:
    !!   - dilation: Array of dilation values.
    !!   - wcoeff: Array of wavelet coefficients.
    !!   - lmax: Alm lmax.
    !!   - mmax: Alm mmax.
    !!   - tr: Returned tr strucutre initialised with the wavelet coefficients
    !!     and dilations passed as input.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_init_computed(dilation, wcoeff, lmax, mmax) result(tr)

      real(s2_sp), intent(in) :: dilation(:,:)
      real(s2_sp),intent(in) :: wcoeff(:,:,:,:)
      integer, intent(in) :: lmax, mmax
      type(cswt_tr) :: tr

      integer :: fail

      ! Check object not already initialised.
      if(tr%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_init_computed')
        return
      end if

      ! Initialise object attributes.
      tr%n_dilation = size(dilation, 1)
      tr%n_alpha = size(wcoeff,2)
      tr%n_beta  = size(wcoeff,3)
      tr%n_gamma = size(wcoeff,4)
      tr%lmax = lmax
      tr%mmax = mmax

      ! Check n_gamma is odd.
      if(mod(tr%n_gamma,2) == 0) then
         call cswt_error(CSWT_ERROR_TR_NGAMMA_INVALID, &
           'cswt_tr_init_computed', &
           comment_add='Number of orientations is even and must be odd')
      end if

      ! Check sizes consistent.
      if(tr%lmax /= (tr%n_alpha - 1) / 2 .or. tr%lmax /= tr%n_beta - 1) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_init_computed', &
           comment_add='Wavelet coefficient array size does not match lmax')
      end if

      ! Allocate space and save dilation.
      allocate(tr%dilation(1:tr%n_dilation, CSWT_SWAV_DILATION_DIM), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_init_computed')
      end if
      tr%dilation = dilation

      ! Allocate space and save wavelet coefficients.
      allocate(tr%wcoeff(1:tr%n_dilation, &
        0:tr%n_alpha-1, 0:tr%n_beta-1, 0:tr%n_gamma-1), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_init_computed')
      end if
      tr%wcoeff = wcoeff
      tr%wcoeff_status = .true.

      tr%swav_mother_status = .false.

      tr%init = .true.

    end function cswt_tr_init_computed


    !--------------------------------------------------------------------------
    ! cswt_tr_init_copy
    !
    !! Initialise a tr structure as a copy of another tr.
    !!
    !! Variables:
    !!   - orig: Original tr to be copied.
    !!   - copy: Returned tr strucutre initialised with a copy of the original.
    !!
    !! Notes:
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function cswt_tr_init_copy(orig) result(copy)

      type(cswt_tr), intent(in) :: orig
      type(cswt_tr) :: copy

      integer :: fail

      ! Check object initialised.
      if(.not. orig%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_init_copy')
      end if 

      ! Copy attributes.
      copy%lmax = orig%lmax
      copy%mmax = orig%mmax
      copy%n_alpha = orig%n_alpha
      copy%n_beta = orig%n_beta
      copy%n_gamma = orig%n_gamma
      copy%n_dilation = orig%n_dilation
      copy%wcoeff_status = orig%wcoeff_status
      copy%swav_mother_status = orig%swav_mother_status

      ! Allocate space for dilation and copy.
      allocate(copy%dilation(1:copy%n_dilation, CSWT_SWAV_DILATION_DIM), &
        stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_init_copy')
      end if  
      copy%dilation = orig%dilation

      ! Allocate space for wavelet coefficients and copy.
      ! (Even if not calculated will be allocated and zero in orig so copy.)
      allocate(copy%wcoeff(1:copy%n_dilation, &
        0:copy%n_alpha-1, 0:copy%n_beta-1, 0:copy%n_gamma-1), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_init_copy')
      end if
      copy%wcoeff = orig%wcoeff

      ! Copy swav_mother if exists.
      if(copy%swav_mother_status) then
         copy%swav_mother = cswt_swav_init(orig%swav_mother)
      end if

      ! Set object as initialised.
      copy%init = .true.

    end function cswt_tr_init_copy


    !--------------------------------------------------------------------------
    ! cswt_tr_init_merg
    !
    !! Initialise a tr structure by merging a number of tr structures that 
    !! have only a single dilation.
    !!
    !! Variables:
    !!   - tr: Array of tr objects to be merged.
    !!   - tr_merg: Returned merged tr strucutre initialised from the 
    !!     merged tr data.
    !!
    !! Notes:
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------
    
    function cswt_tr_init_merg(tr) result(tr_merg)

      type(cswt_tr), intent(in) :: tr(:)
      type(cswt_tr) :: tr_merg

      integer :: n_tr, i_tr, fail

      ! Get number of tr structures to merg.
      n_tr = size(tr)

      ! Check tr_merg not init.
      if(tr_merg%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_init_merg')
        return
      end if

      ! Check each tr structure.
      do i_tr = 1,n_tr

         ! Check each tr init.
         if(.not. tr(i_tr)%init) then
            call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_init_merg')
         end if

         ! Check each tr has only one dilation.
         if(tr(i_tr)%n_dilation /= 1) then
            call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_init_merg', &
              comment_add='Each tr has more than one dilation')
         end if

         ! Check each tr has same size and wcoeff_status attributes.
         if(tr(i_tr)%lmax /= tr(1)%lmax &
              .or. tr(i_tr)%mmax /= tr(1)%mmax &
              .or. tr(i_tr)%n_alpha /= tr(1)%n_alpha &
              .or. tr(i_tr)%n_beta /= tr(1)%n_beta &
              .or. tr(i_tr)%n_gamma /= tr(1)%n_gamma &
              .or. tr(i_tr)%wcoeff_status .neqv. tr(1)%wcoeff_status) then
            call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_init_merg', &
              comment_add='Tr objects have inconsistent attributes')
         end if

      end do

      ! Set attributes for merged object.
      tr_merg%lmax = tr(1)%lmax
      tr_merg%mmax = tr(1)%mmax
      tr_merg%n_alpha = tr(1)%n_alpha
      tr_merg%n_beta = tr(1)%n_beta
      tr_merg%n_gamma = tr(1)%n_gamma
      tr_merg%wcoeff_status = tr(1)%wcoeff_status
      tr_merg%n_dilation = n_tr
      tr_merg%swav_mother_status = .false.  ! Don't save mother wavelets.

      ! Allocate space for dilations.
      allocate(tr_merg%dilation(1:tr_merg%n_dilation, &
        CSWT_SWAV_DILATION_DIM), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_init_merg')
      end if  

      ! Save dilations.
      do i_tr = 1,n_tr
         ! Note, only one dilation in each tr (already checked).
         tr_merg%dilation(i_tr,:) = tr(i_tr)%dilation(1,:)  
      end do

      ! Copy wcoeff if present.
      if(tr_merg%wcoeff_status) then

         ! Allocate space for wcoeff.
         allocate(tr_merg%wcoeff(1:tr_merg%n_dilation, 0:tr_merg%n_alpha-1, &
           0:tr_merg%n_beta-1, 0:tr_merg%n_gamma-1), stat=fail)
         if(fail /= 0) then
           call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_init_merg')
         end if

         ! Save wcoeff.
         do i_tr = 1,n_tr
            ! Note, only one dilation in each tr (already checked).
            tr_merg%wcoeff(i_tr,:,:,:) = tr(i_tr)%wcoeff(1,:,:,:)
         end do

      end if

      ! Set as initialised.
      tr_merg%init = .true.

    end function cswt_tr_init_merg


    !--------------------------------------------------------------------------
    ! cswt_tr_free
    !
    !! Free all data associated with an initialised tr structure and reset 
    !! all other attributes.
    !!
    !! Variables:
    !!   - tr: The tr structure to be freed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_free(tr)

      type(cswt_tr), intent(inout) :: tr

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_free')
      end if 

      ! Free swav.
      if(tr%swav_mother_status) call cswt_swav_free(tr%swav_mother)

      ! Deallocate all space.
      deallocate(tr%wcoeff)
      deallocate(tr%dilation)

      ! Reset other variables.
      tr%lmax = 0
      tr%mmax = 0
      tr%n_alpha = 0
      tr%n_beta = 0
      tr%n_gamma = 0 

      tr%n_dilation = 0
      tr%wcoeff_status = .false.
      tr%swav_mother_status = .false.

      ! Set status to not initialised.
      tr%init = .false.

    end subroutine cswt_tr_free


    !--------------------------------------------------------------------------
    ! Spherical wavelet analysis
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tr_analysis
    !
    !! Perform directional continuous spherical wavelet transform to compute 
    !! wavelet coefficients.  The actual method employed to compute wavelet 
    !! coefficients depends on the method specifier.
    !!
    !! Notes:
    !!   - Data must have intent inout so can write alms to it.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet structure to compute coefficients of.  The 
    !!     computed coefficients are written to the tr%wcoeff array.
    !!   - data: The data sky to take the wavelet transform of.
    !!   - [method]: Optional method specifier to specify the method used to 
    !!     compute the spherical wavelet coefficients.
    !!   - [message_status]: Logical to indicate whether messages are to be 
    !!     displyed during the analysis routine.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_analysis(tr, data, method, message_status, time_it, &
         filename_isotropic)

      type(cswt_tr), intent(inout) :: tr
      type(s2_sky), intent(inout) :: data
      integer, intent(in),optional :: method
      logical, intent(in), optional :: message_status
      logical, intent(in), optional :: time_it
      character(len=*), intent(in), optional :: filename_isotropic

      type(s2_sky) :: sky_temp
      character(len=*), parameter :: msgpfx = 'cswt_tr_analysis> '
      integer :: count_start, count_finish, count_rate, count_max
      logical :: time = .false.
      integer :: method_use
      integer, parameter :: DEFAULT_METHOD = CSWT_TR_METHOD_FAST_FFT_REAL

      ! Set method type.
      if(present(method)) then
         method_use = method
      else
         method_use = DEFAULT_METHOD
      end if

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_analysis')
      end if 

      ! Check data initialised.
      if(.not. s2_sky_get_init(data)) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_analysis', &
          comment_add='Data not initialised')
      end if 

     ! Check wavelet present
      if(.not. tr%swav_mother_status) then
         call cswt_error(CSWT_ERROR_TR_SWAV_NPRES, 'cswt_tr_analysis')
      end if

      ! Check wavelet coefficients not already calculated.
      if(tr%wcoeff_status) then
        call cswt_error(CSWT_ERROR_TR_WCOEFF_COMP, &
          'cswt_tr_analysis')
      end if

      ! Check data sky has the same lmax and mmax (if present) as tr structure.
      if(s2_sky_get_alm_status(data)) then
         ! Check lmax.
         if(s2_sky_get_lmax(data) /= cswt_tr_get_lmax(tr)) then
            call cswt_error(CSWT_ERROR_TR_LMAX_INVALID, &
              'cswt_tr_analysis', &
              comment_add='Data has invalid lmax')
         end if
         ! Check mmax.
         if(s2_sky_get_mmax(data) /= cswt_tr_get_mmax(tr)) then
            call cswt_error(CSWT_ERROR_TR_MMAX_INVALID, &
              'cswt_tr_analysis', &
              comment_add='Data has invalid mmax')
         end if
      end if

      ! Check wavelet sky has the same lmax and mmax (if present) as 
      ! tr structure.
      sky_temp = cswt_swav_get_sky(tr%swav_mother)
      if(s2_sky_get_alm_status(sky_temp)) then
         ! Check lmax.
         if(s2_sky_get_lmax(sky_temp) /= cswt_tr_get_lmax(tr)) then
            call cswt_error(CSWT_ERROR_TR_LMAX_INVALID, &
              'cswt_tr_analysis', &
              comment_add='Spherical wavelet has invalid lmax')
         end if
         ! Check mmax.
         if(s2_sky_get_mmax(sky_temp) /= cswt_tr_get_mmax(tr)) then
            call cswt_error(CSWT_ERROR_TR_MMAX_INVALID, &
              'cswt_tr_analysis', &
              comment_add='Spherical wavelet has invalid mmax')
         end if
      end if
      call s2_sky_free(sky_temp)
  
      ! Set timer status.
      if(present(time_it)) time = time_it

      ! Begin timer if required.
      if(time) then
         call system_clock(count_start, count_rate, count_max)
      end if

      if(method_use == CSWT_TR_METHOD_DIRECT) then
         call cswt_tr_anadirect(tr, data, message=message_status)
      else   
         call cswt_tr_anafast(tr, data, method=method, &
           message=message_status, &
           filename_isotropic=filename_isotropic)
      end if

      ! Stop timer and print time taken if required.
      if(time) then
         call system_clock(count_finish)
         write(*,'(a,a,i10)') msgpfx, 'count_rate:  ', count_rate
         write(*,'(a,a,i10)') msgpfx, 'count_max:   ', count_max
         write(*,'(a,a,i10)') msgpfx, 'count_start: ', count_start
         write(*,'(a,a,i10)') msgpfx, 'count_finish:', count_finish
         write(*,'(a,a,f12.4)') msgpfx, 'time (s) = ', &
           (count_finish - count_start) / real(count_rate,s2_sp)
      end if

      tr%wcoeff_status = .true.

    end subroutine cswt_tr_analysis


    !--------------------------------------------------------------------------
    ! cswt_tr_anadirect
    !
    !! Perform the direct directional spherical wavelet transform.
    !!
    !! Variables:   
    !!   - tr: Spherical wavelet structure to compute coefficients of.  The 
    !!     computed coefficients are written to the tr%wcoeff array.
    !!   - data: The data sky to take the wavelet transform of.
    !!   - [message_status]: Logical to indicate whether messages are to be 
    !!     displyed during the analysis routine.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_anadirect(tr, data, message)

      type(cswt_tr), intent(inout) :: tr
      type(s2_sky), intent(inout) :: data
      logical, intent(in), optional :: message

      type(cswt_swav) :: swav
      type(s2_sky) :: swav_sky, swav_sky_rot
      real(s2_sp), allocatable :: swav_sky_map(:), data_map(:)
      integer :: ia, ib, ig, i_dil, fail
      real(s2_sp) :: alpha, beta, gamma, area_element

      logical :: message_status
      logical, parameter :: DEFAULT_MESSAGE_STATUS = .true.
      logical :: norm_preserve
      logical, parameter :: DEFAULT_NORM_PRESERVE = .true.
      character(len=*), parameter :: msgpfx = 'cswt_tr_anadirect> '

      ! Set message status.
      if(present(message)) then
         message_status = message
      else
         message_status = DEFAULT_MESSAGE_STATUS
      end if

      ! Set dilation norm_preserve status.
      norm_preserve = DEFAULT_NORM_PRESERVE

      ! Perform total spherical convolution for each dilation scale.
      do i_dil = 1,tr%n_dilation
      
         if(message_status) then
            write(*,'(a)') msgpfx
            write(*,'(a,a,i3,a,i3)') msgpfx, 'Considering dilation ', i_dil, &
              ' of ', tr%n_dilation
         end if
       
         ! Create copy of mother wavelet (to dilate and project onto)
         swav = cswt_swav_init(tr%swav_mother)

         ! Dilate mother wavelet.
         call cswt_swav_dilate(swav, tr%dilation(i_dil,:), norm_preserve)

         ! Extract dilated spherical wavelet sky required for spherical
         ! convolution.
         swav_sky = cswt_swav_get_sky(swav)

         ! Allocate space for maps.
         allocate(swav_sky_map(0:s2_sky_get_npix(data)-1), stat=fail)
         allocate(data_map(0:s2_sky_get_npix(data)-1), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_anadirect')
         end if
         swav_sky_map = 0.0e0
         data_map = 0.0e0

         ! Get data map (same for each loop).
         call s2_sky_get_map(data, data_map)

         ! Set area element.
         area_element = 4.0e0 * pi / real(s2_sky_get_npix(data), s2_sp)

         if(message_status) then
            write(*,'(a)') msgpfx
            write(*,'(a,a,i3,a,i3)') msgpfx, 'Considering dilation ', i_dil, &
              ' of ', tr%n_dilation
         end if

         ! Compute convolution.
         do ia = 0,tr%n_alpha-1
            
            ! Print progress.
            if(message_status) then
               write(*,'(a,a,f5.1)') msgpfx, &
                 'Percent of current dilation complete: ', &
                 ia/real(tr%n_alpha,s2_sp)*100.0e0
            end if

            alpha = 2.0e0*pi* ia / real(tr%n_alpha,s2_sp)
            do ib = 0,tr%n_beta-1
               beta = pi* ib / real(tr%n_beta,s2_sp)
               do ig = 0,tr%n_gamma-1
                  gamma = 2.0e0*pi* ig / real(tr%n_gamma,s2_sp)
                 
                  ! Rotate wavelet for current orientation.
                  swav_sky_rot = s2_sky_init(swav_sky)
                  call s2_sky_rotate(swav_sky_rot, alpha, beta, gamma)
                  
                  ! Get map for rotated wavelet.
                  call s2_sky_get_map(swav_sky_rot, swav_sky_map)

                  ! Perfom integration
                  ! Note conjugate disappears since real signals.
                  tr%wcoeff(i_dil,ia,ib,ig) = sum(swav_sky_map * data_map)

                  ! Free rotated wavelet used for current orientation.
                  call s2_sky_free(swav_sky_rot)

                  ! Reset current rotated wavelet map values.
                  ! Extracted again in next loop so no real need to reset
                  ! here but done to flag any bugs.
                  swav_sky_map = 0.0e0

               end do
            end do
         end do

      end do

      ! Scale by area element.
      tr%wcoeff = tr%wcoeff * area_element
      
      ! Free memory.
      call cswt_swav_free(swav)
      call s2_sky_free(swav_sky)
      deallocate(swav_sky_map, data_map)     

    end subroutine cswt_tr_anadirect


    !--------------------------------------------------------------------------
    ! cswt_tr_anafast
    !
    !! Perform the fast directional spherical wavelet transform.
    !!
    !! Notes:
    !!   - For case where alms are already defined (and map isn't) then set 
    !!     dilation to unitary; no dilation will actually be performed and
    !!     analysis can continue.
    !!
    !! Variables:   
    !!   - tr: Spherical wavelet structure to compute coefficients of.  The 
    !!     computed coefficients are written to the tr%wcoeff array.
    !!   - data: The data sky to take the wavelet transform of.
    !!   - [method]: Optional method specifier to specify the method used to 
    !!     compute the spherical wavelet coefficients.
    !!   - [message_status]: Logical to indicate whether messages are to be 
    !!     displyed during the analysis routine.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_anafast(tr, data, method, message, &
         filename_isotropic)

      type(cswt_tr), intent(inout) :: tr
      type(s2_sky), intent(inout) :: data
      integer, intent(in), optional :: method
      logical, intent(in), optional :: message
      character(len=*), intent(in), optional :: filename_isotropic

      integer :: method_use
      logical :: message_status
      logical, parameter :: DEFAULT_MESSAGE_STATUS = .true.
      logical :: norm_preserve
      logical, parameter :: DEFAULT_NORM_PRESERVE = .true.
      character(len=*), parameter :: msgpfx = 'cswt_tr_anafast> '
      real(s2_sp), parameter :: AZBL_CUTOFF_PROP = 0.95
      real(s2_sp) :: azimuthal_bl

      integer, parameter :: DEFAULT_METHOD = CSWT_TR_METHOD_FAST_FFT_REAL
      integer :: i_dil
      type(cswt_swav) :: swav
      type(s2_sky) :: swav_sky
      character(len=S2_STRING_LEN) :: error_comment

      ! Set message status.
      if(present(message)) then
         message_status = message
      else
         message_status = DEFAULT_MESSAGE_STATUS
      end if

      ! Set method type.
      if(present(method)) then
         method_use = method
      else
         method_use = DEFAULT_METHOD
      end if

      ! Set dilation norm_preserve status.
      norm_preserve = DEFAULT_NORM_PRESERVE

      ! Compute data alm.  (Only need to do once.)
      if(message_status) then
         write(*,'(a,a)') msgpfx, 'Starting fast cswt analysis'
         write(*,'(a,a)') msgpfx, 'Computing data alms'
      end if
      call s2_sky_compute_alm(data, tr%lmax, tr%mmax, message_status)

      ! Perform total spherical convolution for each dilation scale.
      do i_dil = 1,tr%n_dilation
      
         if(message_status) then
            write(*,'(a)') msgpfx
            write(*,'(a,a,i3,a,i3)') msgpfx, 'Considering dilation ', i_dil, &
              ' of ', tr%n_dilation
         end if
       
         ! Create copy of mother wavelet (to dilate and project onto)
         swav = cswt_swav_init(tr%swav_mother)

         ! Dilate mother wavelet.
         call cswt_swav_dilate(swav, tr%dilation(i_dil,:), norm_preserve)

         ! Compute alm for dilated wavelet.
         if(message_status) then
            write(*,'(a,a)') msgpfx, 'Computing dilated spherical wavelet alms'
         end if
         call cswt_swav_compute_alm(swav, tr%lmax, tr%mmax, message_status)

         ! Check azimuthal band width sufficient.
         azimuthal_bl = cswt_swav_azimuthal_bl(swav, AZBL_CUTOFF_PROP)
         if(tr%n_gamma < 2*azimuthal_bl) then
            write(error_comment,'(a,f4.0)') 'Azimuthal bandlimit = ', &
              azimuthal_bl 
            call cswt_error(CSWT_ERROR_TR_AZBANDLIM_INVALID, &
              'cswt_tr_anafast', comment_add=trim(error_comment))
         end if

         ! Compute spherical convolution of dilated wavelet and data.

         if(message_status) then
            write(*,'(a,a)') msgpfx, 'Computing spherical convolution'
         end if

         ! Extract dilated spherical wavelet sky required for spherical
         ! convolution.
         swav_sky = cswt_swav_get_sky(swav)

         ! Perform convolution.
         if(method_use == CSWT_TR_METHOD_FAST_ISOTROPIC) then
            call cswt_tr_anafast_isotropic(tr, i_dil, data, swav_sky, &
                 message=message, filename=filename_isotropic)
         else
            call cswt_tr_anafast_spcv(tr, i_dil, data, swav_sky, &
              method_in=method, message=message)
         end if

         ! Free temporary variables.
         call s2_sky_free(swav_sky)
         call cswt_swav_free(swav)

      end do

      if(message_status) then
         write(*,'(a)') msgpfx
         write(*,'(a,a)') msgpfx, 'Fast cswt analysis complete!'
      end if

    end subroutine cswt_tr_anafast


    !--------------------------------------------------------------------------
    ! cswt_tr_anafast_isotropic
    !
    !! Perform the fast isotropic spherical convolution for a particular 
    !! dilation (noting the the convolution is simply given by a product of 
    !! spherical harmonic coefficients).
    !!
    !! Notes:
    !!   - Sizes are checked before here so no checking of size consistency is 
    !!     performed herein.  (All skies must have same lmax and mmax as values
    !!     passed to function.)
    !!   - tr%wcoeff array must be already allocated and must be indexed 
    !!     from 0:N-1 for Euler angle dimensions.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet structure to compute coefficients of.  The 
    !!     computed coefficients are written to the tr%wcoeff array.
    !!   - i_dil: Dilation index to store computed coefficient in.
    !!   - sky: Data to convolve.
    !!   - swav: Dilated wavelet to convolve.
    !!   - [message_status]: Logical to indicate whether messages are to be 
    !!     displyed during the analysis routine.
    !!   - [filename]: Name of fits file to save healpix pixelisation of 
    !!     wavelet coefficients in (if present, else not saved)
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   May 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_anafast_isotropic(tr, i_dil, sky, swav, message, &
         filename)

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: i_dil
      type(s2_sky), intent(in) :: sky, swav
      logical, intent(in), optional :: message
      character(len=*), intent(in), optional :: filename

      logical :: message_status
      logical, parameter :: DEFAULT_MESSAGE_STATUS = .true.
      character(len=*), parameter :: msgpfx = 'cswt_tr_anafast_isotropic> '

      integer :: lmax, mmax, l, m, fail
      complex(s2_spc), allocatable :: sky_alm(:,:)
      complex(s2_spc), allocatable :: swav_alm(:,:)
      complex(s2_spc), allocatable :: wcoeff_alm(:,:)
      real(s2_sp), allocatable :: wcoeff_ab(:,:)
      type(s2_sky) :: wcoeff_sky
      character(len=S2_STRING_LEN) :: filename_use

      ! Check only considering one orientation since isotropic.
      if(tr%n_gamma /= 1) then
         call cswt_error(CSWT_ERROR_TR_NGAMMA_INVALID, &
           'cswt_tr_anafast_isotropic', &
           comment_add='Too many orientations for isotropic cswt algorithm')
      end if

      ! Set message status.
      if(present(message)) then
         message_status = message
      else
         message_status = DEFAULT_MESSAGE_STATUS
      end if

      ! Print method type message.
      if(message_status) then
         write(*, '(a,a,i5)') msgpfx, &
              'Method: ** Fast isotropic **'
      end if

      ! Set sizes.
      lmax = tr%lmax
      mmax = tr%mmax

      ! Get copies of sky alms for manipulation.
      allocate(sky_alm(0:lmax, 0:mmax), stat=fail)
      allocate(swav_alm(0:lmax, 0:mmax), stat=fail)
      allocate(wcoeff_alm(0:lmax, 0:mmax), stat=fail)
      allocate(wcoeff_ab(0:tr%n_alpha-1, 0:tr%n_beta-1), stat=fail)
      sky_alm = 0.0e0
      swav_alm = 0.0e0
      wcoeff_alm = 0.0e0
      wcoeff_ab = 0.0e0
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_anafast_isotropic')
      end if
      call s2_sky_get_alm(sky, sky_alm)
      call s2_sky_get_alm(swav, swav_alm)

      ! Multiply harmonic coefficients.    
      ! Note currently has different normalisation to full 3 rotation
      ! spherical convolution.
      do l = 0,lmax 
         do m = 0,min(l,mmax)
            wcoeff_alm(l,m) = sqrt(4*pi/real(2*l+1,s2_sp)) &
                 * swav_alm(l,0) * sky_alm(l,m)
         end do
      end do

      ! Initialise wavelet coefficient sky object from alms.
      wcoeff_sky = s2_sky_init(wcoeff_alm, lmax, mmax, &
        cswt_swav_get_nside(tr%swav_mother), &
        cswt_swav_get_pix_scheme(tr%swav_mother))

      ! Compute wavelet sky map.
      call s2_sky_compute_map(wcoeff_sky, message=message_status)

      ! Save healpix map of wavelet coefficients if required.
      if(present(filename)) then
         filename_use = filename
         call s2_sky_write_map_file(wcoeff_sky, filename_use)
      end if

      ! Convert healpix map to ECP euler angle pixelisation of tr structure.
      call s2_sky_extract_ab(wcoeff_sky, wcoeff_ab)

      ! Copy coefficients back to tr structure.
      tr%wcoeff(i_dil,:,:,0) = wcoeff_ab

      ! Free memory.
      deallocate(sky_alm, swav_alm, wcoeff_alm)
      deallocate(wcoeff_ab)
      call s2_sky_free(wcoeff_sky)

    end subroutine cswt_tr_anafast_isotropic


    !--------------------------------------------------------------------------
    ! cswt_tr_anafast_spcv
    !
    !! Perform the fast spherical convolution for a particular dilation 
    !! using the fast algorithm of Wandlet and Gorski.
    !!
    !! Notes:
    !!   - Sizes are checked before here so no checking of size consistency is 
    !!     performed herein.  (All skies must have same lmax and mmax as values
    !!     passed to function.)
    !!   - tr%wcoeff array must be already allocated and must be indexed 
    !!     from 0:N-1 for Euler angle dimensions.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet structure to compute coefficients of.  The 
    !!     computed coefficients are written to the tr%wcoeff array.
    !!   - i_dil: Dilation index to store computed coefficient in.
    !!   - sky1: Data to convolve.
    !!   - sky2: Dilated wavelet to convolve.
    !!   - [method_in]: Optional method specifier to specify the method used 
    !!     to compute the spherical wavelet coefficients.
    !!   - [message_status]: Logical to indicate whether messages are to be 
    !!     displyed during the analysis routine.
    !
    !! @author D. J. Mortlock and J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   Originally written by Daniel Mortlock
    !   November 2004 - Modified by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_anafast_spcv(tr, i_dil, sky1, sky2, method_in, message)
      
      use s2_dl_mod, only: s2_dl_beta_operator

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: i_dil
      type(s2_sky), intent(in) :: sky1, sky2
      integer, intent(in), optional :: method_in
      logical, intent(in), optional :: message

      real(s2_dp), parameter :: &
        PIONTWO = 1.570796326794896619231321691639751442099, &
        TWOPI = 6.283185307179586476925286766559005768394

      logical :: message_status
      logical, parameter :: DEFAULT_MESSAGE_STATUS = .true.
      integer :: method
      integer, parameter :: DEFAULT_METHOD = CSWT_TR_METHOD_FAST_FFT_REAL
      character(len=*), parameter :: msgpfx = 'cswt_tr_anafast_spcv> '
      real(s2_dp), parameter :: ALPHA_MAX = TWOPI, BETA_MAX = TWOPI, &
           & GAMMA_MAX = TWOPI

      integer :: n_a, n_b, n_b_out, n_g, n_gon2
      integer :: l, m, mm, mmm, mmmmax_l, a, b, g
      integer :: lmax, mmax
      integer :: fail
      real(s2_dp) :: rm, msign, mmmsign, alpha, beta, gamma, angle
      real(s2_dp), pointer :: dl(:, :) => null() 
      complex(s2_dpc), allocatable :: exppion2(:)
      complex(s2_dpc), allocatable :: tmmm(:,:,:)
      complex(s2_dpc), allocatable :: tmmm_shift(:,:,:)
      complex(s2_dpc) :: sb, expangle
      complex(s2_spc), allocatable :: sky1_alm(:,:), sky2_alm(:,:)
      
      ! FFTW variables
      integer, parameter :: FFTW_ESTIMATE=0,FFTW_MEASURE=1
      integer, parameter :: FFTW_FORWARD=-1,FFTW_BACKWARD=1
      integer, parameter ::  FFTW_REAL_TO_COMPLEX=-1,FFTW_COMPLEX_TO_REAL=1
      integer, parameter ::  FFTW_OUT_OF_PLACE=0,FFTW_IN_PLACE=8
      integer*8 :: plan
      real(s2_dp), allocatable :: fftw_out(:,:,:)

   
      ! Set message status.
      if(present(message)) then
         message_status = message
      else
         message_status = DEFAULT_MESSAGE_STATUS
      end if

      ! Set method type.
      if(present(method_in)) then
         method = method_in
      else
         method = DEFAULT_METHOD
      end if

  
      !------------------------------------------------------------------------

      ! --------------------------------------
      ! CSWT_TR_METHOD_FAST_FFT_REAL
      ! --------------------------------------
      if(method == CSWT_TR_METHOD_FAST_FFT_REAL) then

         ! Print method type message.
         if(message_status) then
           write(*, '(a,a,i5)') msgpfx, &
             'Method: ** Fast real **'
         end if


         ! --------------------------------------
         ! Initialise sizes and data arrays
         ! --------------------------------------

         ! Get copies of sky alms for manipulation.
         allocate(sky1_alm(0:s2_sky_get_lmax(sky1), &
           0:s2_sky_get_mmax(sky1)), stat=fail)
         allocate(sky2_alm(0:s2_sky_get_lmax(sky2), &
           0:s2_sky_get_mmax(sky2)), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_anafast_spcv')
         end if
         call s2_sky_get_alm(sky1, sky1_alm)
         call s2_sky_get_alm(sky2, sky2_alm)

         ! Set sizes.
         lmax = tr%lmax
         mmax = tr%mmax
         n_a = 2 * lmax + 1
         n_b = 2 * lmax + 1      ! Beta ranges over full 2*pi so can use fft.
         n_b_out = lmax + 1      ! Only output beta values up to pi.
         if(mod(tr%n_gamma,2) == 0) then
            ! Even case.
            ! Should already be checked when initialise tr structure 
            ! but check again here.
            call cswt_error(CSWT_ERROR_TR_NGAMMA_INVALID, &
                 'cswt_tr_anafast_spcv')
            ! Should stop here but if continue set n_g odd.
            n_g = tr%n_gamma + 1
         else
            ! Odd case.
            n_g = tr%n_gamma
         end if
         n_gon2 = (n_g - 1) / 2

         ! Print beginning message.
         if(message_status) then
           write(*, '(a,a,i5)') msgpfx, &
             'Creating total convolution to l_max = ', lmax
         end if
         
         ! Allocate space for dl.
         allocate(dl(-lmax:lmax, -lmax:lmax))
         
         ! Precalculate the complex exponentials.
         allocate(exppion2(-lmax:lmax))
         do m = - lmax, lmax
            rm = real(m)
            exppion2(m) = cmplx(cos(rm * PIONTWO), sin(rm * PIONTWO))
         end do
         
         ! Allocate T_mmm, the triple modified Fourier transform 
         ! of T_abg, first setting all its elements to zero. 
         ! Note only store positive m indices since negative m vales
         ! implicitly given by conjugate symmetry relationship.
         allocate(tmmm(0:lmax, 0:2*lmax, 0:2*n_gon2))
         tmmm = 0.0d0
         
         
         ! --------------------------------------
         ! Calculate T_mmm
         ! --------------------------------------
         
         ! Create T_mmm, looping over l, and adding to the elements
         ! as required.
         
         if(message_status) then
            write(*, '(a,a)') msgpfx, 'Calculating T_mmm'
         end if

         do l = 0, lmax
            
            ! For each l value create the plane of the d-matrix.
            call s2_dl_beta_operator(dl, PIONTWO, l)
            
            ! For each l value loop over all the elements of the 3-dimensional
            ! array -- m, mm, and mmm.

            msign = -1.0e0

            do m = 0, l
               msign = - msign

!** TODO: could incorporate additional symmetry to restrict mm>=0 only.               
               do mm = -l, l
                  
                  mmmmax_l = min(n_gon2, l)
                  if (mod(mmmmax_l,2) == 0) then
                     ! Even case.
                     mmmsign = - 1.0
                  else
                     ! Odd case.
                     mmmsign = 1.0
                  end if
                  do mmm = -mmmmax_l, mmmmax_l
                     mmmsign = - mmmsign
                     
                     ! In the product of the beam and sky, it is the complex
                     ! conjugate of the beam multiplied with the sky. Note
                     ! that it is assumed that b_P_lm is set to zero for the
                     ! l < 2 modes that don't exist for the linear polarization
                     ! modes.

                     if ((m < 0) .and. (mmm < 0)) then
                        sb = msign * mmmsign &
                             * sky2_alm(l, - mmm) * conjg(sky1_alm(l, - m))
                     else if (m < 0) then
                        sb = msign &
                             * conjg(sky2_alm(l, mmm) * sky1_alm(l, - m))
                     else if (mmm < 0) then
                        sb = mmmsign &
                             * sky2_alm(l, - mmm) * sky1_alm(l, m)
                     else
                        sb = conjg(sky2_alm(l, mmm)) * sky1_alm(l, m)
                     end if

                     tmmm(m, mm+lmax, mmm+n_gon2) &
                       = tmmm(m, mm+lmax, mmm+n_gon2) &
                       + exppion2(- m)*exppion2(mmm)*dl(mm, m)*dl(mm, mmm) &
                          * sb

                  end do
               end do
            end do
            
         end do

         ! Shift tmmm spectrum so origin at center.

         allocate(tmmm_shift(0:lmax, 0:2*lmax, 0:2*n_gon2))

         do mmm=0,2*n_gon2
            tmmm_shift(:,0:lmax,mmm) = tmmm(:,lmax:2*lmax,mmm)
            tmmm_shift(:,lmax+1:2*lmax,mmm) = tmmm(:,0:lmax-1,mmm)
         end do
         tmmm = tmmm_shift

         if(n_gon2 > 0) then
            tmmm_shift(:,:,0:n_gon2) = tmmm(:,:,n_gon2:2*n_gon2)
            tmmm_shift(:,:,n_gon2+1:2*n_gon2) = tmmm(:,:,0:n_gon2-1)
            tmmm = tmmm_shift
         end if

         deallocate(tmmm_shift)


         ! --------------------------------------
         ! IFFT 
         ! --------------------------------------
    
         if(message_status) then         
            write(*, '(a,a)') msgpfx, &
                 'Fast Fourier transforming from T_mmm to T_abg'
         end if

         ! Allocate space for fftw output 
         allocate(fftw_out(0:2*lmax,0:2*lmax,0:2*n_gon2))
         fftw_out = 0.0d0

         ! Take inverse fft of tmmm (assuming conjugate symmetry 
         ! relationship for values corresponding to negative m index).
         call rfftw3d_f77_create_plan(plan, 2*lmax+1, 2*lmax+1, &
              2*n_gon2+1, FFTW_COMPLEX_TO_REAL, FFTW_ESTIMATE)
         call rfftwnd_f77_one_complex_to_real(plan, &
              tmmm(0:lmax,0:2*lmax,0:2*n_gon2), &
              fftw_out(0:2*lmax,0:2*lmax,0:2*n_gon2))
         call rfftwnd_f77_destroy_plan(plan)
         
         ! Copy wavelet coefficients back (casting as single precision).
         tr%wcoeff(i_dil,0:2*lmax,0:lmax,0:2*n_gon2) &
           = real(fftw_out(0:2*lmax,0:lmax,0:2*n_gon2), kind=s2_sp) 

        
         ! --------------------------------------
         ! Free memory
         ! --------------------------------------

          deallocate(sky1_alm)
          deallocate(sky2_alm)
          deallocate(tmmm)
          deallocate(fftw_out)
          deallocate(exppion2)
          deallocate(dl)


      !------------------------------------------------------------------------

      ! --------------------------------------
      ! CSWT_TR_METHOD_FAST_FFT
      ! --------------------------------------
      elseif(method == CSWT_TR_METHOD_FAST_FFT) then

        ! Print method type message.
         if(message_status) then
           write(*, '(a,a,i5)') msgpfx, &
             'Method: ** Fast (not exploting conjugate symmetry) **'
         end if


         ! --------------------------------------
         ! Initialise sizes and data arrays
         ! --------------------------------------

         ! Get copies of sky alms for manipulation.
         allocate(sky1_alm(0:s2_sky_get_lmax(sky1), &
           0:s2_sky_get_mmax(sky1)), stat=fail)
         allocate(sky2_alm(0:s2_sky_get_lmax(sky2), &
           0:s2_sky_get_mmax(sky2)), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_anafast_spcv')
         end if
         call s2_sky_get_alm(sky1, sky1_alm)
         call s2_sky_get_alm(sky2, sky2_alm)

         ! Set sizes.
         lmax = tr%lmax
         mmax = tr%mmax
         n_a = 2 * lmax + 1
         n_b = 2 * lmax + 1      ! Beta ranges over full 2*pi so can use fft.
         n_b_out = lmax + 1      ! Only output beta values up to pi.
         if(mod(tr%n_gamma,2) == 0) then
            ! Even case.
            ! Should already be checked when initialise tr structure 
            ! but check again here.
            call cswt_error(CSWT_ERROR_TR_NGAMMA_INVALID, &
                 'cswt_tr_anafast_spcv')
            ! Should stop here but if continue set n_g odd.
            n_g = tr%n_gamma + 1
         else
            ! Odd case.
            n_g = tr%n_gamma
         end if
         n_gon2 = (n_g - 1) / 2

         ! Print beginning message.
         if(message_status) then
           write(*, '(a,a,i5)') msgpfx, &
             'Creating total convolution to l_max = ', lmax
         end if
         
         ! Allocate space for dl.
         allocate(dl(-lmax:lmax, -lmax:lmax))
         
         ! Precalculate the complex exponentials.
         allocate(exppion2(-lmax:lmax))
         do m = - lmax, lmax
            rm = real(m)
            exppion2(m) = cmplx(cos(rm * PIONTWO), sin(rm * PIONTWO))
         end do
         
         ! Allocate T_mmm, the triple modified Fourier transform 
         ! of T_abg, first setting all its elements to zero.
         ! (Note that tmmm allocated differently to slow implementation 
         ! above. Indexed from 0 here, whereas indexed from -lmax above.)
         
         allocate(tmmm(0: 2 * lmax, 0: 2 * lmax, 0: 2 * n_gon2))
         tmmm = 0.0d0
         
         
         ! --------------------------------------
         ! Calculate T_mmm
         ! --------------------------------------
         
         ! Create T_mmm, looping over l, and adding to the elements
         ! as required.
         
         if(message_status) then
            write(*, '(a,a)') msgpfx, 'Calculating T_mmm'
         end if

         do l = 0, lmax
            
            ! For each l value create the plane of the d-matrix.
            call s2_dl_beta_operator(dl, PIONTWO, l)
            
            ! For each l value loop over all the elements of the 3-dimensional
            ! array -- m, mm, and mmm.
            
            if (mod(l,2) == 0) then
               ! Even case.
               msign = - 1.0
            else
               ! Odd case.
               msign = 1.0
            end if
            do m = - l, l
               msign = - msign
               
               do mm = - l, l
                  
                  mmmmax_l = min(n_gon2, l)
                  if (mod(mmmmax_l,2) == 0) then
                     ! Even case.
                     mmmsign = - 1.0
                  else
                     ! Odd case.
                     mmmsign = 1.0
                  end if
                  do mmm = - mmmmax_l, mmmmax_l
                     mmmsign = - mmmsign
                     
                     ! In the product of the beam and sky, it is the complex
                     ! conjugate of the beam multiplied with the sky. Note
                     ! that it is assumed that b_P_lm is set to zero for the
                     ! l < 2 modes that don't exist for the linear polarization
                     ! modes.

                     if ((m < 0) .and. (mmm < 0)) then
                        sb = msign * mmmsign &
                             * sky2_alm(l, - mmm) * conjg(sky1_alm(l, - m))
                     else if (m < 0) then
                        sb = msign &
                             * conjg(sky2_alm(l, mmm) * sky1_alm(l, - m))
                     else if (mmm < 0) then
                        sb = mmmsign &
                             * sky2_alm(l, - mmm) * sky1_alm(l, m)
                     else
                        sb = conjg(sky2_alm(l, mmm)) * sky1_alm(l, m)
                     end if
                     
                     tmmm(m + lmax, mm + lmax, mmm + n_gon2) &
                       = tmmm(m + lmax, mm + lmax, mmm + n_gon2) &
                       + exppion2(- m)*exppion2(mmm)*dl(mm, m)*dl(mm, mmm) &
                          * sb

                  end do
               end do
            end do
            
         end do


         ! --------------------------------------
         ! IFFT and modulate
         ! --------------------------------------
    
         if(message_status) then         
            write(*, '(a,a)') msgpfx, &
                 'Fast Fourier transforming from T_mmm to T_abg'
         end if
         
         ! Perform 3d inverse FFT.
         call cswt_tr_anafast_fft3d(tmmm, backward=.true.)
          
         ! Modulate inverse fft calculated to construct the modified 
         ! transform required and write result onto total convolution 
         ! data structure for output beta range only.
         
         do a = 0, n_a - 1
            alpha = ALPHA_MAX * real(a) / real(n_a)
            do b = 0, n_b_out - 1
               beta = BETA_MAX * real(b) / real(n_b)  
               do g = 0, n_g - 1
                  gamma = GAMMA_MAX * real(g) / real(n_g)
                  angle = -( (alpha+beta)*lmax + gamma*n_gon2 )
                  expangle = cmplx(cos(angle), sin(angle))

                  ! Modulated to account for negative indices of required
                  ! transform.
                  tr%wcoeff(i_dil, a, b, g) &
                       = real(expangle * tmmm(a,b,g), kind=s2_sp)  

              end do
            end do
          end do


         ! --------------------------------------
         ! Free memory
         ! --------------------------------------

          deallocate(sky1_alm)
          deallocate(sky2_alm)
          deallocate(tmmm)
          deallocate(exppion2)
          deallocate(dl)


      !------------------------------------------------------------------------

      ! --------------------------------------
      ! CSWT_TR_METHOD_FAST_DFT
      ! --------------------------------------
       else if(method == CSWT_TR_METHOD_FAST_DFT) then

         ! Essentially identical to FFT technique except instead of using fft
         ! followed by modulating the result, a modified inverse DFT is 
         ! applied.  Only present to verify FFT technique.  In general use
         ! FFT technique.

         ! Print method type message.
         if(message_status) then
            write(*, '(a,a,i5)') msgpfx, &
                 'Method: ** Fast DFT (not exploting conjugate symmetry) **'
         end if

         ! --------------------------------------
         ! Initialise sizes and data arrays
         ! --------------------------------------

         ! Get copies of sky alms for manipulation.
         allocate(sky1_alm(0:s2_sky_get_lmax(sky1), &
           0:s2_sky_get_mmax(sky1)), stat=fail)
         allocate(sky2_alm(0:s2_sky_get_lmax(sky2), &
           0:s2_sky_get_mmax(sky2)), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_anafast_spcv')
         end if
         call s2_sky_get_alm(sky1, sky1_alm)
         call s2_sky_get_alm(sky2, sky2_alm)

         ! Set sizes.
         lmax = tr%lmax
         mmax = tr%mmax
         n_a = 2 * lmax + 1
         n_b = 2 * lmax + 1      ! Beta ranges over full 2*pi so can use fft.
         n_b_out = lmax + 1      ! Only output beta values up to pi.
         if(mod(tr%n_gamma,2) == 0) then
            ! Even case.
            ! Should already be checked when initialise tr structure 
            ! but check again here.
            call cswt_error(CSWT_ERROR_TR_NGAMMA_INVALID, &
                 'cswt_tr_anafast_spcv')
            ! Should stop here but if continue set n_g odd.
            n_g = tr%n_gamma + 1
         else
            ! Odd case.
            n_g = tr%n_gamma
         end if
         n_gon2 = (n_g - 1) / 2
 
         ! Print beginning message.
         if(message_status) then
           write(*, '(a,a,i5)') msgpfx, &
             'Creating total convolution to l_max = ', lmax
         end if

         ! Allocate space for dl.
         allocate(dl(-lmax:lmax, -lmax:lmax))
         
         ! Precalculate the complex exponentials.
         allocate(exppion2(-lmax:lmax))
         do m = - lmax, lmax
            rm = real(m)
            exppion2(m) = cmplx(cos(rm * PIONTWO), sin(rm * PIONTWO))
         end do
         
         ! Allocate T_mmm, the triple modified Fourier transform 
         ! of T_abg, first setting all its elements to zero.
         ! (Note that tmmm allocated differently to slow implementation 
         ! above. Indexed from 0 here, whereas indexed from -lmax above.)
         
         allocate(tmmm(0: 2 * lmax, 0: 2 * lmax, 0: 2 * n_gon2))
         tmmm = 0.0d0
         
         
         ! --------------------------------------
         ! Calculate T_mmm
         ! --------------------------------------
         
         ! Create T_mmm, looping over l, and adding to the elements
         ! as required.
         
         if(message_status) then
            write(*, '(a,a)') msgpfx, 'Calculating T_mmm'
         end if

         do l = 0, lmax
            
            ! For each l value create the plane of the d-matrix.
            call s2_dl_beta_operator(dl, PIONTWO, l)
            
            ! For each l value loop over all the elements of the 3-dimensional
            ! array -- m, mm, and mmm.
            
            if (mod(l,2) == 0) then
               ! Even case.
               msign = - 1.0
            else
               ! Odd case.
               msign = 1.0
            end if
            do m = - l, l
               msign = - msign
               
               do mm = - l, l
                  
                  mmmmax_l = min(n_gon2, l)
                  if (mod(mmmmax_l,2) == 0) then
                     ! Even case.
                     mmmsign = - 1.0
                  else
                     ! Odd case.
                     mmmsign = 1.0
                  end if
                  do mmm = - mmmmax_l, mmmmax_l
                     mmmsign = - mmmsign
                     
                     ! In the product of the beam and sky, it is the complex
                     ! conjugate of the beam multiplied with the sky. Note
                     ! that it is assumed that b_P_lm is set to zero for the
                     ! l < 2 modes that don't exist for the linear polarization
                     ! modes.

                     if ((m < 0) .and. (mmm < 0)) then
                        sb = msign * mmmsign &
                             * sky2_alm(l, - mmm) * conjg(sky1_alm(l, - m))
                     else if (m < 0) then
                        sb = msign &
                             * conjg(sky2_alm(l, mmm) * sky1_alm(l, - m))
                     else if (mmm < 0) then
                        sb = mmmsign &
                             * sky2_alm(l, - mmm) * sky1_alm(l, m)
                     else
                        sb = conjg(sky2_alm(l, mmm)) * sky1_alm(l, m)
                     end if
                     
                     tmmm(m + lmax, mm + lmax, mmm + n_gon2) &
                       = tmmm(m + lmax, mm + lmax, mmm + n_gon2) &
                       + exppion2(- m)*exppion2(mmm)*dl(mm, m)*dl(mm, mmm) &
                          * sb

                  end do
               end do
            end do
            
         end do


         ! --------------------------------------
         ! IDFT and modulate
         ! --------------------------------------
    
         if(message_status) then         
            write(*, '(a,a,a)') msgpfx, &
              'Direct discrete Fourier transforming (*not fast*)', &
              ' from T_mmm to T_abg'
         end if

          do a = 0, n_a - 1
            alpha = ALPHA_MAX * real(a) / real(n_a)
            do b = 0, n_b_out - 1
              beta = BETA_MAX * real(b) / real(n_b)
              do g = 0, n_g - 1   
                gamma = GAMMA_MAX * real(g) / real(n_g)
                tr%wcoeff(i_dil,a,b,g) &
                  = cswt_tr_anafast_imdft(tmmm, alpha, beta, gamma)
               end do
            end do
          end do

         ! --------------------------------------
         ! Free memory
         ! --------------------------------------

          deallocate(sky1_alm)
          deallocate(sky2_alm)
          deallocate(tmmm)
          deallocate(exppion2)
          deallocate(dl)


      ! --------------------------------------
      ! Invalid method
      ! --------------------------------------
       else

          call cswt_error(CSWT_ERROR_TR_METHOD_INVALID, &
               'cswt_tr_anafast_spcv')

       end if

    end subroutine cswt_tr_anafast_spcv


    !--------------------------------------------------------------------------
    ! cswt_tr_anafast_fft3d
    !
    !! Perform 3D FFT using healpix fft routines (fftw if linked).
    !!   
    !! Notes:
    !!   - Since the Fourier transform is separable the 3D transform is
    !!     calculated by taking suitable 1D transforms over each dimension
    !!     of the data.
    !!   - No scaling performed - even on reverse transform.
    !!
    !! Variables:
    !!   - x: The 3D input data on input.  The 3D FFT of the inputted data
    !!     on exit.
    !!   - backward: Logical specifying whether the forward FFT (.false.) or
    !!     reverse inverse FFT (.true.) is performed.
    !    
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   February 2004 - Written by Jason McEwen
    !   November 2004 - Incorporated into cswt library by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_anafast_fft3d(x, backward)

      use healpix_fft, only: complex_fft

      complex(s2_dpc), intent(inout) :: x(:,:,:)
      logical, intent(in) :: backward
      
      integer :: i,j,k

      do i = 1,size(x,1)
        do j = 1,size(x,2)
          call complex_fft(x(i,j,:), backward)
        end do 
      end do

      do i = 1,size(x,1)
        do k = 1,size(x,3)
          call complex_fft(x(i,:,k), backward)
        end do 
      end do

      do j = 1,size(x,2)
        do k = 1,size(x,3)
          call complex_fft(x(:,j,k), backward)
        end do 
      end do

    end subroutine cswt_tr_anafast_fft3d


    !--------------------------------------------------------------------------
    ! cswt_tr_anafast_imdft
    !
    !! Perform the modified inverse Fourier transform required to convert
    !! T_mmm to T_abg directly.
    !!
    !! Notes:
    !!   - t may be computed for any continuous values of alpha, beta and
    !!     gamma.
    !!
    !! Variables:
    !!   - tmmm: Input data
    !!     representing the modified 3d Fourier transform of required
    !!     t_abg.  (Now allocated array rather than pointer so indexed from
    !!     1 when accessed herein.)
    !!   - alpha: Alpha Euler angle to calculate t for.
    !!   - beta: Beta Euler angle to calculate t for.
    !!   - gamma: Gamma Euler angle to calculate t for.
    !!   - t: The modified inverse Fourier transform of tmmm for the
    !!     specified angles.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   February 2004 - Written by Jason McEwen
    !   November 2004 - Incorporated into cswt library by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_anafast_imdft(tmmm, alpha, beta, gamma) result(t)

      complex(s2_dpc), intent(in) :: tmmm(:, :, :)
      real(s2_dp), intent(in) :: alpha, beta, gamma
      real(s2_sp)  :: t

      integer :: m, mm, mmm, lmax, n_g, n_gon2
      real(s2_dp) :: angle
      complex(s2_dpc) :: expangle

      lmax = (size(tmmm,1) - 1) / 2

      n_g = size(tmmm,3)
      n_gon2 = (n_g - 1) / 2

      t = 0.0
      do m = - lmax, lmax
         do mm = - lmax, lmax
            do mmm = - n_gon2, n_gon2
               angle = real(m) * alpha + real(mm) * beta + real(mmm) * gamma
               expangle = cmplx(cos(angle), sin(angle))
               t = t &
                 + real(expangle * tmmm(m+lmax+1, mm+lmax+1, mmm+n_gon2+1), &
                   kind = s2_sp)
            end do
         end do
      end do

    end function cswt_tr_anafast_imdft


    !--------------------------------------------------------------------------
    ! Wavelet coefficient operation routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tr_norm
    !
    !! Compute the norm of the wavelet coefficients stored in the tr strucutre 
    !! for each dilation.  The norm is defined as the sum, over all Euler
    !! angles, of the squared coefficient values, i.e. 
    !! norm(a) = 1/N_s * \sum_i ( k_i (W_i)^2 )
    !! where W_i are the wavelet coefficients, k_i the weights (if
    !! applied) and N_s is the number of samples per dilation, 
    !! i.e. N_s = N_alpha * N_beta * N_gamma.  The weights may be optionally
    !! applied.  *Note that N_s=1 now*
    !!  
    !! Notes:
    !!   - Take note that the norm computed is scaled by one over the number
    !!     of samples per dilation.  -- *Not any longer*: must do scaling 
    !!     separately, since in some cases want to scale by number of
    !!     *effective* coefficients, not total number.
    !!
    !! Variables:
    !!   - tr: The tr structure containing the wavelet coefficients to
    !!     compute the dilation of.
    !!   - norm: Computed norm values for each dilation.
    !!   - [weight_in]: Logical to specified whether pixel weights are to be 
    !!     applied.  If not specified then WEIGHT_DEFAULT is assumed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !   August 2005 -   Jason McEwen
    !                   No longer perform any scaling.
    !--------------------------------------------------------------------------

    subroutine cswt_tr_norm(tr, norm, weight_in)

      type(cswt_tr), intent(in) :: tr
      real(s2_sp), intent(out) :: norm(:)
      logical, intent(in), optional :: weight_in

      logical, parameter :: WEIGHT_DEFAULT = .false.
      logical :: weight = WEIGHT_DEFAULT
      integer :: i_dil, i_beta
      real(s2_sp) :: beta, weight_val
      type(cswt_tr) :: tr_temp

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_norm')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_norm')
      end if

      ! Check mfactor array of correct size.
      if(size(norm) /= tr%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_norm', &
           comment_add='Norm array of invalid size')
      end if

      ! Set weight option.
      if(present(weight_in)) then
         weight = weight_in
      end if

      ! Make temporary copy of tr structure to manipulate herein.
      tr_temp = cswt_tr_init(tr)

      ! Square.
      tr_temp%wcoeff = tr_temp%wcoeff**2.0e0

      ! Weight if required, do all dilations at once.
      if(weight) then
         do i_beta = 0,tr_temp%n_beta - 1
            beta = real(i_beta,s2_sp)/real(tr_temp%n_beta, s2_sp) * pi
            weight_val = 2*pi**2.0e0 &
                 / real(tr_temp%n_alpha * tr_temp%n_beta) * sin(beta)               
            tr_temp%wcoeff(:,:,i_beta,:) = tr_temp%wcoeff(:,:,i_beta,:) &
                 * weight_val
         end do
      end if
      
      ! Sum for each dilation.
      do i_dil = 1,tr_temp%n_dilation
          norm(i_dil) = sum(tr_temp%wcoeff(i_dil,:,:,:))       
      end do

      ! Scale by one over n_samples per dilation.
!      norm = norm &
!        / real(tr_temp%n_alpha * tr_temp%n_beta * tr_temp%n_gamma, s2_sp)

      ! Free memory.
      call cswt_tr_free(tr_temp)

    end subroutine cswt_tr_norm


    !--------------------------------------------------------------------------
    ! cswt_tr_multiply
    !
    !! Multiply two wavelet coefficients contained in tr structure arrays. 
    !! Weight the multiplication for each pixel by the size of the pixel on 
    !! the sky if desired.  The product coefficients are written to the 
    !! tr1 structure wavelet coefficient array on exit.
    !!  
    !! Variables:
    !!   - tr1: First tr structure containing the coefficients to multiply.
    !!     Contains product of coefficients on exit.
    !!   - tr2: Second tr structure containing the coefficients to multiply.
    !!     Unchanged on exit.
    !!   - [weight_in]: Logical to specified whether pixel weights are to be 
    !!     applied.  If not specified then WEIGHT_DEFAULT is assumed.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_multiply(tr1, tr2, weight_in)

      type(cswt_tr), intent(inout) :: tr1
      type(cswt_tr), intent(in) :: tr2
      logical, intent(in), optional :: weight_in

      logical, parameter :: WEIGHT_DEFAULT = .false.
      logical :: weight = WEIGHT_DEFAULT
      real(s2_sp) :: beta, weight_val
      integer :: i_beta

      ! Check objects initialised.
      if(.not. tr1%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_multiply')
      end if
      if(.not. tr2%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_multiply')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr1%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_multiply')
      end if
      if(.not. tr2%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_multiply')
      end if

      ! Check sizes of each wavelet coefficient structure are the same.
      if(tr1%n_alpha /= tr2%n_alpha &
         .or. tr1%n_beta /= tr2%n_beta &
         .or. tr1%n_gamma /= tr2%n_gamma &
         .or. tr1%n_dilation /= tr2%n_dilation) then
           call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_multiply', &
             comment_add='Wavelet coefficient sizes inconsistent')
      end if

      ! Set weight option.
      if(present(weight_in)) then
         weight = weight_in
      end if

      ! Perform multiplication.
      if(weight) then
         
         do i_beta = 0,tr1%n_beta - 1

            beta = real(i_beta,s2_sp)/real(tr1%n_beta,s2_sp) * pi

            weight_val = 2*pi**2.0e0 / real(tr1%n_alpha * tr1%n_beta) &
              * sin(beta)

            tr1%wcoeff(:,:,i_beta,:) &
              = tr1%wcoeff(:,:,i_beta,:) &
                * tr2%wcoeff(:,:,i_beta,:) &
                * weight_val

         end do
  
      else

         tr1%wcoeff = tr1%wcoeff * tr2%wcoeff

      end if

    end subroutine cswt_tr_multiply


    !--------------------------------------------------------------------------
    ! cswt_tr_const_multiply_array
    !
    !! Multiply the wavelet coefficients of the tr structure for a given
    !! dilation by the value of mfactor for the corresponding dilation.
    !!  
    !! Variables:
    !!   - tr: tr structure containing the coefficients to multiply.
    !!   - mfactor: Array containing multiplication factor for each dilation.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_const_multiply_array(tr, mfactor)

      type(cswt_tr), intent(inout) :: tr
      real(s2_sp), intent(in) :: mfactor(:)

      integer :: i_dil

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_const_multiply_array')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_const_multiply_array')
      end if

      ! Check mfactor array of correct size.
      if(size(mfactor) /= tr%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_const_multiply_array', &
           comment_add='Multiplication factor array of invalid size')
      end if

      ! Perform multiplication for each dilation.
      do i_dil = 1,tr%n_dilation
         tr%wcoeff(i_dil,:,:,:) = tr%wcoeff(i_dil,:,:,:) * mfactor(i_dil)
      end do

    end subroutine cswt_tr_const_multiply_array


    !--------------------------------------------------------------------------
    ! cswt_tr_const_multiply_val
    !
    !! Multiply the wavelet coefficients of the tr structure by a constant 
    !! value (for all dilations).
    !!  
    !! Variables:
    !!   - tr: tr structure containing the coefficients to multiply.
    !!   - mfactor: Constant value to multiply all coefficients by (for all
    !!     dilations).
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_const_multiply_val(tr, mfactor)
      
      type(cswt_tr), intent(inout) :: tr
      real(s2_sp), intent(in) :: mfactor

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_const_multiply_val')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_const_multiply_val')
      end if

      ! Perform multiplication.
      tr%wcoeff = tr%wcoeff * mfactor

    end subroutine cswt_tr_const_multiply_val


    !--------------------------------------------------------------------------
    ! cswt_tr_const_subtract_array
    !
    !! Subtract the wavelet coefficients of the tr structure for a given
    !! dilation by the value of val for the corresponding dilation.
    !!  
    !! Variables:
    !!   - tr: tr structure containing the coefficients to subtract.
    !!   - val: Array containing subtraction values for each dilation.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_const_subtract_array(tr, val)

      type(cswt_tr), intent(inout) :: tr
      real(s2_sp), intent(in) :: val(:)

      integer :: i_dil

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_const_subtract_array')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_const_subtract_array')
      end if

      ! Check val array of correct size.
      if(size(val) /= tr%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_const_subtract_array', &
           comment_add='Subtract value array of invalid size')
      end if

      ! Perform multiplication for each dilation.
      do i_dil = 1,tr%n_dilation
         tr%wcoeff(i_dil,:,:,:) = tr%wcoeff(i_dil,:,:,:) - val(i_dil)
      end do

    end subroutine cswt_tr_const_subtract_array


    !--------------------------------------------------------------------------
    ! cswt_tr_const_subtract_val
    !
    !! Subtract the wavelet coefficients of the tr structure by a constant 
    !! value (for all dilations).
    !!  
    !! Variables:
    !!   - tr: tr structure containing the coefficients to subtract.
    !!   - val: Constant value to subtract from all coefficients (for all
    !!     dilations).
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_const_subtract_val(tr, val)
      
      type(cswt_tr), intent(inout) :: tr
      real(s2_sp), intent(in) :: val

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_const_subtract_val')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_const_subtract_val')
      end if

      ! Perform multiplication.
      tr%wcoeff = tr%wcoeff - val

    end subroutine cswt_tr_const_subtract_val


    !--------------------------------------------------------------------------
    ! cswt_tr_wcoeff_sigma_each
    !
    !! Compute the mean and standard deviation of the wavelet coefficients of
    !! the tr strucutre *for each dilation*.
    !!  
    !! Notes:
    !!   - If tr wavelet coefficients masked, sigma and mean should be
    !!     computed only over the remaining pixels.  If these pixels are zero,
    !!     then they won't compute to the calculation of mean and sigma, 
    !!     however the number of samples to divide by should be the effective
    !!     number of pixels.  Thus divide by this value, for each dilation, 
    !!     that is passed as an input.
    !!
    !! Variables:
    !!   - tr: tr structure containing the coefficients to calculate the mean
    !!     and standard deviation of.
    !!   - n_effective_samples: Number of effective samples per dilation
    !!     for previously applied mask, if not masked just give full number
    !!     of Euler angle samples (for each dilation - which of course will
    !!     be an array of the same number (n_alpha*n_beta*n_gamma)).
    !!   - mean: Computed mean array for each dilation.
    !!   - sigma: Computed standard deviation array for each dilation.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_wcoeff_sigma_each(tr, n_effective_samples, mean, sigma)

      type(cswt_tr), intent(in) :: tr
      real(s2_sp), intent(out) :: mean(:), sigma(:)
      integer, intent(in) :: n_effective_samples(:)

      integer :: i_dil

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_wcoeff_sigma')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_wcoeff_sigma')
      end if

      ! Check output arrays of correct size.
      if(size(sigma) /= tr%n_dilation .or. size(mean) /= tr%n_dilation &
           .or. size(n_effective_samples) /= tr%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_wcoeff_sigma', &
           comment_add='Output arrays of invalid size')
      end if

      ! Compute mean and sigma for each dilation
      do i_dil = 1,tr%n_dilation
         
         mean(i_dil) = sum(tr%wcoeff(i_dil,:,:,:)) &
            / real(n_effective_samples(i_dil), s2_sp)

         sigma(i_dil) = &
           sqrt( sum( (tr%wcoeff(i_dil,:,:,:) - mean(i_dil))**2.0e0 )  &
                 / real(n_effective_samples(i_dil), s2_sp) )

      end do

    end subroutine cswt_tr_wcoeff_sigma_each

    
    !--------------------------------------------------------------------------
    ! cswt_tr_wcoeff_sigma_all
    !
    !! Compute the mean and standard deviation of the wavelet coefficients of
    !! the tr strucutre *over all dilations*.
    !!  
    !! Notes:
    !!   - If tr wavelet coefficients masked, sigma and mean should be
    !!     computed only over the remaining pixels.  If these pixels are zero,
    !!     then they won't compute to the calculation of mean and sigma, 
    !!     however the number of samples to divide by should be the effective
    !!     number of pixels.  Thus divide by this value that is passed as an 
    !!     input.
    !!
    !! Variables:
    !!   - tr: tr structure containing the coefficients to calculate the mean
    !!     and standard deviation of.
    !!   - n_effective_samples: Number of effective samples, if not masked 
    !!     just give full number of samples: size(tr%wcoeff).
    !!   - mean: Computed mean over all dilations.
    !!   - sigma: Computed standard deviation over all dilations.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_wcoeff_sigma_all(tr, n_effective_samples, mean, sigma)

      type(cswt_tr), intent(in) :: tr
      integer, intent(in) :: n_effective_samples
      real(s2_sp), intent(out) :: mean, sigma

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_wcoeff_sigma_all')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_wcoeff_sigma_all')
      end if

      ! Compute mean and sigma over all coefficients.

      mean = sum(tr%wcoeff(:,:,:,:)) &
            / real(n_effective_samples, s2_sp)

      sigma = sqrt( sum( (tr%wcoeff(:,:,:,:) - mean)**2.0e0 )  &
                 / real(n_effective_samples, s2_sp) )

    end subroutine cswt_tr_wcoeff_sigma_all


    !--------------------------------------------------------------------------
    ! cswt_tr_wcoeff_thres_nsigma
    !
    !! Threshold wavelet coefficients.  Coefficients that remain keep their 
    !! original value and are *not* set to 1.0.  Differs to 
    !! cswt_tr_mask_thres_nsigma routine since the same nsigma is used over
    !! all dilations.
    !!
    !! Notes:
    !!   - If tr wavelet coefficients masked, sigma and mean should be
    !!     computed only over the remaining pixels.  If these pixels are zero,
    !!     then they won't compute to the calculation of mean and sigma, 
    !!     however the number of samples to divide by should be the effective
    !!     number of pixels.  Thus divide by this value that is passed as 
    !!     an input.
    !!
    !! Variables:
    !!   - tr: tr structure containing the coefficients to be thresholded.
    !!     On exit, tr structure containing thresholded coefficients.
    !!   - n_effective_samples: Number of effective samples, if not masked 
    !!     just give full number of samples: size(tr%wcoeff).
    !!   - mode: Mode of thresholding to perform.  
    !!     If mode=CSWT_TR_THRES_NSIGMA_MODE_ABS then the absolute value of
    !!     coefficients are considered and those coefficients below (in
    !!     absolute value)  nsigma * sigma are set to zero.
    !!     If mode=CSWT_TR_THRES_NSIGMA_MODE_ABOVE then those coefficients 
    !!     below nsigma * sigma are set to zero.
    !!     If mode=CSWT_TR_THRES_NSIGMA_MODE_BELOW then those coefficients 
    !!     above -nsigma * sigma (*note minus sign*) are set to zero.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_wcoeff_thres_nsigma(tr, nsigma, mode, &
         n_effective_samples_in)

      type(cswt_tr), intent(inout) :: tr
      real(s2_sp) :: nsigma
      integer, intent(in) :: mode
      integer, intent(in), optional :: n_effective_samples_in
      
      real(s2_sp) :: mean, sigma, thres
      integer :: n_effective_samples

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_wcoeff_thres_nsigma')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_wcoeff_thres_nsigma')
      end if

      ! Set size of effective samples.
      if(present(n_effective_samples_in)) then
         n_effective_samples = n_effective_samples_in
      else
         n_effective_samples = size(tr%wcoeff)
      end if

      ! Compute mean and sigma over all dilations.
      call cswt_tr_wcoeff_sigma_all(tr, n_effective_samples, mean, sigma)

      ! Compute threshold to use.
      thres = nsigma * sigma

      ! Threshold wavelet coefficients depending on mode.
      ! Note difference from mean considered, not simply value 
      ! (although generally the mean will be zero).
      select case (mode)

           case (CSWT_TR_THRES_NSIGMA_MODE_ABS)
              where( abs(tr%wcoeff(:,:,:,:) - mean) < thres )
                 tr%wcoeff(:,:,:,:) = 0.0e0
              end where

             case (CSWT_TR_THRES_NSIGMA_MODE_ABOVE)
                where( (tr%wcoeff(:,:,:,:) - mean) < thres )
                   tr%wcoeff(:,:,:,:) = 0.0e0
                end where

           case (CSWT_TR_THRES_NSIGMA_MODE_BELOW)
              where( (tr%wcoeff(:,:,:,:) - mean) > -thres )
                 tr%wcoeff(:,:,:,:) = 0.0e0
              end where

           case default
              call cswt_error(CSWT_ERROR_TR_THRESNSIG_MODE_INVALID, &
                'cswt_tr_wcoeff_thres_nsigma')
         end select

    end subroutine cswt_tr_wcoeff_thres_nsigma


    !--------------------------------------------------------------------------
    ! Local optimisation routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tr_localmax_ab
    !
    !! Find local maxima in a single ab slice of the wavelet coefficients
    !! (performed for each orientation separately).
    !! The input coefficients should already be thresholded so that only local
    !! non-zero regions remain for which the maximum of each is found.
    !!
    !! Notes:
    !!   - Note the tr structure passed should be thresholded so only need 
    !!     to find local maximum in each connected component region.
    !!   - Data arrays of maximum local values and locations are allocated
    !!     herein, since no prior knowledge of number of objects contained 
    !!     in map.  These arrays *must* be freed by the calling routine.
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!   - The algorithm works by finding the connected components in 
    !!     thresholded wavelet coefficient map.  The maximum of each 
    !!     connected component is taken as the local maximum.
    !!
    !! Variables:
    !!   - tr: *Thresholded* tr structure containing non-zero local regions
    !!     to find maxima of.
    !!   - i_dil: Dilation index of tr%wcoeff map to consider.
    !!   - i_gamma: Gamma Euler angle index of tr%wcoeff to consider 
    !!     (indexed from 1:N).
    !!   - max_val: List of local maximum values found.
    !!   - max_loc: List of locations of local maximum values found.
    !!   - [filename_connected]: If present save cswt file of connected 
    !!     components.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !   April 2006 - Jason McEwen: Updated to search each orientation
    !--------------------------------------------------------------------------

    subroutine cswt_tr_localmax_ab(tr, i_dil, n_regions, max_val, max_loc, &
         max_siz, filename_connected)

      type(cswt_tr), intent(in) :: tr
      integer, intent(in) :: i_dil
      integer, allocatable, intent(out) :: n_regions(:)
      real(s2_sp), allocatable, intent(out) :: max_val(:,:)
      integer, allocatable, intent(out) :: max_loc(:,:,:)
      integer, allocatable, intent(out) :: max_siz(:,:)
      character(len=*), intent(in), optional :: filename_connected

      real(s2_sp), parameter :: TOL  = 1e-10
      type(cswt_tr) :: tr_connected
      integer :: i_gamma, i_gamma0, i_r, fail, max_n_regions

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_localmax_ab')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_localmax_ab')
      end if

      ! Check i_dil and i_gamma are in valid range.
      if(i_dil < 1 .or. i_dil > tr%n_dilation ) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_localmax_ab', &
           comment_add='Dilation index invalid')
      end if

      ! Initialise connected component data structure.
      tr_connected = cswt_tr_init(tr)     

      ! Allocate memory for n_regions for each orientation.
      allocate(n_regions(0:tr_connected%n_gamma-1), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
              'cswt_tr_localmax_ab')
      end if
      
      ! Find connected components for each orientation
      ! and count number of regions in each.
      do i_gamma = 1,tr_connected%n_gamma

         ! Set i_gamma0 indexed from 0, not 1 like i_gamma.
         i_gamma0 = i_gamma - 1
         
         ! Find connected components in ab for this orientation.
         call cswt_tr_connected_ab(tr_connected, i_dil, i_gamma)
         
         ! Get number of connected components, i.e. number of objects detected.
         n_regions(i_gamma0) = maxval(tr_connected%wcoeff(i_dil,:,:,i_gamma0))
         
      end do

      ! Determine max number of regions over all orientations.
      max_n_regions = maxval(n_regions)

      ! Allocate memory for max values and locations.
      allocate(max_val(0:tr_connected%n_gamma-1,1:max_n_regions), stat=fail)
      allocate(max_loc(0:tr_connected%n_gamma-1,1:max_n_regions, 2), stat=fail)
      allocate(max_siz(0:tr_connected%n_gamma-1,1:max_n_regions), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
              'cswt_tr_localmax_ab')
      end if

      ! Find maximum values and locations for each orientation.
      do i_gamma = 1,tr_connected%n_gamma

         ! Set i_gamma0 indexed from 0, not 1 like i_gamma.
         i_gamma0 = i_gamma - 1

         ! Find maximum values and locations in each connected component.
         do i_r = 1,n_regions(i_gamma0)
         
            max_val(i_gamma0,i_r) = maxval(tr%wcoeff(i_dil,:,:,i_gamma0), &
                 mask=abs(tr_connected%wcoeff(i_dil,:,:,i_gamma0) - i_r) < TOL)
            
            max_loc(i_gamma0,i_r,:) = maxloc(tr%wcoeff(i_dil,:,:,i_gamma0), &
                 mask=abs(tr_connected%wcoeff(i_dil,:,:,i_gamma0) - i_r)<TOL) 
            ! Note locations indexed from 1, even though alpha, beta
            ! components indexed from 0 in array.

            max_siz(i_gamma0,i_r) = sum(tr_connected%wcoeff(i_dil,:,:,i_gamma0)-i_r+1, &
                 mask=abs(tr_connected%wcoeff(i_dil,:,:,i_gamma0) - i_r)<TOL)

         end do

      end do

      ! Shift locations so indexed from zero.
      max_loc = max_loc - 1

      ! Save connected components if filename_connected present.
      if(present(filename_connected)) then
         call cswt_tr_io_fits_write_wcoeff(filename_connected, tr_connected)
      end if

      ! Free memory.
      call cswt_tr_free(tr_connected)

    end subroutine cswt_tr_localmax_ab


    !--------------------------------------------------------------------------
    ! cswt_tr_connected_ab
    !
    !! Find wavelet coefficient connected components (defined by non-zero 
    !! pixels) and overwrite with integer label for each connected component.
    !!
    !! Notes:
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!   - Connected components defined by non-zero pixels (non-zero values
    !!     are irrelevant).
    !!
    !! Variables:
    !!   - tr: *Thresholded* tr structure to find connected components of.
    !!     On exit, coefficients contain integer values (although real 
    !!      precision) labelling connected components.
    !!   - i_dil: Dilation index of tr%wcoeff map to consider.
    !!   - i_gamma: Gamma Euler angle index of tr%wcoeff to consider 
    !!     (indexed from 1:N).
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_connected_ab(tr, i_dil, i_gamma)

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: i_dil, i_gamma

      real(s2_sp), parameter :: TOL = 1e-10
      integer :: i, j, ia, ib, isym, i_gamma0, n_equiv_pair, fail
      real(s2_sp) :: l, t
      integer :: current_new_label

      integer, allocatable :: equiv_tab(:,:)
      integer, allocatable :: equiv_class(:,:)
      
      ! Initialise variables (if init above then become 'saved/persistent'
      ! variables).
      n_equiv_pair = 0
      l = 0e0
      t = 0e0

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_connected_ab')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_connected_ab')
      end if

      ! Check i_dil and i_gamma are in valid range.
      if(i_dil < 1 .or. i_dil > tr%n_dilation &
        .or. i_gamma < 1 .or. i_gamma > tr%n_gamma) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_connected_ab', &
           comment_add='Dilation or orientation index invalid')
      end if

      ! Set i_gamma0 (gamma indexed from zero).
      i_gamma0 = i_gamma - 1

      ! Initialise current label.
      current_new_label = 0

      ! Allocate initial space for table of equivalence pairs.
      allocate(equiv_tab(1:200,1:2), stat=fail)      
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
              'cswt_tr_connected_ab')
      end if

      do ib = 0,tr%n_beta-1
         do ia = 0,tr%n_alpha-1

            ! If current point non-zero then set label (either copied or new).
            if(abs(tr%wcoeff(i_dil,ia,ib,i_gamma0)) > TOL) then

               ! Get left and top labels.

               if(ia==0) then
! Set left label as undefined for first column.
!l = 0.0e0
                  ! Wrap left label around sphere and go up one row since current
                  !row not yet considered.
                  if(ib==0) then
                     ! If at very beginning then set l to zero.
                     l = 0.0e0 
                  else
                     l = tr%wcoeff(i_dil, tr%n_alpha-1, ib-1, i_gamma0)
                  end if
               else
                  l = tr%wcoeff(i_dil, ia-1, ib, i_gamma0)
               end if

               if(ib==0) then
                  ! Set top label as undefined for first row.
                  t = 0.0e0
               else
                  t = tr%wcoeff(i_dil, ia, ib-1, i_gamma0)
               end if

               ! Set label of current pixel.
               if( abs(l) > TOL .and. abs(t) > TOL .and. abs(t-l) < TOL ) then
                  
                  ! Both non-zero and the same.

                  ! Copy (either) label.
                  tr%wcoeff(i_dil,ia,ib,i_gamma0) = t

               elseif( abs(l) > TOL .and. abs(t) > TOL .and. abs(t-l) > TOL ) then

                  ! Both non-zero and different.

                  ! Copy top label.
                  tr%wcoeff(i_dil,ia,ib,i_gamma0) = t

                  ! Set labels as equivalent.

                  ! If either l or m different to current row then add, else already in table.
                  if(abs(equiv_tab(n_equiv_pair,1) - l) > TOL &
                       .or. abs(equiv_tab(n_equiv_pair,2) - t) > TOL) then
                     n_equiv_pair = n_equiv_pair + 1
                     if(n_equiv_pair > size(equiv_tab,1)) then
                        call cswt_tr_reallocate_mat(equiv_tab, .true.)
                     end if
                     equiv_tab(n_equiv_pair,1) = int(l)
                     equiv_tab(n_equiv_pair,2) = int(t)
                  end if

               elseif( abs(t) > TOL) then

                  ! If only top label defined then copy.
                  tr%wcoeff(i_dil,ia,ib,i_gamma0) = t

               elseif( abs(l) > TOL) then

                  ! if only left label defined then copt.
                  tr%wcoeff(i_dil,ia,ib,i_gamma0) = l

               else
                  
                  ! If no neighbour labels defined then set new label.
                  current_new_label = current_new_label + 1
                  tr%wcoeff(i_dil,ia,ib,i_gamma0) = real(current_new_label, s2_sp)

               end if

            end if

         end do
      end do

      ! Add symmetric cases to the table to ensure they aren't missed.
      ! (I.e. they will just appear in their own equivalence class 
      ! once classes constructed from table.)
      do isym = 1,current_new_label
         if(isym+n_equiv_pair > size(equiv_tab,1)) then
            call cswt_tr_reallocate_mat(equiv_tab, .true.)
         end if
         equiv_tab(isym+n_equiv_pair,1) = isym
         equiv_tab(isym+n_equiv_pair,2) = isym
      end do

      ! Resolve equivalence pairs into equivalence classes.
      call cswt_tr_equivclass(equiv_tab, n_equiv_pair+current_new_label, &
        equiv_class)

      ! Relabel equivalence classes.
      ! Set new values as negative for now to ensure not replaced in 
      ! subsequent replacement step.  (No negative labels used in 
      ! equivalence pairs.) Overwrite with positive values once 
      ! finished relabeling.
      i = 1
      j = 1
      do while(equiv_class(i,j) /= 0)
         do while(equiv_class(i,j) /= 0)
            
            where( abs(tr%wcoeff(i_dil,:,:,i_gamma0)-equiv_class(i,j)) < TOL )
               tr%wcoeff(i_dil,:,:,i_gamma0) = -real(i,s2_sp)
            end where
            
            j = j + 1            
         end do 
         i = i + 1
         j = 1
      end do

      ! Convert sign of relabeled labels to give positive values.
      tr%wcoeff(i_dil,:,:,i_gamma0) = -tr%wcoeff(i_dil,:,:,i_gamma0)

      ! Free temporary memory used to store equivalence pairs and resolved
      ! class table.
      deallocate(equiv_tab, equiv_class)

    end subroutine cswt_tr_connected_ab


    !--------------------------------------------------------------------------
    ! cswt_tr_equivclass
    !
    !! Resolve a table of equivalence pairs into a table of equivalence 
    !! classes.
    !!
    !! Notes:
    !!   - Note memory for equic is allocated here and *must* be freed be 
    !!     calling routine.
    !!
    !! Variables:
    !!   - tab: Original table of equivalence pairs.  Table may be larger 
    !!     than actual number of pairs - unused cells filled with zeros.
    !!   - n_equiv_pair: Number of equivalence pairs in tab (number of 
    !!     non-zero rows of table).
    !!   - equiv: Table of equivalence classes.  Table may be larger than 
    !!     required - unused cells filled with zeros.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_equivclass(tab, n_equiv_pair, equiv)

      integer, intent(in) :: tab(:,:)
      integer, intent(in) :: n_equiv_pair
      integer, allocatable, intent(out) :: equiv(:,:)

      integer :: M, ieq, ieq_test, jeq, i, j, test_val, fail

      ! Allocate initial space for equivalence class table.
      allocate(equiv(1:50,1:5), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
              'cswt_tr_equivclass')
      end if

      equiv = 0

      M = n_equiv_pair
    
      ieq = 0

      do ieq_test = 1,M

         jeq = 1
         j = 1
    
         test_val = tab(ieq_test,1)
 
         ! Continue if test_val not in equiv array already.
         ! If i value not then j value won't be since already checked.
         if( .not. cswt_tr_present_mat_int(equiv, test_val) ) then

            ieq = ieq + 1

            ! Reallocate more rows if needed.
            if(ieq > size(equiv,1)) then
               call cswt_tr_reallocate_mat(equiv, .true.)
            end if

            equiv(ieq, jeq) = test_val

            do while(j <= jeq) 

               do i = ieq,M
                  
                  if(tab(i,1) == equiv(ieq, j)) then
                    if( .not. cswt_tr_present_vect_int(equiv(ieq,:), tab(i,2)) ) then
                       jeq = jeq + 1
                       ! Reallocate more columns if needed.
                       if(jeq > size(equiv,2)) then                       
                          call cswt_tr_reallocate_mat(equiv, .false.)
                       end if
                       equiv(ieq, jeq) = tab(i,2)
                     end if
                  end if
                  
                  if(tab(i,2) == equiv(ieq, j)) then
                     if( .not. cswt_tr_present_vect_int(equiv(ieq,:), tab(i,1)) ) then
                        jeq = jeq + 1
                        ! Reallocate more columns if needed.
                        if(jeq > size(equiv,2)) then
                           call cswt_tr_reallocate_mat(equiv, .false.)
                        end if
                        equiv(ieq, jeq) = tab(i,1)
                     end if
                  end if
                  
               end do
               
               j = j + 1

            end do

         end if
            
      end do

    end subroutine cswt_tr_equivclass
   

    !--------------------------------------------------------------------------
    ! cswt_tr_reallocate_mat
    !
    !! Dynamically reallocate a 2D array to increase its size, either in the x 
    !! or y direction.
    !!
    !! Notes:
    !!   - Original mat must be allocatable.
    !!   - Current implementation is only for integer mat.
    !!
    !! Variables:
    !!   - mat: Original matrix to increase the size of.  On exit, mat 
    !!     contains the same data but has a larger size.
    !!   - x_axis: Logical indication whether to increase the size of mat in
    !!     the x (down) direction.  If false then size is increased in the y
    !!     (across) direction.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_reallocate_mat(mat, x_axis)

      integer, allocatable :: mat(:,:)
      logical, intent(in) :: x_axis 

      integer, parameter :: INC_X = 50, INC_Y = 5
      integer, allocatable :: temp_mat(:,:)
      integer :: n_x, n_y, fail

      n_x = size(mat,1)
      n_y = size(mat,2)
      
      ! Create temporary copy of values stored in mat to replace later.
      allocate(temp_mat(1:n_x,1:n_y), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
              'cswt_tr_reallocate_mat')
      end if
      temp_mat = mat

      ! Reallocate larger memory for mat and initialise values to zero.
      deallocate(mat)
      if(x_axis) then
         allocate(mat(n_x+INC_X, n_y), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
                 'cswt_tr_reallocate_mat')
         end if
      else
         allocate(mat(n_x, n_y+INC_Y), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
                 'cswt_tr_reallocate_mat')
         end if
      end if
      mat = 0.0e0

      ! Copy values back into resized original matrix.
      mat(1:n_x,1:n_y) = temp_mat

      ! Free temporary memory used.
      deallocate(temp_mat)

    end subroutine cswt_tr_reallocate_mat


    !--------------------------------------------------------------------------
    ! cswt_tr_present_vect_int
    !
    !! Check if a value is present in a vector.
    !!
    !! Notes:
    !!   - Could simply subtract value from vect, then take product of all
    !!     elements.  Result will be zero if value is present.  Errors with 
    !!     this approach however (error with product intrinsic function?), so
    !!     search directly.
    !!
    !! Variables:
    !!   - vect: Vector to search for value in.
    !!   - val: Test value to look for in vector.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_present_vect_int(vect, val) result(pres)

      integer, intent(in) :: vect(:)
      integer, intent(in) :: val
      logical :: pres

      integer :: i

      pres = .false.

      do i = 1,size(vect)

         if(vect(i) == val) then
            pres = .true.
            return
         end if

      end do

    end function cswt_tr_present_vect_int


    !--------------------------------------------------------------------------
    ! cswt_tr_present_mat_int
    !
    !! Check if a value is present in a 2D matrix.
    !!
    !! Notes:
    !!   - Could simply subtract value from mat, then take product of all
    !!     elements.  Result will be zero if value is present.  Errors with 
    !!     this approach however (error with product intrinsic function?), so
    !!     search directly.
    !!
    !! Variables:
    !!   - mat: matrix to search for value in.
    !!   - val: Test value to look for in vector.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - April 2005
    !
    ! Revisions:
    !   April 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_present_mat_int(mat, val) result(pres)

      integer, intent(in) :: mat(:,:)
      integer, intent(in) :: val
      logical :: pres

      integer :: i, j

      pres = .false.

      do i = 1,size(mat,1)
         do j = 1,size(mat,2)

            if(mat(i,j) == val) then
               pres = .true.
               return
            end if

         end do
      end do

    end function cswt_tr_present_mat_int


    !--------------------------------------------------------------------------
    ! Masking routines
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    ! cswt_tr_mask_nonzero
    !
    !! Compute the number of non-zero values (i.e. effective remaining 
    !! coefficients) in a given mask.  (Computed for each dilation and
    !! orientation, i.e. for each alpha-beta `sky' of the mask.)
    !!  
    !! Variables:
    !!   - tr_mask: tr structure containing the mask in the tr_mask%wcoeff
    !!     array.
    !!   - ncoeff_eff: Number of effective coefficients in the mask for each
    !!     dilation and orientation, i.e. for each alpha-beta `sky'.
    !!   - ncoeff_tot: Total number of coefficients in a single alpha-beta
    !!     mask.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_nonzero(tr_mask, ncoeff_eff, ncoeff_tot)

      type(cswt_tr), intent(in) :: tr_mask
      integer, intent(out) :: ncoeff_eff(:,:)
      integer, intent(out) :: ncoeff_tot

      integer :: i_dil, i_alpha, i_beta, i_gamma
      real(s2_sp), parameter :: ZERO_TOL = 0.01

      ! Check object initialised.
      if(.not. tr_mask%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_nonzero')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr_mask%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_mask_nonzero')
      end if

      ! Check input ncoeff_eff is correct size.
      if(size(ncoeff_eff,1) /= tr_mask%n_dilation .or. &
       size(ncoeff_eff,2) /= tr_mask%n_gamma) then
       call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_mask_nonzero', &
        comment_add='Input effective number of coefficient array invalid size')
      end if
      ncoeff_eff = 0

      ! Calculate the effective (non-zero) number of coefficients in mask.
      do i_dil = 1,tr_mask%n_dilation
         do i_gamma = 0,tr_mask%n_gamma-1

            ! Count number of non-zero values.
            do i_alpha = 0,tr_mask%n_alpha-1
               do i_beta = 0,tr_mask%n_beta-1
                  
                  if(abs(tr_mask%wcoeff(i_dil,i_alpha,i_beta,i_gamma)) &
                    > ZERO_TOL) then
                     ncoeff_eff(i_dil, i_gamma+1) = &
                       ncoeff_eff(i_dil, i_gamma+1) + 1
                  end if

               end do
            end do

            ! Previous technique assumes values are ones.
            !            ncoeff_eff(i_dil, i_gamma+1) &
            !              = sum(tr_mask%wcoeff(i_dil,:,:,i_gamma))
            
         end do
      end do
      
      ! Calculate the total number of coefficients in the mask.
      ncoeff_tot = tr_mask%n_alpha * tr_mask%n_beta

    end subroutine cswt_tr_mask_nonzero

 
    !--------------------------------------------------------------------------
    ! cswt_tr_mask_nonzero_weight
    !
    !! Compute the proportion of weighted coefficients in a given mask. 
    !! (Computed for each dilation and orientation, i.e. for each alpha-beta 
    !! `sky' of the mask.)
    !! 
    !! Notes:
    !!   -  The proportion is computed from
    !!      .  \sum_{all i s.t. mask_i > thres} w_i
    !!      .  ------------------------------------
    !!      .            \sum_{all j} w_j
    !!      for the map corresponding to each dilation and orientation.
    !! 
    !! Variables:
    !!   - tr_mask: tr structure containing the mask in the tr_mask%wcoeff
    !!     array.
    !!   - prop: Proportion of effective weighted coefficients for each
    !!     dilation and orientation, i.e. for each alpha-beta `sky'.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_nonzero_weight(tr_mask, prop)

      type(cswt_tr), intent(in) :: tr_mask
      real(s2_sp), intent(out) :: prop(:,:)

      integer :: i_dil, i_alpha, i_beta, i_gamma
      real(s2_sp), parameter :: ZERO_TOL = 0.01
      real(s2_sp) :: beta, weight, mask_sum, weight_sum

      ! Check object initialised.
      if(.not. tr_mask%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_nonzero_weight')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr_mask%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_mask_nonzero_weight')
      end if

     ! Check input ncoeff_eff is correct size.
      if(size(prop,1) /= tr_mask%n_dilation .or. &
       size(prop,2) /= tr_mask%n_gamma) then
       call cswt_error(CSWT_ERROR_SIZE_INVALID, &
         'cswt_tr_mask_nonzero_weight', &
         comment_add='Size of input proportion array invalid')
      end if
      prop = 0.0e0

      ! Calculate the weighted proportion of mask that remains.
      do i_dil = 1,tr_mask%n_dilation
         do i_gamma = 0,tr_mask%n_gamma-1

            weight_sum = 0
            mask_sum = 0

            ! Count number of non-zero values.
            do i_alpha = 0,tr_mask%n_alpha-1
               do i_beta = 0,tr_mask%n_beta-1

                  beta = i_beta/real(tr_mask%n_beta, s2_sp) * pi    
                  weight = 2*pi**2.0e0/real(tr_mask%n_alpha * tr_mask%n_beta) &
                       * sin(beta)
                  
                  weight_sum = weight_sum + weight

                  if(abs(tr_mask%wcoeff(i_dil,i_alpha,i_beta,i_gamma)) &
                    > ZERO_TOL) then
                     mask_sum = mask_sum + weight
                  end if

               end do
            end do

            ! Weight_sum should be the same each time but easiest to computed 
            ! again for each loop. 
            prop(i_dil, i_gamma+1) = mask_sum / weight_sum

         end do
      end do

    end subroutine cswt_tr_mask_nonzero_weight


    !--------------------------------------------------------------------------
    ! cswt_tr_mask_invert
    !
    !! Invert a coefficient mask by converting ones to zeros, and zeros to 
    !! ones.
    !!
    !! Variables:
    !!   - tr_mask: Mask to invert.  On output contains ones in the place of
    !!     original zeros, and zeros in the place of original ones.
    !!
    !! @author J. D. McEwen
    !! @version 0.1 - January 2005
    !
    ! Revisions:
    !   January 2005 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_invert(tr_mask)

      type(cswt_tr), intent(inout) :: tr_mask

      real(s2_sp), parameter :: ZERO_TOL = 1e-4

      ! Check object initialised.
      if(.not. tr_mask%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_apply')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr_mask%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_mask_apply')
      end if

      ! Invert mask (i.e. set 0s to 1s and 1s to 0s).
      where( abs(tr_mask%wcoeff) <= ZERO_TOL )
         tr_mask%wcoeff = 1.0e0
      elsewhere 
         tr_mask%wcoeff = 0.0e0
      end where

    end subroutine cswt_tr_mask_invert


    !--------------------------------------------------------------------------
    ! cswt_tr_mask_apply
    !
    !! Apply an extended coefficient mask by multiplying the wavelet
    !! coefficients of the specified tr structure with the mask (contained 
    !! in the wavelet coefficients of the tr_mask structure).
    !!  
    !! Notes:
    !!   - *Important*: Only specify the display flag if the output may is 
    !!     only used for display purposes.  If subsequent processing is 
    !!     applied erroneous results amy occur.
    !!   - Obviously the sizes of the tr and tr_mask structure must be the
    !!     same.
    !!
    !! Variables:
    !!   - tr: Wavelet transform structure to be masked.
    !!   - tr_mask: Mask to apply.  Mask is stored in the tr_mask%wcoeff array.
    !!   - [display]: Logical to specified whether a display mask is to be 
    !!     applied.  If that is the case then the pixels in the tr structure
    !!     corresponding to the zeros with be replaced with a large negative
    !!     value so that they appear grey in output images.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_apply(tr, tr_mask, display_in)

      type(cswt_tr), intent(inout) :: tr
      type(cswt_tr), intent(in) :: tr_mask
      logical, intent(in), optional :: display_in

      real(s2_sp), parameter :: ZERO_TOl = 0.1
      real(s2_sp), parameter :: FITS_DISPLAY_GREY_MAGIC_NUMBER = -1.6375e30
      logical, parameter :: DISPLAY_DEFAULT = .false.
      logical :: display = DISPLAY_DEFAULT

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_apply')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_mask_apply')
      end if

      ! Check object initialised.
      if(.not. tr_mask%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_apply')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr_mask%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_mask_apply')
      end if

      ! Check sizes consistent.
      if(tr%lmax /= tr_mask%lmax .or. &
         tr%mmax /= tr_mask%mmax .or. &
         tr%n_alpha /= tr_mask%n_alpha .or. &
         tr%n_beta /= tr_mask%n_beta .or. &
         tr%n_gamma /= tr_mask%n_gamma .or. &
         tr%n_dilation /= tr_mask%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_mask_apply', &
           comment_add='Transform structure and mask have inconsistent size')
      end if

      ! Set display status.
      if(present(display_in)) display = display_in

      ! Apply mask.
      tr%wcoeff = tr%wcoeff * tr_mask%wcoeff

      ! Produce display output so that masked regions appear grey, if required.
      if(display) then
         ! Overwrite mask zero region grey.
         where(tr_mask%wcoeff < ZERO_TOL)
            tr%wcoeff = FITS_DISPLAY_GREY_MAGIC_NUMBER
         end where
      end if

    end subroutine cswt_tr_mask_apply


    !--------------------------------------------------------------------------
    ! cswt_tr_mask_thres_nsigma
    !
    !! Produce a mask from the wavelet coefficients contained in the tr 
    !! structure by thresholding the coefficients to binary 1,0 values.
    !! A range of thresholding 'directions' may be performed (see mode
    !! below), nevertheless each is realitive to the nsigma * sigma(a)
    !! (where sigma(a) is computed herein).
    !! 
    !! Notes:
    !!   - Note that on output the binary {0,1} mask is stored in the tr
    !!     wavelet coefficient array. 
    !!   - If tr wavelet coefficients masked, sigma and mean should be
    !!     computed only over the remaining pixels.  If these pixels are zero,
    !!     then they won't compute to the calculation of mean and sigma, 
    !!     however the number of samples to divide by should be the effective
    !!     number of pixels.  Thus divide by this value, for each dilation, 
    !!     that is passed as an input.
    !!
    !! Variables:
    !!   - tr: Strucutre containing wavelet coefficients to be thresholded 
    !!     to construct the mask.  On exit, the binary mask is contained in
    !!     the tr structure wavelet coefficient array.
    !!   - nsigma: Number of sigma threshold is set at.
    !!   - n_effective_samples: Number of effective samples per dilation
    !!     for previously applied mask, if not masked just give full number
    !!     of Euler angle samples (for each dilation - which of course will
    !!     be an array of the same number (n_alpha*n_beta*n_gamma)).
    !!   - mode: Mode of thresholding to perform.  
    !!     If mode=CSWT_TR_THRES_NSIGMA_MODE_ABS then the absolute value of
    !!     coefficients are considered and those coefficients above (in
    !!     absolute value)  nsigma * sigma(a) are set to one, and the
    !!     remainder to zero.
    !!     If mode=CSWT_TR_THRES_NSIGMA_MODE_ABOVE then those coefficients 
    !!     above nsigma * sigma(a) are set to one, and the remainder to zero.
    !!     If mode=CSWT_TR_THRES_NSIGMA_MODE_BELOW then those coefficients 
    !!     below -nsigma * sigma(a) (*note minus sign*)  are set to one,
    !!     and the remainder to zero.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_thres_nsigma(tr, nsigma, n_effective_samples, mode)

      type(cswt_tr), intent(inout) :: tr
      real(s2_sp), intent(in) :: nsigma
      integer, intent(in) :: n_effective_samples(:)
      integer, intent(in) :: mode

      real(s2_sp), allocatable :: sigma(:), mean(:)
      integer :: i_dil, fail

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_xcorr_mask')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_xcorr_mask')
      end if

      ! Check n_effective samples_correct size.
      if(size(n_effective_samples) /= tr%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_wcoeff_sigma', &
           comment_add='Output arrays of invalid size')
      end if

      ! Allocate space.
      allocate(mean(tr%n_dilation), stat=fail)
      allocate(sigma(tr%n_dilation), stat=fail)
      if(fail /= 0) then
        call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
          'cswt_tr_xcorr_mask')
      end if

      ! Compute sigma for each dilation.
      call cswt_tr_wcoeff_sigma(tr, n_effective_samples, mean, sigma)

      ! Threshold wavelet coefficients for each dilation, using 
      ! thres = nsigma * sigma(i_dil).
      ! Note difference from mean considered, not simply value 
      ! (although generally the mean is zero (for xcorr)).
      do i_dil = 1,tr%n_dilation

         select case (mode)

           case (CSWT_TR_THRES_NSIGMA_MODE_ABS)
              where( abs(tr%wcoeff(i_dil,:,:,:) - mean(i_dil)) < nsigma * sigma(i_dil) )
                 tr%wcoeff(i_dil,:,:,:) = 0.0e0
              elsewhere
                 tr%wcoeff(i_dil,:,:,:) = 1.0e0
              end where

             case (CSWT_TR_THRES_NSIGMA_MODE_ABOVE)
                where( (tr%wcoeff(i_dil,:,:,:) - mean(i_dil)) < nsigma * sigma(i_dil) )
                   tr%wcoeff(i_dil,:,:,:) = 0.0e0
                elsewhere
                   tr%wcoeff(i_dil,:,:,:) = 1.0e0
                end where

           case (CSWT_TR_THRES_NSIGMA_MODE_BELOW)
              where( (tr%wcoeff(i_dil,:,:,:) - mean(i_dil)) > -nsigma * sigma(i_dil) )
                 tr%wcoeff(i_dil,:,:,:) = 0.0e0
              elsewhere
                 tr%wcoeff(i_dil,:,:,:) = 1.0e0
              end where

           case default
              call cswt_error(CSWT_ERROR_TR_THRESNSIG_MODE_INVALID, &
                'cswt_tr_mask_thres_nsigma')
         end select

      end do

      ! Free memory used.
      deallocate(mean, sigma)

    end subroutine cswt_tr_mask_thres_nsigma


    !--------------------------------------------------------------------------
    ! cswt_tr_mask_copy
    !
    !! Copy masks defined on the sky (Healpix format) for each dilation to 
    !! a cswt tr structure mask defined in the spherical wavelet coefficient 
    !! domain.  
    !!  
    !! Notes:
    !!   - The wavelet coefficients are not actually required in the 
    !!     algorithm, however the wavelet transform structure of the correct
    !!     size must be initialised.  To overcome this problem the tr_mask 
    !!     structure containing the wavelet coefficient masks is often used.
    !!
    !! Variables:
    !!   - tr_mask: Spherical wavelet transform of sky mask. Overwritten with
    !!     extended coefficient mask on output.  (Note that this strucutre 
    !!     doesn't actually have to have the wavelet coefficients of the
    !!     original sky map mask, however it is easier to just supply that.)
    !!   - sky_mask: Array of sky masks for each dilation that are to be
    !!     convert to the wavelet coefficient domain.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - December 2004
    !
    ! Revisions:
    !   December 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_copy(tr_mask, sky_mask)

      type(cswt_tr), intent(inout) :: tr_mask
      type(s2_sky), intent(in) :: sky_mask(:)
      
      real(s2_sp), allocatable :: sky_mask_ab(:,:)
      integer :: i_dil, i_gamma, fail

      ! Check object initialised.
      if(.not. tr_mask%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_copy')
      end if

      ! Check one sky mask for each dilation.
      if (size(sky_mask) /= tr_mask%n_dilation) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_mask_copy', &
           comment_add='Inconsistent number of dilations')
      end if

      ! Check mask sky initialised for each dilation
      do i_dil = 1,size(sky_mask) 
         if(.not. s2_sky_get_init(sky_mask(i_dil))) then
            call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_copy', &
                 comment_add='Mask sky not initialised')
         end if
      end do

      ! Allocate space for an ab representation of a sky mask.
      ! Overwritten each time through with current sky mask.
      allocate(sky_mask_ab(1:tr_mask%n_alpha, 1:tr_mask%n_beta), stat=fail)
      if(fail /= 0) then
        call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
          'cswt_tr_mask_copy')
      end if

      ! Considered in turn coefficient map corresponding to each dilation
      ! and orientation.
      do i_dil = 1,tr_mask%n_dilation
         do i_gamma = 0,tr_mask%n_gamma-1     ! Note indexed from 0.
            
            call s2_sky_extract_ab(sky_mask(i_dil), sky_mask_ab)           
            tr_mask%wcoeff(i_dil,:,:,i_gamma) = sky_mask_ab

         end do
      end do

      ! Free memory.
      deallocate(sky_mask_ab) 

    end subroutine cswt_tr_mask_copy


    !--------------------------------------------------------------------------
    ! cswt_tr_mask_gen
    !
    !! Generate an extended coefficient exclusion mask from an original sky
    !! mask and the spherical wavelet transform of the mask.
    !!  
    !! Notes:
    !!   - This routine may perform two different methods for generating the
    !!     extended mask: (a) from spherical wavelet coefficients or 
    !!     (b) by using morphological operations.  These are subsequently
    !!     described.
    !!   - (a) The extended coefficient exclusion mask is generated by taking
    !!     the transform of the original sky mask.  The only non-zero wavelet 
    !!     coefficients are those that are distorted by the mask boundary.  
    !!     These distorted coefficients are detected by thresholding.  An 
    !!     intermediate stage of smoothing (by convolving with a Gaussian
    !!     kernel) is performed to smooth out whispy regions that otherwise
    !!     remain when the wavelet support is exactly half contained in the 
    !!     masked region and half not.  A second stage of thresholding is 
    !!     performed to remove these regions once blurred.
    !!   - (b) The extended coefficient exclusion mask is generated by 
    !!     performing morphological operations on the original mask.  The 
    !!     wavelet coefficients are not requried, however the wavelet 
    !!     transform structure of the correct size must be initialised.  To
    !!     overcome this problem the tr_mask structure is actually used (since
    !!     it is already available) and simply overwritten with the new 
    !!     extended mask.  The method itself involves first removing the point
    !!     source mask regions by opening (morphological dialtion followed by
    !!     erosion) the original mask.  The glactic plane exclusion region
    !!     then remains.  This is extended simply by a morphological erosion
    !!     operation where the erosion kernel size is given by the 2*effective 
    !!     size of the wavelet on the sky.  A circular kernel shape is used.
    !!   - *Additional important note*:  Morphological mask generation is only
    !!     valid for the Mexhat, Morlet and Butterfly wavelets (or other 
    !!     wavelets that are defined to have the same effective size on the 
    !!     sky).  Modifications must be made for wavelets with a different
    !!     effective size on the sky.
    !!
    !! Variables:
    !!   - tr_mask: Spherical wavelet transform of sky mask. Overwritten with
    !!     extended coefficient mask on output.
    !!   - sky_mask: Original sky mask to be extended.
    !!   - mode: Integer to specifier what method of generating extended 
    !!     coefficient mask to use.
    !!   - [run_stage1]: Logical to specify whether to perform first stage 
    !!     of thresholding.
    !!   - [run_stage2]: Logical to specify whether to perform second stage 
    !!     of thresholding.
    !!   - [orig_mask_only_in]: Logical to specify if to just convert the 
    !!     original sky mask to the wavelet domain without extending it.
    !!   - [thres1_proportion_in]: Threshold proportion for the first stage
    !!     of dynamic thresholding.
    !!   - [thres2_proportion_in]: Threshold proportion for the second stage
    !!     of dynamic thresholding.
    !!   - [thres1_min_in]: Minimum threshold for the first stage of dynamic
    !!     thresholding.
    !!   - [std_in]: Standard deviation of Gaussian smoothing for each
    !!     dilation index.  If not specified then calculated by default method
    !!     from dilation values.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_mask_gen(tr_mask, sky_mask, mode, &
      run_stage1_in, run_stage2_in, orig_mask_only_in, &
      thres1_proportion_in, thres2_proportion_in, thres1_min_in, std_in)

      real(s2_dp), parameter :: &
        TWOPI = 6.283185307179586476925286766559005768394

      type(cswt_tr), intent(inout) :: tr_mask
      type(s2_sky), intent(in) :: sky_mask
      integer, intent(in) :: mode
      logical, intent(in), optional :: run_stage1_in, run_stage2_in
      logical, intent(in), optional :: orig_mask_only_in
      real(s2_sp), intent(in), optional :: thres1_proportion_in
      real(s2_sp), intent(in), optional :: thres2_proportion_in
      real(s2_sp), intent(in), optional :: thres1_min_in
      real(s2_sp), intent(in), optional :: std_in(:)

      logical :: run_stage1, run_stage2, orig_mask_only
      logical, parameter :: DEFAULT_RUN_STAGE1 = .true.
      logical, parameter :: DEFAULT_RUN_STAGE2 = .true.
      logical, parameter :: DEFAULT_ORIG_MASK_ONLY = .false.
      integer :: i_dil, i_gamma, fail
      integer, parameter :: kernel_size1 = 11
      integer, parameter :: kernel_size2 = 11
      real(s2_sp) :: DILATION_TO_STD_FACTOR = 90.0e0
      real(s2_sp), parameter :: DEFAULT_THRES1_MIN = 1.0e-3
      real(s2_sp) :: thres1_min
      real(s2_sp) :: thres1_proportion, thres2_proportion
      real(s2_sp), parameter :: DEFAULT_THRES1_PROPORTION = 0.1e0
      real(s2_sp), parameter :: DEFAULT_THRES2_PROPORTION = 0.9e0
      real(s2_sp) :: thres1, thres2, min_dilation
      real(s2_sp), allocatable :: sky_mask_ab(:,:)
      real(s2_sp), allocatable :: std(:)
      integer, parameter :: KERNEL_SIZE_POINT_REMOVAL = 5
      real(s2_sp) :: max_dilation, erode_scale_angle
      integer :: erode_scale_pixel

      ! Set optional input threshold parameters.
      if(present(run_stage1_in)) then
         run_stage1 = run_stage1_in
      else
         run_stage1 = DEFAULT_RUN_STAGE1
      end if
      if(present(run_stage2_in)) then
         run_stage2 = run_stage2_in
      else
         run_stage2 = DEFAULT_RUN_STAGE2
      end if
      if(present(orig_mask_only_in)) then
         orig_mask_only = orig_mask_only_in
      else
         orig_mask_only = DEFAULT_ORIG_MASK_ONLY
      end if

      if(present(thres1_min_in)) then
         thres1_min = thres1_min_in
      else
         thres1_min = DEFAULT_THRES1_MIN
      end if
      if(present(thres1_proportion_in)) then
         thres1_proportion = thres1_proportion_in
      else
         thres1_proportion = DEFAULT_THRES1_PROPORTION
      end if
      if(present(thres2_proportion_in)) then
         thres2_proportion = thres2_proportion_in
      else
         thres2_proportion = DEFAULT_THRES2_PROPORTION
      end if

      ! Check object initialised.
      if(.not. tr_mask%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_gen')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr_mask%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_mask_gen')
      end if

      ! Check mask sky initialised.
      if(.not. s2_sky_get_init(sky_mask)) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_mask_gen', &
          comment_add='Mask sky not initialised')
      end if

      ! Set up std for smoothing.
      allocate(std(1:tr_mask%n_dilation), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_mask_gen')
      end if
      if(present(std_in)) then
         if(size(std_in) /= tr_mask%n_dilation) then
            call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_tr_mask_gen', &
              comment_add='Input smoothing std invalid size')
         end if
         std = std_in
      else
         do i_dil = 1,tr_mask%n_dilation
            min_dilation &
              = min(tr_mask%dilation(i_dil,1),  tr_mask%dilation(i_dil,2))
            std(i_dil) = DILATION_TO_STD_FACTOR * min_dilation
         end do
      end if

      ! Get alpha-beta (ecp) map of original mask.
      ! This is same for all masks generated from mask coefficients 
      ! so just get once here before use for all.
      allocate(sky_mask_ab(1:tr_mask%n_alpha, 1:tr_mask%n_beta), stat=fail)
      if(fail /= 0) then
        call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
          'cswt_tr_mask_gen')
      end if
      call s2_sky_extract_ab(sky_mask, sky_mask_ab)
      
      ! Considered in turn coefficient map corresponding to each dilation
      ! and orientation.
      do i_dil = 1,tr_mask%n_dilation
         do i_gamma = 0,tr_mask%n_gamma-1     ! Note indexed from 0.

            ! Generate extended mask by taking transform, thresholding and
            ! smoothing (and thresholding again).
            if(mode == CSWT_TR_MASK_GEN_CSWT) then

               if(run_stage1) then

                  ! Perfom initial thresholding of mask coefficients.
                  ! (If above zero (thres1) then distorted and should become 
                  ! ones of final mask.)
                  thres1 = thres1_proportion &
                       * maxval(abs(tr_mask%wcoeff(i_dil,:,:,i_gamma)))
                  if(thres1 < thres1_min) then 
                     tr_mask%wcoeff(i_dil,:,:,i_gamma) = 1.0e0
                  else
                     where(abs(tr_mask%wcoeff(i_dil,:,:,i_gamma)) < thres1)
                        tr_mask%wcoeff(i_dil,:,:,i_gamma) = 1.0e0
                     elsewhere
                        tr_mask%wcoeff(i_dil,:,:,i_gamma) = 0.0e0
                     end where
                  end if
               end if
               
               if(run_stage2) then
                  
                  ! Smooth thresholded mask transform.
                  ! Whispy regions that remain after first stage of
                  ! thresholding get blurred out.
                  call cswt_tr_smooth_gaussian_ab(tr_mask, &
                       kernel_size1, kernel_size2, &
                       std(i_dil), &
                       i_dil, i_gamma+1)    ! i_gamma plus 1 so indexed from 1.
                  
                  ! Perform second stage thresholding.
                  ! If below thres then set to zero (remove whispy regions
                  ! that remain before smoothing).
                  thres2 = thres2_proportion &
                       * maxval(abs(tr_mask%wcoeff(i_dil,:,:,i_gamma)))
                  where(abs(tr_mask%wcoeff(i_dil,:,:,i_gamma)) > thres2)
                     tr_mask%wcoeff(i_dil,:,:,i_gamma) = 1.0e0
                  elsewhere
                     tr_mask%wcoeff(i_dil,:,:,i_gamma) = 0.0e0
                  end where
                  
               end if
               
               
               if(orig_mask_only) then
                  
                  tr_mask%wcoeff(i_dil,:,:,i_gamma) = sky_mask_ab
                  
               else
                  
                  ! Now have extended regions mask, but need to also include
                  ! original mask.
                  
                  ! Multiply two masks in alpha-beta format to generate final
                  ! extended coefficient mask.
                  tr_mask%wcoeff(i_dil,:,:,i_gamma) &
                       = tr_mask%wcoeff(i_dil,:,:,i_gamma) * sky_mask_ab
                  
               end if

            ! Generate mask by morphological operations.    
            elseif(mode == CSWT_TR_MASK_GEN_MORPH) then

               tr_mask%wcoeff(i_dil,:,:,i_gamma) = sky_mask_ab

               ! Remove point source regions by closing map.
               call cswt_tr_morph_dilate_ab(tr_mask, i_dil, i_gamma+1, &
                 KERNEL_SIZE_POINT_REMOVAL, CSWT_TR_KERNEL_TYPE_CIRC)
               call cswt_tr_morph_erode_ab(tr_mask, i_dil, i_gamma+1, &
                 KERNEL_SIZE_POINT_REMOVAL, CSWT_TR_KERNEL_TYPE_CIRC)

               ! Determine kernel size for morphological erosion (mask
               ! extension).

               max_dilation = max(tr_mask%dilation(i_dil,1),  &
                                  tr_mask%dilation(i_dil,2))

               ! Set erode angle to effective size of wavelet on sky.
! Change here to alter size of entended exclusion region.

! Used to generate extT masks.
               erode_scale_angle = 1.0e0 * &
                 4.0e0 * atan(max_dilation / sqrt(2.0e0))

! Use to generate extT1.5 masks.
!               erode_scale_angle = 1.5e0 * &
!                 4.0e0 * atan(max_dilation / sqrt(2.0e0))

! Use to generate extT2 masks.
!               erode_scale_angle = 2.0e0 * &
!                 4.0e0 * atan(max_dilation / sqrt(2.0e0))

               ! Determine no. of pixels equivalent to angle.
               erode_scale_pixel = erode_scale_angle * tr_mask%n_alpha/TWOPI

               ! Ensure erode_scale_pixel (=kernel_size) is odd.
! Kernel size no longer has to be odd.
!               if(mod(erode_scale_pixel,2) == 0) then
!                  erode_scale_pixel = erode_scale_pixel + 1 !+1=conservative
!               end if

               ! Perform morphological erosion to extend mask.
               call cswt_tr_morph_erode_ab(tr_mask, i_dil, i_gamma+1, &
                 erode_scale_pixel, CSWT_TR_KERNEL_TYPE_CIRC)
               
               ! Multiply by original mask on alpha-beta format.
               tr_mask%wcoeff(i_dil,:,:,i_gamma) &
                    = tr_mask%wcoeff(i_dil,:,:,i_gamma) * sky_mask_ab

            else

               call cswt_error(CSWT_ERROR_TR_MASK_MODE_INVALID, &
                 'cswt_tr_mask_gen')

            end if

         end do
      end do

      ! Free memory.
      deallocate(sky_mask_ab)
      deallocate(std)

    end subroutine cswt_tr_mask_gen


    !--------------------------------------------------------------------------
    ! cswt_tr_morph_erode_ab
    !
    !! Morphologically erode an ecp alpha-beta map with either a square or
    !! circular kernel.
    !!  
    !! Notes:
    !!   - The input tr wavelet coefficient alpha-beta map must be a binary
    !!     map (the ones in the map are eroded, i.e. zero component dilated).
    !!   - Kernel size variable is positive size (so total size of kernel
    !!     is 2*kernel_size + 1).
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!   - Overwrites the tr alpha and beta components with the 
    !!     morphologically dilated values.
    !!
    !! Variables:
    !!   - tr: tr structure containing wavelet coefficient that are 
    !!     morphologically eroded.
    !!   - i_dil: Dilation index of tr%wcoeff map to erode.
    !!   - i_gamma: Gamma Euler angle index of tr%wcoeff to erode.
    !!   - kernel_size: Size of 2D kernel in alpha and beta Euler angle
    !!     dimension (must be odd).
    !!   - kernel_type: Specifier to indicate shape of kernel to apply.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_morph_erode_ab(tr, i_dil, i_gamma, kernel_size, &
      kernel_type)

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: i_dil, i_gamma
      integer, intent(in) :: kernel_size
      integer, intent(in) :: kernel_type

      real(s2_sp), allocatable :: x_ab(:,:), x_ab_morph_erode(:,:)
      real(s2_sp), allocatable :: kernel(:,:)
      integer :: nk1_on2, nk2_on2, ia, ib, ib_use, ik1, ik2, fail
      real(s2_sp) :: cur_val, dist_from_cen_sqd

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_morph_erode_ab')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_morph_erode_ab')
      end if
      
      ! Check i_dil and i_gamma are in valid range.
      if(i_dil < 1 .or. i_dil > tr%n_dilation &
        .or. i_gamma < 1 .or. i_gamma > tr%n_gamma) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_morph_erode_ab', &
           comment_add='Dilation or orientation index invalid')
      end if

      ! Check kernel valid size.
!      if(mod(kernel_size,2) /= 1) then
!         call cswt_error(CSWT_ERROR_TR_KERNEL_INVALID, &
!           'cswt_tr_morph_erode_ab', &
!           comment_add='Kernel size not odd')
!      end if

      ! Allocate space for temporary storage.
      allocate(x_ab(1:tr%n_alpha, 1:tr%n_beta), stat=fail)
      allocate(x_ab_morph_erode(1:tr%n_alpha, 1:tr%n_beta), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_morph_erode_ab')
      end if

      ! Extract out copy of alpha beta wavelet coefficients considered 
      ! so have consistend indices.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      x_ab = tr%wcoeff(i_dil,:,:,i_gamma-1)

      ! Initialise morphologically eroded map with ones.
      x_ab_morph_erode = 1.0e0
      
      ! Set half kernel sizes. 
      nk1_on2 = kernel_size
      nk2_on2 = kernel_size
!     nk1_on2 = kernel_size
!     nk2_on2 = (kernel_size - 1) / 2

      ! Allocate space for kernel and compute.
      allocate(kernel(-nk1_on2:nk1_on2, -nk2_on2:nk2_on2), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_morph_erode_ab')
      end if
      if(kernel_type == CSWT_TR_KERNEL_TYPE_CIRC) then
         ! Calculate circular kernel.
         do ik1 = -nk1_on2,nk1_on2
            do ik2 = -nk2_on2,nk2_on2          
               dist_from_cen_sqd = ik1**2 + ik2**2
               if(dist_from_cen_sqd <= nk1_on2**2) then! Use nk1_on2 since 
                                                       ! kernel size is square.
                  kernel(ik1,ik2) = 1.0e0
               else
                  kernel(ik1,ik2) = 0.0e0
               end if
            end do
         end do
      elseif(kernel_type == CSWT_TR_KERNEL_TYPE_SQR) then
         ! Calculate square kernel, i.e. all values 1.0e0.
         kernel = 1.0e0         
      else
         call cswt_error(CSWT_ERROR_TR_KERNEL_INVALID, &
           'cswt_tr_morph_erode_ab', &
           comment_add='Kernel type invalid')
      end if

      ! Perform morphological erosion.
      do ia = 1,tr%n_alpha
        do ib = 1,tr%n_beta
    
          do ik1 = -nk1_on2,nk1_on2
            do ik2 = -nk2_on2,nk2_on2

              ! Actually reflect beta.
              if(ib + ik2 < 1) then
                ib_use = abs(ib + ik2) + 1
              elseif(ib + ik2 == tr%n_beta + 1) then
                ib_use = tr%n_beta
              elseif(ib + ik2 > tr%n_beta + 1) then
                ib_use = 2*tr%n_beta - (ib + ik2) + 2
              else
                ib_use = ib + ik2
              end if

              ! Wrap alpha (phi) around sphere.
              if(ia + ik1 < 1) then
                cur_val = x_ab(ia + ik1 + tr%n_alpha, ib_use)
              elseif(ia + ik1 > tr%n_alpha) then
                cur_val = x_ab(ia + ik1 - tr%n_alpha, ib_use)               
              else
                cur_val = x_ab(ia + ik1, ib_use)
              end if
                  
              if(kernel(ik1,ik2) == 1.0e0) then
                 x_ab_morph_erode(ia,ib) = cur_val * x_ab_morph_erode(ia,ib)
              end if
              ! Else if kernel element is zero then just ignore.  

            end do
          end do

        end do
      end do

      ! Write smoothed values back to wavelet coefficients.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      tr%wcoeff(i_dil,:,:,i_gamma-1) = x_ab_morph_erode

      ! Free memory used.
      deallocate(x_ab)
      deallocate(x_ab_morph_erode)
      deallocate(kernel)

    end subroutine cswt_tr_morph_erode_ab


    !--------------------------------------------------------------------------
    ! cswt_tr_morph_dilate_ab
    !
    !! Morphologically dilate an ecp alpha-beta map with either a square or
    !! circular kernel.
    !!  
    !! Notes:
    !!   - The input tr wavelet coefficient alpha-beta map must be a binary
    !!     map (the ones in the map are dilated).
    !!   - Kernel size variable is positive size (so total size of kernel
    !!     is 2*kernel_size + 1).
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!   - Overwrites the tr alpha and beta components with the 
    !!     morphologically dilated values.
    !!
    !! Variables:
    !!   - tr: tr structure containing wavelet coefficient that are 
    !!     morphologically dilated.
    !!   - i_dil: Dilation index of tr%wcoeff map to dilate.
    !!   - i_gamma: Gamma Euler angle index of tr%wcoeff to dilate.
    !!   - kernel_size: Size of 2D kernel in alpha and beta Euler angle
    !!     dimension (must be odd).
    !!   - kernel_type: Specifier to indicate shape of kernel to apply.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_morph_dilate_ab(tr, i_dil, i_gamma, kernel_size, &
      kernel_type)

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: i_dil, i_gamma
      integer, intent(in) :: kernel_size
      integer, intent(in) :: kernel_type

      real(s2_sp), allocatable :: x_ab(:,:), x_ab_morph_dilate(:,:)
      real(s2_sp), allocatable :: kernel(:,:)
      integer :: nk1_on2, nk2_on2, ia, ib, ib_use, ik1, ik2, fail
      real(s2_sp) :: cur_val, dist_from_cen_sqd
      
      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_morph_dilate_ab')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_morph_dilate_ab')
      end if
      
      ! Check i_dil and i_gamma are in valid range.
      if(i_dil < 1 .or. i_dil > tr%n_dilation &
        .or. i_gamma < 1 .or. i_gamma > tr%n_gamma) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_morph_dilate_ab', &
           comment_add='Dilation or orientation index invalid')
      end if

      ! Check kernel valid size.
!      if(mod(kernel_size,2) /= 1) then
!         call cswt_error(CSWT_ERROR_TR_KERNEL_INVALID, &
!           'cswt_tr_morph_dilate_ab', &
!           comment_add='Kernel size not odd')
!      end if

      ! Allocate space for temporary storage.
      allocate(x_ab(1:tr%n_alpha, 1:tr%n_beta), stat=fail)
      allocate(x_ab_morph_dilate(1:tr%n_alpha, 1:tr%n_beta), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_morph_dilate_ab')
      end if

      ! Extract out copy of alpha beta wavelet coefficients considered 
      ! so have consistend indices.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      x_ab = tr%wcoeff(i_dil,:,:,i_gamma-1)

      ! Initialise morphologically dilated map with zeros.
      x_ab_morph_dilate = 0.0e0
      
      ! Set half kernel sizes. 
      ! more samples in alpha direction Ensure kernel size is square.
      nk1_on2 = kernel_size
      nk2_on2 = kernel_size
!      nk2_on2 = (kernel_size - 1) / 2
!      nk2_on2 = (kernel_size - 1) / 2

      ! Allocate space for kernel and compute.
      allocate(kernel(-nk1_on2:nk1_on2, -nk2_on2:nk2_on2), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_morph_dilate_ab')
      end if
      if(kernel_type == CSWT_TR_KERNEL_TYPE_CIRC) then
         ! Calculate circular kernel.
         do ik1 = -nk1_on2,nk1_on2
            do ik2 = -nk2_on2,nk2_on2          
               dist_from_cen_sqd = ik1**2 + ik2**2
               if(dist_from_cen_sqd <= nk1_on2**2) then! Use nk1_on2 since 
                                                       ! kernel size is square.
                  kernel(ik1,ik2) = 1.0e0
               else
                  kernel(ik1,ik2) = 0.0e0
               end if
            end do
         end do
      elseif(kernel_type == CSWT_TR_KERNEL_TYPE_SQR) then
         ! Calculate square kernel, i.e. all values 1.0e0.
         kernel = 1.0e0         
      else
         call cswt_error(CSWT_ERROR_TR_KERNEL_INVALID, &
           'cswt_tr_morph_dilate_ab', &
           comment_add='Kernel type invalid')
      end if

      ! Perform morphological dilation.
      do ia = 1,tr%n_alpha
        do ib = 1,tr%n_beta
    
          do ik1 = -nk1_on2,nk1_on2
            do ik2 = -nk2_on2,nk2_on2

              ! Actually reflect beta.
              if(ib + ik2 < 1) then
                ib_use = abs(ib + ik2) + 1
              elseif(ib + ik2 == tr%n_beta + 1) then
                ib_use = tr%n_beta
              elseif(ib + ik2 > tr%n_beta + 1) then
                ib_use = 2*tr%n_beta - (ib + ik2) + 2
              else
                ib_use = ib + ik2
              end if

              ! Wrap alpha (phi) around sphere.
              if(ia + ik1 < 1) then
                cur_val = x_ab(ia + ik1 + tr%n_alpha, ib_use)
              elseif(ia + ik1 > tr%n_alpha) then
                cur_val = x_ab(ia + ik1 - tr%n_alpha, ib_use)               
              else
                cur_val = x_ab(ia + ik1, ib_use)
              end if
              
              if(kernel(ik1,ik2) == 1.0e0) then
                x_ab_morph_dilate(ia,ib) = cur_val + x_ab_morph_dilate(ia,ib)
                x_ab_morph_dilate(ia,ib) = min(x_ab_morph_dilate(ia,ib), 1.0e0)
              end if
              ! Else if kernel element is zero then just ignore.  
                
            end do
          end do

        end do
      end do

      ! Write smoothed values back to wavelet coefficients.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      tr%wcoeff(i_dil,:,:,i_gamma-1) = x_ab_morph_dilate

      ! Free memory used.
      deallocate(x_ab)
      deallocate(x_ab_morph_dilate)
      deallocate(kernel)

    end subroutine cswt_tr_morph_dilate_ab


    !--------------------------------------------------------------------------
    ! cswt_tr_smooth_gaussian
    !
    !! Smooth all tr%wcoeff maps over alpha-beta for each dilation and gamma
    !! Euler orientation.
    !!  
    !! Notes:
    !!   - Overwrites the tr wavelet coefficients with the smoothed values.
    !!
    !! Variables:
    !!   - tr: tr structure containing wavelet coefficient that are smoothed.
    !!   - kernel_size1: Size of 2D kernel in alpha Euler angle dimension 
    !!     (must be odd).
    !!   - kernel_size2: Size of 2D kernel in beta Euler angle dimension
    !!     (must be odd).
    !!   - std: Standard deviation of Gaussian used to smooth.  Array indexed
    !!     over i_dil and i_gamma containing the standard deviation for each
    !!     alpha-beta map.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_smooth_gaussian(tr, kernel_size1, kernel_size2, std)

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: kernel_size1, kernel_size2
      real(s2_sp), intent(in) :: std(:,:)

      integer :: i_dil, i_gamma

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_smooth_gaussian')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_smooth_gaussian')
      end if

      ! Check std size valid.
      if(size(std,1) /= tr%n_dilation .or. size(std,2) /= tr%n_gamma) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_smooth_gaussian', &
           comment_add='Standard deviation size invalid')
      end if

      ! Perform smoothing for each dilation and orientation. 
      do i_dil = 1,tr%n_dilation
         do i_gamma = 1,tr%n_gamma    ! Note indexed from 1.

            call cswt_tr_smooth_gaussian_ab(tr, kernel_size1, kernel_size2, &
              std(i_dil, i_gamma), i_dil, i_gamma) 

         end do
      end do

    end subroutine cswt_tr_smooth_gaussian


    !--------------------------------------------------------------------------
    ! cswt_tr_smooth_gaussian_ab
    !
    !! Smooth an ecp alpha-beta map by convolving it with a Gaussian kernel.
    !!  
    !! Notes:
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!   - Overwrites the tr alpha and beta components with the smoothed 
    !!     values.
    !!
    !! Variables:
    !!   - tr: tr structure containing wavelet coefficient that are smoothed.
    !!   - kernel_size1: Size of 2D kernel in alpha Euler angle dimension 
    !!     (must be odd).
    !!   - kernel_size2: Size of 2D kernel in beta Euler angle dimension
    !!     (must be odd).
    !!   - std: Standard deviation of Gaussian used to smooth.
    !!   - i_dil: Dilation index of tr%wcoeff map to smooth.
    !!   - i_gamma: Gamma Euler angle index of tr%wcoeff to smooth.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_smooth_gaussian_ab(tr, kernel_size1, kernel_size2, &
      std, i_dil, i_gamma) 

      type(cswt_tr), intent(inout) :: tr
      integer, intent(in) :: kernel_size1, kernel_size2
      real(s2_sp), intent(in) :: std
      integer, intent(in) :: i_dil, i_gamma

      real(s2_sp), allocatable :: x_ab(:,:), x_ab_smooth(:,:)
      real(s2_sp), allocatable :: kernel(:,:)
      real(s2_sp) :: normalisation, cur_val
      integer :: nk1_on2, nk2_on2, ik1, ik2, ia, ib, ib_use, fail

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_smooth_gaussian_ab')
      end if

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_smooth_gaussian_ab')
      end if
      
      ! Check i_dil and i_gamma are in valid range.
      if(i_dil < 1 .or. i_dil > tr%n_dilation &
        .or. i_gamma < 1 .or. i_gamma > tr%n_gamma) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_smooth_gaussian_ab', &
           comment_add='Dilation or orientation index invalid')
      end if

      ! Check kernel valid.
      if(mod(kernel_size1, 2) /= 1 .or.mod(kernel_size2, 2) /= 1) then
         call cswt_error(CSWT_ERROR_TR_KERNEL_INVALID, &
           'cswt_tr_smooth_gaussian_ab', &
           comment_add='Kernel size not odd')
      end if
      if(kernel_size1 < 2 .or. kernel_size2 < 2) then
         call cswt_error(CSWT_ERROR_TR_KERNEL_INVALID, &
           'cswt_tr_smooth_gaussian_ab', &
           comment_add='Kernel size less than 2')
      end if

      ! Allocate space for temporary storage.
      allocate(x_ab(1:tr%n_alpha, 1:tr%n_beta), stat=fail)
      allocate(x_ab_smooth(1:tr%n_alpha, 1:tr%n_beta), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_smooth_gaussian_ab')
      end if

      ! Extract out copy of alpha beta wavelet coefficients considered 
      ! so have consistend indices.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      x_ab = tr%wcoeff(i_dil,:,:,i_gamma-1)

      ! Initialised smoothed values to zero.
      x_ab_smooth = 0.0e0

      ! Set half kernel sizes.
      nk1_on2 = (kernel_size1 - 1) / 2
      nk2_on2 = (kernel_size2 - 1) / 2
      
      ! Allocate space for kernel.
      allocate(kernel(-nk1_on2:nk1_on2, -nk2_on2:nk2_on2), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_smooth_gaussian_ab')
      end if

      ! Calculate Gaussian kernel.
      normalisation = 1 / sqrt(2.0e0 * pi * std**2)
      do ik1 = -nk1_on2,nk1_on2
        do ik2 = -nk2_on2,nk2_on2          
          kernel(ik1,ik2) = exp(-(ik1**2 + ik1**2)/ (2.0e0 * std**2) )
        end do
      end do
      kernel = kernel * normalisation

      ! Perform convolution.
      do ia = 1,tr%n_alpha
        do ib = 1,tr%n_beta
    
          do ik1 = -nk1_on2,nk1_on2
            do ik2 = -nk2_on2,nk2_on2
    
              ! Crop beta range at edges.
!              if(ib + ik2 < 1) then
!                ib_use = 1
!              elseif(ib + ik2 > tr%n_beta) then
!                ib_use = tr%n_beta
!              else
!                ib_use = ib + ik2
!              end if

              ! Actually reflect beta.
              if(ib + ik2 < 1) then
                ib_use = abs(ib + ik2) + 1
              elseif(ib + ik2 == tr%n_beta + 1) then
                ib_use = tr%n_beta
              elseif(ib + ik2 > tr%n_beta + 1) then
                ib_use = 2*tr%n_beta - (ib + ik2) + 2
              else
                ib_use = ib + ik2
              end if

              ! Wrap alpha (phi) around sphere.
              if(ia + ik1 < 1) then
                cur_val = x_ab(ia + ik1 + tr%n_alpha, ib_use)
              elseif(ia + ik1 > tr%n_alpha) then
                cur_val = x_ab(ia + ik1 - tr%n_alpha, ib_use)               
              else
                cur_val = x_ab(ia + ik1, ib_use)
              end if
              
              x_ab_smooth(ia, ib) = x_ab_smooth(ia, ib) + &
                & kernel(ik1, ik2) * cur_val
                
            end do
          end do
          
        end do
      end do
    
      ! Write smoothed values back to wavelet coefficients.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      tr%wcoeff(i_dil,:,:,i_gamma-1) = x_ab_smooth

      ! Free memory used.
      deallocate(kernel)
      deallocate(x_ab)
      deallocate(x_ab_smooth)
    
    end subroutine cswt_tr_smooth_gaussian_ab


    !--------------------------------------------------------------------------
    ! Extract and write skies (not just ecp Euler angle samples)
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tr_extract_wcoeff_sky
    !
    !! Extract a sky representation of the ecp (equi-sampled) spherical 
    !! wavelet coefficients for a given dilation and orientation.
    !!  
    !! Notes:
    !!   - Initialises a new sky form the wavelet coefficients of the tr
    !!     structure.
    !!   - The returned sky is subsequently independed of the tr structure
    !!     and should be freed by the calling routine at some point.
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!
    !! Variables:
    !!   - tr: The tr structure containing the computed spherical wavelet 
    !!     coefficients that are to be written to a sky representation.
    !!   - i_dil: The dilation index of the spherical wavelet coefficients
    !!     to write the sky representation of.
    !!   - i_gamma: The orientation index of the spherical wavelet 
    !!     coefficients to write the sky representation of.
    !!   - nside: The Healpix nside resolution parameter of the constructed 
    !!     sky.  
    !!   - interp: Logical to specify whether linear interpolation is to be 
    !!     performed.
    !!   - [pix_scheme_in]:  Pixelisation scheme of initialised sky.  If not 
    !!     specified then the default pixelisation scheme is used.
    !!   - sky: Sky representation of the wavelet coefficients for the 
    !!     specified dilation and orientation.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_extract_wcoeff_sky(tr, i_dil, i_gamma, nside, interp, &
      pix_scheme_in) result(sky)

      type(cswt_tr), intent(in) :: tr
      integer, intent(in) :: i_dil, i_gamma, nside
      logical, intent(in) :: interp
      integer, intent(in), optional :: pix_scheme_in
      type(s2_sky) :: sky

      integer :: pix_scheme
      integer, parameter :: DEFAULT_PIX_SCHEME = S2_SKY_RING

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_extract_wcoeff_sky')
      end if 

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
           'cswt_tr_extract_wcoeff_sky')
      end if
      
      ! Check i_dil and i_gamma are in valid range.
      if(i_dil < 1 .or. i_dil > tr%n_dilation &
        .or. i_gamma < 1 .or. i_gamma > tr%n_gamma) then
         call cswt_error(CSWT_ERROR_SIZE_INVALID, &
           'cswt_tr_extract_wcoeff_sky', &
           comment_add='Dilation or orientation index invalid')
      end if

      ! Set pixelisation scheme.
      if(present(pix_scheme_in)) then
         pix_scheme = pix_scheme_in
      else
         pix_scheme = DEFAULT_PIX_SCHEME
      end if

      ! Initialise sky from alpha-beta array.
      ! Note -1 since i_gamma in range 1:N but data stored from 0:N-1.
      sky = s2_sky_init(tr%wcoeff(i_dil,:,:,i_gamma-1), interp, &
        nside, pix_scheme)

    end function cswt_tr_extract_wcoeff_sky


    !--------------------------------------------------------------------------
    ! cswt_tr_io_fits_write_wcoeff_sky
    !
    !! Extract a sky representation of the ecp (equi-sampled) spherical
    !! wavelet coefficients for a given dilation and orientation and write 
    !! them directly to a fits sky file.
    !!
    !! Notes:
    !!   - Both the dilation (i_dil) and *orientation* (i_gamma) are
    !!     indexed from 1:N.
    !!
    !! Variables:
    !!   - filename: Name of the output fits file to write the sky 
    !!     representation of the wavelet coefficients to.
    !!   - tr: The tr structure containing the computed spherical wavelet 
    !!     coefficients that are to be written fits file sky representation.
    !!   - i_dil: The dilation index of the spherical wavelet coefficients
    !!     to write the sky representation of.
    !!   - i_gamma: The orientation index of the spherical wavelet 
    !!     coefficients to write the sky representation of.
    !!   - nside: The Healpix nside resolution parameter of the constructed 
    !!     sky.  
    !!   - interp: Logical to specify whether linear interpolation is to be 
    !!     performed.
    !!   - [pix_scheme_in]:  Pixelisation scheme of initialised sky.  If not 
    !!     specified then the default pixelisation scheme is used.
    !!   - [comment]: Optional comment to add to the header of the output 
    !!     fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_fits_write_wcoeff_sky(filename, tr, i_dil, i_gamma, &
      nside, interp, pix_scheme_in, comment)

      character(len=*), intent(in) :: filename
      type(cswt_tr), intent(in) :: tr
      integer, intent(in) :: i_dil, i_gamma, nside
      logical, intent(in) :: interp
      integer, intent(in), optional :: pix_scheme_in
      character(len=*), intent(in), optional :: comment

      type(s2_sky) :: sky

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, &
          'cswt_tr_io_fits_write_wcoeff_sky')
      end if 

      ! Extract sky from tr wavelet coefficients.
      sky = cswt_tr_extract_wcoeff_sky(tr, i_dil, i_gamma, nside, interp, &
        pix_scheme_in) 

      ! Write to file.
      call s2_sky_write_map_file(sky, filename, comment)

      ! Free temporary sky extracted.
      call s2_sky_free(sky)

    end subroutine cswt_tr_io_fits_write_wcoeff_sky


    !--------------------------------------------------------------------------
    ! File IO routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tr_io_fits_write_wcoeff
    !
    !! Write a tr structure to an output fits file.  Both the wavelet 
    !! coefficients and the dilation values are written (and other additional
    !! wavelet transform parameters).
    !!
    !! Variables:
    !!   - filename: Name of the output fits file to write the tr structure
    !!     data to.
    !!   - tr: The tr structure containing the data to be written to the
    !!     output fits file.
    !!   - [comment]: Optional comment string to be added to the output fits 
    !!     file header if present.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_fits_write_wcoeff(filename, tr, comment)

      character(len=*), intent(in) :: filename
      type(cswt_tr), intent(in) :: tr
      character(len=*), intent(in), optional :: comment

      integer, parameter :: NUM_WCOEFF_EULER_DIM = 3

      integer :: status,unit,blocksize,bitpix
      integer :: group,dim1,dim2
      logical :: simple, extend, file_exists
      integer :: decimals
      integer :: naxis
      integer :: naxes(1), naxes_wcoeff(NUM_WCOEFF_EULER_DIM)
      integer :: tfields, nrows, varidat
      character(len=32) :: ttype(2), tform(2), tunit(2), extname
      integer :: frow, felem, colnum
      integer :: i_dil

      type(s2_sky) :: temp_sky
      logical :: norm_preserve
      real(s2_sp), allocatable :: admiss(:)
      integer :: fail
      integer :: n_param, i_param
      real(s2_sp), allocatable :: param(:)
      character(len=S2_STRING_LEN) :: param_str

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_io_fits_write_wcoeff')
      end if 

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
        call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, &
          'cswt_tr_io_fits_write_wcoeff')
      end if

      ! Define FITS parameters.

      bitpix=-32 ! Real single precision.
      status=0   ! Initialse error status to zero.

      decimals = 4 ! Number of decimals present in dilation keyword.

      ! Check if file already exists.
      call cswt_tr_io_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call cswt_error(CSWT_ERROR_TR_FILE_EXISTS, &
              'cswt_tr_io_fits_write_wcoeff')
        stop
      end if

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Create the new empty fits file.
      blocksize=1  ! Historical artifact that is ignored.
      call ftinit(unit,filename,blocksize,status)

      ! Write primary header.
      simple=.true.
      extend=.true.
      naxis=0
      naxes(1)=0
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      ! Write additional header keywords.
      call ftpcom(unit, &
        '  Spherical wavelet coefficients created by cswt-0.2.',status)
      call ftpcom(unit, &
        '  Primary extension empty.',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if
      if(tr%swav_mother_status) then
         call ftpkys(unit, 'WAVELET', cswt_swav_get_name(tr%swav_mother), &
              'wavelet type',status)
      end if
      call ftpkyj(unit,'NDIL', tr%n_dilation, &
        'length of dilation array',status)
      call ftpkyj(unit,'NALPHA',tr%n_alpha, &
        'length of alpha Euler angle dimension',status)
      call ftpkyj(unit,'NBETA',tr%n_beta, &
        'length of beta Euler angle dimension',status)
      call ftpkyj(unit,'NGAMMA',tr%n_gamma, &
        'length of gamma Euler angle dimension',status)
      call ftpkyj(unit,'LMAX',tr%lmax, &
        'max spherical harmonic l considered',status)
      call ftpkyj(unit,'MMAX',tr%mmax, &
        'max spherical harmonic m considered',status)

      ! Write wavelet parameters to header if present.
      if(tr%swav_mother_status) then
         temp_sky = cswt_swav_get_sky(tr%swav_mother)
         n_param = s2_sky_get_n_param(temp_sky)
         call ftpkyj(unit,'NPARAM', n_param, &
              'number of wavelet parameters stored',status)
         if(n_param /= 0) then
            allocate(param(1:n_param), stat=fail)
            if(fail /= 0) then
               call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
                    'cswt_tr_io_fits_write_wcoeff')
            end if
            call s2_sky_get_param(temp_sky, param)
            do i_param = 1, n_param
               write(param_str,'(a,i1)') 'PARAM', i_param
               call ftpkye(unit,trim(param_str),param(i_param),decimals, &
                    'parameter value',status)
            end do
            deallocate(param)
         end if
         call s2_sky_free(temp_sky)
      end if

      ! Insert binary table extension for dilations.
      extname='DIL'
      ttype(1)='DIL1'
      ttype(2)='DIL2'
      tform(1)='1E'
      tform(2)='1E'
      tunit(1)='direct (rad)'
      tunit(2)='direct (rad)'
      tfields=2
      nrows=tr%n_dilation
      varidat=0
      call ftibin(unit,nrows,tfields,ttype,tform,tunit,extname,varidat,status)

      ! Write dilations to binary table.
      frow=1
      felem=1
      colnum=1
      call ftpcle(unit,colnum,frow,felem,nrows,tr%dilation(:,1),status)
      colnum=2
      call ftpcle(unit,colnum,frow,felem,nrows,tr%dilation(:,2),status)

      ! Insert additional image extensions for wavelet coefficients at 
      ! each scale.

      ! Define parameters for wavelet coefficient data arrays to save.
      naxis=3
      naxes_wcoeff(1)=tr%n_alpha
      naxes_wcoeff(2)=tr%n_beta
      naxes_wcoeff(3)=tr%n_gamma

      ! Calculate admissibility for all dilations.
      if(tr%swav_mother_status &
           .and. cswt_swav_get_map_status(tr%swav_mother)) then
         temp_sky = cswt_swav_get_sky(tr%swav_mother)
         norm_preserve = .true.
         allocate(admiss(1:tr%n_dilation), stat=fail)
         if(fail /= 0) then
            call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
                 'cswt_tr_io_fits_write_wcoeff')
         end if
         call s2_sky_admiss_dil(temp_sky, &
              tr%dilation(:,1), tr%dilation(:,2), &
              admiss, norm_preserve)
         call s2_sky_free(temp_sky)
      end if

      do i_dil = 1,tr%n_dilation

        ! Insert a new image extension.
        call ftiimg(unit,bitpix,naxis,naxes_wcoeff,status)

        ! Write additional header keywords.
        call ftpkye(unit,'DIL1',tr%dilation(i_dil,1),decimals, &
          'scale of coefficients',status)
        call ftpkye(unit,'DIL2',tr%dilation(i_dil,2),decimals, &
          'scale of coefficients',status)      
        if(tr%swav_mother_status &
             .and. cswt_swav_get_map_status(tr%swav_mother)) then
          call ftpkye(unit,'ADMISS',admiss(i_dil),decimals, &
            'wavelet admissibility',status)
          call ftpkye(unit,'ADMISSPC',(1.0e0-abs(admiss(i_dil)))*100e0, &
            decimals, 'wavelet admissibility percentage',status)
        end if
        call ftpkyj(unit,'NALPHA',tr%n_alpha, &
          'length of alpha Euler angle dimension (=NAXIS1)',status)
        call ftpkyj(unit,'NBETA',tr%n_beta, &
         'length of beta Euler angle dimension (=NAXIS2)',status)
        call ftpkyj(unit,'NGAMMA',tr%n_gamma, &
         'length of gamma Euler angle dimension (=NAXIS3)',status)
        call ftpkys(unit,'EXTNAME','WCOEFF','entension name',status)

        ! Write wavelet coefficients for particular scale as a 3D data cube.
        group = 1
        dim1 = tr%n_alpha
        dim2 = tr%n_beta
        call ftp3de(unit, group, dim1, dim2, &
          tr%n_alpha, tr%n_beta, tr%n_gamma, &
          tr%wcoeff(i_dil,:,:,:), status)

      end do

      ! Free memory associated with admissibility values calculated.
      if(tr%swav_mother_status &
           .and. cswt_swav_get_map_status(tr%swav_mother)) then
        deallocate(admiss)
      end if

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call cswt_tr_io_fits_error_check(status, .true.)

    end subroutine cswt_tr_io_fits_write_wcoeff


    !--------------------------------------------------------------------------
    ! cswt_tr_io_fits_read_wcoeff
    !
    !! Read a fits wavelet coefficient file and allocate a new tr structure
    !! with the data read.
    !!
    !! Notes:
    !!   - The tr data structure is allocated herein by calling 
    !!     cswt_tr_init_computed
    !!   - No wavelet specific attributes are read from headers here since
    !!     the wavelet itself is not contained in the file.  In some cases 
    !!     the file may also be written from a previously read file (after
    !!     manipulation) and wavelet attributes may not be present.
    !!
    !! Variables:
    !!   - filename: Name of the fits input file containing the wavelet
    !!     transform data to be read.
    !!   - tr: The initialised tr structure containing the data read from
    !!     the input fits file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_fits_read_wcoeff(filename, tr)

      character(len=*), intent(in) :: filename
      type(cswt_tr), intent(out) :: tr

      real(s2_sp), parameter :: TOL = 1e-4
      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: nhdu, hdutype, naxis, n_check
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: colnum, frow, felem, nelem
      integer :: group, dim1, dim2
      real(s2_sp) :: nullval, dil_check
      integer :: i_dil, n_dilation, n_alpha, n_beta, n_gamma, fail
      real(s2_sp), allocatable :: dilation(:,:), wcoeff(:,:,:,:)
      integer :: lmax, mmax

      ! Check object not already initialised.
      if(tr%init) then
        call cswt_error(CSWT_ERROR_INIT, 'cswt_tr_io_fits_read_wcoeff')
        return
      end if

      status=0   ! Initialse error status to zero.

      ! Check if file already exists.
      call cswt_tr_io_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
           'cswt_tr_io_fits_read_wcoeff', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)


      ! --------------------------------------
      ! Read primary header.
      ! --------------------------------------

      ! Read size of image from variable size keywords (e.g. NDIL, NALPHA etc.)
      ! (read individually).
      call ftgkyj(unit, 'NDIL', n_dilation, comment, status)
      call ftgkyj(unit, 'NALPHA', n_alpha, comment, status)
      call ftgkyj(unit, 'NBETA', n_beta, comment, status)
      call ftgkyj(unit, 'NGAMMA', n_gamma, comment, status)
      call ftgkyj(unit, 'LMAX', lmax, comment, status)
      call ftgkyj(unit, 'MMAX', mmax, comment, status)

      ! Parameters will not be read since associated with spherical wavelet
      ! which don't have.

      ! Check correct number of HDUs in input file.
      ! (First two extension due to primary and dilation, then have one
      ! extension for coefficients at each scale.)
      hdunum = 1 + 1 + n_dilation
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
         call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
           'cswt_tr_io_fits_read_wcoeff', &
           comment_add='Invalid number of headers')
      end if

      ! Allocate wcoeff.
      allocate(wcoeff(1:n_dilation, 0:n_alpha-1, 0:n_beta-1, 0:n_gamma-1), &
        stat=fail)
      if(fail /= 0) then 
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_io_fits_read_wcoeff')
      end if
      ! Initialise with zeros.
      wcoeff = 0.0e0 

      ! Allocate dilation.
      allocate(dilation(n_dilation, CSWT_SWAV_DILATION_DIM), stat=fail)
      if(fail /= 0) then 
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_io_fits_write_wcoeff')
      end if
      ! Initialise with zeros.
      dilation = 0.0e0 

      ! Don't read admissibility since associated with spherical wavelet which
      ! don't have.


      ! --------------------------------------
      ! Read dilations.
      ! --------------------------------------

      ! Move to second extension (binary table containing dilation values).
      nhdu = 2  ! Dilations stored in second extension.
      call ftmahd(unit, nhdu, hdutype, status)

      ! Check correct hdutype.
      if(hdutype /= 2) then
         call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
           'cswt_tr_io_fits_read_wcoeff', &
           comment_add='Dilations not stored in binary table')
      end if

      ! Read header NAXIS2 and check same as n_a.
      call ftgkyj(unit, 'NAXIS2', naxis, comment, status)
      if(naxis/=n_dilation) then
         call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
           'cswt_tr_io_fits_read_wcoeff', &
           comment_add='Inconsistent number of dilations')
      end if

      ! Read dilation values from binary table.
      frow=1
      felem=1
      nelem=n_dilation
      nullval = -999      ! Arbitrary since will stop and return error 
                          ! if null values detected.
      colnum=1      
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           dilation(:,1),anynull,status)
      colnum = 2
      call ftgcve(unit,colnum,frow,felem,nelem,nullval, &
           dilation(:,2),anynull,status)
      if(anynull) then
         call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
           'cswt_tr_io_fits_read_wcoeff', &
           comment_add='Null dilation values contained in file')
      end if


      ! --------------------------------------
      ! Read wavelet coefficients.
      ! --------------------------------------

      ! Read wavelet coefficients at each dilation scale.

      do i_dil = 1,n_dilation

        ! Move to next extension.
         nhdu = nhdu + 1 
         call ftmahd(unit, nhdu, hdutype, status)  

         ! Read header and check correct sizes.
         call ftgkyj(unit, 'NAXIS1', n_check, comment, status)
         if(n_check/=n_alpha) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid coefficient image size')
         end if
         call ftgkyj(unit, 'NAXIS2', n_check, comment, status)
         if(n_check/=n_beta) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid coefficient image size')
         end if
         call ftgkyj(unit, 'NAXIS3', n_check, comment, status)
         if(n_check/=n_gamma) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid coefficient image size')
         end if
         call ftgkyj(unit, 'NALPHA', n_check, comment, status)
         if(n_check/=n_alpha) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid coefficient image size')
         end if
         call ftgkyj(unit, 'NBETA',  n_check, comment, status)
         if(n_check/=n_beta) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid coefficient image size')
         end if
         call ftgkyj(unit, 'NGAMMA', n_check, comment, status)
         if(n_check/=n_gamma) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid coefficient image size')
         end if
         call ftgkye(unit, 'DIL1', dil_check, comment, status)

         if(abs(dil_check - dilation(i_dil,1)) > TOL) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid dilation value')
         end if
         call ftgkye(unit, 'DIL2', dil_check, comment, status)
         if(abs(dil_check - dilation(i_dil,2)) > TOL) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Invalid dilation value')
         end if

         ! Read coefficients as 3D data cube.
         group = 1
         nullval = -999
         dim1 = n_alpha
         dim2 = n_beta
         call ftg3de(unit, group, nullval, dim1, dim2, &
           n_alpha, n_beta, n_gamma, &
           wcoeff(i_dil,:,:,:), anynull, status)
         if(anynull) then
            call cswt_error(CSWT_ERROR_TR_FILE_INVALID, &
              'cswt_tr_io_fits_read_wcoeff', &
              comment_add='Null dilation values contained in file')
         end if

      end do

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call cswt_tr_io_fits_error_check(status, .true.)

      ! Initialise tr data structure.
      tr = cswt_tr_init_computed(dilation, wcoeff, lmax, mmax)

      ! Free temporary storage used.
      deallocate(dilation)
      deallocate(wcoeff)

    end subroutine cswt_tr_io_fits_read_wcoeff


    !--------------------------------------------------------------------------
    ! cswt_tr_io_fits_error_check
    !
    !! Check if a fits error has occured and print error message.  Halt
    !! program execution if halt flag is set.
    !!
    !! Variables:
    !!   - status: Fits integer status code.
    !!   - halt: Logical to indicate whether to halt program execution if an 
    !!     error is detected.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_fits_error_check(status, halt)

      integer, intent(inout) :: status
      logical, intent(in) :: halt

      character(len=30) :: errtext
      character(len=80) :: errmessage

      !  Check if status is OK (no error); if so, simply return.
      if (status .le. 0) return

      ! The FTGERR subroutine returns a descriptive 30-character text 
      ! string that corresponds to the integer error status number.  
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      ! The FTGMSG subroutine retrieves the oldest message from
      ! the stack and shifts any remaining messages on the stack down one
      ! position.  FTGMSG is called repeatedly until a blank message is
      ! returned, which indicates that the stack is empty.  Each error message
      ! may be up to 80 characters in length. 
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          write(*,*) trim(errmessage)
          call ftgmsg(errmessage)
      end do

      if(halt) stop

    end subroutine cswt_tr_io_fits_error_check


    !--------------------------------------------------------------------------
    ! cswt_tr_io_fits_exists
    !
    !! Check if a fits file exists.
    !!
    !! Variables:
    !!   - filename: Name of fits file to check existence of.
    !!   - status: Fits integer status code.
    !!   - exists: Logical indicating whether the fits file already exists.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_fits_exists(filename, status, exists)

      character(len=*), intent(in) :: filename
      integer, intent(inout) :: status
      logical, intent(out) :: exists

      integer :: unit, blocksize
      logical :: halt

      ! Simply return if status is already greater than zero.
      if (status .gt. 0) return

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      call ftopen(unit, filename, 1, blocksize, status)

      ! Check status of opening file.
      if(status == 0) then

        ! File was opened.  Close it and set exists flag accordingly.
        call ftclos(unit, status)
        exists = .true.

      else if (status == 104) then
        
        ! File does not exist.  Reset status and set exists flag accordingly.
         status = 0
         exists = .false.

      else

        ! Some other error occured while opening file.
        halt = .false.
        call cswt_tr_io_fits_error_check(status, halt)
        call ftclos(unit, status)
        status = 0
        exists = .true.

      end if

      ! Deallocate unit number.
      call ftfiou(unit, status)

    end subroutine cswt_tr_io_fits_exists


    !--------------------------------------------------------------------------
    ! cswt_tr_io_fits_del
    !
    !! Delete a fits file.
    !!
    !! Variables:
    !!   - filename: Name of fits file to detele.
    !!   - status: Fits integer status code.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_fits_del(filename, status)

      character(len=*), intent(in) :: filename
      integer, intent(inout) ::  status

      integer :: unit, blocksize

      ! Simply return if status is greater than zero.
      if (status .gt. 0)return

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Try to open the file, to see if it exists.
      call ftopen(unit,filename,1,blocksize,status)

      if(status .eq. 0) then
         ! File was opened;  so now delete it.
         call ftdelt(unit,status)
      else if(status .eq. 103) then
         ! File doesn't exist, so just reset status to zero and clear errors.
          status=0
          call ftcmsg
      else
         ! There was some other error opening the file; delete the file anyway.
         status=0
         call ftcmsg
         call ftdelt(unit,status)
      end if

      ! Free the unit number for later reuse.
      call ftfiou(unit, status)

    end subroutine cswt_tr_io_fits_del


    !--------------------------------------------------------------------------
    ! cswt_tr_io_txt_dilation_write
    !
    !! Write tr structure dilation values to an output text file.
    !!
    !! Notes:
    !!   - Units stored in tr%dilation is always in radians.  However 
    !!     units read and written to file may be in either radians or arcmins.
    !!
    !! Variables:
    !!   - filename: Name of file to write dilations to.
    !!   - tr: Tr struture containing the dilations to be written.
    !!   - unit: Integer specifying the units to write the dilation values in
    !!     (either in arcmins or radians depending of the global 
    !!     CSWT_TR_DIL_UNIT_* specifier passed).
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_txt_dilation_write(filename, tr, unit)
 
      character(len=*), intent(in) :: filename
      type(cswt_tr), intent(in) :: tr
      integer, intent(in) :: unit

      character(len=1), parameter :: COMMENT_CHAR = '#'
      integer :: fileid = 20, i_dil

      open(unit=fileid, file=filename, status='replace', action='write')

      write(fileid,'(a,a)') COMMENT_CHAR, ' Dilation file'
      write(fileid,'(a,a)') COMMENT_CHAR, ' Written by cswt-0.2'

      if(unit == CSWT_TR_DIL_UNIT_RAD) then
         write(fileid,'(a,a12)') 'unit: ', CSWT_TR_DIL_UNIT_STR_RAD
      else if(unit == CSWT_TR_DIL_UNIT_ARCMIN) then
         write(fileid,'(a,a12)') 'unit: ', CSWT_TR_DIL_UNIT_STR_ARCMIN
      else
         call cswt_error(CSWT_ERROR_TR_DIL_UNIT_INVALID, &
           'cswt_tr_io_txt_dilation_write')
      end if

      write(fileid,'(a,i12)') 'n_dil:', tr%n_dilation

      do i_dil = 1,tr%n_dilation
         
         if(unit == CSWT_TR_DIL_UNIT_RAD) then
            write(fileid,'(i6,f12.4,f12.4)') i_dil, &
              tr%dilation(i_dil,1), tr%dilation(i_dil,2)
         else if(unit == CSWT_TR_DIL_UNIT_ARCMIN) then
            write(fileid,'(i6,f12.4,f12.4)') i_dil, &
              s2_vect_rad_to_arcmin(tr%dilation(i_dil,1)), &
              s2_vect_rad_to_arcmin(tr%dilation(i_dil,2))
         end if

      end do

      close(fileid)

    end subroutine cswt_tr_io_txt_dilation_write


    !--------------------------------------------------------------------------
    ! cswt_tr_io_txt_dilation_read
    !
    !! Read dilation values from an input dilation text file.
    !! 
    !! Notes:
    !!   - Units stored in tr%dilation is always in radians.  However 
    !!     units read and written to file may be in either radians or arcmins.
    !!   - The dilation array is allocated herein and must be freed by the
    !!     calling routine.
    !!
    !! Variables:
    !!   - filename: Name of the input dilation file to be read.
    !!   - dilation: Array of dilation values read from the input file.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_io_txt_dilation_read(filename, dilation)

      character(len=*), intent(in) :: filename
      real(s2_sp), intent(out), allocatable :: dilation(:,:)

      character(len=1), parameter :: COMMENT_CHAR = '#'
      integer :: fileid = 101, fail, n_dilation, i, idum
      character(len=S2_STRING_LEN) :: line, line2, unit_str
      real(s2_sp) :: dilation1_temp, dilation2_temp

      open(fileid, file=filename, form='formatted', status='old')

      ! Ignore leading comment lines.
      line = COMMENT_CHAR
      do while(line(1:1) == COMMENT_CHAR)
         read(fileid,'(a)') line
      end do

      ! Read units from last line read.
      read(line, *) line2, unit_str

      ! Read number of dilations from last line read (that is actually
      ! not a comment line).
      read(fileid, *) line, n_dilation

      ! Allocate space for dilations.
      allocate(dilation(1:n_dilation, CSWT_SWAV_DILATION_DIM), stat=fail)
      if(fail /= 0) then
         call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, &
           'cswt_tr_io_txt_dilation_read')
      end if

      ! Check unit string valid.
      if(trim(unit_str) /= CSWT_TR_DIL_UNIT_STR_RAD .and. &
         trim(unit_str) /= CSWT_TR_DIL_UNIT_STR_ARCMIN) then
         call cswt_error(CSWT_ERROR_TR_DIL_UNIT_INVALID, &
           'cswt_tr_io_txt_dilation_read')
      end if

      ! Read dilations.
      do i = 1,n_dilation
         read(fileid,*) idum, dilation1_temp, dilation2_temp
         if(trim(unit_str) == CSWT_TR_DIL_UNIT_STR_RAD) then
            dilation(i,1) = dilation1_temp
            dilation(i,2) = dilation2_temp
         else if(trim(unit_str) == CSWT_TR_DIL_UNIT_STR_ARCMIN) then
            dilation(i,1) = s2_vect_arcmin_to_rad(dilation1_temp)
            dilation(i,2) = s2_vect_arcmin_to_rad(dilation2_temp)
         end if
      end do

      close(fileid)

    end subroutine  cswt_tr_io_txt_dilation_read


    !--------------------------------------------------------------------------
    ! Get routines
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tr_get_init
    !
    !! Get init variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - init: Copy of object init variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_init(tr) result(init)

      type(cswt_tr), intent(in) :: tr
      logical :: init

      init = tr%init

    end function cswt_tr_get_init


    !--------------------------------------------------------------------------
    ! cswt_tr_get_swav_mother
    !
    !! Get copy of swav_mother variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - swav_mother: Copy of spherical mother wavelet extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    function cswt_tr_get_swav_mother(tr) result(swav_mother)

      type(cswt_tr), intent(in) :: tr
      type(cswt_swav) :: swav_mother

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_swav_mother')
      end if 

     ! Check wavelet present
      if(.not. tr%swav_mother_status) then
         call cswt_error(CSWT_ERROR_TR_SWAV_NPRES, 'cswt_tr_get_swav_mother')
      end if

      ! Make a copy for the returned mother spherical wavelet.
      swav_mother = cswt_swav_init(tr%swav_mother)

    end function cswt_tr_get_swav_mother


    !--------------------------------------------------------------------------
    ! cswt_tr_get_wcoeff
    !
    !! Get copy of the wavelet coefficients.
    !!
    !! Notes:
    !!   - Note that the wcoeff data array is indexed from *one* for the
    !!     dilation dimensions, but from *zero* for the Euler angle
    !!     dimensions.
    !!   - The size of the input wcoeff array is not checked and must adhere 
    !!     to the dimensions of the tr%wcoeff array (and be index similarly 
    !!     as described by the previous note). 
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - swav_mother: Copy of spherical mother wavelet extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_tr_get_wcoeff(tr, wcoeff)

      type(cswt_tr), intent(in) :: tr
      real(s2_sp), intent(out) :: wcoeff(1:,0:,0:,0:)

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_wcoeff')
      end if 

      ! Check wavelet coefficients calculated.
      if(.not. tr%wcoeff_status) then
         call cswt_error(CSWT_ERROR_TR_WCOEFF_NCOMP, 'cswt_tr_get_wcoeff')
      end if

      ! Get copy of wavelet coefficients.
      wcoeff = tr%wcoeff

    end subroutine cswt_tr_get_wcoeff


    !--------------------------------------------------------------------------
    ! cswt_tr_get_dilation
    !
    !! Get copy of dilation variable from the passed tr.
    !!
    !! Notes:
    !!   - The size of the input dilation array is not checked and must adhere 
    !!     to the dimensions of the tr%dilation array.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - dilation: Copy of object dilation variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine cswt_tr_get_dilation(tr, dilation)

      type(cswt_tr), intent(in) :: tr
      real(s2_sp), intent(out) :: dilation(:,:)

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_dilation')
      end if 

      dilation = tr%dilation

    end subroutine cswt_tr_get_dilation


    !--------------------------------------------------------------------------
    ! cswt_tr_get_lmax
    !
    !! Get copy of lmax variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - lmax: Copy of object lmax variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_lmax(tr) result(lmax)

      type(cswt_tr), intent(in) :: tr
      integer :: lmax

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_lmax')
      end if 

      lmax = tr%lmax

    end function  cswt_tr_get_lmax


    !--------------------------------------------------------------------------
    ! cswt_tr_get_mmax
    !
    !! Get copy of mmax variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - mmax: Copy of object mmax variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_mmax(tr) result(mmax)

      type(cswt_tr), intent(in) :: tr
      integer :: mmax

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_mmax')
      end if 

      mmax = tr%mmax

    end function  cswt_tr_get_mmax


    !--------------------------------------------------------------------------
    ! cswt_tr_get_n_alpha
    !
    !! Get copy of n_alpha variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - n_alpha: Copy of object n_alpha variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_n_alpha(tr) result(n_alpha)

      type(cswt_tr), intent(in) :: tr
      integer :: n_alpha

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_n_alpha')
      end if 

      n_alpha = tr%n_alpha

    end function  cswt_tr_get_n_alpha


    !--------------------------------------------------------------------------
    ! cswt_tr_get_n_beta
    !
    !! Get copy of n_beta variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - n_beta: Copy of object n_beta variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_n_beta(tr) result(n_beta)

      type(cswt_tr), intent(in) :: tr
      integer :: n_beta

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_n_beta')
      end if 

      n_beta = tr%n_beta

    end function  cswt_tr_get_n_beta


    !--------------------------------------------------------------------------
    ! cswt_tr_get_n_gamma
    !
    !! Get copy of n_gamma variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - n_gamma: Copy of object n_gamma variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_n_gamma(tr) result(n_gamma)

      type(cswt_tr), intent(in) :: tr
      integer :: n_gamma

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_n_gamma')
      end if 

      n_gamma = tr%n_gamma

    end function  cswt_tr_get_n_gamma


    !--------------------------------------------------------------------------
    ! cswt_tr_get_n_dilation
    !
    !! Get copy of n_dilation variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - n_dilaiton: Copy of object n_dilation variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_n_dilation(tr) result(n_dilation)

      type(cswt_tr), intent(in) :: tr
      integer :: n_dilation

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_n_dilation')
      end if 

      n_dilation = tr%n_dilation

    end function  cswt_tr_get_n_dilation


    !--------------------------------------------------------------------------
    ! cswt_tr_get_wcoeff_status
    !
    !! Get copy of wcoeff_status variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - wcoeff_status: Copy of object wcoeff_status variable extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_wcoeff_status(tr) result(wcoeff_status)

      type(cswt_tr), intent(in) :: tr
      logical :: wcoeff_status

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_wcoeff_status')
      end if 

      wcoeff_status = tr%wcoeff_status

    end function  cswt_tr_get_wcoeff_status


    !--------------------------------------------------------------------------
    ! cswt_tr_get_swav_mother_status
    !
    !! Get copy of swav_mother_status variable from the passed tr.
    !!
    !! Variables:
    !!   - tr: Spherical wavelet transform object to get copy of variable of.
    !!   - swav_mother_status: Copy of object swav_mother_status variable 
    !!     extracted.
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    function cswt_tr_get_swav_mother_status(tr) result(swav_mother_status)

      type(cswt_tr), intent(in) :: tr
      logical :: swav_mother_status

      ! Check object initialised.
      if(.not. tr%init) then
        call cswt_error(CSWT_ERROR_NOT_INIT, 'cswt_tr_get_swav_mother_status')
      end if 

      swav_mother_status = tr%swav_mother_status

    end function  cswt_tr_get_swav_mother_status


end module cswt_tr_mod
