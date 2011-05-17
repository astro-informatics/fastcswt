!------------------------------------------------------------------------------
! cswt_error_mod -- CSWT library error class
!
!! Functionality to handle errors that may occur in the CSWT library.  Public
!! CSWT error codes are defined, with corresponding private error comments and 
!! default halt execution status.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.2 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module cswt_error_mod

  use s2_types_mod, only: S2_STRING_LEN

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: cswt_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: CSWT_ERROR_NUM = 27

  integer, public, parameter :: &
    CSWT_ERROR_NONE = 0, &
    CSWT_ERROR_INIT = 1, &
    CSWT_ERROR_NOT_INIT = 2, &
    CSWT_ERROR_INIT_FAIL = 3, &
    CSWT_ERROR_SIZE_INVALID = 4, &
    CSWT_ERROR_MEM_ALLOC_FAIL = 5, &
    CSWT_ERROR_TMPL_PARAM_INVALID = 6, &
    CSWT_ERROR_SWAV_FTYPE_INVALID = 7, &
    CSWT_ERROR_SWAV_DIL_ERROR = 8, &
    CSWT_ERROR_SWAV_ROT_ERROR = 9, &
    CSWT_ERROR_TR_AZBANDLIM_INVALID = 10, &
    CSWT_ERROR_TR_FILE_EXISTS = 11, &
    CSWT_ERROR_TR_FILE_INVALID = 12, &
    CSWT_ERROR_TR_LMAX_INVALID = 13, &
    CSWT_ERROR_TR_MMAX_INVALID = 14, &
    CSWT_ERROR_TR_NGAMMA_INVALID = 15, &
    CSWT_ERROR_TR_WCOEFF_NCOMP = 16, &
    CSWT_ERROR_TR_WCOEFF_COMP = 17, &
    CSWT_ERROR_TR_SWAV_NPRES = 18, &
    CSWT_ERROR_TR_DIL_UNIT_INVALID = 19, &
    CSWT_ERROR_TR_METHOD_INVALID = 20, &
    CSWT_ERROR_TR_KERNEL_INVALID = 21, &
    CSWT_ERROR_TR_MASK_MODE_INVALID = 22, &
    CSWT_ERROR_TR_THRESNSIG_MODE_INVALID = 23, &
    CSWT_ERROR_ANALYSIS = 24, &
    CSWT_ERROR_MASK_GEN = 25, &
    CSWT_ERROR_PLOT_SWAV = 26

  ! Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.  
  !! Comment associated with each error type.
  character(len=S2_STRING_LEN), parameter :: &
    error_comment(CSWT_ERROR_NUM) = &
      (/ & 
      'No error                                                                 ', &
      'Attempt to initialise object that has already been initialised           ', &
      'Object not initialised                                                   ', &
      'Object initialisation failed                                             ', &
      'Sizes invalid                                                            ', &
      'Memory allocation failed                                                 ', &
      'Invalid number of input parameters in template function                  ', &
      'Invalid function type                                                    ', &
      'Error dilating spherical wavelet                                         ', &
      'Error rotating spherical wavelet                                         ', &
      'Azimuthal band limit not met (increase N_gamma)                          ', &
      'File already exists                                                      ', &
      'File invalid                                                             ', &
      'Lmax invalid                                                             ', &
      'Mmax invalid                                                             ', &
      'N gamma invalid (must be odd)                                            ', &
      'Wavelet coefficients not computed                                        ', &
      'Wavelet coefficients already computed                                    ', &
      'Wavelet not present                                                      ', &
      'Invalid dilation unit                                                    ', &
      'Spherical wavelet transform method type invalid                          ',&
      'Kernel invalid                                                           ',&
      'Invalid mask method specifier                                            ',&
      'Invalid mask thresholding method specifier                               ',&
      'Error in cswt_analysis program                                           ', &
      'Error in cswt_mask_gen program                                           ', &
      'Error in cswt_plot_swav program                                          ' &
      /) 

  !! Default program halt status of each error type.
  logical, parameter :: &
    halt_default(CSWT_ERROR_NUM) = &
      (/ &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .false., &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.,  &
      .true.  /)
  
  
  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! cswt_error
    ! 
    !! Display error message corresponding to error_code and halt program 
    !! execution if required.
    !!
    !! Variables:
    !!   - error_code:
    !!   - [procedure]: Procedure name where s2_error called from.  Displayed 
    !!     when error message printed to screen.
    !!   - [comment_add]: If present, additional comment to append to default 
    !!     error comment.
    !!   - [comment_out]: If present the error comment is copied to comment_out
    !!     on output.
    !!   - [halt_in]: If present overrides default halt value.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine cswt_error(error_code, procedure, comment_add, &
      comment_out, halt_in)

      integer, intent(in) :: error_code
      character(len=*), intent(in), optional :: procedure, comment_add
      character(len=*), intent(inout), optional :: comment_out
      logical, intent(in), optional :: halt_in

      logical :: halt
      character(len=*), parameter :: comment_prefix = 'CSWT_ERROR: '

      !---------------------------------------
      ! Display error message
      !---------------------------------------

      if(present(procedure)) then

        if(present(comment_add)) then
          write(*,'(a,a,a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            '''', &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a,a,a,a,a)') comment_prefix, 'Error ''', &
            trim(error_comment(error_code+1)), &
            ''' occured in procedure ''', &
            trim(procedure), &
            ''''
        end if
 
     else

        if(present(comment_add)) then
          write(*,'(a,a,a,a)') comment_prefix, &
            trim(error_comment(error_code+1)), &
            ' - ', trim(comment_add)
        else
          write(*,'(a,a)') comment_prefix, trim(error_comment(error_code+1))
        end if

      end if

      ! Copy error comment if comment_out present.
      if(present(comment_out)) comment_out = error_comment(error_code+1)

      !---------------------------------------
      ! Halt program execution if required
      !---------------------------------------
      
      if( present(halt_in) ) then
        halt = halt_in
      else
        halt = halt_default(error_code+1)
      end if

      if( halt ) then
        write(*,'(a,a,a,a,a)') comment_prefix, &
          '  Halting program execution ', &
          'due to error ''', trim(error_comment(error_code+1)), ''''
        stop
      end if

    end subroutine cswt_error


end module cswt_error_mod
