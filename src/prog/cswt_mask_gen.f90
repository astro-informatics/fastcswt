!------------------------------------------------------------------------------
! cswt_mask_gen -- CSWT mask_gen program
!
!! Generate an extended coefficient exclusion mask defined in the ecp
!! wavelet coefficient domain from an original mask defined on the sky and 
!! the spherical wavelet transform of the original mask.
!!
!! Notes:
!!   - Morphological operation method is reccommended.  Otherwise parameters 
!!     must be tunned to application.
!!
!! Usage: cswt_mask_gen
!!   - [-help]: Display usage information.
!!   - [-wcoeff filename_wcoeff]: Name of input wavelet coefficient file 
!!     containing coefficients of the original mask considered.
!!   - [-sky filename_sky]: Name of the input mask sky file containing the
!!     original mask defined on the sky
!!   - [-ext file_extension]: File extension number mask sky is stored in.
!!   - [-morph mode_morph]: If true use morphological operator method.
!!   - [-convert mode_convert]: If true just convert sky mask to wavelet
!!     domain.
!!   - [-cswt mode_cswt]: If true use spherical wavelet coefficients and 
!!     thresholding technique.
!!   - [-out filename_out]: Name of output extended coefficient exclusion 
!!     mask file.
!!   - [-wav wavelet_type (only required if method is cswt)]: Wavelet type 
!!      mask corresponds to. 
!!     (Alters parameters used to generate the extended mask.)
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_mask_gen

  use s2_types_mod
  use s2_sky_mod
  use cswt_error_mod
  use cswt_tr_mod

  implicit none
  
  character(len=S2_STRING_LEN), parameter ::  &
    WAV_TYPE_MEXHAT = 'mexhat', &
    WAV_TYPE_MORLET = 'morlet'

  ! *Warning*: These parameters are not generic.  Need to be tuned for 
  ! application.  In general should use morph method, then parameters
  ! not required.
  real(s2_sp), parameter :: &
       THRES1_PROP_MEXHAT = 0.25e0, &
       THRES2_PROP_MEXHAT = 0.9e0, &
       THRES1_PROP_MORLET = 0.8e0, &  !! reduce slowly
       THRES2_PROP_MORLET = 0.01e0, & !! increase slowly
       THRES1_MIN_MEXHAT = 1.0e-3, &
       THRES1_MIN_MORLET = 5.0e-4

  character(len=S2_STRING_LEN) :: &
    filename_wcoeff, &
    filename_sky, &
    filename_out, &
    wavelet_type
  integer :: file_extension_sky
  type(s2_sky) :: sky_mask
  type(cswt_tr) :: tr_mask

  logical :: mode_cswt = .false.
  logical :: mode_morph = .false.
  logical :: mode_convert = .false.
  

  ! Set default parameter values.
  wavelet_type = WAV_TYPE_MEXHAT
  file_extension_sky = 1

  ! Parse options from command line.
  call parse_options()

  ! Check mode not set to more than one type
  if(     (mode_cswt .and. mode_morph) &
     .or. (mode_cswt .and. mode_convert) &
     .or. (mode_convert .and. mode_morph) ) then

     call cswt_error(CSWT_ERROR_MASK_GEN, 'cswt_mask_gen', &
       comment_add='More than one method set')

  end if

  ! Read wavelet coefficients of mask.
  tr_mask = cswt_tr_init(filename_wcoeff)

  ! Read original mask sky.
  sky_mask = s2_sky_init(filename_sky, file_extension_sky)

  ! Generate extended coefficient mask.
  if(mode_morph) then

     call cswt_tr_mask_gen(tr_mask, sky_mask, CSWT_TR_MASK_GEN_MORPH)

  elseif(mode_convert) then

     call cswt_tr_mask_gen(tr_mask, sky_mask, CSWT_TR_MASK_GEN_CSWT, &
       run_stage1_in=.false., &   
       run_stage2_in=.false., &
       orig_mask_only_in=.true.)

  elseif(mode_cswt) then

     ! Generate extended coefficient mask.
     ! (Choose parameters depending on wavelet).

     if(wavelet_type == trim(WAV_TYPE_MEXHAT)) then

        call cswt_tr_mask_gen(tr_mask, sky_mask, CSWT_TR_MASK_GEN_CSWT, &
             run_stage1_in=.true., &   
             run_stage2_in=.true., &
             orig_mask_only_in=.false., &
             thres1_proportion_in=THRES1_PROP_MEXHAT, &
             thres2_proportion_in=THRES2_PROP_MEXHAT, &
             thres1_min_in=THRES1_MIN_MEXHAT)
        
     elseif(wavelet_type == trim(WAV_TYPE_MORLET)) then
        
        call cswt_tr_mask_gen(tr_mask, sky_mask, CSWT_TR_MASK_GEN_CSWT, &
             run_stage1_in=.true., &  
             run_stage2_in=.true., &
             orig_mask_only_in=.false., &
             thres1_proportion_in=THRES1_PROP_MORLET, &
             thres2_proportion_in=THRES2_PROP_MORLET, &
             thres1_min_in=THRES1_MIN_MORLET)

! Add other wavelets here.
        
     else

        call cswt_error(CSWT_ERROR_MASK_GEN, 'cswt_mask_gen', &
             comment_add='Invalid wavelet type')
        
     end if

  else

     call cswt_error(CSWT_ERROR_MASK_GEN, 'cswt_mask_gen', &
          comment_add='No method specified')

  end if

  ! Write extended coefficient mask to output file.
  call cswt_tr_io_fits_write_wcoeff(filename_out, tr_mask)

  ! Free memory.
  call cswt_tr_free(tr_mask)
  call s2_sky_free(sky_mask)


 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    subroutine parse_options()

      use extension, only: getArgument, nArguments
     
      implicit none
      
      integer :: n, i
      character(len=S2_STRING_LEN) :: opt
      character(len=S2_STRING_LEN) :: arg
      
      n = nArguments()
     
      do i = 1,n,2
        
        call getArgument(i,opt)
     
        if (i == n .and. trim(opt) /= '-help') then
          write(*,'(a,a,a)') 'Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: cswt_mask_gen [-wcoeff filename_wcoeff]'
            write(*,'(a)') '                     [-sky filename_sky]'
            write(*,'(a)') '                     [-ext file_extension]'
            write(*,'(a)') '                     [-morph mode_morph]'
            write(*,'(a)') '                     [-convert mode_convert]'
            write(*,'(a)') '                     [-cswt mode_cswt]'
            write(*,'(a)') '                     [-out filename_out]' 
            write(*,'(a,a)') '                     ', &
              '[-wav wavelet_type (only required if method is cswt)]' 
            stop
          
          case ('-wcoeff')
            filename_wcoeff = trim(arg)
         
          case ('-sky')
            filename_sky = trim(arg)

          case ('-ext')
            read(arg,*) file_extension_sky

          case ('-morph')
            read(arg,*) mode_morph

          case ('-convert')
            read(arg,*) mode_convert

          case ('-cswt')
            read(arg,*) mode_cswt

          case ('-out')
            filename_out = trim(arg)

          case ('-wav')
            wavelet_type = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program cswt_mask_gen
