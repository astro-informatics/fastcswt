!------------------------------------------------------------------------------
! cswt_mask_apply -- CSWT mask_apply program
!
!! Apply an extended coefficient exclusion mask to wavelet coefficients.
!!
!! Notes:
!!   - *Important*: Only specify the display flag if the output may is 
!!     only used for display purposes.  If subsequent processing is 
!!     applied erroneous results amy occur.
!!
!! Usage: cswt_mask_apply
!!   - [-help]: Display usage information.
!!   - [-wcoeff filename_wcoeff]: Name of input file containing wavelet 
!!     coefficients.
!!   - [-mask filename_mask]: Name of input file containing extended 
!!     coefficient mask defined over the same domain as the wavelet 
!!     coefficients (hence saved in wavelet coefficient fits file).
!!   - [-display display_mask]: Logical to specified whether a display mask 
!!     is to be applied.  If that is the case then the pixels in the tr
!!     structure corresponding to the zeros with be replaced with a large 
!!     negative  value so that they appear grey in output images.
!!   - [-out filename_out]: Name of output file to write the masked wavelet
!!      coefficient to.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_mask_apply

  use s2_types_mod
  use cswt_tr_mod

  implicit none

  character(len=S2_STRING_LEN) :: &
    filename_wcoeff, &
    filename_mask, &
    filename_out

  type(cswt_tr) :: tr, tr_mask
  logical, parameter :: DISPLAY_MASK_DEFAULT = .false.
  logical :: display_mask = DISPLAY_MASK_DEFAULT

  ! Parse options from command line.
  call parse_options()

  ! Read in wavelet coefficients.
  tr = cswt_tr_init(filename_wcoeff)

  ! Read in extended coefficient mask (in wavelet coefficient structure file).
  tr_mask = cswt_tr_init(filename_mask)

  ! Apply mask.
  call cswt_tr_mask_apply(tr, tr_mask, display_mask)

  ! Write masked wavelet coefficients to output file.
  call cswt_tr_io_fits_write_wcoeff(filename_out, tr)

  ! Free memory.
  call cswt_tr_free(tr)
  call cswt_tr_free(tr_mask)
  

 !----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parse the options passed when program called.
    !
    !! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
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
          write(*,*) 'option ', opt, ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: cswt_mask_gen [-wcoeff filename_wcoeff]'
            write(*,'(a)') '                     [-mask filename_mask]'
            write(*,'(a)') '                     [-display display_mask]'
            write(*,'(a)') '                     [-out filename_out]' 
            stop
          
          case ('-wcoeff')
            filename_wcoeff = trim(arg)
         
          case ('-mask')
            filename_mask = trim(arg)

          case ('-display')
            read(arg,*) display_mask

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program cswt_mask_apply
