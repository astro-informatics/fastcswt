!------------------------------------------------------------------------------
! cswt_mask_invert -- CSWT mask_apply program
!
!! Invert a coefficient mask by converting ones to zeros, and zeros to 
!! ones.
!!
!! Notes:
!!   - Original mask should be in one/zero format, not display format.
!!
!! Usage: cswt_mask_invert
!!   - [-help]: Display usage information.
!!   - [-inp filename_inp]: Name of input file containing original mask
!!     to be inverted.
!!   - [-out filename_out]: Name of output file created that will contain 
!!     the inverted mask.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - January 2005
!
! Revisions:
!   January 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_mask_invert

  use s2_types_mod
  use cswt_tr_mod

  implicit none

  character(len=S2_STRING_LEN) :: &
    filename_in, &
    filename_out

  type(cswt_tr) :: tr_mask

  ! Parse options from command line.
  call parse_options()

  ! Read in wavelet coefficients.
  tr_mask = cswt_tr_init(filename_in)

  ! Invert zeros and ones in mask.
  call cswt_tr_mask_invert(tr_mask)

  ! Write inverted mask to output file.
  call cswt_tr_io_fits_write_wcoeff(filename_out, tr_mask)

  ! Free memory.
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
            write(*,'(a)') 'Usage: cswt_mask_invert [-inp filename_in]'
            write(*,'(a)') '                        [-out filename_out]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)
         
          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program cswt_mask_invert

  
