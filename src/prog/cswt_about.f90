!------------------------------------------------------------------------------
! cswt_about -- CSWT about program
!
!! Display information about the CSWT package.
!!
!! Usage: cswt_about 
!!   - [-help]: Display usage information.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2011
!
! Revisions:
!   November 2011 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_about

  use s2_types_mod

  implicit none

  ! Display info.
  write(*,'(a)') "=========================================================="
  write(*,'(a)') "CSWT package to perform a fast directional continuous "
  write(*,'(a)') "wavelet transform on the sphere"
  write(*,'(a)') "By Jason McEwen"

  write(*,'(a)') "See www.jasonmcewen.org for more information."
  write(*,'(a)') "See LICENSE.txt for license details."

  write(*,'(a,a)') "Version: ", CSWT_VERSION
  write(*,'(a,a)') "Build: ", CSWT_BUILD
  write(*,'(a)') "=========================================================="


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
            write(*,'(a)') 'Usage: cswt_about'
            stop

          case default
            print '("Unknown option ",a," ignored")', trim(opt)            

        end select
      end do

    end subroutine parse_options


end program cswt_about
