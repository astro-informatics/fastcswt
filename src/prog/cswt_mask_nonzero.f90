!------------------------------------------------------------------------------
! cswt_mask_nonzero -- CSWT mask_nonzero program
!
!! Compute number of effective coefficients in an extended coefficient
!! exclusion mask and write results to either the standard output or to a file.
!!
!! Usage: cswt_mask_nonzero
!!   - [-help]: Display usage information.
!!   - [-inp filename_mask]: Name of file containing extended coefficient
!      mask (in wavelet coefficient fits file).
!!   - [-out filename_out]: Name of file to write results to.  If not 
!!     specified then results are written to the screen.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_mask_nonzero

  use s2_types_mod
  use cswt_error_mod
  use cswt_tr_mod

  implicit none

  character(len=*), parameter :: SCREEN = 'screen'
  character(len=S2_STRING_LEN) :: filename_mask, filename_out
  type(cswt_tr) :: tr_mask
  integer, allocatable :: ncoeff_eff(:,:)
  integer :: ncoeff_tot

  integer :: i_dil, i_gamma, fail, fileid

  ! Set default parameters.
  filename_out = SCREEN

  ! Parse options from command line.
  call parse_options() 

  ! Read in mask
  tr_mask = cswt_tr_init(filename_mask)

  ! Allocate space for number of effective coefficient array.
  allocate(ncoeff_eff(1:cswt_tr_get_n_dilation(tr_mask), &
    1:cswt_tr_get_n_gamma(tr_mask)), stat=fail)
  if(fail /= 0) then
     call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_mask_nonzero')
  end if

  ! Compute number of effective coefficients.
  call cswt_tr_mask_nonzero(tr_mask, ncoeff_eff, ncoeff_tot)

  ! Write results.
  if(trim(filename_out) == SCREEN) then
     
     write(*,'(a)') 'EFFECTIVE COEFFICIENTS OF EXTENDED MASK'
     write(*,'(a)') 'Created by cswt_mask_nonzero'
     write(*,*)

     write(*,'(a,a55)') 'Mask file: ', trim(filename_mask)
     write(*,'(a,i27)') 'Number of alpha samples:               ', &
          cswt_tr_get_n_alpha(tr_mask)
     write(*,'(a,i27)') 'Number of beta samples:                ', &
          cswt_tr_get_n_beta(tr_mask)
     write(*,'(a,i27)') 'Number of gamma samples:               ', &
          cswt_tr_get_n_gamma(tr_mask)
     write(*,'(a,i27)') 'Number of dilation samples:            ', &
          cswt_tr_get_n_dilation(tr_mask)
     write(*,'(a,i27)') 'Total number of coefficients (per sky):', &
          ncoeff_tot

     write(*,*)
     write(*,'(a,a)') '-----------------------------------------------', &
          '-------------------'
     write(*,'(a21,a45)') 'Coefficient map', 'Effective number of coefficients'
     write(*,'(a45,a18)') 'Number', 'Proportion'
     write(*,'(a,a)') '-----------------------------------------------', &
          '-------------------'
     
     do i_dil = 1, cswt_tr_get_n_dilation(tr_mask)
        
        write(*,'(a,i3.3)') 'Dilation index i_dil=', i_dil
        
        do i_gamma = 1,cswt_tr_get_n_gamma(tr_mask)
           
           write(*,'(a,i2.2,i16,f16.4)') ' Orientation index i_gamma=', &
                i_gamma, &
                ncoeff_eff(i_dil, i_gamma), &
                ncoeff_eff(i_dil,i_gamma)/real(ncoeff_tot,s2_sp) 
           
        end do
     end do
     
     write(*,'(a,a)') '-----------------------------------------------', &
          '-------------------'
     
  else

     fileid = 200
     open(unit=fileid, file=filename_out, status='replace', action='write')
     
     write(fileid,'(a)') 'EFFECTIVE COEFFICIENTS OF EXTENDED MASK'
     write(fileid,'(a)') 'Created by cswt_mask_nonzero'
     write(fileid,*)
     
     write(fileid,'(a,a55)') 'Mask file: ', trim(filename_mask)
     write(fileid,'(a,i27)') 'Number of alpha samples:               ', &
          cswt_tr_get_n_alpha(tr_mask)
     write(fileid,'(a,i27)') 'Number of beta samples:                ', &
          cswt_tr_get_n_beta(tr_mask)
     write(fileid,'(a,i27)') 'Number of gamma samples:               ', &
          cswt_tr_get_n_gamma(tr_mask)
     write(fileid,'(a,i27)') 'Number of dilation samples:            ', &
          cswt_tr_get_n_dilation(tr_mask)
     write(fileid,'(a,i27)') 'Total number of coefficients (per sky):', &
          ncoeff_tot

     write(fileid,*)
     write(fileid,'(a,a)') '-----------------------------------------------', &
          '-------------------'
     write(fileid,'(a21,a45)') 'Coefficient map', &
          'Effective number of coefficients'
     write(fileid,'(a45,a18)') 'Number', 'Proportion'
     write(fileid,'(a,a)') '-----------------------------------------------', &
          '-------------------'
     
     do i_dil = 1, cswt_tr_get_n_dilation(tr_mask)
        
        write(fileid,'(a,i3.3)') 'Dilation index i_dil=', i_dil
        
        do i_gamma = 1,cswt_tr_get_n_gamma(tr_mask)
           
           write(fileid,'(a,i2.2,i16,f16.4)') ' Orientation index i_gamma=', &
                i_gamma, &
                ncoeff_eff(i_dil, i_gamma), &
                ncoeff_eff(i_dil,i_gamma)/real(ncoeff_tot,s2_sp) 
           
        end do
     end do
     
     write(fileid,'(a,a)') '-----------------------------------------------', &
          '-------------------'

     close(fileid)

  end if

  ! Free memory.
  deallocate(ncoeff_eff)
  call cswt_tr_free(tr_mask)


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
            write(*,'(a)') 'Usage: cswt_mask_nonzero [-inp filename_mask]'
            write(*,'(a)') '                         [-out filename_out]' 
            stop
          
          case ('-inp')
            filename_mask = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("Unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program cswt_mask_nonzero
