!------------------------------------------------------------------------------
! cswt_tr_nsigma -- CSWT tr_nsigma program
!
!! Convert wavelet coefficients to nsigma values.  The coefficient map for 
!! each dilation is considered; the mean for each dilation is subtracted, 
!! before dividing by the sigma value computed for that dilation.
!! Note that the mean and sigma are only computed over the non-masked pixels.
!!
!! Usage: cswt_analysis 
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of input fits wavelet coefficient file.
!!   - [-mask filename_mask]: Name of input mask file.
!!   - [-out filename_out]: Name of output wavelet coeffifficient file 
!!     containing the nsigma values.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - December 2004
!
! Revisions:
!   December 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_tr_nsigma

  use s2_types_mod
  use cswt_tr_mod
  use cswt_error_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_in, filename_out, filename_mask
  type(cswt_tr) :: tr, tr_mask
  integer :: n_dilation, fail
  real(s2_sp), allocatable :: sigma(:), mean(:)
  integer, allocatable :: ncoeff_eff(:,:), n_effective_samples(:)
  integer :: ncoeff_tot

  ! Parse options from command line.
  call parse_options()

  ! Read wavelet coefficient file.
  tr = cswt_tr_init(filename_in)
  n_dilation = cswt_tr_get_n_dilation(tr)

  ! Apply mask here.
  tr_mask = cswt_tr_init(filename_mask)
  call cswt_tr_mask_apply(tr, tr_mask)

  ! Compute effective number of coefficients for each dilation.
  ! Note that ncoeff_eff gives effective number for each dilation 
  ! and orientation, sum over orientations to get total number of 
  ! effective coefficients for each dilation.
  allocate(ncoeff_eff(n_dilation, cswt_tr_get_n_gamma(tr)), stat=fail)
  allocate(n_effective_samples(n_dilation), stat=fail)
  if(fail /= 0) then
     call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_nsigma')
  end if
  call cswt_tr_mask_nonzero(tr_mask, ncoeff_eff, ncoeff_tot)
  n_effective_samples = sum(ncoeff_eff,2)

  ! Compute mean and sigma for each dilation.
  allocate(sigma(n_dilation), stat=fail)
  allocate(mean(n_dilation), stat=fail)
  if(fail /= 0) then
     call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_tr_nsigma')
  end if
  call cswt_tr_wcoeff_sigma(tr, n_effective_samples, mean, sigma)

  ! Subtract mean of wavelet coefficients for each dilation.
  call cswt_tr_const_subtract(tr, mean)

  ! Scale wavelet coefficients by sigma to get nsigma maps
  ! for each dilation.
  call cswt_tr_const_multiply(tr, 1.0e0/sigma)

  ! Write nsigma wavelet coefficient maps.
  call cswt_tr_io_fits_write_wcoeff(filename_out, tr)

  ! Free memory.
  call cswt_tr_free(tr)
  call cswt_tr_free(tr_mask)
  deallocate(sigma, mean)
  deallocate(ncoeff_eff, n_effective_samples)


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
          write(*,*) 'option ', opt, ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: cswt_tr_nsigma [-inp filename_in]'
            write(*,'(a)') '                      [-mask filename_mask]' 
            write(*,'(a)') '                      [-out filename_out]' 
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-mask')
            filename_mask = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program cswt_tr_nsigma
