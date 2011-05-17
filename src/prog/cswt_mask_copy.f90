!------------------------------------------------------------------------------
! cswt_mask_copy -- CSWT mask_copy program
!
!! Generate a coefficient exclusion mask defined in the ecp wavelet 
!! coefficient domain by simply converting a list of Healpix sky mask to 
!! the wavelet domain.
!!
!! Notes:
!!   - The wavelet coefficients are not actually required in the 
!!     algorithm, however the wavelet transform structure of the correct
!!     size must be initialised.  To overcome this problem the tr_mask 
!!     structure containing the wavelet coefficient masks is often used.
!!
!! Usage: cswt_mask_copy
!!   - [-help]: Display usage information.
!!   - [-wcoeff filename_wcoeff]: Name of input wavelet coefficient file 
!!     initialised with wavelet coefficient tr structure of correct size
!!     (convenient to just use saved wavelet transform coefficient file
!!      here, see note above).
!!   - [-ndil n_dilation]: Number of dilations.  (Must be the same as the 
!!     number of sky files read and the same size as defined in the 
!!     spherical wavelet coefficient tr stucture.
!!   - [-sky filename_sky(i) (repeat n_dil times)]: Sky mask filenames.
!!     Parameters are read sequentially on order appear in command line.
!!   - [-ext file_extension]: File extension number sky masks are stored in.
!!   - [-out filename_out]: Name of output extended coefficient exclusion 
!!     mask file.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - December 2004
!
! Revisions:
!   December 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_mask_copy

  use s2_types_mod
  use s2_sky_mod
  use cswt_error_mod
  use cswt_tr_mod

  implicit none
  
  character(len=S2_STRING_LEN) :: filename_wcoeff, filename_out
  character(len=S2_STRING_LEN), allocatable :: filename_sky(:)

  integer :: file_extension_sky = 1
  integer :: i_dil = 0, n_dilation, fail
  type(s2_sky), allocatable :: sky_mask(:)
  type(cswt_tr) :: tr_mask
  
  ! Parse options from command line.
  call parse_options()

  ! Read wavelet coefficients of mask to get initialised tr structure of
  ! correct size.
  tr_mask = cswt_tr_init(filename_wcoeff)

  ! Check tr structure has correct size, i.e. same number of dilations as 
  ! read from command line.
  if(n_dilation /= cswt_tr_get_n_dilation(tr_mask)) then
     call cswt_error(CSWT_ERROR_SIZE_INVALID, 'cswt_mask_copy', &
       comment_add='Dilation sizes inconsistent')
  end if

  ! Allocate space for sky masks.
  allocate(sky_mask(n_dilation), stat=fail)
  if(fail /= 0) then
     call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_mask_copy')
  end if

  ! Read original sky masks for each dilation.
  do i_dil = 1,n_dilation
     sky_mask(i_dil) = s2_sky_init(filename_sky(i_dil), file_extension_sky)
  end do

  ! Convert original masks to wavelet domain.
  call cswt_tr_mask_copy(tr_mask, sky_mask)

  ! Write extended coefficient mask to output file.
  call cswt_tr_io_fits_write_wcoeff(filename_out, tr_mask)

  ! Free memory.
  call cswt_tr_free(tr_mask)
  do i_dil = 1,n_dilation
     call s2_sky_free(sky_mask(i_dil)) 
  end do
  deallocate(sky_mask)
  deallocate(filename_sky)


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
          write(*,*) 'option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: cswt_mask_copy [-wcoeff filename_wcoeff]'
            write(*,'(a)') '                      [-ndil n_dilation]'
            write(*,'(a,a)') '                      ', &
              '[-sky filename_sky(i) (repeat n_dil times)]'
            write(*,'(a)') '                      [-ext file_extension]'
            write(*,'(a)') '                      [-out filename_out]' 
            stop
          
          case ('-wcoeff')
            filename_wcoeff = trim(arg)
         
         case ('-ndil')
            read(arg,*) n_dilation
            ! Allocate space for parameters.
            allocate(filename_sky(n_dilation), stat=fail)
            if(fail /= 0) then
               call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_mask_copy')
            end if

          case ('-sky')
            ! Parameters read in order appear in command line.
            i_dil = i_dil + 1
            filename_sky(i_dil) = trim(arg)

          case ('-ext')
            read(arg,*) file_extension_sky

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a," ignored")', trim(opt)

        end select
      end do

    end subroutine parse_options


end program cswt_mask_copy
