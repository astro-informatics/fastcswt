!------------------------------------------------------------------------------
! cswt_tr2sky -- CSWT tr2sky program
!
!! Convert a spherical wavelet transform fits file (ecp sampled) to a sky 
!! Healpix fits file.  Just the sky corresponding to a specified dilation 
!! and orientation index may be written, or skies corresponding to all
!! dilations and orientations may be written.
!!
!! Usage: cswt_tr2sky
!!   - [-help]: Display usage information.
!!   - [-inp filename_data]: Name of input wavelet coefficient fits file.
!!   - [-out filename_out]: Name of output sky coefficient map (only required 
!!     if all=.false.)
!!   - [-idil dilation_index]: Dilation index of coefficients to write (only 
!!     required if all=.false.)
!!   - [-igamma orientation_index]: Orientation index of coefficients to write
!!     (only required if all=.false.)
!!   - [-nside healpix_nside]: Healpix nside resolution of written sky. 
!!     (Independent of resolution of coefficients since interpolation
!!      performed [or closest pixel taken].)
!!   - [-interp interpolation_status]: Logical to specify whether to perform
!!     linear interpolation (if true) or else simply take nearest pixel.
!!   - [-all all_status]: Logical to specify whether skies corresponding to 
!!     all dilations and orientations are to be written.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_tr2sky

  use s2_types_mod
  use cswt_tr_mod

  implicit none

  character(len=S2_STRING_LEN) :: filename_in, filename_out
  integer :: i_dil, i_gamma, nside
  logical :: interp, all

  integer, parameter :: FITS_FILENAME_EXT_LEN = 5
  type(cswt_tr) :: tr

  ! Set default parameters.
  nside = 128
  interp = .true.
  all = .false.

  ! Parse options from command line.
  call parse_options()

  ! Read file.
  tr = cswt_tr_init(filename_in)

  ! Write specified map(s).
  if(all) then

     ! Write all maps.
     do i_dil = 1,cswt_tr_get_n_dilation(tr)
       do i_gamma = 1,cswt_tr_get_n_gamma(tr)

         if(cswt_tr_get_n_dilation(tr) < 100) then
           write(filename_out, '(a,a,i2.2,a,i2.2,a)') &
             trim(filename_in(1:len(trim(filename_in)) &
                              - FITS_FILENAME_EXT_LEN)), &
             '_sky_id', i_dil, '_ig', i_gamma, '.fits'
          else
            write(filename_out, '(a,a,i3.3,a,i2.2,a)') &
              trim(filename_in(1:len(trim(filename_in)) &
                                - FITS_FILENAME_EXT_LEN)), &
              '_sky_id', i_dil, '_ig', i_gamma, '.fits'
         end if

         call cswt_tr_io_fits_write_wcoeff_sky(filename_out, tr, &
           i_dil, i_gamma, nside, interp)

        end do
     end do

  else

     ! Write single map specified.
     call cswt_tr_io_fits_write_wcoeff_sky(filename_out, tr, &
          i_dil, i_gamma, nside, interp)

  end if

  ! Free memory.
  call cswt_tr_free(tr)


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
            write(*,'(a)') 'Usage: cswt_tr2sky [-inp filename_data]'
            write(*,'(a)') '                   [-out filename_out]' 
            write(*,'(a)') '                   [-idil dilation_index]' 
            write(*,'(a)') '                   [-igamma orientation_index]' 
            write(*,'(a)') '                   [-nside healpix_nside]' 
            write(*,'(a)') &
              '                   [-interp interpolation_status]' 
            write(*,'(a)') &
              '                   [-all all_status]' 
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-idil')
            read(arg,*) i_dil

          case ('-igamma')
            read(arg,*) i_gamma

         case ('-nside')
            read(arg,*) nside

         case ('-interp')
            read(arg,*) interp

          case ('-all')
            read(arg,*) all

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program cswt_tr2sky
