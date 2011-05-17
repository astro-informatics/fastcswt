!------------------------------------------------------------------------------
! cswt_plot_swav -- CSWT plot_swav program
!
!! Plot spherical wavelet for a range of dilations.
!!
!! Usage: cswt_analysis 
!!   - [-help]: Display usage information.
!!   - [-dil filename_dilation]: Name of dilation file containing dilations
!!     to consider.
!!   - [-out filename_out_prefix]: Name prefix of output spherical wavelet
!!     fits file. (Actual output file is appended with '_id**.fits' 
!!     where ** is the dilation index.)
!!   - [-wav wavelet_type]: String specifying wavelet type.
!!   - [-nside nside]: Healpix nside resolution parameter.
!!   - [-cen center_on_equator]: Logical to specify whether to rotate
!!     spherical wavelet to equator for visualisation purposes.
!!   - [-norm norm_preserve]: Logical specifying whether to perform a
!!      norm preserving dilation.
!!   - [-npar n_param]: Number of wavelet parameters specified.
!!   - [-par param(i) (repeat n_param times)]: Wavelet parameters.  
!!     Parameters are read sequentially on order appear in command line.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_plot_swav

  use s2_types_mod
  use s2_sky_mod
  use cswt_error_mod
  use cswt_tmpl_mod
  use cswt_swav_mod
  use cswt_tr_mod

  implicit none

  character(len=S2_STRING_LEN), parameter ::  &
    WAV_TYPE_FILE = 'file', &
    WAV_TYPE_MEXHAT = 'mexhat', &
    WAV_TYPE_BUTTERFLY = 'butterfly', &
    WAV_TYPE_MORLET = 'morlet'
  character(len=S2_STRING_LEN) :: &
    filename_dilation, &
    filename_out_prefix, &
    filename_out, &
    wavelet_type

  logical :: center_on_equator, norm_preserve
  real(s2_sp), parameter :: ALPHA_CENTER = 0.0e0
  real(s2_sp), parameter :: BETA_CENTER = pi/2.0e0
  real(s2_sp), parameter :: GAMMA_CENTER = pi

  integer :: i_param = 0, n_param = 0
  integer :: nside, i_dil, fail
  real(s2_sp), allocatable :: param(:), dilation(:,:)
  type(cswt_swav) :: swav_mother, swav

  ! Set default parameters.
  wavelet_type = WAV_TYPE_MEXHAT
  nside = 128
  center_on_equator = .true.
  norm_preserve = .true.

  ! Parse options from command line.
  call parse_options()

  ! Check correct number of parameters were read.
  if(i_param /= n_param) then
     call cswt_error(CSWT_ERROR_PLOT_SWAV, 'cswt_plot_swav', &
       comment_add='Inconsistent number of parameters in command line')
  end if

  ! Initialise specified wavelet.
  if(wavelet_type == trim(WAV_TYPE_MEXHAT)) then

     if(n_param == 0) then
        swav_mother = cswt_swav_init(cswt_tmpl_s2_mexhat, nside, &
             name=trim(wavelet_type), fun_type=S2_SKY_FUN_TYPE_SPHERE)
     else
        swav_mother = cswt_swav_init(cswt_tmpl_s2_mexhat, nside, &
             param=param, name=trim(wavelet_type), &
             fun_type=S2_SKY_FUN_TYPE_SPHERE)
     end if
     
  elseif(wavelet_type == trim(WAV_TYPE_MORLET)) then
     
     if(n_param == 0) then
        swav_mother = cswt_swav_init(cswt_tmpl_s2_morlet, nside, &
             name=trim(wavelet_type), fun_type=S2_SKY_FUN_TYPE_SPHERE)
     else
        swav_mother = cswt_swav_init(cswt_tmpl_s2_morlet, nside, &
             param=param, name=trim(wavelet_type), &
             fun_type=S2_SKY_FUN_TYPE_SPHERE)
     end if
     
  elseif(wavelet_type == trim(WAV_TYPE_BUTTERFLY)) then
     
     if(n_param == 0) then
        swav_mother = cswt_swav_init(cswt_tmpl_s2_butterfly, nside, &
             name=trim(wavelet_type), fun_type=S2_SKY_FUN_TYPE_SPHERE)
     else
        swav_mother = cswt_swav_init(cswt_tmpl_s2_butterfly, nside, &
             param=param, name=trim(wavelet_type), &
             fun_type=S2_SKY_FUN_TYPE_SPHERE)
     end if
     
  else
     
     call cswt_error(CSWT_ERROR_ANALYSIS, 'cswt_analysis', &
          comment_add='Invalid wavelet type')
 
  end if
  
  ! Read dilations.
  call cswt_tr_io_txt_dilation_read(filename_dilation, dilation)

  ! For each dilation, perform dilation and write dilated wavelet to file.
  do i_dil = 1,size(dilation,1)

     swav = cswt_swav_init(swav_mother)

     call cswt_swav_dilate(swav, dilation(i_dil,:), norm_preserve)
     if(center_on_equator) then
        call cswt_swav_rotate(swav, ALPHA_CENTER, BETA_CENTER, GAMMA_CENTER)
     end if

     write(filename_out, '(a,a,i2.2,a)') trim(filename_out_prefix), &
       '_id', i_dil, '.fits'
     call cswt_swav_write_map_file(swav, filename_out)

     call cswt_swav_free(swav)

  end do

  ! Free data.
  call cswt_swav_free(swav_mother)
  deallocate(dilation)
  if(allocated(param)) deallocate(param)


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
            write(*,'(a)') 'Usage: cswt_plot_swav [-dil filename_dilation]'
            write(*,'(a)') '                      [-out filename_out_prefix]' 
            write(*,'(a)') '                      [-wav wavelet_type]' 
            write(*,'(a)') '                      [-nside nside]' 
            write(*,'(a)') '                      [-cen center_on_equator]' 
            write(*,'(a)') '                      [-norm norm_preserve]' 
            write(*,'(a)') '                      [-npar n_param]' 
            write(*,'(a)') &
              '                      [-par param(i) (repeat n_param times)]' 
            stop
          
          case ('-dil')
            filename_dilation = trim(arg)
         
          case ('-out')
            filename_out_prefix = trim(arg)

          case ('-wav')
            wavelet_type = trim(arg)

          case ('-nside')
            read(arg,*) nside

          case ('-cen')
            read(arg,*) center_on_equator

          case ('-norm')
            read(arg,*) norm_preserve

          case ('-npar')
            read(arg,*) n_param
            ! Allocate space for parameters.
            allocate(param(n_param), stat=fail)
            if(fail /= 0) then
               call cswt_error(CSWT_ERROR_MEM_ALLOC_FAIL, 'cswt_analysis')
            end if

          case ('-par')
            ! Parameters read in order appear in command line.
            i_param = i_param + 1
            read(arg,*) param(i_param)  

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program cswt_plot_swav
