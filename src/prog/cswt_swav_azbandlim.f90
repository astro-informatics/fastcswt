!------------------------------------------------------------------------------
! cswt_swav_azbandlim -- CSWT swav_azbandlim program
!
!! Determine azimuthal band limit of wavelet.
!!
!! Usage: cswt_swav_azbandlim
!!   - [-help]: Display usage information.
!!   - [-wav wavelet_type]: String specifying wavelet type.
!!   - [-nside nside]: healpix nside of wavelet to comsider.
!!   - [-cutoff cutoff_prop]: Cutoff proporiton for energy of cms.
!!   - [-lmax lmax]: Healpix lmax.
!!   - [-mmax mmax]: Healpix mmax.
!!   - [-dil1 dil1]: Component a of dilation to analyse wavelet at.
!!   - [-dil2 dil2]: Component b of dilation to anaylse wavelet at.
!!   - [-norm norm_preserve]: Logical specifying whether to perform a
!!      norm preserving dilation.
!!   - [-npar n_param]: Number of wavelet parameters specified.
!!   - [-par param(i) (repeat n_param times)]: Wavelet parameters.  
!!     Parameters are read sequentially on order appear in command line.    
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - May 2005
!
! Revisions:
!   May 2005 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program cswt_swav_azbandlim

  use s2_types_mod
  use s2_sky_mod
  use cswt_error_mod
  use cswt_tmpl_mod
  use cswt_swav_mod

  implicit none

  character(len=S2_STRING_LEN), parameter ::  &
    WAV_TYPE_FILE = 'file', &
    WAV_TYPE_MEXHAT = 'mexhat', &
    WAV_TYPE_BUTTERFLY = 'butterfly', &
    WAV_TYPE_MORLET = 'morlet'
  character(len=S2_STRING_LEN) :: wavelet_type

  logical :: norm_preserve
  integer :: i_param = 0, n_param = 0
  integer :: nside, lmax, mmax, fail, mmax_min
  real(s2_sp), allocatable :: param(:)
  real(s2_sp) :: dilation(1,2) = 1.0e0
  real(s2_sp) :: dil1 = 1.0e0, dil2 = 1.0e0
  real(s2_sp) :: cutoff_prop
  type(cswt_swav) :: swav

  ! Set default parameters.
  wavelet_type = WAV_TYPE_MEXHAT
  nside = 128
  lmax = 256
  mmax = 256
  cutoff_prop = 0.95e0
  norm_preserve = .true.

  ! Parse options from command line.
  call parse_options()

  dilation(1,1) = dil1
  dilation(1,2) = dil2

  ! Check correct number of parameters were read.
  if(i_param /= n_param) then
     call cswt_error(CSWT_ERROR_PLOT_SWAV, 'cswt_plot_swav', &
       comment_add='Inconsistent number of parameters in command line')
  end if

  ! Initialise specified wavelet.
  if(wavelet_type == trim(WAV_TYPE_MEXHAT)) then

     if(n_param == 0) then
        swav = cswt_swav_init(cswt_tmpl_s2_mexhat, nside, &
             name=trim(wavelet_type), fun_type=S2_SKY_FUN_TYPE_SPHERE)
     else
        swav = cswt_swav_init(cswt_tmpl_s2_mexhat, nside, &
             param=param, name=trim(wavelet_type), &
             fun_type=S2_SKY_FUN_TYPE_SPHERE)
     end if
     
  elseif(wavelet_type == trim(WAV_TYPE_MORLET)) then
     
     if(n_param == 0) then
        swav = cswt_swav_init(cswt_tmpl_s2_morlet, nside, &
             name=trim(wavelet_type), fun_type=S2_SKY_FUN_TYPE_SPHERE)
     else
        swav = cswt_swav_init(cswt_tmpl_s2_morlet, nside, &
             param=param, name=trim(wavelet_type), &
             fun_type=S2_SKY_FUN_TYPE_SPHERE)
     end if
     
  elseif(wavelet_type == trim(WAV_TYPE_BUTTERFLY)) then
     
     if(n_param == 0) then
        swav = cswt_swav_init(cswt_tmpl_s2_butterfly, nside, &
             name=trim(wavelet_type), fun_type=S2_SKY_FUN_TYPE_SPHERE)
     else
        swav = cswt_swav_init(cswt_tmpl_s2_butterfly, nside, &
             param=param, name=trim(wavelet_type), &
             fun_type=S2_SKY_FUN_TYPE_SPHERE)
     end if
     
  else
     
     call cswt_error(CSWT_ERROR_ANALYSIS, 'cswt_analysis', &
          comment_add='Invalid wavelet type')
 
  end if
  
  ! Dilate if required
  call cswt_swav_dilate(swav, dilation, norm_preserve)

  ! Compute alms.
  call cswt_swav_compute_alm(swav, lmax, mmax)

  ! Compute wavelet azimuthal band limit.
  mmax_min = cswt_swav_azimuthal_bl(swav, cutoff_prop)

  write(*,*) 'Azimuthal band limit: ', mmax_min

  ! Free data.
  call cswt_swav_free(swav)
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
            write(*,'(a,a)') 'Usage: cswt_swav_azimuthal_bl ', &
                 '[-wav wavelet_type]'
            write(*,'(a,a)') '                              ', &
                 '[-nside nside]' 
            write(*,'(a,a)') '                              ', &
                  '[-cutoff cutoff_prop]' 
            write(*,'(a,a)') '                              ', &
                  '[-lmax lmax]' 
            write(*,'(a,a)') '                              ', &
                 '[-mmax mmax]' 
            write(*,'(a,a)') '                              ', &
                 '[-dil1 dil1]' 
            write(*,'(a,a)') '                              ', &
                 '[-dil2 dil2]' 
            write(*,'(a,a)') '                              ', &
                 '[-norm norm_preserve]' 
            write(*,'(a,a)') '                              ', &
                 '[-npar n_param]' 
            write(*,'(a,a)') '                              ', &
                 '[-par param(i) (repeat n_param times)]' 
            stop
          
          case ('-wav')
            wavelet_type = trim(arg)

          case ('-nside')
            read(arg,*) nside

          case ('-cutoff')
            read(arg,*) cutoff_prop

          case ('-lmax')
            read(arg,*) lmax

          case ('-mmax')
            read(arg,*) mmax

          case ('-dil1')
            read(arg,*) dil1

          case ('-dil2')
            read(arg,*) dil2

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


end program cswt_swav_azbandlim
