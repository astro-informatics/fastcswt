!------------------------------------------------------------------------------
! cswt_analysis -- CSWT analysis program
!
!! Perform fast directional continuous spherical wavelet transform.
!!
!! Usage: cswt_analysis 
!!   - [-help]: Display usage information.
!!   - [-dil filename_dilation]: Name of dilation file containing dilations
!!     to consider.
!!   - [-inp filename_data]: Name of data file containing sky to take 
!!     transform of.
!!   - [-ext file_extension]: Extension of filename_data to read sky from.
!!   - [-out filename_out]: Name of output wavelet coefficient file.
!!   - [-method method_str]: String specifying cswt algorithm to use.
!!   - [-message message]: Logical specifying whether to display messages.
!!   - [-filename_isotropic filename_isotropic]: Optional name of fits output
!!     file of wavelet coefficients in healpix format for isotropic algorithm.
!!   - [-wav wavelet_type]: String specifying wavelet type.
!!   - [-wavfile filename_wavelet]: Name of wavelet file if wavelet is to be
!!     read from file (i.e. -wav file)
!!   - [-lmax lmax]: lmax of wavelet transform.
!!   - [-mmax mmax]: mmax of wavelet transform.
!!   - [-ngamma n_gamma]: Number of orientations to consider in transform.
!!   - [-time time_it]: Logical specifying whether to time algorithm.
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

program cswt_analysis

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
    filename_in, &
    filename_out, &
    filename_wavelet, &
    filename_isotropic, &
    wavelet_type
  integer :: i_param = 0, n_param = 0, n_gamma = 1, lmax = 0, mmax = 0
  integer :: nside, fail
  integer :: file_extension = 1
  real(s2_sp), allocatable :: param(:)

  integer :: method
  character(len=S2_STRING_LEN), parameter ::  &
    METHOD_FFT_REAL = 'fft_real', &
    METHOD_FFT = 'fft', &
    METHOD_DFT = 'dft', &
    METHOD_ISOTROPIC = 'isotropic', &
    METHOD_DIRECT = 'direct'
  character(len=S2_STRING_LEN) :: method_str = METHOD_FFT_REAL
  logical :: message = .true.

  type(s2_sky) :: data
  type(cswt_swav) :: swav_mother
  type(cswt_tr) :: tr

  logical :: time_it = .false.
  logical :: isotropic_out = .false.

  ! Set default parameters.
  wavelet_type = WAV_TYPE_MEXHAT
!**
! For ISW work use this lmax.  
! Set as default status here so don't have to continually enter 
! at command line (where mistake could occur).
  lmax = 128
  mmax = 128
!**

  ! Parse options from command line.
  call parse_options()

  ! Check correct number of parameters were read.
  if(i_param /= n_param) then
     call cswt_error(CSWT_ERROR_ANALYSIS, 'cswt_analysis', &
       comment_add='Inconsistent number of parameters in command line')
  end if

  ! Read in data sky.
  data = s2_sky_init(filename_in, file_extension)
  nside = s2_sky_get_nside(data)

  ! Set up tr data structure.
  if(wavelet_type == trim(WAV_TYPE_FILE)) then

     ! If loading mother wavelet from file construct tr structure directly.
     tr = cswt_tr_init(filename_wavelet, filename_dilation, lmax, mmax, &
       n_gamma, trim(wavelet_type))

  else 

     ! Otherwise construct mother wavelet analytically before constructing tr.

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

     ! Construct tr from mother wavelet.
     tr = cswt_tr_init(swav_mother, filename_dilation, lmax, mmax, n_gamma)
     call cswt_swav_free(swav_mother)

  end if

  ! Set method.
   select case(trim(method_str))

     case(trim(METHOD_FFT_REAL))
        method = CSWT_TR_METHOD_FAST_FFT_REAL

     case(trim(METHOD_FFT))
        method = CSWT_TR_METHOD_FAST_FFT

     case(trim(METHOD_DFT))
        method = CSWT_TR_METHOD_FAST_DFT
        
     case(trim(METHOD_ISOTROPIC))
        method = CSWT_TR_METHOD_FAST_ISOTROPIC

     case(trim(METHOD_DIRECT))
        method = CSWT_TR_METHOD_DIRECT

     case default
        call cswt_error(CSWT_ERROR_TR_METHOD_INVALID, 'cswt_analysis')

  end select

  ! Perform analysis.
  if(isotropic_out) then
     call cswt_tr_analysis(tr, data, method, message, time_it=time_it, &
       filename_isotropic=filename_isotropic)
  else
     call cswt_tr_analysis(tr, data, method, message, time_it=time_it)
  end if

  ! Write wavelet coefficients.
  call cswt_tr_io_fits_write_wcoeff(filename_out, tr)

  ! Free data.
  call s2_sky_free(data)
  call cswt_tr_free(tr)
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
            write(*,'(a)') 'Usage: cswt_analysis [-dil filename_dilation]'
            write(*,'(a)') '                     [-inp filename_data]'
            write(*,'(a)') '                     [-ext file_extension]' 
            write(*,'(a)') '                     [-out filename_out]' 
            write(*,'(a)') '                     [-method method_str]' 
            write(*,'(a)') '                     [-message message]' 
            write(*,'(a,a)') '                     ', &
              '[-filename_isotropic filename_isotropic]' 
            write(*,'(a)') '                     [-wav wavelet_type]' 
            write(*,'(a)') '                     [-wavfile filename_wavelet]' 
            write(*,'(a)') '                     [-lmax lmax]' 
            write(*,'(a)') '                     [-mmax mmax]' 
            write(*,'(a)') '                     [-ngamma n_gamma]' 
            write(*,'(a)') '                     [-time time_it]' 
            write(*,'(a)') '                     [-npar n_param]' 
            write(*,'(a)') &
              '                     [-par param(i) (repeat n_param times)]' 
            stop

          case ('-method')
            method_str = trim(arg)

          case ('-message')
             read(arg,*) message

         case ('-filename_isotropic')
            filename_isotropic = trim(arg)
            isotropic_out = .true.
          
          case ('-dil')
            filename_dilation = trim(arg)
         
          case ('-inp')
            filename_in = trim(arg)

          case ('-ext')
            read(arg,*) file_extension

          case ('-out')
            filename_out = trim(arg)

          case ('-wav')
            wavelet_type = trim(arg)

          case ('-wavfile')
            filename_wavelet = trim(arg)

          case ('-lmax')
            read(arg,*) lmax

          case ('-mmax')
            read(arg,*) mmax

          case ('-ngamma')
            read(arg,*) n_gamma

          case ('-time')
            read(arg,*) time_it

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


end program cswt_analysis
