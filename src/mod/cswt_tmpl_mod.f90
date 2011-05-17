!------------------------------------------------------------------------------
! cswt_tmpl_mod -- CSWT library template class
!
!! Contains definitions of template functions defined on both the sky and the
!! plane used to initialise spherical wavelets.  The functions defined on the
!! plane are projected onto the sphere numerically, whereas those functions 
!! defined on the sphere have already been analytically projected from
!! the plane.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.2 - November 2004
!
! Revisions:
!   November 2004 - Written by Jason McEwen
!------------------------------------------------------------------------------

module cswt_tmpl_mod

  use s2_types_mod
  use cswt_error_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    cswt_tmpl_s2_costhetaon2, &
    cswt_tmpl_s2_gaussian, &
    cswt_tmpl_s2_butterfly, &
    cswt_tmpl_s2_mexhat, &
    cswt_tmpl_s2_morlet, &
    cswt_tmpl_pln_gaussian, &
    cswt_tmpl_pln_butterfly, &
    cswt_tmpl_pln_mexhat, &
    cswt_tmpl_pln_morlet


  !----------------------------------------------------------------------------

  contains

    
    !--------------------------------------------------------------------------
    ! Functions analytically defined directly on the sphere
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    ! cswt_tmpl_s2_costhetaon2
    !
    !! Template function defined on the sphere.  
    !!   f(theta,phi) = cos(theta/2.0e0)
    !!
    !! Notes:
    !!   - Param array not used.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_s2_costhetaon2(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val
      
      val = cos(theta/2.0e0)

    end function cswt_tmpl_s2_costhetaon2


    !--------------------------------------------------------------------------
    ! cswt_tmpl_s2_gaussian
    !
    !! Template function defined on the sphere.  
    !! 2d Gaussian.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as eccentricity.
    !!   - Two parameters in param array: taken as standard deviations.
    !!   - Error occurs if param array has length greater than 2.
    !!   - If no parameter array is given then sigma_x=1 and sigma_y=1.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_s2_gaussian(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: r 

      r = 2e0 * tan(theta/2.0e0)

      val = 2.0e0 / (1.0e0 + cos(theta)) * cswt_tmpl_pln_gaussian(r, phi, param)

    end function cswt_tmpl_s2_gaussian


    !--------------------------------------------------------------------------
    ! cswt_tmpl_s2_butterfly
    !
    !! Template function defined on the sphere.  
    !! Butterfly is a Gaussian in y direction and first derrivative of 
    !! Gaussian in x direction.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as eccentricity.
    !!   - Two parameters in param array: taken as standard deviations.
    !!   - Error occurs if param array has length greater than 2.
    !!   - If no parameter array is given then sigma_x=1 and sigma_y=1.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_s2_butterfly(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: r

      r = 2e0 * tan(theta/2.0e0)

      val = 2.0e0 / (1.0e0 + cos(theta)) * cswt_tmpl_pln_butterfly(r, phi, param)

    end function cswt_tmpl_s2_butterfly


    !--------------------------------------------------------------------------
    ! cswt_tmpl_s2_mexhat
    !
    !! Template function defined on the sphere.  
    !! Mexican hat is negative of Laplacian of 2D Gaussian.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as eccentricity.
    !!   - Two parameters in param array: taken as standard deviations.
    !!   - Error occurs if param array has length greater than 2.
    !!   - If no parameter array is given then sigma_x=1 and sigma_y=1.
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_s2_mexhat(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: r

      r = 2e0 * tan(theta/2.0e0)

      val = 2.0e0 / (1.0e0 + cos(theta)) * cswt_tmpl_pln_mexhat(r, phi, param)

    end function cswt_tmpl_s2_mexhat


    !--------------------------------------------------------------------------
    ! cswt_tmpl_s2_morlet
    !
    !! Template function defined on the sphere.  
    !! Real Morlet wavelet.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as k0 of wave vector, where
    !!     k=(k0,0).
    !!   - Error occurs if param array has length greater than 1.
    !!   - If no parameter array is given then k0=10
    !!
    !! Variables:
    !!   - theta: Theta spherical coordiante in range [0,pi].
    !!   - phi: Phi spherical coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_s2_morlet(theta, phi, param) result(val)

      real(s2_sp), intent(in) :: theta, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: r

      r = 2e0 * tan(theta/2.0e0)

      val = 2.0e0 / (1.0e0 + cos(theta)) * cswt_tmpl_pln_morlet(r, phi, param)

    end function cswt_tmpl_s2_morlet


    !--------------------------------------------------------------------------
    ! Functions analytically defined on the plane
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! cswt_tmpl_pln_gaussian
    !
    !! Template function defined on the plane.  Later projected onto sphere
    !! by numerical stereographic projection.  
    !! 2d Gaussian.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as eccentricity.
    !!   - Two parameters in param array: taken as standard deviations.
    !!   - Error occurs if param array has length greater than 2.
    !!   - If no parameter array is given then sigma_x=1 and sigma_y=1.
    !!
    !! Variables:
    !!   - r: R polar coordinate in 2d.
    !!   - phi: Phi polar coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_pln_gaussian(r, phi, param) result(val)

      real(s2_sp),intent(in) :: r, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: sigma_x = 1.0e0
      real(s2_sp) :: sigma_y = 1.0e0
      real(s2_sp) :: x, y

      if(present(param)) then
 
         ! If only one parameter then taken as eccentricity.
         if(size(param) == 1) then

            if(param(1) >= 1.0e0 .or. param(1) < 0.0e0) then
               call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
                 'cswt_tmpl_pln_gaussian', &
                 comment_add='Eccentricity not in range [0,1)')
            end if

            sigma_y = 1.0
            sigma_x = sigma_y * (1 - param(1)**2)**(1/4.0d0)

         ! If two parameters then taken as standard deviations.
         else if(size(param) == 2) then

            sigma_x = param(1)
            sigma_y = param(2)

         else

            call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
              'cswt_tmpl_pln_gaussian')

         end if

      end if

      x = r * cos(phi)
      y = r * sin(phi)

      ! Normalised.
!      val = 1 / (2.0e0 * pi * sigma_x * sigma_y) &
!            * exp(- (x**2/sigma_x**2 + y**2/sigma_y**2) / 2.0e0)

      ! Unit amplitude.
      val = exp(- (x**2/sigma_x**2 + y**2/sigma_y**2) / 2.0e0)
      
    end function cswt_tmpl_pln_gaussian


    !--------------------------------------------------------------------------
    ! cswt_tmpl_pln_butterfly
    !
    !! Template function defined on the plane.  Later projected onto sphere
    !! by numerical stereographic projection.  
    !! Butterfly is a Gaussian in y direction and first derrivative of 
    !! Gaussian in x direction.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as eccentricity.
    !!   - Two parameters in param array: taken as standard deviations.
    !!   - Error occurs if param array has length greater than 2.
    !!   - If no parameter array is given then sigma_x=1 and sigma_y=1.
    !!
    !! Variables:
    !!   - r: R polar coordinate in 2d.
    !!   - phi: Phi polar coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_pln_butterfly(r, phi, param) result(val)

      real(s2_sp),intent(in) :: r, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: sigma_x = 1.0e0
      real(s2_sp) :: sigma_y = 1.0e0
      real(s2_sp) :: x, y

      if(present(param)) then
 
         ! If only one parameter then taken as eccentricity.
         if(size(param) == 1) then

            if(param(1) >= 1.0e0 .or. param(1) < 0.0e0) then
               call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
                 'cswt_tmpl_pln_butterfly', &
                 comment_add='Eccentricity not in range [0,1)')
            end if

            sigma_y = 1.0
            sigma_x = sigma_y * (1 - param(1)**2)**(1/4.0d0)

         ! If two parameters then taken as standard deviations.
         else if(size(param) == 2) then

            sigma_x = param(1)
            sigma_y = param(2)

         else

            call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
              'cswt_tmpl_pln_butterfly')

         end if

      end if

      x = r * cos(phi)
      y = r * sin(phi)

      ! Normalised.
!      val = -x / (2.0e0 * pi * sigma_x**3.0e0 * sigma_y) &
!            * exp(- (x**2/sigma_x**2 + y**2/sigma_y**2) / 2.0e0)

      ! Unscaled amplitude.
      val = -x * exp(- (x**2/sigma_x**2 + y**2/sigma_y**2) / 2.0e0)
            
    end function cswt_tmpl_pln_butterfly


    !--------------------------------------------------------------------------
    ! cswt_tmpl_pln_mexhat
    !
    !! Template function defined on the plane.  Later projected onto sphere
    !! by numerical stereographic projection.  
    !! Mexican hat is negative of Laplacian of 2D Gaussian.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as eccentricity.
    !!   - Two parameters in param array: taken as standard deviations.
    !!   - Error occurs if param array has length greater than 2.
    !!   - If no parameter array is given then sigma_x=1 and sigma_y=1.
    !!
    !! Variables:
    !!   - r: R polar coordinate in 2d.
    !!   - phi: Phi polar coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_pln_mexhat(r, phi, param) result(val)

      real(s2_sp),intent(in) :: r, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: sigma_x = 1.0e0
      real(s2_sp) :: sigma_y = 1.0e0
      real(s2_sp) :: x, y

      if(present(param)) then
 
         ! If only one parameter then taken as eccentricity.
         if(size(param) == 1) then

            if(param(1) >= 1.0e0 .or. param(1) < 0.0e0) then
               call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
                 'cswt_tmpl_pln_mexhat', &
                 comment_add='Eccentricity not in range [0,1)')
            end if

            sigma_y = 1.0
            sigma_x = sigma_y * (1 - param(1)**2)**(1/4.0d0)

         ! If two parameters then taken as standard deviations.
         else if(size(param) == 2) then

            sigma_x = param(1)
            sigma_y = param(2)

         else

            call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
              'cswt_tmpl_pln_mexhat')

         end if

      end if

      x = r * cos(phi)
      y = r * sin(phi)

      ! Normalised.
!      val = 1 / (2.0e0 * pi * sigma_x**3 * sigma_y**3) &
!        * ( sigma_x**2 + sigma_y**2 &
!            - x**2/(sigma_x/sigma_y)**2 -y**2/(sigma_y/sigma_x)**2 ) &
!        * exp(- (x**2/sigma_x**2 + y**2/sigma_y**2) / 2.0e0)

      ! Unit amplitude.
      val = 1e0 / (sigma_x**2 + sigma_y**2) &
        * ( sigma_x**2 + sigma_y**2 &
            - x**2/(sigma_x/sigma_y)**2 -y**2/(sigma_y/sigma_x)**2 ) &
        * exp(- (x**2/sigma_x**2 + y**2/sigma_y**2) / 2.0e0)

    end function cswt_tmpl_pln_mexhat

    
    !--------------------------------------------------------------------------
    ! cswt_tmpl_pln_morlet
    !
    !! Template function defined on the plane.  Later projected onto sphere
    !! by numerical stereographic projection.  
    !! Real Morlet wavelet.
    !!
    !! Notes:
    !!   - One parameter in param array: taken as k0 of wave vector, where
    !!     k=(k0,0).
    !!   - Error occurs if param array has length greater than 1.
    !!   - If no parameter array is given then k0=10
    !!
    !! Variables:
    !!   - r: R polar coordinate in 2d.
    !!   - phi: Phi polar coordinate in range [0,2pi).
    !!   - [param]: Optional parameter array.
    !
    !! @author J. D. McEwen
    !! @version 0.2 - November 2004
    !
    ! Revisions:
    !   November 2004 - Written by Jason McEwen
    !--------------------------------------------------------------------------
   
    function cswt_tmpl_pln_morlet(r, phi, param) result(val)

      real(s2_sp),intent(in) :: r, phi
      real(s2_sp), intent(in), optional :: param(:)
      real(s2_sp) :: val

      real(s2_sp) :: k(2) = (/10.0e0, 0.0e0/)
      real(s2_sp) :: r_scale, x, y, angle

      if(present(param)) then
 
         ! If one parameter then taken as k0, where k=(k0,0).
         if(size(param) == 1) then

            k(1) = param(1)
            k(2) = 0.0e0

         else

            call cswt_error(CSWT_ERROR_TMPL_PARAM_INVALID, &
              'cswt_tmpl_pln_morlet')

         end if

      end if

      ! Radius scaled so effective size on the sky is the same as the Mexican
      ! Hat wavelet of equivalent dilation.
      r_scale = r / sqrt(2.0d0)

      x = r_scale * cos(phi)
      y = r_scale * sin(phi)      

      angle = k(1)*x + k(2)*y

      val = real(cmplx(cos(angle), sin(angle)) * exp(-(r_scale**2)), &
            kind=s2_sp)

    end function cswt_tmpl_pln_morlet


end module cswt_tmpl_mod
