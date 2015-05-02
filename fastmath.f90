!f90

!This fortran module contains subroutines for calculating the intensity at each point.
!There are subroutines both with and without both symmetry and the small-angle approximation.
!The first number in the title indicates if it's symmetric (1=symmetric,0=asymmetric)
!The second number in the title indicates if the small-angle approxmation is used.

!OpenMP is used in the outermost loop of each as the problem is embarassingly parallel.

!Copyright ANU/Scott Medling, 2015.

module fastmath
contains

subroutine sumintensity00(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*8, dimension(npts,4), intent(in) :: points
   real*8, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity,temp_intensity_2,total_intensity,QdotR
   real*8, dimension(3) :: Q
   integer*4 :: p
   !'asymmetry'; no small angle approximation
   total_intensity = 0
   !$OMP PARALLEL DO PRIVATE(Q,QdotR,temp_intensity,temp_intensity_2) SHARED(mask,intensity) REDUCTION(+:total_intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(x_pixels*2*ehc))**2 /)
                !!POSSIBLE FORMULA ERROR, CHECK X_PIXELS IN DENOMINATOR!!
            !intensity(i,j)= SUM(density(p)*COS(DOT_PRODUCT(Q,R(p))))**2 + SUM(density(p)*SIN(DOT_PRODUCT(Q,R(p))))**2
            temp_intensity = 0
            temp_intensity_2 = 0
            do p=1,npts
               QdotR = DOT_PRODUCT(Q,points(p,1:3))
               temp_intensity = temp_intensity + points(p,4)*COS(QdotR)
               temp_intensity_2 = temp_intensity_2 + points(p,4)*SIN(QdotR)
            end do
            intensity(i,j) = temp_intensity**2 + temp_intensity_2**2
            total_intensity = total_intensity + intensity(i,j)
            !total_intensity = total_intensity + temp_intensity**2 + temp_intensity_2**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   intensity = intensity / total_intensity
   return
end subroutine sumintensity00

subroutine sumintensity01(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*8, dimension(npts,4), intent(in) :: points
   real*8, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity,temp_intensity_2,total_intensity
   real*8, dimension(3) :: Q,QP
   integer*4 :: p
   !'asymmetry'; small angle approximation
   total_intensity = 0
   !$OMP PARALLEL DO PRIVATE(Q,temp_intensity,temp_intensity_2) SHARED(mask,intensity) REDUCTION(+:total_intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, 0.0d0 /)
            QP = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                 2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(x_pixels*2*ehc))**2 /)
                !!POSSIBLE FORMULA ERROR, CHECK X_PIXELS IN DENOMINATOR!!
            !intensity(i,j)= SUM(density(p)*COS(DOT_PRODUCT(Q,R(p))))**2 + SUM(density(p)*SIN(DOT_PRODUCT(QP,R(p))))**2
            temp_intensity = 0
            temp_intensity_2 = 0
            do p=1,npts
               temp_intensity = temp_intensity + (points(p,4)*COS(DOT_PRODUCT(Q,points(p,1:3))))
               temp_intensity_2 = temp_intensity_2 + (points(p,4)*SIN(DOT_PRODUCT(QP,points(p,1:3))))
            end do
            intensity(i,j) = temp_intensity**2 + temp_intensity_2**2
            total_intensity = total_intensity + intensity(i,j)
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   intensity = intensity / total_intensity
   return
end subroutine sumintensity01

subroutine sumintensity10(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*8, dimension(npts,4), intent(in) :: points
   real*8, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity,total_intensity
   real*8, dimension(3) :: Q
   integer*4 :: p
   !'symmetry'; no small angle approximation
   total_intensity = 0
   !$OMP PARALLEL DO PRIVATE(Q,temp_intensity,temp_intensity_2) SHARED(mask,intensity) REDUCTION(+:total_intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(x_pixels*2*ehc))**2 /)
                !!POSSIBLE FORMULA ERROR, CHECK X_PIXELS IN DENOMINATOR!!
            temp_intensity = 0
            do p=1,npts
               !temp_intensity = temp_intensity + (density(p)*COS(DOT_PRODUCT(Q,R)))**2
               temp_intensity = temp_intensity + (points(p,4)*COS(DOT_PRODUCT(Q,points(p,1:3))))
            end do
            intensity(i,j) = temp_intensity**2
            total_intensity = total_intensity + temp_intensity**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   intensity = intensity / total_intensity
   return
end subroutine sumintensity10

subroutine sumintensity11(qsize,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   integer*4, intent(in) :: npts
   real*8, dimension(npts,4), intent(in) :: points
   real*8, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity,total_intensity
   real*8, dimension(3) :: Q
   integer*4 :: p
   !'symmetry'; small angle approximation
   total_intensity = 0
   !$OMP PARALLEL DO PRIVATE(Q,temp_intensity,temp_intensity_2) SHARED(mask,intensity) REDUCTION(+:total_intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, 0.0d0 /)
            temp_intensity = 0   !or 0.0d0?
            do p=1,npts
               !temp_intensity = temp_intensity + (density(p)*COS(DOT_PRODUCT(Q,R)))
               temp_intensity = temp_intensity + (points(p,4)*COS(DOT_PRODUCT(Q,points(p,1:3))))
            end do
            intensity(i,j) = temp_intensity**2
            total_intensity = total_intensity + temp_intensity**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   intensity = intensity / total_intensity
   return
end subroutine sumintensity11

end module fastmath
