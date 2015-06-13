!f90

!The sumint module contains subroutines for calculating the intensity at each point.
!There are subroutines both with and without both symmetry and the small-angle approximation.
!The first number in the title indicates if it's symmetric (1=symmetric,0=asymmetric)
!The second number in the title indicates if the small-angle approxmation is used.
!Single-precision is used for the input array to keep the size down. Output array is double.

!The density module contains vectorized subroutines to calculate the density at every point in real space.

!OpenMP is used in the outermost loop of each as the problem is embarassingly parallel.

!Copyright ANU/Scott Medling, 2015.

!TODO:
!Tell OpenMP to synchronize access to shared variables (points) if possible.
!!$OMP PARALLEL DO ... COLLAPSE(2) !tells it to collapse the double do loop
!If running with more than 100,000 points (>3MB) and grid no larger than 100x100, it may be faster to move the points loop from innermost to outermost.  Note, this will require storing a temporary Q(3,x,y) array in order to still avoid calculating Q 100,000 times.
!Divide up points into smaller chuncks, make temp_intensity into arrays so can sqaure at end.
!Change mask implmentation so x,y coordinates of relevant pixels are passed in so there's no need for an if statement.


module sumint
contains

subroutine sumintensity00(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*4, dimension(4,npts), intent(in) :: points
   real*4, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8, dimension(3) :: Q
   real*8 :: temp_intensity,temp_intensity_2,QdotR
   integer*4 :: p
   !'asymmetry'; no small angle approximation
   !$OMP PARALLEL DO PRIVATE(Q,QdotR,temp_intensity,temp_intensity_2) SHARED(mask,points,intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(y_pixels*2*ehc))**2 /)
                !!TODO: POSSIBLE FORMULA ERROR, CHECK PIXELS IN DENOMINATOR!!
            !intensity(i,j)= SUM(density(p)*COS(DOT_PRODUCT(Q,R(p))))**2 + SUM(density(p)*SIN(DOT_PRODUCT(Q,R(p))))**2
            temp_intensity = 0
            temp_intensity_2 = 0
            do p=1,npts
               QdotR = DOT_PRODUCT(Q,points(1:3,p))
               temp_intensity = temp_intensity + points(4,p)*COS(QdotR)
               temp_intensity_2 = temp_intensity_2 + points(4,p)*SIN(QdotR)
            end do
            intensity(i,j) = temp_intensity**2 + temp_intensity_2**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   return
end subroutine sumintensity00

subroutine test00(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*4, dimension(4,npts), intent(in) :: points
   real*4, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8, dimension(3) :: Q
   real*8 :: QdotR
   real*8, dimension(x_pixels,y_pixels) :: temp_intensity,temp_intensity_2
   integer*4 :: p
   !'asymmetry'; no small angle approximation
   temp_intensity(:,:) = 0
   temp_intensity_2(:,:) = 0
   !$OMP PARALLEL DO PRIVATE(Q,QdotR) SHARED(mask,points) REDUCTION(+:temp_intensity,temp_intensity_2)
   do p=1,npts
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(y_pixels*2*ehc))**2 /)
                !!TODO: POSSIBLE FORMULA ERROR, CHECK PIXELS IN DENOMINATOR!!
            !intensity(i,j)= SUM(density(p)*COS(DOT_PRODUCT(Q,R(p))))**2 + SUM(density(p)*SIN(DOT_PRODUCT(Q,R(p))))**2
            QdotR = DOT_PRODUCT(Q,points(1:3,p))
            temp_intensity(i,j) = temp_intensity(i,j) + points(4,p)*COS(QdotR)
            temp_intensity_2(i,j) = temp_intensity_2(i,j) + points(4,p)*SIN(QdotR)
         end if
      end do
   end do
   end do
   !$OMP END PARALLEL DO
   intensity(:,:) = temp_intensity(:,:)**2 + temp_intensity_2(:,:)**2
   return
end subroutine test00

!Computationally, there is no reason to ever use this; it takes (slightly) longer than sumint00.
subroutine sumintensity01(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*4, dimension(4,npts), intent(in) :: points
   real*4, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity,temp_intensity_2
   real*8, dimension(3) :: Q,QP
   integer*4 :: p
   !'asymmetry'; small angle approximation
   !$OMP PARALLEL DO PRIVATE(Q,temp_intensity,temp_intensity_2) SHARED(mask,points,intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, 0.0d0 /)
            QP = (/ Q(1), Q(2), &
                 2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(y_pixels*2*ehc))**2 /)
                !!TODO: POSSIBLE FORMULA ERROR, CHECK PIXELS IN DENOMINATOR!!
            !intensity(i,j)= SUM(density(p)*COS(DOT_PRODUCT(Q,R(p))))**2 + SUM(density(p)*SIN(DOT_PRODUCT(QP,R(p))))**2
            temp_intensity = 0
            temp_intensity_2 = 0
            do p=1,npts
               temp_intensity = temp_intensity + (points(4,p)*COS(DOT_PRODUCT(Q,points(1:3,p))))
               temp_intensity_2 = temp_intensity_2 + (points(4,p)*SIN(DOT_PRODUCT(QP,points(1:3,p))))
            end do
            intensity(i,j) = temp_intensity**2 + temp_intensity_2**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   return
end subroutine sumintensity01

subroutine sumintensity10(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*4, dimension(4,npts), intent(in) :: points
   real*4, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity
   real*8, dimension(3) :: Q
   integer*4 :: p
   !'symmetry'; no small angle approximation
   !$OMP PARALLEL DO PRIVATE(Q,temp_intensity) SHARED(mask,points,intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(y_pixels*2*ehc))**2 /)
                !!TODO: POSSIBLE FORMULA ERROR, CHECK PIXELS IN DENOMINATOR!!
            temp_intensity = 0
            do p=1,npts
               !temp_intensity = temp_intensity + (density(p)*COS(DOT_PRODUCT(Q,R)))**2
               temp_intensity = temp_intensity + (points(4,p)*COS(DOT_PRODUCT(Q,points(1:3,p))))
            end do
            intensity(i,j) = temp_intensity**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   return
end subroutine sumintensity10

subroutine sumintensity11(qsize,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   integer*4, intent(in) :: npts
   real*4, dimension(4,npts), intent(in) :: points
   real*4, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(out) :: intensity
   real*8 :: temp_intensity
   real*8, dimension(2) :: Q
   integer*4 :: p
   !'symmetry'; small angle approximation
   !$OMP PARALLEL DO PRIVATE(Q,temp_intensity) SHARED(mask,points,intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize /)
            temp_intensity = 0   !or 0.0d0?
            do p=1,npts
               !temp_intensity = temp_intensity + (density(p)*COS(DOT_PRODUCT(Q,R)))
               temp_intensity = temp_intensity + (points(4,p)*COS(DOT_PRODUCT(Q,points(1:2,p))))
            end do
            intensity(i,j) = temp_intensity**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   return
end subroutine sumintensity11


subroutine sumintensity00add(qsize,ehc,mask,x_pixels,y_pixels,points,npts,intensity)
   real*8, intent(in) :: qsize
   integer*4, intent(in) :: x_pixels,y_pixels
   real*8, intent(in) :: ehc
   integer*4, intent(in) :: npts
   real*4, dimension(4,npts), intent(in) :: points
   real*4, dimension(x_pixels,y_pixels), intent(in) :: mask
   real*8, dimension(x_pixels,y_pixels), intent(inout) :: intensity
   real*8, dimension(3) :: Q
   real*8 :: temp_intensity,temp_intensity_2,QdotR
   integer*4 :: p
   !'asymmetry'; no small angle approximation
   !$OMP PARALLEL DO PRIVATE(Q,QdotR,temp_intensity,temp_intensity_2) SHARED(mask,points,intensity)
   do j=1,y_pixels
      do i=1,x_pixels
         if (mask(i,j) > 0) then
            Q = (/ i*qsize/x_pixels-0.5*qsize, j*qsize/y_pixels-0.5*qsize, &
                2*ehc*sin(sqrt((i-0.5*x_pixels)**2+(j-0.5*y_pixels)**2)*qsize/(y_pixels*2*ehc))**2 /)
                !!TODO: POSSIBLE FORMULA ERROR, CHECK PIXELS IN DENOMINATOR!!
            !intensity(i,j)= SUM(density(p)*COS(DOT_PRODUCT(Q,R(p))))**2 + SUM(density(p)*SIN(DOT_PRODUCT(Q,R(p))))**2
            temp_intensity = 0
            temp_intensity_2 = 0
            do p=1,npts
               QdotR = DOT_PRODUCT(Q,points(1:3,p))
               temp_intensity = temp_intensity + points(4,p)*COS(QdotR)
               temp_intensity_2 = temp_intensity_2 + points(4,p)*SIN(QdotR)
            end do
            intensity(i,j) = intensity(i,j) + temp_intensity**2 + temp_intensity_2**2
         else
            intensity(i,j) = 0
         end if
      end do
   end do
   !$OMP END PARALLEL DO
   return
end subroutine sumintensity00add


end module sumint

module density
contains

subroutine d1sphere(radius_1,rho_1,points,npts)
   real*4, intent(in) :: radius_1,rho_1
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (points(1,i)**2+points(2,i)**2+points(3,i)**2 < radius_1**2) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !OMP END PARALLEL DO
   return
end subroutine d1sphere

subroutine d2cylinder(radius_1,rho_1,points,npts)
   real*4, intent(in) :: radius_1,rho_1
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (points(1,i)**2+points(2,i)**2 < radius_1**2) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !OMP END PARALLEL DO
   return
end subroutine d2cylinder

subroutine d3coreshell(radius_1,radius_2,rho_1,rho_2,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1,rho_2
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   real*4 :: dist
   !$OMP PARALLEL DO
   do i=1,npts
      dist = points(1,i)**2+points(2,i)**2
      if (dist < radius_2**2) then
         points(4,i) = rho_1
      else if (dist < radius_1**2) then
         points(4,i) = rho_2
      else
         points(4,i) = 0
      end if
   end do
   !OMP END PARALLEL DO
   return
end subroutine d3coreshell

subroutine d4gaussian(radius_2,points,npts)
   real*4, intent(in) :: radius_2
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      dist = points(1,i)**2+points(2,i)**2
      points(4,i) = EXP(-dist/(radius_2)**2)
   end do
   !$OMP END PARALLEL DO
end subroutine d4gaussian

subroutine d5choppedcone(radius_1,radius_2,rho_1,z_dim,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1,z_dim
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (SQRT(SUM(points(1:2,i)**2)) < points(3,i)*(radius_2-radius_1)/z_dim+(radius_1+radius_2)/2) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !$OMP END PARALLEL DO
end subroutine d5choppedcone

subroutine d6hexprism(radius_1,rho_1,points,npts)
   real*4, intent(in) :: radius_1,rho_1
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   real*4, dimension(3) :: coords
   real*4, parameter :: sqrt3 = SQRT(3.)
   !$OMP PARALLEL DO
   do i=1,npts
      coords = points(1:3,i)/radius_1
      if ((coords(2)**2>0.75) .OR. (coords(2)+(coords(1)-1)*sqrt3 > 0) .OR. (coords(2)+(coords(1)+1)*sqrt3 < 0) &
          .OR. (coords(2)-(coords(1)-1)*sqrt3 < 0) .OR. (coords(2)-(coords(1)+1)*sqrt3 > 0)) then
         points(4,i) = 0
      else
         points(4,i) = rho_1
      end if
   end do
   !$OMP END PARALLEL DO
end subroutine d6hexprism

subroutine d7rectprism(radius_2,rho_1,points,npts)
   real*4, intent(in) :: radius_2,rho_1
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (points(1,i) < radius_2) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !$OMP END PARALLEL DO
end subroutine d7rectprism

subroutine d11doubleslit(radius_1,radius_2,rho_1,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (((-radius_1/2 < points(1,i)) .AND. (points(1,i) < -radius_2/2)) &
          .OR. ((radius_2/2 < points(1,i)) .AND. (points(1,i) < radius_1/2))) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !$OMP END PARALLEL DO
end subroutine d11doubleslit

subroutine d13sine(radius_1,radius_2,rho_1,rho_2,z_dim,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1,rho_2,z_dim
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   real*4, parameter :: pi = 4.0*ATAN(1.0)
   !$OMP PARALLEL DO
   do i=1,npts
      if (SQRT(SUM(points(1:2,i)**2)) < (radius_1+radius_2)/2 + (radius_1-radius_2)*SIN(points(3,i)*rho_2*2*pi/z_dim)/2) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !$OMP END PARALLEL DO
end subroutine d13sine

subroutine d14doublecone(radius_1,radius_2,rho_1,z_dim,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1,z_dim
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (SQRT(SUM(points(1:2,i)**2)) < radius_2+ABS(points(3,i))*(radius_1-radius_2)/z_dim/2) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !$OMP END PARALLEL DO
end subroutine d14doublecone

subroutine d15elipticalcylinder(radius_1,radius_2,rho_1,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   !$OMP PARALLEL DO
   do i=1,npts
      if (points(1,i)**2/radius_1**2+points(2,i)**2/radius_2**2 < 1) then
         points(4,i) = rho_1
      else
         points(4,i) = 0
      end if
   end do
   !OMP END PARALLEL DO
   return
end subroutine d15elipticalcylinder

subroutine d16asymmhexpyr(radius_1,radius_2,rho_1,z_dim,points,npts)
   real*4, intent(in) :: radius_1,radius_2,rho_1,z_dim
   real*4, dimension(4,npts), intent(inout) :: points
   integer*4, intent(in) :: npts
   real*4 :: x,y,scale_by
   real*4, parameter :: sqrt3 = SQRT(3.)
   !$OMP PARALLEL DO
   do i=1,npts
      scale_by = 2.*z_dim/(z_dim-2.*points(3,i))*MAX(1.,(radius_1+radius_2)/radius_1)/radius_1
      x = points(1,i)*scale_by
      y = points(2,i)*scale_by
      if ((ABS(y) > sqrt3/2.) .OR. (ABS(y) > (-sqrt3*(ABS(x)-radius_2/radius_1)+sqrt3))) then
         points(4,i) = 0
      else
         points(4,i) = rho_1
      end if
   end do
   !OMP END PARALLEL DO
   return
end subroutine d16asymmhexpyr


!Template:
!subroutine d#name(...,points,npts)
!   real*4, intent(in) :: ...
!   real*4, dimension(4,npts), intent(inout) :: points
!   integer*4, intent(in) :: npts
!   !$OMP PARALLEL DO
!   do i=1,npts
!      points(4,i) = ???
!   end do
!   !$OMP END PARALLEL DO
!end subroutine d#name

!Note, unlike Python, fortran:
!array indices start at 1, not 0.
!array slices include the last numbered element.
!multidimensional array order is swapped.

end module density

