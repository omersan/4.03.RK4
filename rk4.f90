!-----------------------------------------------------------------------------!
! MAE 5093 DOCUMENTATION | Engineering Numerical Analysis
!-----------------------------------------------------------------------------!
! >>> Basic ODE solvers
!     Runke-Kutta 4
!-----------------------------------------------------------------------------!
! References: 
! * Fundamentals of Engineering Numerical Analysis by P. Moin (2012) -- Ch.1
! * Numerical Recipes: The Art of Scientific Computing, 2nd Edition (1992) 
!-----------------------------------------------------------------------------!
! Written by Omer San
!            CFDLab, Oklahoma State University, cfdlab.osu@gmail.com
!            www.cfdlab.org
! 
! Last updated: Sep. 8, 2015
!-----------------------------------------------------------------------------!


program rk4
implicit none
real*8 ::h,tmax,t,w
integer::j,k,n,np,ne
real*8,allocatable ::y(:),r(:),k1(:),k2(:),k3(:),k4(:)

!Solve y'' + w*y = 0
w = 4.0d0

h = 0.10d0
tmax = 6.0d0

n = nint(tmax/h)

allocate(y(2))
allocate(r(2))
allocate(k1(2))
allocate(k2(2))
allocate(k3(2))
allocate(k4(2))

!Initial condition
y(1) = 1.0d0
y(2) = 0.0d0

ne = 2 !number of equation

t = 0.0d0

open(12, file="numerical.plt")
write(12,*)'variables ="t","f","df"'
write(12,*) t,y(1),y(2)


!Time integration with RK4 scheme
do j=1,n

t = dfloat(j)*h

    call RHS(ne,w,y,t,r)
	do k=1,ne
		k1(k) = h*r(k)      
    	end do
	
    call RHS(ne,w,y+k1/2.0d0,t+h/2.0d0, r)
	do k=1,ne
		k2(k) = h*r(k)      
    	end do

    call RHS(ne,w,y+k2/2.0d0,t+h/2.0d0,r)
	do k=1,ne
		k3(k) = h*r(k)      
    	end do

    call RHS(ne,w,y+k3,t+h,r)
	do k=1,ne
		k4(k) = h*r(k)      
    	end do
    
 	do k=1,ne
		y(k) = y(k) + (k1(k)+2.0d0*(k2(k)+k3(k))+k4(k))/6.0d0     
    end do

    write(12,*) t,y(1),y(2)      
end do

close(12)
   
! Writing exact solution using np points
np = 2000 !use 2000 points to plot curve
open(12, file="exact_curve.plt")
write(12,*)'variables ="t","f","df"'
	do j=0,np
		t = dfloat(j)*(tmax)/(dfloat(np))
		write(12,*) t,dcos(w*t),-w*dsin(w*t)
	end do
close(12)

end

!-----------------------------------------------------------------------------!
!Right Hand Side (RHS)
!-----------------------------------------------------------------------------!
subroutine RHS(ne,w,y,t,r)
implicit none
integer::ne
real*8 ::y(ne),r(ne),w,t

r(1) =  y(2)
r(2) = -w*w*y(1)

end 








