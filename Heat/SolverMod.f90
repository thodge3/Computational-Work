module SolverMod
implicit none
contains


subroutine RK2(f,nx,ny,dx,dy,c,TimeStep,NumSteps,w)
use TimeMod

implicit none
double precision,dimension(:,:,:) :: f
double precision,dimension(201,201) :: dt
double precision :: dx,dy,c,w,TimeStep

integer :: nx,ny,NumSteps,Lboundx,Rboundx,Lboundy,Rboundy

integer :: i,j,k




Lboundx = (nx+1)/2.0d0 - floor(.5d0/dx)
Rboundx = (nx+1)/2.0d0 + floor(.5d0/dx)
Lboundy = (ny+1)/2.0d0 - floor(.5d0/dy)
Rboundy = (ny+1)/2.0d0 + floor(.5d0/dy)


do i = 1,NumSteps -1

do j = Lboundx,Rboundx
do k = Lboundy,Rboundy

f(k,j,1) = 300 + 150*sin(w*TimeStep)

end do 
end do
call TimeDerivative(f(:,:,i),dx,dy,c,nx,ny,dt)

do j = 1,nx
do k = 1,ny

f(k,j,i+1) = f(k,j,i) + dt(k,j) * TimeStep/2.0d0
end do
end do


do j = Lboundx,Rboundx
do k = Lboundy,Rboundy

f(k,j,i+1) = 300 + 150*sin(w*(i-0.5d0)*TimeStep)

end do
end do

call TimeDerivative(f(:,:,i+1),dx,dy,c,nx,ny,dt)

do j = 1,nx
do k = 1,ny

f(k,j,i+1) = f(k,j,i) + dt(k,j) * TimeStep

end do 
end do



do j = Lboundx,Rboundx
do k = Lboundy,Rboundy

f(k,j,i+1) = 300 + 150*sin(w*i* TimeStep)

end do 
end do

end do
return
end subroutine

subroutine RK4(f,nx,ny,dx,dy,c,TimeStep,NumSteps,w)
use TimeMod

implicit none
double precision,dimension(:,:,:) :: f
double precision,dimension(201,201) :: dt1,dt2,dt3,dt4
double precision :: dx,dy,c,w,TimeStep

integer :: nx,ny,NumSteps,Lboundx,Rboundx,Lboundy,Rboundy

integer :: i,j,k


Lboundx = (nx+1)/2.0d0 - floor(.5d0/dx)
Rboundx = (nx+1)/2.0d0 + floor(.5d0/dx)
Lboundy = (ny+1)/2.0d0 - floor(.5d0/dy)
Rboundy = (ny+1)/2.0d0 + floor(.5d0/dy)



do j = Lboundx,Rboundx
do k = Lboundy,Rboundy

f(k,j,1) = 300 + 150*sin(w*TimeStep)

end do 
end do

do i = 1,NumSteps-1
	

	
	call TimeDerivative(f(:,:,i),dx,dy,c,nx,ny,dt1)


	do j = 1,nx
		do k = 1,ny
			f(k,j,i+1) = f(k,j,i) + dt1(k,j) * TimeStep/2.0d0
		end do
	end do

	do j = Lboundx,Rboundx
		do k = Lboundy,Rboundy
			f(k,j,i+1) = 300 + 150*sin(w*(i-0.5d0) * TimeStep)
		end do
	end do

	call TimeDerivative(f(:,:,i+1),dx,dy,c,nx,ny,dt2)
	
	do j = 1,nx
		do k = 1,ny
			f(k,j,i+1) = f(k,j,i) + dt2(k,j) * TimeStep/2.0d0
		end do
	end do

	do j = Lboundx,Rboundx
		do k = Lboundy,Rboundy

			f(k,j,i+1) = 300 + 150*sin(w*(i-0.5d0) * TimeStep)
		end do
	end do

	call TimeDerivative(f(:,:,i+1),dx,dy,c,nx,ny,dt3)

	do j = 1,nx
		do k = 1,ny
			f(k,j,i+1) = f(k,j,i) + dt3(k,j)*TimeStep
		end do
	end do

	do j = Lboundx,Rboundx
		do k = Lboundy,Rboundy
			f(k,j,i+1) = 300 + 150*sin(w*i*TimeStep)
		end do
	end do


	call TimeDerivative(f(:,:,i),dx,dy,c,nx,ny,dt4)
	
	do j = 1,nx
		do k = 1,ny
			f(k,j,i+1) = f(k,j,i) + (1/6.0d0)*(dt1(k,j) + (2.0d0*dt2(k,j)) + (2.0d0*dt3(k,j)) + dt4(k,j))
		end do
	end do

	do j = Lboundx,Rboundx
		do k = Lboundy,Rboundy
			
			f(k,j,i+1) = 300 + 150*sin(w*i*TimeStep)
		end do
	end do
end do

return
end subroutine























end module

