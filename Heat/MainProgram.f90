program MainProgram
use TimeMod
use SolverMod



implicit none



double precision , dimension (201,201,4001) :: f
double precision , dimension (201,201) :: dt
double precision :: dx,dy,c,w,TimeStep

integer :: nx,ny,NumSteps,Lboundx,Lboundy,Rboundx,Rboundy,i,j,k

f = 300
NumSteps = 4001
nx = 201
ny = 201
dx = 0.28d0
dy = 0.28d0
c = 22
w = 1.571
TimeStep = 5E-4

Lboundx = (nx+1)/2.0d0 - floor(.5d0/dx)
Rboundx = (nx+1)/2.0d0 + floor(.5d0/dx)
Lboundy = (ny+1)/2.0d0 - floor(.5d0/dy)
Rboundy = (ny+1)/2.0d0 + floor(.5d0/dy)


call RK2(f,nx,ny,dx,dy,c,TimeStep,NumSteps,w)

open (unit  = 12, file = 'output.csv')

do i = 1,ny

write(12,*) (f(i,j,1000),j=1,nx) 

end do

close (unit = 12)

open (unit  = 13, file = 'output1.csv')

do i = 1,ny

write(13,*) (f(i,j,2000),j=1,nx) 

end do

close (unit = 13)

open (unit  = 14, file = 'output2.csv')

do i = 1,ny

write(14,*) (f(i,j,3000),j=1,nx) 

end do

close (unit = 14)


open (unit  = 15, file = 'output3.csv')

do i = 1,ny

write(15,*) (f(i,j,4001),j=1,nx) 

end do

close (unit = 15)
end program


