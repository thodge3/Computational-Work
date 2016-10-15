module SolverMod
implicit none
contains
subroutine Euler(x,y,z,vx,vy,vz,ax,ay,az,h,numBodies,numSteps,m,v0,delta,r0,vmag2)
use ForceEnergyMod
implicit none

double precision, dimension(:,:) :: x,y,z,vx,vy,vz,ax,ay,az,vmag2
double precision :: m,v0,delta,r0,h,fx,fy,fz,f
integer :: numBodies,numSteps
integer :: i,j,k

ax = 0.0
ay = 0.0
az = 0.0

do i = 1,numSteps -1
do j = 1,numBodies
do k = j+1,numBodies

call Force(x(i,j),y(i,j),z(i,j),x(i,k),y(i,k),z(i,k),fx,fy,fz,f,r0,delta,v0)
ax(i,j) = ax(i,j) + fx/m 
ay(i,j) = ay(i,j) + fy/m
az(i,j) = az(i,j) + fz/m

ax(i,k) = ax(i,k) - (fx/m)
ay(i,k) = ay(i,k) - (fy/m)
az(i,k) = az(i,k) - (fz/m)

end do

vx(i+1,j) = vx(i,j) + h*ax(i,j)
vy(i+1,j) = vy(i,j) + h*ay(i,j)
vz(i+1,j) = vz(i,j) + h*az(i,j)

vmag2(i+1,j) = sqrt(vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)

x(i+1,j) = x(i,j) + h*((vx(i,j) + vx(i+1,j))/2)
y(i+1,j) = y(i,j) + h*((vy(i,j) + vy(i+1,j))/2)
z(i+1,j) = z(i,j) + h*((vz(i,j) + vz(i+1,j))/2) 

end do
if (  mod(i,100) == 0 ) then
print*,i,"Euler"
end if
end do


return
end subroutine
end module
