module ForceEnergyMod
implicit none
contains
subroutine Force(x1,y1,z1,x2,y2,z2,fx,fy,fz,f,r0,delta,v0)
implicit none
double precision  :: x1,y1,z1,x2,y2,z2,r0,v0,delta
double precision  :: fx,fy,fz
double precision :: f,r

r = sqrt((x2-x1)**2 + ( y2-y1)**2 + (z2-z1)**2)

f = ((2*v0)/delta)* (1-exp(-(r-r0)/delta)) * (exp(-(r-r0)/delta))


fx = f * ((x2-x1)/r)
fy = f * ((y2-y1)/r)
fz = f * ((z2-z1)/r)


return
end subroutine

subroutine Kinetic(vx,vy,vz,m,numSteps,numBodies,KE,TKE)
implicit none

double precision , dimension (:,:), intent (in) :: vx,vy,vz
double precision , dimension (:,:), intent (out):: KE
double precision , dimension (:), intent (out):: TKE

integer :: numSteps, numBodies, i , j
double precision :: m, Vmag

TKE = 0

open(unit = 17 , file = 'kinetic.csv')
do i = 1,numSteps
do j = 1,numBodies

Vmag = (vx(i,j)**2 + vy(i,j)**2 + vz(i,j)**2)

KE(i,j) = (1.0/2.0) * m * Vmag

TKE(i) = TKE(i) + KE(i,j)
end do
write(17,*) TKE(i)
if ( mod(i,100) == 0) then
print*,i,"kinetic"
end if
end do
close(unit = 17)


end subroutine 

subroutine Potential(x,y,z,numSteps,numBodies,v0,delta,r0,Vp)

implicit none

double precision, dimension (:), intent(out) :: Vp
double precision, dimension (:,:),intent(in) :: x,y,z
integer :: numSteps, numBodies, i, j, k
double precision :: v0, delta, r0, dx,dy,dz,r

open(unit = 16,file = 'potential.csv')

do i = 1,numSteps
do j = 1,numBodies
do k = j+1,numBodies

dx = x(i,j) - x(i,k)
dy = y(i,j) - y(i,k)
dz = z(i,j) - z(i,k)


r = sqrt(dx**2 + dy**2 + dz**2)

Vp(i) = v0*(1-exp(-(r-r0)/delta))**2 - v0 + Vp(i) 

end do 
end do
write(16,*) Vp(i) 
if ( mod(i,100) == 0) then
print*,i,"potential"
end if
end do
close(unit = 16)
end subroutine 










end module ForceEnergyMod



