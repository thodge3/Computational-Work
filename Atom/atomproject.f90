Program atomproject
use SolverMod
use ForceEnergyMod
use StatisticsMod

implicit none


double precision,dimension(2500,513) :: x,y,z,vx,vy,vz,ax,ay,az,KE,vmag2
double precision,dimension(2500) :: TKE, Vp, Val,bound
double precision :: r0,delta,v0,m,h,Vmag, space,starx,stary,starz,fx,fy,fz,f
integer, dimension (1) :: count1
integer :: numSteps,numBodies,numBins,t, i, j,a,b,c,step, maxp,kkk,hhh, timeslice, jj,timeslice2,timeslice3


maxp = 8

bound = 1
r0 = 1E-13
!r0 = 224E-12
!r0 = 2E-10
delta = 2.5E-11
v0 = 4E-36
m =  4.485E-26
h =  8E-9


timeslice = 430
timeslice2 = 702
timeslice3 = 2400

numSteps = 2500
numBodies = maxp**3
numBins = 10





space = .4E-10

starx = maxp * space
stary = maxp * space
starz = maxp * space

open(unit = 13, file = 'newoutput.csv')
do a = 1,maxp
do b = 1,maxp
do c = 1,maxp

step = ((a - 1) * maxp**2) + ((b - 1) * maxp) + (c)

x(1,step) = starx - ((a - 1) * space)
y(1,step) = stary - ((b - 1) * space)
z(1,step) = starz - ((c - 1) * space)

vx(1,step) = 0
vy(1,step) = 0
vz(1,step) = 0

write(13,*) x(1,step) , " , " , y(1,step) , " , " , z(1,step) 
end do 
end do
end do

close(unit = 13)






call Euler(x,y,z,vx,vy,vz,ax,ay,az,h,numBodies,numSteps,m,v0,delta,r0,vmag2)
call Kinetic(vx,vy,vz,m,numSteps,numBodies,KE,TKE)
call Potential(x,y,z,numSteps,numBodies,v0,delta,r0,Vp)

open(unit = 14, file = '2here.csv')
do hhh = 1,numbodies
write(14,*) x(timeslice,hhh) , ",",y(timeslice,hhh),",",z(timeslice,hhh)
end do
close(unit = 14)

open(unit = 15, file = 'ke.csv')
do jj = 1,numbodies
write(15,*) vmag2(timeslice,jj)
end do
close(unit = 15)

open(unit = 18, file = '3here.csv')
do hhh = 1,numbodies
write(18,*) x(timeslice2,hhh) , ",",y(timeslice2,hhh),",",z(timeslice2,hhh)
end do
close(unit = 18)

open(unit = 19, file = '4here.csv')
do hhh = 1,numbodies
write(19,*) x(timeslice3,hhh) , ",",y(timeslice3,hhh),",",z(timeslice3,hhh)
end do
close(unit = 19)

do t = 1,numSteps,10
	call Counter(KE(49,:), bound,count1,numBodies,numBins )

end do



end program atomproject

