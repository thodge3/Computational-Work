program projectile
implicit none

double precision :: coefA, coefB, rho, cD,cL, m, A, vMag, wMag, wxvx, wxvy, wxvz, h, g, ax, ay, az, vMag22, pi, theta
integer :: i,n

double precision, allocatable, dimension(:) :: Vx, Vy, Vz, Rx, Ry, Rz, Wx, Wy, Wz


n = 100000
h = .1
rho = 1.225
cD = .33
cL = .12
!m = .145
m = 12805
A = ((1.775)**2)*(3.14159)

allocate(Vx(n) , Vy(n) , Vz(n) , Rx(n) , Ry(n) , Rz(n) , Wx(n) , Wy(n) , Wz(n) )
g = 9.81

vMag22 = 600
pi = acos(-1.0)
theta = 45 * (pi/180.0)

Vx(1) = vMag22 * cos(theta)
Vy(1) = 0
Vz(1) = vMag22 * sin(theta)
Rx(1) = 0
Ry(1) = 0
Rz(1) = 0
Wx(1) = 0
Wy(1) = 0
Wz(1) = 0

coefA = -1*((rho * A) / (2.0*m))*(cD)
coefB = (rho*cL*A) / (2.0*m)


do i = 1,n-1

vMag = sqrt( (Vx(i)**2) + (Vy(i)**2) + (Vz(i)**2) )
wMag = sqrt( (Wx(1)**2) + (Wy(1)**2) + (Wz(1)**2) )


wxvx = ( (Wy(1) * Vz(i)) - (Wz(1) * Vy(i)) )
wxvy = -( (Wx(1) * Vz(i)) - (Wz(1) * Vx(i)) )
wxvz = ( (Wx(1) * Vy(i)) - (Wy(1) * Vx(i)) )


!ax = ( (coefA * vMag * Vx(i) ) + ((coefB * vMag) ) * (wxvx) )
ax =  ((coefB * vMag) ) * (wxvx) 
!ay = ( (coefA * vMag * Vy(i) ) + ((coefB * vMag) ) * (wxvy) )
ay =  ((coefB * vMag) ) * (wxvy) 
!az = ( (coefA * vMag * Vz(i) ) - g + ((coefB * vMag)) * (wxvz) )
az = - g + ((coefB * vMag)) * (wxvz) 


Vx(i+1) = Vx(i) + ax * h
Vy(i+1) = Vy(i) + ay * h
Vz(i+1) = Vz(i) + az * h


Rx(i+1) = Rx(i) + ( h * ( (Vx(i+1) + Vx(i)) / 2.0 ) )
Ry(i+1) = Ry(i) + ( h * ( (Vy(i+1) + Vy(i)) / 2.0 ) )
Rz(i+1) = Rz(i) + ( h * ( (Vz(i+1) + Vz(i)) / 2.0 ) )

if(Rz(i+1) .lt. 0 )then
exit
end if


if(Rz(i) .gt. 0) then
!print*,Rx(i),Rz(i)
end if
end do


print*, Rx(i), Rz(i)



open (unit = 14, file = 'text1.csv')
open (unit = 15, file = 'text2.csv')
open(unit = 13, file = 'text.csv')
do i = 1,n-1
write(13,*) Rx(i)
write(14,*) Ry(i)
write(15,*) Rz(i)
end do
close(unit = 14)
close (unit = 15)
close(unit = 13)


deallocate(Vx,Vy,Vz,Rx,Ry,Rz,Wx,Wy,Wz)
end program projectile 

