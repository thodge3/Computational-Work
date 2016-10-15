program projectile
implicit none

real*16 ::  cD,cL, A, vMag, wMag, wxvx, wxvy, wxvz, h, g, ax, ay, az, T0, P0, gasc, L, molar, qqq, CxVx,CxVy,CxVz&
,vM2, chem, sI, vMag1, vMag2, dryMass, spefImpulse, resul, resul3, theta,pi
integer :: i,n

real*16, allocatable, dimension(:) :: Vx, Vy, Vz, Rx, Ry, Rz, Wx, Wy, Wz, coefA, coefB, rho, pressure, abtemp, Cx, Cy,Cz&
, coefC,m, resul2

logical, allocatable, dimension (:) :: mask



n = 100000
h = 0.1
!rho = 1.2
cD = .33
cL = .12
!m=12500
A = ((1.775)**2)*(3.14159)
gasc = 8.314
L = 0.0065
molar = 0.0289
!molar = 0.055
allocate(Vx(n) , Vy(n) , Vz(n) , Rx(n) , Ry(n) , Rz(n) , Wx(n) , Wy(n) , Wz(n), rho(n), coefA(n), coefB(n), pressure(n),&
	abtemp(n), Cx(n), Cy(n), Cz(n), coefC(n), m(n), mask(n), resul2(n) )
g = 9.81
vM2 = 1715
pi = acos(-1.0d0)
theta = 45 * (pi/180.0)

!spefImpulse = 239
!spefImpulse = 125
dryMass = 4008


!m(1) = 2970000
m(1) = 12805
Vx(1) = vM2 * cos(theta)
Vy(1) = 0
Vz(1) = vM2 * sin(theta)
Rx(1) = 0
Ry(1) = 0
Rz(1) = 0
!!!!!!!!!!!!!!!!! Rotation of Rocket !!!!!!!!!!!!!!!!!!!!!
Wx(1) = 0
Wy(1) = 0
Wz(1) = 0
!!!!!!!!!!!!!!!!! Rotation of Earth (Corliolis Effect) !!!!!!!!!!!!!!!!!!!!!
Cx(1) = 0
Cy(1) = 0
!Cz(1) = 0.001
Cz(1) = 7.29E-5
!Cz(1) = 0.01
!!!!!!!!!!!!!!!!! Initial Density of Air at Ground level !!!!!!!!!!!!!!!!!!!!!
rho(1) = 1.225

coefA(1) = -1*((rho(1) * A) / (2.0*m(1)))*(cD)
coefB(1) = (rho(1)*cL*A) / (2.0*m(1))
coefC(1) = 1.0/m(1)

abtemp(1) = 288.15
pressure(1) = 101.325
T0 = 288.15
P0 = 101.325
!sI = 239
!chem = spefImpulse * g

!!!!!!!!!!!!!!!!! Exponent of Changing Pressure !!!!!!!!!!!!!!!!!!!!!
qqq = (g*molar) / (gasc * L)
print*, qqq



!do i = 2,n-1
do i = 1, n-1
if (Rz(i) .gt. 30000 ) then
exit
end if


if (Rz(i) .lt. -100) then
exit
end if



!write(12, *),Vx(i) , Vy(i) , Vz(i) , Rx(i) , Ry(i) , Rz(i) , Wx(i) , Wy(i) , Wz(i), rho(i), coefA(i), coefB(i), pressure(i),&
!	abtemp(i), Cx(i), Cy(i), Cz(i), coefC(i) 


vMag = sqrt( (Vx(i)**2) + (Vy(i)**2) + (Vz(i)**2) )
wMag = sqrt( (Wx(1)**2) + (Wy(1)**2) + (Wz(1)**2) )

vMag1 = sqrt( (Vx(i)**2) + (Vy(i)**2) + (Vz(i)**2) )

wxvx = ( (Wy(1) * Vz(i)) - (Wz(1) * Vy(i)) )
wxvy = -( (Wx(1) * Vz(i)) - (Wz(1) * Vx(i)) )
wxvz = ( (Wx(1) * Vy(i)) - (Wy(1) * Vx(i)) )


!!!!!!!!!!!!!!!!! Coriolis Effect !!!!!!!!!!!!!!!!!!!!!

CxVx = -2*  ( (Cy(1) * Vz(i)) - (Cz(1) * Vy(i)))
CxVy = -2*  (-((Cz(1) * Vx(i)) - (Cx(1) * Vz(i))))
CxVz = -2*  ((Cx(1) * Vy(i)) - (Cy(1) * Vx(i)))



ax = ( (coefA(i) * vMag * Vx(i) ) + ((coefB(i) * vMag) ) * (wxvx)  + (coefC(1)*CxVx*vMag) )
ay = ( (coefA(i) * vMag * Vy(i) ) + ((coefB(i) * vMag) ) * (wxvy)  + (coefC(1)*CxVy*vMag) )
az = ( (coefA(i) * vMag * Vz(i) ) - g + ((coefB(i) * vMag)) * (wxvz) + ( coefC(1)*CxVz*vMag) )

!Rx(i+1) = 2*Rx(i) - Rx(i-1) + ax(i) * h**2 
!Ry(i+1) = 2*Ry(i) - Rx(i-1) + ay(i) * h**2
!Rz(i+1) = 2*Rz(i) - Rz(i-1) + az(i) * h**2

!Vx(i) = Vx(i-1) + (((ax(i-1) + ax(i)) / 2.0d0) * h)
!Vy(i) = Vy(i-1) + (((ay(i-1) + ay(i)) / 2.0d0) * h)
!Vz(i) = Vz(i-1) + (((az(i-1) + az(i)) / 2.0d0) * h)


Vx(i+1) = Vx(i) + ax * h
Vy(i+1) = Vy(i) + ay * h
Vz(i+1) = Vz(i) + az * h

vMag2 = sqrt( (Vx(i+1)**2) + (Vy(i+1)**2) + (Vz(i+1)**2) )

Rx(i+1) = Rx(i) + ( h * ( (Vx(i+1) + Vx(i)) / 2.0 ) )
Ry(i+1) = Ry(i) + ( h * ( (Vy(i+1) + Vy(i)) / 2.0 ) )
Rz(i+1) = Rz(i) + ( h * ( (Vz(i+1) + Vz(i)) / 2.0 ) )


!Vx(n) = Vx(n-1) + (((ax(n-1) + ax(n)) / 2.0d0) * h)
!Vy(n) = Vy(n-1) + (((ay(n-1) + ay(n)) / 2.0d0) * h)
!Vz(n) = Vz(n-1) + (((az(n-1) + az(n)) / 2.0d0) * h)

!!!!!!!!!!!!!!!!! Equations for the Varying Density as a function of Height !!!!!!!!!!!!!!!!!!!!!


abtemp(i+1) = T0 - (L*Rz(i+1) )

pressure(i+1) = p0* (1.0d0- (L*Rz(i+1)) / T0)**(qqq)

rho(i+1) = ((pressure(i+1)*1000) * molar)  / (gasc * abtemp(i+1))

!m(i+1) =  m(i) *  exp(((vMag2-vMag1) )/(chem))
m(i+1) = m(i) - (130.0*h)
!m(i+1) = m(i) - ((spefImpulse*g*(dble(i)*h)))
!m(i+1) = m(1)
if (m(i+1) .lt. dryMass)then
m(i+1) = dryMass
end if


coefA(i+1) = -1*((rho(i+1) * A) / (2.0*m(i+1)))*(cD)
coefB(i+1) = (rho(i+1)*cL*A) / (2.0*m(i+1))

if(abtemp(i+1) .lt. 1E-10) then
abtemp(i+1) = 1E-10
end if
print*, Rx(i), Ry(i), Rz(i)


end do

resul = maxval(Rz)
resul2 = maxloc(Rz)
resul3 = maxval(resul2)

print*, resul, Rx(277),resul3

print*, '------------------------'

!print*, wxvx, wxvy, wxvz
!print*, coefC(1)
!print*, CxVx, CxVy, CxVz
print*,'------------'
!print*, Vx(1), Vz(1), Rx(940)



open (unit = 14, file = 'text1.csv')
open (unit = 15, file = 'text2.csv')
open(unit = 13, file = 'text.csv')
do i = 1,19201
!do i = 1,n-1
write(13,*) Rx(i)
write(14,*) Ry(i)
write(15,*) Rz(i)
end do
close(unit = 14)
close (unit = 15)
close(unit = 13)


deallocate(Vx,Vy,Vz,Rx,Ry,Rz,Wx,Wy,Wz,rho,coefA,coefB,pressure,abtemp)
end program projectile 

