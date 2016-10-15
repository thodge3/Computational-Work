program NewMain



implicit none

double precision, dimension (10,10,10) :: H
double precision, dimension (10000,10000,10000) :: psi
double precision, dimension (10) :: x,W
double precision, dimension (1000) :: work
double precision :: L,dx,dxsquared,V0,hbar,mass,a,pi,b
integer :: i,j,k,n,m, ierr, lwork,ii,jj,kk,numpoints, num, nee, nums,hh,g
character :: JOBZ = 'V', UPLO = 'L'
double precision :: start, finish

num = 10
lwork = 3*num-1
hbar = 1.0545E-34
pi = acos(-1.0)
mass = 9.1E-31
numpoints = 100;
L = 1.0E-9
a = 0.0d0
b = L
nums = 10000
dx = L / (dble(numpoints) )

V0 = ( (hbar**2) * (pi**2) * 5) / (2.0 * mass * (L**2) )
dxsquared = dx**2


do n = 1, num
	do m = 1, num
	do g = 1, num
		
		H(m,g,n) = (trap(f2,a,b,nums,L,n,m) - trap(f,a,b,nums,L,n,m))*2*(V0/L)
	end do
	end do
H(n,n,n) = H(n,n,n) + ((( (hbar**2) * (n**2) * (pi**2)  ) / (2*mass * (L**2))) - ((3.0d0 / 4.0d0) * V0))


end do


call cpu_time(start)

call DSYEV(JOBZ,UPLO,num,H,num,W,work,lwork,ierr)

call cpu_time(finish)

!print '("Time = ",f6.5," Seconds.")', finish-start
!print*, ierr



!do i = 1,num
!	do j = 1,num

!	print*, H(j,i)
!	end do
!end do




psi = 0
do ii = 1, numpoints
do hh = 1, num
	do jj = 1, num
		do kk = 1, num

	psi(hh,jj,ii) = psi(hh,jj,ii) + H(hh,kk,jj) * sqrt(2.0d0/(L**2)) &
	*sin((kk*pi*(dble(ii)-1.0d0)*dx)/L)*sin((kk*pi*(dble(ii) - 1.0d0)*dx)/L)
!print*, psi(jj,ii)

		end do 
	end do
	end do
end do




!do i = 1,3
	do k = 1, numpoints
!	do j = 1, numpoints

	!write(9 + i , *) (psi(1,k,j))
	write(13,*) (psi(i,k,1),j=1,numpoints) 

	end do
!	end do
!end do







contains

double precision function f2(x,n,m,L) 

double precision :: x,pi,L
integer :: n,m

pi = 4*atan(1.0d0)

	f2 = x**2 * sin(n*pi*x / L) * sin(m*pi*x/L)

return
end function



double precision function f(x,n,m,L) 

double precision :: x,pi,L
integer :: n,m

pi = 4*atan(1.0d0)

	f = x * sin(n*pi*x / L) * sin(m*pi*x/L)

return
end function f



double precision function trap(f,a,b,numsplit,L,n,m)
implicit none
double precision :: h, Tot, xj, f,Funct, pi,L,a,b
integer :: j,numsplit,k,n,m

!numsplit = 2168
pi = 4.0d0*atan(1.0d0)
a = 0.0d0
!L = 2.0E-9
b = L
h = (b-a)/numsplit

tot  = 0.0d0
xj = 0.0d0

do j = 1,n-1
	xj = a + (j  * h)
	Tot = Tot + f(xj,n,m,L) 
end do

Tot = ((2.0d0*Tot) + f(a,n,m,L) + f(b,n,m,L)) * (h/2.0d0)

!print*,n,  Tot
end function trap 




!double precision function Trapezoid(Funct,n,dx,m,L)

!implicit none
!double precision :: Funct
!double precision :: dx,L
!integer :: n,m,i
!double precision :: summ

!nee = floor(L / dx)

!summ = 0

!do i = 2 , nee-1

!	summ = dx * Funct((i-1)*dx,n,m,L) + summ
!end do

!summ = summ + dx*(Funct(0,n,m,L) + Funct((n-1)*dx,n,m,L) / 2.0d0)

!return

!end function




end program
