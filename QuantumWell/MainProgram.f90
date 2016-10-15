program MainProgram
use BuildMatrix
use GaussianElimination


implicit none

double precision :: hbar,pi,m,l, V0, E, one, dE,Emax,Emin,eps
double precision, dimension (1:1000,1:1000) :: a
double precision, dimension (1:1000) :: b,x,V,Aq
integer :: MaxNodes,n,nodes, i,j,ii,jj,neednodes, ierr
logical :: choice
integer, dimension (1000) :: pivot



n = 1000
hbar = 1.0545E-34
pi = 4*atan(1.0)
m = 9.1E-31
L = 1.0E-9
V0 = ( (hbar**2) * (pi**2) * 5) / (2.0 * m * (L**2) )

MaxNodes = 10
choice = .true.


E = -V0
eps = 1.0E-8

do neednodes = 0,MaxNodes,1
	
	Emin = E
	Emax = ( (hbar**2) * (pi**2) * ( (neednodes + 1)**2)  ) / (2.0d0*m*(L**2))
	E = (Emax-Emin) /2.0d0
	choice = .true.


	do  while(choice)
		call Matrix(n,L,V0,a,b,hbar,m,E)

		call DGESV(n,1,a,n,pivot,b,n,ierr)


		!call Gaussian(a,b,n)

		call checknodes(b,n,nodes,L,V0,pi,m,hbar)

		one = 1.0d0

		nodes = 0

	do i = 1, n-1
		if (sign(one,b(i))*sign(one,b(i+1)) < 0)then
		nodes = nodes + 1
		end if
	end do



	if( (nodes .eq. neednodes) .and.(b(2) - abs(b(n)/b(2))  .lt. eps) ) then

		do i = 1,n
		write(10+ neednodes,*)b(i)
		end do
		Aq = 0
		do i = 1,n

		Aq(neednodes + 1) = Aq(neednodes + 1) + (b(i)**2 * (L/n)) 
		
		
		end do
		
		write(neednodes +100,*)Aq(neednodes+1)
		choice = .false.

	else if(nodes .gt. neednodes) then
		dE = E - Emin
		Emax = E
		E = E - dE/2.0d0
	
	else if( (nodes .lt. neednodes) .or. ( b(2) - abs(b(n) / b(2)) .gt. eps) ) then
	
		dE = E - Emax
		Emin = E
		E = E - dE/2.0d0
	
	end if

	end do
end do



end program
