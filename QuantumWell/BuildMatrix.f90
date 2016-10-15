module BuildMatrix


contains
subroutine Matrix(n,L,V0,a,b,hbar,m,E)

implicit none

double precision, dimension (1:n,1:n) :: a
double precision, dimension (1:n) :: b,x,V
double precision :: L,dx,dxsquared,V0,hbar,m, E
integer :: i,j,k,n




L = 1.0E-9
!a(1:n,1:n)
a = 0.0d0
dx = L / (dble(n) + 1.0E-9)

dxsquared = dx**2

do i = 1,n

	x(i) = dble(i) * dx
	V(i) = V0 * ((x(i) - 0.5E-9)**2 - L)
	a(i,i) = dxsquared * (V(i) - E) + (hbar**2 / m)

end do


do j = 1,n-1

	a(j,j+1) = -hbar**2 / (2.0d0 * m)
	a(j+1,j) = -hbar**2 / (2.0d0 * m)

end do


b = 0.0d0
b(2) = 1.0d0



a(1,1) = 1.0d0
a(2,2) = 1.0d0
a(1,2) = 0
a(2,1) = 0
a(2,3) = 0



return

end subroutine

end module


