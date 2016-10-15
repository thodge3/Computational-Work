module GaussianElimination

implicit none

contains
subroutine Gaussian(a,b,n)

double precision, dimension (:,:) :: a
double precision, dimension (:) :: b
double precision :: mult
integer :: i,j,k, n


do i = 1,n
	do j = i+1, n
			mult = (-a(j,i)/a(i,i))
		do k = 1,n
			a(j,k) = mult * a(i,k) + a(j,k)
		end do
			b(j) = mult*b(i) + b(j)
	end do
end do

do j = n,1,-1
	do i = j-1,1,-1
		b(i) = (-a(i,j)/a(j,j) * b(j) + b(i))
	end do
end do


do i = 1,n
	b(i) = b(i) / a(i,i)
end do

return
end subroutine 

subroutine checknodes (b,n,nodes,L,V0,pi,m,hbar)
implicit none
integer :: maxnodes,neednodes,n,nodes
double precision :: hbar,pi,m,L,Emax,Emin,E,V0
double precision, dimension(5) :: b

E = -V0
do neednodes = 0,maxnodes 
	Emin = E
	Emax = (hbar**2 * pi**2) / ((2.0*m*L**2) * (neednodes + 1)**2) 
	E = (Emax - Emin) / 2.0

end do
return
end subroutine 
end module
