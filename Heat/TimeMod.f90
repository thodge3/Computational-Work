module TimeMod

contains
subroutine TimeDerivative(f,dx,dy,c,nx,ny,dt)

double precision , dimension (:,:) :: f
double precision :: dx,dy,c
integer :: nx,ny
double precision , dimension (:,:) :: dt
integer :: i,j

do i = 1,nx-1
do j = 2,ny-1

	dt(i,j) = (c*(f(i+1,j) - 2*f(i,j) + f(i-1,j)) / dx**2) + (c*(f(i,j+1) - 2*f(i,j) + f(i,j-1)) / dy**2)

end do
end do
return
end subroutine 
end module
