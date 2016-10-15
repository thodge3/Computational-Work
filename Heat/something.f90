
module Heat

implicit none

subroutine RK2(f,nx,ny,
use timeDerivative,
do i = 1.numsteps

call timeDerivative(f,nx,ny,dx,dy,dt,b)

do j = 1,nx
do k = 1,ny

f(k,j,i+1) = f(k,j,i) + dt(k,j) * h/2.0

end do
end do
call timeDerivative(



end program





