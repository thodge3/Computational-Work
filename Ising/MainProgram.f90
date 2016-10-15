program MainProgram


double precision, dimension (15,15) :: cfg
double precision :: p, dE, mu, B, J,boltz, E, SummE, SummE2, SummM, SummM2
double precision :: rng1, rng2, A,FF

integer :: x,y,N,x1,x2,y1,y2,i,t,jj,initial,kkk,jjj, hh, gg


call init_random_seed()


N = 15
boltz = 1.3806E-23
mu = 9.27E-24
J = 0.9E-26
!J = 0
B = 0.001
initial = 0

!do i = 1,N
!do jj = 1,N

!	initial  = initial + cfg(i,jj)
!end do
!end do

!E = initial

!print* , E
!print*, '-------------------'



E = 1E-21 * N**2
M = 2* N**2
SummE = 0
SummE2 = 0
SummM = 0
SummM2 = 0

!do hh = 1,N
!	do gg = 1,N

		
!call random_number (FF)
!if (FF<0.5)then
!	cfg(gg,hh) = -1
!end if
!if (FF > 0.5) then
!	cfg(gg,hh) = 1
!end if

!end do
!end do



cfg = -1


open (unit  = 14, file = 'output2.csv')

do kkk = 1,N

write(14,*) (cfg(kkk,jjj),jjj=1,N) 

end do

close (unit = 12)



open (unit  = 15, file = 'output3.csv')
open (unit = 16, file = 'output5.csv')


do t = 1, 1100

if (t == 500)then

open (unit  = 14, file = 'output4.csv')

do kkk = 1,N

write(14,*) (cfg(kkk,jjj),jjj=1,N) 

end do

close (unit = 14)
end if


do i = 1, 90000

call random_number(rng1)
call random_number(rng2)

x = floor(rng1 * N + 1 )
y = floor(rng2 * N + 1 )




cfg(x,y) = cfg(x,y) * (-1)



x1 = x-1
x2 = x+1
y1 = y-1
y2 = y+1

if (x1 == 0) x1 = n
if (x2 == n+1) x2 = 0
if (y1 == 0) y1 = n
if (y2 == n+1) y2 = 0

dE = -2.0d0  * cfg(x,y) * J * (cfg(x1,y) + cfg(x2,y) + cfg(x,y1) + cfg(x,y2) ) - 2.0d0*B*cfg(x,y) * mu




p = exp(-dE / boltz*dble(t))

call random_number(A)

if ( A< p ) then

	E = E + dE
	!M = M + 2 * cfg(x,y)
	M = M - 2 * cfg(x,y)
else 

	cfg(x,y) = cfg(x,y) * (-1)
end if

SummE = SummE + E
SummM = SummM + M
SummE2 = SummE2 + E**2
SummM2 = SummM2 + M**2





end do

write(15,*) (SummM/(N*1E5))
write(16,*) ((-SummE/(N*1E5*boltz*1044)))
!write(15,*) (-SummE/(N*1E5*boltz*1044))
SummE = 0
SummM = 0

if (MOD(t,100) == 0) then
!print*, t
end if

end do

close (unit = 15)
close(unit = 16)

print*, '-------------------'


!print*, (SummE / N)
!print*, (SummM / N)
!print*, (SummE2 / N)
!print*, (SummM2 / N)

open (unit  = 13, file = 'output1.csv')

do kkk = 1,N

write(13,*) (cfg(kkk,jjj),jjj=1,N) 

end do

close (unit = 13)



contains
SUBROUTINE init_random_seed()
     INTEGER :: i, n, clock
     INTEGER, DIMENSION(:), ALLOCATABLE :: seed

     CALL RANDOM_SEED(size = n)
     ALLOCATE(seed(n))

     CALL SYSTEM_CLOCK(COUNT=clock)

     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     CALL RANDOM_SEED(PUT = seed)

     DEALLOCATE(seed)
   END SUBROUTINE

end program
