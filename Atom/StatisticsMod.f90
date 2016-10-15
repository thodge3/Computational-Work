module StatisticsMod
implicit none

contains
subroutine Counter(Val,bound,count1,numBodies,numBins)
implicit none

double precision, dimension (:) :: Val, bound
integer, dimension (:) :: count1
integer :: numBodies, numBins, i,j


do i = 1,numBodies
do j = 1, numBins

if ( Val(i) < bound(j) ) then

count1(j) = count1(j) + 1

exit

end if
end do
end do

end subroutine
end module
