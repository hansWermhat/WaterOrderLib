! Fortran library for basic sorting algorithms

subroutine genRandInt(m, n, randInt)
  implicit none
  integer, intent(in) :: m, n
  integer, intent(out) :: randInt
  integer :: s(13), i, j
  real :: r

  ! simple (but annoying) method of assuring randomness                            
  do j = 1, 13
     call system_clock(count=s(j))
  enddo
  i = 1
  call random_seed(size=i)
  call random_seed(put=s(1:13))

  ! call random number r                                                           
  call random_number(r)

  ! grab the corresponding random integer
  randInt = (m) + floor(((n-1)+1-(m))*r)
end subroutine genRandInt

! Adding a depth first sorting algorithm, should be O(|n_verts| + |n_edges|)!
recursive subroutine depthFirstSort(vertex, array, visited, M, N, final_visited)
     implicit none
     integer, intent(in) :: N
     !integer, intent(in) :: vertex
     integer, dimension(N,N), intent(in) :: array
     integer*4 :: vertex
     integer*4 :: visited(:)
     integer, dimension(N), intent(out) :: final_visited
     integer*4 :: M
!     integer, dimension(N), intent(inout) :: visited
     integer, dimension(N) :: arrayRow
     integer, dimension(M) :: trueInds
     integer :: i, count, dummy, iVert
     
     count = 1

     arrayRow = array(vertex, :)

     ! do loop to find all 1-indices in "arrayRow"
     do i = 1, N
        if (arrayRow(i).eq.1) then
           trueInds(count) = i
           count = count + 1
        else
           cycle
        endif 
     enddo

     visited(vertex) = 1

     ! set up a dummy variable 
     dummy = M
     ! do loop through trueInds to perform recursion
     do i = 1, dummy
        iVert = trueInds(i)
        if ( visited(iVert).eq.1 ) then
           cycle
        else
           M = sum(array(iVert, :))
           final_visited = visited
           call depthFirstSort(iVert,array,visited,M,N, final_visited)
           final_visited = visited
        endif

     enddo
        
end subroutine depthFirstSort

! Adding a depth first sorting algorithm, should be O(|n_verts| + |n_edges|)!
recursive subroutine quickSort(array, first, last)
  implicit none
  integer first, last
  real array(*)
  integer pivotInd, i, j, N
  real :: pivot
  real :: temp

  N = last - first + 1
  pivotInd = first + floor(N/2.0)
  pivot = array(pivotInd)

  i = first; j = last
  do 
     do while (array(i)<pivot)
        i = i + 1
     enddo

     do while (pivot < array(j))
        j = j - 1
     enddo

     if ( i >= j ) then 
        exit
     endif

     temp = array(i)
     array(i) = array(j)
     array(j) = temp
     i = i+1; j=j-1
  enddo

  if (first==last) then
     print *, first, last
     print *, i, j
     print *, array(first), array(last)
     print *, array(i), array(j)
     stop 123
  endif
!    print *, first, i-1
!    print *, last, j+1
!    stop 123
  if (first < i-1) then 
     !       print *, first, i-1
     call quickSort(array, first, i-1)
  endif
  if (j+1 < last) call quickSort(array, j+1, last) 
end subroutine quickSort

subroutine qSort(array, array_output, N)
  implicit none
  integer, intent(in) :: N
  real, dimension(N), intent(in) ::  array
  real, dimension(N), intent(out) ::  array_output
  
  array_output = array
  call quickSort(array_output, 1, N)

end subroutine qSort
!program test
!  use quicksort_mod
!  integer :: N
!  real, dimension(100000000) :: array

!  N = 100000000

!  call random_number(array)
  
!  call quickSort(array, 1, N)
!  print *, array(1:4)

!end program
