! Fortran library for investigating structure and properties of liquid water
! Works for both bulk water and nearby solutes
! SpherePoints and SphereSurfaceAreas are copied from sim.geom.py
! with minor modifications made.


! Adding a quicksort algorithm, should be O(N*logN)!
recursive subroutine quickSort(array)
     implicit none
!     integer, intent(in) :: Dim
     real(kind(1d0)), intent(inout) :: array(:)
!     real(kind(1d0)), dimension(size(array)), intent(out) :: new_array
!     real(kind(1d0)), dimension(:), intent(in) :: array
!     real(kind(1d0)), dimension(size(array)), intent(out) :: new_array
     real(kind(1d0)) :: temp, pivot
     integer :: i, j, last, left, right

!     do i = 1, size(array) 
!        new_array(i) = array(i)
!     enddo
     
     last = size(array)
     if (last.lt.50) then ! use insersion sort on small arrays
        do i=2,last
           temp = array(i)
           do j=i-1,1,-1
              if (array(j).le.temp) exit
              array(j+1) = array(j)
           enddo
           array(j+1) = temp
        enddo
        return
     endif

     ! find median of three pivot
     ! and place sentinels at first and last elements
     temp = array(last/2)
     array(last/2) = array(2)
     if (temp.gt.array(last)) then
        array(2) = array(last)
        array(last) = temp
     else
        array(2) = temp
     endif
     if (array(1).gt.array(last)) then
        temp = array(1)
        array(1) = array(2) 
        array(2) = temp
     endif
     pivot = array(2)

     left = 3 
     right = last - 1
     do 
        do while(array(left).lt.pivot)
           left = left + 1
        enddo
        if (left.ge.right) exit
        temp = array(left)
        array(left) = array(right)
        array(right) = temp
        left = left + 1
        right = right - 1
     enddo
     if (left.eq.right) left = left + 1
     call quicksort(array(1:left-1))
     call quicksort(array(left:))

end subroutine quickSort

! Same as in geom.py
! Finds centroid of given atom set
subroutine Centroid(Pos, Ret, N, Dim)
    implicit none
    integer, intent(in) :: N, Dim
    real(8), dimension(N,Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(out) :: Ret
    Ret = sum(Pos, 1) / real(N)
end subroutine

! Taken from geom.py
subroutine crossProd3(r1, r2, Ret, Dim)
    integer, intent(in) :: Dim
    real(8), dimension(Dim), intent(in) :: r1, r2
    real(8), dimension(Dim), intent(out) :: Ret
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    Ret(1) = r1(2)*r2(3) - r1(3)*r2(2)
    Ret(2) = r1(3)*r2(1) - r1(1)*r2(3)
    Ret(3) = r1(1)*r2(2) - r1(2)*r2(1)
end subroutine

! Same as in geom.py
subroutine reimage(Pos, RefPos, BoxL, ReimagedPos, NPos, Dim)
    implicit none
    integer, intent(in) :: NPos, Dim
    real(8), dimension(NPos, Dim), intent(in) :: Pos
    real(8), dimension(Dim), intent(in) :: RefPos
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NPos, Dim), intent(out) :: ReimagedPos
    integer :: i
    real(8), dimension(Dim) :: distvec, iBoxL
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0)
    do i = 1, NPos
        distvec = Pos(i,:) - RefPos
        distvec = distvec - BoxL * anint(distvec * iBoxL)
        ReimagedPos(i,:) = RefPos + distvec
    enddo            
end subroutine

! Same as in proteinlib.f90
real(8) function RgWeights(Pos, Weights, NAtom, Dim)
  implicit none
  integer, intent(in) :: NAtom, Dim
  real(8), dimension(NAtom,Dim), intent(in) :: Pos
  real(8), dimension(NAtom), intent(in) :: Weights
  real(8), dimension(Dim) :: Center
  integer :: i
  Center = sum(Pos, 1) / real(NAtom)
  RgWeights = 0.
  do i = 1, NAtom
    RgWeights = RgWeights + Weights(i) * sum((Pos(i,:) - Center)**2)
  enddo
  RgWeights = RgWeights / sum(Weights)
  RgWeights = sqrt(RgWeights)
end function

! Same as in geom.py
! Places points on a sphere to later define SASA
subroutine SpherePoints(N, Points)
    implicit none
    integer, intent(in) :: N
    real(8), dimension(N,3), intent(out) :: Points
    real(8) :: off, y, phi, r
    real(8) :: inc 
    real(8), parameter :: pi = 3.1415926535897931D0 
    integer :: k
    inc = pi * (3. - sqrt(5.))
    Points = 0.
    off = 2. / real(N)
    do k = 1, N
        y = real(k-1) * off - 1. + (off * 0.5)
        r = sqrt(max(1. - y*y, 0.))
        phi = real(k-1) * inc
        Points(k,1) = cos(phi)*r
        Points(k,2) = y
        Points(k,3) = sin(phi)*r
    enddo
end subroutine

! Modified to also return Exposed, so can find which points (then atoms) are exposed
subroutine SphereSurfaceAreas(Pos, Radii, Points, nExp, BoxL, Areas, Exposed, NSphere, NPoints, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8), dimension(NPoints, Dim), intent(in) :: Points
    integer, intent(in) :: nExp
    real(8), dimension(Dim), intent(in) :: BoxL
    real(8), dimension(NSphere), intent(out) :: Areas
    logical, dimension(NSphere), intent(out) :: Exposed
    real(8), parameter :: pi = 3.141592653589D0 
    integer, intent(in) :: NSphere, NPoints
    integer :: i, j, k
    real(8), dimension(NPoints,Dim) :: ThisPoints
    real(8) :: AreaPerPoint
    real(8), dimension(NSphere) :: RadiiSq
    real(8), dimension(Dim) :: iPos, jPos
    real(8), dimension(Dim) :: distvec, iBoxL
    logical, dimension(NPoints) :: tempExposed
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    iBoxL = merge(1.d0/BoxL, 0.d0, BoxL >= 0.d0) 
    Areas = 0.
    Exposed = .false.
    RadiiSq = Radii*Radii
    do i = 1, NSphere
        iPos = Pos(i,:)
        AreaPerPoint = 4.*pi*Radii(i)**2 / real(NPoints)
        tempExposed = .true.
        do k = 1, NPoints
            ThisPoints(k,:) = Points(k,:) * Radii(i) + iPos
        enddo
        do j = 1, NSphere
            if (i == j) cycle
            jPos = Pos(j,:)
            distvec = jPos - iPos 
            distvec = distvec - BoxL * anint(distvec * iBoxL)
            jPos = iPos + distvec
            if (.not. any(tempExposed)) exit
            !first check if spheres are far from each other
            if (sum((jPos-iPos)**2) > (Radii(i) + Radii(j))**2) cycle
            do k = 1, NPoints
                if (.not. tempExposed(k)) cycle
                if (sum((ThisPoints(k,:) - jPos)**2) < RadiiSq(j)) tempExposed(k) = .false.
            enddo
        enddo
        Areas(i) = AreaPerPoint * real(count(tempExposed))
        if (count(tempExposed) >= nExp) Exposed(i) = .true.
    enddo
end subroutine

! Same as in geom.py
subroutine SphereVolumes(Pos, Radii, dx, Volumes, NSphere, Dim)
    implicit none
    integer, intent(in) :: Dim
    real(8), dimension(NSphere, Dim), intent(in) :: Pos
    real(8), dimension(NSphere), intent(in) :: Radii
    real(8),  intent(in) :: dx
    real(8), dimension(NSphere), intent(out) :: Volumes
    integer, intent(in) :: NSphere
    real(8), dimension(NSphere) :: RadiiSq
    real(8) :: minDistSq, DistSq, dV
    integer :: i,j
    real(8), dimension(Dim) :: Pos2, minPos, maxPos
    if (Dim /= 3) then
        print *, 'Expecting three dimensions and found', Dim
        stop
    endif
    RadiiSq = Radii*Radii
    Volumes = 0.
    dV = dx*dx*dx
    minPos = (/minval(Pos(:,1) - Radii), minval(Pos(:,2) - Radii), minval(Pos(:,3) - Radii)/)
    maxPos = (/maxval(Pos(:,1) + Radii), maxval(Pos(:,2) + Radii), maxval(Pos(:,3) + Radii)/)
    maxPos = maxPos + dx * 0.5    
    !first do a coarse grid check to see which spheres are where
    Pos2 = minPos
    do while (all(Pos2 < maxPos))
        j = 0
        minDistSq = huge(1.d0)
        do i = 1, NSphere
            DistSq = sum((Pos(i,:) - Pos2)**2)
            if (DistSq < minDistSq .and. DistSq < RadiiSq(i)) then
                minDistSq = DistSq
                j = i
            endif
        enddo
        if (j > 0) Volumes(j) = Volumes(j) + dV
        Pos2(1) = Pos2(1) + dx
        do i = 1, 2
            if (Pos2(i) >= maxPos(i)) then
                Pos2(i) = minPos(i)
                Pos2(i+1) = Pos2(i+1) + dx
            endif
        enddo
    enddo   
end subroutine

! Finds the area of a triangle in 3 dimentions
subroutine triangleArea(Pos, area)
    implicit none
    real(8), dimension(3,3), intent(in) :: Pos
    real(8), intent(out) :: area
    real(8), dimension(3) :: v1, v2
    real(8) :: cosTheta, sinTheta, v1Sq, v2Sq
    v1 = Pos(2,:)-Pos(1,:)
    v2 = Pos(3,:)-Pos(1,:)
    v1Sq = sum(v1**2.0)
    v2Sq = sum(v2**2.0)
    cosTheta = sum(v1*v2)/(v1Sq*v2Sq)**0.5
    sinTheta = (1.0-cosTheta**2.0)**0.5
    area = (v1Sq*v2Sq)**0.5*sinTheta
end subroutine

! Map 3D triangle to 2D
subroutine transformTriangle(Pos, NPos, xyPos)
    implicit none
    integer, intent(in) :: NPos
    real(8), dimension(NPos,3,3), intent(in) :: Pos
    real(8), dimension(NPos,3,2), intent(out) :: xyPos
    real(8), dimension(NPos,3,3) :: newPos
    real(8), dimension(3) :: v1, v2, n, xhat, yhat
    real(8), dimension(3,3) :: R, tempPos, iPos
    integer :: i, k
    external :: crossProd3    

    do i = 1, NPos
       iPos = Pos(i,:,:)
       v1 = iPos(2,:) - iPos(1,:)
       v2 = iPos(3,:) - iPos(1,:)

       call crossProd3(v1,v2,n,3)
       n = n/(sum(n**2.0))**0.5
       xhat = v1/(sum(v1**2.0))**0.5
       call crossProd3(n,xhat,yhat,3)

       R(1,:) = xhat
       R(2,:) = yhat
       R(3,:) = n

       do k = 1, 3 
          tempPos(k,:) = iPos(k,:)-iPos(1,:)
          newPos(i,k,:) = matmul(R, tempPos(k,:))
       enddo
    enddo
    xyPos = newPos(:,:,:2)
end subroutine

! Interpolating value at the centroid position within a triangle
! (Should generalize for any position if this is necessary...)
subroutine propertyBarycentric(Pos, Prop, NPos, newProp)
  implicit None
  integer, intent(in) :: NPos
  real(8), dimension(NPos,3,3), intent(in) :: Pos
  real(8), dimension(NPos,3), intent(in) :: Prop
  real(8), dimension(NPos), intent(out) :: newProp
  real(8), dimension(NPos,3,2) :: newPos
  integer :: i
  external :: transformTriangle

  call transformTriangle(Pos, NPos, newPos)

  do i = 1, NPos
     newProp(i) = sum( Prop(i,:) )/3.0
  enddo 
end subroutine
