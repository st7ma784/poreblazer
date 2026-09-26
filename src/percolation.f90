!----------------------------------------------------------
! This module performs percolation analisys of lattice: returns number and sizes of lattice clusters
! and wheather the system is percolated
!----------------------------------------------------------

Module percolation

    Implicit None
    Save

    Private
    Public :: percolation_calc, percolation_calc_simple
    Public :: clusteranalysis, clusteranalysis_exact     ! for tests/test_percolation.f90

    ! Cluster labelling: 0 is Poreblazer 3.0.5's (the default, so results match upstream),
    ! 1 is exact union-find. Poreblazer's labelling can split one connected cluster into
    ! several; see FORK.md.
    Integer, Public :: percolation_labelling = 0


Contains

!-----------------------------------------------------------------------------------------------
!!subroutine that initiates percolation analysis
!-----------------------------------------------------------------------------------------------

    Subroutine percolation_calc(lattice_in, n_sites, nn_sites, cl_summary)

        Integer*1, Dimension(:,:,:), Intent(InOut)                                     :: lattice_in
        Integer, Dimension(:), Intent(InOut), Optional                               :: n_sites
        Integer, Intent(InOut), Optional                                             :: nn_sites
        Integer, Dimension(:,:), Intent(InOut)                                       :: cl_summary
        Integer, Dimension(:,:,:), allocatable                                       :: cluster
        Integer, Dimension(:), allocatable                                           :: cl
        Integer, Dimension(:), allocatable                                           :: trcl
        Integer                                                                      :: nc

        allocate(cluster(size(lattice_in, 1), size(lattice_in, 2), size(lattice_in, 3)))
        allocate(cl(5000000),trcl(5000000))
        nc = 0
        cluster = 0
        cl = 0
        trcl = 0

        If(percolation_labelling == 1) Then
            Call clusteranalysis_exact(lattice_in,cluster,cl,nc)
        Else
            Call clusteranalysis(lattice_in,cluster,cl,trcl,nc)
        End If
        Call span(lattice_in, n_sites, nn_sites, cluster, cl, nc, cl_summary)
        deallocate(cluster)
        deallocate(cl)
        deallocate(trcl)
    End Subroutine percolation_calc

!-----------------------------------------------------------------------------------------------
!!subroutine that initiates percolation analysis
!-----------------------------------------------------------------------------------------------

    Subroutine percolation_calc_simple(lattice_in, cl_summary, spanning)

        Integer*1, Dimension(:,:,:), Intent(InOut)                                       :: lattice_in
        Integer, Dimension(:,:), Intent(InOut)                                         :: cl_summary
        Integer, Dimension(:,:,:), allocatable                                         :: cluster
        Integer, Dimension(:), allocatable                                             :: cl,trcl
        Integer, Intent(InOut)                                                         :: spanning
        Integer                                                                        :: nc

        allocate(cluster(size(lattice_in, 1), size(lattice_in, 2), size(lattice_in, 3)))
        allocate(cl(100000),trcl(100000))
        nc = 0
        cluster = 0
        cl = 0
        trcl = 0
        spanning = 0

        If(percolation_labelling == 1) Then
            Call clusteranalysis_exact(lattice_in,cluster,cl,nc)
        Else
            Call clusteranalysis(lattice_in,cluster,cl,trcl,nc)
        End If

        Call span_simple(lattice_in, cluster, cl, nc, cl_summary, spanning)
        deallocate(cluster, cl, trcl)
    End Subroutine percolation_calc_simple



!---------------------------------------------------------------------
! Subroutine that takes a 3D grid of sites (0 or 1) and performs cluster
! analysis
!---------------------------------------------------------------------

    Subroutine clusteranalysis(ngrid,cluster,cl,trcl,nc)

        Integer*1, Dimension(:,:,:), Intent(InOut)     :: ngrid
        Integer, Dimension(:,:,:), Intent(InOut)  :: cluster
        Integer, Dimension(:), Intent(InOut)      :: cl,trcl
        Integer, Intent(InOut)                       :: nc
        Integer, Dimension(26)                       :: local
        Integer                                      :: i, j, k, l, LX, LY, LZ, scenario, ncold, trcli, icount, i1
        Logical :: loop
        Real :: xr, yr, zr
        Integer, Dimension(:), allocatable                                             :: newlabel


        ! newlabel(t) is the final label of the clusters whose root label is t (0 until the
        ! first site of that cluster is met). Upstream searched a list of the labels so far
        ! for every site, which is quadratic in the number of clusters and overflowed after
        ! 100000 clusters; the labels are the same.
        allocate(newlabel(size(trcl)))
        newlabel = 0

        LX = size(ngrid,1)
        LY = size(ngrid,2)
        LZ = size(ngrid,3)

        do k=1, LZ
            do j=1, LY
                do i=1, LX
                   If(ngrid(i,j,k) == 1) Then
                       If(cluster(i,j,k) /= 0) cycle

                       Call reveal_local(i,j,k,cluster,ngrid,local, scenario)

                       If(scenario == 1) Then
                           Call scenario_1(i,j,k, cluster, ngrid, cl, trcl, nc)
                       Else
                           Call scenario_2(i,j,k, cluster, ngrid, cl, trcl, nc)
                       End If

                    End If

                end do
            end do
        end do
     
        icount = 0 
        nc = 0

        do k=1, LZ
            do j=1, LY
                do i=1, LX

                        If(ngrid(i,j,k) == 0) cycle
                        icount =  icount + 1
                        trcli = trcl(cluster(i, j, k))
                        if(newlabel(trcli) == 0) then
                            nc = nc + 1
                            newlabel(trcli) = nc
                        end if
                        cluster(i, j, k) = newlabel(trcli)
                        cl(newlabel(trcli)) = cl(newlabel(trcli)) + 1
                end do
            end do
        end do
    
     do i=1, nc
     end do

    ! Several additional loops may be required to reconnect the neighbouring clusters

!        do l=1,100
!            ncold=nc
!            do k=1, LZ
!                do j=1, LY
!                    do i=1, LX
!                     If(ngrid(i,j,k) == 0) cycle
!                     write(20,*) i, j, k, cluster(i,j,k), cl(cluster(i,j,k)), nc
!
!                        If(ngrid(i,j,k) == 0) cycle
!                        Call update(i,j,k,cluster,ngrid,cl,nc)
!
!                    end do
!                end do
!            end do
!
!            If(nc==ncold) exit
!         end do
!     print*, "hello2"
     deallocate(newlabel)
     End Subroutine clusteranalysis

!---------------------------------------------------------------------
! Exact cluster labelling: union-find over the occupied sites with the same periodic
! six-neighbour connectivity, clusters numbered in order of their first site in scan
! order (the order clusteranalysis numbers them in). Two sites share a label exactly
! when a path of occupied neighbours joins them.
!---------------------------------------------------------------------

    Subroutine clusteranalysis_exact(ngrid, cluster, cl, nc)
        Integer*1, Dimension(:,:,:), Intent(In)   :: ngrid
        Integer, Dimension(:,:,:), Intent(InOut)  :: cluster
        Integer, Dimension(:), Intent(InOut)      :: cl
        Integer, Intent(InOut)                    :: nc
        Integer, Dimension(:), allocatable        :: parent
        Integer                                   :: i, j, k, LX, LY, LZ, s, r

        LX = size(ngrid,1)
        LY = size(ngrid,2)
        LZ = size(ngrid,3)
        allocate(parent(LX*LY*LZ))
        do s=1, LX*LY*LZ
            parent(s) = s
        end do
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    if(ngrid(i,j,k) /= 1) cycle
                    s = i + LX*(j-1) + LX*LY*(k-1)
                    if(ngrid(modulo(i,LX)+1,j,k) == 1) call unite(s, modulo(i,LX)+1 + LX*(j-1) + LX*LY*(k-1))
                    if(ngrid(i,modulo(j,LY)+1,k) == 1) call unite(s, i + LX*modulo(j,LY) + LX*LY*(k-1))
                    if(ngrid(i,j,modulo(k,LZ)+1) == 1) call unite(s, i + LX*(j-1) + LX*LY*modulo(k,LZ))
                end do
            end do
        end do
        ! Number the clusters in order of their first site; a numbered root holds -number
        cluster = 0
        nc = 0
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    if(ngrid(i,j,k) /= 1) cycle
                    r = find(i + LX*(j-1) + LX*LY*(k-1))
                    if(parent(r) > 0) then
                        nc = nc + 1
                        if(nc > size(cl)) stop "clusteranalysis_exact: more clusters than the cluster size array holds"
                        parent(r) = -nc
                    end if
                    cluster(i,j,k) = -parent(r)
                    cl(-parent(r)) = cl(-parent(r)) + 1
                end do
            end do
        end do
        deallocate(parent)

    contains

        ! The root of x, halving the path on the way (roots point to themselves, or hold
        ! -number once numbered)
        integer function find(x)
            integer, intent(in) :: x
            integer :: y
            y = x
            do while(parent(y) > 0 .and. parent(y) /= y)
                if(parent(parent(y)) > 0) parent(y) = parent(parent(y))
                y = parent(y)
            end do
            find = y
        end function find

        subroutine unite(a, b)
            integer, intent(in) :: a, b
            integer :: ra, rb
            ra = find(a)
            rb = find(b)
            if(ra /= rb) parent(max(ra, rb)) = min(ra, rb)
        end subroutine unite

    End Subroutine clusteranalysis_exact

!---------------------------------------------------------------------
! Subroutine which reveals status of the neighbouring sites
!---------------------------------------------------------------------

    Subroutine reveal_local(i0,j0,k0,cluster,ngrid,local, scenario)
        Integer, Intent(In)                      :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(In)     :: cluster
        Integer*1, Dimension(:,:,:), Intent(In)     ::ngrid
        Integer, Dimension(:), Intent(InOut)        :: local
        Integer, Intent(InOut)                     :: scenario
        Integer                                  :: i1, j1, k1, i, j, k, ic

    ! Default scenario 1: neighbouring sites are not occupied
    ! or occupied but not assigned to any clusters

     scenario = 1
     ic = 1

    ! Six neighbouring sites

     k1=k0-1; j1=j0; i1=i0    
     If(k1<1) k1=size(cluster,3)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then
           local(ic) = cluster(i1,j1,k1)

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
        End If

     k1=k0+1; j1=j0; i1=i0
     If(k1>size(cluster,3)) k1= 1
        If(ngrid(i1,j1,k1)==1) Then                 
           If(cluster(i1,j1,k1)/=0) Then             
           local(ic) = cluster(i1,j1,k1)                 

           ! Scenario 2 where some of the sites are occupied and                 
           ! assigned to some clusters                 

           scenario = 2                 
           Return                 
           End If             
        End If

     k1=k0; j1=j0-1; i1=i0
     If(j1<1) j1=size(cluster,2)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then                 
           local(ic) = cluster(i1,j1,k1)                 

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
        End If
      
     k1=k0; j1=j0+1; i1=i0
     If(j1>size(cluster,2)) j1= 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then                 
           local(ic) = cluster(i1,j1,k1)                 

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
        End If

     k1=k0; j1=j0; i1=i0-1
     If(i1<1) i1=size(cluster,1)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then
           local(ic) = cluster(i1,j1,k1)

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters  

           scenario = 2
           Return
           End If
        End If

     k1=k0; j1=j0; i1=i0+1
     If(i1>size(cluster,1)) i1 = 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then
           local(ic) = cluster(i1,j1,k1)

           ! Scenario 2 where some of the sites are occupied and
           ! assigned to some clusters

           scenario = 2
           Return
           End If
       End If

     Return

    End Subroutine reveal_local

!---------------------------------------------------------------------
! Subroutine which considers scenario 1 where neighbouring sites are not occupied
! or occupied but not assigned to any clusters
!---------------------------------------------------------------------

    Subroutine scenario_1(i0,j0,k0, cluster,ngrid, cl, trcl, nc)
        Integer, Intent(In)                         :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer*1, Dimension(:,:,:), Intent(InOut) :: ngrid
        Integer, Dimension(:), Intent(InOut)      :: cl, trcl
        Integer, Intent(InOut)                      :: nc

        nc = nc + 1
        trcl(nc) = nc
        cluster(i0,j0,k0) = nc
    End Subroutine scenario_1

!---------------------------------------------------------------------
! Subroutine which considers scenario 2 where some of the sites areoccupied and
! assigned to some clusters
!---------------------------------------------------------------------

    Subroutine scenario_2(i0,j0,k0, cluster,ngrid, cl, trcl, nc)
        Integer, Intent(In)                         :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer*1, Dimension(:,:,:), Intent(InOut) :: ngrid
        Integer, Dimension(:), Intent(InOut)      :: cl, trcl
        Integer, Intent(InOut)                      :: nc
        Integer                                     :: i,j,k,i1,j1,k1, lowest, current, trlowest

        lowest = huge(0); trlowest = huge(0)
        ! Six neighbouring sites

        k1=k0-1; j1=j0; i1=i0            
        If(k1<1) k1=size(cluster,3)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        k1=k0+1; j1=j0; i1=i0            
        If(k1>size(cluster,3)) k1= 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)               
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        k1=k0; j1=j0-1; i1=i0
        If(j1<1) j1=size(cluster,2)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)               
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))           

           End If
        End If

        k1=k0; j1=j0+1; i1=i0
        If(j1>size(cluster,2)) j1= 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        k1=k0; j1=j0; i1=i0-1
        If(i1<1) i1=size(cluster,1)
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)               
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))           

           End If
        End If

        k1=k0; j1=j0; i1=i0+1
        If(i1>size(cluster,1)) i1 = 1
        If(ngrid(i1,j1,k1)==1) Then
           If(cluster(i1,j1,k1)/=0) Then

            If(lowest>cluster(i1,j1,k1))  lowest = cluster(i1,j1,k1)
            If(trlowest>trcl(cluster(i1,j1,k1))) trlowest = trcl(cluster(i1,j1,k1))

           End If
        End If

        If(lowest< huge(0)) then
        cluster(i0, j0, k0) = lowest
        trcl(lowest) = trlowest
        else
        nc = nc + 1
        trcl(nc) = nc
        cluster(i0,j0,k0) = nc
        return
        end if

        k1=k0-1; j1=j0; i1=i0
        If(k1<1) k1=size(cluster,3)
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If         

        k1=k0+1; j1=j0; i1=i0 
        If(k1>size(cluster,3)) k1 = 1
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest 
        End If

        k1=k0; j1=j0-1; i1=i0 
        If(j1<1) j1=size(cluster,2)
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If

        k1=k0; j1=j0+1; i1=i0        
        If(j1>size(cluster,2)) j1 = 1
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End IF


        k1=k0; j1=j0; i1=i0-1           
        If(i1<1) i1=size(cluster,1)
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If  

        k1=k0; j1=j0; i1=i0+1           
        If(i1>size(cluster,1)) i1 = 1
        If(cluster(i1,j1,k1)/=0) Then
        trcl(cluster(i1,j1,k1)) = trlowest
        End If

        Return


    End Subroutine scenario_2

!---------------------------------------------------------------------
! Subroutine to update the cluster structure
!---------------------------------------------------------------------

    subroutine permutate(lowest, current, cluster, nc, cl)
        Integer, Intent(In)                         :: lowest, current
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer, Dimension(:), Intent(InOut)      :: cl
        Integer, Intent(InOut)                      :: nc
        Integer                                     :: i,j,k, LX, LY, LZ

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)

        do k=1, LZ
            do j=1, LY
                do i=1, LX

                    If(cluster(i,j,k)==current) Then
                        cl(lowest) = cl(lowest) + 1
                        cluster(i,j,k) = lowest
                    End If

                end do
            end do
        end do

        ! Eliminate current cluster by copying last cluster (nc) into it

        cl(current) = cl(nc)

        do k=1, LZ
            do j=1, LY
                do i=1, LX

                    If(cluster(i,j,k)==nc) Then
                        cluster(i,j,k) = current
                    End If

                end do
            end do
        end do
        cl(nc) = 0
        nc = nc - 1
    End Subroutine permutate

!----------------------------------------------------------------------------
!Subroutine to determine if a percolating cluster exists in x y or z-direction
!----------------------------------------------------------------------------

!There could be more than one spanning cluster in each direction
!Useful output: if there exists a spanning cluster

!Look at cluster information and find out if any one of the clusters has an available site in every x y or z-axis position.

    subroutine span(lattice_in, n_sites, nn_sites, cluster, cl, nc, cl_summary)

        Integer*1, Dimension(:,:,:), Intent(InOut) :: lattice_in
        Integer, Dimension(:), Intent(InOut), Optional :: n_sites     ! the sites of the spanning clusters, if wanted
        Integer, Intent(InOut), Optional          :: nn_sites
        Integer, Dimension(:,:,:), Intent(In)     :: cluster
        Integer, Dimension(:), Intent(In)         :: cl
        Integer, Intent(In)                       :: nc
        Integer, Dimension(:,:), Intent(Out)      :: cl_summary
        Integer                                   :: x_span, y_span, z_span,potentialspan !x_span=0 if there is no spanning cluster
        Integer                                   :: n,i,j,k, LX, LY, LZ, ic, icount, spanning, i1, nn
        Integer, Dimension(:), allocatable        :: cand
        Logical(kind=1), Dimension(:,:), allocatable :: xo, yo, zo
        Integer, Dimension(:), allocatable :: x_array, y_array, z_array
        Logical :: attempt

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)
        allocate(x_array(LX), y_array(LY), z_array(LZ))
        attempt = .False.
        ic = 0
        cl_summary = 0

        ! Which planes each candidate cluster occupies, found in one pass over the grid
        ! (upstream scanned the grid up to three times per candidate); the test per
        ! candidate, in label order, is unchanged
        Call plane_occupancy(cluster, cl, nc, cand, xo, yo, zo)
        do n=1, nc
            if(cand(n) == 0) cycle
            attempt = .True.
            potentialspan = n
            spanning = 0
            if(all(xo(:, cand(n)))) spanning = spanning + 1
            if(all(yo(:, cand(n)))) spanning = spanning + 1
            if(all(zo(:, cand(n)))) spanning = spanning + 1
            if(spanning>0) then
                ic = ic + 1
                cl_summary(1, 1)      = ic
                cl_summary(ic+1, 1)   = potentialspan
                cl_summary(ic+1, 2)   = spanning
            end if
        end do
        deallocate(cand, xo, yo, zo)

        if(attempt.eqv..False.) then
            write(*,*) " The system is NOT percolated in ANY direction "
            spanning = 0

            deallocate(x_array, y_array, z_array)
            return
        end if

        lattice_in = 0
        nn = 0
        if(present(n_sites)) n_sites = 0
        icount = 0
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    icount = icount + 1
                    do i1=1, ic
                        if(cluster(i,j,k) ==  cl_summary(i1+1, 1)) then
                            lattice_in(i,j,k) = 1
                            nn = nn + 1
                            if(present(n_sites)) n_sites(nn) = icount
                        end if
                    end do
                end do
            end do
        end do
        if(present(nn_sites)) nn_sites = nn
        deallocate(x_array, y_array, z_array)
        
    end subroutine span

!----------------------------------------------------------------------------
!Subroutine to determine if a percolating cluster exists in x y or z-direction
!----------------------------------------------------------------------------

!There could be more than one spanning cluster in each direction
!Useful output: if there exists a spanning cluster

!Look at cluster information and find out if any one of the clusters has an available site in every x y or z-axis position.

    subroutine span_simple(lattice_in, cluster, cl, nc, cl_summary, spanning)

        Integer*1, Dimension(:,:,:), Intent(InOut)  :: lattice_in
        Integer, Dimension(:,:,:), Intent(InOut)  :: cluster(:,:,:)
        Integer, Dimension(:), Intent(In)         :: cl
        Integer, Intent(In)                       :: nc
        Integer, Dimension(:,:), Intent(Out)      :: cl_summary
        Integer, Intent(InOut)                    :: spanning
        Integer                                   :: x_span, y_span, z_span,potentialspan !x_span=0 if there is no spanning cluster
        Integer                                   :: n,i,j,k, LX, LY, LZ, ic
        Integer, Dimension(:), allocatable        :: cand
        Logical(kind=1), Dimension(:,:), allocatable :: xo, yo, zo
        Integer, Dimension(:), allocatable  :: x_array, y_array, z_array
        Logical :: attempt

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)
        allocate(x_array(LX), y_array(LY), z_array(LZ))
        attempt = .False.
        ic = 0
        cl_summary = 0

        ! Which planes each candidate cluster occupies, found in one pass over the grid
        ! (upstream scanned the grid up to three times per candidate); the test per
        ! candidate, in label order, is unchanged
        Call plane_occupancy(cluster, cl, nc, cand, xo, yo, zo)
        do n=1, nc
            if(cand(n) == 0) cycle
            attempt = .True.
            potentialspan = n
            spanning = 0
            if(all(xo(:, cand(n)))) spanning = spanning + 1
            if(all(yo(:, cand(n)))) spanning = spanning + 1
            if(all(zo(:, cand(n)))) spanning = spanning + 1
            if(spanning>0) then
                ic = ic + 1
                cl_summary(1, 1)      = ic
                cl_summary(ic+1, 1)   = potentialspan
                cl_summary(ic+1, 2)   = spanning
                deallocate(x_array, y_array, z_array, cand, xo, yo, zo)
                return
            end if
        end do
        deallocate(cand, xo, yo, zo)

        if(attempt.eqv..False.) then
        !    write(*,*) " The system is NOT percolated in ANY direction "
            spanning = 0
        end if

        deallocate(x_array, y_array, z_array)

    end subroutine span_simple

!----------------------------------------------------------------------------
! For each cluster big enough to span (at least LX, LY or LZ sites), cand(n) numbers it
! and xo(i, cand(n)), yo(j, ...) and zo(k, ...) say whether it has a site in plane
! i, j or k. One pass over the grid.
!----------------------------------------------------------------------------

    subroutine plane_occupancy(cluster, cl, nc, cand, xo, yo, zo)
        Integer, Dimension(:,:,:), Intent(In)                  :: cluster
        Integer, Dimension(:), Intent(In)                      :: cl
        Integer, Intent(In)                                    :: nc
        Integer, Dimension(:), allocatable, Intent(Out)        :: cand
        Logical(kind=1), Dimension(:,:), allocatable, Intent(Out) :: xo, yo, zo
        Integer                                                :: n, m, i, j, k, c, LX, LY, LZ

        LX = size(cluster,1)
        LY = size(cluster,2)
        LZ = size(cluster,3)
        allocate(cand(max(1, nc)))
        cand = 0
        m = 0
        do n=1, nc
            if(cl(n)>=LX.or.cl(n)>=LY.or.cl(n)>=LZ) then
                m = m + 1
                cand(n) = m
            end if
        end do
        allocate(xo(LX, max(1, m)), yo(LY, max(1, m)), zo(LZ, max(1, m)))
        xo = .False.
        yo = .False.
        zo = .False.
        if(m == 0) return
        do k=1, LZ
            do j=1, LY
                do i=1, LX
                    c = cluster(i,j,k)
                    if(c < 1 .or. c > nc) cycle
                    c = cand(c)
                    if(c == 0) cycle
                    xo(i, c) = .True.
                    yo(j, c) = .True.
                    zo(k, c) = .True.
                end do
            end do
        end do
    end subroutine plane_occupancy


!---------------------------------------------------------------------
! Subroutine which searches for two clusters sharing a site and
! reconnects them into one large cluster
!---------------------------------------------------------------------

    subroutine update(i0,j0,k0, cluster,ngrid, cl, nc)
        Integer, Intent(In)                         :: i0,j0,k0
        Integer, Dimension(:,:,:), Intent(InOut) :: cluster
        Integer*1, Dimension(:,:,:), Intent(InOut) :: ngrid
        Integer, Dimension(:), Intent(InOut)      :: cl
        Integer, Intent(InOut)                      :: nc
        Integer                                     :: i,j,k,i1,j1,k1, lowest

        lowest = huge(0)
        !First we find the lowest cluster observed among the neighbours
        do k=k0-1, k0+1, 1
            do j=j0-1, j0+1, 1
                do i=i0-1, i0+1, 1

                    if (k/=k0.and.(i/=i0.or.j/=j0)) cycle
                    if (j/=j0.and.(k/=k0.or.i/=i0)) cycle
                    if (i/=i0.and.(j/=j0.or.k/=k0)) cycle

                    i1 = i
                    j1 = j
                    k1 = k

                    If(i1<1) i1=size(cluster,1)
                    If(i1>size(cluster,1)) i1= 1
                    If(j1<1) j1=size(cluster,2)
                    If(j1>size(cluster,2)) j1= 1
                    If(k1<1) k1=size(cluster,3)
                    If(k1>size(cluster,3)) k1= 1

                    If(ngrid(i1,j1,k1)==1) Then
                        If(cluster(i1,j1,k1)/=0) Then
                            If(lowest>cluster(i1,j1,k1)) Then
                                lowest = cluster(i1,j1,k1)
                            End If
                        End If
                    End If

                end do
            end do
        end do

        If(lowest==0) Then
            print*, 'something is wrong, it is not scenario 2'
            stop
        End If

        do k=k0-1, k0+1, 1
            do j=j0-1, j0+1, 1
                do i=i0-1, i0+1, 1

                    if (k/=k0.and.(i/=i0.or.j/=j0)) cycle
                    if (j/=j0.and.(k/=k0.or.i/=i0)) cycle
                    if (i/=i0.and.(j/=j0.or.k/=k0)) cycle

                    i1 = i
                    j1 = j
                    k1 = k

                    If(i1<1) i1=size(cluster,1)
                    If(i1>size(cluster,1)) i1= 1
                    If(j1<1) j1=size(cluster,2)
                    If(j1>size(cluster,2)) j1= 1
                    If(k1<1) k1=size(cluster,3)
                    If(k1>size(cluster,3)) k1= 1

                    If((i1/=i0.or.j1/=j0.or.k1/=k0).and.cluster(i1,j1,k1)==0) cycle

                    If(ngrid(i1,j1,k1)==1) Then
                        If(cluster(i1,j1,k1)==0) Then
                            cl(lowest) = cl(lowest) + 1
                            cluster(i1,j1,k1) = lowest
                        ElseIf(cluster(i1,j1,k1)==lowest)Then
                            cycle
                        Else
                            Call permutate(lowest, cluster(i1,j1,k1), cluster, nc, cl)
                        End If
                    End If

                end do
            end do
        end do
    End Subroutine update

End Module percolation
