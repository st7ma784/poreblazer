! Unit tests of the percolation analysis's cluster labelling.
!
!   gfortran -O2 src/percolation.f90 tests/test_percolation.f90 -o test_percolation && ./test_percolation
!
! 1. The smallest lattice we found on which Poreblazer 3.0.5's labelling (clusteranalysis)
!    splits a connected cluster: 8 sites, one cluster, labelled as two. The default
!    labelling must keep doing this (the fork matches upstream); the exact labelling
!    must find one cluster.
! 2. Exact labelling (clusteranalysis_exact) against an independent flood fill on random
!    periodic lattices: the same clusters, numbered in order of their first site.
program test_percolation
    use percolation, only: clusteranalysis, clusteranalysis_exact
    implicit none
    integer :: failures, trial, L, nc, nref
    integer*1, allocatable :: grid(:,:,:)
    integer, allocatable :: cluster(:,:,:), ref(:,:,:), cl(:), trcl(:)
    real, allocatable :: r(:,:,:)
    real :: p
    integer :: seed(64)

    failures = 0
    allocate(cl(5000000), trcl(5000000))

    ! 1. The 8-site counterexample (a 6 x 6 x 1 lattice; nothing touches the edges)
    allocate(grid(6,6,1), cluster(6,6,1))
    grid = 0
    grid(2,3,1) = 1; grid(5,3,1) = 1
    grid(2,4,1) = 1; grid(4,4,1) = 1; grid(5,4,1) = 1
    grid(2,5,1) = 1; grid(3,5,1) = 1; grid(4,5,1) = 1
    cluster = 0; cl = 0; trcl = 0; nc = 0
    call clusteranalysis(grid, cluster, cl, trcl, nc)
    call check(nc == 2, "Poreblazer's labelling splits the 8-site cluster into 2 (as upstream)")
    cluster = 0; cl = 0; nc = 0
    call clusteranalysis_exact(grid, cluster, cl, nc)
    call check(nc == 1, "exact labelling finds 1 cluster in the 8-site lattice")
    call check(cl(1) == 8, "exact labelling puts all 8 sites in it")
    deallocate(grid, cluster)

    ! 2. Random periodic lattices
    seed = 12345
    call random_seed(put=seed(1:size_of_seed()))
    do trial = 1, 200
        L = 4 + mod(trial * 7, 27)
        p = 0.2 + 0.5 * mod(trial * 13, 100) / 100.0
        allocate(grid(L, L + 1, L + 2), cluster(L, L + 1, L + 2), ref(L, L + 1, L + 2), r(L, L + 1, L + 2))
        call random_number(r)
        grid = 0
        where (r < p) grid = 1
        cluster = 0; cl = 0; nc = 0
        call clusteranalysis_exact(grid, cluster, cl, nc)
        call flood_fill(grid, ref, nref)
        if (nc /= nref .or. any(cluster /= ref)) then
            call check(.false., "exact labelling matches a flood fill on a random lattice")
            print '(a, i4, a, i4, f6.2, a, 2i8)', '    trial', trial, ' L', L, p, ' clusters', nc, nref
        end if
        if (sum(cl(1:nc)) /= count(grid == 1)) call check(.false., "exact cluster sizes add up to the sites")
        deallocate(grid, cluster, ref, r)
    end do
    call check(failures == 0, "exact labelling matches a flood fill on 200 random lattices")

    if (failures > 0) then
        print '(i0, a)', failures, ' failure(s)'
        stop 1
    end if
    print '(a)', 'All percolation tests passed'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (ok) then
            print '(a, a)', 'ok    ', what
        else
            print '(a, a)', 'FAIL  ', what
            failures = failures + 1
        end if
    end subroutine check

    integer function size_of_seed()
        call random_seed(size=size_of_seed)
        size_of_seed = min(size_of_seed, 64)
    end function size_of_seed

    ! Clusters of occupied sites with periodic six-neighbour connectivity, by breadth-first
    ! flood fill from each unlabelled site in scan order (so numbered by first site)
    subroutine flood_fill(grid, label, n)
        integer*1, intent(in) :: grid(:,:,:)
        integer, intent(out) :: label(:,:,:)
        integer, intent(out) :: n
        integer, allocatable :: queue(:,:)
        integer :: LX, LY, LZ, i, j, k, head, tail, d, x, y, z, nb(3)
        integer, parameter :: step(3, 6) = reshape([1,0,0, -1,0,0, 0,1,0, 0,-1,0, 0,0,1, 0,0,-1], [3, 6])
        LX = size(grid, 1); LY = size(grid, 2); LZ = size(grid, 3)
        allocate(queue(3, LX * LY * LZ))
        label = 0
        n = 0
        do k = 1, LZ; do j = 1, LY; do i = 1, LX
            if (grid(i,j,k) /= 1 .or. label(i,j,k) /= 0) cycle
            n = n + 1
            label(i,j,k) = n
            head = 1; tail = 1
            queue(:, 1) = [i, j, k]
            do while (head <= tail)
                x = queue(1, head); y = queue(2, head); z = queue(3, head)
                head = head + 1
                do d = 1, 6
                    nb(1) = modulo(x - 1 + step(1, d), LX) + 1
                    nb(2) = modulo(y - 1 + step(2, d), LY) + 1
                    nb(3) = modulo(z - 1 + step(3, d), LZ) + 1
                    if (grid(nb(1), nb(2), nb(3)) == 1 .and. label(nb(1), nb(2), nb(3)) == 0) then
                        label(nb(1), nb(2), nb(3)) = n
                        tail = tail + 1
                        queue(:, tail) = nb
                    end if
                end do
            end do
        end do; end do; end do
    end subroutine flood_fill

end program test_percolation
