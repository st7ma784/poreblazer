!--------------------------------------------------------------------------------------
! This is the main program for characterization of porous solids
!
! Poreblazer 3.0.5 / Windows XP
!
! Lev Sarkisov 2018
!--------------------------------------------------------------------------------------

module parameters
        use fundcell, only: Fundamental_Cell
        type(Fundamental_Cell) :: fcell                  ! contains information about the strucutre of the simulation cell
        real*8  :: hicut, hicut2                          ! cutoff potential distance (A) and its squared value
        real*8  :: coeff_surface, coeff_surface2         ! additional coefficient to switch between hard sphere (1) or potential minima (1.122) accessible surface
        real*8  :: temp                                  ! temperature (K) required for Helium porosimetry calculations
        integer :: nsample                               ! number of tests for surface area
        integer :: iseed                                 ! random seed number
        integer :: vis_option                            ! visualization options
end module parameters

module atoms
        integer*4, allocatable :: atype(:)               ! this array is natoms long, for each atom i containing the number corresponding to the type of atom
        character*10, allocatable :: aname(:)            ! this array is ntypes long, for each atom i containing the name corresponding to the type of atom
        real*8, allocatable :: asigma(:), asigma2(:)     ! these arrays are ntypes long, for each atom type i containing the collision diameter for that type and squared value of it
        real*8, allocatable :: aeps(:)                   ! this array is ntypes long, for each atom type i containing the LJ epsilon, K
        real*8, allocatable :: amass(:)                  ! this array is ntypes long, for each atom type i containing molar mass
        real*8, allocatable :: asigma2_he(:), aeps_he(:) ! these arrays contain mixed values for sigma and epsilon for atom type i interacting with atom of helium
        real*8, allocatable :: asigma_n(:), asigma2_n(:) ! these arrays contain mixed values for sigma and squared value of it for atom type i interacting with atom of nitrogen
end module atoms

module adsorbent
        use  vector, only:       vectype
        integer*4 :: natypes, natoms                                       ! number of atom types, and atoms in adsorbent structure
        real*8, allocatable :: coords(:, :), coords_temp(:,:)              ! coordinates of the adsorbent (natoms long array)
        type(vectype), allocatable :: matvec(:)                            ! 
        character*10, allocatable :: adsname(:)                            ! names of all atoms in the input coordinates file
end module adsorbent

module lattice
        real*8                   :: cube_size                                         ! size of the lattice cubelet
        integer                  :: ncubesx, ncubesy, ncubesz                         ! number of lattice cubes in each direction and other parameters
        integer*8                :: ntot                                              ! total number of lattice cubes
        integer                  :: spanning                                          ! logical variable for percolation
        integer*8                :: nhe_cubes, nn_cubes, ng_cubes                     ! number of lattice cubes accessible to helium, nitrogen and a point, respectively (whole grid)
        ! The grid is split over MPI ranks in slabs of z planes (module parallel): this
        ! rank holds planes z0..z1, and every grid array below is dimensioned (:, :, z0:z1).
        ! With one rank (and in the serial build) z0 = 1 and z1 = ncubesz.
        integer                  :: z0 = 1, z1 = 0
        integer                  :: ng_local = 0, nn_local = 0                        ! the counts in this rank's slab
        integer*1, allocatable   :: lattice_space(:,:,:)                              ! lattice cubes: cube (i,j,k) is 1 if occupied, 0 if empty (inaccessible)
        integer*1, allocatable   :: lattice_temp(:,:,:)                               ! temp array for various uses: cube (i,j,k) is 1 if occupied, 0 if empty (inaccessible)
        integer, allocatable     :: cl_summary(:,:)                                   ! cluster summary: stores information about percolating clusters and the number of percolated dimensions
        real*8, allocatable      :: lattice_rdist2(:,:,:)                             ! squared distances between the center of the cubelet and nearest atom of the structure
        integer*1, allocatable   :: lattice_space_he(:,:,:), lattice_space_n(:,:,:)   ! lattice cubes: cube (i,j,k) is 1 if occupied, 0 if empty (inaccessible) for helim and nitrogen
        integer, allocatable     :: n_cubes(:)                                        ! nitrogen-accessible cubes of this slab, by local cube number (sized to their count)
        real*8, allocatable      :: lattice_lj_he(:,:,:)                              ! for each cube, summary of the LJ interaction between the He atom in the center of the cube and the structure

contains

    ! The indices (j, k, l) of this slab's local cube number icount (numbered with j
    ! fastest, then k, then l from plane z0); l is the plane in the whole grid
    subroutine cube_ijk(icount, j, k, l)
        integer, intent(in)  :: icount
        integer, intent(out) :: j, k, l
        j = mod(icount - 1, ncubesx) + 1
        k = mod((icount - 1)/ncubesx, ncubesy) + 1
        l = (icount - 1)/(ncubesx*ncubesy) + z0
    end subroutine cube_ijk

    ! The local cube number of cube (j, k, l) of this slab
    integer function cube_local(j, k, l)
        integer, intent(in) :: j, k, l
        cube_local = ((l-z0) * ncubesx * ncubesy) + ((k-1) * ncubesx) + j
    end function cube_local

end module lattice


module distributions
        real*4, allocatable :: psd_cumul(:), psd(:)                     ! cumulative pore volume Vp function and negative derivative of that function, i.e. pore size distribution
        real*4              :: binsize, maxpore                                       ! size of the bin for PSD analysis (A), maximum anticipated pore diameter (A)
        integer             :: nbins                                                  ! number of bins for PSD analysis
end module distributions

module results
        real*8 :: volume, smass, density, rho, pore_v_he, pore_geom, stotal, pore_lim, pore_max
        integer :: sys_perc
end module results

module atomcells
        ! Atoms binned into cells, for searches that need only nearby atoms. The cells are
        ! boxes in the cell's slanted coordinates, in which Poreblazer's minimum image wraps
        ! each component, so the gap g between two cells bounds the slanted separation of
        ! any points in them. The Cartesian distance is |U d| for the unslant matrix U, and
        ! |U d|^2 >= metric |d|^2, where metric is a lower bound on the smallest eigenvalue
        ! of U^T U (atomcells_metric): so a search to reach/sqrt(metric) in slanted
        ! coordinates finds every atom within reach. metric is 1 in orthorhombic cells.
        implicit none
        integer                  :: ncell(3) = 1, ncells = 1                          ! cells along each side, in total
        real*8                   :: cellw(3) = 0.0d0                                  ! cell widths (A)
        integer, allocatable     :: cell_start(:), cell_atoms(:)                      ! atoms of cell n: cell_atoms(cell_start(n):cell_start(n+1)-1), ascending

contains

    ! Bin the atoms (matvec) into cells about target A wide, or into one cell
    subroutine atomcells_build(target, use_cells)
        use parameters, only: fcell
        use adsorbent, only:  natoms, matvec
        real*8, intent(in)       :: target
        logical, intent(in)      :: use_cells
        integer                  :: i, d, n, c(3)
        integer, allocatable     :: atom_cell(:), fill(:)

        if(allocated(cell_start)) deallocate(cell_start, cell_atoms)
        ncell = 1
        cellw = fcell%ell
        if(use_cells) then
            do d=1, 3
                ncell(d) = max(1, int(fcell%ell(d)/target))
                cellw(d) = fcell%ell(d)/dble(ncell(d))
            end do
        end if
        ncells = ncell(1)*ncell(2)*ncell(3)
        allocate(atom_cell(natoms), fill(ncells), cell_start(ncells+1), cell_atoms(natoms))
        cell_start = 0
        do i=1, natoms
            do d=1, 3
                c(d) = atomcells_index(matvec(i)%comp(d), d)
            end do
            atom_cell(i) = 1 + c(1) + ncell(1)*(c(2) + ncell(2)*c(3))
            cell_start(atom_cell(i)+1) = cell_start(atom_cell(i)+1) + 1
        end do
        cell_start(1) = 1
        do n=1, ncells
            cell_start(n+1) = cell_start(n+1) + cell_start(n)
        end do
        fill = cell_start(1:ncells)
        do i=1, natoms
            cell_atoms(fill(atom_cell(i))) = i
            fill(atom_cell(i)) = fill(atom_cell(i)) + 1
        end do
        deallocate(atom_cell, fill)
    end subroutine atomcells_build

    ! A lower bound on the smallest eigenvalue of U^T U, U the unslant matrix (Gershgorin's
    ! circles, less a margin for rounding): 1 in orthorhombic cells, 0.5 for a hexagonal
    ! cell. Zero or less when the cell is too oblique for the bound to help.
    real*8 function atomcells_metric()
        use parameters, only: fcell
        use fundcell, only:   fundcell_isortho
        real*8                   :: u(3,3), g(3,3)
        integer                  :: i, j
        if(fundcell_isortho(fcell)) then
            atomcells_metric = 1.0d0
            return
        end if
        u = fcell%unslantmatrix%comp
        g = matmul(transpose(u), u)
        atomcells_metric = huge(0.0d0)
        do i=1, 3
            atomcells_metric = min(atomcells_metric, g(i,i) - sum(abs(g(i,:))) + abs(g(i,i)))
        end do
        atomcells_metric = atomcells_metric*(1.0d0 - 1.0d-9)
    end function atomcells_metric

    ! The cell (from 0) along side d of coordinate x, wrapped into the cell
    integer function atomcells_index(x, d)
        use parameters, only: fcell
        real*8, intent(in)       :: x
        integer, intent(in)      :: d
        real*8                   :: s
        s = x - fcell%ell(d)*floor(x/fcell%ell(d))
        atomcells_index = max(0, min(int(s/cellw(d)), ncell(d)-1))
    end function atomcells_index

    ! The cells along side d that can lie within reach of a point in cell a, each once
    ! (periodic), with the smallest gap between a point in cell a and one in each
    subroutine atomcells_neighbours(a, d, reach, nb, gap, n)
        integer, intent(in)      :: a, d
        real*8, intent(in)       :: reach
        integer, intent(out)     :: nb(:), n
        real*8, intent(out)      :: gap(:)
        integer                  :: m, b, o
        n = 0
        m = int(reach/cellw(d)) + 1
        if(2*m+1 >= ncell(d)) then
            do b=0, ncell(d)-1
                o = abs(b - a)
                o = min(o, ncell(d) - o)
                n = n + 1
                nb(n) = b
                gap(n) = max(0, o-1)*cellw(d)
            end do
        else
            do o=-m, m
                n = n + 1
                nb(n) = modulo(a + o, ncell(d))
                gap(n) = max(0, abs(o)-1)*cellw(d)
            end do
        end if
    end subroutine atomcells_neighbours

end module atomcells

program poreblazer
use parallel, only: par_init, par_finalize, par_rank
implicit none
interface initialize
    subroutine initialize(arg)
        character(len=*), optional :: arg
    end subroutine initialize
end interface initialize


character(len=256) :: arg
character(len=16) :: timing_env
logical :: timing
integer*8 :: clock0, clock1, rate

! With MPI, rank 0 writes the log and every output file; the other ranks only compute
call par_init()
if(par_rank > 0) open(6, file='/dev/null', status='old', action='write')
call get_environment_variable("POREBLAZER_TIMING", timing_env)
timing = len_trim(timing_env) > 0

write(*,*) "!-------------------------------------------------------!"
write(*,*) "! Welcome to Poreblazer_v3.0.5/Linux/Windows XP         !"
write(*,*) "! Developed by: Lev Sarkisov, 2012/18                   !"
write(*,*) "! Cite: Sarkisov, Harrison, Molecular Simulation, 2011  !"
write(*,*) "!-------------------------------------------------------!"
write(*,*)

!--------------------------------------------------------------------------------------
!       The program first goes through several preliminary initializations, associated
!       with reading parameters in.
!--------------------------------------------------------------------------------------

write(*,*) "Step 1/8"
if (command_argument_count() > 0) then
    call get_command_argument(1, arg)
    call initialize(arg)
else
    call initialize
end if

!--------------------------------------------------------------------------------------
!       Now we will perform a number of preliminary calculations:
!       - distance from each lattice cube center to all probes
!       - lattice accessible to Helium and cluster accessible to Helium, Kh(He), pore volume
!       - lattice accessible to Nitrogen
!       - Limiting probe size
!--------------------------------------------------------------------------------------

write(*,*) "Step 2/8"
call system_clock(clock0, rate)
call lattice_calculations
call report_time("lattice")

!--------------------------------------------------------------------------------------
!      Now we sequentially call several routines where the actual
!      analysis of the structure is performed
!--------------------------------------------------------------------------------------

write(*,*) "Step 3/8"
call system_clock(clock0, rate)
call helium_lattice              ! returns lattice accessible to helium
call report_time("helium lattice")
write(*,*) "Step 4/8"
call system_clock(clock0, rate)
call nitrogen_lattice            ! returns lattice accessible to nitrogen
call report_time("nitrogen lattice")
write(*,*) "Step 5/8"
call system_clock(clock0, rate)
call helium_volume               ! helium volume
call report_time("helium volume")
write(*,*) "Step 6/8"
call system_clock(clock0, rate)
call surface_area                ! accessible surface area
call report_time("surface area")
write(*,*) "Step 7/8"
call system_clock(clock0, rate)
call pore_distribution           ! geometric pore size distribution
call report_time("pore size distribution")
write(*,*) "Step 8/8"
call system_clock(clock0, rate)
call limiting_diameter           ! calculation of the pore limiting diamater
call report_time("limiting diameter")

call finalize

write(*,*) "!-------------------------------------------------------!"
write(*,*) "! Simulation complete                                   !"
write(*,*) "!-------------------------------------------------------!"

call par_finalize()

contains

    ! With POREBLAZER_TIMING set, rank 0 writes each step's wall time to stderr (not the log)
    subroutine report_time(step)
        character(len=*), intent(in) :: step
        call system_clock(clock1)
        if(timing .and. par_rank == 0) write(0, '(a, a, f10.3, a)') " timing: ", step, dble(clock1 - clock0)/dble(rate), " s"
    end subroutine report_time

end program

!--------------------------------------------------------------------------------------
! THE END OF THE PROGRAM
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
! SUBROUTINES
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
! Initialization: the main purpose of this subroutine is to read in various types of
! simulation parameters and initialize arrays
!--------------------------------------------------------------------------------------

subroutine initialize(filename)

    use parameters
    use atoms
    use adsorbent
    use lattice
    use distributions

    use defaults, only: rdbl
    use fundcell, only: fundcell_init, fundamental_cell,  fundcell_slant, fundcell_unslant
    use  vector, only:  vectype
    use random, only:   random_init
    use results
    use percolation, only: percolation_labelling, percolation_set_slab
    use parallel, only: par_rank, par_size, par_slab, par_abort, par_bcast_line

    implicit none
    character(len=1024)                   :: stdin_lines(3)
    character(len=*), optional :: filename                                            ! input filename, defaults stdin
    integer*4                             :: i, j, k, l                                 ! cycle indicies
    logical                               :: check                                      ! logical variable, required for input analysis
    character(256)                        :: filename1, filename2, filename3            ! data file names
    real*8                                :: sigma_he, eps_he, sigma_n                  ! sigma (A) and epsilon (K) of helium; sigma (A) of nitrogen atom
    type(vectype)                         ::  atvec1
    character(256)                        :: visline                                    ! the visualisation line of defaults.dat
    integer                               :: ios
    character(256)                        :: offset_env                                 ! POREBLAZER_GRID_OFFSET
    real*8                                :: grid_offset(3)                             ! its fractions of a cubelet

    ! Initialization cycle

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Proceed with initialization and data input            !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    if (present(filename)) then
        open(4, file=trim(adjustl(filename)), status='old')
        read(4,'(a)') filename1
        read(4,'(a)') filename2
        read(4,'(a)') filename3

        read(4,*) fcell%ell(1), fcell%ell(2), fcell%ell(3)    ! non-orthorhombic dimensions of the simulation cell
        read(4,*) fcell%anglebc, fcell%angleac, fcell%angleab ! cell angles  (RDbl)  (alpha,beta,gamma)
        close(4)
    else
        filename3 = "defaults.dat"
        ! Only rank 0 can read stdin under mpirun: it reads the three lines and shares them
        stdin_lines = ' '
        if(par_rank == 0) then
            do i=1, 3
                read(*,'(a)') stdin_lines(i)
            end do
        end if
        do i=1, 3
            call par_bcast_line(stdin_lines(i))
        end do
        read(stdin_lines(1),*) filename2  !  name of the xyz coordinates file
        read(stdin_lines(2),*) fcell%ell(1), fcell%ell(2), fcell%ell(3)    ! non-orthorhombic dimensions of the simulation cell
        read(stdin_lines(3),*) fcell%anglebc, fcell%angleac, fcell%angleab ! cell angles  (RDbl)  (alpha,beta,gamma)
    end if

    open(2, file=Trim(adjustl(filename2)), status='old')
    open(3, file=Trim(adjustl(filename3)), status='old')

    read(3,*) filename1  !  name of the atom types file
    read(3,*) sigma_he, eps_he, sigma_n, temp
    read(3,*) hicut, nsample

    open(1, file=Trim(adjustl(filename1)), status='old')

    hicut2 = hicut*hicut

    coeff_surface = 1.122 ! hard-coded now
    coeff_surface2 = coeff_surface*coeff_surface

    ! Here we read in the number of atom types and their collision diameters

    read(1,*) natypes !  number of atom types
    allocate (aname(natypes),asigma(natypes),aeps(natypes), &
    amass(natypes), asigma2(natypes), asigma2_he(natypes), aeps_he(natypes), asigma_n(natypes), asigma2_n(natypes))

    do i=1, natypes  ! reading in characteristics of atoms
        read(1,*) aname(i), asigma(i), aeps(i), amass(i)
        asigma2(i)    = asigma(i)*asigma(i)
        asigma2_he(i) = (0.5*(asigma(i) + sigma_he))**2.0
        aeps_he(i)    = sqrt(aeps(i)*eps_he)
        asigma_n(i)   = (0.5*(asigma(i) + sigma_n))
        asigma2_n(i)  = (0.5*(asigma(i) + sigma_n))**2.0
    end do

    ! Now we initialize the adsrobent structure

    read(2,*) natoms !  number of atoms in the strucutre
    allocate (adsname(natoms), coords(3, natoms), coords_temp(3, natoms), atype(natoms), matvec(natoms))

    ! Here we read in x, y, z, coordinates of the atoms from the input file
    ! The format of the file is xyz  file

    read(2,*)

    smass = 0.0 ! mass of the unit cell g/mol

    do i=1, natoms
        read(2,*) adsname(i), coords(1, i), coords(2, i), coords(3, i)

        check=.False.

        do j=1, natypes
            if(aname(j)==trim(adjustl(adsname(i)))) then
                atype(i)=j
                smass = smass + amass(j)
                check=.True.
                exit
            end if
        end do
        if(.Not.check) then
            print*, 'No such atom type:   ', trim(adjustl(adsname(i))), aname(i)
            stop
        end if
    end do

    ! Parameters of the simulations



    read(3,*) cube_size            ! lattice cube_size size in A and number of rotations per molecule

    ! Parameters of the simulation box

    fcell%origin%comp=0.0                                 ! coordinates of the origin of the cell (0,0,0)

    call fundcell_init(fcell)

    do i=1, natoms
    matvec(i)%comp = coords(:,i)
    end do

    rho = (smass/fcell%volume)*1.660539

    fcell%eff(1) = fcell%ell(1)
    fcell%eff(2) = fcell%ell(2)
    fcell%eff(3) = fcell%ell(3)
  
      do i=1, natoms
      matvec(i) = fundcell_slant(fcell, matvec(i))
      end do

    do i=1, natoms
    coords_temp(:,i) =  matvec(i)%comp 
    end do

    fcell%mx =  minval(coords_temp(1, :))
    fcell%my =  minval(coords_temp(2, :))
    fcell%mz =  minval(coords_temp(3, :))

    coords_temp(1, :) = coords_temp(1, :) - fcell%mx
    coords_temp(2, :) = coords_temp(2, :) - fcell%my
    coords_temp(3, :) = coords_temp(3, :) - fcell%mz

    do i=1, natoms
    matvec(i)%comp = coords_temp(:,i)
    end do

    atvec1%comp(1) = fcell%mx
    atvec1%comp(2) = fcell%my
    atvec1%comp(3) = fcell%mz

    atvec1 = fundcell_unslant(fcell, atvec1)

    ! Adjusting the position of the system in case it is seriously outside the simulation cell

    coords(1, :) = coords(1, :) - atvec1%comp(1)
    coords(2, :) = coords(2, :) - atvec1%comp(2)
    coords(3, :) = coords(3, :) - atvec1%comp(3)


    ncubesx = int(fcell%eff(1)/cube_size)                ! number of cubelets on x side
    cube_size = fcell%eff(1)/dble(ncubesx)               ! corrected size of the cubelet
    ncubesx = int(fcell%eff(1)/cube_size)                ! number of cubelets on x side
    ncubesy = int(fcell%eff(2)/cube_size)                ! number of cubelets on y side
    ncubesz = int(fcell%eff(3)/cube_size)                ! number of cubelets on z side

    ! With POREBLAZER_GRID_OFFSET="fx fy fz" (fork only), the structure moves against the grid by
    ! these fractions of a cubelet along the cell's sides, to measure how much the results depend
    ! on where the grid falls; without it the grid starts at the lowest atom, as in 3.0.5
    call get_environment_variable("POREBLAZER_GRID_OFFSET", offset_env)
    if(len_trim(offset_env) > 0) then
        read(offset_env,*,iostat=ios) grid_offset
        if(ios /= 0) call par_abort("POREBLAZER_GRID_OFFSET needs three numbers, fractions of a cubelet")
        write(*,'(a, 3f9.5)') "  Grid offset (cubelets):", grid_offset
        atvec1%comp = grid_offset*cube_size
        do i=1, natoms
            matvec(i)%comp = matvec(i)%comp + atvec1%comp
        end do
        atvec1 = fundcell_unslant(fcell, atvec1)
        do i=1, natoms
            coords(:, i) = coords(:, i) + atvec1%comp
        end do
    end if

    if(par_rank == 0) then
    open(40, file = 'new.xyz')
    write(40,*) natoms     !  number of atoms in the strucutre

    write(40,*)

    do i=1, natoms
        write(40,'(a, 3f10.5)') adsname(i), coords(1, i), coords(2, i), coords(3, i)
    end do

    close(40)
    end if

    ntot = int(ncubesx, 8)*ncubesy*ncubesz               ! total number of lattice cubelets in the system

    ! This rank's slab of z planes (all of them on one rank)
    if(ncubesz < par_size) call par_abort("More MPI ranks than grid planes along z; use fewer ranks")
    call par_slab(ncubesz, par_rank, z0, z1)
    if(int(ncubesx, 8)*ncubesy*(z1 - z0 + 1) > huge(0)) &
        call par_abort("More than 2**31 grid cubes on one rank; use more MPI ranks or a larger cubelet size")
    call percolation_set_slab(ncubesz, z0)

    allocate(lattice_space(ncubesx, ncubesy, z0:z1))   ! lattice_temp: in limiting_diameter
    allocate(lattice_rdist2(ncubesx, ncubesy, z0:z1))
    allocate(lattice_space_he(ncubesx, ncubesy, z0:z1))
    allocate(lattice_space_n(ncubesx, ncubesy, z0:z1))
    allocate(lattice_lj_he(ncubesx, ncubesy, z0:z1))
    allocate(cl_summary(20,2))

    lattice_space = 0
    lattice_rdist2 = 0.0
    lattice_space_he = 0
    lattice_space_n = 0
    nhe_cubes = 0;  nn_cubes = 0; ng_cubes = 0
    lattice_lj_he = 0.0d0
    cl_summary = 0

    read(3,*) maxpore, binsize    ! information related to pore size distrubution
    nbins = int(maxpore/binsize)
    binsize = maxpore/real(nbins)

    allocate(psd_cumul(nbins+100), psd(nbins+100))

    read(3,*) iseed               ! random seed

    if(iseed>0) then
        iseed = -iseed
    end if

    call random_init(iseed)


    ! The visualisation option, optionally followed by the cluster labelling (fork only):
    ! 0 (or absent) is Poreblazer 3.0.5's labelling, 1 is exact union-find labelling
    read(3,'(a)') visline
    read(visline,*) vis_option
    read(visline,*,iostat=ios) vis_option, percolation_labelling
    if(ios /= 0) percolation_labelling = 0
    if(percolation_labelling == 1) then
        write(*,*) "Percolation labelling: exact"
    else
        percolation_labelling = 0
        write(*,*) "Percolation labelling: poreblazer"
    end if
    if(par_size > 1) then
        if(percolation_labelling /= 1) call par_abort("Poreblazer's percolation labelling needs one MPI rank: "// &
            "use exact labelling (1 after the visualisation option in defaults.dat) to split the grid")
        write(*,'(a,i0,a,i0,a)') " MPI ranks: ", par_size, " (slabs of about ", ncubesz/par_size, " z planes)"
    end if

    close(1)
    close(2)
    close(3)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Initialization complete                               !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine initialize

!--------------------------------------------------------------------------------------
! Finalize: formally deallocate all arrays
!--------------------------------------------------------------------------------------                

subroutine finalize
    use atoms
    use adsorbent
    use lattice
    use distributions

    implicit none

    ! These arrays are allocated in initialize
    deallocate(aname,asigma,aeps, amass, asigma2, asigma2_he, aeps_he, asigma_n, asigma2_n)
    deallocate(adsname, coords, atype, matvec)
    deallocate(lattice_space)
    if(allocated(lattice_temp)) deallocate(lattice_temp)
    deallocate(lattice_rdist2)
    if(allocated(lattice_space_he)) deallocate(lattice_space_he)
    if(allocated(lattice_space_n)) deallocate(lattice_space_n)
    if(allocated(lattice_lj_he)) deallocate(lattice_lj_he)
    deallocate(cl_summary)
    deallocate(psd_cumul, psd)

    ! Allocated in lattice_calculations and pore_distribution
    if(allocated(n_cubes)) deallocate(n_cubes)

end subroutine finalize

!--------------------------------------------------------------------------------------
!       Now we will perform a number of preliminary calculations:
!       - distance from each lattice cube center to all atoms, detecting the shortest distance
!       - lattice accessible to a helium atom and cluster accessible to a helium atom
!       - preliminary information for helium and geometric pore volume
!       - lattice accessible to a nitrogen atom
!--------------------------------------------------------------------------------------

subroutine lattice_calculations

    use parameters
    use atoms
    use adsorbent
    use lattice
    use atomcells

    use fundcell, only:     fundamental_cell, fundcell_snglminimage, fundcell_isortho
    use vector, only:       vectype
    use parallel, only:     par_sum_i8

    implicit none

    real*8, parameter                  :: cell_target = 2.0d0       ! width of the cell-list cells (A)
    real*8, parameter                  :: margin = 1.0d-3           ! cells are searched to hicut + margin (A)
    real*8                             :: maxsigma, reach, x, xl, yl, zl, metric, u(3,3), erc
    integer*4                          :: i, j, k, l, icount, a, b, d, n, ix, iy, iz, ncand, maxnb, necand
    integer                            :: c(3), ncubes(3), nnb(3), lfirst, llast, nhe_local
    logical                            :: exact, ortho
    integer, allocatable               :: all_atoms(:), pfirst(:,:), plast(:,:), nbr(:,:), cand(:), ecand(:)
    real*8, allocatable                :: gap(:,:), d2(:)
    real*8, allocatable                :: ax(:), ay(:), az(:), aov(:), ahalf(:), as2he(:), alj(:)
    real*8, allocatable                :: cx(:), cy(:), cz(:), cov(:), chalf(:), cs2he(:), ceps(:)
    real*8, allocatable                :: ex(:), ey(:), ez(:), eov(:), ehalf(:), es2he(:), eeps(:)


    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting preliminary lattice calculations             !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    ! Each cubelet writes only its own elements of the lattice arrays, so the cubelets are
    ! independent and run in parallel. In upstream 3.0.5 the Lennard-Jones temporaries were
    ! shared between threads (a race, and slow, as every thread wrote the same variables),
    ! and the cubelet lists were filled with a counter read after other threads could have
    ! incremented it. The lists are now built afterwards, in cubelet order, so the results
    ! do not depend on the number of threads.
    !
    ! Cell list: the atoms are binned into cells about cell_target A wide and the cubelets
    ! are processed cell by cell. Each cubelet checks only the atoms in cells that can lie
    ! within the cutoff, in ascending atom order, so the Lennard-Jones sum, the overlap test
    ! and the nearest atom are exactly those found by checking every atom. When the nearest
    ! atom or the nearest surface could lie beyond the cutoff (a pore wider than about twice
    ! the cutoff), the cubelet checks every atom instead, as upstream does. The list is used
    ! in slanted coordinates, to the cutoff scaled for the cell's shape (module atomcells).
    !
    ! The distances are computed in a loop of their own over contiguous copies of the
    ! positions, with the arithmetic of fundcell_snglMinImage (and, in non-orthorhombic
    ! cells, of its multiplication by the unslant matrix), so the compiler can keep it
    ! tight; the loop that uses them runs in atom order as before.
    maxsigma = maxval(asigma)
    ortho = fundcell_isortho(fcell)
    metric = atomcells_metric()
    reach = (hicut + margin)/sqrt(max(metric, tiny(0.0d0)))   ! in slanted coordinates
    call atomcells_build(cell_target, metric > 0.05d0 .and. 0.5*maxsigma < hicut)
    u = fcell%unslantmatrix%comp
    ncubes = (/ ncubesx, ncubesy, ncubesz /)
    maxnb = maxval(ncell)
    xl = fcell%ell(1)
    yl = fcell%ell(2)
    zl = fcell%ell(3)

    ! Contiguous copies of the positions and of the type parameters the cubelet loop uses
    allocate(all_atoms(natoms), ax(natoms), ay(natoms), az(natoms), aov(natoms), ahalf(natoms), &
             as2he(natoms), alj(natoms))
    do i=1, natoms
        all_atoms(i) = i
        ax(i) = matvec(i)%comp(1)
        ay(i) = matvec(i)%comp(2)
        az(i) = matvec(i)%comp(3)
        aov(i) = 0.25*asigma2(atype(i))
        ahalf(i) = 0.5*asigma(atype(i))
        as2he(i) = asigma2_he(atype(i))
        alj(i) = aeps_he(atype(i))
    end do

    ! Range of cubelets in each cell along each side
    allocate(pfirst(0:maxnb-1, 3), plast(0:maxnb-1, 3))
    pfirst = 1
    plast = 0
    do d=1, 3
        do j=ncubes(d), 1, -1
            x = dble(j-1)*cube_size + (0.5*cube_size)
            a = max(0, min(int(x/cellw(d)), ncell(d)-1))
            pfirst(a, d) = j
            if(plast(a, d) < j) plast(a, d) = j
        end do
    end do

    ! Only this rank's planes z0..z1 (all of them on one rank). A cubelet whose answer the
    ! cutoff list cannot settle (the nearest atom or surface could lie beyond the cutoff: a
    ! pore wider than about twice the cutoff) is retried with the atoms within twice the
    ! radius, then four times, ..., until it is settled or the list holds every atom.
    ! Upstream checks every atom at once; the result is the same (lattice_cube), but in a
    ! large cell with wide pores the cost stays local. The wider list is built once per
    ! cell and kept for the cell's other cubelets.
    !$omp parallel default(shared) private(a, b, c, d, n, ix, iy, iz, j, k, l, ncand, nnb, exact, &
    !$omp nbr, gap, cand, cx, cy, cz, cov, chalf, cs2he, ceps, d2, lfirst, llast, &
    !$omp necand, ecand, ex, ey, ez, eov, ehalf, es2he, eeps, erc)
    allocate(nbr(maxnb, 3), gap(maxnb, 3), cand(natoms), d2(natoms), &
             cx(natoms), cy(natoms), cz(natoms), cov(natoms), chalf(natoms), cs2he(natoms), ceps(natoms))
    allocate(ecand(natoms), ex(natoms), ey(natoms), ez(natoms), eov(natoms), ehalf(natoms), es2he(natoms), &
             eeps(natoms))
    !$omp do collapse(3) schedule(dynamic)
    do iz=0, ncell(3)-1
        do iy=0, ncell(2)-1
            do ix=0, ncell(1)-1
                lfirst = max(pfirst(iz,3), z0)
                llast = min(plast(iz,3), z1)
                if(lfirst > llast) cycle
                c = (/ ix, iy, iz /)

                ! Candidate atoms: those in cells closer than the cutoff, in ascending order
                call cell_atoms_within(c, reach, ncand, cand, cx, cy, cz, cov, chalf, cs2he, ceps)
                erc = 0.0d0

                do l=lfirst, llast
                    do k=pfirst(iy,2), plast(iy,2)
                        do j=pfirst(ix,1), plast(ix,1)
                            call lattice_cube(j, k, l, ncand, cand, cx, cy, cz, cov, chalf, cs2he, ceps, d2, hicut, exact)
                            if(.not. exact .and. erc > 0.0d0) &
                                call lattice_cube(j, k, l, necand, ecand, ex, ey, ez, eov, ehalf, es2he, eeps, d2, erc, exact)
                            do while(.not. exact)
                                if(erc == 0.0d0) erc = hicut
                                erc = 2.0d0*erc
                                call cell_atoms_within(c, (erc + margin)/sqrt(max(metric, tiny(0.0d0))), necand, ecand, &
                                                       ex, ey, ez, eov, ehalf, es2he, eeps)
                                call lattice_cube(j, k, l, necand, ecand, ex, ey, ez, eov, ehalf, es2he, eeps, d2, erc, exact)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
    !$omp end do
    deallocate(nbr, gap, cand, d2, cx, cy, cz, cov, chalf, cs2he, ceps)
    deallocate(ecand, ex, ey, ez, eov, ehalf, es2he, eeps)
    !$omp end parallel
    deallocate(all_atoms, ax, ay, az, aov, ahalf, as2he, alj, pfirst, plast)

    ! Counts of geometrically, helium- and nitrogen-accessible cubelets (in this slab, then
    ! the whole grid), and the list of this slab's nitrogen-accessible cubelets in cubelet
    ! order (the PSD samples from it). The helium volume reads the helium mask directly.
    nhe_local = 0
    do l=z0, z1
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space(j,k,l) < 1) cycle
                ng_local = ng_local + 1
                if(lattice_space_he(j,k,l) == 1) nhe_local = nhe_local + 1
                if(lattice_space_n(j,k,l) == 1) nn_local = nn_local + 1
            end do
        end do
    end do
    ng_cubes = par_sum_i8(int(ng_local, 8))
    nhe_cubes = par_sum_i8(int(nhe_local, 8))
    nn_cubes = par_sum_i8(int(nn_local, 8))
    ! (the list of nitrogen-accessible cubelets is made after their percolation analysis)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Preliminary lattice complete                          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

contains

    ! The atoms in cells that can lie within reach (slanted) of a point in cell c, in
    ! ascending order, with their positions and type parameters: the atoms within reach
    ! times sqrt(metric) (Cartesian) of it, and perhaps some farther. Every atom when
    ! there are no cells. Gathered from the cells and sorted (upstream-order sums need
    ! ascending atoms), so the cost follows the atoms nearby, not every atom.
    subroutine cell_atoms_within(c, reach, n, list, x, y, z, ov, half, s2he, eps)
        integer, intent(in)                :: c(3)
        real*8, intent(in)                 :: reach
        integer, intent(out)               :: n
        integer, intent(out)               :: list(:)
        real*8, intent(out)                :: x(:), y(:), z(:), ov(:), half(:), s2he(:), eps(:)
        integer                            :: d, i1, i2, i3, b, m, nb(maxnb, 3), nn(3)
        real*8                             :: g(maxnb, 3)
        do d=1, 3
            call atomcells_neighbours(c(d), d, reach, nb(:,d), g(:,d), nn(d))
        end do
        n = 0
        do i3=1, nn(3)
            do i2=1, nn(2)
                do i1=1, nn(1)
                    if(g(i1,1)**2 + g(i2,2)**2 + g(i3,3)**2 >= reach*reach) cycle
                    b = 1 + nb(i1,1) + ncell(1)*(nb(i2,2) + ncell(2)*nb(i3,3))
                    do m=cell_start(b), cell_start(b+1)-1
                        n = n + 1
                        list(n) = cell_atoms(m)
                    end do
                end do
            end do
        end do
        call sort_ints(list, n)
        do m=1, n
            b = list(m)
            x(m) = ax(b)
            y(m) = ay(b)
            z(m) = az(b)
            ov(m) = aov(b)
            half(m) = ahalf(b)
            s2he(m) = as2he(b)
            eps(m) = alj(b)
        end do
    end subroutine cell_atoms_within

    ! In-place heap sort of v(1:n), ascending
    subroutine sort_ints(v, n)
        integer, intent(inout)             :: v(:)
        integer, intent(in)                :: n
        integer                            :: i, last, t, root, child, m
        do i=n, 1, -1
            ! sift down from i (building the heap), then from 1 (taking the largest off)
            if(i > n/2) cycle
            root = i
            do while(2*root <= n)
                child = 2*root
                if(child < n) then
                    if(v(child) < v(child+1)) child = child + 1
                end if
                if(v(root) >= v(child)) exit
                t = v(root); v(root) = v(child); v(child) = t
                root = child
            end do
        end do
        do last=n, 2, -1
            t = v(1); v(1) = v(last); v(last) = t
            m = last - 1
            root = 1
            do while(2*root <= m)
                child = 2*root
                if(child < m) then
                    if(v(child) < v(child+1)) child = child + 1
                end if
                if(v(root) >= v(child)) exit
                t = v(root); v(root) = v(child); v(child) = t
                root = child
            end do
        end do
    end subroutine sort_ints

    ! One cubelet, checking the atoms list(1:nlist) (ascending), whose positions and type
    ! parameters are x, y, z, ov, half, s2he and eps(1:nlist); d2 is workspace. The list
    ! holds every atom within rc (A, at least the cutoff), so the overlap test and the
    ! Lennard-Jones sum are exact. exact is false, and nothing is stored, when an atom
    ! outside the list could still be the nearest atom or give the nearest surface; with
    ! a longer list (a larger rc) the answer is the same, as the atoms added are beyond
    ! the cutoff and farther than those found.
    subroutine lattice_cube(j, k, l, nlist, list, x, y, z, ov, half, s2he, eps, d2, rc, exact)
        integer, intent(in)                :: j, k, l, nlist
        integer, intent(in)                :: list(nlist)
        real*8, intent(in)                 :: x(nlist), y(nlist), z(nlist), ov(nlist), half(nlist), s2he(nlist), eps(nlist)
        real*8, intent(inout)              :: d2(nlist)
        real*8, intent(in)                 :: rc
        logical, intent(out)               :: exact
        type(vectype)                      :: atvec1, sepvec
        real*8                             :: px, py, pz, d1, dd2, d3, q1, q2, q3, r1, r2, r3, s1, s2, s3
        real*8                             :: rdist2, rdist6, rdist12, rdist2_ref, rdist_surface, rdist_surface_ref
        real*8                             :: sig2_rdist2, lj_energy, lj_sum
        integer                            :: n, amin
        logical                            :: overlap

        ! calculate the coordinates of the center of each cubelet
        atvec1%comp = (/dble(j-1), dble(k-1), dble(l-1)/) * cube_size + (0.5 * cube_size)

        ! squared minimum-image distances, as fundcell_snglMinImage computes them; the
        ! wrap is anint(q), written so that it vectorises: truncate |q|, add 1 when the part
        ! dropped is at least 1/2, restore the sign. Exact for |q| < 2**31 (separations are
        ! within a few cell lengths), so the distances are unchanged. One loop per cell
        ! shape, so that each vectorises.
        px = atvec1%comp(1)
        py = atvec1%comp(2)
        pz = atvec1%comp(3)
        if(ortho) then
            do n=1, nlist
                d1 = px - x(n)
                dd2 = py - y(n)
                d3 = pz - z(n)
                q1 = d1/xl
                q2 = dd2/yl
                q3 = d3/zl
                r1 = dble(int(abs(q1)))
                r2 = dble(int(abs(q2)))
                r3 = dble(int(abs(q3)))
                r1 = r1 + merge(1.0d0, 0.0d0, abs(q1) - r1 >= 0.5d0)
                r2 = r2 + merge(1.0d0, 0.0d0, abs(q2) - r2 >= 0.5d0)
                r3 = r3 + merge(1.0d0, 0.0d0, abs(q3) - r3 >= 0.5d0)
                d1 = d1 - xl*(sign(r1, q1))
                dd2 = dd2 - yl*(sign(r2, q2))
                d3 = d3 - zl*(sign(r3, q3))
                d2(n) = d1*d1+dd2*dd2+d3*d3
            end do
        else
            do n=1, nlist
                d1 = px - x(n)
                dd2 = py - y(n)
                d3 = pz - z(n)
                q1 = d1/xl
                q2 = dd2/yl
                q3 = d3/zl
                r1 = dble(int(abs(q1)))
                r2 = dble(int(abs(q2)))
                r3 = dble(int(abs(q3)))
                r1 = r1 + merge(1.0d0, 0.0d0, abs(q1) - r1 >= 0.5d0)
                r2 = r2 + merge(1.0d0, 0.0d0, abs(q2) - r2 >= 0.5d0)
                r3 = r3 + merge(1.0d0, 0.0d0, abs(q3) - r3 >= 0.5d0)
                d1 = d1 - xl*(sign(r1, q1))
                dd2 = dd2 - yl*(sign(r2, q2))
                d3 = d3 - zl*(sign(r3, q3))
                ! the unslant matrix times (d1, dd2, d3), and its squared norm, summed in the
                ! order of matrix_m_mult_v and vector_getnormsq
                s1 = u(1,1)*d1 + u(1,2)*dd2 + u(1,3)*d3
                s2 = u(2,1)*d1 + u(2,2)*dd2 + u(2,3)*d3
                s3 = u(3,1)*d1 + u(3,2)*dd2 + u(3,3)*d3
                d2(n) = s1*s1+s2*s2+s3*s3
            end do
        end if

        overlap = .False.

        amin = 1
        rdist2_ref = huge(0.0d0)
        rdist_surface_ref = huge(0.0d0)
        lj_sum = 0.0d0

        do n=1, nlist                                      ! for each cubelet go through the atoms of the adsorbent structure
            rdist2 = d2(n)

            ! First check if the point is "inside" a sphere
            if(rdist2<ov(n)) then                              ! ignore the cubelet icount if it is "inside" an atom
                overlap=.True.
                exit
            end if

            if(rdist2<rdist2_ref) then                         ! in the next few lines we identify if the returned distance squared rdist2 is the smallest so far between
                rdist2_ref = rdist2                            ! cubelet icount and an atom of the structure, without an overlap
                amin = list(n)
            end if

            rdist_surface = sqrt(rdist2)- half(n)

            if(rdist_surface<rdist_surface_ref) then
                rdist_surface_ref = rdist_surface
            end if

            if(rdist2<hicut2) then                            ! ignore the atom if it is beyond the cutoff radius
                sig2_rdist2 = s2he(n)/rdist2                      ! if an atom is within the cut-off, this is a convinient place to calculate its Lennard-Jones interaction
                rdist6 = sig2_rdist2*sig2_rdist2*sig2_rdist2                         ! with a helium atom placed in the center of the cubelet icount for later use in the helium volume calculation
                rdist12 = rdist6*rdist6                           ! based in the second virial approach
                lj_energy  = eps(n)*(rdist12-rdist6)
                lj_sum = lj_sum + lj_energy
            end if
        end do

        exact = overlap .or. nlist == natoms .or. &
                (rdist2_ref < rc*rc .and. rdist_surface_ref + 0.5*maxsigma < rc)
        if(.not. exact) return

        lattice_lj_he(j,k,l) = lj_sum

        if(overlap.eqv..True.) then                          ! if an overlap was detected, this whole cubelet is ignored
            return
        end if

        lattice_space(j,k,l) = 1                          ! otherwise the cubelet is geometrically accessible

        lattice_rdist2(j,k,l) = rdist_surface_ref*rdist_surface_ref ! lattice_rdist2(j,k,l) stores the shortest squared distance between cubelet j, k, l and nearest atom (without overlap)

        if(rdist2_ref>0.25*asigma2_he(atype(amin))) then  ! accessible to a helium atom
            lattice_space_he(j,k,l) = 1
        end if

        if(rdist2_ref>asigma2_n(atype(amin))) then        ! accessible to a nitrogen atom
            lattice_space_n(j,k,l) = 1
        end if

    end subroutine lattice_cube

end subroutine lattice_calculations

!--------------------------------------------------------------------------------------
!       Here we use the information obtained in lattice_calculations
!       to generate percolated lattice accessible to helium
!--------------------------------------------------------------------------------------

subroutine helium_lattice
    use lattice
    use percolation, only: percolation_calc

    implicit none

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of helium-accessible lattice               !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    call percolation_calc(lattice_space_he, cl_summary=cl_summary)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of helium-accessible lattice complete      !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine helium_lattice

!--------------------------------------------------------------------------------------
!       Here we use the information obtained in lattice_calculations
!       to generate percolated lattice accessible to nitrogen atom/probe
!--------------------------------------------------------------------------------------

subroutine nitrogen_lattice
    use parameters
    use lattice
    use percolation, only:  percolation_calc
    use parallel, only:     par_sum_i8
    use fundcell, only:     fundcell_unslant
    use vector, only:       vectype
    use defaults, only:     radTodeg

    implicit none
    integer :: j, k, l, n

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of nitrogen-accessible lattice             !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    call percolation_calc(lattice_space_n, cl_summary=cl_summary)

    ! The list of this slab's nitrogen-accessible cubelets that percolate, in cubelet order,
    ! which the PSD samples from: the cubelets the analysis left marked (all of them if no
    ! cluster could span), as it listed them. Made only now, so that it is not held while
    ! the grid is labelled.
    nn_local = 0
    do l=z0, z1
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space_n(j,k,l) == 1) nn_local = nn_local + 1
            end do
        end do
    end do
    nn_cubes = par_sum_i8(int(nn_local, 8))
    allocate(n_cubes(max(1, nn_local)))
    n = 0
    do l=z0, z1
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space_n(j,k,l) /= 1) cycle
                n = n + 1
                n_cubes(n) = cube_local(j, k, l)
            end do
        end do
    end do

    ! vusalization options: 1 -xyz, 2 - grd, 3 - both; 0 - none
    if(vis_option == 1) then 
    call nitrogen_lattice_vis(1)
    elseif(vis_option == 2) then
    call nitrogen_lattice_vis(2)
    elseif(vis_option == 3) then
    call nitrogen_lattice_vis(1)
    call nitrogen_lattice_vis(2)
    endif

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Generation of nitrogen-accessible lattice complete    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine nitrogen_lattice

!--------------------------------------------------------------------------------------
!       In this subroutine we calculate helium-accessible volume,
!       using the second virial approach and geometric volume
!       accessible to a point probe
!--------------------------------------------------------------------------------------

subroutine helium_volume

    use parameters
    use lattice
    use adsorbent
    use fundcell, only: fundamental_cell, fundcell_getvolume
    use results
    use parallel, only: par_rank, par_size, par_send_d, par_recv_d

    implicit none
    integer ::                            i, j, k, l, cube_number
    real*8  ::                            buf(1)
    real*8  ::                            lj_energy, bf !, pore_v_he, pore_geom, volume  ! Moved these variables to the results module.

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting pore volume calculations                     !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    pore_v_he = 0.0

    ! The helium-accessible cubelets that percolate, in cubelet order: the percolation
    ! analysis leaves exactly these marked in lattice_space_he (upstream kept a list).
    ! The Boltzmann factors are computed in parallel (in place of the energies, which
    ! nothing reads afterwards), then summed in cubelet order, so the sum is upstream's:
    ! with MPI, each rank continues the running sum of the ranks below it.
    !$omp parallel do default(shared) private(j, k, l, lj_energy)
    do l=z0, z1
    do k=1, ncubesy
    do j=1, ncubesx
        if(lattice_space_he(j,k,l) /= 1) cycle
        lj_energy = lattice_lj_he(j,k,l)            ! total LJ energy of the cubelet and its environment within cut-off distance
        lj_energy = 4.0 * lj_energy
        lattice_lj_he(j,k,l) = exp(-lj_energy/temp) ! Boltzmann factor for the cubelet
    end do
    end do
    end do
    !$omp end parallel do

    buf = 0.0d0
    if(par_rank > 0) call par_recv_d(buf, 1, par_rank - 1, 21)
    pore_v_he = buf(1)
    do l=z0, z1
    do k=1, ncubesy
    do j=1, ncubesx
        if(lattice_space_he(j,k,l) /= 1) cycle
        bf = lattice_lj_he(j,k,l)
        pore_v_he = pore_v_he + bf                  ! second virial as a sum of Boltzmann factors over all cubelets
    end do
    end do
    end do
    if(par_size > 1) then
        buf(1) = pore_v_he
        if(par_rank < par_size - 1) call par_send_d(buf, 1, par_rank + 1, 21)
        if(par_rank == par_size - 1) call par_send_d(buf, 1, 0, 22)
        if(par_rank == 0) call par_recv_d(buf, 1, par_size - 1, 22)
        pore_v_he = buf(1)
    end if
    ! nothing reads the helium energies or the helium mask again
    deallocate(lattice_lj_he, lattice_space_he)

    pore_v_he = pore_v_he/dble(ntot)            ! averaging over the whole sample

    pore_geom = dble(ng_cubes)/dble(ntot)       ! geometric volume is simply proportional to all cubelets accessible to the point probe

    volume = fundcell_getvolume(fcell)          ! volume of the simulation cell

    pore_v_he = volume*pore_v_he
    pore_geom = volume*pore_geom

    write(*,'(a,f12.3)') " System volume in A^3:    ",  volume
    write(*,'(a,f12.3)') " System mass, g/mol:      ",  smass
    write(*,'(a,f12.3)') " System density, g/cm^3:  ",  rho
    write(*,*) "-------------------------------------------------------"
    write(*,*) "If calculated density is  different from the  reported "
    write(*,*) "crystallagraphic one,  make sure you remove unresolved "
    write(*,*) "atoms from the adsrobent structure file"
    write(*,*) "------------------------------------------------------"
    write(*,*)


    write(*,'(a,f12.3)') " Helium volume in A^3:    ",  pore_v_he
    write(*,'(a,f12.3)') " Helium volume in cm^3/g: ",  pore_v_he*0.602214129/smass
    write(*,*)
    write(*,'(a,f12.3)') " Geometric (point accessible) volume in A^3:    ",  pore_geom
    write(*,'(a,f12.3)') " Geometric (point accessible) volume in cm^3/g: ",  pore_geom*0.602214129/smass
    write(*,*)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Pore volume calculations complete                     !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

end subroutine helium_volume

!--------------------------------------------------------------------------------------
!       Subroutine for calculation of the accessible surface area
!--------------------------------------------------------------------------------------

subroutine surface_area

    use parameters
    use atoms
    use adsorbent
    use lattice

    use defaults, only:   rdbl, pi
    use fundcell, only:   fundamental_cell, fundcell_snglminimage, fundcell_getvolume, fundcell_maptocell, fundcell_slant, &
                          fundcell_isortho
    use atomcells
    use vector, only:     vectype
    Use random, only:     rranf
    use results
    use parallel, only:   par_sum_i4_array

    implicit none
    real*8                                :: phi, costheta, theta
    !real*8                                :: stotal, volume  ! Moved these variables to the results module.
    real*8                                :: sjreal, sfrac, stotalreduced, rdist2
    integer                               :: nx, ny, nz, ncount, i, j, k, nx_temp, ny_temp, nz_temp
    type(VecType)                         :: atvec1, atvec2, sepvec, atvec_temp
    logical                               :: deny
    integer, parameter                    :: block = 256        ! atoms whose random numbers are drawn at a time
    real*8                                :: reach, xl, yl, zl, d1, d2, d3
    integer                               :: i0, i1, ii, m, n1, n2, n3, b, maxnb, c(3), nnb(3)
    logical                               :: ortho, usecells
    real*8                                :: metric
    real*8, allocatable                   :: rnd(:,:,:), ax(:), ay(:), az(:), gap(:,:)
    integer, allocatable                  :: counts(:), nbr(:,:)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting surface area calculations                    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

 !   open(50, file="surface.xyz", status='unknown')

    if (nn_cubes == 0) then
        ! If there are no nitrogen-accessible points, the surface area is zero.
        stotal = 0.0
        volume = fundcell_getvolume(fcell)
        return
    end if

    stotal = 0.0      ! initialize cumulative accessible surface area

    ! The trial points are the same as upstream's: every trial draws its two random numbers
    ! before any test, so the numbers for a block of atoms are drawn first, in upstream's
    ! order, and the atoms then run in parallel. The surface area is summed in atom order.
    ! The overlap test checks only atoms in nearby cells (module atomcells), as no farther
    ! atom can overlap the point.
    ortho = fundcell_isortho(fcell)
    metric = atomcells_metric()
    usecells = metric > 0.05d0
    reach = (sqrt(coeff_surface2*maxval(asigma2_n)) + 1.0d-3)/sqrt(max(metric, tiny(0.0d0)))
    call atomcells_build(reach, usecells)
    maxnb = maxval(ncell)
    xl = fcell%ell(1)
    yl = fcell%ell(2)
    zl = fcell%ell(3)
    allocate(rnd(2, nsample, block), counts(block), ax(natoms), ay(natoms), az(natoms))
    do k=1, natoms
        ax(k) = matvec(k)%comp(1)
        ay(k) = matvec(k)%comp(2)
        az(k) = matvec(k)%comp(3)
    end do

    do i0=1, natoms, block
    i1 = min(natoms, i0 + block - 1)
    do i=i0, i1
        do j=1, nsample
            rnd(1, j, i-i0+1) = rranf()
            rnd(2, j, i-i0+1) = rranf()
        end do
    end do

    !$omp parallel default(shared) private(i, ii, j, k, m, n1, n2, n3, b, c, nnb, nbr, gap, ncount, phi, &
    !$omp costheta, theta, atvec1, atvec_temp, sepvec, nx, ny, nz, rdist2, deny, d1, d2, d3)
    allocate(nbr(maxnb, 3), gap(maxnb, 3))
    !$omp do schedule(dynamic)
    do i=i0, i1    ! atoms of the structure
        ii = i - i0 + 1
        ncount = 0
        do j=1, nsample   ! number of sample points for each atom
            ! generate random vector of length 1
            ! first generate phi -pi pi
            phi=pi - rnd(1, j, ii)*2.0*pi
            ! generate theta -pi:pi
            costheta = 1 - rnd(2, j, ii) * 2.0
            theta = acos(costheta)
            atvec1%comp(1) = sin(theta) * cos(phi)
            atvec1%comp(2) = sin(theta) * sin(phi)
            atvec1%comp(3) = costheta
            ! make this vector of (sigma+probe_diameter)/2.0 length
            atvec1%comp = atvec1%comp * (coeff_surface * asigma_n(atype(i))) +  coords(:,i)
            ! translate the center of the coordinate to the particle i center and apply PBC
            ! apply PBCs to ensure that the selected point is within the simulation cell
            atvec_temp = atvec1
            if(fcell%orthoflag.eqv..True.) then
            atvec1 = fundcell_maptocell(fcell,atvec1)
            else
            atvec1 = fundcell_slant(fcell, atvec1)
            end if
            ! locate the cubelet in which the point sits
            nx = int(atvec1%comp(1)/cube_size) + 1
            ny = int(atvec1%comp(2)/cube_size) + 1
            nz = int(atvec1%comp(3)/cube_size) + 1
            if(nx>ncubesx) nx = nx - (nx/ncubesx)*ncubesx
            if(nx<1) nx = nx + (1-(nx/ncubesx))*ncubesx
            if(ny>ncubesy) ny = ny - (ny/ncubesy)*ncubesy
            if(ny<1) ny = ny +  (1-(ny/ncubesy))*ncubesy
            if(nz>ncubesz) nz = nz - (nz/ncubesz)*ncubesz
            if(nz<1) nz = nz +  (1-(nz/ncubesz))*ncubesz
            if(nz < z0 .or. nz > z1) cycle        ! another rank's planes: it counts this point
            if(lattice_space_n(nx,ny,nz)<1) cycle ! reject the point if it is within an atom
            !----------------------
            ! now we perform the overlap test, i.e. we ensure that the selected point is not inside
            ! any other surf_coeff*(atom_diamter+probe_diameter)/2 distance in the system
            deny=.False.
            if(usecells) then
                do m=1, 3
                    c(m) = atomcells_index(atvec1%comp(m), m)
                    call atomcells_neighbours(c(m), m, reach, nbr(:,m), gap(:,m), nnb(m))
                end do
                do n3=1, nnb(3)
                    do n2=1, nnb(2)
                        do n1=1, nnb(1)
                            if(gap(n1,1)**2 + gap(n2,2)**2 + gap(n3,3)**2 >= reach*reach) cycle
                            b = 1 + nbr(n1,1) + ncell(1)*(nbr(n2,2) + ncell(2)*nbr(n3,3))
                            do m=cell_start(b), cell_start(b+1)-1
                                k = cell_atoms(m)
                                if(k==i) cycle
                                if(ortho) then
                                    ! as fundcell_snglMinImage computes it
                                    d1 = atvec1%comp(1) - ax(k)
                                    d2 = atvec1%comp(2) - ay(k)
                                    d3 = atvec1%comp(3) - az(k)
                                    d1 = d1 - xl*(anint(d1/xl))
                                    d2 = d2 - yl*(anint(d2/yl))
                                    d3 = d3 - zl*(anint(d3/zl))
                                    rdist2 = d1*d1+d2*d2+d3*d3
                                else
                                    Call fundcell_snglMinImage(fcell,atvec1,matvec(k),sepvec,rdist2)
                                end if
                                if(rdist2<coeff_surface2*asigma2_n(atype(k))) then
                                    deny=.True.
                                    go to 100
                                end if
                            end do
                        end do
                    end do
                end do
100             continue
            else
                do k=1, natoms
                    if(k==i) cycle
                    Call fundcell_snglMinImage(fcell,atvec1,matvec(k),sepvec,rdist2)
                    if(rdist2<coeff_surface2*asigma2_n(atype(k))) then
                        deny=.True.
                        exit
                    end if
                end do
            end if
            !----------------------
            if(deny) cycle
            ncount=ncount+1
        end do
        counts(ii) = ncount
    end do
    !$omp end do
    deallocate(nbr, gap)
    !$omp end parallel
    ! every trial point was counted by the rank holding its plane
    call par_sum_i4_array(counts, i1 - i0 + 1)

    do i=i0, i1
        ! fraction of the accessible surface area for sphere i
        sfrac=dble(counts(i-i0+1))/dble(nsample)
        ! surface area for sphere i in real units (A^2)
        sjreal=4.0*pi*coeff_surface2*asigma2_n(atype(i))*sfrac
        stotal=stotal+sjreal
    end do
    end do
    deallocate(rnd, counts, ax, ay, az)

    ! converting stotal on Surface per Volume

    volume = fundcell_getvolume(fcell)
    stotalreduced=stotal/(volume)*1.E4

    ! report results

    write(*,'(a,f12.2)') ' Accessible surface area in A^2:                 ', stotal
    write(*,'(a,f12.2)') ' Accessible surface area per volume in m^2/cm^3: ', stotalreduced
    write(*,'(a,f12.2)') ' Accessible surface area per mass in m^2/g:      ', &
    stotalreduced / ( smass/(0.6022141*volume))
    write(*,*)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Surface area calculations complete                    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

end subroutine surface_area


!--------------------------------------------------------------------------------------
!       Subroutine for calculation of the PSD
!--------------------------------------------------------------------------------------

subroutine pore_distribution

    use parameters
    use atoms
    use adsorbent
    use lattice
    use distributions
    use atomcells, only: atomcells_metric


    use fundcell, only: fundcell_init, fundamental_cell, fundcell_snglminimage, fundcell_slant, fundcell_isortho
    Use vector, only: vectype
    Use random, Only: rranf
    use parallel, only: par_rank, par_size, par_sum_i4_array, par_max_d, par_max_d_array, par_exscan_i8


    implicit none
    type(vectype)                         :: atvec1, atvec2, sepvec
    integer                               :: i, j, k, l, nx, ny, nz, nx1, ny1, nz1, icount, m, bin
    integer*8                             :: isite, nn_offset
    integer, allocatable                  :: sxyz(:,:)
    real*8, allocatable                   :: sbest(:)
    integer                               :: bz0, bz1, pass
    real*8                                :: sigma_ref, sigma2_ref, rdist2
    real*8                                :: deldis1, deldis2, deldis
    integer, parameter                    :: nsamples = 10000
    integer, allocatable                  :: bins(:)
    real*8                                :: rmax2, best, v, t, g2, xq(3), lo, hi, metric
    integer                               :: bsz, nb(3), nblocks, n, d, bx, by, bz, j1, k1, l1
    logical                               :: ortho
    real*8, allocatable                   :: blockmax(:,:,:), bkey(:)
    integer, allocatable                  :: bxs(:), bys(:), bzs(:)


    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting pore size distribution calculations          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    psd_cumul=0.0
    psd=0.0
    icount = 0

    ! Upstream sorted these radii here, for the search below and the limiting diameter; the
    ! search works by blocks of the grid and needs only the largest, and the limiting
    ! diameter reads the grid, so no sorted copy is kept (8 bytes per accessible cubelet)
    rmax2 = -huge(0.0d0)
    !$omp parallel do default(shared) private(j, k, l) reduction(max:rmax2)
    do l=z0, z1
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space(j,k,l) >= 1) rmax2 = max(rmax2, lattice_rdist2(j,k,l))
            end do
        end do
    end do
    !$omp end parallel do
    rmax2 = par_max_d(rmax2)

    if (nn_cubes == 0) then
        ! If there are no nitrogen-accessible cubes, the pore size is zero.
        ! psd() and psd_cumul() have already been initialised to zero.
        return
    end if

    ! Draw every sample site first, in the order the serial code drew them, so that the
    ! samples can run in parallel and give the same distribution as before. A site is a
    ! position in the list of nitrogen-accessible cubelets of the whole grid, which with
    ! MPI is the ranks' lists in rank order: the rank holding a site gives its cubelet.
    allocate(sxyz(3, nsamples), sbest(nsamples), bins(nsamples))
    sxyz = 0
    nn_offset = par_exscan_i8(int(nn_local, 8))
    do i = 1, nsamples
        isite = int(rranf()*dble(nn_cubes), 8) + 1 ! randomly select an available cubelet
        if(isite > nn_cubes) isite = nn_cubes
        if(isite > nn_offset .and. isite <= nn_offset + nn_local) &
            call cube_ijk(n_cubes(int(isite - nn_offset)), sxyz(1, i), sxyz(2, i), sxyz(3, i))
    end do
    call par_sum_i4_array(sxyz, 3*nsamples)
    ! The radius of each sample's own cubelet, which contains it: a first answer, shared so that every
    ! rank starts its search from it
    sbest = -1.0d0
    do i = 1, nsamples
        if(sxyz(3, i) >= z0 .and. sxyz(3, i) <= z1) sbest(i) = lattice_rdist2(sxyz(1, i), sxyz(2, i), sxyz(3, i))
    end do
    call par_max_d_array(sbest, nsamples)

    ! Each sample takes the largest sphere, among the geometrically accessible cubelets,
    ! that contains its point: upstream scans the cubelets from the largest radius down and
    ! stops at the first that contains it, which is that maximum. Here the cubelets are
    ! grouped in blocks of bsz^3, and the blocks are visited from the largest radius they
    ! hold down, stopping at the first block that cannot beat the best sphere so far. A
    ! block farther from the point than its largest radius is skipped: the distance is at
    ! least sqrt(metric) times the slanted gap (module atomcells; metric is 1 in
    ! orthorhombic cells).
    ! The containment test is upstream's, so the same sphere radius is found. With MPI,
    ! each rank finds the largest sphere among its own slab's cubelets (the blocks hold
    ! only those), and the sample's answer is the largest over the ranks.
    ortho = fundcell_isortho(fcell)
    metric = atomcells_metric()
    bsz = max(8, int(sqrt(rmax2)/(6.0d0*cube_size)) + 1)
    nb = (/ (ncubesx+bsz-1)/bsz, (ncubesy+bsz-1)/bsz, (ncubesz+bsz-1)/bsz /)
    bz0 = (z0-1)/bsz
    bz1 = (z1-1)/bsz
    allocate(blockmax(0:nb(1)-1, 0:nb(2)-1, bz0:bz1))
    blockmax = -1.0d0
    do l=z0, z1
        do k=1, ncubesy
            do j=1, ncubesx
                if(lattice_space(j,k,l)<1) cycle
                blockmax((j-1)/bsz, (k-1)/bsz, (l-1)/bsz) = max(blockmax((j-1)/bsz, (k-1)/bsz, (l-1)/bsz), &
                                                               lattice_rdist2(j,k,l))
            end do
        end do
    end do
    ! The blocks, sorted by the largest radius they hold (ascending; visited from the end)
    nblocks = nb(1)*nb(2)*(bz1-bz0+1)
    allocate(bkey(nblocks), bxs(nblocks), bys(nblocks), bzs(nblocks))
    n = 0
    do l=bz0, bz1
        do k=0, nb(2)-1
            do j=0, nb(1)-1
                n = n + 1
                bkey(n) = blockmax(j, k, l)
                bxs(n) = j
                bys(n) = k
                bzs(n) = l
            end do
        end do
    end do
    !$omp parallel default(shared)
    !$omp single
    call parallel_sort(1, nblocks, bkey, bxs, bys, bzs)
    !$omp end single
    !$omp end parallel

    ! Two passes with MPI. First, the rank holding each sample's cubelet searches its own
    ! slab; the largest sphere there is usually the answer, and the largest over the ranks is
    ! then every rank's first answer. Second, the other ranks search theirs, needing only
    ! spheres larger than that: without it, a rank whose slab cannot contain a better sphere
    ! would visit every block within reach (in wide pores, most of them). One rank does only
    ! the first pass, which is the serial search.
    do pass = 1, min(2, par_size)
    !$omp parallel do default(shared) schedule(dynamic, 16) private(i, n, d, nx, ny, nz, j1, k1, l1, &
    !$omp atvec1, atvec2, sepvec, rdist2, best, v, t, g2, xq, lo, hi, bx, by, bz)
    do i = 1, nsamples
        nx = sxyz(1, i)
        ny = sxyz(2, i)
        nz = sxyz(3, i)
        if((nz >= z0 .and. nz <= z1) .neqv. (pass == 1)) cycle

        atvec1%comp(1) = dble(nx-1)*cube_size+0.5*cube_size ! this is the center of the selected cubelet
        atvec1%comp(2) = dble(ny-1)*cube_size+0.5*cube_size
        atvec1%comp(3) = dble(nz-1)*cube_size+0.5*cube_size
        xq = atvec1%comp

        ! start from the sample's own cubelet (it contains the point)
        best = sbest(i)

        do n = nblocks, 1, -1
            if(bkey(n) <= best) exit                      ! no block left can beat it
            bx = bxs(n)
            by = bys(n)
            bz = bzs(n)
            if(metric > 0.0d0) then
                ! smallest (periodic, slanted) separation of the point and the block's cubelet
                ! centres; times metric, a lower bound on the squared distance
                g2 = 0.0d0
                do d = 1, 3
                    if(d == 1) j1 = bx
                    if(d == 2) j1 = by
                    if(d == 3) j1 = bz
                    lo = dble(j1*bsz)*cube_size + 0.5*cube_size
                    hi = dble(min((j1+1)*bsz, ncubes_of(d)) - 1)*cube_size + 0.5*cube_size
                    t = xq(d) - 0.5d0*(lo + hi)
                    t = t - fcell%ell(d)*anint(t/fcell%ell(d))
                    t = max(0.0d0, abs(t) - 0.5d0*(hi - lo))
                    g2 = g2 + t*t
                end do
                if(metric*g2 > bkey(n)*(1.0d0 + 1.0d-9) + 1.0d-9) cycle
            end if
            do l1 = max(bz*bsz+1, z0), min((bz+1)*bsz, ncubesz, z1)
                do k1 = by*bsz+1, min((by+1)*bsz, ncubesy)
                    do j1 = bx*bsz+1, min((bx+1)*bsz, ncubesx)
                        if(lattice_space(j1,k1,l1)<1) cycle
                        v = lattice_rdist2(j1,k1,l1)
                        if(v <= best) cycle
                        atvec2%comp(1) = dble(j1-1)*cube_size+0.5*cube_size
                        atvec2%comp(2) = dble(k1-1)*cube_size+0.5*cube_size
                        atvec2%comp(3) = dble(l1-1)*cube_size+0.5*cube_size
                        call fundcell_snglminimage(fcell,atvec1,atvec2,sepvec,rdist2)
                        if(rdist2>v) cycle                ! the point is outside this sphere
                        best = v
                    end do
                end do
            end do
        end do
        sbest(i) = best
    end do
    !$omp end parallel do
    call par_max_d_array(sbest, nsamples)
    end do
    deallocate(blockmax, bkey, bxs, bys, bzs)

    do i = 1, nsamples
        sigma2_ref = sbest(i)
        sigma_ref = sqrt(sigma2_ref)                   ! sigma here is distance, not diameter
        bin=int(2.0*sigma_ref/binsize)+1               ! distribution bin (2 is needed to make sigma proper diameter)
        if (bin > ubound(psd_cumul, 1)) bin = ubound(psd_cumul, 1)
        bins(i) = bin
    end do

    ! Update the cumulative distribution in sample order
    do i = 1, nsamples
        do m=1, bins(i)
            psd_cumul(m)=psd_cumul(m)+1
        end do
    end do
    deallocate(sxyz, sbest, bins)
    deallocate(n_cubes, lattice_space_n)   ! nothing reads them again

    ! rank 0 writes the files (every rank has the same distribution)
    if(par_rank == 0) then
    open(13, file='psd_cumulative.txt', status='unknown') ! The file for cumulative Vp function
    open(14, file='psd.txt', status='unknown')            ! File containing pore size distribution
    else
    open(13, file='/dev/null', status='old', action='write')
    open(14, file='/dev/null', status='old', action='write')
    end if

    write(13, *) '# Cumulative accessible volume distribution as a function of probe diameter'
    write(13, *) '# '
    write(13, *) '# d(probe)                Volume Fraction'

    psd_cumul=psd_cumul/psd_cumul(1)

    do i=1, nbins
        write(13,*) binsize*real(i-1)-binsize/2.0, psd_cumul(i)
    end do


    psd(1)=0.0
    psd(2)=0.0
    psd(nbins)=0.0

    do i=2, nbins-1               ! numerical differentiation where for calculation in point i, we use points i-1 and i+1
        deldis1=psd_cumul(i+1)
        deldis2=psd_cumul(i-1)
        deldis=deldis1-deldis2
        psd(i)=-1.0*(deldis)/(binsize*2.0) ! note that this way of taking derivatives is very crude
    end do

    write(14,*) '# Derivative distribution function -dV(r)/dr (or -dV(d)/dd) vs d'

    do i=2, nbins-1
        write(14,*) binsize*real(i-1)-binsize/2.0, psd(i)
    end do

    close(13)
    close(14)

    write(*,*)  "Cumulative PSD and differential PSD "
    write(*,*)  "have been stored in files psd_cumulative.txt"
    write(*,*)  "and psd.txt"

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Pore size distribution calculations complete          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

contains

    integer function ncubes_of(d)
        integer, intent(in) :: d
        if(d == 1) then
            ncubes_of = ncubesx
        else if(d == 2) then
            ncubes_of = ncubesy
        else
            ncubes_of = ncubesz
        end if
    end function ncubes_of

end subroutine pore_distribution

!--------------------------------------------------------------------------------------
!       Subroutine for for limiting pore diameter and for maximum pore size
!--------------------------------------------------------------------------------------

subroutine limiting_diameter

    use parameters
    use atoms
    use lattice

    use percolation, only: percolation_calc_simple
    use results
    use parallel, only: par_max_d, par_min_d, par_allgather_d

    implicit none
    real*8  :: rhigh, rlow, rmiddle, rmax, rdiff, rdiff_old
    integer :: i, j, percolation_type, i1, k1, l1
    integer :: nlow
    real*8  :: lmax, lmin, v, low2(2), top1, top2, t1, t2, vabove
    real*8, allocatable :: lows(:)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting limiting pore diameter and maximum pore size !"
    write(*,*) "! analysis                                              !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    ! Upstream bisected on a sorted list of the radii. The bisection needs only the largest
    ! and smallest radius, and at each probe radius rmiddle the radii next to it in that list:
    ! the second largest below rmiddle (PA1(j-1)) and the smallest not below (PA1(j+1)). These
    ! are read from the grid as each step marks it, on every rank, with no sorted copy.
    allocate(lattice_temp(ncubesx, ncubesy, z0:z1))
    lmax = -huge(0.0d0)
    lmin = huge(0.0d0)
    !$omp parallel do default(shared) private(i1, k1, l1) reduction(max:lmax) reduction(min:lmin)
    do l1=z0, z1
        do k1=1, ncubesy
            do i1=1, ncubesx
                if(lattice_space(i1,k1,l1) < 1) cycle
                lmax = max(lmax, lattice_rdist2(i1,k1,l1))
                lmin = min(lmin, lattice_rdist2(i1,k1,l1))
            end do
        end do
    end do
    !$omp end parallel do
    rmax  = par_max_d(lmax)
    rhigh = rmax
    rlow  = par_min_d(lmin)
    rdiff = 0.0

    do i=1, 100000 ! here we employ bisection method + percolation analysis to find the largest probe
        ! with respect to which the system remains lattice-percolated

        rdiff_old = rdiff
        rdiff = rhigh-rlow
        if(abs(rdiff)<0.25.or.(rdiff==rdiff_old)) exit ! convergence is achieved
        rmiddle = 0.5*(rhigh+rlow)

        ! the geometrically accessible cubelets with a radius of at least rmiddle (upstream
        ! marked the sorted entries above j, the last radius below rmiddle: the same cubelets),
        ! with the two largest radii below rmiddle and the smallest not below it
        top1 = -huge(0.0d0)
        top2 = -huge(0.0d0)
        vabove = huge(0.0d0)
        !$omp parallel default(shared) private(i1, k1, l1, v, t1, t2)
        t1 = -huge(0.0d0)
        t2 = -huge(0.0d0)
        !$omp do reduction(min:vabove)
        do l1=z0, z1
            do k1=1, ncubesy
                do i1=1, ncubesx
                    lattice_temp(i1,k1,l1) = 0
                    if(lattice_space(i1,k1,l1) < 1) cycle
                    v = lattice_rdist2(i1,k1,l1)
                    if(.not. (v < rmiddle)) then
                        lattice_temp(i1,k1,l1) = 1
                        vabove = min(vabove, v)
                    else if(v > t1) then
                        t2 = t1
                        t1 = v
                    else if(v > t2) then
                        t2 = v
                    end if
                end do
            end do
        end do
        !$omp end do
        !$omp critical
        call keep_two_largest(t1, top1, top2)
        call keep_two_largest(t2, top1, top2)
        !$omp end critical
        !$omp end parallel

        call percolation_calc_simple(lattice_temp, cl_summary, spanning)

        if(spanning>0) then
            ! PA1(j-1): the second largest radius below rmiddle over the ranks, counting
            ! repeats (the smallest radius if there is only one below)
            nlow = 0
            if(top1 > -huge(0.0d0)) then
                nlow = nlow + 1
                low2(nlow) = top1
            end if
            if(top2 > -huge(0.0d0)) then
                nlow = nlow + 1
                low2(nlow) = top2
            end if
            call par_allgather_d(low2(1:nlow), lows)
            if(size(lows) >= 2) then
                rlow = second_largest(lows)
            else
                rlow = par_min_d(lmin)
            end if
            deallocate(lows)
            percolation_type = spanning
            cycle
        else
            ! PA1(j+1): the smallest radius not below rmiddle (the largest if there is none)
            rhigh = par_min_d(vabove)
            if(rhigh == huge(0.0d0)) rhigh = rmax
            cycle
        end if

    end do
    deallocate(lattice_temp)

    write(*,'(a,f12.2)') ' Pore limiting diameter in A: ', 2.0*sqrt(rhigh)
    write(*,'(a,f12.2)') ' Maximum pore diameter in A:  ', 2.0*sqrt(rmax)
    if(percolation_type==1) then
        write(*, *) ' The system is percolated in ', percolation_type, ' dimension (channels)'
    elseif(percolation_type==2) then
        write(*, *) ' The system is percolated in ', percolation_type, ' dimensions (slits)'
    elseif(percolation_type==3) then
        write(*, *) 'The system is percolated in ', percolation_type, ' dimensions (3D pores)'
    else
    end if

    ! Store the output in the results module.
    pore_lim = 2.0*sqrt(rhigh)
    pore_max = 2.0*sqrt(rmax)
    sys_perc = percolation_type

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Limiting pore diameter and maximum pore size          !"
    write(*,*) "! analysis complete                                     !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return

contains

    ! Keep a and b the largest and second largest of the values offered so far (with repeats)
    subroutine keep_two_largest(x, a, b)
        real*8, intent(in)    :: x
        real*8, intent(inout) :: a, b
        if(x > a) then
            b = a
            a = x
        else if(x > b) then
            b = x
        end if
    end subroutine keep_two_largest

    ! The second largest of v, counting repeats
    real*8 function second_largest(v)
        real*8, intent(in) :: v(:)
        real*8             :: a, b
        integer            :: n
        a = -huge(0.0d0)
        b = -huge(0.0d0)
        do n=1, size(v)
            if(v(n) > a) then
                b = a
                a = v(n)
            else if(v(n) > b) then
                b = v(n)
            end if
        end do
        second_largest = b
    end function second_largest

end subroutine limiting_diameter

!--------------------------------------------------------------------------------------
!       Here we visualize the nitrogen accessible network
!       using two options: xyz format or grd format
!--------------------------------------------------------------------------------------

subroutine nitrogen_lattice_vis(option)
    use parameters
    use lattice
    use adsorbent
    use percolation, only:  percolation_calc
    use fundcell, only:     fundcell_unslant
    use vector, only:       vectype
    use defaults, only:     radTodeg
    use parallel, only:     par_rank, par_size, par_owner, par_send_i8, par_recv_i8, par_send_r4, par_recv_r4

    implicit none
    integer                            :: option
    integer                            :: nx, ny, nz, nxlow, nxup, nylow, nyup, nzlow, nzup, i, j, k, r
    integer                            :: ncx, ncy, ncz
    real                               :: x, y, z, field
    type(vectype)                      :: atvec1
    real, allocatable                  :: plane(:,:)
    integer*8                          :: cnt(1)
    integer*8, allocatable             :: cubes(:)

    ! Rank 0 writes both files; with MPI the other ranks send it their parts in order
    if(option==1) then
    ! Creating a grid of lattice sites accessible to nitrogen probe in xyz fornat
    if(par_rank > 0) then
        cnt(1) = nn_local
        call par_send_i8(cnt, 1, 0, 31)
        allocate(cubes(max(1, nn_local)))
        do i=1, nn_local
            cubes(i) = n_cubes(i)
        end do
        if(nn_local > 0) call par_send_i8(cubes, nn_local, 0, 32)
        deallocate(cubes)
        return
    end if
    open(110, file='nitrogen_network.xyz', status='unknown')

    if(nn_cubes <= huge(0)) then
        write(110,*) int(nn_cubes)
    else
        write(110,*) nn_cubes
    end if
    write(110,*)

    do r=0, par_size-1
        if(r == 0) then
            allocate(cubes(max(1, nn_local)))
            cnt(1) = nn_local
            do i=1, nn_local
                cubes(i) = n_cubes(i)
            end do
        else
            call par_recv_i8(cnt, 1, r, 31)
            allocate(cubes(max(1_8, cnt(1))))
            if(cnt(1) > 0) call par_recv_i8(cubes, int(cnt(1)), r, 32)
        end if
        do i=1, int(cnt(1))
        ! the local cube number on rank r: j fastest, then k, then the plane from its first
        nx = int(mod(cubes(i) - 1, int(ncubesx, 8))) + 1
        ny = int(mod((cubes(i) - 1)/ncubesx, int(ncubesy, 8))) + 1
        nz = int((cubes(i) - 1)/(int(ncubesx, 8)*ncubesy)) + slab_first(r)

        x = dble(nx-1)*cube_size+0.5*cube_size ! this is the center of the selected cubelet
        y = dble(ny-1)*cube_size+0.5*cube_size
        z = dble(nz-1)*cube_size+0.5*cube_size

        x = x +  fcell%mx
        y = y +  fcell%my
        z = z +  fcell%mz

        atvec1%comp(1) = x
        atvec1%comp(2) = y
        atvec1%comp(3) = z

        atvec1 = fundcell_unslant(fcell, atvec1)

        write(110,*) "N ", atvec1%comp(1), atvec1%comp(2), atvec1%comp(3)
        end do
        deallocate(cubes)
    end do
    close(110)
    return
    endif

    if(option==2) then
   allocate(plane(ncubesx, ncubesy))
   if(par_rank > 0) then
       do k=z0, z1
           call fill_plane(k)
           call par_send_r4(plane, ncubesx*ncubesy, 0, 33)
       end do
       deallocate(plane)
       return
   end if
    open(111, file='nitrogen_network.grd', status='unknown')

    write(111,*) "Poreblazer nitrogen accessible grid: field = r"
    write(111,*) "(1p,e12.5)"
    write(111,'(6f12.3)') fcell%ell(1), fcell%ell(2), fcell%ell(3), radTodeg*fcell%anglebc, radTodeg*fcell%angleac, &
    radTodeg*fcell%angleab

    ncx = ncubesx-1
    ncy = ncubesy-1
    ncz = ncubesz-1  

    write(111,*) ncx, ncy, ncz   

    write(111,*) "1 ",  " 0 ", ncx, " 0 ", ncy, " 0 ", ncz

   ! One write statement per plane (one record per value by format reversion) instead of
   ! one per cubelet: the same file, about 30% faster. The formatting itself is most of
   ! what is left, and libgfortran serialises it, so formatting on threads is slower.
   do k=1,  ncubesz
    r = 0
    if(par_size > 1) r = par_owner(ncubesz, k)
    if(r == 0) then
        call fill_plane(k)
    else
        call par_recv_r4(plane, ncubesx*ncubesy, r, 33)
    end if
    write(111,'(e12.5)') plane
   end do
   deallocate(plane)
    close(111)
    return
    end if

    return

contains

    ! The first plane of rank r's slab
    integer function slab_first(r)
        use parallel, only: par_slab
        integer, intent(in) :: r
        integer             :: a, b
        call par_slab(ncubesz, r, a, b)
        slab_first = a
    end function slab_first

    ! plane: the pore radius of every cubelet of plane k (one of this rank's)
    subroutine fill_plane(k)
        integer, intent(in) :: k
        integer             :: i, j
        do j=1, ncubesy
            do i=1, ncubesx
                if(lattice_rdist2(i,j,k) == 0.0) then
                    plane(i,j) = 0.0
                else
                    plane(i,j) = sqrt(lattice_rdist2(i,j,k))
                endif
            end do
        end do
    end subroutine fill_plane

end subroutine nitrogen_lattice_vis

!==============================================================

!==============================================================
! Quicksort of arr(lo:hi) in ascending order, carrying brr, crr and drr along, with the
! two halves of each partition sorted as OpenMP tasks. Ranges below cutoff use the serial
! sort. Call from inside a parallel region's single construct.

recursive subroutine parallel_sort(lo, hi, arr, brr, crr, drr)
    implicit none
    integer, intent(in)       :: lo, hi
    real*8                    :: arr(*)
    integer                   :: brr(*), crr(*), drr(*)
    integer, parameter        :: cutoff = 100000
    integer                   :: i, j, tempi
    real*8                    :: pivot, temp

    if(hi - lo < cutoff) then
        if(hi > lo) call sort(hi-lo+1, arr(lo), brr(lo), crr(lo), drr(lo))
        return
    end if

    ! Hoare partition around the median of the first, middle and last values
    pivot = max(min(arr(lo), arr(hi)), min(max(arr(lo), arr(hi)), arr((lo+hi)/2)))
    i = lo - 1
    j = hi + 1
    do
        do
            i = i + 1
            if(arr(i) >= pivot) exit
        end do
        do
            j = j - 1
            if(arr(j) <= pivot) exit
        end do
        if(i >= j) exit
        temp = arr(i);  arr(i) = arr(j);  arr(j) = temp
        tempi = brr(i); brr(i) = brr(j); brr(j) = tempi
        tempi = crr(i); crr(i) = crr(j); crr(j) = tempi
        tempi = drr(i); drr(i) = drr(j); drr(j) = tempi
    end do

    !$omp task default(shared) firstprivate(lo, j)
    call parallel_sort(lo, j, arr, brr, crr, drr)
    !$omp end task
    !$omp task default(shared) firstprivate(hi, j)
    call parallel_sort(j+1, hi, arr, brr, crr, drr)
    !$omp end task
    !$omp taskwait

end subroutine parallel_sort

!==============================================================
! Ascending sort of arr(lo:hi) (values only): quicksort with the two halves of each
! partition as OpenMP tasks, insertion sort for short ranges. Call from inside a parallel
! region's single construct.

recursive subroutine parallel_sort_values(lo, hi, arr)
    implicit none
    integer, intent(in)       :: lo, hi
    real*8                    :: arr(*)
    integer, parameter        :: cutoff = 100000, small = 16
    integer                   :: i, j
    real*8                    :: pivot, temp

    if(hi - lo < small) then
        do i=lo+1, hi
            temp = arr(i)
            j = i - 1
            do while(j >= lo)
                if(arr(j) <= temp) exit
                arr(j+1) = arr(j)
                j = j - 1
            end do
            arr(j+1) = temp
        end do
        return
    end if

    ! Hoare partition around the median of the first, middle and last values
    pivot = max(min(arr(lo), arr(hi)), min(max(arr(lo), arr(hi)), arr((lo+hi)/2)))
    i = lo - 1
    j = hi + 1
    do
        do
            i = i + 1
            if(arr(i) >= pivot) exit
        end do
        do
            j = j - 1
            if(arr(j) <= pivot) exit
        end do
        if(i >= j) exit
        temp = arr(i); arr(i) = arr(j); arr(j) = temp
    end do

    if(hi - lo < cutoff) then
        call parallel_sort_values(lo, j, arr)
        call parallel_sort_values(j+1, hi, arr)
    else
        !$omp task default(shared) firstprivate(lo, j)
        call parallel_sort_values(lo, j, arr)
        !$omp end task
        !$omp task default(shared) firstprivate(hi, j)
        call parallel_sort_values(j+1, hi, arr)
        !$omp end task
        !$omp taskwait
    end if

end subroutine parallel_sort_values

subroutine sort(n,arr,brr,crr,drr)
    integer n,m,nstack
    real*8 arr(n)
    integer brr(n),crr(n), drr(n)
    parameter (M=7,NSTACK=50)
    integer i,ir,j,jstack,k,l,istack(NSTACK)
    real*8 a,temp
    integer b, c, d, tempi
    jstack=0
    l=1
    ir=n
    1     if(ir-l.lt.M)then
    do 12 j=l+1,ir
    a=arr(j)
    b=brr(j)
    c=crr(j)
    d=drr(j)
    do 11 i=j-1,1,-1
    if(arr(i).le.a)goto 2
    arr(i+1)=arr(i)
    brr(i+1)=brr(i)
    crr(i+1)=crr(i)
    drr(i+1)=drr(i)
    11        continue
    i=0
    2         arr(i+1)=a
    brr(i+1)=b
    crr(i+1)=c
    drr(i+1)=d
    12      continue
    if(jstack.eq.0)return
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
    else
    k=(l+ir)/2
    temp=arr(k)
    arr(k)=arr(l+1)
    arr(l+1)=temp
    tempi=brr(k)
    brr(k)=brr(l+1)
    brr(l+1)=tempi
    tempi=crr(k)
    crr(k)=crr(l+1)
    crr(l+1)=tempi
    tempi=drr(k)
    drr(k)=drr(l+1)
    drr(l+1)=tempi
    if(arr(l+1).gt.arr(ir))then
    temp=arr(l+1)
    arr(l+1)=arr(ir)
    arr(ir)=temp
    tempi=brr(l+1)
    brr(l+1)=brr(ir)
    brr(ir)=tempi
    tempi=crr(l+1)
    crr(l+1)=crr(ir)
    crr(ir)=tempi
    tempi=drr(l+1)
    drr(l+1)=drr(ir)
    drr(ir)=tempi
    endif
    if(arr(l).gt.arr(ir))then
    temp=arr(l)
    arr(l)=arr(ir)
    arr(ir)=temp
    tempi=brr(l)
    brr(l)=brr(ir)
    brr(ir)=tempi
    tempi=crr(l)
    crr(l)=crr(ir)
    crr(ir)=tempi
    tempi=drr(l)
    drr(l)=drr(ir)
    drr(ir)=tempi
    endif
    if(arr(l+1).gt.arr(l))then
    temp=arr(l+1)
    arr(l+1)=arr(l)
    arr(l)=temp
    tempi=brr(l+1)
    brr(l+1)=brr(l)
    brr(l)=tempi
    tempi=crr(l+1)
    crr(l+1)=crr(l)
    crr(l)=tempi
    tempi=drr(l+1)
    drr(l+1)=drr(l)
    drr(l)=tempi
    endif
    i=l+1
    j=ir
    a=arr(l)
    b=brr(l)
    c=crr(l)
    d=drr(l)
    3       continue
    i=i+1
    if(arr(i).lt.a)goto 3
    4       continue
    j=j-1
    if(arr(j).gt.a)goto 4
    if(j.lt.i)goto 5
    temp=arr(i)
    arr(i)=arr(j)
    arr(j)=temp
    tempi=brr(i)
    brr(i)=brr(j)
    brr(j)=tempi
    tempi=crr(i)
    crr(i)=crr(j)
    crr(j)=tempi
    tempi=drr(i)
    drr(i)=drr(j)
    drr(j)=tempi
    goto 3
    5       arr(l)=arr(j)
    arr(j)=a
    brr(l)=brr(j)
    brr(j)=b
    crr(l)=crr(j)
    crr(j)=c
    drr(l)=drr(j)
    drr(j)=d
    jstack=jstack+2
    if(jstack.gt.NSTACK)pause 'NSTACK too small in sort2'
    if(ir-i+1.ge.j-l)then
    istack(jstack)=ir
    istack(jstack-1)=i
    ir=j-1
    else
    istack(jstack)=j-1
    istack(jstack-1)=l
    l=i
    endif
    endif
    goto 1
end subroutine sort
!  (C) Copr. 1986-92 Numerical Recipes Software .2.
