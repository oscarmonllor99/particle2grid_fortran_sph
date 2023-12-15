! Created on Mon Nov 14 2022
!author: ÓSCAR MONLLOR BERBEGAL
!Basic idea: turn particle field to continuum with a SPH kernel

module SPH
    implicit none

contains
    ! RECURSIVE SORTING by https://gist.github.com/1AdAstra1/6f7785373efe5bb6c254d2e20c78ccc4
    recursive subroutine quicksort(a)
    implicit none
    real :: a(:)
    real x, t
    integer :: first = 1, last
    integer i, j

    last = size(a, 1)
    x = a( (first+last) / 2 )
    i = first
    j = last
    
    do
        do while (a(i) < x)
            i=i+1
        end do
        do while (x < a(j))
            j=j-1
        end do
        if (i >= j) exit
        t = a(i);  a(i) = a(j);  a(j) = t
        i=i+1
        j=j-1
    end do
    
    if (first < i - 1) call quicksort(a(first : i - 1))
    if (j + 1 < last)  call quicksort(a(j + 1 : last))
    end subroutine quicksort



    subroutine sorting_sph(id, x_pos, y_pos, z_pos, hpart, kneigh, clean_neighbours, num_neigh)
        !Para la particula "id", calculamos la distancia a sus vecinos (clean_neighbours)
        !Después ordenamos el array de distancias y tomamos la posición "kneigh", cuyo valor es "h"
        implicit none
        integer :: kneigh
        real, dimension(:) :: x_pos, y_pos, z_pos, hpart
        integer, dimension(:) :: clean_neighbours
        integer :: num_neigh

        real, dimension(num_neigh) :: distance
        integer :: id, ip0, ip2
        do ip0=1,num_neigh
            ip2 = clean_neighbours(ip0)
            distance(ip0) = sqrt( (x_pos(id) - x_pos(ip2))**2 &
                                + (y_pos(id) - y_pos(ip2))**2 &
                                + (z_pos(id) - z_pos(ip2))**2 )
        enddo

        call quicksort(distance)
        hpart(id) = distance(kneigh)
    end subroutine


    subroutine sph_lenght(partNum, x_pos, y_pos, z_pos, Lx, Ly, Lz, kneigh, hpart)
        use omp_lib
        implicit none
        
        integer :: partNum, kneigh
        real, dimension(partNum) :: x_pos, y_pos, z_pos, hpart
        real :: Lx, Ly, Lz

        integer :: max_num_part  !maximum number of particles per cell
        integer :: max_neigh = 4096
        integer :: N_aux = 128 ! spacing of auxiliar grid to find neighbours quickly
        integer :: x_cell, y_cell, z_cell
        integer, dimension(:,:,:,:), allocatable :: grid_aux
        integer, dimension(:,:,:), allocatable :: cell_part_num
        integer, dimension(:,:), allocatable :: part_pos
        integer, dimension(:), allocatable :: neighbours
        integer :: zero_at, num_neigh
        real :: resx , resy, resz
        integer :: ip, ipp, ix, iy, iz
        integer :: thread_id
        integer :: next_cell
        integer :: warning_done = 0
        integer, dimension(:), allocatable :: cell_particles
        integer :: cell_space

        !f2py intent(in) partNum, x_pos, y_pos, z_pos, Lx, Ly, Lz, kneigh, 
        !f2py intent(out) hpart
        !f2py depend(partNum) x_pos, y_pos, z_pos, hpart

        max_num_part = 2*kneigh !maximum number of particles to take into account inside a cell to calculate kneigh
                                ! I assume that more than 2*kneigh is not necessary, as thus h would already be less than the grid resolution, 
                                ! so h smaller (putting more particles in the cell) would give the same result
                                ! this is to avoid grid_aux too large.

        !allocate subroutine arrays
        allocate(grid_aux(N_aux, N_aux, N_aux, max_num_part))
        allocate(cell_part_num(N_aux, N_aux, N_aux))
        allocate(part_pos(partNum, 3))
        allocate(neighbours(max_neigh))
        allocate(cell_particles(max_num_part))
        grid_aux(:,:,:,:) = 0
        cell_part_num(:,:,:) = 0
        part_pos(:,:)= 0

        ! Putting particles in an auxiliar grid
        ! Se crea la malla auxiliar para poder localizar de forma rápida los vecinos de cada partícula, 
        ! sobre los cuales calcular h
        ! En grid_aux se guarda qué partículas hay en qué celda, hasta un máximo de max_num_part,
        ! para evitar arrays demasiado grandes
        ! En cell_part_num el número de partículas de cada celda

        ! write(*,*) ' --> Creating auxiliar grid to find neighbours'
        resx = Lx/N_aux
        resy = Ly/N_aux
        resz = Lz/N_aux
        do ip=1,partNum
            x_cell = int(x_pos(ip)/resx) + 1
            y_cell = int(y_pos(ip)/resy) + 1
            z_cell = int(z_pos(ip)/resz) + 1
            part_pos(ip,1) = x_cell
            part_pos(ip,2) = y_cell
            part_pos(ip,3) = z_cell
            zero_at = cell_part_num(x_cell, y_cell, z_cell) + 1
            if (zero_at <= max_num_part) then
                grid_aux(x_cell, y_cell, z_cell, zero_at) = ip
                cell_part_num(x_cell, y_cell, z_cell) = cell_part_num(x_cell, y_cell, z_cell) + 1
            else 
                if (warning_done == 0) then
                ! print*,'WARNING! La partícula',ip,'no cabe en su celda, se necesita una malla más fina o un max_num_part mayor'
                ! print*,'Warning en el', real(ip)/real(partNum)*100,'%'
                warning_done = 1
                endif
            endif
        enddo

        ! write(*,*) '... done'
        ! write(*,*) 'CHECK: max in cell_part_num', maxval(cell_part_num)

        warning_done = 0
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! calculating h for each cell (I assume all particle sin a cell share the same h)
        ! Ahora se calculan primero los vecinos de cada partícula, que se van colocando en
        ! el array neighbours. Si hay menos de "kneigh" vecinos en la propia celda de la partícula,
        ! se mira en las celdas adyacentes hasta que en "neighbours" hay más partículas que "kneigh"
        ! luego con las partículas de "neighbours" se calcula la "kneigh" (normalmente la 32) más lejana.

        ! write(*,*) ' --> Calculating h for each particle'
        
        !$OMP PARALLEL SHARED(partNum, part_pos, grid_aux, cell_part_num, max_neigh, kneigh, hpart, N_aux), &
        !$OMP          PRIVATE(ip, ipp, thread_id, x_cell, y_cell, z_cell, next_cell, zero_at, &
        !$OMP                    cell_particles, cell_space, warning_done, neighbours, &
        !$OMP                    ix, iy, iz, num_neigh)
        !$OMP DO REDUCTION(+:hpart)
        do ip=1,partNum
            thread_id = OMP_get_thread_num()
            x_cell = part_pos(ip, 1)
            y_cell = part_pos(ip, 2)
            z_cell = part_pos(ip, 3)
            next_cell = 0 !see the 27 cells sorrounding x_cell, y_cell, z_cell. If there is more than 32 neigh
                                ! I stop, else next_cell += 1 and continue searching
            neighbours(:) = 0 !in each position, a ip neighbour
            zero_at = 1 !first 0 in neighbours, that is, where to put the new neighbours found
            cell_particles(:) = grid_aux(x_cell, y_cell, z_cell, :) !particles in this cell
            cell_space = cell_part_num(x_cell, y_cell, z_cell) !space needed to put cell particles inside neighbours
            if (cell_space > max_neigh-zero_at) then !max_neigh-zero_at = remaining space
                cell_space = max_neigh-zero_at
                if (warning_done == 0) then
                    warning_done = 1
                endif
            endif

            do ipp=1,cell_space
                neighbours(ipp+zero_at-1) = cell_particles(ipp)
            enddo
            zero_at = zero_at + cell_space
                
            outer: do while (zero_at <= kneigh) 
            !look at sorrounding cells for neighbours
                do iz=-1-next_cell, 1+next_cell
                if (z_cell+iz > 0 .and. iz+z_cell < N_aux+1) then
                    do iy=-1-next_cell, 1+next_cell
                    if (iy+y_cell > 0 .and. iy+y_cell < N_aux+1) then
                        do ix=-1-next_cell, 1+next_cell
                        if (ix+x_cell > 0 .and. ix+x_cell < N_aux+1) then
                            if ( (abs(ix) > next_cell) .or. &
                                 (abs(iy) > next_cell) .or. &
                                 (abs(iz) > next_cell)) then
                                cell_particles(:) = grid_aux(x_cell+ix, y_cell+iy, z_cell+iz, :) !repeat the process made for the particle cell
                                cell_space = cell_part_num(x_cell+ix, y_cell+iy, z_cell+iz)
                                if (cell_space > max_neigh-zero_at) then
                                    cell_space = max_neigh-zero_at
                                endif
                                do ipp=1,cell_space
                                    neighbours(ipp+zero_at-1) = cell_particles(ipp)
                                enddo
                                zero_at = zero_at + cell_space !move zero_at

                                !si ya hay más partículas que "kneigh" podria no mirar más, aunque esto 
                                !puede hacer h más grande de forma artificial, pero acelera el proceso
                                if (zero_at>kneigh) exit outer

                            endif
                        endif
                        enddo
                    endif
                    enddo
                endif
                enddo
                next_cell = next_cell + 1  !look at the next cells if neighbours found are less than kneigh
            enddo outer
                
            num_neigh = max(kneigh, zero_at - 1) !number of neighbours collected to find h
            ! Now sort neighbours for ip particle and find h
            call sorting_sph(ip, x_pos, y_pos, z_pos, hpart, kneigh, neighbours(:num_neigh), num_neigh)

        enddo ! h calculated for all particles
        !$OMP END DO
        !$OMP END PARALLEL

        ! if (warning_done == 1) then
        !     write(*,*) 'WARNING: Maximum number of neighbours exceeded by a particle in its own cell'
        ! endif
        ! write(*,*) '... done'

        !deallocation
        deallocate(grid_aux)
        deallocate(cell_part_num)
        deallocate(part_pos)
        deallocate(neighbours)
        deallocate(cell_particles)
    end subroutine

    subroutine sph_kernel(r, h, W)
        implicit none
        real :: q,r,h,W
        q = r/h
        if(q<1) then
            W = 1 - 1.5*q**2*(1-0.5*q)
        else if ( q >= 1 .and. q < 2 ) then
            W = 0.25*(2-q)**3
        endif
    end subroutine

    subroutine sph_density(partNum, x_pos, y_pos, z_pos, field, nx, ny, nz, Lx, Ly, Lz, hpart, grid_field)
        use omp_lib
        implicit none
        integer :: partNum, nx, ny, nz
        real, dimension(:) :: x_pos, y_pos, z_pos, field
        real :: Lx, Ly, Lz
        real :: resx , resy, resz
        real :: dx2 , dy2, dz2
        real :: r, W, norm
        real, dimension(:) :: hpart
        real, dimension(:,:,:) :: grid_field
        real, dimension(:,:,:), allocatable :: contribution
        integer :: ip, ix, iy, iz, ix2, iy2, iz2
        integer :: thread_id
        integer :: x_cell, y_cell, z_cell
        integer :: until_cell_x, until_cell_y, until_cell_z
        integer, dimension(:,:), allocatable :: part_pos

        allocate(part_pos(partNum,3))
        part_pos(:,:) = 0

        resx = Lx/nx
        resy = Ly/ny
        resz = Lz/nz

        !$OMP PARALLEL SHARED(partNum, x_pos, y_pos, z_pos, resx, resy, resz, part_pos), &
        !$OMP          PRIVATE(ip, x_cell, y_cell, z_cell)
        !$OMP DO REDUCTION(+:part_pos)
        do ip=1,partNum 
            x_cell = int(x_pos(ip)/resx) + 1
            y_cell = int(y_pos(ip)/resy) + 1
            z_cell = int(z_pos(ip)/resz) + 1
            part_pos(ip,1) = x_cell
            part_pos(ip,2) = y_cell
            part_pos(ip,3) = z_cell
        enddo
        !$OMP END DO
        !$OMP END PARALLEL

        ! write(*,*) '--> Maximum number of cells covered by a particle (x,y,z)', int(2*maxval(hpart)/resx), &
        !                                                             int(2*maxval(hpart)/resy), &
        !                                                             int(2*maxval(hpart)/resz)
        ! write(*,*) '--> From particle field to density field using h and SPH kernel'

        !en este bucle recorremos las partículas y dependiendo de su "h" se recorren más o menos celdas alrededor
        !de la celda en la cual vive esa partícula, dando su contribución a la malla en el campo "field".
        
        !$OMP PARALLEL SHARED(partNum, part_pos, hpart, field), &
        !$OMP          PRIVATE(ip, thread_id, x_cell, y_cell, z_cell, contribution, until_cell_x, &
        !$OMP                  until_cell_y, until_cell_z, dx2, dy2, dz2, ix, iy, iz, ix2, iy2, iz2, W, r, norm)
        !$OMP DO REDUCTION(+:grid_field)
        do ip=1,partNum
            thread_id = OMP_get_thread_num()
            x_cell = part_pos(ip, 1)
            y_cell = part_pos(ip, 2)
            z_cell = part_pos(ip, 3)
            until_cell_x = int(2*hpart(ip)/resx) + 1
            until_cell_y = int(2*hpart(ip)/resy) + 1
            until_cell_z = int(2*hpart(ip)/resz) + 1
            allocate(contribution(2*until_cell_x+1, 2*until_cell_y+1, 2*until_cell_z+1)) !this is weird, should allocate contribution outside DO with same dimension as grid_field
            norm = 0.
            do iz=z_cell-until_cell_z, z_cell+until_cell_z, 1
             if (iz>0 .and. iz < nz+1) then
              do iy=y_cell-until_cell_y, y_cell+until_cell_y, 1
               if (iy>0 .and. iy < ny+1) then
                do ix=x_cell-until_cell_x, x_cell+until_cell_x, 1
                 if (ix>0 .and. ix < nx+1) then
                    dx2 = ((ix*resx-resx/2) - x_pos(ip))**2
                    dy2 = ((iy*resy-resy/2) - y_pos(ip))**2
                    dz2 = ((iz*resz-resz/2) - z_pos(ip))**2
                    r = sqrt(dx2 + dy2 + dz2)
                    W = 0.
                    call sph_kernel(r, hpart(ip), W)
                    norm = norm + W
                    !From whole grid to little grid over which the particle contributes
                    ix2 = ix-x_cell+until_cell_x + 1
                    iy2 = iy-y_cell+until_cell_y + 1
                    iz2 = iz-z_cell+until_cell_z + 1
                    contribution(ix2, iy2, iz2) = field(ip)*W
                 endif
                enddo
               endif
              enddo
             endif
            enddo
            
            !normalization
            if (norm > 0.) then
                contribution(:,:,:) = contribution(:,:,:)/norm
                !contributing
                do iz=z_cell-until_cell_z, z_cell+until_cell_z, 1
                 if (iz>0 .and. iz < nz+1) then
                  do iy=y_cell-until_cell_y, y_cell+until_cell_y, 1
                   if (iy>0 .and. iy < ny+1) then
                    do ix=x_cell-until_cell_x, x_cell+until_cell_x, 1
                     if (ix>0 .and. ix < nx+1) then
                        ix2 = ix-x_cell+until_cell_x + 1
                        iy2 = iy-y_cell+until_cell_y + 1
                        iz2 = iz-z_cell+until_cell_z + 1
                        grid_field(ix, iy, iz) = grid_field(ix, iy, iz) + contribution(ix2, iy2, iz2)
                     endif
                    enddo
                   endif
                  enddo
                 endif
                enddo
            else if (norm == 0.) then !if h is such small that don't contributes -> Nearest Grid Point
                grid_field(x_cell, y_cell, z_cell) = grid_field(x_cell, y_cell, z_cell) + field(ip)
            endif

            deallocate(contribution)
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
        ! write(*,*) '... done'
        deallocate(part_pos)
    end subroutine
    
    subroutine main(ncores, x_pos, y_pos, z_pos, Lx, Ly, Lz, field, kneigh, nx, ny, nz, partNum, grid_field, hpart)
        use omp_lib
        implicit none
        integer :: ncores, partNum, nx, ny, nz, kneigh
        real, dimension(partNum) :: x_pos, y_pos, z_pos, field, hpart
        real :: Lx, Ly, Lz
        real, dimension(nx,ny,nz) :: grid_field
        ! integer :: i,j,k

        !f2py intent(in) ncores, x_pos, y_pos, z_pos, Lx, Ly, Lz, field, kneigh, nx, ny, nz, partNum
        !f2py intent(out) grid_field, hpart
        !f2py depend(nx, ny, nz) grid_field
        !f2py depend(partNum) x_pos, y_pos, z_pos, field, hpart

        !control
        if (partNum<kneigh) then
            write(*,*) 'STOP: decrease kneigh !!!!!!!!!! ----> partNum, kneigh', partNum, kneigh
            stop
        endif
        
        !SET THE NUMBER OF CORES
        call OMP_SET_NUM_THREADS(ncores)

        !Start:
        !allocate h lenghts and initialize
        ! write(*,*) '                   '
        ! write(*,*) '*******************'
        ! write(*,*) '******** SPH ******'
        ! write(*,*) '*******************' 
        ! write(*,*) 'CHECK: partNum', partNum
        !Calculate h for each particle
        call sph_lenght(partNum, x_pos(:), y_pos(:), z_pos(:), Lx, Ly, Lz, kneigh, hpart(:))
        !From particle field to density field
        call sph_density(partnum, x_pos(:), y_pos(:), z_pos(:), field(:), nx, ny, nz, Lx, Ly, Lz, hpart(:), grid_field(:,:,:))
        !Finished, now deallocate h lenghts
        !End and save
        ! open (31,FILE='grid_field', STATUS='unknown',ACTION='write', FORM='unformatted')
        ! write(31) (((grid_field(i,j,k),i=1,nx),j=1,ny),k=1,nz)
        ! close(31)
        ! write(*,*) '****** END SPH ****'
        ! write(*,*) '                   '
    end subroutine

end module


!A LITTLE PROGRAM TO TEST THE MODULE

! program try
!     use SPH
!     implicit none
!     !INPUT:
!     integer :: partNum, nx, ny, nz, kneigh
!     real, dimension(:), allocatable :: x_pos, y_pos, z_pos, field
!     real :: Lx, Ly, Lz
!     integer :: ip
!     integer :: i,j,k
!     !OUPUT:
!     real, dimension(:), allocatable :: hpart
!     real, dimension(:,:,:), allocatable :: grid_field
!     partNum = int(1e6)
!     nx = 50; ny = 50 ; nz = 50
!     kneigh = 32
!     Lx = 40.; Ly = 40.; Lz = 40.

!     allocate(x_pos(partNum))
!     allocate(y_pos(partNum))
!     allocate(z_pos(partNum))
!     allocate(field(partNum))
!     allocate(hpart(partNum))
!     allocate(grid_field(nx, ny, nz))
!     field(:) = 1.
!     grid_field(:,:,:) = 0.
!     hpart(:) = 0.

!     do ip=1,partNum
!         call random_number(x_pos(ip))
!         call random_number(y_pos(ip))
!         call random_number(z_pos(ip))
!         x_pos(ip) = x_pos(ip)*Lx
!         y_pos(ip) = y_pos(ip)*Ly
!         z_pos(ip) = z_pos(ip)*Lz
!     enddo

!     call main(x_pos, y_pos, z_pos, Lx, Ly, Lz, field, kneigh, nx, ny, nz, partNum, grid_field)

!     open (31,FILE='grid_field', STATUS='unknown',ACTION='write', FORM='unformatted')
!     write(31) (((grid_field(i,j,k),i=1,nx),j=1,ny),k=1,nz)
!     close(31)

!     deallocate(x_pos)
!     deallocate(y_pos)
!     deallocate(z_pos)
!     deallocate(field)
!     deallocate(hpart)
!     deallocate(grid_field)

! end program