!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HOW TO USE THE LINKED CELL LIST:                                  !
! === == === === ====== ==== ====                                   !
! 1) call the routine "count_cells", with the proper arguments.     !
!    This routine should only be called once.                       !
!                                                                   !
! 2) call the routine "neigh_list", to generate linked lists:       !
!    part_in_cell,  r_nei, l_nei.                                   !
!    This routine should only be called once.                       !
!                                                                   !
! 3) call get_cell_neigh_list to generate lists of all neighboring  !
!    cells for each cell.                                           !
!    This routine should only be called once.                       !
!                                                                   !
! 4) to loop particles in a cell:                                   !
!   i_part = part_in_cell(i_cell)                                   !
!   do while(i_part .ne. 0)                                         !
!       !Do whatever with r0(:,i_part)                              !
!       i_part = r_nei(i_part) !gets next particle in this cell     !
!   end do                                                          !
!                                                                   !
! 5) To loop between possible interacting particles of              !
!       a particle (i_part) :                                       !
!   a) call get_cell, to obtain the cell where i_part               !
!        is located (i_cell)                                        !
!   b) loop between neighboring cells of i_cell neig_cells(i_cell,:)!
!   c) loop over all particles in those cells (see 4), above)       !
!   do j=1,3**n_dim-1                                               !
!       j_cell = neig_cells(i_cell,j)                               !
!       j_part = part_in_cell(j_cell)                               !
!       do while(j_part .ne. 0)                                     ! 
!           !Do whatever with r0(:,j_part) and r0(:,i_part)         !  
!           j_part = r_nei(j_part) !gets next particle in this cell !    
!       end do                                                      !
!   end do                                                          !
!   d) loop over particles in the same cell also (i_cell)           !
!   j_part = part_in_cell(i_cell)                                   !
!   do while(j_part .ne. 0)                                         !
!       if(j_part .eq. i_part) then                                 !
!           j_part = r_nei(i_part)                                  !
!           cycle    ! don't count auto-interactions                !
!       end if                                                      !
!       !Do whatever with r0(:,i_part) and r0(:,j_part)             !
!       j_part = r_nei(j_part) !gets next particle in this cell     !
!   end do                                                          !
!                                                                   !
! 6) Check if particle i_part leaves a cell and update list:        !
!    a) Before updating position, get it's cell (old_cell)          !
!    b) call subroutine  "update_part_cell"                         !
!    This will automatically check if particle changes cell,        !
!    and will update linked lists  and cell list                    !
!                                                                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine neigh_list(r0, n_part, n_dim, n_cells_tot, boundary, l_cell, part_in_cell, r_nei, l_nei)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Routine to make cell lists in a linked way                        !
! part_in_cell(i_cell) gives the index of a                         !
! particle in the cell.                                             !
! r_nei(i_part) gives the next particle in the same cell.           !
! If there are no more particles in this cell,                      !
! then r_nei(i_part) = 0                                            !
! l_nei(j_part) gives the previous particle in the cell. If this is !
! the first particle in the cell, then l_nei(j_part) = 0            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
real(kind=8), intent(in) :: l_cell, r0(n_dim,n_part), boundary(n_dim)
integer, intent(in) :: n_part, n_dim, n_cells_tot
integer, intent(out) :: part_in_cell(n_cells_tot), r_nei(n_part), l_nei(n_part)
integer :: i_part, j_part, i_cell, i_dim, n_cells(n_dim)

!!!!!!!!!!! DO LINKED LIST FROM SCRATCH !!!!!!!!!!!!!!!!
part_in_cell(:) = 0
l_nei(:) = 0
r_nei(:) = 0
do i_part = 1, n_part

    !First get the cell index of the particle. Cell index is between 1 and
    !n_cells_tot
    call get_cell(r0(:,i_part), n_dim, n_cells, l_cell, i_cell)
    
    ! If Cell is empty, link cell with particle i_part
    if( part_in_cell(i_cell) .eq. 0 ) then
        part_in_cell(i_cell) = i_part
    else !If cell has already a particle in 
        j_part = part_in_cell(i_cell)
        part_in_cell(i_cell) = i_part
        r_nei(i_part) = j_part
        l_nei(i_part) = 0 !Done already
        l_nei(j_part) = i_part
    end if
end do

end subroutine

subroutine count_cells(n_dim, boundary, l_cell, n_cells, n_cells_tot)
implicit none
!This routine does the binning, and output is number of cells in each direction
!and total number of cells
integer, intent(in) :: n_dim ! dimension of space 
real(kind=8), intent(in) :: boundary(n_dim), l_cell ! length of simulation box
                                        !in each direction, linear length of a
                                        !cell
integer, intent(out) :: n_cells(n_dim), n_cells_tot !number of cells in each
                                        !direction, total number fo cells
!Cells are cubes of length l_cell
! l_cell should be equal to r_cut
!Set number of cells in each direction 
n_cells(:) = int( boundary(:) / l_cell ) + 1

!Set  total number of  cells in the system
n_cells_tot = 1
do i_dim = 1, n_dim
    n_cells_tot = n_cells_tot * n_cells(i_dim)
end do

end subroutine

subroutine get_cell(r0_part, n_dim, n_cells, l_cell, i_cell)
implicit none
real(kind=8), intent(in) :: r0_part(n_dim), l_cell
integer, intent(in) :: n_dim, n_cells(n_dim)
integer, intent(out) :: i_cell
integer :: pf, i_dim
pf = 1
i_cell = 1 
    do i_dim = 1, n_dim 
        i_cell = i_cell + pf * int( r0_part(i_dim) / l_cell )
        pf = pf * n_cells(i_dim) 
    end do
end subroutine

subroutine update_part_cell(r0_part, i_part, n_dim, n_cells, n_cells_tot, l_cell, old_cell, n_part, part_in_cell, r_nei, l_nei)
!This routine updates the linked lists of the cells, if a particle leaves a cell

!! IN VARIABLES
!r0_part  is the coordinates of the position of a particle
! i_part  is the index of the particle
! n_dim   is the number of space dimnesions of the system
! n_cells is an array which contains the number of cells in each of the n_dim directions
! n_cells_tot is the total number of cells in the system
! l_cell  is the length of the cell box
! old_cell is the cell where the particle (i_part) was before updating positions

!! IN OUT VARIABLES
! part_in_cell is an array that stores a particle index for each cell or 0 if the cell is empty
! r_nei is the right cell neighbour of each particle
! l_nei is the left cell neighbour of each particle

implicit none
real(kind=8), intent(in) :: r0_part(n_dim), l_cell
integer, intent(in) :: n_dim, old_cell, n_cells(n_dim), i_part, n_cells_tot, n_part
integer, intent(inout) :: part_in_cell(n_cells_tot), r_nei(n_part), l_nei(n_part)
integer ::  new_cell, j_part

!get new cell
call get_cell(r0_part, n_dim, n_cells, l_cell, new_cell)

if( new_cell .ne. old_cell ) then !Update cell list

    !Remove particle from old cell, and relink old cell.
    if( part_in_cell(old_cell) .eq. i_part ) then ! iif l_nei(i_part) = 0
        j_part = r_nei(i_part)
        part_in_cell(old_cell) = j_part
        if( j_part .ne. 0 ) l_nei(j_part) = 0
    else !l_nei(i_part) > 0
        j_part = r_nei(i_part)
        r_nei( l_nei(i_part) ) = j_part
        if( j_part .ne. 0 ) l_nei(j_part) = l_nei(i_part)
    end if
   
    !Add particle to new cell 
    if( part_in_cell(new_cell) .eq. 0 ) then !If there are no particles in this cell
        part_in_cell(new_cell) = i_part !Put particle in new cell
        l_nei(i_part) = 0
        r_nei(i_part) = 0
    else !If cell has already a particle in 
        j_part = part_in_cell(new_cell)
        part_in_cell(new_cell) = i_part
        r_nei(i_part) = j_part
        l_nei(i_part) = 0 !Done already
        l_nei(j_part) = i_part
    end if
end if

end subroutine

subroutine get_cell_neigh_list(n_cells_tot, n_cells, n_dim, cell_neigh_ls)
!Generate a list of neighboring cells for each cell
implicit none
integer, intent(in) :: n_dim, n_cells_tot, n_cells(n_dim)
integer, intent(out) :: cell_neigh_ls(n_cells_tot,3**n_dim-1)
integer :: i_cell

do i_cell = 1, n_cells_tot
    call make_cell_list(i_cell, n_cells, n_dim, cell_neigh_ls(i_cell,:))
end do

end subroutine get_cell_neigh_list

subroutine make_cell_list(in_cell, n_cells, n_dim, cell_list)
!This routine gives an array with the index of all the neighbouring cells of the
!input cell. Already corrected fo boundary conditions, assuming PBC in all 
!directions
!
! POSSIBLE OPTIMIZATION : Loop between neighboring cells in one direction only,
! to count interactions only once. ( do i=2,3; do j=2,3; bla; end do; end do)
! instead of 3**2-1 neighbors, (3**2 -1 )/2 neighbors.
implicit none
integer, intent(in) :: in_cell, n_dim, n_cells(n_dim)
integer, intent(out) :: cell_list( 3**n_dim - 1 )
integer :: i_dim, j_dim, i_cell, j_cell, i_list, x_cell, y_cell, z_cell, i, j, k, l

select case(n_dim)

case(1)
    i_list = 1
    do i = 1, 3
        j_cell = in_cell + i - 2
        if ( j_cell .eq. in_cell ) cycle ! Dont count in_cell as neighbour
    
        !HERE CHECK PBC
        if ( j_cell .eq. 0 ) j_cell = n_cells(1)
        if ( j_cell .eq. ( n_cells(1) + 1 ) ) j_cell = 1

        cell_list(i_list) = j_cell
        i_list = i_list + 1
    end do
    
     
case(2)
    y_cell = (in_cell - 1) / n_cells(1) + 1  ! Get y coordinate of in_cell
    x_cell = in_cell - n_cells(1) * (y_cell - 1)
    i_list = 1
    do i = 1, 3
        do j = 1, 3
            j_cell = in_cell + i-2 + (j-2) * n_cells(1)
            if ( j_cell .eq. in_cell ) cycle ! don't count same cell

            !HERE CHECK FOR PBC IN ALL DIRECTIONS 
            if ( (x_cell .eq. n_cells(1)) .and. (i.eq.3) ) j_cell = j_cell - n_cells(1)
            if ( (x_cell .eq. 1) .and. (i.eq.1) ) j_cell = j_cell + n_cells(1)
            if ( (y_cell .eq. n_cells(2) ) .and. (j.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2)    
            if ( (y_cell .eq. 1) .and. (j.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)           
            
            cell_list(i_list) = j_cell
            i_list = i_list + 1
        end do
    end do
    
case(3)
    z_cell = (in_cell - 1) / (n_cells(1)*n_cells(2)) + 1
    y_cell = ( in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - 1 ) / n_cells(1) + 1 
    x_cell = in_cell - n_cells(1) * n_cells(2) * (z_cell - 1) - n_cells(1) * (y_cell - 1) 
    i_list = 1
    do i = 1, 3
        do j = 1, 3 
            do k = 1, 3
                j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) * n_cells(1)
                if ( j_cell .eq. in_cell ) cycle

                !HERE CHECK PBC
                if ( (x_cell .eq. 1) .and. (i.eq.1) ) j_cell = j_cell + n_cells(1)
                if ( (x_cell .eq. n_cells(1)) .and. (i.eq.3) ) j_cell = j_cell - n_cells(1)
                if ( (y_cell .eq. 1) .and. (j.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)
                if ( (y_cell .eq. n_cells(2) ) .and. (j.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2) 
                if ( (z_cell .eq. 1) .and. (k.eq.1) ) j_cell = j_cell + n_cells(1)*n_cells(2)*n_cells(3)
                if ( (z_cell .eq. n_cells(3) ) .and. (k.eq.3) ) j_cell = j_cell - n_cells(1)*n_cells(2)*n_cells(3)   

                cell_list(i_list) = j_cell
                i_list = i_list + 1
            end do
        end do
    end do
    
case(4)
    i_list = 1
    do i = 1, 3
        do j = 1, 3 
            do k = 1, 3
                do l = 1, 3
                    j_cell = in_cell + i-2 + (j-2) * n_cells(1) + (k-2) * n_cells(2) + (l-2) * n_cells(3)
                    if ( j_cell .eq. in_cell ) cycle
                    !HERE CHECK PBC
                    cell_list(i_list) = j_cell
                    i_list = i_list + 1
                end do
            end do
        end do
    end do

case default
    print*, "ERROR in routine make_cell_list. "
    print*, "This routine works fine for 3 dimensions or less "
    call exit()

end select 

end subroutine
