!
! © 2025. Triad National Security, LLC. All rights reserved.
!
! This program was produced under U.S. Government contract 89233218CNA000001
! for Los Alamos National Laboratory (LANL), which is operated by
! Triad National Security, LLC for the U.S. Department of Energy/National Nuclear
! Security Administration. All rights in the program are reserved by
! Triad National Security, LLC, and the U.S. Department of Energy/National
! Nuclear Security Administration. The Government is granted for itself and
! others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide
! license in this material to reproduce, prepare derivative works,
! distribute copies to the public, perform publicly and display publicly,
! and to permit others to do so.
!
! Author:
!    Kai Gao, kaigao@lanl.gov
!


module acoustic_iso_3d_boundary_saving

    use libflit
    use acoustic_iso_3d_vars

    implicit none

    integer :: xbwbeg1
    integer :: xbwbeg2
    integer :: xbwend1
    integer :: xbwend2
    integer :: xbwbegy
    integer :: xbwendy
    integer :: xbwbegz
    integer :: xbwendz

    integer :: ybwbeg1
    integer :: ybwbeg2
    integer :: ybwend1
    integer :: ybwend2
    integer :: ybwbegx
    integer :: ybwendx
    integer :: ybwbegz
    integer :: ybwendz

    integer :: zbwbegx
    integer :: zbwendx
    integer :: zbwbegy
    integer :: zbwendy
    integer :: zbwbeg1
    integer :: zbwbeg2
    integer :: zbwend1
    integer :: zbwend2

    integer :: bwiounit
    integer :: bwrecl

contains

    !
    !> Save final step wavefields for elastic media
    !
    subroutine output_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='stream', status='replace', action='write')

        write (funit) vx, vy, vz, p

        close (funit)

    end subroutine output_final_step_wavefield

    !
    !> Input final step wavefields for elastic media
    !
    !
    subroutine input_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='stream', status='old', action='read')

        read (funit) vx, vy, vz, p

        close (funit, status='delete')

    end subroutine input_final_step_wavefield

    !
    !> Prepare boundary saving
    !
    subroutine prepare_boundary_saving

        ! x
        xbwbeg1 = max(nx1, 1 - fdhalf + 1)
        xbwbeg2 = min(nx2, 1 - 1 + 1)
        if (xbwbeg2 < xbwbeg1) then
            xbwbeg1 = nx1 - 1
            xbwbeg2 = nx1 - 2
        end if

        xbwend1 = max(nx1, nx + 1)
        xbwend2 = min(nx2, nx + fdhalf)
        if (xbwend2 < xbwend1) then
            xbwend1 = nx2 - 1
            xbwend2 = nx2 - 2
        end if

        xbwbegy = max(ny1, 1 - fdhalf)
        xbwendy = min(ny2, ny + fdhalf)
        xbwbegz = max(nz1, 1 - fdhalf)
        xbwendz = min(nz2, nz + fdhalf)

        ! y
        ybwbegx = max(nx1, 1 - fdhalf)
        ybwendx = min(nx2, nx + fdhalf)

        ybwbeg1 = max(ny1, 1 - fdhalf + 1)
        ybwbeg2 = min(ny2, 1 - 1 + 1)
        if (ybwbeg2 < ybwbeg1) then
            ybwbeg1 = ny1 - 1
            ybwbeg2 = ny1 - 2
        end if

        ybwend1 = max(ny1, ny + 1)
        ybwend2 = min(ny2, ny + fdhalf)
        if (ybwend2 < ybwend1) then
            ybwend1 = ny2 - 1
            ybwend2 = ny2 - 2
        end if

        ybwbegz = max(nz1, 1 - fdhalf)
        ybwendz = min(nz2, nz + fdhalf)

        ! z
        zbwbegx = max(nx1, 1 - fdhalf)
        zbwendx = min(nx2, nx + fdhalf)
        zbwbegy = max(ny1, 1 - fdhalf)
        zbwendy = min(ny2, ny + fdhalf)

        zbwbeg1 = max(nz1, 1 - fdhalf + 1)
        zbwbeg2 = min(nz2, 1 - 1 + 1)
        if (zbwbeg2 < zbwbeg1) then
            zbwbeg1 = nz1 - 1
            zbwbeg2 = nz1 - 2
        end if

        zbwend1 = max(nz1, nz + 1)
        zbwend2 = min(nz2, nz + fdhalf)
        if (zbwend2 < zbwend1) then
            zbwend1 = nz2 - 1
            zbwend2 = nz2 - 2
        end if

        bwrecl = 0
        ! x
        bwrecl = bwrecl + (xbwbeg2 - xbwbeg1 + 1)*(xbwendy - xbwbegy + 1)*(xbwendz - xbwbegz + 1)
        bwrecl = bwrecl + (xbwend2 - xbwend1 + 1)*(xbwendy - xbwbegy + 1)*(xbwendz - xbwbegz + 1)
        ! y
        bwrecl = bwrecl + (ybwendx - ybwbegx + 1)*(ybwbeg2 - ybwbeg1 + 1)*(ybwendz - ybwbegz + 1)
        bwrecl = bwrecl + (ybwendx - ybwbegx + 1)*(ybwend2 - ybwend1 + 1)*(ybwendz - ybwbegz + 1)
        ! z
        bwrecl = bwrecl + (zbwendx - zbwbegx + 1)*(zbwendy - zbwbegy + 1)*(zbwbeg2 - zbwbeg1 + 1)
        bwrecl = bwrecl + (zbwendx - zbwbegx + 1)*(zbwendy - zbwbegy + 1)*(zbwend2 - zbwend1 + 1)

    end subroutine prepare_boundary_saving

    !
    !> Save boundary wavefield
    !
    subroutine save_boundary_wavefield(t)

        integer, intent(in) :: t

        if (bwrecl == 0) then
            return
        end if

        write (bwiounit, rec=t) &
            ! left boundary
        vx(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

    end subroutine save_boundary_wavefield

    !
    !> Inject boundary wavefield as boundary condition:
    !>        wavefield = R, not wavefield=wavefield + R
    !
    subroutine inject_boundary_wavefield(t)

        integer, intent(in) :: t

        if (bwrecl == 0) then
            return
        end if

        read (bwiounit, rec=t) &
            ! left boundary
        vx(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vx(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

    end subroutine inject_boundary_wavefield

    !
    !> Open boundary-saving file
    !
    subroutine open_boundary_saving

        if (bwrecl == 0) then
            return
        end if

        open (newunit=bwiounit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_boundary_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='direct', recl=4*bwrecl)

    end subroutine open_boundary_saving

    !
    !> Close boundary-saving file
    !
    subroutine close_boundary_saving(delete)

        logical, intent(in), optional :: delete

        if (bwrecl == 0) then
            return
        end if

        if (present(delete)) then
            if (delete) then
                close (bwiounit, status='delete')
            else
                close (bwiounit)
            end if
        else
            close (bwiounit)
        end if

    end subroutine close_boundary_saving

end module

! The boundary saving can be also implemented in the following point-by-point way,
! which is slightly memory efficient but requires intermediate arrays

!
!    integer, allocatable, dimension(:) :: ii, jj, kk
!    real, allocatable, dimension(:) :: bwx, bwy, bwz
!
!integer :: i, j, k
!
!ii = zeros(1)
!jj = zeros(1)
!kk = zeros(1)
!
!bwrecl = 0
!do k = nz1, nz2
!    do j = ny1, ny2
!        do i = nx1, nx2
!
!            if ((i >= 1 - fdhalf + 1 .and. i <= 1 - 1 + 1 &
    !                    .and. j >= 1 - fdhalf .and. j <= ny + fdhalf &
    !                    .and. k >= 1 - fdhalf .and. k <= nz + fdhalf) &
    !                    .or. &
    !                    (i >= nx + 1 .and. i <= nx + fdhalf &
    !                    .and. j >= 1 - fdhalf .and. j <= ny + fdhalf &
    !                    .and. k >= 1 - fdhalf .and. k <= nz + fdhalf) &
    !                    .or. &
    !                    (i >= 1 - fdhalf .and. i <= nx + fdhalf &
    !                    .and. j >= 1 - fdhalf + 1 .and. j <= 1 - 1 + 1 &
    !                    .and. k >= 1 - fdhalf .and. k <= nz + fdhalf) &
    !                    .or. &
    !                    (i >= 1 - fdhalf .and. i <= nx + fdhalf &
    !                    .and. j >= ny + 1 .and. j <= ny + fdhalf &
    !                    .and. k >= 1 - fdhalf .and. k <= nz + fdhalf) &
    !                    .or. &
    !                    (i >= 1 - fdhalf .and. i <= nx + fdhalf &
    !                    .and. j >= 1 - fdhalf .and. j <= ny + fdhalf &
    !                    .and. k >= 1 - fdhalf + 1 .and. k <= 1 - 1 + 1) &
    !                    .or. &
    !                    (i >= 1 - fdhalf .and. i <= nx + fdhalf &
    !                    .and. j >= 1 - fdhalf .and. j <= ny + fdhalf &
    !                    .and. k >= nz + 1 .and. k <= nz + fdhalf)) then
!
!                if (bwrecl == 0) then
!                    ii = [i]
!                    jj = [j]
!                    kk = [k]
!                else
!                    ii = [ii, i]
!                    jj = [jj, j]
!                    kk = [kk, k]
!                end if
!
!                bwrecl = bwrecl + 1
!
!            end if
!        end do
!    end do
!end do
!
!bwx = zeros(bwrecl)
!bwy = zeros(bwrecl)
!bwz = zeros(bwrecl)
!
!
!!$omp parallel do private(i)
!do i = 1, size(bwx)
!    bwx(i) = vx(ii(i), jj(i), kk(i))
!    bwy(i) = vy(ii(i), jj(i), kk(i))
!    bwz(i) = vz(ii(i), jj(i), kk(i))
!end do
!!$omp end parallel do
!
!write (bwiounit, rec=t) bwx, bwy, bwz
!
!
!read (bwiounit, rec=t) bwx, bwy, bwz
!
!!$omp parallel do private(i)
!do i = 1, size(bwx)
!    vx(ii(i), jj(i), kk(i)) = bwx(i)
!    vy(ii(i), jj(i), kk(i)) = bwy(i)
!    vz(ii(i), jj(i), kk(i)) = bwz(i)
!end do
!!$omp end parallel do
