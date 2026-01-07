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


module elastic_vhtiort_3d_boundary_saving

    use libflit
    use elastic_vhtiort_3d_vars

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

        write (funit) vx, vy, vz, stressxx, stressyy, stresszz, stressyz, stressxz, stressxy

        close (funit)

    end subroutine output_final_step_wavefield

    !
    !> Input final step wavefields for elastic media
    !
    subroutine input_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin.'//num2str(rankid_group), &
            form='unformatted', access='stream', status='old', action='read')

        read (funit) vx, vy, vz, stressxx, stressyy, stresszz, stressyz, stressxz, stressxy

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

        ! For elastic wavefield reconstruction, must save all three particle velocity wavefields
        bwrecl = bwrecl*3

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
        vx(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2)

    end subroutine save_boundary_wavefield

    !
    !> Inject boundary wavefield as boundary condition:
    !>        wavefield=R, not wavefield=wavefield+R
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
        vx(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vx(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vy(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vy(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vy(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vy(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
            ! top boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vy(zbwbegx:zbwendx, zbwbegy:zbwendy, zbwend1:zbwend2), &
            ! left boundary
        vz(xbwbeg1:xbwbeg2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! right boundary
        vz(xbwend1:xbwend2, xbwbegy:xbwendy, xbwbegz:xbwendz), &
            ! front boundary
        vz(ybwbegx:ybwendx, ybwbeg1:ybwbeg2, ybwbegz:ybwendz), &
            ! back boundary
        vz(ybwbegx:ybwendx, ybwend1:ybwend2, ybwbegz:ybwendz), &
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
    !> Close boundary saving file
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

    end subroutine

end module
