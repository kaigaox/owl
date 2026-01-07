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


module elastic_tti_2d_boundary_saving

    use libflit
    use elastic_tti_2d_vars

    implicit none

    integer :: xbwbeg1
    integer :: xbwbeg2
    integer :: xbwend1
    integer :: xbwend2
    integer :: xbwbegy
    integer :: xbwendy
    integer :: xbwbegz
    integer :: xbwendz

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
    !> Save final step wavefields
    !
    subroutine output_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin', &
            form='unformatted', access='stream', status='replace', action='write')

        write (funit) &
            vx_hxiz, vx_ixhz, &
            vz_hxiz, vz_ixhz, &
            stressxx_ixiz, stressxx_hxhz, &
            stresszz_ixiz, stresszz_hxhz, &
            stressxz_ixiz, stressxz_hxhz

        close (funit)

    end subroutine output_final_step_wavefield

    !
    !> Input final step wavefields for elastic media
    !
    subroutine input_final_step_wavefield

        integer :: funit

        open (newunit=funit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_final_step_wavefield.bin', &
            form='unformatted', access='stream', status='old', action='read')

        read (funit) &
            vx_hxiz, vx_ixhz, &
            vz_hxiz, vz_hxiz, &
            stressxx_ixiz, stressxx_hxhz, &
            stresszz_ixiz, stresszz_hxhz, &
            stressxz_ixiz, stressxz_hxhz

        close (funit, status='delete')

    end subroutine input_final_step_wavefield

    !
    !> Prepare boundary saving
    !
    subroutine prepare_boundary_saving

        ! x
        xbwbeg1 = 1 - fdhalf + 1
        xbwbeg2 = 1
        xbwend1 = nx + 1
        xbwend2 = nx + fdhalf
        xbwbegz = 1 - fdhalf
        xbwendz = nz + fdhalf

        ! z
        zbwbegx = 1 - fdhalf
        zbwendx = nx + fdhalf
        zbwbeg1 = 1 - fdhalf + 1
        zbwbeg2 = 1
        zbwend1 = nz + 1
        zbwend2 = nz + fdhalf

        bwrecl = 0
        ! x
        bwrecl = bwrecl + (xbwbeg2 - xbwbeg1 + 1)*(xbwendz - xbwbegz + 1)
        bwrecl = bwrecl + (xbwend2 - xbwend1 + 1)*(xbwendz - xbwbegz + 1)
        ! z
        bwrecl = bwrecl + (zbwendx - zbwbegx + 1)*(zbwbeg2 - zbwbeg1 + 1)
        bwrecl = bwrecl + (zbwendx - zbwbegx + 1)*(zbwend2 - zbwend1 + 1)

        bwrecl = bwrecl*2*2

    end subroutine prepare_boundary_saving

    !
    !> Save boundary wavefield
    !
    subroutine save_boundary_wavefield(t)

        integer, intent(in) :: t

        write (bwiounit, rec=t) &
            ! left boundary
        vx_hxiz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vx_hxiz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vx_hxiz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_hxiz(zbwbegx:zbwendx, zbwend1:zbwend2), &
            ! left boundary
        vz_hxiz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vz_hxiz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vz_hxiz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_hxiz(zbwbegx:zbwendx, zbwend1:zbwend2), &
            ! left boundary
        vx_ixhz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vx_ixhz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vx_ixhz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_ixhz(zbwbegx:zbwendx, zbwend1:zbwend2), &
            ! left boundary
        vz_ixhz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vz_ixhz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vz_ixhz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_ixhz(zbwbegx:zbwendx, zbwend1:zbwend2)

    end subroutine save_boundary_wavefield

    !
    !> Inject boundary wavefield as boundary condition:
    !>        wavefield=R, not wavefield=wavefield+R
    !
    subroutine inject_boundary_wavefield(t)

        integer, intent(in) :: t

        read (bwiounit, rec=t) &
            ! left boundary
        vx_hxiz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vx_hxiz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vx_hxiz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_hxiz(zbwbegx:zbwendx, zbwend1:zbwend2), &
            ! left boundary
        vz_hxiz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vz_hxiz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vz_hxiz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_hxiz(zbwbegx:zbwendx, zbwend1:zbwend2), &
            ! left boundary
        vx_ixhz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vx_ixhz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vx_ixhz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vx_ixhz(zbwbegx:zbwendx, zbwend1:zbwend2), &
            ! left boundary
        vz_ixhz(xbwbeg1:xbwbeg2, xbwbegz:xbwendz), &
            ! right boundary
        vz_ixhz(xbwend1:xbwend2, xbwbegz:xbwendz), &
            ! top boundary
        vz_ixhz(zbwbegx:zbwendx, zbwbeg1:zbwbeg2), &
            ! bottom boundary
        vz_ixhz(zbwbegx:zbwendx, zbwend1:zbwend2)

    end subroutine inject_boundary_wavefield

    !
    !> Open boundary saving file
    !
    subroutine open_boundary_saving

        open (newunit=bwiounit, file=tidy(dir_working)//'/shot_' &
            //num2str(sgmtr%id) &
            //'_boundary_wavefield.bin', &
            form='unformatted', access='direct', recl=4*bwrecl)

    end subroutine open_boundary_saving

    !
    !> Close boundary saving file
    !
    subroutine close_boundary_saving(delete)

        logical, intent(in), optional :: delete

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
