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


#ifdef _dim3_

!
!> Compute Hilbert transform of a wavefield along some axis
!
function compute_hilbert_transform(w, dim) result(v)

    real, dimension(:, :, :) :: w
    integer :: dim
    real, allocatable, dimension(:, :, :) :: v

    integer :: i, j, k

    v = w

    select case (dim)

        case (1)
            !$omp parallel do private(i, j, k) collapse(2) schedule(auto)
            do k = nz1_interior, nz2_interior
                do j = ny1_interior, ny2_interior
                    call hilbert_transform(v(:, j, k))
                end do
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i, j, k) collapse(2) schedule(auto)
            do k = nz1_interior, nz2_interior
                do i = nx1_interior, nx2_interior
                    call hilbert_transform(v(i, :, k))
                end do
            end do
            !$omp end parallel do

        case (3)
            !$omp parallel do private(i, j, k) collapse(2) schedule(auto)
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior
                    call hilbert_transform(v(i, j, :))
                end do
            end do
            !$omp end parallel do

    end select

end function

!
!> Compute FWI gradient based on directional cross-correlation of wavefields
!
function compute_directional_gradient(ws, wr, wsh, wrh, sgnh, dim) result(v)

    real, dimension(:, :, :) :: ws, wr, wsh, wrh
    integer :: sgnh
    integer :: dim
    real, allocatable, dimension(:, :, :) :: v

    integer :: i, j, k

    allocate(v(interior_region))

    select case (dim)

        case (1)
            !$omp parallel do private(i, j, k) collapse(2) schedule(auto)
            do k = nz1_interior, nz2_interior
                do j = ny1_interior, ny2_interior
                    v(nx1_interior:nx2_interior, j, k) = &
                        (ws(nx1_interior:nx2_interior, j, k)*wr(nx1_interior:nx2_interior, j, k) &
                        + sgnh*wsh(nx1_interior:nx2_interior, j, k)*wrh(nx1_interior:nx2_interior, j, k))
                end do
            end do
            !$omp end parallel do

        case (2)
            !$omp parallel do private(i, j, k) collapse(2) schedule(auto)
            do k = nz1_interior, nz2_interior
                do i = nx1_interior, nx2_interior
                    v(i, ny1_interior:ny2_interior, k) = &
                        (ws(i, ny1_interior:ny2_interior, k)*wr(i, ny1_interior:ny2_interior, k) &
                        + sgnh*wsh(i, ny1_interior:ny2_interior, k)*wrh(i, ny1_interior:ny2_interior, k))
                end do
            end do
            !$omp end parallel do

        case (3)
            !$omp parallel do private(i, j, k) collapse(2) schedule(auto)
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior
                    v(i, j, nz1_interior:nz2_interior) = &
                        (ws(i, j, nz1_interior:nz2_interior)*wr(i, j, nz1_interior:nz2_interior) &
                        + sgnh*wsh(i, j, nz1_interior:nz2_interior)*wrh(i, j, nz1_interior:nz2_interior))
                end do
            end do
            !$omp end parallel do

    end select

end function

#endif
