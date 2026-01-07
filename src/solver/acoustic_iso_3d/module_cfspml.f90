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


module acoustic_iso_3d_cfspml

    use libflit
    use acoustic_iso_3d_vars

    implicit none

    integer, parameter :: npower_k = 2
    integer, parameter :: npower_d = 2
    integer, parameter :: npower_a = 1
    real :: R0, kmax, alphamax

    real, allocatable, dimension(:) :: axi, ayi, azi, bxi, byi, bzi, kxi, kyi, kzi
    real, allocatable, dimension(:) :: axh, ayh, azh, bxh, byh, bzh, kxh, kyh, kzh

contains

    subroutine damp_coef(bdvp, d, dist, a, b, k)

        real, intent(in) :: bdvp, d, dist
        real, intent(inout) :: a, b, k

        real :: damp, alpha
        real :: nd

        nd = dist/(pml*d)

        ! Diminishing alpha (Roden and Gedney, 2000); residual alpha to ensure long-time stability
        alpha = max(alphamax*(1.0 - nd**npower_a), 0.25*alphamax)
        damp = bdvp*(npower_d + 1.0)*log(1.0/R0)/(2.0*pml*d)*nd**npower_d
        k = 1.0 + (kmax - 1.0)*nd**npower_k

        a = -2.0*abs(dt)*damp/k/(2.0 + abs(dt)*(alpha + damp/k))
        b = (2.0 - abs(dt)*(alpha + damp/k))/(2.0 + abs(dt)*(alpha + damp/k))

    end subroutine damp_coef

    !
    !> Compute damping coefficients for isotropic media
    !
    subroutine compute_cfspml_damping_coef

        integer :: i, j, k
        real :: xdisti, xdisth, ydisti, ydisth, zdisti, zdisth

        ! Compute CFS-PML coefficients
        alphamax = 1.0*maxval(sgmtr%srcr(:)%f0)*const_pi
        R0 = 1.0e-5
        kmax = 1.0

        call alloc_array(axi, [nx1, nx2])
        call alloc_array(axh, [nx1, nx2])
        call alloc_array(ayi, [ny1, ny2])
        call alloc_array(ayh, [ny1, ny2])
        call alloc_array(azi, [nz1, nz2])
        call alloc_array(azh, [nz1, nz2])

        call alloc_array(bxi, [nx1, nx2])
        call alloc_array(bxh, [nx1, nx2])
        call alloc_array(byi, [ny1, ny2])
        call alloc_array(byh, [ny1, ny2])
        call alloc_array(bzi, [nz1, nz2])
        call alloc_array(bzh, [nz1, nz2])

        call alloc_array(kxi, [nx1, nx2])
        call alloc_array(kxh, [nx1, nx2])
        call alloc_array(kyi, [ny1, ny2])
        call alloc_array(kyh, [ny1, ny2])
        call alloc_array(kzi, [nz1, nz2])
        call alloc_array(kzh, [nz1, nz2])

        kxi = 1.0
        kxh = 1.0
        kyi = 1.0
        kyh = 1.0
        kzi = 1.0
        kzh = 1.0

        !$omp parallel do private(i, xdisti, xdisth)
        do i = nx1, nx2

            if (i <= 1) then
                xdisti = abs((i - 1)*dx)
                xdisth = abs((i - 1 - 0.5)*dx)
            else if (i >= nx + 1) then
                xdisti = abs((i - nx)*dx)
                xdisth = abs((i - nx - 0.5)*dx)
            end if

            if (i <= 1 .or. i >= nx + 1) then
                call damp_coef(pmlvp, dx, xdisti, axi(i), bxi(i), kxi(i))
                call damp_coef(pmlvp, dx, xdisth, axh(i), bxh(i), kxh(i))
            end if

        end do
        !$omp end parallel do

        !$omp parallel do private(j, ydisti, ydisth)
        do j = ny1, ny2

            if (j <= 1) then
                ydisti = abs((j - 1)*dy)
                ydisth = abs((j - 1 - 0.5)*dy)
            else if (j >= ny + 1) then
                ydisti = abs((j - ny)*dy)
                ydisth = abs((j - ny - 0.5)*dy)
            end if

            ! compute only for boundaries
            if (j <= 1 .or. j >= ny + 1) then
                call damp_coef(pmlvp, dy, ydisti, ayi(j), byi(j), kyi(j))
                call damp_coef(pmlvp, dy, ydisth, ayh(j), byh(j), kyh(j))
            end if

        end do
        !$omp end parallel do

        !$omp parallel do private(k, zdisti, zdisth)
        do k = nz1, nz2

            if (k <= 1) then
                zdisti = abs((k - 1)*dz)
                zdisth = abs((k - 1 - 0.5)*dz)
            else if (k >= nz + 1) then
                zdisti = abs((k - nz)*dz)
                zdisth = abs((k - nz - 0.5)*dz)
            end if

            ! compute only for boundaries
            if (k <= 1 .or. k >= nz + 1) then
                call damp_coef(pmlvp, dz, zdisti, azi(k), bzi(k), kzi(k))
                call damp_coef(pmlvp, dz, zdisth, azh(k), bzh(k), kzh(k))
            end if

        end do
        !$omp end parallel do

    end subroutine

end module
