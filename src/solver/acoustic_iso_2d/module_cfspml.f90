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


module acoustic_iso_2d_cfspml

    use libflit
    use acoustic_iso_2d_vars

    implicit none

    integer, parameter :: npower_k = 2
    integer, parameter :: npower_d = 2
    integer, parameter :: npower_a = 1
    real :: R0, kmax, alphamax

    real, allocatable, dimension(:) :: axi, azi, bxi, bzi, kxi, kzi
    real, allocatable, dimension(:) :: axh, azh, bxh, bzh, kxh, kzh

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

    end subroutine

    !
    !> Compute damping coefficients for isotropic media
    !
    subroutine compute_cfspml_damping_coef

        integer :: i, j
        real :: xdisti, xdisth, zdisti, zdisth
        real :: pmlvp

        ! Compute CSF-PML coefficients
        alphamax = 1.0*maxval(sgmtr%srcr(:)%f0)*const_pi
        R0 = 1.0e-5
        kmax = 1.0

        ! allocate memory for coefficient arrays
        call alloc_array(axi, [1, nx], pad=pml)
        call alloc_array(axh, [1, nx], pad=pml)
        call alloc_array(azi, [1, nz], pad=pml)
        call alloc_array(azh, [1, nz], pad=pml)

        call alloc_array(bxi, [1, nx], pad=pml)
        call alloc_array(bxh, [1, nx], pad=pml)
        call alloc_array(bzi, [1, nz], pad=pml)
        call alloc_array(bzh, [1, nz], pad=pml)

        call alloc_array(kxi, [1, nx], pad=pml)
        call alloc_array(kxh, [1, nx], pad=pml)
        call alloc_array(kzi, [1, nz], pad=pml)
        call alloc_array(kzh, [1, nz], pad=pml)

        kxi = 1.0
        kxh = 1.0
        kzi = 1.0
        kzh = 1.0

        pmlvp = max(maxval(vp(1, :)), maxval(vp(nx, :)), maxval(vp(:, 1)), maxval(vp(:, nz)))

        !$omp parallel do private(i, xdisti, xdisth)
        do i = -pml + 1, nx + pml

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

        !$omp parallel do private(j, zdisti, zdisth)
        do j = -pml + 1, nz + pml

            if (j <= 1) then
                zdisti = abs((j - 1)*dz)
                zdisth = abs((j - 1 - 0.5)*dz)
            else if (j >= nz + 1) then
                zdisti = abs((j - nz)*dz)
                zdisth = abs((j - nz - 0.5)*dz)
            end if

            if (j <= 1 .or. j >= nz + 1) then
                call damp_coef(pmlvp, dz, zdisti, azi(j), bzi(j), kzi(j))
                call damp_coef(pmlvp, dz, zdisth, azh(j), bzh(j), kzh(j))
            end if

        end do
        !$omp end parallel do

    end subroutine

end module
