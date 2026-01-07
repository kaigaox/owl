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

module elastic_tti_2d_wavefield

    use libflit
    use elastic_tti_2d_vars
    use elastic_tti_2d_cfspml

    implicit none

    ! Finite-difference stencils
#include 'macro_fd_stencil.f90'

contains

    !
    !> Update wavefields
    !
    subroutine update_wavefield(dt, &
            stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
            stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
            vx_ixhz, vz_ixhz, &
            vx_hxiz, vz_hxiz, &
            memory_pdxvx_hxhz, &
            memory_pdxvz_hxhz, &
            memory_pdxvx_ixiz, &
            memory_pdxvz_ixiz, &
            memory_pdxxx_hxiz, &
            memory_pdxxz_hxiz, &
            memory_pdxxx_ixhz, &
            memory_pdxxz_ixhz, &
            memory_pdzvx_hxhz, &
            memory_pdzvz_hxhz, &
            memory_pdzvx_ixiz, &
            memory_pdzvz_ixiz, &
            memory_pdzxz_hxiz, &
            memory_pdzzz_hxiz, &
            memory_pdzxz_ixhz, &
            memory_pdzzz_ixhz)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :), intent(inout) :: &
            stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
            stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
            vx_ixhz, vz_ixhz, &
            vx_hxiz, vz_hxiz
        real, allocatable, dimension(:, :), intent(inout) :: &
            memory_pdxvx_hxhz, &
            memory_pdxvz_hxhz, &
            memory_pdxvx_ixiz, &
            memory_pdxvz_ixiz, &
            memory_pdxxx_hxiz, &
            memory_pdxxz_hxiz, &
            memory_pdxxx_ixhz, &
            memory_pdxxz_ixhz, &
            memory_pdzvx_hxhz, &
            memory_pdzvz_hxhz, &
            memory_pdzvx_ixiz, &
            memory_pdzvz_ixiz, &
            memory_pdzxz_hxiz, &
            memory_pdzzz_hxiz, &
            memory_pdzxz_ixhz, &
            memory_pdzzz_ixhz

        integer :: i, j
        real :: pdxvx, pdzvx, pdxvz, pdzvz
        real :: pdxxx, pdzxz, pdxxz, pdzzz
        real :: rho_hx, rho_hz
        real :: c11_ixiz, c13_ixiz, c15_ixiz, c33_ixiz, c35_ixiz, c55_ixiz
        real :: c11_hxhz, c13_hxhz, c15_hxhz, c33_hxhz, c35_hxhz, c55_hxhz

        ! Stress in the center of a cell (i, j)
        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz, &
            !$omp c11_ixiz, c13_ixiz, c15_ixiz, c33_ixiz, c35_ixiz, c55_ixiz) collapse(2) schedule(auto)
        do j = nz1, nz2
            do i = nx1, nx2

                c11_ixiz = c11(i, j)
                c13_ixiz = c13(i, j)
                c15_ixiz = c15(i, j)
                c33_ixiz = c33(i, j)
                c35_ixiz = c35(i, j)
                c55_ixiz = c55(i, j)

                pdxvx = pdxw_stencil(vx_hxiz, i, j)/dx
                pdxvz = pdxw_stencil(vz_hxiz, i, j)/dx
                pdzvx = pdzw_stencil(vx_ixhz, i, j)/dz
                pdzvz = pdzw_stencil(vz_ixhz, i, j)/dz

                if (i <= 0 .or. i >= nx + 1 .or. j <= 0 .or. j >= nz + 1) then
                    memory_pdxvx_ixiz(i, j) = axii(i, j)*pdxvx + bxii(i, j)*memory_pdxvx_ixiz(i, j)
                    pdxvx = (pdxvx + memory_pdxvx_ixiz(i, j))/kxii(i, j)
                    memory_pdxvz_ixiz(i, j) = axii(i, j)*pdxvz + bxii(i, j)*memory_pdxvz_ixiz(i, j)
                    pdxvz = (pdxvz + memory_pdxvz_ixiz(i, j))/kxii(i, j)
                    memory_pdzvx_ixiz(i, j) = azii(i, j)*pdzvx + bzii(i, j)*memory_pdzvx_ixiz(i, j)
                    pdzvx = (pdzvx + memory_pdzvx_ixiz(i, j))/kzii(i, j)
                    memory_pdzvz_ixiz(i, j) = azii(i, j)*pdzvz + bzii(i, j)*memory_pdzvz_ixiz(i, j)
                    pdzvz = (pdzvz + memory_pdzvz_ixiz(i, j))/kzii(i, j)
                end if

                stressxx_ixiz(i, j) = stressxx_ixiz(i, j) + dt*( &
                    c11_ixiz*pdxvx &
                    + c13_ixiz*pdzvz &
                    + c15_ixiz*(pdxvz + pdzvx))

                stresszz_ixiz(i, j) = stresszz_ixiz(i, j) + dt*( &
                    c13_ixiz*pdxvx &
                    + c33_ixiz*pdzvz &
                    + c35_ixiz*(pdxvz + pdzvx))

                stressxz_ixiz(i, j) = stressxz_ixiz(i, j) + dt*( &
                    c15_ixiz*pdxvx &
                    + c35_ixiz*pdzvz &
                    + c55_ixiz*(pdxvz + pdzvx))

            end do
        end do
        !$omp end parallel do

        ! Stress on the corners of a cell (i+-1/2, j+-1/2)
        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz, &
            !$omp c11_hxhz, c13_hxhz, c15_hxhz, c33_hxhz, c35_hxhz, c55_hxhz) collapse(2) schedule(auto)
        do j = nz1 - 1, nz2 - 1
            do i = nx1 - 1, nx2 - 1

                pdxvx = pdxw_stencil(vx_ixhz, i, j + 1)/dx
                pdxvz = pdxw_stencil(vz_ixhz, i, j + 1)/dx
                pdzvx = pdzw_stencil(vx_hxiz, i + 1, j)/dz
                pdzvz = pdzw_stencil(vz_hxiz, i + 1, j)/dz

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j + 1 <= 1 .or. j + 1 >= nz + 1) then
                    memory_pdxvx_hxhz(i + 1, j + 1) = axhh(i + 1, j + 1)*pdxvx + bxhh(i + 1, j + 1)*memory_pdxvx_hxhz(i + 1, j + 1)
                    pdxvx = (pdxvx + memory_pdxvx_hxhz(i + 1, j + 1))/kxhh(i + 1, j + 1)
                    memory_pdxvz_hxhz(i + 1, j + 1) = axhh(i + 1, j + 1)*pdxvz + bxhh(i + 1, j + 1)*memory_pdxvz_hxhz(i + 1, j + 1)
                    pdxvz = (pdxvz + memory_pdxvz_hxhz(i + 1, j + 1))/kxhh(i + 1, j + 1)
                    memory_pdzvx_hxhz(i + 1, j + 1) = azhh(i + 1, j + 1)*pdzvx + bzhh(i + 1, j + 1)*memory_pdzvx_hxhz(i + 1, j + 1)
                    pdzvx = (pdzvx + memory_pdzvx_hxhz(i + 1, j + 1))/kzhh(i + 1, j + 1)
                    memory_pdzvz_hxhz(i + 1, j + 1) = azhh(i + 1, j + 1)*pdzvz + bzhh(i + 1, j + 1)*memory_pdzvz_hxhz(i + 1, j + 1)
                    pdzvz = (pdzvz + memory_pdzvz_hxhz(i + 1, j + 1))/kzhh(i + 1, j + 1)
                end if

                ! ! Harmonic average may result in instability
                ! c11_hxhz = 4.0/sum(1.0/c11(i:i + 1, j:j + 1))
                ! c13_hxhz = 4.0/sum(1.0/c13(i:i + 1, j:j + 1))
                ! c15_hxhz = 4.0/sum(1.0/c15(i:i + 1, j:j + 1))
                ! c33_hxhz = 4.0/sum(1.0/c33(i:i + 1, j:j + 1))
                ! c35_hxhz = 4.0/sum(1.0/c35(i:i + 1, j:j + 1))
                ! c55_hxhz = 4.0/sum(1.0/c55(i:i + 1, j:j + 1))

                c11_hxhz = 0.25*sum(c11(i:i + 1, j:j + 1))
                c13_hxhz = 0.25*sum(c13(i:i + 1, j:j + 1))
                c15_hxhz = 0.25*sum(c15(i:i + 1, j:j + 1))
                c33_hxhz = 0.25*sum(c33(i:i + 1, j:j + 1))
                c35_hxhz = 0.25*sum(c35(i:i + 1, j:j + 1))
                c55_hxhz = 0.25*sum(c55(i:i + 1, j:j + 1))

                stressxx_hxhz(i + 1, j + 1) = stressxx_hxhz(i + 1, j + 1) + dt*( &
                    c11_hxhz*pdxvx &
                    + c13_hxhz*pdzvz &
                    + c15_hxhz*(pdxvz + pdzvx))

                stresszz_hxhz(i + 1, j + 1) = stresszz_hxhz(i + 1, j + 1) + dt*( &
                    c13_hxhz*pdxvx &
                    + c33_hxhz*pdzvz &
                    + c35_hxhz*(pdxvz + pdzvx))

                stressxz_hxhz(i + 1, j + 1) = stressxz_hxhz(i + 1, j + 1) + dt*( &
                    c15_hxhz*pdxvx &
                    + c35_hxhz*pdzvz &
                    + c55_hxhz*(pdxvz + pdzvx))

            end do
        end do
        !$omp end parallel do

        ! Particle velocity on the vertical edges of a cell (i+-1/2, j)
        !$omp parallel do private(i, j, pdxxx, pdzxz, pdxxz, pdzzz, rho_hx) collapse(2) schedule(auto)
        do j = nz1, nz2
            do i = nx1 - 1, nx2 - 1

                pdxxx = pdxw_stencil(stressxx_ixiz, i, j)/dx
                pdxxz = pdxw_stencil(stressxz_ixiz, i, j)/dx
                pdzxz = pdzw_stencil(stressxz_hxhz, i + 1, j)/dz
                pdzzz = pdzw_stencil(stresszz_hxhz, i + 1, j)/dz

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j <= 0 .or. j >= nz + 1) then
                    memory_pdxxx_hxiz(i + 1, j) = axhi(i + 1, j)*pdxxx + bxhi(i + 1, j)*memory_pdxxx_hxiz(i + 1, j)
                    pdxxx = (pdxxx + memory_pdxxx_hxiz(i + 1, j))/kxhi(i + 1, j)
                    memory_pdxxz_hxiz(i + 1, j) = axhi(i + 1, j)*pdxxz + bxhi(i + 1, j)*memory_pdxxz_hxiz(i + 1, j)
                    pdxxz = (pdxxz + memory_pdxxz_hxiz(i + 1, j))/kxhi(i + 1, j)
                    memory_pdzxz_hxiz(i + 1, j) = azhi(i + 1, j)*pdzxz + bzhi(i + 1, j)*memory_pdzxz_hxiz(i + 1, j)
                    pdzxz = (pdzxz + memory_pdzxz_hxiz(i + 1, j))/kzhi(i + 1, j)
                    memory_pdzzz_hxiz(i + 1, j) = azhi(i + 1, j)*pdzzz + bzhi(i + 1, j)*memory_pdzzz_hxiz(i + 1, j)
                    pdzzz = (pdzzz + memory_pdzzz_hxiz(i + 1, j))/kzhi(i + 1, j)
                end if

                rho_hx = 0.5*(rho(i, j) + rho(i + 1, j))

                vx_hxiz(i + 1, j) = vx_hxiz(i + 1, j) + dt/rho_hx*(pdxxx + pdzxz)
                vz_hxiz(i + 1, j) = vz_hxiz(i + 1, j) + dt/rho_hx*(pdxxz + pdzzz)

            end do
        end do
        !$omp end parallel do

        ! Particle velocity on the horizontal edges of a cell (i, j+-1/2)
        !$omp parallel do private(i, j, pdxxx, pdzxz, pdxxz, pdzzz, rho_hz) collapse(2) schedule(auto)
        do j = nz1 - 1, nz2 - 1
            do i = nx1, nx2

                pdxxx = pdxw_stencil(stressxx_hxhz, i, j + 1)/dx
                pdxxz = pdxw_stencil(stressxz_hxhz, i, j + 1)/dx
                pdzxz = pdzw_stencil(stressxz_ixiz, i, j)/dz
                pdzzz = pdzw_stencil(stresszz_ixiz, i, j)/dz

                if (i <= 0 .or. i >= nx + 1 .or. j + 1 <= 1 .or. j + 1 >= nz + 1) then
                    memory_pdxxx_ixhz(i, j + 1) = axih(i, j + 1)*pdxxx + bxih(i, j + 1)*memory_pdxxx_ixhz(i, j + 1)
                    pdxxx = (pdxxx + memory_pdxxx_ixhz(i, j + 1))/kxih(i, j + 1)
                    memory_pdxxz_ixhz(i, j + 1) = axih(i, j + 1)*pdxxz + bxih(i, j + 1)*memory_pdxxz_ixhz(i, j + 1)
                    pdxxz = (pdxxz + memory_pdxxz_ixhz(i, j + 1))/kxih(i, j + 1)
                    memory_pdzxz_ixhz(i, j + 1) = azih(i, j + 1)*pdzxz + bzih(i, j + 1)*memory_pdzxz_ixhz(i, j + 1)
                    pdzxz = (pdzxz + memory_pdzxz_ixhz(i, j + 1))/kzih(i, j + 1)
                    memory_pdzzz_ixhz(i, j + 1) = azih(i, j + 1)*pdzzz + bzih(i, j + 1)*memory_pdzzz_ixhz(i, j + 1)
                    pdzzz = (pdzzz + memory_pdzzz_ixhz(i, j + 1))/kzih(i, j + 1)
                end if

                rho_hz = 0.5*(rho(i, j) + rho(i, j + 1))

                vx_ixhz(i, j + 1) = vx_ixhz(i, j + 1) + dt/rho_hz*(pdxxx + pdzxz)
                vz_ixhz(i, j + 1) = vz_ixhz(i, j + 1) + dt/rho_hz*(pdxxz + pdzzz)

            end do
        end do
        !$omp end parallel do

    end subroutine update_wavefield

    !
    !> Update wavefields in model with topographic free surface
    !
    subroutine update_wavefield_free_surface(dt, &
            stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
            stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
            vx_ixhz, vz_ixhz, &
            vx_hxiz, vz_hxiz, &
            memory_pdxvx_hxhz, &
            memory_pdxvz_hxhz, &
            memory_pdxvx_ixiz, &
            memory_pdxvz_ixiz, &
            memory_pdxxx_hxiz, &
            memory_pdxxz_hxiz, &
            memory_pdxxx_ixhz, &
            memory_pdxxz_ixhz, &
            memory_pdzvx_hxhz, &
            memory_pdzvz_hxhz, &
            memory_pdzvx_ixiz, &
            memory_pdzvz_ixiz, &
            memory_pdzxx_hxiz, &
            memory_pdzxz_hxiz, &
            memory_pdzzz_hxiz, &
            memory_pdzxx_ixhz, &
            memory_pdzxz_ixhz, &
            memory_pdzzz_ixhz)

        real, intent(in) :: dt
        real, allocatable, dimension(:, :), intent(inout) :: &
            stressxx_ixiz, stresszz_ixiz, stressxz_ixiz, &
            stressxx_hxhz, stresszz_hxhz, stressxz_hxhz, &
            vx_ixhz, vz_ixhz, &
            vx_hxiz, vz_hxiz
        real, allocatable, dimension(:, :), intent(inout) :: &
            memory_pdxvx_hxhz, &
            memory_pdxvz_hxhz, &
            memory_pdxvx_ixiz, &
            memory_pdxvz_ixiz, &
            memory_pdxxx_hxiz, &
            memory_pdxxz_hxiz, &
            memory_pdxxx_ixhz, &
            memory_pdxxz_ixhz, &
            memory_pdzvx_hxhz, &
            memory_pdzvz_hxhz, &
            memory_pdzvx_ixiz, &
            memory_pdzvz_ixiz, &
            memory_pdzxx_hxiz, &
            memory_pdzxz_hxiz, &
            memory_pdzzz_hxiz, &
            memory_pdzxx_ixhz, &
            memory_pdzxz_ixhz, &
            memory_pdzzz_ixhz

        integer :: i, j
        real :: pdxvx, pdzvx, pdxvz, pdzvz
        real :: pdxxx, pdzxx, pdzxz, pdxxz, pdzzz
        real :: pd1, pd2, pd3
        real :: tpaz, tpcz
        real :: d11, d12, d21, d22, dd, f1, f2
        real :: hx, eta
        real :: rho_hx, rho_hz
        real :: c11_ixiz, c13_ixiz, c15_ixiz, c33_ixiz, c35_ixiz, c55_ixiz
        real :: c11_hxhz, c13_hxhz, c15_hxhz, c33_hxhz, c35_hxhz, c55_hxhz

        ! Particle velocity above the free surface (j = 1) are zero
        !$omp parallel do private(j) schedule(auto)
        do j = 1, fdhalf
            vx_hxiz(:, 1 - j) = 0.0
            vz_hxiz(:, 1 - j) = 0.0
            vx_ixhz(:, 2 - j) = 0.0
            vz_ixhz(:, 2 - j) = 0.0
        end do
        !$omp end parallel do

        ! Stress in the center of a cell (i, j)
        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz, pd1, pd2, pd3, tpaz, tpcz, &
            !$omp d11, d12, d21, d22, dd, f1, f2, hx, eta, &
            !$omp c11_ixiz, c13_ixiz, c15_ixiz, c33_ixiz, c35_ixiz, c55_ixiz) collapse(2) schedule(auto)
        do j = 1, nz + pml
            do i = -pml + 1, nx + pml

                eta = eta_zz_i(j)

                tpaz = max(0.0, eta_max - eta)/(topo_i(i) + depth_max)*slopex_i(i)
                tpcz = eta_max/(topo_i(i) + depth_max)

                c11_ixiz = c11(i, j)
                c13_ixiz = c13(i, j)
                c15_ixiz = c15(i, j)
                c33_ixiz = c33(i, j)
                c35_ixiz = c35(i, j)
                c55_ixiz = c55(i, j)

                ! x derivatives
                pdxvx = pdxw_stencil(vx_hxiz, i, j)/dx
                pdxvz = pdxw_stencil(vz_hxiz, i, j)/dx

                if (j == 1) then
                    ! At the free surface j = 1, i.e., z = (j-1)*dz = 0

                    ! In the PML region
                    if (i <= 0 .or. i >= nx + 1) then
                        memory_pdxvx_ixiz(i, j) = axii(i, j)*pdxvx + bxii(i, j)*memory_pdxvx_ixiz(i, j)
                        pdxvx = (pdxvx + memory_pdxvx_ixiz(i, j))/kxii(i, j)
                        memory_pdxvz_ixiz(i, j) = axii(i, j)*pdxvz + bxii(i, j)*memory_pdxvz_ixiz(i, j)
                        pdxvz = (pdxvz + memory_pdxvz_ixiz(i, j))/kxii(i, j)
                    end if

                    ! z derivatives are solved from x
                    hx = slopex_i(i)
                    d11 = tpcz*(hx**2*c11_ixiz + hx*2.0*c15_ixiz + c55_ixiz)
                    d12 = tpcz*(hx**2*c15_ixiz + hx*(c13_ixiz + c55_ixiz) + c35_ixiz)
                    d22 = tpcz*(hx**2*c55_ixiz + hx*2.0*c35_ixiz + c33_ixiz)
                    dd = d12**2 - d11*d22
                    f1 = (-hx*c11_ixiz - c15_ixiz)*pdxvx + (-hx*c15_ixiz - c55_ixiz)*pdxvz
                    f2 = (-hx*c15_ixiz - c13_ixiz)*pdxvx + (-hx*c55_ixiz - c35_ixiz)*pdxvz
                    pdzvx = (-d22*f1 + d12*f2)/dd
                    pdzvz = (-d11*f2 + d12*f1)/dd

                else
                    ! Below the free surface, conventional derivatives

                    pdzvx = pdzw_stencil(vx_ixhz, i, j)/dz*eta_dz_scaling_i(j)
                    pdzvz = pdzw_stencil(vz_ixhz, i, j)/dz*eta_dz_scaling_i(j)

                    if (i <= 0 .or. i >= nx + 1 .or. j >= nz + 1) then
                        memory_pdxvx_ixiz(i, j) = axii(i, j)*pdxvx + bxii(i, j)*memory_pdxvx_ixiz(i, j)
                        pdxvx = (pdxvx + memory_pdxvx_ixiz(i, j))/kxii(i, j)
                        memory_pdxvz_ixiz(i, j) = axii(i, j)*pdxvz + bxii(i, j)*memory_pdxvz_ixiz(i, j)
                        pdxvz = (pdxvz + memory_pdxvz_ixiz(i, j))/kxii(i, j)
                        memory_pdzvx_ixiz(i, j) = azii(i, j)*pdzvx + bzii(i, j)*memory_pdzvx_ixiz(i, j)
                        pdzvx = (pdzvx + memory_pdzvx_ixiz(i, j))/kzii(i, j)
                        memory_pdzvz_ixiz(i, j) = azii(i, j)*pdzvz + bzii(i, j)*memory_pdzvz_ixiz(i, j)
                        pdzvz = (pdzvz + memory_pdzvz_ixiz(i, j))/kzii(i, j)
                    end if

                end if

                pd1 = pdxvx + pdzvx*tpaz
                pd2 = pdzvz*tpcz
                pd3 = pdzvx*tpcz + pdxvz + pdzvz*tpaz

                stressxx_ixiz(i, j) = stressxx_ixiz(i, j) + dt*( &
                    c11_ixiz*pd1 &
                    + c13_ixiz*pd2 &
                    + c15_ixiz*pd3)

                stresszz_ixiz(i, j) = stresszz_ixiz(i, j) + dt*( &
                    c13_ixiz*pd1 &
                    + c33_ixiz*pd2 &
                    + c35_ixiz*pd3)

                stressxz_ixiz(i, j) = stressxz_ixiz(i, j) + dt*( &
                    c15_ixiz*pd1 &
                    + c35_ixiz*pd2 &
                    + c55_ixiz*pd3)

            end do
        end do
        !$omp end parallel do

        ! Stress on the corners of a cell (i+-1/2, j+-1/2)
        !$omp parallel do private(i, j, pdxvx, pdzvx, pdxvz, pdzvz, pd1, pd2, pd3, tpaz, tpcz, &
            !$omp d11, d12, d22, dd, f1, f2, hx, eta, &
            !$omp c11_hxhz, c13_hxhz, c15_hxhz, c33_hxhz, c35_hxhz, c55_hxhz) collapse(2) schedule(auto)
        do j = 1, nz + pml - 1
            do i = -pml, nx + pml - 1

                eta = eta_zz_h(j + 1)

                tpaz = max(0.0, eta_max - eta)/(topo_h(i + 1) + depth_max)*slopex_h(i + 1)
                tpcz = eta_max/(topo_h(i + 1) + depth_max)

                c11_hxhz = 0.25*sum(c11(i:i + 1, j:j + 1))
                c13_hxhz = 0.25*sum(c13(i:i + 1, j:j + 1))
                c15_hxhz = 0.25*sum(c15(i:i + 1, j:j + 1))
                c33_hxhz = 0.25*sum(c33(i:i + 1, j:j + 1))
                c35_hxhz = 0.25*sum(c35(i:i + 1, j:j + 1))
                c55_hxhz = 0.25*sum(c55(i:i + 1, j:j + 1))

                pdxvx = pdxw_stencil(vx_ixhz, i, j + 1)/dx
                pdxvz = pdxw_stencil(vz_ixhz, i, j + 1)/dx

                pdzvx = pdzw_stencil(vx_hxiz, i + 1, j)/dz*eta_dz_scaling_h(j + 1)
                pdzvz = pdzw_stencil(vz_hxiz, i + 1, j)/dz*eta_dz_scaling_h(j + 1)

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j + 1 >= nz + 1) then
                    memory_pdxvx_hxhz(i + 1, j + 1) = axhh(i + 1, j + 1)*pdxvx + bxhh(i + 1, j + 1)*memory_pdxvx_hxhz(i + 1, j + 1)
                    pdxvx = (pdxvx + memory_pdxvx_hxhz(i + 1, j + 1))/kxhh(i + 1, j + 1)
                    memory_pdxvz_hxhz(i + 1, j + 1) = axhh(i + 1, j + 1)*pdxvz + bxhh(i + 1, j + 1)*memory_pdxvz_hxhz(i + 1, j + 1)
                    pdxvz = (pdxvz + memory_pdxvz_hxhz(i + 1, j + 1))/kxhh(i + 1, j + 1)
                    memory_pdzvx_hxhz(i + 1, j + 1) = azhh(i + 1, j + 1)*pdzvx + bzhh(i + 1, j + 1)*memory_pdzvx_hxhz(i + 1, j + 1)
                    pdzvx = (pdzvx + memory_pdzvx_hxhz(i + 1, j + 1))/kzhh(i + 1, j + 1)
                    memory_pdzvz_hxhz(i + 1, j + 1) = azhh(i + 1, j + 1)*pdzvz + bzhh(i + 1, j + 1)*memory_pdzvz_hxhz(i + 1, j + 1)
                    pdzvz = (pdzvz + memory_pdzvz_hxhz(i + 1, j + 1))/kzhh(i + 1, j + 1)
                end if

                pd1 = pdxvx + pdzvx*tpaz
                pd2 = pdzvz*tpcz
                pd3 = pdzvx*tpcz + pdxvz + pdzvz*tpaz

                stressxx_hxhz(i + 1, j + 1) = stressxx_hxhz(i + 1, j + 1) + dt*( &
                    c11_hxhz*pd1 &
                    + c13_hxhz*pd2 &
                    + c15_hxhz*pd3)

                stresszz_hxhz(i + 1, j + 1) = stresszz_hxhz(i + 1, j + 1) + dt*( &
                    c13_hxhz*pd1 &
                    + c33_hxhz*pd2 &
                    + c35_hxhz*pd3)

                stressxz_hxhz(i + 1, j + 1) = stressxz_hxhz(i + 1, j + 1) + dt*( &
                    c15_hxhz*pd1 &
                    + c35_hxhz*pd2 &
                    + c55_hxhz*pd3)

            end do
        end do
        !$omp end parallel do

        ! Apply mirror boundary condition to mimic free surface
        !$omp parallel do private(j) schedule(auto)
        do j = 1, fdhalf

            stressxx_ixiz(:, 1 - j) = -stressxx_ixiz(:, 1 + j)
            stresszz_ixiz(:, 1 - j) = -stresszz_ixiz(:, 1 + j)
            stressxz_ixiz(:, 1 - j) = -stressxz_ixiz(:, 1 + j)

            ! When the free surface is topographic, the boundary condition is sigma \cdot n = 0,
            ! but zz and xz are not necessarily zero. Therefore,
            ! the following should not be imposed here unless it is horizontal free surface:
            ! stresszz_ixiz(i, 1) = 0.0
            ! stressxz_ixiz(i, 1) = 0.0

            stressxx_hxhz(:, 2 - j) = -stressxx_hxhz(:, 1 + j)
            stresszz_hxhz(:, 2 - j) = -stresszz_hxhz(:, 1 + j)
            stressxz_hxhz(:, 2 - j) = -stressxz_hxhz(:, 1 + j)

        end do
        !$omp end parallel do

        ! Particle velocity on the vertical edges of a cell (i+-1/2, j)
        !$omp parallel do private(i, j, pdxxx, pdzxx, pdzxz, pdxxz, pdzzz, pd1, pd2, tpaz, tpcz, eta, rho_hx) collapse(2) schedule(auto)
        do j = 1, nz + pml
            do i = -pml, nx + pml - 1

                eta = eta_zz_i(j)

                tpaz = max(0.0, eta_max - eta)/(topo_h(i + 1) + depth_max)*slopex_h(i + 1)
                tpcz = eta_max/(topo_h(i + 1) + depth_max)

                pdxxx = pdxw_stencil(stressxx_ixiz, i, j)/dx
                pdxxz = pdxw_stencil(stressxz_ixiz, i, j)/dx

                pdzxx = pdzw_stencil(stressxx_hxhz, i + 1, j)/dz*eta_dz_scaling_i(j)
                pdzxz = pdzw_stencil(stressxz_hxhz, i + 1, j)/dz*eta_dz_scaling_i(j)
                pdzzz = pdzw_stencil(stresszz_hxhz, i + 1, j)/dz*eta_dz_scaling_i(j)

                if (i + 1 <= 1 .or. i + 1 >= nx + 1 .or. j >= nz + 1) then
                    memory_pdxxx_hxiz(i + 1, j) = axhi(i + 1, j)*pdxxx + bxhi(i + 1, j)*memory_pdxxx_hxiz(i + 1, j)
                    pdxxx = (pdxxx + memory_pdxxx_hxiz(i + 1, j))/kxhi(i + 1, j)
                    memory_pdxxz_hxiz(i + 1, j) = axhi(i + 1, j)*pdxxz + bxhi(i + 1, j)*memory_pdxxz_hxiz(i + 1, j)
                    pdxxz = (pdxxz + memory_pdxxz_hxiz(i + 1, j))/kxhi(i + 1, j)
                    memory_pdzxx_hxiz(i + 1, j) = azhi(i + 1, j)*pdzxx + bzhi(i + 1, j)*memory_pdzxx_hxiz(i + 1, j)
                    pdzxx = (pdzxx + memory_pdzxx_hxiz(i + 1, j))/kzhi(i + 1, j)
                    memory_pdzxz_hxiz(i + 1, j) = azhi(i + 1, j)*pdzxz + bzhi(i + 1, j)*memory_pdzxz_hxiz(i + 1, j)
                    pdzxz = (pdzxz + memory_pdzxz_hxiz(i + 1, j))/kzhi(i + 1, j)
                    memory_pdzzz_hxiz(i + 1, j) = azhi(i + 1, j)*pdzzz + bzhi(i + 1, j)*memory_pdzzz_hxiz(i + 1, j)
                    pdzzz = (pdzzz + memory_pdzzz_hxiz(i + 1, j))/kzhi(i + 1, j)
                end if

                pd1 = pdxxx &
                    + tpaz*pdzxx &
                    + tpcz*pdzxz
                pd2 = pdxxz &
                    + tpaz*pdzxz &
                    + tpcz*pdzzz

                rho_hx = 0.5*(rho(i, j) + rho(i + 1, j))

                vx_hxiz(i + 1, j) = vx_hxiz(i + 1, j) + dt/rho_hx*pd1
                vz_hxiz(i + 1, j) = vz_hxiz(i + 1, j) + dt/rho_hx*pd2

            end do
        end do
        !$omp end parallel do

        ! Particle velocity on the horizontal edges of a cell (i, j+-1/2)
        !$omp parallel do private(i, j, pdxxx, pdzxx, pdzxz, pdxxz, pdzzz, pd1, pd2, tpaz, tpcz, eta, rho_hz) collapse(2) schedule(auto)
        do j = 1, nz + pml - 1
            do i = -pml + 1, nx + pml

                eta = eta_zz_h(j + 1)

                tpaz = max(0.0, eta_max - eta)/(topo_i(i) + depth_max)*slopex_i(i)
                tpcz = eta_max/(topo_i(i) + depth_max)

                pdxxx = pdxw_stencil(stressxx_hxhz, i, j + 1)/dx
                pdxxz = pdxw_stencil(stressxz_hxhz, i, j + 1)/dx

                pdzxx = pdzw_stencil(stressxx_ixiz, i, j)/dz*eta_dz_scaling_h(j + 1)
                pdzxz = pdzw_stencil(stressxz_ixiz, i, j)/dz*eta_dz_scaling_h(j + 1)
                pdzzz = pdzw_stencil(stresszz_ixiz, i, j)/dz*eta_dz_scaling_h(j + 1)

                if (i <= 0 .or. i >= nx + 1 .or. j + 1 >= nz + 1) then
                    memory_pdxxx_ixhz(i, j + 1) = axih(i, j + 1)*pdxxx + bxih(i, j + 1)*memory_pdxxx_ixhz(i, j + 1)
                    pdxxx = (pdxxx + memory_pdxxx_ixhz(i, j + 1))/kxih(i, j + 1)
                    memory_pdxxz_ixhz(i, j + 1) = axih(i, j + 1)*pdxxz + bxih(i, j + 1)*memory_pdxxz_ixhz(i, j + 1)
                    pdxxz = (pdxxz + memory_pdxxz_ixhz(i, j + 1))/kxih(i, j + 1)
                    memory_pdzxx_ixhz(i, j + 1) = azih(i, j + 1)*pdzxx + bzih(i, j + 1)*memory_pdzxx_ixhz(i, j + 1)
                    pdzxx = (pdzxx + memory_pdzxx_ixhz(i, j + 1))/kzih(i, j + 1)
                    memory_pdzxz_ixhz(i, j + 1) = azih(i, j + 1)*pdzxz + bzih(i, j + 1)*memory_pdzxz_ixhz(i, j + 1)
                    pdzxz = (pdzxz + memory_pdzxz_ixhz(i, j + 1))/kzih(i, j + 1)
                    memory_pdzzz_ixhz(i, j + 1) = azih(i, j + 1)*pdzzz + bzih(i, j + 1)*memory_pdzzz_ixhz(i, j + 1)
                    pdzzz = (pdzzz + memory_pdzzz_ixhz(i, j + 1))/kzih(i, j + 1)
                end if

                pd1 = pdxxx &
                    + tpaz*pdzxx &
                    + tpcz*pdzxz
                pd2 = pdxxz &
                    + tpaz*pdzxz &
                    + tpcz*pdzzz

                rho_hz = 0.5*(rho(i, j) + rho(i, j + 1))

                vx_ixhz(i, j + 1) = vx_ixhz(i, j + 1) + dt/rho_hz*pd1
                vz_ixhz(i, j + 1) = vz_ixhz(i, j + 1) + dt/rho_hz*pd2

            end do
        end do
        !$omp end parallel do

    end subroutine

    !
    !> Add source
    !> For modeling in free surface, only add source on points below the free surface.
    !
    subroutine add_source(t)

        integer, intent(in) :: t

        integer :: k
        integer :: nbeg, nend
        integer :: irx, irz, sgx, sgz, shx, shz
        real :: amp, polar
        real :: dz_s

        do k = 1, sgmtr%ns

            if (yn_free_surface) then
                dz_s = eta_dz_i(sgmtr%srcr(k)%gz)*(topo_max + depth_max)/eta_max
            else
                dz_s = dz
            end if

            nbeg = nint(sgmtr%srcr(k)%t0/dt) + 1
            nend = nbeg + sgmtr%srcr(k)%nt - 1

            sgx = sgmtr%srcr(k)%gx
            sgz = sgmtr%srcr(k)%gz
            shx = sgmtr%srcr(k)%hx
            shz = sgmtr%srcr(k)%hz

            if (t >= nbeg .and. t <= nend) then

                amp = sgmtr%srcr(k)%stf(t - nbeg + 1)*sgmtr%srcr(k)%amp*dt

                select case (sgmtr%srcr(k)%mechanism)

                    case ('force')
                        ! Force vector
                        polar = sgmtr%srcr(k)%polar

                        amp = amp/rho(sgx, sgz)

                        !$omp parallel do private(irx, irz) collapse(2) schedule(auto)
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                    vx_hxiz(shx + irx, sgz + irz) = &
                                        vx_hxiz(shx + irx, sgz + irz) + sin(polar)*amp &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                    vz_hxiz(shx + irx, sgz + irz) = &
                                        vz_hxiz(shx + irx, sgz + irz) + cos(polar)*amp &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                end if
                                if (ifelse(yn_free_surface, shz + irz >= 2, .true.)) then
                                    vx_ixhz(sgx + irx, shz + irz) = &
                                        vx_ixhz(sgx + irx, shz + irz) + sin(polar)*amp &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                    vz_ixhz(sgx + irx, shz + irz) = &
                                        vz_ixhz(sgx + irx, shz + irz) + cos(polar)*amp &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                end if
                            end do
                        end do
                        !$omp end parallel do

                    case ('explosion')
                        ! Explosive source

                        amp = amp/(dx*dz_s)

                        !$omp parallel do private(irx, irz) collapse(2) schedule(auto)
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                    stressxx_ixiz(sgx + irx, sgz + irz) = &
                                        stressxx_ixiz(sgx + irx, sgz + irz) - amp &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                    stresszz_ixiz(sgx + irx, sgz + irz) = &
                                        stresszz_ixiz(sgx + irx, sgz + irz) - amp &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                end if
                                if (ifelse(yn_free_surface, shz + irz >= 2, .true.)) then
                                    stressxx_hxhz(shx + irx, shz + irz) = &
                                        stressxx_hxhz(shx + irx, shz + irz) - amp &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                    stresszz_hxhz(shx + irx, shz + irz) = &
                                        stresszz_hxhz(shx + irx, shz + irz) - amp &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                end if
                            end do
                        end do
                        !$omp end parallel do

                    case ('mt')
                        ! Moment tensor

                        amp = amp/(dx*dz_s)

                        !$omp parallel do private(irx, irz) collapse(2) schedule(auto)
                        do irz = -nkw, nkw
                            do irx = -nkw, nkw
                                if (ifelse(yn_free_surface, sgz + irz >= 2, .true.)) then
                                    stressxx_ixiz(sgx + irx, sgz + irz) = &
                                        stressxx_ixiz(sgx + irx, sgz + irz) - amp*sgmtr%srcr(k)%moment_tensor(1, 1) &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                    stresszz_ixiz(sgx + irx, sgz + irz) = &
                                        stresszz_ixiz(sgx + irx, sgz + irz) - amp*sgmtr%srcr(k)%moment_tensor(3, 3) &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                    stressxz_ixiz(sgx + irx, sgz + irz) = &
                                        stressxz_ixiz(sgx + irx, sgz + irz) - amp*sgmtr%srcr(k)%moment_tensor(1, 3) &
                                        *sgmtr%srcr(k)%interp_ix(irx) &
                                        *sgmtr%srcr(k)%interp_iz(irz)
                                end if
                                if (ifelse(yn_free_surface, shz + irz >= 2, .true.)) then
                                    stressxx_hxhz(shx + irx, shz + irz) = &
                                        stressxx_hxhz(shx + irx, shz + irz) - amp*sgmtr%srcr(k)%moment_tensor(1, 1) &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                    stresszz_hxhz(shx + irx, shz + irz) = &
                                        stresszz_hxhz(shx + irx, shz + irz) - amp*sgmtr%srcr(k)%moment_tensor(3, 3) &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                    stressxz_hxhz(shx + irx, shz + irz) = &
                                        stressxz_hxhz(shx + irx, shz + irz) - amp*sgmtr%srcr(k)%moment_tensor(1, 3) &
                                        *sgmtr%srcr(k)%interp_hx(irx) &
                                        *sgmtr%srcr(k)%interp_hz(irz)
                                end if
                            end do
                        end do
                        !$omp end parallel do

                end select
            end if
        end do

    end subroutine add_source

end module
