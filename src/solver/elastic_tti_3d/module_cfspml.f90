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


module elastic_tti_3d_cfspml

    use libflit
    use elastic_tti_3d_vars
    use blas95
    use lapack95

    implicit none

    integer, parameter :: npower_k = 2
    integer, parameter :: npower_d = 2
    integer, parameter :: npower_a = 1
    real :: eigenderiv_max = -0.01
    real :: alphamax
    integer :: npower
    real :: R0
    real :: kmax

    real, allocatable, dimension(:, :, :) :: idax, iday, idaz, hdax, hday, hdaz
    real, allocatable, dimension(:, :, :) :: alphaxi, alphayi, alphazi, alphaxh, alphayh, alphazh
    real, allocatable, dimension(:, :, :) :: kappaxi, kappayi, kappazi, kappaxh, kappayh, kappazh

    real, allocatable, dimension(:, :) :: dampratio_left, dampratio_right
    real, allocatable, dimension(:, :) :: dampratio_front, dampratio_back
    real, allocatable, dimension(:, :) :: dampratio_top, dampratio_bottom

contains

    !
    !> Calculate the eigenvalue derivative for M-PML
    !
    function mpml_eigenvalue_derivative(c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, &
            direction, ratio, k1, k2, k3) result(deriv)

        real, intent(in) :: c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66
        character(len=*), intent(in) :: direction
        real, intent(in) :: ratio, k1, k2, k3
        real :: deriv

        real, dimension(1:27, 1:27) :: a, vl, vr
        real, dimension(1:27, 1:27) :: b
        real, dimension(1:27) :: wr, wi
        integer, dimension(1:3) :: psn
        integer :: i, l
        real, dimension(1:3) :: eigd
        real :: temp(1, 1)
        real :: evl(1:1, 1:27), evr(1:27, 1:1)
        real :: diag(1:27)
        integer :: ind

        ! Undamped system and its left/rigth eigenvectors
        ! For efficiency purpose, assign blocks instead one by one
        a = 0.0
        a(1:18, 19:27) = reshape([ &
            c11*k1, c16*k1, c15*k1, c11*k1, c16*k1, c15*k1, c11*k1, c16*k1, c15*k1, &
            c12*k1, c26*k1, c25*k1, c12*k1, c26*k1, c25*k1, c12*k1, c26*k1, c25*k1, &
            c13*k1, c36*k1, c35*k1, c13*k1, c36*k1, c35*k1, c13*k1, c36*k1, c35*k1, &
            c14*k1, c46*k1, c45*k1, c14*k1, c46*k1, c45*k1, c14*k1, c46*k1, c45*k1, &
            c15*k1, c56*k1, c55*k1, c15*k1, c56*k1, c55*k1, c15*k1, c56*k1, c55*k1, &
            c16*k1, c66*k1, c56*k1, c16*k1, c66*k1, c56*k1, c16*k1, c66*k1, c56*k1, &
            c16*k2, c12*k2, c14*k2, c16*k2, c12*k2, c14*k2, c16*k2, c12*k2, c14*k2, &
            c26*k2, c22*k2, c24*k2, c26*k2, c22*k2, c24*k2, c26*k2, c22*k2, c24*k2, &
            c36*k2, c23*k2, c34*k2, c36*k2, c23*k2, c34*k2, c36*k2, c23*k2, c34*k2, &
            c46*k2, c24*k2, c44*k2, c46*k2, c24*k2, c44*k2, c46*k2, c24*k2, c44*k2, &
            c56*k2, c25*k2, c45*k2, c56*k2, c25*k2, c45*k2, c56*k2, c25*k2, c45*k2, &
            c66*k2, c26*k2, c46*k2, c66*k2, c26*k2, c46*k2, c66*k2, c26*k2, c46*k2, &
            c15*k3, c14*k3, c13*k3, c15*k3, c14*k3, c13*k3, c15*k3, c14*k3, c13*k3, &
            c25*k3, c24*k3, c23*k3, c25*k3, c24*k3, c23*k3, c25*k3, c24*k3, c23*k3, &
            c35*k3, c34*k3, c33*k3, c35*k3, c34*k3, c33*k3, c35*k3, c34*k3, c33*k3, &
            c45*k3, c44*k3, c34*k3, c45*k3, c44*k3, c34*k3, c45*k3, c44*k3, c34*k3, &
            c55*k3, c45*k3, c35*k3, c55*k3, c45*k3, c35*k3, c55*k3, c45*k3, c35*k3, &
            c56*k3, c46*k3, c36*k3, c56*k3, c46*k3, c36*k3, c56*k3, c46*k3, c36*k3 &
            ], [18, 9], order=[2, 1])/1.0e9
        a(19:27, 1:18) = reshape([ &
            k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, &
            0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, 0.0, 0.0, 0.0, 0.0, k1, 0.0, &
            0.0, 0.0, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, 0.0, 0.0, k2, &
            0.0, k2, 0.0, 0.0, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, 0.0, &
            0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, 0.0, 0.0, k2, 0.0, 0.0, 0.0, 0.0, 0.0, k2, 0.0, 0.0, &
            0.0, 0.0, 0.0, 0.0, k3, 0.0, 0.0, 0.0, 0.0, 0.0, k3, 0.0, 0.0, 0.0, 0.0, 0.0, k3, 0.0, &
            0.0, 0.0, 0.0, k3, 0.0, 0.0, 0.0, 0.0, 0.0, k3, 0.0, 0.0, 0.0, 0.0, 0.0, k3, 0.0, 0.0, &
            0.0, 0.0, k3, 0.0, 0.0, 0.0, 0.0, 0.0, k3, 0.0, 0.0, 0.0, 0.0, 0.0, k3, 0.0, 0.0, 0.0 &
            ], [9, 18], order=[2, 1])

        ! Calcualte the eigevalue derivative
        select case (direction)
            case ('x')
                diag = -([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
                    ratio, ratio, ratio, ratio, ratio, ratio, &
                    ratio, ratio, ratio, ratio, ratio, ratio, &
                    1.0, 1.0, 1.0, &
                    ratio, ratio, ratio, &
                    ratio, ratio, ratio])

            case ('y')
                diag = -([ratio, ratio, ratio, ratio, ratio, ratio, &
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
                    ratio, ratio, ratio, ratio, ratio, ratio, &
                    ratio, ratio, ratio, &
                    1.0, 1.0, 1.0, &
                    ratio, ratio, ratio])

            case ('z')
                diag = -([ratio, ratio, ratio, ratio, ratio, ratio, &
                    ratio, ratio, ratio, ratio, ratio, ratio, &
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, &
                    ratio, ratio, ratio, &
                    ratio, ratio, ratio, &
                    1.0, 1.0, 1.0])
        end select
        b = 0.0
        do i = 1, 27
            b(i, i) = diag(i)
        end do

        ! Solve for left and right eigenvectors using LAPACK:
        ! The reason for computing the eigenvectors and eigenvalues of A rather
        ! than tilde(A) is that it is only necessary to compute the eigenvectors around d = 0
        ! When d = 0, tilde(A) = A
        wr = 0
        wi = 0
        vl = 0
        vr = 0
        call geev(a, wr, wi, vl, vr)
        where (abs(wr) < 1.0e-2*maxval(abs(wr)))
            wr = 0
        end where
        where (abs(wi) < 1.0e-2*maxval(abs(wr)))
            wi = 0
        end where

        ! Sort and use the smallest eigenvalues
        l = 1
        do i = 1, 27
            ! Here, matrix A_0 is stripped off imaginary unit i,
            ! therefore its eigenvalues should be strictly real numbers and 0s
            ! and here wr (real part of the eigevalue) are its eigenvalues
            if (l <= 3 .and. wr(i) < 0 .and. wi(i) == 0) then
                psn(l) = i
                l = l + 1
            end if
        end do
        psn = sort(psn)

        ! Compute the derivative of eigenvalues
        temp = 0
        do i = 1, 3
            ind = psn(i)
            if (ind >= 1 .and. ind <= 27) then
                evl(1, :) = vl(:, ind)
                evr(:, 1) = vr(:, ind)
                temp = matmul(evl, evr)
                if (temp(1, 1) == 0) then
                    temp = temp + 1.0e-15
                else
                    temp = temp + sign(1.0e-15, temp(1, 1))
                end if
                temp = matmul(matmul(evl, b), evr)/temp(1, 1)
                eigd(i) = temp(1, 1)
            end if
        end do

        deriv = maxval(pack(eigd, mask=(.not. isnan(eigd))))

    end function mpml_eigenvalue_derivative

    !
    !> Get the damping ratio for MPML
    !
    !  Refs: Method adopted from Meza-Fajardo and Papageorgiou (2008)
    !        with some extension to deal with general anisotropic media
    !
    function damp_ratio(c11, c12, c13, c14, c15, c16, &
            c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, direction) &
            result(ratio)

        real, intent(in) :: c11, c12, c13, c14, c15, c16
        real, intent(in) ::      c22, c23, c24, c25, c26
        real, intent(in) ::           c33, c34, c35, c36
        real, intent(in) ::                c44, c45, c46
        real, intent(in) ::                     c55, c56
        real, intent(in) ::                          c66
        character(len=*), intent(in) :: direction
        real :: ratio

        real :: kx, ky, kz, dtheta, dphi, derig, theta, phi
        integer :: i, j, ng1, ng2
        real :: alpha, dratio

        ! Devide pi, can be coarser
        ng1 = 60
        ng2 = 60
        dtheta = const_pi/ng1
        dphi = const_pi/ng2

        ! Critical value
        alpha = eigenderiv_max

        ! Update step size
        dratio = 2.0e-3

        ! Calculate the ratio
        ratio = dratio
        do i = 1, ng1
            do j = 1, ng2

                ! To avoid singularity, add a deviation to regular positions
                theta = (i + 0.5)*dtheta
                phi = (j + 0.5)*dphi

                ! Wavenumber vector
                kx = sin(theta)*cos(phi)
                ky = sin(theta)*sin(phi)
                kz = cos(theta)

                derig = mpml_eigenvalue_derivative(c11, c12, c13, c14, c15, c16, &
                    c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, &
                    direction, ratio, kx, ky, kz)

                do while (derig >= alpha .and. ratio < 1.0)

                    ! Increase ratio until derivative of eigenvalue smaller than alpha
                    ratio = ratio + dratio

                    ! Compute for the derivative of eigenvalues
                    derig = mpml_eigenvalue_derivative(c11, c12, c13, c14, c15, c16, &
                        c22, c23, c24, c25, c26, c33, c34, c35, c36, c44, c45, c46, c55, c56, c66, &
                        direction, ratio, kx, ky, kz)

                end do

            end do
        end do

    end function damp_ratio

    !
    !> Compute MPML damping profile ratio for a face
    !
    subroutine compute_damping_profile_ratio_face(ratio, axis, position)

        real, allocatable, dimension(:, :), intent(inout) :: ratio
        character(len=*), intent(in) :: axis, position

        real, allocatable, dimension(:, :) :: r
        integer :: n1, n2, m1, m2, i, j, ii, jj
        integer :: interval, bindex
        integer :: n1a, n1b, n2a, n2b

        bindex = 0

        select case (axis)
            case ('x')
                n1 = ny2 - ny1 + 1
                n2 = nz2 - nz1 + 1
                n1a = ny1
                n1b = ny2
                n2a = nz1
                n2b = nz2
                if (position == 'left' .and. ((nx1 <= 1 .and. 1 <= nx2) .or. (nx2 <= 1))) then
                    bindex = min(1, nx2)
                else if (position == 'right' .and. ((nx1 <= nx .and. nx <= nx2) .or. (nx1 >= nx))) then
                    bindex = max(nx, nx1)
                end if
            case ('y')
                n1 = nx2 - nx1 + 1
                n2 = nz2 - nz1 + 1
                n1a = nx1
                n1b = nx2
                n2a = nz1
                n2b = nz2
                if (position == 'front' .and. ((ny1 <= 1 .and. 1 <= ny2) .or. (ny2 <= 1))) then
                    bindex = min(1, ny2)
                else if (position == 'back' .and. ((ny1 <= ny .and. ny <= ny2) .or. (ny1 >= ny))) then
                    bindex = max(ny, ny1)
                end if
            case ('z')
                n1 = nx2 - nx1 + 1
                n2 = ny2 - ny1 + 1
                n1a = nx1
                n1b = nx2
                n2a = ny1
                n2b = ny2
                if (position == 'top' .and. ((nz1 <= 1 .and. 1 <= nz2) .or. (nz2 <= 1))) then
                    bindex = min(1, nz2)
                else if (position == 'bottom' .and. ((nz1 <= nz .and. nz <= nz2) .or. (nz1 >= nz))) then
                    bindex = max(nz, nz1)
                end if
        end select

        if (bindex /= 0) then

            ! Compute damping profile ratio only for selected points
            ! to save computational time
            interval = min(10, n1, n2)
            m1 = 1
            m2 = 1
            do i = 1, n1, interval
                m1 = m1 + 1
            end do
            do j = 1, n2, interval
                m2 = m2 + 1
            end do

            ! Here, m1/m2 can be overstepping outside of the block boundaries,
            ! and therefore we must check.
            if (n1a + (m1 - 1)*interval > n1b) then
                m1 = m1 - 1
            end if
            if (n2a + (m2 - 1)*interval > n2b) then
                m2 = m2 - 1
            end if

            call alloc_array(r, [1, m1, 1, m2])

            call warn(date_time_compact() &
                //' >> Computing MPML damping profile ratio for ' &
                //num2str(m1)//' x '//num2str(m2)//' points ')

            select case (axis)
                case ('x')
                    !$omp parallel do private(i, j, ii, jj) collapse(2) schedule(dynamic)
                    do jj = 1, m2
                        do ii = 1, m1
                            i = (ii - 1)*interval + n1a
                            j = (jj - 1)*interval + n2a
                            r(ii, jj) = damp_ratio( &
                                c11(bindex, i, j), c12(bindex, i, j), c13(bindex, i, j), c14(bindex, i, j), c15(bindex, i, j), c16(bindex, i, j), &
                                c22(bindex, i, j), c23(bindex, i, j), c24(bindex, i, j), c25(bindex, i, j), c26(bindex, i, j), &
                                c33(bindex, i, j), c34(bindex, i, j), c35(bindex, i, j), c36(bindex, i, j), &
                                c44(bindex, i, j), c45(bindex, i, j), c46(bindex, i, j), &
                                c55(bindex, i, j), c56(bindex, i, j), &
                                c66(bindex, i, j), 'x')
                        end do
                    end do
                    !$omp end parallel do
                case ('y')
                    !$omp parallel do private(i, j, ii, jj) collapse(2) schedule(dynamic)
                    do jj = 1, m2
                        do ii = 1, m1
                            i = (ii - 1)*interval + n1a
                            j = (jj - 1)*interval + n2a
                            r(ii, jj) = damp_ratio( &
                                c11(i, bindex, j), c12(i, bindex, j), c13(i, bindex, j), c14(i, bindex, j), c15(i, bindex, j), c16(i, bindex, j), &
                                c22(i, bindex, j), c23(i, bindex, j), c24(i, bindex, j), c25(i, bindex, j), c26(i, bindex, j), &
                                c33(i, bindex, j), c34(i, bindex, j), c35(i, bindex, j), c36(i, bindex, j), &
                                c44(i, bindex, j), c45(i, bindex, j), c46(i, bindex, j), &
                                c55(i, bindex, j), c56(i, bindex, j), &
                                c66(i, bindex, j), 'y')
                        end do
                    end do
                    !$omp end parallel do
                case ('z')
                    !$omp parallel do private(i, j, ii, jj) collapse(2) schedule(dynamic)
                    do jj = 1, m2
                        do ii = 1, m1
                            i = (ii - 1)*interval + n1a
                            j = (jj - 1)*interval + n2a
                            r(ii, jj) = damp_ratio( &
                                c11(i, j, bindex), c12(i, j, bindex), c13(i, j, bindex), c14(i, j, bindex), c15(i, j, bindex), c16(i, j, bindex), &
                                c22(i, j, bindex), c23(i, j, bindex), c24(i, j, bindex), c25(i, j, bindex), c26(i, j, bindex), &
                                c33(i, j, bindex), c34(i, j, bindex), c35(i, j, bindex), c36(i, j, bindex), &
                                c44(i, j, bindex), c45(i, j, bindex), c46(i, j, bindex), &
                                c55(i, j, bindex), c56(i, j, bindex), &
                                c66(i, j, bindex), 'z')
                        end do
                    end do
                    !$omp end parallel do
            end select

            ! Linear interpolate to the entire central region
            call alloc_array(ratio, [n1a, n1b, n2a, n2b], &
                source=interp(r, [m1, m2], [interval*1.0, interval*1.0], [0.0, 0.0], &
                [n1, n2], [1.0, 1.0], [0.0, 0.0], &
                ['linear', 'linear']))

            call warn(date_time_compact() &
                //' >> MPML damping ratio at '//tidy(position)//' boundary: ' &
                //num2str(minval(ratio), '(f6.3)')//' -- '//num2str(maxval(ratio), '(f6.3)'))

        else

            call alloc_array(ratio, [n1a, n1b, n2a, n2b])

        end if

    end subroutine compute_damping_profile_ratio_face

    !
    !> Compute damping coefficients for MPML
    !
    subroutine compute_damping_profile_ratio

        call compute_damping_profile_ratio_face(dampratio_left, 'x', 'left')
        call compute_damping_profile_ratio_face(dampratio_right, 'x', 'right')
        call compute_damping_profile_ratio_face(dampratio_front, 'y', 'front')
        call compute_damping_profile_ratio_face(dampratio_back, 'y', 'back')
        call compute_damping_profile_ratio_face(dampratio_top, 'z', 'top')
        call compute_damping_profile_ratio_face(dampratio_bottom, 'z', 'bottom')

    end subroutine compute_damping_profile_ratio

    !
    !> Compute CPML damping coefficients for 3D anisotropic media
    !
    subroutine compute_cfspml_damping_coef

        integer :: i, j, k
        real :: nd, ndh
        integer :: nxmin, nxmax, nymin, nymax, nzmin, nzmax
        real :: alpha_extra

        nxmin = 1
        nxmax = nx
        nymin = 1
        nymax = ny
        nzmin = 1
        nzmax = nz

        ! Compute CFS-MPML coefficients
        alphamax = 1.0*maxval(sgmtr%srcr(:)%f0)*const_pi
        R0 = 1.0e-5
        kmax = 1.0

        ! Diminishing alpha (Roden and Gedney, 2000); residual alpha to ensure long-time stability
        alpha_extra = 0.25*alphamax

        ! First Compute boundary layer P-wave velocity and damping ratio
        if (aniso_param == 'iso') then
            call alloc_array(dampratio_left, [ny1 - 1, ny2 + 1, nz1 - 1, nz2 + 1])
            call alloc_array(dampratio_right, [ny1 - 1, ny2 + 1, nz1 - 1, nz2 + 1])
            call alloc_array(dampratio_front, [nx1 - 1, nx2 + 1, nz1 - 1, nz2 + 1])
            call alloc_array(dampratio_back, [nx1 - 1, nx2 + 1, nz1 - 1, nz2 + 1])
            call alloc_array(dampratio_top, [nx1 - 1, nx2 + 1, ny1 - 1, ny2 + 1])
            call alloc_array(dampratio_bottom, [nx1 - 1, nx2 + 1, ny1 - 1, ny2 + 1])
            dampratio_left = 0
            dampratio_right = 0
            dampratio_front = 0
            dampratio_back = 0
            dampratio_top = 0
            dampratio_bottom = 0
        else
            call compute_damping_profile_ratio
        end if

        ! Second Compute damping coefficients for each CPML point
        call alloc_array(idax, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(iday, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(idaz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(hdax, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(hday, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(hdaz, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(alphaxi, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(alphayi, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(alphazi, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(alphaxh, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(alphayh, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(alphazh, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(kappaxi, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(kappayi, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(kappazi, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(kappaxh, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(kappayh, [nx1, nx2, ny1, ny2, nz1, nz2])
        call alloc_array(kappazh, [nx1, nx2, ny1, ny2, nz1, nz2])

        alphaxi = alphamax
        alphayi = alphamax
        alphazi = alphamax
        alphaxh = alphamax
        alphayh = alphamax
        alphazh = alphamax
        kappaxi = 1.0
        kappayi = 1.0
        kappazi = 1.0
        kappaxh = 1.0
        kappayh = 1.0
        kappazh = 1.0

        !$omp parallel do private(i, j, k) collapse(3)
        do k = nz1_interior, nz2_interior
            do j = ny1_interior, ny2_interior
                do i = nx1_interior, nx2_interior
                    alphaxi(i, j, k) = 0.0
                    alphayi(i, j, k) = 0.0
                    alphazi(i, j, k) = 0.0
                    alphaxh(i, j, k) = 0.0
                    alphayh(i, j, k) = 0.0
                    alphazh(i, j, k) = 0.0
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do private(i, j, k, nd, ndh) collapse(3) schedule(auto)
        do k = nz1, nz2
            do j = ny1, ny2
                do i = nx1, nx2

                    if (i <= nxmin) then
                        ! Left boundary

                        nd = abs(i - nxmin)*1.0/pml
                        ndh = abs(i - nxmin - 0.5)*1.0/pml
                        idax(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dx)*nd**npower_d
                        hdax(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dx)*ndh**npower_d
                        alphaxi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                        alphaxh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                        kappaxi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                        kappaxh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                    else if (i >= nxmax) then
                        ! Right boundary

                        nd = abs(i - nxmax)*1.0/pml
                        ndh = abs(i - nxmax - 0.5)*1.0/pml
                        idax(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dx)*nd**npower_d
                        hdax(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dx)*ndh**npower_d
                        alphaxi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                        alphaxh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                        kappaxi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                        kappaxh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                    end if

                    if (j <= nymin) then
                        ! Front boundary

                        nd = abs(j - nymin)*1.0/pml
                        ndh = abs(j - nymin - 0.5)*1.0/pml
                        iday(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dy)*nd**npower_d
                        hday(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dy)*ndh**npower_d
                        alphayi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                        alphayh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                        kappayi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                        kappayh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                    else if (j >= nymax) then
                        ! Back boundary

                        nd = abs(j - nymax)*1.0/pml
                        ndh = abs(j - nymax - 0.5)*1.0/pml
                        iday(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dy)*nd**npower_d
                        hday(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dy)*ndh**npower_d
                        alphayi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                        alphayh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                        kappayi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                        kappayh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                    end if

                    if (yn_free_surface) then
                        ! Free-surface modeling uses depth-varying mesh

                        if (k >= nzmax) then

                            ! Bottom boundary
                            nd = abs(k - nzmax)*1.0/pml
                            ndh = abs(k - nzmax - 0.5)*1.0/pml
                            idaz(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*eta_dz_i(k))*nd**npower_d
                            hdaz(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*eta_dz_i(k))*ndh**npower_d
                            alphazi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                            alphazh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                            kappazi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                            kappazh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                        end if

                    else

                        if (k <= nzmin) then
                            ! Top boundary

                            nd = abs(k - nzmin)*1.0/pml
                            ndh = abs(k - nzmin - 0.5)*1.0/pml
                            idaz(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dz)*nd**npower_d
                            hdaz(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dz)*ndh**npower_d
                            alphazi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                            alphazh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                            kappazi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                            kappazh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                        else if (k >= nzmax) then

                            ! Bottom boundary
                            nd = abs(k - nzmax)*1.0/pml
                            ndh = abs(k - nzmax - 0.5)*1.0/pml
                            idaz(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dz)*nd**npower_d
                            hdaz(i, j, k) = log(1.0/R0)*(npower_d + 1.0)*pmlvp/(2.0*pml*dz)*ndh**npower_d
                            alphazi(i, j, k) = max(alphamax*(1.0 - nd**npower_a), alpha_extra)
                            alphazh(i, j, k) = max(alphamax*(1.0 - ndh**npower_a), alpha_extra)
                            kappazi(i, j, k) = 1.0 + (kmax - 1.0)*nd**npower_k
                            kappazh(i, j, k) = 1.0 + (kmax - 1.0)*ndh**npower_k

                        end if

                    end if

                end do
            end do
        end do
        !$omp end parallel do

    end subroutine

end module
