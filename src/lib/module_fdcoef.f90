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


module mod_fdcoef

    implicit none

    ! L=8
    real, parameter :: coef81 = 0.1259312e+1
    real, parameter :: coef82 = -0.1280347e+0
    real, parameter :: coef83 = 0.3841945e-1
    real, parameter :: coef84 = -0.1473229e-1
    real, parameter :: coef85 = 0.5924913e-2
    real, parameter :: coef86 = -0.2248618e-2
    real, parameter :: coef87 = 0.7179226e-3
    real, parameter :: coef88 = -0.1400855e-3

    ! L=7
    real, parameter :: coef71 = 0.1254799e+1
    real, parameter :: coef72 = -0.1238928e+0
    real, parameter :: coef73 = 0.3494371e-1
    real, parameter :: coef74 = -0.1208897e-1
    real, parameter :: coef75 = 0.4132531e-2
    real, parameter :: coef76 = -0.1197110e-2
    real, parameter :: coef77 = 0.2122227e-3

    ! L=6
    real, parameter :: coef61 = 0.1247662e+1
    real, parameter :: coef62 = -0.1175538e+0
    real, parameter :: coef63 = 0.2997970e-1
    real, parameter :: coef64 = -0.8719078e-2
    real, parameter :: coef65 = 0.2215897e-2
    real, parameter :: coef66 = -0.3462075e-3

    ! L=5
    real, parameter :: coef51 = 0.1236607e+1
    real, parameter :: coef52 = -0.1082265e+0
    real, parameter :: coef53 = 0.2343440e-1
    real, parameter :: coef54 = -0.5033546e-2
    real, parameter :: coef55 = 0.6817483e-3

    ! L=4
    real, parameter :: coef41 = 0.1217990e+1
    real, parameter :: coef42 = -0.9382142e-1
    real, parameter :: coef43 = 0.1507536e-1
    real, parameter :: coef44 = -0.1700324e-2

    ! L=3
    real, parameter :: coef31 = 0.1185991e+1
    real, parameter :: coef32 = -0.7249965e-1
    real, parameter :: coef33 = 0.6301572e-2

    ! L=2
    real, parameter :: coef21 = 0.1129042e+1
    real, parameter :: coef22 = -0.4301412e-1

    ! FD coefficients in the array form
#ifdef _fdorder16_
    integer, parameter :: fdhalf = 8
    real, parameter, dimension(1:8) :: fdcoefs = &
        [0.1259312e+1, -0.1280347e+0, 0.3841945e-1, -0.1473229e-1, &
        0.5924913e-2, -0.2248618e-2, 0.7179226e-3, -0.1400855e-3]
#endif

#ifdef _fdorder14_
    integer, parameter :: fdhalf = 7
    real, parameter, dimension(1:7) :: fdcoefs = &
        [0.1254799e+1, -0.1238928e+0, 0.3494371e-1, -0.1208897e-1, &
        0.4132531e-2, -0.1197110e-2, 0.2122227e-3]
#endif

#ifdef _fdorder12_
    integer, parameter :: fdhalf = 6
    real, parameter, dimension(1:6) :: fdcoefs = &
        [0.1247662e+1, -0.1175538e+0, 0.2997970e-1, -0.8719078e-2, &
        0.2215897e-2, -0.3462075e-3]
#endif

#ifdef _fdorder10_
    integer, parameter :: fdhalf = 5
    real, parameter, dimension(1:5) :: fdcoefs = &
        [0.1236607e+1, -0.1082265e+0, 0.2343440e-1, -0.5033546e-2, 0.6817483e-3]
#endif

#ifdef _fdorder8_
    integer, parameter :: fdhalf = 4
    real, parameter, dimension(1:4) :: fdcoefs = &
        [0.1217990e+1, -0.9382142e-1, 0.1507536e-1, -0.1700324e-2]
#endif

#ifdef _fdorder6_
    integer, parameter :: fdhalf = 3
    real, parameter, dimension(1:3) :: fdcoefs = &
        [0.1185991e+1, -0.7249965e-1, 0.6301572e-2]
#endif

#ifdef _fdorder4_
    integer, parameter :: fdhalf = 2
    real, parameter, dimension(1:2) :: fdcoefs = &
        [0.1129042e+1, -0.4301412e-1]
#endif

#ifdef _fdorder2_
    integer, parameter :: fdhalf = 1
    real, parameter, dimension(1:1) :: fcoefs = [0.5]
#endif

end module
