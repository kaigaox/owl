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


! L=8
#ifdef _fdorder16_

#define pdxvx_stencil (coef81*(vx(i + 1, j) - vx(i, j)) \
+ coef82*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef83*(vx(i + 3, j) - vx(i - 2, j)) \
+ coef84*(vx(i + 4, j) - vx(i - 3, j)) \
+ coef85*(vx(i + 5, j) - vx(i - 4, j)) \
+ coef86*(vx(i + 6, j) - vx(i - 5, j)) \
+ coef87*(vx(i + 7, j) - vx(i - 6, j)) \
+ coef88*(vx(i + 8, j) - vx(i - 7, j)))

#define pdzvz_stencil (coef81*(vz(i, j + 1) - vz(i, j)) \
+ coef82*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef83*(vz(i, j + 3) - vz(i, j - 2)) \
+ coef84*(vz(i, j + 4) - vz(i, j - 3)) \
+ coef85*(vz(i, j + 5) - vz(i, j - 4)) \
+ coef86*(vz(i, j + 6) - vz(i, j - 5)) \
+ coef87*(vz(i, j + 7) - vz(i, j - 6)) \
+ coef88*(vz(i, j + 8) - vz(i, j - 7)))

#define pdxp_stencil (coef81*(p(i + 1, j) - p(i, j)) \
+ coef82*(p(i + 2, j) - p(i - 1, j)) \
+ coef83*(p(i + 3, j) - p(i - 2, j)) \
+ coef84*(p(i + 4, j) - p(i - 3, j)) \
+ coef85*(p(i + 5, j) - p(i - 4, j)) \
+ coef86*(p(i + 6, j) - p(i - 5, j)) \
+ coef87*(p(i + 7, j) - p(i - 6, j)) \
+ coef88*(p(i + 8, j) - p(i - 7, j)))

#define pdzp_stencil (coef81*(p(i, j + 1) - p(i, j)) \
+ coef82*(p(i, j + 2) - p(i, j - 1)) \
+ coef83*(p(i, j + 3) - p(i, j - 2)) \
+ coef84*(p(i, j + 4) - p(i, j - 3)) \
+ coef85*(p(i, j + 5) - p(i, j - 4)) \
+ coef86*(p(i, j + 6) - p(i, j - 5)) \
+ coef87*(p(i, j + 7) - p(i, j - 6)) \
+ coef88*(p(i, j + 8) - p(i, j - 7)))

#endif

! L=7
#ifdef _fdorder14_

#define pdxp_stencil (coef71*(p(i + 1, j) - p(i, j)) \
+ coef72*(p(i + 2, j) - p(i - 1, j)) \
+ coef73*(p(i + 3, j) - p(i - 2, j)) \
+ coef74*(p(i + 4, j) - p(i - 3, j)) \
+ coef75*(p(i + 5, j) - p(i - 4, j)) \
+ coef76*(p(i + 6, j) - p(i - 5, j)) \
+ coef77*(p(i + 7, j) - p(i - 6, j)))

#define pdzp_stencil (coef71*(p(i, j + 1) - p(i, j)) \
+ coef72*(p(i, j + 2) - p(i, j - 1)) \
+ coef73*(p(i, j + 3) - p(i, j - 2)) \
+ coef74*(p(i, j + 4) - p(i, j - 3)) \
+ coef75*(p(i, j + 5) - p(i, j - 4)) \
+ coef76*(p(i, j + 6) - p(i, j - 5)) \
+ coef77*(p(i, j + 7) - p(i, j - 6)))

#define pdxvx_stencil (coef71*(vx(i + 1, j) - vx(i, j)) \
+ coef72*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef73*(vx(i + 3, j) - vx(i - 2, j)) \
+ coef74*(vx(i + 4, j) - vx(i - 3, j)) \
+ coef75*(vx(i + 5, j) - vx(i - 4, j)) \
+ coef76*(vx(i + 6, j) - vx(i - 5, j)) \
+ coef77*(vx(i + 7, j) - vx(i - 6, j)))

#define pdzvz_stencil (coef71*(vz(i, j + 1) - vz(i, j)) \
+ coef72*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef73*(vz(i, j + 3) - vz(i, j - 2)) \
+ coef74*(vz(i, j + 4) - vz(i, j - 3)) \
+ coef75*(vz(i, j + 5) - vz(i, j - 4)) \
+ coef76*(vz(i, j + 6) - vz(i, j - 5)) \
+ coef77*(vz(i, j + 7) - vz(i, j - 6)))

#endif

! L=6
#ifdef _fdorder12_

#define pdxp_stencil (coef61*(p(i + 1, j) - p(i, j)) \
+ coef62*(p(i + 2, j) - p(i - 1, j)) \
+ coef63*(p(i + 3, j) - p(i - 2, j)) \
+ coef64*(p(i + 4, j) - p(i - 3, j)) \
+ coef65*(p(i + 5, j) - p(i - 4, j)) \
+ coef66*(p(i + 6, j) - p(i - 5, j)))

#define pdzp_stencil (coef61*(p(i, j + 1) - p(i, j)) \
+ coef62*(p(i, j + 2) - p(i, j - 1)) \
+ coef63*(p(i, j + 3) - p(i, j - 2)) \
+ coef64*(p(i, j + 4) - p(i, j - 3)) \
+ coef65*(p(i, j + 5) - p(i, j - 4)) \
+ coef66*(p(i, j + 6) - p(i, j - 5)))

#define pdxvx_stencil (coef61*(vx(i + 1, j) - vx(i, j)) \
+ coef62*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef63*(vx(i + 3, j) - vx(i - 2, j)) \
+ coef64*(vx(i + 4, j) - vx(i - 3, j)) \
+ coef65*(vx(i + 5, j) - vx(i - 4, j)) \
+ coef66*(vx(i + 6, j) - vx(i - 5, j)))

#define pdzvz_stencil (coef61*(vz(i, j + 1) - vz(i, j)) \
+ coef62*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef63*(vz(i, j + 3) - vz(i, j - 2)) \
+ coef64*(vz(i, j + 4) - vz(i, j - 3)) \
+ coef65*(vz(i, j + 5) - vz(i, j - 4)) \
+ coef66*(vz(i, j + 6) - vz(i, j - 5)))

#endif

! L=5
#ifdef _fdorder10_

#define pdxp_stencil (coef51*(p(i + 1, j) - p(i, j)) \
+ coef52*(p(i + 2, j) - p(i - 1, j)) \
+ coef53*(p(i + 3, j) - p(i - 2, j)) \
+ coef54*(p(i + 4, j) - p(i - 3, j)) \
+ coef55*(p(i + 5, j) - p(i - 4, j)))

#define pdzp_stencil (coef51*(p(i, j + 1) - p(i, j)) \
+ coef52*(p(i, j + 2) - p(i, j - 1)) \
+ coef53*(p(i, j + 3) - p(i, j - 2)) \
+ coef54*(p(i, j + 4) - p(i, j - 3)) \
+ coef55*(p(i, j + 5) - p(i, j - 4)))

#define pdxvx_stencil (coef51*(vx(i + 1, j) - vx(i, j)) \
+ coef52*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef53*(vx(i + 3, j) - vx(i - 2, j)) \
+ coef54*(vx(i + 4, j) - vx(i - 3, j)) \
+ coef55*(vx(i + 5, j) - vx(i - 4, j)))

#define pdzvz_stencil (coef51*(vz(i, j + 1) - vz(i, j)) \
+ coef52*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef53*(vz(i, j + 3) - vz(i, j - 2)) \
+ coef54*(vz(i, j + 4) - vz(i, j - 3)) \
+ coef55*(vz(i, j + 5) - vz(i, j - 4)))

#endif

! L=4
#ifdef _fdorder8_

#define pdxp_stencil (coef41*(p(i + 1, j) - p(i, j)) \
+ coef42*(p(i + 2, j) - p(i - 1, j)) \
+ coef43*(p(i + 3, j) - p(i - 2, j)) \
+ coef44*(p(i + 4, j) - p(i - 3, j)))

#define pdzp_stencil (coef41*(p(i, j + 1) - p(i, j)) \
+ coef42*(p(i, j + 2) - p(i, j - 1)) \
+ coef43*(p(i, j + 3) - p(i, j - 2)) \
+ coef44*(p(i, j + 4) - p(i, j - 3)))

#define pdxvx_stencil (coef41*(vx(i + 1, j) - vx(i, j)) \
+ coef42*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef43*(vx(i + 3, j) - vx(i - 2, j)) \
+ coef44*(vx(i + 4, j) - vx(i - 3, j)))

#define pdzvz_stencil (coef41*(vz(i, j + 1) - vz(i, j)) \
+ coef42*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef43*(vz(i, j + 3) - vz(i, j - 2)) \
+ coef44*(vz(i, j + 4) - vz(i, j - 3)))

#endif

! L=3
#ifdef _fdorder6_

#define pdxp_stencil (coef31*(p(i + 1, j) - p(i, j)) \
+ coef32*(p(i + 2, j) - p(i - 1, j)) \
+ coef33*(p(i + 3, j) - p(i - 2, j)))

#define pdzp_stencil (coef31*(p(i, j + 1) - p(i, j)) \
+ coef32*(p(i, j + 2) - p(i, j - 1)) \
+ coef33*(p(i, j + 3) - p(i, j - 2)))

#define pdxvx_stencil (coef31*(vx(i + 1, j) - vx(i, j)) \
+ coef32*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef33*(vx(i + 3, j) - vx(i - 2, j)))

#define pdzvz_stencil (coef31*(vz(i, j + 1) - vz(i, j)) \
+ coef32*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef33*(vz(i, j + 3) - vz(i, j - 2)))

#endif

! L=2
#ifdef _fdorder4_

#define pdxp_stencil (coef21*(p(i + 1, j) - p(i, j)) \
+ coef22*(p(i + 2, j) - p(i - 1, j)))

#define pdzp_stencil (coef21*(p(i, j + 1) - p(i, j)) \
+ coef22*(p(i, j + 2) - p(i, j - 1)))

#define pdxvx_stencil (coef21*(vx(i + 1, j) - vx(i, j)) \
+ coef22*(vx(i + 2, j) - vx(i - 1, j)))

#define pdzvz_stencil (coef21*(vz(i, j + 1) - vz(i, j)) \
+ coef22*(vz(i, j + 2) - vz(i, j - 1)))

#endif

! L=1
#ifdef _fdorder2_

#define pdxp_stencil (p(i + 1, j) - p(i, j))

#define pdzp_stencil (p(i, j + 1) - p(i, j))

#define pdxvx_stencil (vx(i + 1, j) - vx(i, j))

#define pdzvz_stencil (vz(i, j + 1) - vz(i, j))

#endif
