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


! L = 8
#ifdef _fdorder16_

#define pdxw_stencil(w, i, j, k) (coef81*(w(i + 1, j, k) - w(i, j, k)) \
+ coef82*(w(i + 2, j, k) - w(i - 1, j, k)) \
+ coef83*(w(i + 3, j, k) - w(i - 2, j, k)) \
+ coef84*(w(i + 4, j, k) - w(i - 3, j, k)) \
+ coef85*(w(i + 5, j, k) - w(i - 4, j, k)) \
+ coef86*(w(i + 6, j, k) - w(i - 5, j, k)) \
+ coef87*(w(i + 7, j, k) - w(i - 6, j, k)) \
+ coef88*(w(i + 8, j, k) - w(i - 7, j, k)))

#define pdyw_stencil(w, i, j, k) (coef81*(w(i, j + 1, k) - w(i, j, k)) \
+ coef82*(w(i, j + 2, k) - w(i, j - 1, k)) \
+ coef83*(w(i, j + 3, k) - w(i, j - 2, k)) \
+ coef84*(w(i, j + 4, k) - w(i, j - 3, k)) \
+ coef85*(w(i, j + 5, k) - w(i, j - 4, k)) \
+ coef86*(w(i, j + 6, k) - w(i, j - 5, k)) \
+ coef87*(w(i, j + 7, k) - w(i, j - 6, k)) \
+ coef88*(w(i, j + 8, k) - w(i, j - 7, k)))

#define pdzw_stencil(w, i, j, k) (coef81*(w(i, j, k + 1) - w(i, j, k)) \
+ coef82*(w(i, j, k + 2) - w(i, j, k - 1)) \
+ coef83*(w(i, j, k + 3) - w(i, j, k - 2)) \
+ coef84*(w(i, j, k + 4) - w(i, j, k - 3)) \
+ coef85*(w(i, j, k + 5) - w(i, j, k - 4)) \
+ coef86*(w(i, j, k + 6) - w(i, j, k - 5)) \
+ coef87*(w(i, j, k + 7) - w(i, j, k - 6)) \
+ coef88*(w(i, j, k + 8) - w(i, j, k - 7)))

#endif

! L = 7
#ifdef _fdorder14_

#define pdxw_stencil(w, i, j, k) (coef71*(w(i + 1, j, k) - w(i, j, k)) \
+ coef72*(w(i + 2, j, k) - w(i - 1, j, k)) \
+ coef73*(w(i + 3, j, k) - w(i - 2, j, k)) \
+ coef74*(w(i + 4, j, k) - w(i - 3, j, k)) \
+ coef75*(w(i + 5, j, k) - w(i - 4, j, k)) \
+ coef76*(w(i + 6, j, k) - w(i - 5, j, k)) \
+ coef77*(w(i + 7, j, k) - w(i - 6, j, k)))

#define pdyw_stencil(w, i, j, k) (coef71*(w(i, j + 1, k) - w(i, j, k)) \
+ coef72*(w(i, j + 2, k) - w(i, j - 1, k)) \
+ coef73*(w(i, j + 3, k) - w(i, j - 2, k)) \
+ coef74*(w(i, j + 4, k) - w(i, j - 3, k)) \
+ coef75*(w(i, j + 5, k) - w(i, j - 4, k)) \
+ coef76*(w(i, j + 6, k) - w(i, j - 5, k)) \
+ coef77*(w(i, j + 7, k) - w(i, j - 6, k)))

#define pdzw_stencil(w, i, j, k) (coef71*(w(i, j, k + 1) - w(i, j, k)) \
+ coef72*(w(i, j, k + 2) - w(i, j, k - 1)) \
+ coef73*(w(i, j, k + 3) - w(i, j, k - 2)) \
+ coef74*(w(i, j, k + 4) - w(i, j, k - 3)) \
+ coef75*(w(i, j, k + 5) - w(i, j, k - 4)) \
+ coef76*(w(i, j, k + 6) - w(i, j, k - 5)) \
+ coef77*(w(i, j, k + 7) - w(i, j, k - 6)))

#endif

! L = 6
#ifdef _fdorder12_

#define pdxw_stencil(w, i, j, k) (coef61*(w(i + 1, j, k) - w(i, j, k)) \
+ coef62*(w(i + 2, j, k) - w(i - 1, j, k)) \
+ coef63*(w(i + 3, j, k) - w(i - 2, j, k)) \
+ coef64*(w(i + 4, j, k) - w(i - 3, j, k)) \
+ coef65*(w(i + 5, j, k) - w(i - 4, j, k)) \
+ coef66*(w(i + 6, j, k) - w(i - 5, j, k)))

#define pdyw_stencil(w, i, j, k) (coef61*(w(i, j + 1, k) - w(i, j, k)) \
+ coef62*(w(i, j + 2, k) - w(i, j - 1, k)) \
+ coef63*(w(i, j + 3, k) - w(i, j - 2, k)) \
+ coef64*(w(i, j + 4, k) - w(i, j - 3, k)) \
+ coef65*(w(i, j + 5, k) - w(i, j - 4, k)) \
+ coef66*(w(i, j + 6, k) - w(i, j - 5, k)))

#define pdzw_stencil(w, i, j, k) (coef61*(w(i, j, k + 1) - w(i, j, k)) \
+ coef62*(w(i, j, k + 2) - w(i, j, k - 1)) \
+ coef63*(w(i, j, k + 3) - w(i, j, k - 2)) \
+ coef64*(w(i, j, k + 4) - w(i, j, k - 3)) \
+ coef65*(w(i, j, k + 5) - w(i, j, k - 4)) \
+ coef66*(w(i, j, k + 6) - w(i, j, k - 5)))

#endif

! L = 5
#ifdef _fdorder10_

#define pdxw_stencil(w, i, j, k) (coef51*(w(i + 1, j, k) - w(i, j, k)) \
+ coef52*(w(i + 2, j, k) - w(i - 1, j, k)) \
+ coef53*(w(i + 3, j, k) - w(i - 2, j, k)) \
+ coef54*(w(i + 4, j, k) - w(i - 3, j, k)) \
+ coef55*(w(i + 5, j, k) - w(i - 4, j, k)))

#define pdyw_stencil(w, i, j, k) (coef51*(w(i, j + 1, k) - w(i, j, k)) \
+ coef52*(w(i, j + 2, k) - w(i, j - 1, k)) \
+ coef53*(w(i, j + 3, k) - w(i, j - 2, k)) \
+ coef54*(w(i, j + 4, k) - w(i, j - 3, k)) \
+ coef55*(w(i, j + 5, k) - w(i, j - 4, k)))

#define pdzw_stencil(w, i, j, k) (coef51*(w(i, j, k + 1) - w(i, j, k)) \
+ coef52*(w(i, j, k + 2) - w(i, j, k - 1)) \
+ coef53*(w(i, j, k + 3) - w(i, j, k - 2)) \
+ coef54*(w(i, j, k + 4) - w(i, j, k - 3)) \
+ coef55*(w(i, j, k + 5) - w(i, j, k - 4)))

#endif

! L = 4
#ifdef _fdorder8_

#define pdxw_stencil(w, i, j, k) (coef41*(w(i + 1, j, k) - w(i, j, k)) \
+ coef42*(w(i + 2, j, k) - w(i - 1, j, k)) \
+ coef43*(w(i + 3, j, k) - w(i - 2, j, k)) \
+ coef44*(w(i + 4, j, k) - w(i - 3, j, k)))

#define pdyw_stencil(w, i, j, k) (coef41*(w(i, j + 1, k) - w(i, j, k)) \
+ coef42*(w(i, j + 2, k) - w(i, j - 1, k)) \
+ coef43*(w(i, j + 3, k) - w(i, j - 2, k)) \
+ coef44*(w(i, j + 4, k) - w(i, j - 3, k)))

#define pdzw_stencil(w, i, j, k) (coef41*(w(i, j, k + 1) - w(i, j, k)) \
+ coef42*(w(i, j, k + 2) - w(i, j, k - 1)) \
+ coef43*(w(i, j, k + 3) - w(i, j, k - 2)) \
+ coef44*(w(i, j, k + 4) - w(i, j, k - 3)))

#endif

! L = 3
#ifdef _fdorder6_

#define pdxw_3(w, i, j, k) (coef31*(w(i + 1, j, k) - w(i, j, k)) \
+ coef32*(w(i + 2, j, k) - w(i - 1, j, k)) \
+ coef33*(w(i + 3, j, k) - w(i - 2, j, k)))

#define pdyw_3(w, i, j, k) (coef31*(w(i, j + 1, k) - w(i, j, k)) \
+ coef32*(w(i, j + 2, k) - w(i, j - 1, k)) \
+ coef33*(w(i, j + 3, k) - w(i, j - 2, k)))

#define pdzw_3(w, i, j, k) (coef31*(w(i, j, k + 1) - w(i, j, k)) \
+ coef32*(w(i, j, k + 2) - w(i, j, k - 1)) \
+ coef33*(w(i, j, k + 3) - w(i, j, k - 2)))

#endif

! L = 2
#ifdef _fdorder4_

#define pdxw_2(w, i, j, k) (coef21*(w(i + 1, j, k) - w(i, j, k)) \
+ coef22*(w(i + 2, j, k) - w(i - 1, j, k)))

#define pdyw_2(w, i, j, k) (coef21*(w(i, j + 1, k) - w(i, j, k)) \
+ coef22*(w(i, j + 2, k) - w(i, j - 1, k)))

#define pdzw_2(w, i, j, k) (coef21*(w(i, j, k + 1) - w(i, j, k)) \
+ coef22*(w(i, j, k + 2) - w(i, j, k - 1)))

#endif

! L = 1
#ifdef _fdorder2_

#define pdxw_1(w, i, j, k) (w(i + 1, j, k) - w(i, j, k))

#define pdyw_1(w, i, j, k) (w(i, j + 1, k) - w(i, j, k))

#define pdzw_1(w, i, j, k) (w(i, j, k + 1) - w(i, j, k))

#endif
