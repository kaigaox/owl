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

#define pdxw_stencil(w, i, j) (coef81*(w(i + 1, j) - w(i, j)) \
+coef82*(w(i + 2, j) - w(i - 1, j)) \
+coef83*(w(i + 3, j) - w(i - 2, j)) \
+coef84*(w(i + 4, j) - w(i - 3, j)) \
+coef85*(w(i + 5, j) - w(i - 4, j)) \
+coef86*(w(i + 6, j) - w(i - 5, j)) \
+coef87*(w(i + 7, j) - w(i - 6, j)) \
+coef88*(w(i + 8, j) - w(i - 7, j)))

#define pdzw_stencil(w, i, j) (coef81*(w(i, j + 1) - w(i, j)) \
+coef82*(w(i, j + 2) - w(i, j - 1)) \
+coef83*(w(i, j + 3) - w(i, j - 2)) \
+coef84*(w(i, j + 4) - w(i, j - 3)) \
+coef85*(w(i, j + 5) - w(i, j - 4)) \
+coef86*(w(i, j + 6) - w(i, j - 5)) \
+coef87*(w(i, j + 7) - w(i, j - 6)) \
+coef88*(w(i, j + 8) - w(i, j - 7)))

#endif

! L=7
#ifdef _fdorder14_

#define pdxw_stencil(w, i, j) (coef71*(w(i + 1, j) - w(i, j)) \
+coef72*(w(i + 2, j) - w(i - 1, j)) \
+coef73*(w(i + 3, j) - w(i - 2, j)) \
+coef74*(w(i + 4, j) - w(i - 3, j)) \
+coef75*(w(i + 5, j) - w(i - 4, j)) \
+coef76*(w(i + 6, j) - w(i - 5, j)) \
+coef77*(w(i + 7, j) - w(i - 6, j)))

#define pdzw_stencil(w, i, j) (coef71*(w(i, j + 1) - w(i, j)) \
+coef72*(w(i, j + 2) - w(i, j - 1)) \
+coef73*(w(i, j + 3) - w(i, j - 2)) \
+coef74*(w(i, j + 4) - w(i, j - 3)) \
+coef75*(w(i, j + 5) - w(i, j - 4)) \
+coef76*(w(i, j + 6) - w(i, j - 5)) \
+coef77*(w(i, j + 7) - w(i, j - 6)))

#endif

! L=6
#ifdef _fdorder12_

#define pdxw_stencil(w, i, j) (coef61*(w(i + 1, j) - w(i, j)) \
+coef62*(w(i + 2, j) - w(i - 1, j)) \
+coef63*(w(i + 3, j) - w(i - 2, j)) \
+coef64*(w(i + 4, j) - w(i - 3, j)) \
+coef65*(w(i + 5, j) - w(i - 4, j)) \
+coef66*(w(i + 6, j) - w(i - 5, j)))

#define pdzw_stencil(w, i, j) (coef61*(w(i, j + 1) - w(i, j)) \
+coef62*(w(i, j + 2) - w(i, j - 1)) \
+coef63*(w(i, j + 3) - w(i, j - 2)) \
+coef64*(w(i, j + 4) - w(i, j - 3)) \
+coef65*(w(i, j + 5) - w(i, j - 4)) \
+coef66*(w(i, j + 6) - w(i, j - 5)))

#endif

! L=5
#ifdef _fdorder10_

#define pdxw_stencil(w, i, j) (coef51*(w(i + 1, j) - w(i, j)) \
+coef52*(w(i + 2, j) - w(i - 1, j)) \
+coef53*(w(i + 3, j) - w(i - 2, j)) \
+coef54*(w(i + 4, j) - w(i - 3, j)) \
+coef55*(w(i + 5, j) - w(i - 4, j)))

#define pdzw_stencil(w, i, j) (coef51*(w(i, j + 1) - w(i, j)) \
+coef52*(w(i, j + 2) - w(i, j - 1)) \
+coef53*(w(i, j + 3) - w(i, j - 2)) \
+coef54*(w(i, j + 4) - w(i, j - 3)) \
+coef55*(w(i, j + 5) - w(i, j - 4)))

#endif

! L=4
#ifdef _fdorder8_

#define pdxw_stencil(w, i, j) (coef41*(w(i + 1, j) - w(i, j)) \
+coef42*(w(i + 2, j) - w(i - 1, j)) \
+coef43*(w(i + 3, j) - w(i - 2, j)) \
+coef44*(w(i + 4, j) - w(i - 3, j)))

#define pdzw_stencil(w, i, j) (coef41*(w(i, j + 1) - w(i, j)) \
+coef42*(w(i, j + 2) - w(i, j - 1)) \
+coef43*(w(i, j + 3) - w(i, j - 2)) \
+coef44*(w(i, j + 4) - w(i, j - 3)))

#endif

! L=3
#ifdef _fdorder6_

#define pdxw_3(w, i, j) (coef31*(w(i + 1, j) - w(i, j)) \
+coef32*(w(i + 2, j) - w(i - 1, j)) \
+coef33*(w(i + 3, j) - w(i - 2, j)))

#define pdzw_3(w, i, j) (coef31*(w(i, j + 1) - w(i, j)) \
+coef32*(w(i, j + 2) - w(i, j - 1)) \
+coef33*(w(i, j + 3) - w(i, j - 2)))

#endif

! L=2
#ifdef _fdorder4_

#define pdxw_2(w, i, j) (coef21*(w(i + 1, j) - w(i, j)) \
+coef22*(w(i + 2, j) - w(i - 1, j)))

#define pdzw_2(w, i, j) (coef21*(w(i, j + 1) - w(i, j)) \
+coef22*(w(i, j + 2) - w(i, j - 1)))

#endif

! L=1
#ifdef _fdorder2_

#define pdxw_1(w, i, j) (w(i + 1, j) - w(i, j))

#define pdzw_1(w, i, j) (w(i, j + 1) - w(i, j))

#endif
