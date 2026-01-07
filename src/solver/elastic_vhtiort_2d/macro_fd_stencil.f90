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

#define pdzvx_stencil (coef81*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef82*(vx(i + 1, j + 2) - vx(i + 1, j - 1)) \
+ coef83*(vx(i + 1, j + 3) - vx(i + 1, j - 2)) \
+ coef84*(vx(i + 1, j + 4) - vx(i + 1, j - 3)) \
+ coef85*(vx(i + 1, j + 5) - vx(i + 1, j - 4)) \
+ coef86*(vx(i + 1, j + 6) - vx(i + 1, j - 5)) \
+ coef87*(vx(i + 1, j + 7) - vx(i + 1, j - 6)) \
+ coef88*(vx(i + 1, j + 8) - vx(i + 1, j - 7)))

#define pdxvz_stencil (coef81*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef82*(vz(i + 2, j + 1) - vz(i - 1, j + 1)) \
+ coef83*(vz(i + 3, j + 1) - vz(i - 2, j + 1)) \
+ coef84*(vz(i + 4, j + 1) - vz(i - 3, j + 1)) \
+ coef85*(vz(i + 5, j + 1) - vz(i - 4, j + 1)) \
+ coef86*(vz(i + 6, j + 1) - vz(i - 5, j + 1)) \
+ coef87*(vz(i + 7, j + 1) - vz(i - 6, j + 1)) \
+ coef88*(vz(i + 8, j + 1) - vz(i - 7, j + 1)))

#define pdxxx_stencil (coef81*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef82*(stressxx(i + 2, j) - stressxx(i - 1, j)) \
+ coef83*(stressxx(i + 3, j) - stressxx(i - 2, j)) \
+ coef84*(stressxx(i + 4, j) - stressxx(i - 3, j)) \
+ coef85*(stressxx(i + 5, j) - stressxx(i - 4, j)) \
+ coef86*(stressxx(i + 6, j) - stressxx(i - 5, j)) \
+ coef87*(stressxx(i + 7, j) - stressxx(i - 6, j)) \
+ coef88*(stressxx(i + 8, j) - stressxx(i - 7, j)))

#define pdxxz_stencil (coef81*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef82*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)) \
+ coef83*(stressxz(i + 3, j + 1) - stressxz(i - 2, j + 1)) \
+ coef84*(stressxz(i + 4, j + 1) - stressxz(i - 3, j + 1)) \
+ coef85*(stressxz(i + 5, j + 1) - stressxz(i - 4, j + 1)) \
+ coef86*(stressxz(i + 6, j + 1) - stressxz(i - 5, j + 1)) \
+ coef87*(stressxz(i + 7, j + 1) - stressxz(i - 6, j + 1)) \
+ coef88*(stressxz(i + 8, j + 1) - stressxz(i - 7, j + 1)))

#define pdzzz_stencil (coef81*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef82*(stresszz(i, j + 2) - stresszz(i, j - 1)) \
+ coef83*(stresszz(i, j + 3) - stresszz(i, j - 2)) \
+ coef84*(stresszz(i, j + 4) - stresszz(i, j - 3)) \
+ coef85*(stresszz(i, j + 5) - stresszz(i, j - 4)) \
+ coef86*(stresszz(i, j + 6) - stresszz(i, j - 5)) \
+ coef87*(stresszz(i, j + 7) - stresszz(i, j - 6)) \
+ coef88*(stresszz(i, j + 8) - stresszz(i, j - 7)))

#define pdzxz_stencil (coef81*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef82*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)) \
+ coef83*(stressxz(i + 1, j + 3) - stressxz(i + 1, j - 2)) \
+ coef84*(stressxz(i + 1, j + 4) - stressxz(i + 1, j - 3)) \
+ coef85*(stressxz(i + 1, j + 5) - stressxz(i + 1, j - 4)) \
+ coef86*(stressxz(i + 1, j + 6) - stressxz(i + 1, j - 5)) \
+ coef87*(stressxz(i + 1, j + 7) - stressxz(i + 1, j - 6)) \
+ coef88*(stressxz(i + 1, j + 8) - stressxz(i + 1, j - 7)))

#endif

! L=7
#ifdef _fdorder14_

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

#define pdzvx_stencil (coef71*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef72*(vx(i + 1, j + 2) - vx(i + 1, j - 1)) \
+ coef73*(vx(i + 1, j + 3) - vx(i + 1, j - 2)) \
+ coef74*(vx(i + 1, j + 4) - vx(i + 1, j - 3)) \
+ coef75*(vx(i + 1, j + 5) - vx(i + 1, j - 4)) \
+ coef76*(vx(i + 1, j + 6) - vx(i + 1, j - 5)) \
+ coef77*(vx(i + 1, j + 7) - vx(i + 1, j - 6)))

#define pdxvz_stencil (coef71*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef72*(vz(i + 2, j + 1) - vz(i - 1, j + 1)) \
+ coef73*(vz(i + 3, j + 1) - vz(i - 2, j + 1)) \
+ coef74*(vz(i + 4, j + 1) - vz(i - 3, j + 1)) \
+ coef75*(vz(i + 5, j + 1) - vz(i - 4, j + 1)) \
+ coef76*(vz(i + 6, j + 1) - vz(i - 5, j + 1)) \
+ coef77*(vz(i + 7, j + 1) - vz(i - 6, j + 1)))

#define pdxxx_stencil (coef71*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef72*(stressxx(i + 2, j) - stressxx(i - 1, j)) \
+ coef73*(stressxx(i + 3, j) - stressxx(i - 2, j)) \
+ coef74*(stressxx(i + 4, j) - stressxx(i - 3, j)) \
+ coef75*(stressxx(i + 5, j) - stressxx(i - 4, j)) \
+ coef76*(stressxx(i + 6, j) - stressxx(i - 5, j)) \
+ coef77*(stressxx(i + 7, j) - stressxx(i - 6, j)))

#define pdxxz_stencil (coef71*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef72*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)) \
+ coef73*(stressxz(i + 3, j + 1) - stressxz(i - 2, j + 1)) \
+ coef74*(stressxz(i + 4, j + 1) - stressxz(i - 3, j + 1)) \
+ coef75*(stressxz(i + 5, j + 1) - stressxz(i - 4, j + 1)) \
+ coef76*(stressxz(i + 6, j + 1) - stressxz(i - 5, j + 1)) \
+ coef77*(stressxz(i + 7, j + 1) - stressxz(i - 6, j + 1)))

#define pdzzz_stencil (coef71*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef72*(stresszz(i, j + 2) - stresszz(i, j - 1)) \
+ coef73*(stresszz(i, j + 3) - stresszz(i, j - 2)) \
+ coef74*(stresszz(i, j + 4) - stresszz(i, j - 3)) \
+ coef75*(stresszz(i, j + 5) - stresszz(i, j - 4)) \
+ coef76*(stresszz(i, j + 6) - stresszz(i, j - 5)) \
+ coef77*(stresszz(i, j + 7) - stresszz(i, j - 6)))

#define pdzxz_stencil (coef71*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef72*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)) \
+ coef73*(stressxz(i + 1, j + 3) - stressxz(i + 1, j - 2)) \
+ coef74*(stressxz(i + 1, j + 4) - stressxz(i + 1, j - 3)) \
+ coef75*(stressxz(i + 1, j + 5) - stressxz(i + 1, j - 4)) \
+ coef76*(stressxz(i + 1, j + 6) - stressxz(i + 1, j - 5)) \
+ coef77*(stressxz(i + 1, j + 7) - stressxz(i + 1, j - 6)))

#endif

! L=6
#ifdef _fdorder12_

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

#define pdzvx_stencil (coef61*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef62*(vx(i + 1, j + 2) - vx(i + 1, j - 1)) \
+ coef63*(vx(i + 1, j + 3) - vx(i + 1, j - 2)) \
+ coef64*(vx(i + 1, j + 4) - vx(i + 1, j - 3)) \
+ coef65*(vx(i + 1, j + 5) - vx(i + 1, j - 4)) \
+ coef66*(vx(i + 1, j + 6) - vx(i + 1, j - 5)))

#define pdxvz_stencil (coef61*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef62*(vz(i + 2, j + 1) - vz(i - 1, j + 1)) \
+ coef63*(vz(i + 3, j + 1) - vz(i - 2, j + 1)) \
+ coef64*(vz(i + 4, j + 1) - vz(i - 3, j + 1)) \
+ coef65*(vz(i + 5, j + 1) - vz(i - 4, j + 1)) \
+ coef66*(vz(i + 6, j + 1) - vz(i - 5, j + 1)))

#define pdxxx_stencil (coef61*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef62*(stressxx(i + 2, j) - stressxx(i - 1, j)) \
+ coef63*(stressxx(i + 3, j) - stressxx(i - 2, j)) \
+ coef64*(stressxx(i + 4, j) - stressxx(i - 3, j)) \
+ coef65*(stressxx(i + 5, j) - stressxx(i - 4, j)) \
+ coef66*(stressxx(i + 6, j) - stressxx(i - 5, j)))

#define pdxxz_stencil (coef61*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef62*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)) \
+ coef63*(stressxz(i + 3, j + 1) - stressxz(i - 2, j + 1)) \
+ coef64*(stressxz(i + 4, j + 1) - stressxz(i - 3, j + 1)) \
+ coef65*(stressxz(i + 5, j + 1) - stressxz(i - 4, j + 1)) \
+ coef66*(stressxz(i + 6, j + 1) - stressxz(i - 5, j + 1)))

#define pdzzz_stencil (coef61*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef62*(stresszz(i, j + 2) - stresszz(i, j - 1)) \
+ coef63*(stresszz(i, j + 3) - stresszz(i, j - 2)) \
+ coef64*(stresszz(i, j + 4) - stresszz(i, j - 3)) \
+ coef65*(stresszz(i, j + 5) - stresszz(i, j - 4)) \
+ coef66*(stresszz(i, j + 6) - stresszz(i, j - 5)))

#define pdzxz_stencil (coef61*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef62*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)) \
+ coef63*(stressxz(i + 1, j + 3) - stressxz(i + 1, j - 2)) \
+ coef64*(stressxz(i + 1, j + 4) - stressxz(i + 1, j - 3)) \
+ coef65*(stressxz(i + 1, j + 5) - stressxz(i + 1, j - 4)) \
+ coef66*(stressxz(i + 1, j + 6) - stressxz(i + 1, j - 5)))

#endif

! L=5
#ifdef _fdorder10_

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

#define pdzvx_stencil (coef51*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef52*(vx(i + 1, j + 2) - vx(i + 1, j - 1)) \
+ coef53*(vx(i + 1, j + 3) - vx(i + 1, j - 2)) \
+ coef54*(vx(i + 1, j + 4) - vx(i + 1, j - 3)) \
+ coef55*(vx(i + 1, j + 5) - vx(i + 1, j - 4)))

#define pdxvz_stencil (coef51*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef52*(vz(i + 2, j + 1) - vz(i - 1, j + 1)) \
+ coef53*(vz(i + 3, j + 1) - vz(i - 2, j + 1)) \
+ coef54*(vz(i + 4, j + 1) - vz(i - 3, j + 1)) \
+ coef55*(vz(i + 5, j + 1) - vz(i - 4, j + 1)))

#define pdxxx_stencil (coef51*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef52*(stressxx(i + 2, j) - stressxx(i - 1, j)) \
+ coef53*(stressxx(i + 3, j) - stressxx(i - 2, j)) \
+ coef54*(stressxx(i + 4, j) - stressxx(i - 3, j)) \
+ coef55*(stressxx(i + 5, j) - stressxx(i - 4, j)))

#define pdxxz_stencil (coef51*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef52*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)) \
+ coef53*(stressxz(i + 3, j + 1) - stressxz(i - 2, j + 1)) \
+ coef54*(stressxz(i + 4, j + 1) - stressxz(i - 3, j + 1)) \
+ coef55*(stressxz(i + 5, j + 1) - stressxz(i - 4, j + 1)))

#define pdzzz_stencil (coef51*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef52*(stresszz(i, j + 2) - stresszz(i, j - 1)) \
+ coef53*(stresszz(i, j + 3) - stresszz(i, j - 2)) \
+ coef54*(stresszz(i, j + 4) - stresszz(i, j - 3)) \
+ coef55*(stresszz(i, j + 5) - stresszz(i, j - 4)))

#define pdzxz_stencil (coef51*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef52*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)) \
+ coef53*(stressxz(i + 1, j + 3) - stressxz(i + 1, j - 2)) \
+ coef54*(stressxz(i + 1, j + 4) - stressxz(i + 1, j - 3)) \
+ coef55*(stressxz(i + 1, j + 5) - stressxz(i + 1, j - 4)))

#endif

! L=4
#ifdef _fdorder8_

#define pdxvx_stencil (coef41*(vx(i + 1, j) - vx(i, j)) \
+ coef42*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef43*(vx(i + 3, j) - vx(i - 2, j)) \
+ coef44*(vx(i + 4, j) - vx(i - 3, j)))

#define pdzvz_stencil (coef41*(vz(i, j + 1) - vz(i, j)) \
+ coef42*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef43*(vz(i, j + 3) - vz(i, j - 2)) \
+ coef44*(vz(i, j + 4) - vz(i, j - 3)))

#define pdzvx_stencil (coef41*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef42*(vx(i + 1, j + 2) - vx(i + 1, j - 1)) \
+ coef43*(vx(i + 1, j + 3) - vx(i + 1, j - 2)) \
+ coef44*(vx(i + 1, j + 4) - vx(i + 1, j - 3)))

#define pdxvz_stencil (coef41*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef42*(vz(i + 2, j + 1) - vz(i - 1, j + 1)) \
+ coef43*(vz(i + 3, j + 1) - vz(i - 2, j + 1)) \
+ coef44*(vz(i + 4, j + 1) - vz(i - 3, j + 1)))

#define pdxxx_stencil (coef41*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef42*(stressxx(i + 2, j) - stressxx(i - 1, j)) \
+ coef43*(stressxx(i + 3, j) - stressxx(i - 2, j)) \
+ coef44*(stressxx(i + 4, j) - stressxx(i - 3, j)))

#define pdxxz_stencil (coef41*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef42*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)) \
+ coef43*(stressxz(i + 3, j + 1) - stressxz(i - 2, j + 1)) \
+ coef44*(stressxz(i + 4, j + 1) - stressxz(i - 3, j + 1)))

#define pdzzz_stencil (coef41*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef42*(stresszz(i, j + 2) - stresszz(i, j - 1)) \
+ coef43*(stresszz(i, j + 3) - stresszz(i, j - 2)) \
+ coef44*(stresszz(i, j + 4) - stresszz(i, j - 3)))

#define pdzxz_stencil (coef41*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef42*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)) \
+ coef43*(stressxz(i + 1, j + 3) - stressxz(i + 1, j - 2)) \
+ coef44*(stressxz(i + 1, j + 4) - stressxz(i + 1, j - 3)))

#endif

! L=3
#ifdef _fdorder6_

#define pdxvx_3 (coef31*(vx(i + 1, j) - vx(i, j)) \
+ coef32*(vx(i + 2, j) - vx(i - 1, j)) \
+ coef33*(vx(i + 3, j) - vx(i - 2, j)))

#define pdzvz_3 (coef31*(vz(i, j + 1) - vz(i, j)) \
+ coef32*(vz(i, j + 2) - vz(i, j - 1)) \
+ coef33*(vz(i, j + 3) - vz(i, j - 2)))

#define pdzvx_3 (coef31*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef32*(vx(i + 1, j + 2) - vx(i + 1, j - 1)) \
+ coef33*(vx(i + 1, j + 3) - vx(i + 1, j - 2)))

#define pdxvz_3 (coef31*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef32*(vz(i + 2, j + 1) - vz(i - 1, j + 1)) \
+ coef33*(vz(i + 3, j + 1) - vz(i - 2, j + 1)))

#define pdxxx_3 (coef31*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef32*(stressxx(i + 2, j) - stressxx(i - 1, j)) \
+ coef33*(stressxx(i + 3, j) - stressxx(i - 2, j)))

#define pdxxz_3 (coef31*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef32*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)) \
+ coef33*(stressxz(i + 3, j + 1) - stressxz(i - 2, j + 1)))

#define pdzzz_3 (coef31*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef32*(stresszz(i, j + 2) - stresszz(i, j - 1)) \
+ coef33*(stresszz(i, j + 3) - stresszz(i, j - 2)))

#define pdzxz_3 (coef31*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef32*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)) \
+ coef33*(stressxz(i + 1, j + 3) - stressxz(i + 1, j - 2)))

#endif

! L=2
#ifdef _fdorder4_

#define pdxvx_2 (coef21*(vx(i + 1, j) - vx(i, j)) \
+ coef22*(vx(i + 2, j) - vx(i - 1, j)))

#define pdzvz_2 (coef21*(vz(i, j + 1) - vz(i, j)) \
+ coef22*(vz(i, j + 2) - vz(i, j - 1)))

#define pdzvx_2 (coef21*(vx(i + 1, j + 1) - vx(i + 1, j)) \
+ coef22*(vx(i + 1, j + 2) - vx(i + 1, j - 1)))

#define pdxvz_2 (coef21*(vz(i + 1, j + 1) - vz(i, j + 1)) \
+ coef22*(vz(i + 2, j + 1) - vz(i - 1, j + 1)))

#define pdxxx_2 (coef21*(stressxx(i + 1, j) - stressxx(i, j)) \
+ coef22*(stressxx(i + 2, j) - stressxx(i - 1, j)))

#define pdxxz_2 (coef21*(stressxz(i + 1, j + 1) - stressxz(i, j + 1)) \
+ coef22*(stressxz(i + 2, j + 1) - stressxz(i - 1, j + 1)))

#define pdzzz_2 (coef21*(stresszz(i, j + 1) - stresszz(i, j)) \
+ coef22*(stresszz(i, j + 2) - stresszz(i, j - 1)))

#define pdzxz_2 (coef21*(stressxz(i + 1, j + 1) - stressxz(i + 1, j)) \
+ coef22*(stressxz(i + 1, j + 2) - stressxz(i + 1, j - 1)))

#endif

! L=1
#ifdef _fdorder2_

#define pdxvx_1 (vx(i + 1, j) - vx(i, j))

#define pdzvz_1 (vz(i, j + 1) - vz(i, j))

#define pdzvx_1 (vx(i + 1, j + 1) - vx(i + 1, j))

#define pdxvz_1 (vz(i + 1, j + 1) - vz(i, j + 1))

#define pdxxx_1 (stressxx(i + 1, j) - stressxx(i, j))

#define pdxxz_1 (stressxz(i + 1, j + 1) - stressxz(i, j + 1))

#define pdzzz_1 (stresszz(i, j + 1) - stresszz(i, j))

#define pdzxz_1 (stressxz(i + 1, j + 1) - stressxz(i + 1, j))

#endif
