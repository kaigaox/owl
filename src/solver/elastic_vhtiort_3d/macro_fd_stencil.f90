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

#define pdxxx_stencil ( \
coef81*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef82*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) + \
coef83*(stressxx(i + 3, j, k) - stressxx(i - 2, j, k)) + \
coef84*(stressxx(i + 4, j, k) - stressxx(i - 3, j, k)) + \
coef85*(stressxx(i + 5, j, k) - stressxx(i - 4, j, k)) + \
coef86*(stressxx(i + 6, j, k) - stressxx(i - 5, j, k)) + \
coef87*(stressxx(i + 7, j, k) - stressxx(i - 6, j, k)) + \
coef88*(stressxx(i + 8, j, k) - stressxx(i - 7, j, k)) \
)

#define pdxxy_stencil ( \
coef81*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef82*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) + \
coef83*(stressxy(i + 3, j + 1, k) - stressxy(i - 2, j + 1, k)) + \
coef84*(stressxy(i + 4, j + 1, k) - stressxy(i - 3, j + 1, k)) + \
coef85*(stressxy(i + 5, j + 1, k) - stressxy(i - 4, j + 1, k)) + \
coef86*(stressxy(i + 6, j + 1, k) - stressxy(i - 5, j + 1, k)) + \
coef87*(stressxy(i + 7, j + 1, k) - stressxy(i - 6, j + 1, k)) + \
coef88*(stressxy(i + 8, j + 1, k) - stressxy(i - 7, j + 1, k)) \
)

#define pdxxz_stencil ( \
coef81*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef82*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) + \
coef83*(stressxz(i + 3, j, k + 1) - stressxz(i - 2, j, k + 1)) + \
coef84*(stressxz(i + 4, j, k + 1) - stressxz(i - 3, j, k + 1)) + \
coef85*(stressxz(i + 5, j, k + 1) - stressxz(i - 4, j, k + 1)) + \
coef86*(stressxz(i + 6, j, k + 1) - stressxz(i - 5, j, k + 1)) + \
coef87*(stressxz(i + 7, j, k + 1) - stressxz(i - 6, j, k + 1)) + \
coef88*(stressxz(i + 8, j, k + 1) - stressxz(i - 7, j, k + 1)) \
)

#define pdyxy_stencil ( \
coef81*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef82*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) + \
coef83*(stressxy(i + 1, j + 3, k) - stressxy(i + 1, j - 2, k)) + \
coef84*(stressxy(i + 1, j + 4, k) - stressxy(i + 1, j - 3, k)) + \
coef85*(stressxy(i + 1, j + 5, k) - stressxy(i + 1, j - 4, k)) + \
coef86*(stressxy(i + 1, j + 6, k) - stressxy(i + 1, j - 5, k)) + \
coef87*(stressxy(i + 1, j + 7, k) - stressxy(i + 1, j - 6, k)) + \
coef88*(stressxy(i + 1, j + 8, k) - stressxy(i + 1, j - 7, k)) \
)

#define pdyyy_stencil ( \
coef81*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef82*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) + \
coef83*(stressyy(i, j + 3, k) - stressyy(i, j - 2, k)) + \
coef84*(stressyy(i, j + 4, k) - stressyy(i, j - 3, k)) + \
coef85*(stressyy(i, j + 5, k) - stressyy(i, j - 4, k)) + \
coef86*(stressyy(i, j + 6, k) - stressyy(i, j - 5, k)) + \
coef87*(stressyy(i, j + 7, k) - stressyy(i, j - 6, k)) + \
coef88*(stressyy(i, j + 8, k) - stressyy(i, j - 7, k)))

#define pdyyz_stencil ( \
coef81*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef82*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) + \
coef83*(stressyz(i, j + 3, k + 1) - stressyz(i, j - 2, k + 1)) + \
coef84*(stressyz(i, j + 4, k + 1) - stressyz(i, j - 3, k + 1)) + \
coef85*(stressyz(i, j + 5, k + 1) - stressyz(i, j - 4, k + 1)) + \
coef86*(stressyz(i, j + 6, k + 1) - stressyz(i, j - 5, k + 1)) + \
coef87*(stressyz(i, j + 7, k + 1) - stressyz(i, j - 6, k + 1)) + \
coef88*(stressyz(i, j + 8, k + 1) - stressyz(i, j - 7, k + 1)))

#define pdzxz_stencil ( \
coef81*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef82*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) + \
coef83*(stressxz(i + 1, j, k + 3) - stressxz(i + 1, j, k - 2)) + \
coef84*(stressxz(i + 1, j, k + 4) - stressxz(i + 1, j, k - 3)) + \
coef85*(stressxz(i + 1, j, k + 5) - stressxz(i + 1, j, k - 4)) + \
coef86*(stressxz(i + 1, j, k + 6) - stressxz(i + 1, j, k - 5)) + \
coef87*(stressxz(i + 1, j, k + 7) - stressxz(i + 1, j, k - 6)) + \
coef88*(stressxz(i + 1, j, k + 8) - stressxz(i + 1, j, k - 7)) \
)

#define pdzyz_stencil ( \
coef81*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef82*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) + \
coef83*(stressyz(i, j + 1, k + 3) - stressyz(i, j + 1, k - 2)) + \
coef84*(stressyz(i, j + 1, k + 4) - stressyz(i, j + 1, k - 3)) + \
coef85*(stressyz(i, j + 1, k + 5) - stressyz(i, j + 1, k - 4)) + \
coef86*(stressyz(i, j + 1, k + 6) - stressyz(i, j + 1, k - 5)) + \
coef87*(stressyz(i, j + 1, k + 7) - stressyz(i, j + 1, k - 6)) + \
coef88*(stressyz(i, j + 1, k + 8) - stressyz(i, j + 1, k - 7)) \
)

#define pdzzz_stencil ( \
coef81*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef82*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) + \
coef83*(stresszz(i, j, k + 3) - stresszz(i, j, k - 2)) + \
coef84*(stresszz(i, j, k + 4) - stresszz(i, j, k - 3)) + \
coef85*(stresszz(i, j, k + 5) - stresszz(i, j, k - 4)) + \
coef86*(stresszz(i, j, k + 6) - stresszz(i, j, k - 5)) + \
coef87*(stresszz(i, j, k + 7) - stresszz(i, j, k - 6)) + \
coef88*(stresszz(i, j, k + 8) - stresszz(i, j, k - 7)) \
)

#define pdxvx_stencil ( \
coef81*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef82*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef83*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef84*(vx(i + 4, j, k) - vx(i - 3, j, k)) + \
coef85*(vx(i + 5, j, k) - vx(i - 4, j, k)) + \
coef86*(vx(i + 6, j, k) - vx(i - 5, j, k)) + \
coef87*(vx(i + 7, j, k) - vx(i - 6, j, k)) + \
coef88*(vx(i + 8, j, k) - vx(i - 7, j, k)) \
)

#define pdyvx_stencil ( \
coef81*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef82*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) + \
coef83*(vx(i + 1, j + 3, k) - vx(i + 1, j - 2, k)) + \
coef84*(vx(i + 1, j + 4, k) - vx(i + 1, j - 3, k)) + \
coef85*(vx(i + 1, j + 5, k) - vx(i + 1, j - 4, k)) + \
coef86*(vx(i + 1, j + 6, k) - vx(i + 1, j - 5, k)) + \
coef87*(vx(i + 1, j + 7, k) - vx(i + 1, j - 6, k)) + \
coef88*(vx(i + 1, j + 8, k) - vx(i + 1, j - 7, k)) \
)

#define pdzvx_stencil ( \
coef81*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef82*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) + \
coef83*(vx(i + 1, j, k + 3) - vx(i + 1, j, k - 2)) + \
coef84*(vx(i + 1, j, k + 4) - vx(i + 1, j, k - 3)) + \
coef85*(vx(i + 1, j, k + 5) - vx(i + 1, j, k - 4)) + \
coef86*(vx(i + 1, j, k + 6) - vx(i + 1, j, k - 5)) + \
coef87*(vx(i + 1, j, k + 7) - vx(i + 1, j, k - 6)) + \
coef88*(vx(i + 1, j, k + 8) - vx(i + 1, j, k - 7)) \
)

#define pdxvy_stencil ( \
coef81*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef82*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) + \
coef83*(vy(i + 3, j + 1, k) - vy(i - 2, j + 1, k)) + \
coef84*(vy(i + 4, j + 1, k) - vy(i - 3, j + 1, k)) + \
coef85*(vy(i + 5, j + 1, k) - vy(i - 4, j + 1, k)) + \
coef86*(vy(i + 6, j + 1, k) - vy(i - 5, j + 1, k)) + \
coef87*(vy(i + 7, j + 1, k) - vy(i - 6, j + 1, k)) + \
coef88*(vy(i + 8, j + 1, k) - vy(i - 7, j + 1, k)) \
)

#define pdyvy_stencil ( \
coef81*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef82*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef83*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef84*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef85*(vy(i, j + 5, k) - vy(i, j - 4, k)) + \
coef86*(vy(i, j + 6, k) - vy(i, j - 5, k)) + \
coef87*(vy(i, j + 7, k) - vy(i, j - 6, k)) + \
coef88*(vy(i, j + 8, k) - vy(i, j - 7, k)) \
)

#define pdzvy_stencil ( \
coef81*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef82*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) + \
coef83*(vy(i, j + 1, k + 3) - vy(i, j + 1, k - 2)) + \
coef84*(vy(i, j + 1, k + 4) - vy(i, j + 1, k - 3)) + \
coef85*(vy(i, j + 1, k + 5) - vy(i, j + 1, k - 4)) + \
coef86*(vy(i, j + 1, k + 6) - vy(i, j + 1, k - 5)) + \
coef87*(vy(i, j + 1, k + 7) - vy(i, j + 1, k - 6)) + \
coef88*(vy(i, j + 1, k + 8) - vy(i, j + 1, k - 7)) \
)

#define pdxvz_stencil ( \
coef81*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef82*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) + \
coef83*(vz(i + 3, j, k + 1) - vz(i - 2, j, k + 1)) + \
coef84*(vz(i + 4, j, k + 1) - vz(i - 3, j, k + 1)) + \
coef85*(vz(i + 5, j, k + 1) - vz(i - 4, j, k + 1)) + \
coef86*(vz(i + 6, j, k + 1) - vz(i - 5, j, k + 1)) + \
coef87*(vz(i + 7, j, k + 1) - vz(i - 6, j, k + 1)) + \
coef88*(vz(i + 8, j, k + 1) - vz(i - 7, j, k + 1)) \
)

#define pdyvz_stencil ( \
coef81*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef82*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) + \
coef83*(vz(i, j + 3, k + 1) - vz(i, j - 2, k + 1)) + \
coef84*(vz(i, j + 4, k + 1) - vz(i, j - 3, k + 1)) + \
coef85*(vz(i, j + 5, k + 1) - vz(i, j - 4, k + 1)) + \
coef86*(vz(i, j + 6, k + 1) - vz(i, j - 5, k + 1)) + \
coef87*(vz(i, j + 7, k + 1) - vz(i, j - 6, k + 1)) + \
coef88*(vz(i, j + 8, k + 1) - vz(i, j - 7, k + 1)) \
)

#define pdzvz_stencil ( \
coef81*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef82*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef83*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef84*(vz(i, j, k + 4) - vz(i, j, k - 3)) + \
coef85*(vz(i, j, k + 5) - vz(i, j, k - 4)) + \
coef86*(vz(i, j, k + 6) - vz(i, j, k - 5)) + \
coef87*(vz(i, j, k + 7) - vz(i, j, k - 6)) + \
coef88*(vz(i, j, k + 8) - vz(i, j, k - 7)) \
)

#endif

! L = 7
#ifdef _fdorder14_

#define pdxxx_stencil ( \
coef71*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef72*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) + \
coef73*(stressxx(i + 3, j, k) - stressxx(i - 2, j, k)) + \
coef74*(stressxx(i + 4, j, k) - stressxx(i - 3, j, k)) + \
coef75*(stressxx(i + 5, j, k) - stressxx(i - 4, j, k)) + \
coef76*(stressxx(i + 6, j, k) - stressxx(i - 5, j, k)) + \
coef77*(stressxx(i + 7, j, k) - stressxx(i - 6, j, k)) \
)

#define pdyyy_stencil ( \
coef71*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef72*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) + \
coef73*(stressyy(i, j + 3, k) - stressyy(i, j - 2, k)) + \
coef74*(stressyy(i, j + 4, k) - stressyy(i, j - 3, k)) + \
coef75*(stressyy(i, j + 5, k) - stressyy(i, j - 4, k)) + \
coef76*(stressyy(i, j + 6, k) - stressyy(i, j - 5, k)) + \
coef77*(stressyy(i, j + 7, k) - stressyy(i, j - 6, k)) \
)

#define pdzzz_stencil ( \
coef71*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef72*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) + \
coef73*(stresszz(i, j, k + 3) - stresszz(i, j, k - 2)) + \
coef74*(stresszz(i, j, k + 4) - stresszz(i, j, k - 3)) + \
coef75*(stresszz(i, j, k + 5) - stresszz(i, j, k - 4)) + \
coef76*(stresszz(i, j, k + 6) - stresszz(i, j, k - 5)) + \
coef77*(stresszz(i, j, k + 7) - stresszz(i, j, k - 6)) \
)

#define pdxxy_stencil ( \
coef71*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef72*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) + \
coef73*(stressxy(i + 3, j + 1, k) - stressxy(i - 2, j + 1, k)) + \
coef74*(stressxy(i + 4, j + 1, k) - stressxy(i - 3, j + 1, k)) + \
coef75*(stressxy(i + 5, j + 1, k) - stressxy(i - 4, j + 1, k)) + \
coef76*(stressxy(i + 6, j + 1, k) - stressxy(i - 5, j + 1, k)) + \
coef77*(stressxy(i + 7, j + 1, k) - stressxy(i - 6, j + 1, k)) \
)

#define pdxxz_stencil ( \
coef71*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef72*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) + \
coef73*(stressxz(i + 3, j, k + 1) - stressxz(i - 2, j, k + 1)) + \
coef74*(stressxz(i + 4, j, k + 1) - stressxz(i - 3, j, k + 1)) + \
coef75*(stressxz(i + 5, j, k + 1) - stressxz(i - 4, j, k + 1)) + \
coef76*(stressxz(i + 6, j, k + 1) - stressxz(i - 5, j, k + 1)) + \
coef77*(stressxz(i + 7, j, k + 1) - stressxz(i - 6, j, k + 1)) \
)

#define pdyxy_stencil ( \
coef71*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef72*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) + \
coef73*(stressxy(i + 1, j + 3, k) - stressxy(i + 1, j - 2, k)) + \
coef74*(stressxy(i + 1, j + 4, k) - stressxy(i + 1, j - 3, k)) + \
coef75*(stressxy(i + 1, j + 5, k) - stressxy(i + 1, j - 4, k)) + \
coef76*(stressxy(i + 1, j + 6, k) - stressxy(i + 1, j - 5, k)) + \
coef77*(stressxy(i + 1, j + 7, k) - stressxy(i + 1, j - 6, k)) \
)

#define pdyyz_stencil ( \
coef71*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef72*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) + \
coef73*(stressyz(i, j + 3, k + 1) - stressyz(i, j - 2, k + 1)) + \
coef74*(stressyz(i, j + 4, k + 1) - stressyz(i, j - 3, k + 1)) + \
coef75*(stressyz(i, j + 5, k + 1) - stressyz(i, j - 4, k + 1)) + \
coef76*(stressyz(i, j + 6, k + 1) - stressyz(i, j - 5, k + 1)) + \
coef77*(stressyz(i, j + 7, k + 1) - stressyz(i, j - 6, k + 1)) \
)

#define pdzxz_stencil ( \
coef71*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef72*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) + \
coef73*(stressxz(i + 1, j, k + 3) - stressxz(i + 1, j, k - 2)) + \
coef74*(stressxz(i + 1, j, k + 4) - stressxz(i + 1, j, k - 3)) + \
coef75*(stressxz(i + 1, j, k + 5) - stressxz(i + 1, j, k - 4)) + \
coef76*(stressxz(i + 1, j, k + 6) - stressxz(i + 1, j, k - 5)) + \
coef77*(stressxz(i + 1, j, k + 7) - stressxz(i + 1, j, k - 6)) \
)

#define pdzyz_stencil ( \
coef71*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef72*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) + \
coef73*(stressyz(i, j + 1, k + 3) - stressyz(i, j + 1, k - 2)) + \
coef74*(stressyz(i, j + 1, k + 4) - stressyz(i, j + 1, k - 3)) + \
coef75*(stressyz(i, j + 1, k + 5) - stressyz(i, j + 1, k - 4)) + \
coef76*(stressyz(i, j + 1, k + 6) - stressyz(i, j + 1, k - 5)) + \
coef77*(stressyz(i, j + 1, k + 7) - stressyz(i, j + 1, k - 6)) \
)

#define pdxvx_stencil ( \
coef71*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef72*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef73*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef74*(vx(i + 4, j, k) - vx(i - 3, j, k)) + \
coef75*(vx(i + 5, j, k) - vx(i - 4, j, k)) + \
coef76*(vx(i + 6, j, k) - vx(i - 5, j, k)) + \
coef77*(vx(i + 7, j, k) - vx(i - 6, j, k)) \
)

#define pdyvx_stencil ( \
coef71*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef72*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) + \
coef73*(vx(i + 1, j + 3, k) - vx(i + 1, j - 2, k)) + \
coef74*(vx(i + 1, j + 4, k) - vx(i + 1, j - 3, k)) + \
coef75*(vx(i + 1, j + 5, k) - vx(i + 1, j - 4, k)) + \
coef76*(vx(i + 1, j + 6, k) - vx(i + 1, j - 5, k)) + \
coef77*(vx(i + 1, j + 7, k) - vx(i + 1, j - 6, k)) \
)

#define pdzvx_stencil ( \
coef71*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef72*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) + \
coef73*(vx(i + 1, j, k + 3) - vx(i + 1, j, k - 2)) + \
coef74*(vx(i + 1, j, k + 4) - vx(i + 1, j, k - 3)) + \
coef75*(vx(i + 1, j, k + 5) - vx(i + 1, j, k - 4)) + \
coef76*(vx(i + 1, j, k + 6) - vx(i + 1, j, k - 5)) + \
coef77*(vx(i + 1, j, k + 7) - vx(i + 1, j, k - 6)) \
)

#define pdxvy_stencil ( \
coef71*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef72*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) + \
coef73*(vy(i + 3, j + 1, k) - vy(i - 2, j + 1, k)) + \
coef74*(vy(i + 4, j + 1, k) - vy(i - 3, j + 1, k)) + \
coef75*(vy(i + 5, j + 1, k) - vy(i - 4, j + 1, k)) + \
coef76*(vy(i + 6, j + 1, k) - vy(i - 5, j + 1, k)) + \
coef77*(vy(i + 7, j + 1, k) - vy(i - 6, j + 1, k)) \
)

#define pdyvy_stencil ( \
coef71*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef72*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef73*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef74*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef75*(vy(i, j + 5, k) - vy(i, j - 4, k)) + \
coef76*(vy(i, j + 6, k) - vy(i, j - 5, k)) + \
coef77*(vy(i, j + 7, k) - vy(i, j - 6, k)) \
)

#define pdzvy_stencil ( \
coef71*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef72*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) + \
coef73*(vy(i, j + 1, k + 3) - vy(i, j + 1, k - 2)) + \
coef74*(vy(i, j + 1, k + 4) - vy(i, j + 1, k - 3)) + \
coef75*(vy(i, j + 1, k + 5) - vy(i, j + 1, k - 4)) + \
coef76*(vy(i, j + 1, k + 6) - vy(i, j + 1, k - 5)) + \
coef77*(vy(i, j + 1, k + 7) - vy(i, j + 1, k - 6)) \
)

#define pdxvz_stencil ( \
coef71*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef72*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) + \
coef73*(vz(i + 3, j, k + 1) - vz(i - 2, j, k + 1)) + \
coef74*(vz(i + 4, j, k + 1) - vz(i - 3, j, k + 1)) + \
coef75*(vz(i + 5, j, k + 1) - vz(i - 4, j, k + 1)) + \
coef76*(vz(i + 6, j, k + 1) - vz(i - 5, j, k + 1)) + \
coef77*(vz(i + 7, j, k + 1) - vz(i - 6, j, k + 1)) \
)

#define pdyvz_stencil ( \
coef71*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef72*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) + \
coef73*(vz(i, j + 3, k + 1) - vz(i, j - 2, k + 1)) + \
coef74*(vz(i, j + 4, k + 1) - vz(i, j - 3, k + 1)) + \
coef75*(vz(i, j + 5, k + 1) - vz(i, j - 4, k + 1)) + \
coef76*(vz(i, j + 6, k + 1) - vz(i, j - 5, k + 1)) + \
coef77*(vz(i, j + 7, k + 1) - vz(i, j - 6, k + 1)) \
)

#define pdzvz_stencil ( \
coef71*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef72*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef73*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef74*(vz(i, j, k + 4) - vz(i, j, k - 3)) + \
coef75*(vz(i, j, k + 5) - vz(i, j, k - 4)) + \
coef76*(vz(i, j, k + 6) - vz(i, j, k - 5)) + \
coef77*(vz(i, j, k + 7) - vz(i, j, k - 6)) \
)

#endif

! L = 6
#ifdef _fdorder12_

#define pdxxx_stencil ( \
coef61*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef62*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) + \
coef63*(stressxx(i + 3, j, k) - stressxx(i - 2, j, k)) + \
coef64*(stressxx(i + 4, j, k) - stressxx(i - 3, j, k)) + \
coef65*(stressxx(i + 5, j, k) - stressxx(i - 4, j, k)) + \
coef66*(stressxx(i + 6, j, k) - stressxx(i - 5, j, k)) \
)

#define pdyyy_stencil ( \
coef61*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef62*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) + \
coef63*(stressyy(i, j + 3, k) - stressyy(i, j - 2, k)) + \
coef64*(stressyy(i, j + 4, k) - stressyy(i, j - 3, k)) + \
coef65*(stressyy(i, j + 5, k) - stressyy(i, j - 4, k)) + \
coef66*(stressyy(i, j + 6, k) - stressyy(i, j - 5, k)) \
)

#define pdzzz_stencil ( \
coef61*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef62*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) + \
coef63*(stresszz(i, j, k + 3) - stresszz(i, j, k - 2)) + \
coef64*(stresszz(i, j, k + 4) - stresszz(i, j, k - 3)) + \
coef65*(stresszz(i, j, k + 5) - stresszz(i, j, k - 4)) + \
coef66*(stresszz(i, j, k + 6) - stresszz(i, j, k - 5)) \
)

#define pdxxy_stencil ( \
coef61*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef62*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) + \
coef63*(stressxy(i + 3, j + 1, k) - stressxy(i - 2, j + 1, k)) + \
coef64*(stressxy(i + 4, j + 1, k) - stressxy(i - 3, j + 1, k)) + \
coef65*(stressxy(i + 5, j + 1, k) - stressxy(i - 4, j + 1, k)) + \
coef66*(stressxy(i + 6, j + 1, k) - stressxy(i - 5, j + 1, k)) \
)

#define pdxxz_stencil ( \
coef61*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef62*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) + \
coef63*(stressxz(i + 3, j, k + 1) - stressxz(i - 2, j, k + 1)) + \
coef64*(stressxz(i + 4, j, k + 1) - stressxz(i - 3, j, k + 1)) + \
coef65*(stressxz(i + 5, j, k + 1) - stressxz(i - 4, j, k + 1)) + \
coef66*(stressxz(i + 6, j, k + 1) - stressxz(i - 5, j, k + 1)) \
)

#define pdyxy_stencil ( \
coef61*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef62*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) + \
coef63*(stressxy(i + 1, j + 3, k) - stressxy(i + 1, j - 2, k)) + \
coef64*(stressxy(i + 1, j + 4, k) - stressxy(i + 1, j - 3, k)) + \
coef65*(stressxy(i + 1, j + 5, k) - stressxy(i + 1, j - 4, k)) + \
coef66*(stressxy(i + 1, j + 6, k) - stressxy(i + 1, j - 5, k)) \
)

#define pdyyz_stencil ( \
coef61*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef62*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) + \
coef63*(stressyz(i, j + 3, k + 1) - stressyz(i, j - 2, k + 1)) + \
coef64*(stressyz(i, j + 4, k + 1) - stressyz(i, j - 3, k + 1)) + \
coef65*(stressyz(i, j + 5, k + 1) - stressyz(i, j - 4, k + 1)) + \
coef66*(stressyz(i, j + 6, k + 1) - stressyz(i, j - 5, k + 1)) \
)

#define pdzxz_stencil ( \
coef61*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef62*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) + \
coef63*(stressxz(i + 1, j, k + 3) - stressxz(i + 1, j, k - 2)) + \
coef64*(stressxz(i + 1, j, k + 4) - stressxz(i + 1, j, k - 3)) + \
coef65*(stressxz(i + 1, j, k + 5) - stressxz(i + 1, j, k - 4)) + \
coef66*(stressxz(i + 1, j, k + 6) - stressxz(i + 1, j, k - 5)) \
)

#define pdzyz_stencil ( \
coef61*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef62*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) + \
coef63*(stressyz(i, j + 1, k + 3) - stressyz(i, j + 1, k - 2)) + \
coef64*(stressyz(i, j + 1, k + 4) - stressyz(i, j + 1, k - 3)) + \
coef65*(stressyz(i, j + 1, k + 5) - stressyz(i, j + 1, k - 4)) + \
coef66*(stressyz(i, j + 1, k + 6) - stressyz(i, j + 1, k - 5)) \
)

#define pdxvx_stencil ( \
coef61*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef62*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef63*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef64*(vx(i + 4, j, k) - vx(i - 3, j, k)) + \
coef65*(vx(i + 5, j, k) - vx(i - 4, j, k)) + \
coef66*(vx(i + 6, j, k) - vx(i - 5, j, k)) \
)

#define pdyvx_stencil ( \
coef61*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef62*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) + \
coef63*(vx(i + 1, j + 3, k) - vx(i + 1, j - 2, k)) + \
coef64*(vx(i + 1, j + 4, k) - vx(i + 1, j - 3, k)) + \
coef65*(vx(i + 1, j + 5, k) - vx(i + 1, j - 4, k)) + \
coef66*(vx(i + 1, j + 6, k) - vx(i + 1, j - 5, k)) \
)

#define pdzvx_stencil ( \
coef61*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef62*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) + \
coef63*(vx(i + 1, j, k + 3) - vx(i + 1, j, k - 2)) + \
coef64*(vx(i + 1, j, k + 4) - vx(i + 1, j, k - 3)) + \
coef65*(vx(i + 1, j, k + 5) - vx(i + 1, j, k - 4)) + \
coef66*(vx(i + 1, j, k + 6) - vx(i + 1, j, k - 5)) \
)

#define pdxvy_stencil ( \
coef61*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef62*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) + \
coef63*(vy(i + 3, j + 1, k) - vy(i - 2, j + 1, k)) + \
coef64*(vy(i + 4, j + 1, k) - vy(i - 3, j + 1, k)) + \
coef65*(vy(i + 5, j + 1, k) - vy(i - 4, j + 1, k)) + \
coef66*(vy(i + 6, j + 1, k) - vy(i - 5, j + 1, k)) \
)

#define pdyvy_stencil ( \
coef61*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef62*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef63*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef64*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef65*(vy(i, j + 5, k) - vy(i, j - 4, k)) + \
coef66*(vy(i, j + 6, k) - vy(i, j - 5, k)) \
)

#define pdzvy_stencil ( \
coef61*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef62*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) + \
coef63*(vy(i, j + 1, k + 3) - vy(i, j + 1, k - 2)) + \
coef64*(vy(i, j + 1, k + 4) - vy(i, j + 1, k - 3)) + \
coef65*(vy(i, j + 1, k + 5) - vy(i, j + 1, k - 4)) + \
coef66*(vy(i, j + 1, k + 6) - vy(i, j + 1, k - 5)) \
)

#define pdxvz_stencil ( \
coef61*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef62*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) + \
coef63*(vz(i + 3, j, k + 1) - vz(i - 2, j, k + 1)) + \
coef64*(vz(i + 4, j, k + 1) - vz(i - 3, j, k + 1)) + \
coef65*(vz(i + 5, j, k + 1) - vz(i - 4, j, k + 1)) + \
coef66*(vz(i + 6, j, k + 1) - vz(i - 5, j, k + 1)) \
)

#define pdyvz_stencil ( \
coef61*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef62*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) + \
coef63*(vz(i, j + 3, k + 1) - vz(i, j - 2, k + 1)) + \
coef64*(vz(i, j + 4, k + 1) - vz(i, j - 3, k + 1)) + \
coef65*(vz(i, j + 5, k + 1) - vz(i, j - 4, k + 1)) + \
coef66*(vz(i, j + 6, k + 1) - vz(i, j - 5, k + 1)) \
)

#define pdzvz_stencil ( \
coef61*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef62*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef63*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef64*(vz(i, j, k + 4) - vz(i, j, k - 3)) + \
coef65*(vz(i, j, k + 5) - vz(i, j, k - 4)) + \
coef66*(vz(i, j, k + 6) - vz(i, j, k - 5)) \
)

#endif

! L = 5
#ifdef _fdorder10_

#define pdxxx_stencil ( \
coef51*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef52*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) + \
coef53*(stressxx(i + 3, j, k) - stressxx(i - 2, j, k)) + \
coef54*(stressxx(i + 4, j, k) - stressxx(i - 3, j, k)) + \
coef55*(stressxx(i + 5, j, k) - stressxx(i - 4, j, k)) \
)

#define pdyyy_stencil ( \
coef51*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef52*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) + \
coef53*(stressyy(i, j + 3, k) - stressyy(i, j - 2, k)) + \
coef54*(stressyy(i, j + 4, k) - stressyy(i, j - 3, k)) + \
coef55*(stressyy(i, j + 5, k) - stressyy(i, j - 4, k)) \
)

#define pdzzz_stencil ( \
coef51*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef52*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) + \
coef53*(stresszz(i, j, k + 3) - stresszz(i, j, k - 2)) + \
coef54*(stresszz(i, j, k + 4) - stresszz(i, j, k - 3)) + \
coef55*(stresszz(i, j, k + 5) - stresszz(i, j, k - 4)) \
)

#define pdxxy_stencil ( \
coef51*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef52*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) + \
coef53*(stressxy(i + 3, j + 1, k) - stressxy(i - 2, j + 1, k)) + \
coef54*(stressxy(i + 4, j + 1, k) - stressxy(i - 3, j + 1, k)) + \
coef55*(stressxy(i + 5, j + 1, k) - stressxy(i - 4, j + 1, k)) \
)

#define pdxxz_stencil ( \
coef51*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef52*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) + \
coef53*(stressxz(i + 3, j, k + 1) - stressxz(i - 2, j, k + 1)) + \
coef54*(stressxz(i + 4, j, k + 1) - stressxz(i - 3, j, k + 1)) + \
coef55*(stressxz(i + 5, j, k + 1) - stressxz(i - 4, j, k + 1)) \
)

#define pdyxy_stencil ( \
coef51*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef52*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) + \
coef53*(stressxy(i + 1, j + 3, k) - stressxy(i + 1, j - 2, k)) + \
coef54*(stressxy(i + 1, j + 4, k) - stressxy(i + 1, j - 3, k)) + \
coef55*(stressxy(i + 1, j + 5, k) - stressxy(i + 1, j - 4, k)) \
)

#define pdyyz_stencil ( \
coef51*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef52*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) + \
coef53*(stressyz(i, j + 3, k + 1) - stressyz(i, j - 2, k + 1)) + \
coef54*(stressyz(i, j + 4, k + 1) - stressyz(i, j - 3, k + 1)) + \
coef55*(stressyz(i, j + 5, k + 1) - stressyz(i, j - 4, k + 1)) \
)

#define pdzxz_stencil ( \
coef51*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef52*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) + \
coef53*(stressxz(i + 1, j, k + 3) - stressxz(i + 1, j, k - 2)) + \
coef54*(stressxz(i + 1, j, k + 4) - stressxz(i + 1, j, k - 3)) + \
coef55*(stressxz(i + 1, j, k + 5) - stressxz(i + 1, j, k - 4)) \
)

#define pdzyz_stencil ( \
coef51*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef52*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) + \
coef53*(stressyz(i, j + 1, k + 3) - stressyz(i, j + 1, k - 2)) + \
coef54*(stressyz(i, j + 1, k + 4) - stressyz(i, j + 1, k - 3)) + \
coef55*(stressyz(i, j + 1, k + 5) - stressyz(i, j + 1, k - 4)) \
)

#define pdxvx_stencil ( \
coef51*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef52*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef53*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef54*(vx(i + 4, j, k) - vx(i - 3, j, k)) + \
coef55*(vx(i + 5, j, k) - vx(i - 4, j, k)) \
)

#define pdyvx_stencil ( \
coef51*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef52*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) + \
coef53*(vx(i + 1, j + 3, k) - vx(i + 1, j - 2, k)) + \
coef54*(vx(i + 1, j + 4, k) - vx(i + 1, j - 3, k)) + \
coef55*(vx(i + 1, j + 5, k) - vx(i + 1, j - 4, k)) \
)

#define pdzvx_stencil ( \
coef51*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef52*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) + \
coef53*(vx(i + 1, j, k + 3) - vx(i + 1, j, k - 2)) + \
coef54*(vx(i + 1, j, k + 4) - vx(i + 1, j, k - 3)) + \
coef55*(vx(i + 1, j, k + 5) - vx(i + 1, j, k - 4)) \
)

#define pdxvy_stencil ( \
coef51*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef52*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) + \
coef53*(vy(i + 3, j + 1, k) - vy(i - 2, j + 1, k)) + \
coef54*(vy(i + 4, j + 1, k) - vy(i - 3, j + 1, k)) + \
coef55*(vy(i + 5, j + 1, k) - vy(i - 4, j + 1, k)) \
)

#define pdyvy_stencil ( \
coef51*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef52*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef53*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef54*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef55*(vy(i, j + 5, k) - vy(i, j - 4, k)) \
)

#define pdzvy_stencil ( \
coef51*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef52*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) + \
coef53*(vy(i, j + 1, k + 3) - vy(i, j + 1, k - 2)) + \
coef54*(vy(i, j + 1, k + 4) - vy(i, j + 1, k - 3)) + \
coef55*(vy(i, j + 1, k + 5) - vy(i, j + 1, k - 4)) \
)

#define pdxvz_stencil ( \
coef51*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef52*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) + \
coef53*(vz(i + 3, j, k + 1) - vz(i - 2, j, k + 1)) + \
coef54*(vz(i + 4, j, k + 1) - vz(i - 3, j, k + 1)) + \
coef55*(vz(i + 5, j, k + 1) - vz(i - 4, j, k + 1)) \
)

#define pdyvz_stencil ( \
coef51*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef52*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) + \
coef53*(vz(i, j + 3, k + 1) - vz(i, j - 2, k + 1)) + \
coef54*(vz(i, j + 4, k + 1) - vz(i, j - 3, k + 1)) + \
coef55*(vz(i, j + 5, k + 1) - vz(i, j - 4, k + 1)) \
)

#define pdzvz_stencil ( \
coef51*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef52*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef53*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef54*(vz(i, j, k + 4) - vz(i, j, k - 3)) + \
coef55*(vz(i, j, k + 5) - vz(i, j, k - 4)) \
)

#endif

! L = 4
#ifdef _fdorder8_

#define pdxxx_stencil ( \
coef41*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef42*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) + \
coef43*(stressxx(i + 3, j, k) - stressxx(i - 2, j, k)) + \
coef44*(stressxx(i + 4, j, k) - stressxx(i - 3, j, k)) \
)

#define pdyyy_stencil ( \
coef41*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef42*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) + \
coef43*(stressyy(i, j + 3, k) - stressyy(i, j - 2, k)) + \
coef44*(stressyy(i, j + 4, k) - stressyy(i, j - 3, k)) \
)

#define pdzzz_stencil ( \
coef41*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef42*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) + \
coef43*(stresszz(i, j, k + 3) - stresszz(i, j, k - 2)) + \
coef44*(stresszz(i, j, k + 4) - stresszz(i, j, k - 3)) \
)

#define pdxxy_stencil ( \
coef41*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef42*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) + \
coef43*(stressxy(i + 3, j + 1, k) - stressxy(i - 2, j + 1, k)) + \
coef44*(stressxy(i + 4, j + 1, k) - stressxy(i - 3, j + 1, k)) \
)

#define pdxxz_stencil ( \
coef41*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef42*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) + \
coef43*(stressxz(i + 3, j, k + 1) - stressxz(i - 2, j, k + 1)) + \
coef44*(stressxz(i + 4, j, k + 1) - stressxz(i - 3, j, k + 1)) \
)

#define pdyxy_stencil ( \
coef41*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef42*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) + \
coef43*(stressxy(i + 1, j + 3, k) - stressxy(i + 1, j - 2, k)) + \
coef44*(stressxy(i + 1, j + 4, k) - stressxy(i + 1, j - 3, k)) \
)

#define pdyyz_stencil ( \
coef41*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef42*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) + \
coef43*(stressyz(i, j + 3, k + 1) - stressyz(i, j - 2, k + 1)) + \
coef44*(stressyz(i, j + 4, k + 1) - stressyz(i, j - 3, k + 1)) \
)

#define pdzxz_stencil ( \
coef41*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef42*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) + \
coef43*(stressxz(i + 1, j, k + 3) - stressxz(i + 1, j, k - 2)) + \
coef44*(stressxz(i + 1, j, k + 4) - stressxz(i + 1, j, k - 3)) \
)

#define pdzyz_stencil ( \
coef41*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef42*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) + \
coef43*(stressyz(i, j + 1, k + 3) - stressyz(i, j + 1, k - 2)) + \
coef44*(stressyz(i, j + 1, k + 4) - stressyz(i, j + 1, k - 3)) \
)

#define pdxvx_stencil ( \
coef41*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef42*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef43*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef44*(vx(i + 4, j, k) - vx(i - 3, j, k)) \
)

#define pdyvx_stencil ( \
coef41*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef42*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) + \
coef43*(vx(i + 1, j + 3, k) - vx(i + 1, j - 2, k)) + \
coef44*(vx(i + 1, j + 4, k) - vx(i + 1, j - 3, k)) \
)

#define pdzvx_stencil ( \
coef41*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef42*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) + \
coef43*(vx(i + 1, j, k + 3) - vx(i + 1, j, k - 2)) + \
coef44*(vx(i + 1, j, k + 4) - vx(i + 1, j, k - 3)) \
)

#define pdxvy_stencil ( \
coef41*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef42*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) + \
coef43*(vy(i + 3, j + 1, k) - vy(i - 2, j + 1, k)) + \
coef44*(vy(i + 4, j + 1, k) - vy(i - 3, j + 1, k)) \
)

#define pdyvy_stencil ( \
coef41*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef42*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef43*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef44*(vy(i, j + 4, k) - vy(i, j - 3, k)) \
)

#define pdzvy_stencil ( \
coef41*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef42*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) + \
coef43*(vy(i, j + 1, k + 3) - vy(i, j + 1, k - 2)) + \
coef44*(vy(i, j + 1, k + 4) - vy(i, j + 1, k - 3)) \
)

#define pdxvz_stencil ( \
coef41*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef42*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) + \
coef43*(vz(i + 3, j, k + 1) - vz(i - 2, j, k + 1)) + \
coef44*(vz(i + 4, j, k + 1) - vz(i - 3, j, k + 1)) \
)

#define pdyvz_stencil ( \
coef41*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef42*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) + \
coef43*(vz(i, j + 3, k + 1) - vz(i, j - 2, k + 1)) + \
coef44*(vz(i, j + 4, k + 1) - vz(i, j - 3, k + 1)) \
)

#define pdzvz_stencil ( \
coef41*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef42*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef43*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef44*(vz(i, j, k + 4) - vz(i, j, k - 3)) \
)

#endif

! L = 3
#ifdef _fdorder6_

#define pdxxx_3 ( \
coef31*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef32*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) + \
coef33*(stressxx(i + 3, j, k) - stressxx(i - 2, j, k)) \
)

#define pdyyy_3 ( \
coef31*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef32*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) + \
coef33*(stressyy(i, j + 3, k) - stressyy(i, j - 2, k)) \
)

#define pdzzz_3 ( \
coef31*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef32*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) + \
coef33*(stresszz(i, j, k + 3) - stresszz(i, j, k - 2)) \
)

#define pdxxy_3 ( \
coef31*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef32*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) + \
coef33*(stressxy(i + 3, j + 1, k) - stressxy(i - 2, j + 1, k)) \
)

#define pdxxz_3 ( \
coef31*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef32*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) + \
coef33*(stressxz(i + 3, j, k + 1) - stressxz(i - 2, j, k + 1)) \
)

#define pdyxy_3 ( \
coef31*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef32*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) + \
coef33*(stressxy(i + 1, j + 3, k) - stressxy(i + 1, j - 2, k)) \
)

#define pdyyz_3 ( \
coef31*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef32*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) + \
coef33*(stressyz(i, j + 3, k + 1) - stressyz(i, j - 2, k + 1)) \
)

#define pdzxz_3 ( \
coef31*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef32*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) + \
coef33*(stressxz(i + 1, j, k + 3) - stressxz(i + 1, j, k - 2)) \
)

#define pdzyz_3 ( \
coef31*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef32*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) + \
coef33*(stressyz(i, j + 1, k + 3) - stressyz(i, j + 1, k - 2)) \
)

#define pdxvx_3 ( \
coef31*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef32*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef33*(vx(i + 3, j, k) - vx(i - 2, j, k)) \
)

#define pdyvx_3 ( \
coef31*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef32*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) + \
coef33*(vx(i + 1, j + 3, k) - vx(i + 1, j - 2, k)) \
)

#define pdzvx_3 ( \
coef31*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef32*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) + \
coef33*(vx(i + 1, j, k + 3) - vx(i + 1, j, k - 2)) \
)

#define pdxvy_3 ( \
coef31*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef32*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) + \
coef33*(vy(i + 3, j + 1, k) - vy(i - 2, j + 1, k)) \
)

#define pdyvy_3 ( \
coef31*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef32*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef33*(vy(i, j + 3, k) - vy(i, j - 2, k)) \
)

#define pdzvy_3 ( \
coef31*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef32*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) + \
coef33*(vy(i, j + 1, k + 3) - vy(i, j + 1, k - 2)) \
)

#define pdxvz_3 ( \
coef31*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef32*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) + \
coef33*(vz(i + 3, j, k + 1) - vz(i - 2, j, k + 1)) \
)

#define pdyvz_3 ( \
coef31*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef32*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) + \
coef33*(vz(i, j + 3, k + 1) - vz(i, j - 2, k + 1)) \
)

#define pdzvz_3 ( \
coef31*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef32*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef33*(vz(i, j, k + 3) - vz(i, j, k - 2)) \
)

#endif

! L = 2
#ifdef _fdorder4_

#define pdxxx_2 ( \
coef21*(stressxx(i + 1, j, k) - stressxx(i, j, k)) + \
coef22*(stressxx(i + 2, j, k) - stressxx(i - 1, j, k)) \
)

#define pdyyy_2 ( \
coef21*(stressyy(i, j + 1, k) - stressyy(i, j, k)) + \
coef22*(stressyy(i, j + 2, k) - stressyy(i, j - 1, k)) \
)

#define pdzzz_2 ( \
coef21*(stresszz(i, j, k + 1) - stresszz(i, j, k)) + \
coef22*(stresszz(i, j, k + 2) - stresszz(i, j, k - 1)) \
)

#define pdxxy_2 ( \
coef21*(stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k)) + \
coef22*(stressxy(i + 2, j + 1, k) - stressxy(i - 1, j + 1, k)) \
)

#define pdxxz_2 ( \
coef21*(stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1)) + \
coef22*(stressxz(i + 2, j, k + 1) - stressxz(i - 1, j, k + 1)) \
)

#define pdyxy_2 ( \
coef21*(stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k)) + \
coef22*(stressxy(i + 1, j + 2, k) - stressxy(i + 1, j - 1, k)) \
)

#define pdyyz_2 ( \
coef21*(stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1)) + \
coef22*(stressyz(i, j + 2, k + 1) - stressyz(i, j - 1, k + 1)) \
)

#define pdzxz_2 ( \
coef21*(stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k)) + \
coef22*(stressxz(i + 1, j, k + 2) - stressxz(i + 1, j, k - 1)) \
)

#define pdzyz_2 ( \
coef21*(stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k)) + \
coef22*(stressyz(i, j + 1, k + 2) - stressyz(i, j + 1, k - 1)) \
)

#define pdxvx_2 ( \
coef21*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef22*(vx(i + 2, j, k) - vx(i - 1, j, k)) \
)

#define pdyvx_2 ( \
coef21*(vx(i + 1, j + 1, k) - vx(i + 1, j, k)) + \
coef22*(vx(i + 1, j + 2, k) - vx(i + 1, j - 1, k)) \
)

#define pdzvx_2 ( \
coef21*(vx(i + 1, j, k + 1) - vx(i + 1, j, k)) + \
coef22*(vx(i + 1, j, k + 2) - vx(i + 1, j, k - 1)) \
)

#define pdxvy_2 ( \
coef21*(vy(i + 1, j + 1, k) - vy(i, j + 1, k)) + \
coef22*(vy(i + 2, j + 1, k) - vy(i - 1, j + 1, k)) \
)

#define pdyvy_2 ( \
coef21*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef22*(vy(i, j + 2, k) - vy(i, j - 1, k)) \
)

#define pdzvy_2 ( \
coef21*(vy(i, j + 1, k + 1) - vy(i, j + 1, k)) + \
coef22*(vy(i, j + 1, k + 2) - vy(i, j + 1, k - 1)) \
)

#define pdxvz_2 ( \
coef21*(vz(i + 1, j, k + 1) - vz(i, j, k + 1)) + \
coef22*(vz(i + 2, j, k + 1) - vz(i - 1, j, k + 1)) \
)

#define pdyvz_2 ( \
coef21*(vz(i, j + 1, k + 1) - vz(i, j, k + 1)) + \
coef22*(vz(i, j + 2, k + 1) - vz(i, j - 1, k + 1)) \
)

#define pdzvz_2 ( \
coef21*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef22*(vz(i, j, k + 2) - vz(i, j, k - 1)) \
)

#endif

! L = 1
#ifdef _fdorder2_

#define pdxxx_1 (stressxx(i + 1, j, k) - stressxx(i, j, k))

#define pdyyy_1 (stressyy(i, j + 1, k) - stressyy(i, j, k))

#define pdzzz_1 (stresszz(i, j, k + 1) - stresszz(i, j, k))

#define pdxxy_1 (stressxy(i + 1, j + 1, k) - stressxy(i, j + 1, k))

#define pdxxz_1 (stressxz(i + 1, j, k + 1) - stressxz(i, j, k + 1))

#define pdyxy_1 (stressxy(i + 1, j + 1, k) - stressxy(i + 1, j, k))

#define pdyyz_1 (stressyz(i, j + 1, k + 1) - stressyz(i, j, k + 1))

#define pdzxz_1 (stressxz(i + 1, j, k + 1) - stressxz(i + 1, j, k))

#define pdzyz_1 (stressyz(i, j + 1, k + 1) - stressyz(i, j + 1, k))

#define pdxvx_1 (vx(i + 1, j, k) - vx(i, j, k))

#define pdyvx_1 (vx(i + 1, j + 1, k) - vx(i + 1, j, k))

#define pdzvx_1 (vx(i + 1, j, k + 1) - vx(i + 1, j, k))

#define pdxvy_1 (vy(i + 1, j + 1, k) - vy(i, j + 1, k))

#define pdyvy_1 (vy(i, j + 1, k) - vy(i, j, k))

#define pdzvy_1 (vy(i, j + 1, k + 1) - vy(i, j + 1, k))

#define pdxvz_1 (vz(i + 1, j, k + 1) - vz(i, j, k + 1))

#define pdyvz_1 (vz(i, j + 1, k + 1) - vz(i, j, k + 1))

#define pdzvz_1 (vz(i, j, k + 1) - vz(i, j, k))

#endif
