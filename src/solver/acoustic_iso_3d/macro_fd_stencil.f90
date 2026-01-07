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

#define pdxp_stencil ( \
coef81*(p(i + 1, j, k) - p(i, j, k)) + \
coef82*(p(i + 2, j, k) - p(i - 1, j, k)) + \
coef83*(p(i + 3, j, k) - p(i - 2, j, k)) + \
coef84*(p(i + 4, j, k) - p(i - 3, j, k)) + \
coef85*(p(i + 5, j, k) - p(i - 4, j, k)) + \
coef86*(p(i + 6, j, k) - p(i - 5, j, k)) + \
coef87*(p(i + 7, j, k) - p(i - 6, j, k)) + \
coef88*(p(i + 8, j, k) - p(i - 7, j, k)) \
)

#define pdyp_stencil ( \
coef81*(p(i, j + 1, k) - p(i, j, k)) + \
coef82*(p(i, j + 2, k) - p(i, j - 1, k)) + \
coef83*(p(i, j + 3, k) - p(i, j - 2, k)) + \
coef84*(p(i, j + 4, k) - p(i, j - 3, k)) + \
coef85*(p(i, j + 5, k) - p(i, j - 4, k)) + \
coef86*(p(i, j + 6, k) - p(i, j - 5, k)) + \
coef87*(p(i, j + 7, k) - p(i, j - 6, k)) + \
coef88*(p(i, j + 8, k) - p(i, j - 7, k)))

#define pdzp_stencil ( \
coef81*(p(i, j, k + 1) - p(i, j, k)) + \
coef82*(p(i, j, k + 2) - p(i, j, k - 1)) + \
coef83*(p(i, j, k + 3) - p(i, j, k - 2)) + \
coef84*(p(i, j, k + 4) - p(i, j, k - 3)) + \
coef85*(p(i, j, k + 5) - p(i, j, k - 4)) + \
coef86*(p(i, j, k + 6) - p(i, j, k - 5)) + \
coef87*(p(i, j, k + 7) - p(i, j, k - 6)) + \
coef88*(p(i, j, k + 8) - p(i, j, k - 7)) \
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

! L=7
#ifdef _fdorder14_

#define pdxp_stencil ( \
coef71*(p(i + 1, j, k) - p(i, j, k)) + \
coef72*(p(i + 2, j, k) - p(i - 1, j, k)) + \
coef73*(p(i + 3, j, k) - p(i - 2, j, k)) + \
coef74*(p(i + 4, j, k) - p(i - 3, j, k)) + \
coef75*(p(i + 5, j, k) - p(i - 4, j, k)) + \
coef76*(p(i + 6, j, k) - p(i - 5, j, k)) + \
coef77*(p(i + 7, j, k) - p(i - 6, j, k)) \
)

#define pdyp_stencil ( \
coef71*(p(i, j + 1, k) - p(i, j, k)) + \
coef72*(p(i, j + 2, k) - p(i, j - 1, k)) + \
coef73*(p(i, j + 3, k) - p(i, j - 2, k)) + \
coef74*(p(i, j + 4, k) - p(i, j - 3, k)) + \
coef75*(p(i, j + 5, k) - p(i, j - 4, k)) + \
coef76*(p(i, j + 6, k) - p(i, j - 5, k)) + \
coef77*(p(i, j + 7, k) - p(i, j - 6, k)) \
)

#define pdzp_stencil ( \
coef71*(p(i, j, k + 1) - p(i, j, k)) + \
coef72*(p(i, j, k + 2) - p(i, j, k - 1)) + \
coef73*(p(i, j, k + 3) - p(i, j, k - 2)) + \
coef74*(p(i, j, k + 4) - p(i, j, k - 3)) + \
coef75*(p(i, j, k + 5) - p(i, j, k - 4)) + \
coef76*(p(i, j, k + 6) - p(i, j, k - 5)) + \
coef77*(p(i, j, k + 7) - p(i, j, k - 6)) \
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

#define pdyvy_stencil ( \
coef71*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef72*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef73*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef74*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef75*(vy(i, j + 5, k) - vy(i, j - 4, k)) + \
coef76*(vy(i, j + 6, k) - vy(i, j - 5, k)) + \
coef77*(vy(i, j + 7, k) - vy(i, j - 6, k)) \
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

! L=6
#ifdef _fdorder12_

#define pdxp_stencil ( \
coef61*(p(i + 1, j, k) - p(i, j, k)) + \
coef62*(p(i + 2, j, k) - p(i - 1, j, k)) + \
coef63*(p(i + 3, j, k) - p(i - 2, j, k)) + \
coef64*(p(i + 4, j, k) - p(i - 3, j, k)) + \
coef65*(p(i + 5, j, k) - p(i - 4, j, k)) + \
coef66*(p(i + 6, j, k) - p(i - 5, j, k)) \
)

#define pdyp_stencil ( \
coef61*(p(i, j + 1, k) - p(i, j, k)) + \
coef62*(p(i, j + 2, k) - p(i, j - 1, k)) + \
coef63*(p(i, j + 3, k) - p(i, j - 2, k)) + \
coef64*(p(i, j + 4, k) - p(i, j - 3, k)) + \
coef65*(p(i, j + 5, k) - p(i, j - 4, k)) + \
coef66*(p(i, j + 6, k) - p(i, j - 5, k)) \
)

#define pdzp_stencil ( \
coef61*(p(i, j, k + 1) - p(i, j, k)) + \
coef62*(p(i, j, k + 2) - p(i, j, k - 1)) + \
coef63*(p(i, j, k + 3) - p(i, j, k - 2)) + \
coef64*(p(i, j, k + 4) - p(i, j, k - 3)) + \
coef65*(p(i, j, k + 5) - p(i, j, k - 4)) + \
coef66*(p(i, j, k + 6) - p(i, j, k - 5)) \
)

#define pdxvx_stencil ( \
coef61*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef62*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef63*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef64*(vx(i + 4, j, k) - vx(i - 3, j, k)) + \
coef65*(vx(i + 5, j, k) - vx(i - 4, j, k)) + \
coef66*(vx(i + 6, j, k) - vx(i - 5, j, k)) \
)

#define pdyvy_stencil ( \
coef61*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef62*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef63*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef64*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef65*(vy(i, j + 5, k) - vy(i, j - 4, k)) + \
coef66*(vy(i, j + 6, k) - vy(i, j - 5, k)) \
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

! L=5
#ifdef _fdorder10_

#define pdxp_stencil ( \
coef51*(p(i + 1, j, k) - p(i, j, k)) + \
coef52*(p(i + 2, j, k) - p(i - 1, j, k)) + \
coef53*(p(i + 3, j, k) - p(i - 2, j, k)) + \
coef54*(p(i + 4, j, k) - p(i - 3, j, k)) + \
coef55*(p(i + 5, j, k) - p(i - 4, j, k)) \
)

#define pdyp_stencil ( \
coef51*(p(i, j + 1, k) - p(i, j, k)) + \
coef52*(p(i, j + 2, k) - p(i, j - 1, k)) + \
coef53*(p(i, j + 3, k) - p(i, j - 2, k)) + \
coef54*(p(i, j + 4, k) - p(i, j - 3, k)) + \
coef55*(p(i, j + 5, k) - p(i, j - 4, k)) \
)

#define pdzp_stencil ( \
coef51*(p(i, j, k + 1) - p(i, j, k)) + \
coef52*(p(i, j, k + 2) - p(i, j, k - 1)) + \
coef53*(p(i, j, k + 3) - p(i, j, k - 2)) + \
coef54*(p(i, j, k + 4) - p(i, j, k - 3)) + \
coef55*(p(i, j, k + 5) - p(i, j, k - 4)) \
)

#define pdxvx_stencil ( \
coef51*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef52*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef53*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef54*(vx(i + 4, j, k) - vx(i - 3, j, k)) + \
coef55*(vx(i + 5, j, k) - vx(i - 4, j, k)) \
)

#define pdyvy_stencil ( \
coef51*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef52*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef53*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef54*(vy(i, j + 4, k) - vy(i, j - 3, k)) + \
coef55*(vy(i, j + 5, k) - vy(i, j - 4, k)) \
)

#define pdzvz_stencil ( \
coef51*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef52*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef53*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef54*(vz(i, j, k + 4) - vz(i, j, k - 3)) + \
coef55*(vz(i, j, k + 5) - vz(i, j, k - 4)) \
)

#endif

! L=4
#ifdef _fdorder8_

#define pdxp_stencil ( \
coef41*(p(i + 1, j, k) - p(i, j, k)) + \
coef42*(p(i + 2, j, k) - p(i - 1, j, k)) + \
coef43*(p(i + 3, j, k) - p(i - 2, j, k)) + \
coef44*(p(i + 4, j, k) - p(i - 3, j, k)) \
)

#define pdyp_stencil ( \
coef41*(p(i, j + 1, k) - p(i, j, k)) + \
coef42*(p(i, j + 2, k) - p(i, j - 1, k)) + \
coef43*(p(i, j + 3, k) - p(i, j - 2, k)) + \
coef44*(p(i, j + 4, k) - p(i, j - 3, k)) \
)

#define pdzp_stencil ( \
coef41*(p(i, j, k + 1) - p(i, j, k)) + \
coef42*(p(i, j, k + 2) - p(i, j, k - 1)) + \
coef43*(p(i, j, k + 3) - p(i, j, k - 2)) + \
coef44*(p(i, j, k + 4) - p(i, j, k - 3)) \
)

#define pdxvx_stencil ( \
coef41*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef42*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef43*(vx(i + 3, j, k) - vx(i - 2, j, k)) + \
coef44*(vx(i + 4, j, k) - vx(i - 3, j, k)) \
)

#define pdyvy_stencil ( \
coef41*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef42*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef43*(vy(i, j + 3, k) - vy(i, j - 2, k)) + \
coef44*(vy(i, j + 4, k) - vy(i, j - 3, k)) \
)

#define pdzvz_stencil ( \
coef41*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef42*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef43*(vz(i, j, k + 3) - vz(i, j, k - 2)) + \
coef44*(vz(i, j, k + 4) - vz(i, j, k - 3)) \
)

#endif

! L=3
#ifdef _fdorder6_

#define pdxp_stencil ( \
coef31*(p(i + 1, j, k) - p(i, j, k)) + \
coef32*(p(i + 2, j, k) - p(i - 1, j, k)) + \
coef33*(p(i + 3, j, k) - p(i - 2, j, k)) \
)

#define pdyp_stencil ( \
coef31*(p(i, j + 1, k) - p(i, j, k)) + \
coef32*(p(i, j + 2, k) - p(i, j - 1, k)) + \
coef33*(p(i, j + 3, k) - p(i, j - 2, k)) \
)

#define pdzp_stencil ( \
coef31*(p(i, j, k + 1) - p(i, j, k)) + \
coef32*(p(i, j, k + 2) - p(i, j, k - 1)) + \
coef33*(p(i, j, k + 3) - p(i, j, k - 2)) \
)

#define pdxvx_stencil ( \
coef31*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef32*(vx(i + 2, j, k) - vx(i - 1, j, k)) + \
coef33*(vx(i + 3, j, k) - vx(i - 2, j, k)) \
)

#define pdyvy_stencil ( \
coef31*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef32*(vy(i, j + 2, k) - vy(i, j - 1, k)) + \
coef33*(vy(i, j + 3, k) - vy(i, j - 2, k)) \
)

#define pdzvz_stencil ( \
coef31*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef32*(vz(i, j, k + 2) - vz(i, j, k - 1)) + \
coef33*(vz(i, j, k + 3) - vz(i, j, k - 2)) \
)

#endif

! L=2
#ifdef _fdorder4_

#define pdxp_stencil ( \
coef21*(p(i + 1, j, k) - p(i, j, k)) + \
coef22*(p(i + 2, j, k) - p(i - 1, j, k)) \
)

#define pdyp_stencil ( \
coef21*(p(i, j + 1, k) - p(i, j, k)) + \
coef22*(p(i, j + 2, k) - p(i, j - 1, k)) \
)

#define pdzp_stencil ( \
coef21*(p(i, j, k + 1) - p(i, j, k)) + \
coef22*(p(i, j, k + 2) - p(i, j, k - 1)) \
)

#define pdxvx_stencil ( \
coef21*(vx(i + 1, j, k) - vx(i, j, k)) + \
coef22*(vx(i + 2, j, k) - vx(i - 1, j, k)) \
)

#define pdyvy_stencil ( \
coef21*(vy(i, j + 1, k) - vy(i, j, k)) + \
coef22*(vy(i, j + 2, k) - vy(i, j - 1, k)) \
)

#define pdzvz_stencil ( \
coef21*(vz(i, j, k + 1) - vz(i, j, k)) + \
coef22*(vz(i, j, k + 2) - vz(i, j, k - 1)) \
)

#endif

! L=1
#ifdef _fdorder2_

#define pdxp_stencil (p(i + 1, j, k) - p(i, j, k))

#define pdyp_stencil (p(i, j + 1, k) - p(i, j, k))

#define pdzp_stencil (p(i, j, k + 1) - p(i, j, k))

#define pdxvx_stencil (vx(i + 1, j, k) - vx(i, j, k))

#define pdyvy_stencil (vy(i, j + 1, k) - vy(i, j, k))

#define pdzvz_stencil (vz(i, j, k + 1) - vz(i, j, k))

#endif
