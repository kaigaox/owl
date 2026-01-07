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


#define thomsen_dc11_dvp ( rho*vp*(2.0*(sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(sin(tiphi)**2 + cos(tiphi)**2*cos(tithe)**2) + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2)*sin(tiphi)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2) + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2)*cos(tiphi)**2*cos(tithe)**2 + (2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2*cos(tithe)**2)*sin(tithe)**2*cos(tiphi)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc11_dvs ( 4.0*rho*vs*(-tidel*vp**2*sin(tiphi)**2 - tidel*vp**2*cos(tiphi)**2*cos(tithe)**2 - vp**2*sin(tiphi)**2 - vp**2*cos(tiphi)**2*cos(tithe)**2 + vs**2*sin(tiphi)**2 + vs**2*cos(tiphi)**2*cos(tithe)**2 + sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2 + sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*cos(tiphi)**2*cos(tithe)**2)*sin(tithe)**2*cos(tiphi)**2/sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) )

#define thomsen_dc11_deps ( 2.0*rho*vp**2*(sin(tiphi)**2*sin(tithe)**2 - sin(tithe)**2 + 1)**2 )

#define thomsen_dc11_ddel ( rho*vp**2*(vp**2 - vs**2)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2)*sin(tithe)**2*cos(tiphi)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc11_dgam ( 0 )

#define thomsen_dc12_dvp ( rho*vp*(2.0*(sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(sin(tiphi)**2 + cos(tiphi)**2*cos(tithe)**2) + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2)*cos(tiphi)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2) + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2)*sin(tiphi)**2*cos(tithe)**2 + (2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2*cos(tithe)**2)*sin(tiphi)**2*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc12_dvs ( -rho*vs*(1.0*(4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam*cos(tithe)**2 + 1)*sin(tiphi)**2 + 2.0*(2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam + 1)*cos(tithe)**2 + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2)*cos(tiphi)**2)*cos(tiphi)**2 + (4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + (2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2)*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))))*sin(tiphi)**2*sin(tithe)**2 + (4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam + 1)*sin(tiphi)**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(8.0*tigam + 4.0*cos(tithe)**2)*cos(tiphi)**2 + 2*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2*cos(tiphi)**2)*sin(tiphi)**2*cos(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc12_deps ( 2.0*rho*vp**2*(-sin(tiphi)**4*sin(tithe)**4 + sin(tiphi)**2*sin(tithe)**4 - sin(tithe)**2 + 1) )

#define thomsen_dc12_ddel ( rho*vp**2*(vp**2 - vs**2)*(2.0*sin(tiphi)**4*sin(tithe)**2 - 2.0*sin(tiphi)**2*sin(tithe)**2 + 1.0)*sin(tithe)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc12_dgam ( -4.0*rho*vs**2*cos(tithe)**2 )

#define thomsen_dc13_dvp ( 2.0*rho*vp*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(sin(tiphi)**2 + cos(tiphi)**2*cos(tithe)**2) + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2)*sin(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2*cos(tithe)**2)*cos(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc13_dvs ( -1.0*rho*vs*((4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + 1.0*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2)*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))))*cos(tithe)**2 + (4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam + 1)*sin(tiphi)**2 + 4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tiphi)**2*cos(tithe)**2 + 2.0*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2*cos(tiphi)**2)*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc13_deps ( 2.0*rho*vp**2*(-(1 - cos(tithe)**2)**2*cos(tiphi)**2 - cos(tithe)**2 + 1) )

#define thomsen_dc13_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*((sin(tiphi)**2 + cos(tiphi)**2*cos(tithe)**2)*cos(tithe)**2 + sin(tithe)**4*cos(tiphi)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc13_dgam ( -4.0*rho*vs**2*sin(tiphi)**2*sin(tithe)**2 )

#define thomsen_dc14_dvp ( 1.0*rho*vp*(-sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2) + 2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + 2*(-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*sin(tithe)**2*cos(tiphi)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2*cos(tithe)**2)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc14_dvs ( 1.0*rho*vs*(sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(8.0*tigam + 4.0)*sin(tiphi)**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(8.0*tigam + 4.0*cos(tithe)**2)*cos(tiphi)**2 - 4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + (2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2)*(-tidel*vp**2 - vp**2 + vs**2 - sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))) + 2*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2*cos(tiphi)**2)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc14_deps ( 2.0*rho*vp**2*(-sin(tiphi)**2*sin(tithe)**2 + sin(tithe)**2 - 1)*sin(tiphi)*sin(tithe)*cos(tithe) )

#define thomsen_dc14_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*(2.0*sin(tiphi)**2*sin(tithe)**2 - 2.0*sin(tithe)**2 + 1.0)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc14_dgam ( 1.0*rho*vs**2*(cos(tiphi - 2*tithe) - cos(tiphi + 2*tithe)) )

#define thomsen_dc15_dvp ( 1.0*rho*vp*(-sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2) + 2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + 2*(-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*sin(tithe)**2*cos(tiphi)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2*cos(tithe)**2)*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc15_dvs ( 1.0*rho*vs*(-4.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 + 4.0*tidel*vp**2*sin(tithe)**2 - 2.0*tidel*vp**2 - 4.0*vp**2*sin(tiphi)**2*sin(tithe)**2 + 4.0*vp**2*sin(tithe)**2 - 2.0*vp**2 + 4.0*vs**2*sin(tiphi)**2*sin(tithe)**2 - 4.0*vs**2*sin(tithe)**2 + 2.0*vs**2 + 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 - 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**2 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc15_deps ( 2.0*rho*vp**2*(-sin(tiphi)**2*sin(tithe)**2 + sin(tithe)**2 - 1)*sin(tithe)*cos(tiphi)*cos(tithe) )

#define thomsen_dc15_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*(2.0*sin(tiphi)**2*sin(tithe)**2 - 2.0*sin(tithe)**2 + 1.0)*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc15_dgam ( 0 )

#define thomsen_dc16_dvp ( rho*vp*(-1.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2) + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(2.0*sin(tiphi)**2 + 2*cos(tiphi)**2*cos(tithe)**2) + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2)*cos(tithe)**2 + (2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2*cos(tiphi)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2*cos(tithe)**2)*sin(tithe)**2 + 2.0*(-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*sin(tithe)**2*cos(tiphi)**2)*sin(tiphi)*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc16_dvs ( 1.0*rho*vs*(-4.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 + 4.0*tidel*vp**2*sin(tithe)**2 - 2.0*tidel*vp**2 - 4.0*vp**2*sin(tiphi)**2*sin(tithe)**2 + 4.0*vp**2*sin(tithe)**2 - 2.0*vp**2 + 4.0*vs**2*sin(tiphi)**2*sin(tithe)**2 - 4.0*vs**2*sin(tithe)**2 + 2.0*vs**2 + 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 - 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**2 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) )

#define thomsen_dc16_deps ( 2.0*rho*vp**2*(-sin(tiphi)**2*sin(tithe)**2 + sin(tithe)**2 - 1)*sin(tiphi)*sin(tithe)**2*cos(tiphi) )

#define thomsen_dc16_ddel ( rho*vp**2*(vp**2 - vs**2)*(2.0*sin(tiphi)**2*sin(tithe)**2 - 2.0*sin(tithe)**2 + 1.0)*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc16_dgam ( 0 )

#define thomsen_dc22_dvp ( rho*vp*(2.0*(sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(-sin(tiphi)**2*sin(tithe)**2 + 1) + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*sin(tithe)**2)*cos(tiphi)**2 + (2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(-sin(tiphi)**2*sin(tithe)**2 + 1) + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*sin(tithe)**2)*sin(tiphi)**2*cos(tithe)**2 + (2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tiphi)**2*sin(tithe)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*cos(tithe)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2)*sin(tiphi)**2*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc22_dvs ( 4.0*rho*vs*(tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 - tidel*vp**2 + vp**2*sin(tiphi)**2*sin(tithe)**2 - vp**2 - vs**2*sin(tiphi)**2*sin(tithe)**2 + vs**2 - sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 + sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)**2*sin(tithe)**2/sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) )

#define thomsen_dc22_deps ( 2.0*rho*vp**2*(sin(tiphi)*sin(tithe) - 1)**2*(sin(tiphi)*sin(tithe) + 1)**2 )

#define thomsen_dc22_ddel ( rho*vp**2*(vp**2 - vs**2)*(-2.0*sin(tiphi)**4*sin(tithe)**4 + 2.0*sin(tiphi)**2*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc22_dgam ( 0 )

#define thomsen_dc23_dvp ( 2.0*rho*vp*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(-sin(tiphi)**2*sin(tithe)**2 + 1) + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*sin(tithe)**2)*sin(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tiphi)**2*sin(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2)*cos(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc23_dvs ( -1.0*rho*vs*((4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tiphi)**2*sin(tithe)**2 - 2.0*(sin(tiphi)**2*sin(tithe)**2 - 1)*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))))*cos(tithe)**2 + (4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam + 1)*cos(tiphi)**2 + 4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tiphi)**2*cos(tithe)**2 + 2.0*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tiphi)**2*sin(tithe)**2)*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc23_deps ( 2.0*rho*vp**2*(-sin(tiphi)**2*sin(tithe)**4 + sin(tithe)**2) )

#define thomsen_dc23_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*((-sin(tiphi)**2*sin(tithe)**2 + 1)*cos(tithe)**2 + sin(tiphi)**2*sin(tithe)**4)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc23_dgam ( -4.0*rho*vs**2*sin(tithe)**2*cos(tiphi)**2 )

#define thomsen_dc24_dvp ( 1.0*rho*vp*(-2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(-sin(tiphi)**2*sin(tithe)**2 + 1) + 2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tiphi)**2*sin(tithe)**2 + 2*(-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*sin(tiphi)**2*sin(tithe)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*cos(tithe)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc24_dvs ( 1.0*rho*vs*(4.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*tidel*vp**2 + 4.0*vp**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*vp**2 - 4.0*vs**2*sin(tiphi)**2*sin(tithe)**2 + 2.0*vs**2 - 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc24_deps ( 2.0*rho*vp**2*(sin(tiphi)**2*sin(tithe)**2 - 1)*sin(tiphi)*sin(tithe)*cos(tithe) )

#define thomsen_dc24_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*(-2.0*sin(tiphi)**2*sin(tithe)**2 + 1.0)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc24_dgam ( 0 )

#define thomsen_dc25_dvp ( 1.0*rho*vp*(-2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*(-sin(tiphi)**2*sin(tithe)**2 + 1) + 2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tiphi)**2*sin(tithe)**2 + 2*(-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*sin(tiphi)**2*sin(tithe)**2 + 2*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tiphi)**2*cos(tithe)**2 + 2.0*(tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tiphi)**2)*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc25_dvs ( 1.0*rho*vs*(4.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*tidel*vp**2 + 8.0*tigam*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) + 4.0*vp**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*vp**2 - 4.0*vs**2*sin(tiphi)**2*sin(tithe)**2 + 2.0*vs**2 - 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc25_deps ( 2.0*rho*vp**2*(sin(tiphi)**2*sin(tithe)**2 - 1)*sin(tithe)*cos(tiphi)*cos(tithe) )

#define thomsen_dc25_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*(-2.0*sin(tiphi)**2*sin(tithe)**2 + 1.0)*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc25_dgam ( 1.0*rho*vs**2*(-sin(tiphi - 2*tithe) + sin(tiphi + 2*tithe)) )

#define thomsen_dc26_dvp ( rho*vp*(-8.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 + 4.0*tidel*vp**2 + 4.0*tidel*vs**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*tidel*vs**2 + 4.0*tieps*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 - 4.0*tieps*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) - 4.0*vp**2*sin(tiphi)**2*sin(tithe)**2 + 2.0*vp**2 + 4.0*vs**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*vs**2 + 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 - 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc26_dvs ( 1.0*rho*vs*(4.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*tidel*vp**2 + 4.0*vp**2*sin(tiphi)**2*sin(tithe)**2 - 2.0*vp**2 - 4.0*vs**2*sin(tiphi)**2*sin(tithe)**2 + 2.0*vs**2 - 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**2 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) )

#define thomsen_dc26_deps ( 2.0*rho*vp**2*(sin(tiphi)**2*sin(tithe)**2 - 1)*sin(tiphi)*sin(tithe)**2*cos(tiphi) )

#define thomsen_dc26_ddel ( rho*vp**2*(vp**2 - vs**2)*(-2.0*sin(tiphi)**2*sin(tithe)**2 + 1.0)*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc26_dgam ( 0 )

#define thomsen_dc33_dvp ( 2.0*rho*vp*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*cos(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*sin(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tithe)**2)*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc33_dvs ( 0.5*rho*vs*(1 - cos(4*tithe))*(-tidel*vp**2 - vp**2 + vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc33_deps ( 2.0*rho*vp**2*sin(tithe)**4 )

#define thomsen_dc33_ddel ( 0.25*rho*vp**2*(1 - cos(4*tithe))*(vp**2 - vs**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc33_dgam ( 0 )

#define thomsen_dc34_dvp ( 2.0*rho*vp*(-sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*sin(tithe)**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tithe)**2 + (-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc34_dvs ( 0.25*rho*vs*(cos(tiphi - 4*tithe) - cos(tiphi + 4*tithe))*(tidel*vp**2 + vp**2 - vs**2 - sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc34_deps ( -2.0*rho*vp**2*sin(tiphi)*sin(tithe)**3*cos(tithe) )

#define thomsen_dc34_ddel ( -0.125*rho*vp**2*(vp**2 - vs**2)*(cos(tiphi - 4*tithe) - cos(tiphi + 4*tithe))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc34_dgam ( 0 )

#define thomsen_dc35_dvp ( 2.0*rho*vp*(-sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*sin(tithe)**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tithe)**2 + (-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc35_dvs ( 0.25*rho*vs*(-sin(tiphi - 4*tithe) + sin(tiphi + 4*tithe))*(tidel*vp**2 + vp**2 - vs**2 - sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc35_deps ( -2.0*rho*vp**2*sin(tithe)**3*cos(tiphi)*cos(tithe) )

#define thomsen_dc35_ddel ( 0.125*rho*vp**2*(vp**2 - vs**2)*(sin(tiphi - 4*tithe) - sin(tiphi + 4*tithe))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc35_dgam ( 0 )

#define thomsen_dc36_dvp ( 2.0*rho*vp*(4*tidel*vp**2*sin(tithe)**2 - 2*tidel*vp**2 - 2*tidel*vs**2*sin(tithe)**2 + tidel*vs**2 - 2*tieps*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**2 + 2*vp**2*sin(tithe)**2 - vp**2 - 2*vs**2*sin(tithe)**2 + vs**2 - 2*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**2 + sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc36_dvs ( rho*vs*(-4.0*tidel*vp**2*sin(tithe)**2 + 2.0*tidel*vp**2 + 8.0*tigam*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) - 4.0*vp**2*sin(tithe)**2 + 2.0*vp**2 + 4.0*vs**2*sin(tithe)**2 - 2.0*vs**2 + 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**2 - 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)**2*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc36_deps ( -0.25*rho*vp**2*(cos(2*tithe) - 1)**2*sin(2*tiphi) )

#define thomsen_dc36_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*(2*sin(tithe)**4 - sin(tithe)**2)*sin(tiphi)*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc36_dgam ( 4.0*rho*vs**2*sin(tiphi)*sin(tithe)**2*cos(tiphi) )

#define thomsen_dc44_dvp ( 2.0*rho*vp*(-2*tidel*vp**2 - 2*tidel*(vp**2 - vs**2) - 2*vp**2 + 2*vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1) + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tiphi)**2*sin(tithe)**2*cos(tithe)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc44_dvs ( 2.0*rho*vs*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*((2*tigam + 1)*cos(tiphi)**2 - sin(tiphi)**2*cos(2*tithe)) + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tiphi)**2*cos(tithe)**2)*sin(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(-2*sin(tiphi)**2*sin(tithe)**2 + 1) + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tiphi)**2*sin(tithe)**2)*cos(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc44_deps ( 2.0*rho*vp**2*sin(tiphi)**2*sin(tithe)**2*cos(tithe)**2 )

#define thomsen_dc44_ddel ( -2.0*rho*vp**2*(vp**2 - vs**2)*sin(tiphi)**2*sin(tithe)**2*cos(tithe)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc44_dgam ( 2.0*rho*vs**2*sin(tithe)**2*cos(tiphi)**2 )

#define thomsen_dc45_dvp ( 0.0625*rho*vp*(2*sin(2*tiphi) - sin(2*tiphi - 4*tithe) - sin(2*tiphi + 4*tithe))*(-2*tidel*vp**2 - 2*tidel*(vp**2 - vs**2) - 2*vp**2 + 2*vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1) + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc45_dvs ( 2.0*rho*vs*((1 - cos(4*tithe))*(tidel*vp**2 + vp**2 - vs**2 - sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))/8 + (-2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(tigam + cos(tithe)**2) + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*cos(tithe)**2)*sin(tithe)**2)*sin(tiphi)*cos(tiphi)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc45_deps ( 0.0625*rho*vp**2*(2*sin(2*tiphi) - sin(2*tiphi - 4*tithe) - sin(2*tiphi + 4*tithe)) )

#define thomsen_dc45_ddel ( 0.0625*rho*vp**2*(vp**2 - vs**2)*(-2*sin(2*tiphi) + sin(2*tiphi - 4*tithe) + sin(2*tiphi + 4*tithe))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc45_dgam ( -2.0*rho*vs**2*sin(tiphi)*sin(tithe)**2*cos(tiphi) )

#define thomsen_dc46_dvp ( 4.0*rho*vp*(-2*tidel*vp**2 + tidel*vs**2 + tieps*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) - vp**2 + vs**2 + sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)**2*sin(tithe)**3*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc46_dvs ( rho*vs*(-4.0*tigam*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tiphi)**2 - 2.0*(-2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(tigam + cos(tithe)**2) + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*cos(tithe)**2)*sin(tiphi)**2 + 2.0*(tidel*vp**2 + vp**2 - vs**2 - sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)**2*sin(tithe)**2 + (2.0*tidel*vp**2 + 2.0*vp**2 - 2.0*vs**2 - 4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam + 1) + 2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tiphi)**2)*sin(tithe)*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc46_deps ( 2.0*rho*vp**2*sin(tiphi)**2*sin(tithe)**3*cos(tiphi)*cos(tithe) )

#define thomsen_dc46_ddel ( -2.0*rho*vp**2*(vp**2 - vs**2)*sin(tiphi)**2*sin(tithe)**3*cos(tiphi)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc46_dgam ( 0.5*rho*vs**2*(sin(tiphi - 2*tithe) - sin(tiphi + 2*tithe)) )

#define thomsen_dc55_dvp ( 2.0*rho*vp*(-2*tidel*vp**2 - 2*tidel*(vp**2 - vs**2) - 2*vp**2 + 2*vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1) + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2*cos(tiphi)**2*cos(tithe)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc55_dvs ( 2.0*rho*vs*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*((2*tigam + 1)*sin(tiphi)**2 - cos(tiphi)**2*cos(2*tithe)) + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*cos(tiphi)**2*cos(tithe)**2)*sin(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(sin(tiphi)**2 + cos(tiphi)**2*cos(2*tithe)) + (tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2*cos(tiphi)**2)*cos(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc55_deps ( 2.0*rho*vp**2*sin(tithe)**2*cos(tiphi)**2*cos(tithe)**2 )

#define thomsen_dc55_ddel ( -2.0*rho*vp**2*(vp**2 - vs**2)*sin(tithe)**2*cos(tiphi)**2*cos(tithe)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc55_dgam ( 2.0*rho*vs**2*sin(tiphi)**2*sin(tithe)**2 )

#define thomsen_dc56_dvp ( 4.0*rho*vp*(-2*tidel*vp**2 + tidel*vs**2 + tieps*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) - vp**2 + vs**2 + sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(tiphi)*sin(tithe)**3*cos(tiphi)**2*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc56_dvs ( rho*vs*(-2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*((2*tigam + 1)*sin(tiphi)**2 - cos(tiphi)**2*cos(2*tithe)) + 2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(sin(tiphi)**2 + cos(tiphi)**2*cos(2*tithe)) + 2.0*(-tidel*vp**2 - vp**2 + vs**2 - sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*cos(tiphi)**2*cos(tithe)**2 + 2.0*(tidel*vp**2 + vp**2 - vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2*cos(tiphi)**2 + (2.0*tidel*vp**2 + 4.0*tigam*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) + 2.0*vp**2 - 2.0*vs**2 - 4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tigam + 1) + 2.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*cos(tiphi)**2)*sin(tiphi)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc56_deps ( 2.0*rho*vp**2*sin(tiphi)*sin(tithe)**3*cos(tiphi)**2*cos(tithe) )

#define thomsen_dc56_ddel ( -2.0*rho*vp**2*(vp**2 - vs**2)*sin(tiphi)*sin(tithe)**3*cos(tiphi)**2*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc56_dgam ( -0.5*rho*vs**2*(cos(tiphi - 2*tithe) - cos(tiphi + 2*tithe)) )

#define thomsen_dc66_dvp ( rho*vp*(1 - cos(4*tiphi))*(1 - cos(2*tithe))**2*(-8.0*tidel*vp**2 + 4.0*tidel*vs**2 + 4.0*tieps*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) - 4.0*vp**2 + 4.0*vs**2 + 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))/(32*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))) )

#define thomsen_dc66_dvs ( rho*vs*(-4.0*tidel*vp**2*sin(tiphi)**4*sin(tithe)**4 + 4.0*tidel*vp**2*sin(tiphi)**2*sin(tithe)**4 - 4.0*tigam*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**2 + 4.0*tigam*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) - 4.0*vp**2*sin(tiphi)**4*sin(tithe)**4 + 4.0*vp**2*sin(tiphi)**2*sin(tithe)**4 + 4.0*vs**2*sin(tiphi)**4*sin(tithe)**4 - 4.0*vs**2*sin(tiphi)**2*sin(tithe)**4 + 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**4*sin(tithe)**4 - 4.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tiphi)**2*sin(tithe)**4 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc66_deps ( -0.0625*rho*vp**2*(cos(4*tiphi) - 1)*(cos(2*tithe) - 1)**2 )

#define thomsen_dc66_ddel ( 0.0625*rho*vp**2*(vp**2 - vs**2)*(cos(4*tiphi) - 1)*(cos(2*tithe) - 1)**2/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc66_dgam ( 2.0*rho*vs**2*cos(tithe)**2 )
