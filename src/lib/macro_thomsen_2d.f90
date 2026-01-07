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


#define thomsen_dc11_dvp ( 2*rho*vp*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tithe)**2)*sin(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*cos(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc11_dvs ( rho*vs*(1 - cos(4*tithe))*(-tidel*vp**2 - vp**2 + vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))/(2*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))) )

#define thomsen_dc11_deps ( 2*rho*vp**2*cos(tithe)**4 )

#define thomsen_dc11_ddel ( rho*vp**2*(1 - cos(4*tithe))*(vp**2 - vs**2)/(4*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))) )

#define thomsen_dc13_dvp ( 2.0*rho*vp*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tithe)**2)*cos(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc13_dvs ( 1.0*rho*vs*((-4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2 + 2.0*(-tidel*vp**2 - vp**2 + vs**2 - sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*cos(tithe)**2)*cos(tithe)**2 + (-4.0*sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tithe)**2 + 2.0*(-tidel*vp**2 - vp**2 + vs**2 - sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(tithe)**2)*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc13_deps ( 0.25*rho*vp**2*(1 - cos(4*tithe)) )

#define thomsen_dc13_ddel ( 1.0*rho*vp**2*(vp**2 - vs**2)*(sin(tithe)**4 + cos(tithe)**4)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc15_dvp ( 2.0*rho*vp*(-sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*cos(tithe)**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*sin(tithe)**2 + (-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*sin(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tithe)**2)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc15_dvs ( 0.5*rho*vs*(-tidel*vp**2 - vp**2 + vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))*sin(4*tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc15_deps ( -2.0*rho*vp**2*sin(tithe)*cos(tithe)**3 )

#define thomsen_dc15_ddel ( 0.25*rho*vp**2*(vp**2 - vs**2)*sin(4*tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc33_dvp ( 2.0*rho*vp*((sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*cos(tithe)**2 + (sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*sin(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*cos(tithe)**2)*sin(tithe)**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc33_dvs ( 0.5*rho*vs*(1 - cos(4*tithe))*(-tidel*vp**2 - vp**2 + vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc33_deps ( 2.0*rho*vp**2*sin(tithe)**4 )

#define thomsen_dc33_ddel ( 0.25*rho*vp**2*(1 - cos(4*tithe))*(vp**2 - vs**2)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc35_dvp ( 2.0*rho*vp*(-sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1)*sin(tithe)**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*cos(tithe)**2 + (-tidel*vp**2 - tidel*(vp**2 - vs**2) - vp**2 + vs**2)*cos(tithe)**2 + (tidel*vp**2 + tidel*(vp**2 - vs**2) + vp**2 - vs**2)*sin(tithe)**2)*sin(tithe)*cos(tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc35_dvs ( 0.5*rho*vs*(tidel*vp**2 + vp**2 - vs**2 - sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4))*sin(4*tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc35_deps ( -2.0*rho*vp**2*sin(tithe)**3*cos(tithe) )

#define thomsen_dc35_ddel ( -0.25*rho*vp**2*(vp**2 - vs**2)*sin(4*tithe)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc55_dvp ( 0.25*rho*vp*(1 - cos(4*tithe))*(-2*tidel*vp**2 - 2*tidel*(vp**2 - vs**2) - 2*vp**2 + 2*vs**2 + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2))*(2*tieps + 1) + sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)))/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )

#define thomsen_dc55_dvs ( 1.0*rho*vs*(0.5*tidel*vp**2*(1 - cos(4*tithe)) + 0.5*vp**2*(1 - cos(4*tithe)) - 0.5*vs**2*(1 - cos(4*tithe)) + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*sin(tithe)**4 + 2.0*sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4)*cos(tithe)**4)/sqrt(2*tidel*vp**4 - 2*tidel*vp**2*vs**2 + vp**4 - 2*vp**2*vs**2 + vs**4) )

#define thomsen_dc55_deps ( 0.25*rho*vp**2*(1 - cos(4*tithe)) )

#define thomsen_dc55_ddel ( 0.25*rho*vp**2*(vp**2 - vs**2)*(cos(4*tithe) - 1)/sqrt((vp**2 - vs**2)*(2*tidel*vp**2 + vp**2 - vs**2)) )
