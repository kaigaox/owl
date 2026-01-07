import sympy as sp
from mod_cij import *

#=======================================================================================
# Anisotropy and Cji related functions

"""
	Rule for elasticity tensor <--> Voigt elasticity matrix conversion
"""
pairs = [(0, 0), (1, 1), (2, 2), (1, 2), (0, 2), (0, 1)]

"""
    Rotation that maps the +z axis onto the direction:
        m = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
    Convention: R = Rz(phi) @ Ry(theta)
"""
def R_from_theta_phi(theta, phi):

    cth, sth = sp.cos(theta), sp.sin(theta)
    cph, sph = sp.cos(phi), sp.sin(phi)

    Rz = sp.Matrix([[cph, -sph, 0.0], 
                    [sph, cph, 0.0], 
                    [0.0, 0.0, 1.0]])
    Ry = sp.Matrix([[cth, 0.0, sth], 
                    [0.0, 1.0, 0.0], 
                    [-sth, 0.0, cth]])
    
    return Rz @ Ry

"""
	Convert 6x6 elasticity matrix in Voigt notation to 3x3x3x3 elasticity tensor
"""
def voigt_to_tensor(Cv):

    C = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3)
    
    for I, (i, j) in enumerate(pairs):
        for J, (k, l) in enumerate(pairs):
            val = Cv[I, J]
            C[i, j, k, l] = val
            C[j, i, k, l] = val
            C[i, j, l, k] = val
            C[j, i, l, k] = val
    
    return C

"""
	Convert 3x3x3x3 elasticity tensor to 6x6 elasticity matrix in Voigt notation 
"""
def tensor_to_voigt(C):

    Cv = sp.MutableDenseMatrix.zeros(6, 6)
    
    for I, (i, j) in enumerate(pairs):
        for J, (k, l) in enumerate(pairs):
            Cv[I, J] = C[i, j, k, l]
    
    return sp.Matrix(Cv)


"""
	Rotate elasticity matrix 
"""
def rotate_C_voigt(Cv, th, phi):

    A = R_from_theta_phi(th, phi)
    C = voigt_to_tensor(Cv)
    
    Crot = sp.MutableDenseNDimArray.zeros(3, 3, 3, 3)
    Crot = np.einsum('pi,qj,rk,sl,ijkl->pqrs', A, A, A, A, C, optimize=True)
    
    return tensor_to_voigt(Crot)

"""
	Symbols used for gradient computation
"""
vp, vs, rho = sp.symbols('vp vs rho', positive=True, real=True)
tieps, tidel, tigam, tieta, tithe, tiphi = sp.symbols('tieps tidel tigam tieta tithe tiphi', real=True)

#=======================================================================================
# 2D Thomsen: vp, vs, rho, eps, delta

# VTI stiffnesses (axis = z)
C33 = rho*vp**2
C44 = rho*vs**2
C11 = C33*(1 + 2*tieps)

term = (C33 - C44)
C13 = -C44 + sp.sqrt(term**2 + 2*tidel*C33*term)

C66 = C44 * (1 + 2*tigam)
C12 = C11 - 2*C66

Cv_VTI = sp.Matrix([
    [C11, C12, C13, 0,   0,   0],
    [C12, C11, C13, 0,   0,   0],
    [C13, C13, C33, 0,   0,   0],
    [0,   0,   0,   C44, 0,   0],
    [0,   0,   0,   0,   C44, 0],
    [0,   0,   0,   0,   0,   C66]
])

# TTI stiffness
Cv_TTI = rotate_C_voigt(Cv_VTI, tithe, 0)

sel = []
for i in [0, 2, 4]:
    for j in [0, 2, 4]:
        if j >= i:
            sel.append((i, j))

C2D = sp.Matrix([Cv_TTI[i,j] for (i,j) in sel])

# Jacobian: dC/dm
p = sp.Matrix([vp, vs, tieps, tidel])
J2D = C2D.jacobian(p)

# Print out
param_names = ['vp', 'vs', 'eps', 'del']

for l, ij in enumerate(sel):
    for k, pname in enumerate(param_names):
        
        i = int(ij[0])
        j = int(ij[1])
        
        print('#define thomsen_dc' + str(i + 1) + str(j + 1) + '_d' + pname + ' (', sp.simplify(J2D[l, k]), ')')

        print(' ')


#=======================================================================================
# 2D Alkhalifah-Tsvankin (Oh et al., 2023): vp, vs, rho, eps, eta,  
# where eta = (eps - delta)/(1 + 2*delta)

# # To procedure to get c13 is:

# c11, c13, c33, c44 = sp.symbols('c11 c13 c33 c44', real=True)
# vp, vs, rho = sp.symbols('vp vs rho', positive=True, real=True)
# eps, delta, eta = sp.symbols('eps delta eta', real=True)

# sol = sp.solve(eta - (eps - delta)/(1 + 2*delta), delta)
# delta = sol[0]

# c44 = rho*vs**2
# c33 = rho*vp**2/(1 + 2*eps)
# sol = sp.solve(delta - ((c13 + c44)**2 - (c33 - c44)**2)/(2*c33*(c33 - c44)), c13)

# print('c13 =', sp.simplify(sol[1]))


C11 = rho*vp**2
C44 = rho*vs**2
C33 = C11/(1 + 2*tieps)

C13 = rho*(sp.sqrt((vp**2/(1 + 2*tieps) - vs**2)*(vp**2/(1 + 2*tieta) - vs**2)) - vs**2)

C66 = C44 * (1 + 2*tigam)
C12 = C11 - 2*C66

Cv_VTI = sp.Matrix([
    [C11, C12, C13, 0,   0,   0],
    [C12, C11, C13, 0,   0,   0],
    [C13, C13, C33, 0,   0,   0],
    [0,   0,   0,   C44, 0,   0],
    [0,   0,   0,   0,   C44, 0],
    [0,   0,   0,   0,   0,   C66]
])

# TTI stiffness
Cv_TTI = rotate_C_voigt(Cv_VTI, tithe, 0)

sel = []
for i in [0, 2, 4]:
    for j in [0, 2, 4]:
        if j >= i:
            sel.append((i, j))

C2D = sp.Matrix([Cv_TTI[i,j] for (i,j) in sel])

# Jacobian: ∂C/∂m
p = sp.Matrix([vp, vs, tieps, tieta])
J2D = C2D.jacobian(p)

# Print out
param_names = ['vp', 'vs', 'eps', 'eta']

for l, ij in enumerate(sel):
    for k, pname in enumerate(param_names):
        
        i = int(ij[0])
        j = int(ij[1])
        
        print('#define alkhalifah_tsvankin_dc' + str(i + 1) + str(j + 1) + '_d' + pname + ' (', sp.simplify(J2D[l, k]), ')')

        print(' ')



#=======================================================================================
# 3D Thomsen: vp, vs, rho, eps, delta, gamma


# VTI stiffnesses (axis = z)
C33 = rho*vp**2
C44 = rho*vs**2
C11 = C33*(1 + 2*tieps)

term = (C33 - C44)
C13 = -C44 + sp.sqrt(term**2 + 2*tidel*C33*term)

C66 = C44 * (1 + 2*tigam)
C12 = C11 - 2*C66

Cv_VTI = sp.Matrix([
    [C11, C12, C13, 0,   0,   0],
    [C12, C11, C13, 0,   0,   0],
    [C13, C13, C33, 0,   0,   0],
    [0,   0,   0,   C44, 0,   0],
    [0,   0,   0,   0,   C44, 0],
    [0,   0,   0,   0,   0,   C66]
])

# TTI stiffness
Cv_TTI = rotate_C_voigt(Cv_VTI, tithe, tiphi)

sel = []
for i in range(6):
    for j in range(6):
        if j >= i:
            sel.append((i, j))

C2D = sp.Matrix([Cv_TTI[i,j] for (i,j) in sel])

# Jacobian: dC/dm
p = sp.Matrix([vp, vs, tieps, tidel, tigam])
J2D = C2D.jacobian(p)

# Print out
param_names = ['vp', 'vs', 'eps', 'del', 'gam']

for l, ij in enumerate(sel):
    for k, pname in enumerate(param_names):
        
        i = int(ij[0])
        j = int(ij[1])
        
        print('#define thomsen_dc' + str(i + 1) + str(j + 1) + '_d' + pname + ' (', sp.simplify(J2D[l, k]), ')')

        print(' ')
        
        
#=======================================================================================
# 2D Alkhalifah-Tsvankin (Oh et al., 2023): vp, vs, rho, eps, eta, gamma, 
# where eta = (eps - delta)/(1 + 2*delta)

# VTI stiffnesses (axis = z)
C11 = rho*vp**2
C44 = rho*vs**2
C33 = C11/(1 + 2*tieps)

C13 = rho*(sp.sqrt((vp**2/(1 + 2*tieps) - vs**2)*(vp**2/(1 + 2*tieta) - vs**2)) - vs**2)

C66 = C44 * (1 + 2*tigam)
C12 = C11 - 2*C66

Cv_VTI = sp.Matrix([
    [C11, C12, C13, 0,   0,   0],
    [C12, C11, C13, 0,   0,   0],
    [C13, C13, C33, 0,   0,   0],
    [0,   0,   0,   C44, 0,   0],
    [0,   0,   0,   0,   C44, 0],
    [0,   0,   0,   0,   0,   C66]
])

# TTI stiffness
Cv_TTI = rotate_C_voigt(Cv_VTI, tithe, tiphi)

sel = []
for i in range(6):
    for j in range(6):
        if j >= i:
            sel.append((i, j))

C2D = sp.Matrix([Cv_TTI[i,j] for (i,j) in sel])

# Jacobian: dC/dm
p = sp.Matrix([vp, vs, tieps, tieta, tigam])
J2D = C2D.jacobian(p)

# Print out
param_names = ['vp', 'vs', 'eps', 'eta', 'gam']

for l, ij in enumerate(sel):
    for k, pname in enumerate(param_names):
        
        i = int(ij[0])
        j = int(ij[1])
        
        print('#define alkhalifah_tsvankin_dc' + str(i + 1) + str(j + 1) + '_d' + pname + ' (', sp.simplify(J2D[l, k]), ')')

        print(' ')
