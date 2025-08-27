import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import solve, norm

# =====================================================
# 1. Problem parameters and mesh
# =====================================================
L = 1.0           # Beam length
EI = 1.0          # Flexural rigidity
p_ext = 10.0       # Distributed load (applied in y direction)
r_penalty = 1e2   # Penalty parameter (and used in multiplier update)

N_elem = 20              # Number of finite elements
N_node = N_elem + 1      # Number of nodes
h = L / N_elem           # Element length

# For (x,y): using cubic Hermite elements => 2 dofs per node for x, 2 for y.
# So total degrees for (x,y) is 2*N_node + 2*N_node = 4*N_node.
total_xy = 4 * N_node

# For lambda: linear element => 1 dof per node.
ndof_lambda = N_node

# =====================================================
# 2. Cubic Hermite shape functions on the reference element [0,1]
# =====================================================
def H1(xi):
    return 1 - 3*xi**2 + 2*xi**3

def H2(xi):
    return h * (xi - 2*xi**2 + xi**3)

def H3(xi):
    return 3*xi**2 - 2*xi**3

def H4(xi):
    return h * (xi**3 - xi**2)

def dH1_dxi(xi):
    return -6*xi + 6*xi**2

def dH2_dxi(xi):
    return h * (1 - 4*xi + 3*xi**2)

def dH3_dxi(xi):
    return 6*xi - 6*xi**2

def dH4_dxi(xi):
    return h * (3*xi**2 - 2*xi)

def d2H1_dxi2(xi):
    return -6 + 12*xi

def d2H2_dxi2(xi):
    return h * (-4 + 6*xi)

def d2H3_dxi2(xi):
    return 6 - 12*xi

def d2H4_dxi2(xi):
    return h * (6*xi - 2)

# Map to physical derivatives:
def H_vals(xi):
    return np.array([H1(xi), H2(xi), H3(xi), H4(xi)])

def H_d_phys(xi):
    return np.array([dH1_dxi(xi), dH2_dxi(xi), dH3_dxi(xi), dH4_dxi(xi)]) / h

def H_dd_phys(xi):
    return np.array([d2H1_dxi2(xi), d2H2_dxi2(xi), d2H3_dxi2(xi), d2H4_dxi2(xi)]) / (h**2)

# =====================================================
# 3. Linear shape functions for lambda on [0,1]
# =====================================================
def M_vals(xi):
    return np.array([1 - xi, xi])

# =====================================================
# 4. Quadrature: 3-point Gauss-Legendre on [0,1]
# =====================================================
xi_q = np.array([0.1127016654, 0.5, 0.8872983346])
w_q  = np.array([0.2777777778, 0.4444444444, 0.2777777778])

# =====================================================
# 5. Assembly routines for the (x,y) subproblem
# =====================================================
def assemble_residual_xy(u_xy, lam_arr):
    """
    Assemble the global residual for the (x,y) subproblem.
    u_xy: vector of length total_xy (which is 4*N_node) containing the displacement DOFs.
           The first 2*N_node entries are for x and the next 2*N_node are for y.
    lam_arr: vector of length N_node for lambda.
    Returns R_xy of length total_xy.
    """
    R_xy = np.zeros(total_xy)
    # Loop over elements (each element has two nodes)
    for e in range(N_elem):
        nodes = [e, e+1]
        # Global indices for x (2 dofs per node)
        idx_x = [2*nodes[0], 2*nodes[0] + 1,
                 2*nodes[1], 2*nodes[1] + 1]
        # Global indices for y: offset by 2*N_node
        idx_y = [2*N_node + 2*nodes[0], 2*N_node + 2*nodes[0] + 1,
                 2*N_node + 2*nodes[1], 2*N_node + 2*nodes[1] + 1]
        ux = u_xy[idx_x]
        uy = u_xy[idx_y]
        # For lambda (linear element), indices are just the node numbers
        idx_lam = [e, e+1]
        lam_local = lam_arr[idx_lam]
        
        R_loc_x = np.zeros(4)
        R_loc_y = np.zeros(4)
        
        for xi, wq in zip(xi_q, w_q):
            jac = h  # ds = h dxi
            N = H_vals(xi)
            N_d = H_d_phys(xi)
            N_dd = H_dd_phys(xi)
            M_lin = M_vals(xi)
            
            x_val = np.dot(N, ux)
            xp = np.dot(N_d, ux)
            xpp = np.dot(N_dd, ux)
            
            y_val = np.dot(N, uy)
            yp = np.dot(N_d, uy)
            ypp = np.dot(N_dd, uy)
            
            lam_val = np.dot(M_lin, lam_local)
            
            # Bending energy contributions
            for a in range(4):
                R_loc_x[a] += EI * xpp * N_dd[a] * wq * jac
                R_loc_y[a] += EI * ypp * N_dd[a] * wq * jac
            # External load (only in y-equation)
            for a in range(4):
                R_loc_y[a] -= p_ext * N[a] * wq * jac
            # Constraint contributions
            S_val = xp**2 + yp**2 - 1.0
            for a in range(4):
                R_loc_x[a] += 2 * (lam_val + r_penalty * S_val) * xp * N_d[a] * wq * jac
                R_loc_y[a] += 2 * (lam_val + r_penalty * S_val) * yp * N_d[a] * wq * jac
        
        for i_local, i_global in enumerate(idx_x):
            R_xy[i_global] += R_loc_x[i_local]
        for i_local, i_global in enumerate(idx_y):
            R_xy[i_global] += R_loc_y[i_local]
    return R_xy

def assemble_jacobian_fd_xy(u_xy, lam_arr, eps=1e-6):
    R0 = assemble_residual_xy(u_xy, lam_arr)
    n = len(u_xy)
    J_fd = np.zeros((n, n))
    for i in range(n):
        u_pert = u_xy.copy()
        u_pert[i] += eps
        R1 = assemble_residual_xy(u_pert, lam_arr)
        J_fd[:, i] = (R1 - R0) / eps
    return J_fd

# =====================================================
# 6b. Compute constraint residual S at nodes
# =====================================================
def compute_constraint(u_xy):
    """
    Compute S = x'(s)^2 + y'(s)^2 - 1 at each node.
    For cubic Hermite, the derivative at a node is stored as the second dof.
    Returns an array of length N_node.
    """
    S_node = np.zeros(N_node)
    # x: stored in u_xy[0:2*N_node], y: in u_xy[2*N_node:4*N_node]
    for i in range(N_node):
        x_prime = u_xy[2*i + 1]
        y_prime = u_xy[2*N_node + 2*i + 1]
        S_node[i] = x_prime**2 + y_prime**2 - 1.0
    return S_node

# =====================================================
# 7. Alternating Direction Iteration (ADMM-style)
# =====================================================
# u_xy: vector of length total_xy = 4*N_node
# lam_arr: vector of length N_node

# Initialize u_xy with a straight beam:
u_xy = np.zeros(total_xy)
s_nodes = np.linspace(0, L, N_node)
for i in range(N_node):
    # For x: displacement = s, slope = 1
    u_xy[2*i] = s_nodes[i]
    u_xy[2*i + 1] = 1.0
    # For y: displacement = 0, slope = 0 (stored in second half)
    u_xy[2*N_node + 2*i] = 0.0
    u_xy[2*N_node + 2*i + 1] = 0.0

# Initialize lambda:
lam_arr = np.zeros(N_node)

# ADMM iteration parameters
max_outer = 20      # max outer iterations
tol_outer = 1e-5    # tolerance on constraint residual (L2 norm)
max_inner = 20      # max inner (Newton) iterations for (x,y) subproblem
tol_inner = 1e-6    # tolerance for inner Newton solve

for outer in range(max_outer):
    # Step 1: Solve (x,y) subproblem with fixed lambda via Newton's method
    for inner in range(max_inner):
        R_xy = assemble_residual_xy(u_xy, lam_arr)
        J_xy = assemble_jacobian_fd_xy(u_xy, lam_arr, eps=1e-6)
        # Apply essential BCs at s=0 for u_xy:
        bc_indices = [0, 1, 2*N_node, 2*N_node+1]  # x(0)=0, x'(0)=1, y(0)=0, y'(0)=0
        bc_values = {0: 0.0, 1: 1.0, 2*N_node: 0.0, 2*N_node+1: 0.0}
        for idx in bc_indices:
            R_xy[idx] = u_xy[idx] - bc_values[idx]
            J_xy[idx, :] = 0.0
            J_xy[idx, idx] = 1.0
        norm_inner = norm(R_xy)
        print(f"Outer {outer}, Inner {inner}: |R_xy| = {norm_inner:.3e}")
        if norm_inner < tol_inner:
            break
        delta = solve(J_xy, -R_xy)
        u_xy = u_xy + delta
    # Step 2: Update lambda at nodes
    S_node = compute_constraint(u_xy)
    lam_arr = lam_arr + r_penalty * S_node
    norm_constraint = norm(S_node)
    print(f"Outer {outer}: Constraint norm = {norm_constraint:.3e}")
    if norm_constraint < tol_outer:
        print("Convergence of alternating direction method achieved.")
        break

# =====================================================
# 8. Plot the beam in physical space
# =====================================================
# Extract nodal x and y: for each node, x is stored at u_xy[2*i] and y at u_xy[2*N_node + 2*i]
x_sol = np.array([u_xy[2*i] for i in range(N_node)])
y_sol = np.array([u_xy[2*N_node + 2*i] for i in range(N_node)])
lam_sol = lam_arr.copy()

plt.figure(figsize=(8,5))
plt.plot(x_sol, y_sol, 'bo-', linewidth=2, markersize=6, label='Beam shape')
plt.xlabel('x (physical coordinate)')
plt.ylabel('y (physical coordinate)')
plt.title('Beam Configuration in Physical Space (ADMM)')
plt.legend()
plt.grid(True)
plt.axis('equal')
plt.show()

plt.figure(figsize=(8,4))
plt.plot(s_nodes, lam_sol, 'k.-', label='lambda(s)')
plt.xlabel('s')
plt.ylabel('lambda')
plt.title('Lagrange Multiplier')
plt.legend()
plt.grid(True)
plt.show()

# Plot constraint error at nodes:
S_error = compute_constraint(u_xy)
plt.figure(figsize=(8,4))
plt.plot(s_nodes, S_error, 'r.-', label='Constraint error')
plt.xlabel('s')
plt.ylabel("x'^2+y'^2 - 1")
plt.title('Constraint Residual at Nodes')
plt.legend()
plt.grid(True)
plt.show()
