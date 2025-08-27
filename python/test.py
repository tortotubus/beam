import numpy as np
import matplotlib.pyplot as plt

# Problem parameters
L = 1.0            # Length of the rod
N = 100            # Number of nodes
EI = 1.0           # Bending stiffness
rho = 1.0          # Linear density
g = 9.81           # Gravity
r = 1e3            # ADMM penalty parameter
max_iters = 500    # Number of ADMM iterations
gd_steps = 50      # Gradient descent steps for x,y-subproblem
gd_step_size = 1e-3

# Discretization
s = np.linspace(0, L, N)
h = s[1] - s[0]

# Finite-difference operators for first and fourth derivatives
def Dx(u):
    # Forward difference, assume last derivative equals previous
    du = np.empty_like(u)
    du[:-1] = (u[1:] - u[:-1]) / h
    du[-1] = du[-2]
    return du

def D4(u):
    # Central finite-difference 4th derivative interior, zero at boundaries
    d4 = np.zeros_like(u)
    d4[2:-2] = (u[:-4] - 4*u[1:-3] + 6*u[2:-2] - 4*u[3:-1] + u[4:]) / h**4
    return d4

# Initialization: straight horizontal rod
x = s.copy()
y = np.zeros_like(s)
p = Dx(x)
q = Dx(y)
lam_p = np.zeros_like(s)
lam_q = np.zeros_like(s)

# Clamped boundary function
def enforce_bc(x, y):
    x[0], y[0] = 0.0, 0.0
    x[1], y[1] = h * 1.0, 0.0  # x'(0)=1, y'(0)=0
    x[-1], y[-1] = L, 0.0
    x[-2], y[-2] = L - h * 1.0, 0.0  # x'(L)=1, y'(L)=0

# ADMM iterations
for it in range(max_iters):
    # 1) x,y subproblem: gradient descent on J + penalty terms
    for _ in range(gd_steps):
        # bending gradient
        grad_x = EI * D4(x)
        grad_y = EI * D4(y) - rho * g  # include gravity in y-gradient

        # penalty gradient
        dp = Dx(x) - p + lam_p
        dq = Dx(y) - q + lam_q
        dp_pre = np.zeros_like(dp); dp_pre[1:] = dp[:-1]
        dq_pre = np.zeros_like(dq); dq_pre[1:] = dq[:-1]
        grad_x += -r * (dp - dp_pre) / h
        grad_y += -r * (dq - dq_pre) / h

        # gradient descent
        x -= gd_step_size * grad_x
        y -= gd_step_size * grad_y
        enforce_bc(x, y)

    # 2) p,q subproblem: pointwise projection onto unit circle
    w_p = Dx(x) + lam_p
    w_q = Dx(y) + lam_q
    norm_w = np.sqrt(w_p**2 + w_q**2) + 1e-12
    p = w_p / norm_w
    q = w_q / norm_w

    # 3) dual updates
    lam_p += Dx(x) - p
    lam_q += Dx(y) - q

# Plotting the final configuration
plt.figure(figsize=(6, 4))
plt.plot(x, y, '-o', markevery=10)
plt.title('Static Inextensible Rod via ADMM')
plt.xlabel('x(s)')
plt.ylabel('y(s)')
plt.grid(True)
plt.show()
