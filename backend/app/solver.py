# backend/app/solver.py

import numpy as np
from typing import Dict, Literal
from pydantic import BaseModel

# 在这里定义 TempBC
class TempBC(BaseModel):
    type: Literal['Dirichlet','Neumann']
    value: float

# class SimParams(BaseModel):
#     Lx:       float
#     Ly:       float
#     U:        float
#     Re:       float
#     Nx:       int
#     Ny:       int
#     dt:       float         # ← 新增
#     nt:       int           # ← 新增
#     rho:      float
#     alpha: float
#     method:   Literal['projection','vorticity']
#     problemType: Literal['channel','cavity']
#     tempBC:   Dict[str, TempBC]   # keys: 'top','bottom','left','right'

# def apply_bcs(u: np.ndarray, v: np.ndarray, U: float) -> None:
#     """
#     Apply no-slip and inlet/outlet boundary conditions:
#       - bottom (y=0): u=v=0
#       - left (x=0):  u=U, v=0 (driven lid)
#       - right (x=Lx): du/dx=0, v=0
#       - top (y=Ly):  u=v=0
#     """
#     # bottom
#     u[0, :] = 0.0
#     v[0, :] = 0.0
#     # left
#     u[:, 0] = U
#     v[:, 0] = 0.0
#     # right (zero-gradient for u, v=0)
#     u[:, -1] = u[:, -2]
#     v[:, -1] = 0.0
#     # top
#     u[-1, :] = 0.0
#     v[-1, :] = 0.0

def apply_bcs_channel(u, v, U):
    # 通道流：左入口 u=U，其余壁面无滑
    u[0 ,:] = 0;    v[0 ,:] = 0    # bottom
    u[: ,0] = U;    v[: ,0] = 0    # left inlet
    u[:,-1] = u[:,-2]; v[:,-1] = 0 # right zero‐grad
    u[-1,:] = 0;    v[-1,:] = 0    # top

def apply_bcs_cavity(u, v, U):
    # 驱动腔：top 壁驱动 u=U，其余三边无滑
    # u[0 ,:] = 0;    v[0 ,:] = 0   # bottom
    # u[: ,0] = 0;    v[: ,0] = 0   # left
    # u[:,-1] = 0;   v[:,-1] = 0    # right
    u[-1,:] = U;   v[-1,:] = 0    # top lid

def pressure_poisson(p: np.ndarray, dx: float, dy: float, b: np.ndarray,problemType: str,  nit: int = 50) -> np.ndarray:
    """
    Solve the Pressure Poisson equation:
      ∇²p = b
    with Neumann/Dirichlet BCs applied each iteration.
    """
    pn = np.empty_like(p)
    for _ in range(nit):
        pn[:] = p
        p[1:-1,1:-1] = (
            ((pn[1:-1,2:] + pn[1:-1,0:-2]) * dy**2 +
             (pn[2:,1:-1] + pn[0:-2,1:-1]) * dx**2 -
             b[1:-1,1:-1] * dx**2 * dy**2)
            / (2 * (dx**2 + dy**2))
        )
        # Boundary conditions
        if problemType=='channel':
            p[:, -1]  = 0.0          # right
            p[:,  0]  = p[:,  1]     # left
            p[ 0, :]  = p[ 1, :]     # bottom
            p[-1, :]  = p[-2, :]     # top
        else:
            p[:, -1]  = p[:, -2]     # right
            p[:,  0]  = p[:,  1]     # left
            p[ 0, :]  = p[ 1, :]     # bottom
            p[-1, :]  = 0.0          # top
    return p

def update_velocity_upwind(
    u: np.ndarray,
    v: np.ndarray,
    p: np.ndarray,
    un: np.ndarray,
    vn: np.ndarray,
    dx: float,
    dy: float,
    dt: float,
    rho: float,
    nu: float
) -> None:
    """
    Update the interior of u and v in-place using:
      - 1st-order upwind for the convective terms
      - central differences for pressure gradient and diffusion
    
    Parameters
    ----------
    u, v : 2D arrays of shape (Ny, Nx)
        Velocity fields to be updated.
    p : 2D array of shape (Ny, Nx)
        Pressure field at current time step.
    un, vn : 2D arrays of shape (Ny, Nx)
        u and v at previous time step.
    dx, dy : float
        Grid spacing in x and y.
    dt : float
        Time step.
    rho : float
        Density.
    nu : float
        Kinematic viscosity.
    """
    # alias center and neighbor slices
    uc = un[1:-1,1:-1]; ue = un[1:-1,2:]; uw = un[1:-1,0:-2]
    un_ = un[2:,1:-1]; us = un[0:-2,1:-1]
    vc = vn[1:-1,1:-1]; ve = vn[1:-1,2:]; vw = vn[1:-1,0:-2]
    vn_ = vn[2:,1:-1]; vs = vn[0:-2,1:-1]

    # upwind convective derivatives
    du_dx = np.where(uc > 0, (uc - uw)/dx, (ue - uc)/dx)
    du_dy = np.where(vc > 0, (uc - us)/dy, (un_ - uc)/dy)
    dv_dx = np.where(uc > 0, (vc - vw)/dx, (ve - vc)/dx)
    dv_dy = np.where(vc > 0, (vc - vs)/dy, (vn_ - vc)/dy)

    # pressure gradients (central)
    dp_dx = (p[1:-1,2:] - p[1:-1,0:-2]) / (2*rho*dx)
    dp_dy = (p[2:,1:-1] - p[0:-2,1:-1]) / (2*rho*dy)

    # diffusion terms (central)
    diff_u = nu * (
        (un[1:-1,2:] - 2*uc + un[1:-1,0:-2]) / dx**2 +
        (un[2:,1:-1] - 2*uc + un[0:-2,1:-1]) / dy**2
    )
    diff_v = nu * (
        (vn[1:-1,2:] - 2*vc + vn[1:-1,0:-2]) / dx**2 +
        (vn[2:,1:-1] - 2*vc + vn[0:-2,1:-1]) / dy**2
    )

    # update interior points in-place
    u[1:-1,1:-1] = uc \
        - dt * (uc * du_dx + vc * du_dy) \
        - dt * dp_dx \
        + dt * diff_u

    v[1:-1,1:-1] = vc \
        - dt * (uc * dv_dx + vc * dv_dy) \
        - dt * dp_dy \
        + dt * diff_v
#----------------------------------以上是需要的辅助函数--------------------------------
#---------------------------------- 以下是数值求解方法  -------------------------------

def run_projection(
    Lx: float,
    Ly: float,
    U:  float,
    Re: float,
    Nx: int,
    Ny: int,
    dt: float, nt: int,
    rho: float,
    problemType: str
):
    bc_func = apply_bcs_channel if problemType=='channel' else apply_bcs_cavity
    """
    Run the 2D incompressible projection-method solver.

    Parameters
    ----------
    Lx, Ly : float
        Domain dimensions in x and y.
    U : float
        Inlet velocity (left boundary).
    Re : float
        Reynolds number.
    Nx, Ny : int
        Number of grid points in x and y.

    Returns
    -------
    x, y : 1D np.ndarray
        Grid coordinates.
    u_hist, v_hist, p_hist : list of 2D np.ndarray
        Time-histories of the u, v velocity fields and pressure field,
        sampled at regular intervals.
    """
    # Physical parameters
    # rho = 1.0   # TODO: 1. 把rho设为可以从网页输入的参数
    nu  = U * Ly / Re

    # Grid spacing
    dx = Lx / (Nx - 1)
    dy = Ly / (Ny - 1)

    # Time-stepping parameters
    # dt = 0.005
    # nt = 2000

    # How often to save a snapshot (20 frames total)
    skip = max(1, nt // 50)

    # Initialize fields
    u = np.zeros((Ny, Nx))
    v = np.zeros((Ny, Nx))
    p = np.zeros((Ny, Nx))
    b = np.zeros((Ny, Nx))

    # Storage for snapshots
    u_hist = []
    v_hist = []
    p_hist = []
    t_hist = []

    # Prepare mesh coordinates
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)

    # Time-stepping loop
    for n in range(nt):
        # 1) Apply BCs to current velocity
        bc_func(u, v, U)

        # 2) Copy for calculations
        un = u.copy()
        vn = v.copy()

        # 3) Build source term for pressure Poisson
        b[1:-1,1:-1] = (
            rho * (
                ( (un[1:-1,2:] - un[1:-1,0:-2]) / (2*dx) +
                  (vn[2:,1:-1] - vn[0:-2,1:-1]) / (2*dy)
                ) / dt
                - ((un[1:-1,2:] - un[1:-1,0:-2]) / (2*dx))**2
                - 2 * ( (un[2:,1:-1] - un[0:-2,1:-1]) / (2*dy) *
                        (vn[1:-1,2:] - vn[1:-1,0:-2]) / (2*dx)
                      )
                - ((vn[2:,1:-1] - vn[0:-2,1:-1]) / (2*dy))**2
            )
        )

        # 4) Solve for pressure
        p = pressure_poisson(p, dx, dy, b, problemType)

        # 5) Update velocities
        if problemType=='channel':
            u[1:-1,1:-1] = (
                un[1:-1,1:-1]
                - un[1:-1,1:-1] * dt/dx * (un[1:-1,1:-1] - un[1:-1,0:-2])
                - vn[1:-1,1:-1] * dt/dy * (un[1:-1,1:-1] - un[0:-2,1:-1])
                - dt/(2*rho*dx) * (p[1:-1,2:] - p[1:-1,0:-2])
                + nu * (
                    dt/dx**2 * (un[1:-1,2:] - 2*un[1:-1,1:-1] + un[1:-1,0:-2]) +
                    dt/dy**2 * (un[2:,1:-1] - 2*un[1:-1,1:-1] + un[0:-2,1:-1])
                )
            )
            v[1:-1,1:-1] = (
                vn[1:-1,1:-1]
                - un[1:-1,1:-1] * dt/dx * (vn[1:-1,1:-1] - vn[1:-1,0:-2])
                - vn[1:-1,1:-1] * dt/dy * (vn[1:-1,1:-1] - vn[0:-2,1:-1])
                - dt/(2*rho*dy) * (p[2:,1:-1] - p[0:-2,1:-1])
                + nu * (
                    dt/dx**2 * (vn[1:-1,2:] - 2*vn[1:-1,1:-1] + vn[1:-1,0:-2]) +
                    dt/dy**2 * (vn[2:,1:-1] - 2*vn[1:-1,1:-1] + vn[0:-2,1:-1])
                )
            )
        else:
            update_velocity_upwind(u, v, p, un, vn, dx, dy, dt, rho, nu)

        # 6) Re-apply BCs after update
        bc_func(u, v, U)

        # 7) Save snapshots
        if n % skip == 0:
            u_hist.append(u.copy())
            v_hist.append(v.copy())
            p_hist.append(p.copy())
            t_hist.append(n * dt)

    return x, y, u_hist, v_hist, p_hist, t_hist

def run_vorticity_stream(
    Lx: float, Ly: float, U: float, Re: float,
    Nx: int, Ny: int, dt: float, nt: int, rho: float, 
    tempBC: Dict[str, TempBC],problemType: str      # 直接从 SimParams 拿过来
):
    bc_func = apply_bcs_channel if problemType=='channel' else apply_bcs_cavity
    """
    2D vorticity–stream function 方法。
    返回 x, y, u_hist, v_hist, psi_hist, omega_hist, T_hist, t_hist
    """
    # 1. 初始化参数
    nu = U * Ly / Re
    alpha = nu  # 热扩散系数
    dx, dy = Lx/(Nx-1), Ly/(Ny-1)
    # dt, nt = 0.01, 1000
    skip = max(1, nt // 100)

    # 2. 字段初始化
    omega = np.zeros((Ny, Nx))
    psi   = np.zeros((Ny, Nx))
    T     = np.zeros((Ny, Nx)) + 20  # 初始 T

    # 3. 网格
    x = np.linspace(0, Lx, Nx)
    y = np.linspace(0, Ly, Ny)
    X, Y = np.meshgrid(x, y)

    # 4. 初始流函数：带入口驱动
    if problemType=='channel':
        psi[:, :] = U * Y

    # 5. 历史数据容器
    psi_hist, omega_hist = [], []
    u_hist, v_hist = [], []
    T_hist, t_hist = [], []

    # Precompute constants for SOR
    coefdx = dx**2 / (2.0*(dx**2 + dy**2))
    coefdy = dy**2 / (2.0*(dx**2 + dy**2))
    coefdxdy = dx**2 * dy**2 / (2.0*(dx**2 + dy**2))

    # 6. 计算主循环
    for n in range(nt):
        # 边界条件 for ω
        omega[:, 0]   = (psi[:, 0] - psi[:, 1])*2.0/dx**2
        omega[0, :]   = (psi[0, :] - psi[1, :])*2.0/dy**2
        if problemType=='channel':
            omega[-1, :]  = (psi[-1, :] - psi[-2, :])*2.0/dy**2
            omega[:, -1]  = omega[:, -2]
        else:
            omega[-1,:] = (psi[-1,:] - psi[-2,:])*2.0/dy**2 - 2.0*U/dy
            omega[:, -1]  = (psi[:, -1] - psi[:, -2]) * 2.0 / dx**2

        # --- Step 2: Advect and diffuse ω (explicit) ---
        px = -(psi[2:,1:-1] - psi[:-2,1:-1]) * (omega[1:-1,2:] - omega[1:-1,:-2]) / (4*dy*dx)
        py =  (omega[2:,1:-1] - omega[:-2,1:-1]) * (psi[1:-1,2:] - psi[1:-1,:-2]) / (4*dy*dx)
        pxy = (omega[1:-1,2:] + omega[1:-1,:-2] - 2*omega[1:-1,1:-1]) / dx**2 \
            + (omega[2:,1:-1] + omega[:-2,1:-1] - 2*omega[1:-1,1:-1]) / dy**2
        omega[1:-1,1:-1] += dt * (px + py + nu * pxy)

        # --- 温度场边界条件设置 ---

        # Top boundary (y = Ly, last row)
        if tempBC['top'].type == 'Dirichlet':
            T[-1, :] = tempBC['top'].value
        else:
            # Neumann: ∂T/∂y = value → T[-1,:] = T[-2,:] + value*dy
            T[-1, :] = T[-2, :] + tempBC['top'].value * dy

        # Bottom boundary (y = 0, first row)
        if tempBC['bottom'].type == 'Dirichlet':
            T[0, :] = tempBC['bottom'].value
        else:
            # Neumann: ∂T/∂y = value → T[0,:] = T[1,:] - value*dy
            T[0, :] = T[1, :] - tempBC['bottom'].value * dy

        # Left boundary (x = 0, first column)
        if tempBC['left'].type == 'Dirichlet':
            T[:, 0] = tempBC['left'].value
        else:
            # Neumann: ∂T/∂x = value → T[:,0] = T[:,1] - value*dx
            T[:, 0] = T[:, 1] - tempBC['left'].value * dx

        # Right boundary (x = Lx, last column)
        if tempBC['right'].type == 'Dirichlet':
            T[:, -1] = tempBC['right'].value
        else:
            # Neumann: ∂T/∂x = value → T[:,-1] = T[:,-2] + value*dx
            T[:, -1] = T[:, -2] + tempBC['right'].value * dx

        px_T = -(psi[2:,1:-1] - psi[:-2,1:-1]) * (T[1:-1,2:] - T[1:-1,:-2]) / (4*dy*dx)
        py_T =  (T[2:,1:-1] - T[:-2,1:-1]) * (psi[1:-1,2:] - psi[1:-1,:-2]) / (4*dy*dx)
        pxy_T = (T[1:-1,2:] + T[1:-1,:-2] - 2*T[1:-1,1:-1]) / dx**2 \
            + (T[2:,1:-1] + T[:-2,1:-1] - 2*T[1:-1,1:-1]) / dy**2
        T[1:-1,1:-1] += dt * (px_T + py_T + alpha * pxy_T)

        # --- Step 3: Solve ∇²ψ = -ω using SOR ---
        err = 1.0
        counter = 0
        while err > 1e-5 and counter < 200:
            psi_old = psi.copy()
            psi[1:-1,1:-1] = coefdxdy * omega[1:-1,1:-1] \
                           + coefdx * (psi[2:,1:-1] + psi[:-2,1:-1]) \
                           + coefdy * (psi[1:-1,2:] + psi[1:-1,:-2])
            # psi[-1, :] = psi[-2, :]#??????-------------
            if problemType=='channel':
                psi[:, -1] = psi[:, -2]
            # psi[:, 0] = U * y              # Inlet(Neumann)
            # psi[0, :] = 0                  # Bottom wall
            # psi[-1, :] = U * Ly            # Top wall
            err = np.linalg.norm(psi - psi_old, ord=np.inf)
            counter += 1

        # 保存快照
        if n % skip == 0:
            # 从 ψ 计算 u,v
            u = np.zeros_like(psi)
            v = np.zeros_like(psi)
            u[1:-1,1:-1] =  (psi[2:,1:-1] - psi[:-2,1:-1])/(2*dy)
            v[1:-1,1:-1] = -(psi[1:-1,2:] - psi[1:-1,:-2])/(2*dx)
            bc_func(u, v, U)

            psi_hist.append(psi.copy())
            omega_hist.append(omega.copy())
            u_hist.append(u.copy())
            v_hist.append(v.copy())
            T_hist.append(T.copy())
            t_hist.append(n * dt)

    return x, y, u_hist, v_hist, psi_hist, omega_hist, T_hist, t_hist