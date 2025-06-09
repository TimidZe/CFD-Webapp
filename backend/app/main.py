# backend/app/main.py
from fastapi import FastAPI
from pydantic import BaseModel
from typing import List, Literal, Dict, Optional
import numpy as np
from .solver import run_projection, run_vorticity_stream  # 把数值求解方法封装成函数
from fastapi.middleware.cors import CORSMiddleware
from fastapi.encoders import jsonable_encoder
from fastapi.responses import JSONResponse
import math

class TempBC(BaseModel):
    type: Literal['Dirichlet', 'Neumann']
    value: float

class SimParams(BaseModel):
    Nx: int
    Ny: int
    Lx: float
    Ly: float
    dt: float         # ← 新增
    nt: int           # ← 新增
    Re: float
    U: float
    rho: float
    method: Literal['projection','vorticity']
    problemType: Literal['channel','cavity']
    tempBC: Dict[str, TempBC]    # keys: 'top','bottom','left','right'

class SimResult(BaseModel):
    x: List[float]
    y: List[float]
    u_frames: List[List[List[float]]]
    v_frames: List[List[List[float]]]
    # projection 独有：
    p_frames: List[List[List[float]]]
    # vorticity–stream function 独有：
    psi_frames:   Optional[list[list[list[float]]]] = None
    omega_frames: Optional[list[list[list[float]]]] = None
    T_frames:     Optional[list[list[list[float]]]] = None
    t_frames: List[float]

app = FastAPI()
# 在现有 app = FastAPI() 之后、所有路由注册之前，添加：
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],            # 或者只写你的前端地址 ["http://localhost:3000"]
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

def sanitize(obj):
    """
    递归地将所有非有限 float (NaN, inf) 替换为 None，其它值原样返回。
    """
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, list):
        return [sanitize(v) for v in obj]
    if isinstance(obj, dict):
        return {k: sanitize(v) for k, v in obj.items()}
    return obj

# @app.post("/simulate", response_model=SimResult)
@app.post("/simulate")
def simulate(params: SimParams):
    if params.method == 'projection':
        x, y, u_hist, v_hist, p_hist, t_hist = run_projection(
            params.Lx,
            params.Ly,
            params.U,
            params.Re,
            params.Nx,
            params.Ny,
            params.dt, 
            params.nt,
            params.rho,
            params.problemType
        )
        result = {
            "x": x.tolist(),
            "y": y.tolist(),
            "u_frames": [u.tolist() for u in u_hist],
            "v_frames": [v.tolist() for v in v_hist],
            "p_frames": [p.tolist() for p in p_hist],
            "psi_frames": None,
            "omega_frames": None,
            "T_frames": None,
            "t_frames": t_hist,
        }
        # return {
        #     "x": x.tolist(), "y": y.tolist(),
        #     "u_frames": [u.tolist() for u in u_hist],
        #     "v_frames": [v.tolist() for v in v_hist],
        #     "p_frames": [p.tolist() for p in p_hist],
        #     "t_frames": t_hist
        # }
    else:  # vorticity–stream function
        # alpha = params.k / (params.rho * params.cp)
        x, y, u_hist, v_hist, psi_hist, omega_hist, T_hist, t_hist = run_vorticity_stream(
            params.Lx,
            params.Ly,
            params.U,
            params.Re,
            params.Nx,
            params.Ny,
            params.dt, params.nt,
            params.rho,
            params.tempBC,
            params.problemType
        )
        result = {
            "x": x.tolist(),
            "y": y.tolist(),
            "u_frames": [u.tolist() for u in u_hist],
            "v_frames": [v.tolist() for v in v_hist],
            "p_frames": [],  # projection-only
            "psi_frames": [s.tolist() for s in psi_hist],
            "omega_frames": [o.tolist() for o in omega_hist],
            "T_frames": [T.tolist() for T in T_hist],
            "t_frames": t_hist,
        }
        # return {
        #     "x": x.tolist(), "y": y.tolist(),
        #     "u_frames": [u.tolist() for u in u_hist],
        #     "v_frames": [v.tolist() for v in v_hist],
        #     "psi_frames": [s.tolist() for s in psi_hist],
        #     "omega_frames": [o.tolist() for o in omega_hist],
        #     "T_frames": [T.tolist() for T in T_hist],
        #     "p_frames": [],                           # ← 占位，确保字段存在
        #     "t_frames": t_hist
        # }
    # 递归清理 NaN/Inf，再返回
    safe = sanitize(result)
    return JSONResponse(content=safe)
    




