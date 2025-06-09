import axios from 'axios';

export type Edge = 'top' | 'bottom' | 'left' | 'right';

export interface TempBC {
  type: 'Dirichlet' | 'Neumann';
  value: number;
}

export type Field =
  | 'speed'
  | 'u'
  | 'v'
  | 'p'
  | 'psi'
  | 'omega'
  | 'T';

export interface SimParams {
  Nx: number;
  Ny: number;
  Lx: number;
  Ly: number;
  // 时间步长与总步数
  dt: number;
  nt: number;
  Re: number;
  U: number;
  rho: number;                 // 新增：密度
  // k:    number;  // 热导率 W/(m·K)
  // cp:   number;  // 比热    J/(kg·K)
  method: 'projection' | 'vorticity';
  problemType: 'channel' | 'cavity';
  // —— 新增：温度场边界条件 —— //
  tempBC: Record<Edge, TempBC>;
}

export interface SimResult {
  x:number[]; y:number[];
  u_frames:number[][][];
  v_frames:number[][][];
  p_frames?:number[][][];   // projection 独有
  psi_frames?:number[][][]; // vorticity–stream 独有
  omega_frames?:number[][][];
  T_frames?:number[][][];
  t_frames:number[];        // 新增：每帧对应的模拟时间（秒）
}

export const simulate = async (params: SimParams): Promise<SimResult> => {
  const resp = await axios.post<SimResult>('http://localhost:8000/simulate', params);
  return resp.data;
};
