import React, { useState, useEffect, useRef } from 'react';
import ContourPlot from './components/ContourPlot';
import { SimParams, Edge, TempBC } from './api';
 
function App() {
  // ——— 新增：最大允许 dt 的状态 ———
  const [maxDt, setMaxDt] = useState<number>(0);
  // ① 四条边的初始 TempBC
  const initialTempBC: Record<Edge, TempBC> = {
    top:    { type: 'Dirichlet', value: 0 },
    bottom: { type: 'Dirichlet', value: 0 },
    left:   { type: 'Dirichlet', value: 0 },
    right:  { type: 'Dirichlet', value: 0 },
  };
  // ② params：实时由表单更新
  const [params, setParams] = useState<SimParams>({
      Nx: 101,
      Ny: 21,
      Lx: 5,
      Ly: 1,
      dt: 0.005, nt: 2000,       // ← 新增
      Re: 20,
      U: 1,
      rho: 1.0,
      method: 'projection',
      tempBC: initialTempBC,  
  });

  // ③ runParams：只有点击 Run 时才同步一次
  const [runParams, setRunParams] = useState<SimParams>(params);
  const [runKey, setRunKey] = useState(0);
  const [field, setField] = useState<'u'|'v'|'p'|'psi'|'omega'|'T'>('u');
  // 四条边的枚举，确保 TS 能识别
  const edges: Edge[] = ['top','bottom','left','right'];

  // 用来记录上一次提示的 dt_max，避免重复弹窗
  const lastPromptedMaxDt = useRef<number | null>(null);

  // 用来显示拒绝自动修改后的警告
  const [manualWarning, setManualWarning] = useState<string | null>(null);

  // 稳定性检查 & 弹窗
  useEffect(() => {
    const { Lx, Ly, Nx, Ny, U, Re, dt } = params;
    const dx = Lx / (Nx - 1);
    const dy = Ly / (Ny - 1);
    const nu = U * Ly / Re;

    const dtConvX  = dx / U;
    const dtConvY  = dy / U;
    const dtVisc   = (dx*dx * dy*dy) / (2 * nu * (dx*dx + dy*dy));
    const dtMax    = Math.min(dtConvX, dtConvY, dtVisc);

    // 更新界面上长期显示的最大 dt
    setMaxDt(dtMax);

    // 只有当 dt > dtMax 且还没对这个 dtMax 弹过窗时，才弹
    if (dt > dtMax && lastPromptedMaxDt.current !== dtMax) {
      lastPromptedMaxDt.current = dtMax;
      // 建议值 = 0.1 * dtMax，并保留两位有效数字
      const suggestedRaw = 0.1 * dtMax;
      const suggested = Number(suggestedRaw.toPrecision(2));
      const msg = 
        `当前 dt = ${dt.toExponential(2)}，超过稳定上限 dt_max = ${dtMax.toExponential(2)}。\n` +
        `建议将 dt 调整为 ${suggested.toExponential(2)}。\n` +
        `是否自动调整？`;
      if (window.confirm(msg)) {
        // 用户接受，自动修改 dt
        setParams(p => ({ ...p, dt: suggested }));
        setManualWarning(null);
      } else {
        // 用户拒绝，保留原 dt，并显示警告
        setManualWarning(
          `请手动将 dt ≤ ${dtMax.toExponential(2)}，否则模拟可能失稳。`
        );
      }
    }
    // 如果 dt <= dtMax，重置提示标志和警告
    else if (dt <= dtMax) {
      lastPromptedMaxDt.current = null;
      setManualWarning(null);
    }
  }, [
    params.Lx, params.Ly,
    params.Nx, params.Ny,
    params.U, params.Re,
    params.dt,
  ]);

  // 2. 简单的表单栏
  return (
    <div style={{ display: 'flex', height: '100vh' }}>
      {/* —— 侧边栏 —— */}
      <aside style={{ width: 300, padding: 16, borderRight: '1px solid #ccc' }}>
        <h3>Simulation Parameters</h3>

        {/* 网格大小 */}
        <label>
          Nx:&nbsp;
          <input
            type="number" 
            min={10}
            value={params.Nx}
            onChange={e => setParams(p => ({ ...p, Nx: +e.target.value }))}
          />
        </label>
        <br/>
        <label>
          Ny:&nbsp;
          <input
            type="number"
            min={10}
            value={params.Ny}
            onChange={e => setParams(p => ({ ...p, Ny: +e.target.value }))}
          />
        </label>
        <br/>

        {/* —— 域尺寸 —— */}
        <label>
          Lx:&nbsp;
          <input
            type="number" step="0.1"
            value={params.Lx}
            onChange={e => setParams(p => ({ ...p, Lx: +e.target.value }))}
          />
        </label>
        <br/>
        <label>
          Ly:&nbsp;
          <input
            type="number" step="0.1"
            value={params.Ly}
            onChange={e => setParams(p => ({ ...p, Ly: +e.target.value }))}
          />
        </label>
        <br/>

        {/* —— 时间参数 —— */}
        <label>
          dt:&nbsp;
          <input
            type="number" step="0.001"
            value={params.dt}
            onChange={e => setParams(p => ({ ...p, dt: +e.target.value }))}
          />
        </label>
        <div>最大允许 dt: <strong>{maxDt.toExponential(2)}</strong></div>
        {manualWarning && (
          <div style={{ color:'red', marginTop:4 }}>
            {manualWarning}
          </div>
        )}
        <br/>
        <label>
          nt:&nbsp;
          <input
            type="number" min={1}
            value={params.nt}
            onChange={e => setParams(p => ({ ...p, nt: +e.target.value }))}
          />
        </label>
        <hr/>

        {/* 数值方法选择 */}
        <label>
          Numerical Method:&nbsp;
          <select
            value={params.method}
            onChange={e => setParams(p => ({
              ...p,
              method: e.target.value as SimParams['method']
            }))}
          >
            <option value="projection">Projection</option>
            <option value="vorticity">Vorticity-Stream</option>
          </select>
        </label>
        <hr/>

        {/* 物理参数 */}
        <label>
          Density rho:&nbsp;
          <input
            type="number" step="0.1"
            value={params.rho}
            onChange={e => setParams(p => ({ ...p, rho: +e.target.value }))}
          />
        </label>
        <br/>
        <label>
          Inlet Speed U:&nbsp;
          <input
            type="number" step="0.1"
            value={params.U}
            onChange={e => setParams(p => ({ ...p, U: +e.target.value }))}
          />
        </label>
        <br/>
        <label>
          Reynolds Re:&nbsp;
          <input
            type="number" step="1"
            value={params.Re}
            onChange={e => setParams(p => ({ ...p, Re: +e.target.value }))}
          />
        </label>
        <hr/>

        {/* 展示场 */}
        <label>
          Simulation Results:&nbsp;
          <select value={field} onChange={e=>setField(e.target.value as any)}>
            {params.method === 'projection' ? (
              <>
                <option value="u">u</option>
                <option value="v">v</option>
                <option value="p">p</option>
              </>
            ) : (
              <>
                <option value="u">u</option>
                <option value="v">v</option>
                <option value="psi">ψ (stream)</option>
                <option value="omega">ω (vorticity)</option>
                <option value="T">T (temperature)</option>
              </>
            )}
          </select>
        </label>
        <hr/>

        {/* 温度场边界：每条边独立设置 */}
        <h4>Temperature BC (per edge)</h4>
        {edges.map(edge => (
          <div key={edge} style={{ marginBottom: 8 }}>
            <strong>{edge}</strong>&nbsp;
            <select
              value={params.tempBC[edge].type}
              onChange={e => {
                const type = e.target.value as TempBC['type'];
                setParams(p => ({
                  ...p,
                  tempBC: {
                    ...p.tempBC,
                    [edge]: { ...p.tempBC[edge], type }
                  }
                }));
              }}
            >
              <option value="Dirichlet">Dirichlet</option>
              <option value="Neumann">Neumann</option>
            </select>
            &nbsp;
            <input
              type="number" step="0.1"
              value={params.tempBC[edge].value}
              onChange={e => {
                const value = Number(e.target.value);
                setParams(p => ({
                  ...p,
                  tempBC: {
                    ...p.tempBC,
                    [edge]: { ...p.tempBC[edge], value }
                  }
                }));
              }}
              style={{ width: 60 }}
            />
          </div>
        ))}
        <hr/>

        {/* Run Simulation 按钮 */}
        <button
          style={{
            marginTop: '16px',
            padding: '8px 16px',
            fontSize: '1rem',
            cursor: 'pointer'
          }}
          onClick={() => {
            // 点击 Run 时，把当前输入 snapshot 给 runParams，再触发图表重载
            setRunParams(params);
            setRunKey(k => k + 1);
          }}
        >
          Run Simulation
        </button>
      </aside>

      {/* 主视图区 */}
      <main style={{ flex: 1, padding: 8 }}>
         <ContourPlot
          key={runKey}
          params={runParams}
          field={field}
          interval={250}
        />
      </main>
    </div>
  );
}

export default App;
