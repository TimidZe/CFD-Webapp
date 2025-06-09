import React, { useState, useEffect, useRef } from 'react';
import ContourPlot from './components/ContourPlot';
import ProfilePlot from './components/ProfilePlot';  // 新增
import { SimParams, Edge, TempBC, SimResult, Field } from './api';
import { simulate } from './api';

function App() {
  //  —— 最大允许 dt —— 
  const [maxDt, setMaxDt] = useState<number>(0);

  // —— 四条边的初始 TempBC —— 
  const initialTempBC: Record<Edge, TempBC> = {
    top:    { type: 'Dirichlet', value: 0 },
    bottom: { type: 'Dirichlet', value: 0 },
    left:   { type: 'Dirichlet', value: 0 },
    right:  { type: 'Dirichlet', value: 0 },
  };

  // —— 实时表单输入 params —— 
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
      problemType: 'channel',
      tempBC: initialTempBC,  
  });

  // —— 点击 Run 时真正运行的参数 runParams + key ——
  const [runParams, setRunParams] = useState<SimParams>(params);
  const [runKey, setRunKey] = useState(0);

  // —— 选择展示的场 field —— 
  // const [field, setField] = useState<'u'|'v'|'p'|'psi'|'omega'|'T'>('u');
  const [field, setField] = useState<Field>('speed');

  // 四条边的枚举，确保 TS 能识别
  const edges: Edge[] = ['top','bottom','left','right'];

  // —— 稳定性检查相关 —— 
  // 用来记录上一次提示的 dt_max，避免重复弹窗
  const lastPromptedMaxDt = useRef<number | null>(null);
  // 用来显示拒绝自动修改后的警告
  const [manualWarning, setManualWarning] = useState<string | null>(null);


  // —— 新增：加载与完成提示 —— 
  const [isLoading, setIsLoading] = useState(false);
  const [simMsg, setSimMsg]       = useState<string | null>(null);

  // 1D 曲线剖面的方向：x 方向或 y 方向
  const [dir, setDir] = useState<'x'|'y'>('x');

  // 稳定性检查 & 弹窗
  useEffect(() => {
    const { Lx, Ly, Nx, Ny, U, Re, dt } = params;
    const dx = Lx / (Nx - 1);
    const dy = Ly / (Ny - 1);
    const nu = U * Ly / Re;

    const dtConvX  = dx / 2*U;
    const dtConvY  = dy / 2*U;
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
        `Current dt = ${dt.toExponential(2)}, exceeding the upper limit of stability dt_max = ${dtMax.toExponential(2)}。\n` +
        `It is recommended to adjust dt to ${suggested.toExponential(2)}。\n` +
        `Automatically adjust?`;
      if (window.confirm(msg)) {
        // 用户接受，自动修改 dt
        setParams(p => ({ ...p, dt: suggested }));
        setManualWarning(null);
      } else {
        // 用户拒绝，保留原 dt，并显示警告
        setManualWarning(
          `Please adjust dt ≤ ${dtMax.toExponential(2)}, otherwise the simulation may become unstable.`
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

  // —— 1D 剖面位置列表 —— 
  const [positions, setPositions] = useState<number[]>([0.1, 0.2, 0.5, 1.0]);
  const addPosition = () => {
    if (positions.length < 4) setPositions(ps => [...ps, 0]);
  };
  const updatePosition = (i: number, v: number) =>
    setPositions(ps => ps.map((p, idx) => idx === i ? v : p));
  const removePosition = (i: number) =>
    setPositions(ps => ps.filter((_, idx) => idx !== i));

  // —— 模拟结果 & 动画控制状态 —— 
  const [data, setData]             = useState<SimResult | null>(null);
  const [frameIndex, setFrameIndex] = useState(0);
  const [isPlaying, setIsPlaying]   = useState(false);
  const timerRef                    = useRef<number | null>(null);

  // —— 执行模拟 —— 
  const runSim = () => {
    setIsLoading(true);
    setSimMsg(null);
    simulate(runParams).then(res => {
      setData(res);
      setFrameIndex(0);
      setIsLoading(false);
      setSimMsg('Simulation Completed!');
      // 3 秒后自动隐藏
      setTimeout(() => setSimMsg(null), 3000);
      setIsPlaying(false);
    });
  };

  // 首次挂载或 runParams 改变时
  useEffect(runSim, [runParams]);

  // —— 播放/暂停 & 自动帧推进 —— 
  useEffect(() => {
    if (!data) return;
    if (isPlaying) {
      timerRef.current = window.setInterval(() => {
        setFrameIndex(i => (i + 1) % data.u_frames.length);
      }, 500);
    } else if (timerRef.current !== null) {
      clearInterval(timerRef.current);
      timerRef.current = null;
    }
    return () => {
      if (timerRef.current !== null) {
        clearInterval(timerRef.current);
        timerRef.current = null;
      }
    };
  }, [isPlaying, data]);


  // 2. 简单的表单栏
  return (
    <div style={{ display: 'flex', height: '100vh' }}>
      {/* —— 侧边栏 —— */}
      <aside style={{ width: 300, padding: 16, borderRight: '1px solid #ccc', overflowY:'auto' }}>
        <h3>Simulation Parameters</h3>

        {/* 网格大小Nx, Ny */}
        <div style={{ display: 'flex', gap: '16px', alignItems: 'center' }}>
          <label style={{ display:'flex', alignItems:'center', margin:0 }}>
            <span style={{ fontFamily: 'Times New Roman, serif' }}>N</span>
            <sub style={{ fontFamily: 'Times New Roman, serif' }}>x</sub>:
            <input
              type="number" min={10}
              value={params.Nx}
              onChange={e => setParams(p => ({ ...p, Nx: +e.target.value }))}
              style={{ width: 60, marginLeft: 8 }}
            />
          </label>
          <label style={{ display:'flex', alignItems:'center', margin:0 }}>
            <span style={{ fontFamily: 'Times New Roman, serif' }}>N</span>
            <sub style={{ fontFamily: 'Times New Roman, serif' }}>y</sub>:
            <input
              type="number" min={10}
              value={params.Ny}
              onChange={e => setParams(p => ({ ...p, Ny: +e.target.value }))}
              style={{ width: 60, marginLeft: 8 }}
            />
          </label>
        </div>

        {/* —— 域尺寸Lx, Ly —— */}
        <div style={{ display: 'flex', gap: '16px', alignItems: 'center' }}>
          <label style={{ display:'flex', alignItems:'center', margin:0 }}>
            <span style={{ fontFamily: 'Times New Roman, serif' }}>L</span>
            <sub style={{ fontFamily: 'Times New Roman, serif' }}>x</sub>:
            <input
              type="number" step="0.1"
              value={params.Lx}
              onChange={e => setParams(p => ({ ...p, Lx: +e.target.value }))}
              style={{ width: 43, marginLeft: 10 }}
            />
            <span style={{ marginLeft: 4, fontFamily: 'Times New Roman, serif' }}>m</span>
          </label>
          <label style={{ margin: 0 }}>
            <span style={{ fontFamily: 'Times New Roman, serif' }}>L</span>
            <sub style={{ fontFamily: 'Times New Roman, serif' }}>y</sub>:
            <input
              type="number" step="0.1"
              value={params.Ly}
              onChange={e => setParams(p => ({ ...p, Ly: +e.target.value }))}
              style={{ width: 43, marginLeft: 10 }}
            />
            <span style={{ marginLeft: 4, fontFamily: 'Times New Roman, serif' }}>m</span>
          </label>
        </div>

        {/* —— 时间参数dt, nt —— */}
        <label>
            <span style={{ fontFamily: 'Times New Roman, serif' }}>Δt</span>:
          <input
            type="number" step="0.001"
            value={params.dt}
            onChange={e => setParams(p => ({ ...p, dt: +e.target.value }))}
            style={{ width: 174, marginLeft: 12 }}
          />
        </label>
        <div>Maximum Δt: <strong>{maxDt.toExponential(2)}</strong></div>
        {manualWarning && (
          <div style={{ color:'red', marginTop:4 }}>
            {manualWarning}
          </div>
        )}
        <br/>
        <label>
          Nt:&nbsp;
          <input
            type="number" min={1}
            value={params.nt}
            onChange={e => setParams(p => ({ ...p, nt: +e.target.value }))}
            style={{ width: 173, marginLeft: 6 }}
          />
        </label>
        <hr/>

        <label>
          Problem Type:&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
          <select
            value={params.problemType}
            onChange={e=>setParams(p=>({
              ...p,
              problemType: e.target.value as 'channel'|'cavity'
            }))}
            style={{ width: 145  }}
          >
            <option value="channel">2D Channel Flow</option>
            <option value="cavity">2D Lid-Driven Cavity</option>
          </select>
        </label>

        {/* 数值方法选择 */}
        <label>Numerical Method:&nbsp;
          <select
            value={params.method}
            onChange={e => setParams(p => ({...p,method: e.target.value as SimParams['method']
            }))}
            style={{ width: 144  }}
          >
            <option value="projection">Projection</option>
            <option value="vorticity">Vorticity-Stream</option>
          </select>
        </label><hr/>

        {/* 物理参数 */}
        <label> Density <span style={{ fontFamily: 'Times New Roman, serif' }}>ρ</span>:
          <input
            type="number" step="0.1"
            value={params.rho}
            onChange={e => setParams(p => ({ ...p, rho: +e.target.value }))}
            style={{ width: 80, marginLeft: 36 }}
          />
          <span style={{ marginLeft: 4, fontFamily: 'Times New Roman, serif' }}>kg/m³</span>
        </label><br/>
        <label>
          Inlet speed U:
          <input
            type="number" step="0.1"
            value={params.U}
            onChange={e => setParams(p => ({ ...p, U: +e.target.value }))}
            style={{ width: 80, marginLeft: 8 }}
          />
          <span style={{ marginLeft: 4, fontFamily: 'Times New Roman, serif' }}>m/s</span>
        </label>
        <br/>
        <label>
          Reynolds Re:
          <input
            type="number" step="1"
            value={params.Re}
            onChange={e => setParams(p => ({ ...p, Re: +e.target.value }))}
            style={{ width: 80, marginLeft: 15 }}
          />
        </label>
        <hr/>


        {/* 展示场 */}
        <label>
          Simulation Results:&nbsp;
          <select value={field} onChange={e=>setField(e.target.value as any)}>
            {params.method === 'projection' ? (
              <>
                <option value="speed">speed</option>
                <option value="u">u</option>
                <option value="v">v</option>
                <option value="p">p</option>
              </>
            ) : (
              <>
                <option value="speed">speed</option>
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
        <h4>Temperature BC</h4>
        {edges.map(edge => (
          <div key={edge} style={{ display: 'flex', alignItems: 'center',marginBottom: 8 }}>
            {/* <strong>{edge}</strong>&nbsp; */}
            <strong style={{ width: 70 }}>{edge}</strong>
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
            {/* 根据类型显示不同单位 */}
            <span style={{ fontStyle: 'italic' }}>
              {params.tempBC[edge].type === 'Dirichlet' ? 'K' : 'K/m'}
            </span>
          </div>
        ))}
        <hr/>

        {/* Run & Animation Controls */}
        <button onClick={()=>{
          setRunParams(params);
          setRunKey(k=>k+1);
        }}disabled={isLoading}
          style={{ opacity: isLoading ? 0.6 : 1 }}
        >
          {isLoading ? 'Simulating…' : 'Run Simulation'}
        </button>

         {/* 模拟完成提示 */}
          {simMsg && (
            <div style={{
              marginTop: 8,
              padding: '4px 8px',
              background: '#e0ffe0',
              border: '1px solid #0a0',
              borderRadius: 4,
              color: '#060'
            }}>
              {simMsg}
            </div>
          )}

        <div style={{marginTop:16}}>
          <button onClick={()=>setIsPlaying(p=>!p)}>
            {isPlaying ? 'Pause' : 'Play'}
          </button>
          {data && (
            <input
              type="range"
              min={0}
              max={data.u_frames.length-1}
              value={frameIndex}
              onChange={e=>setFrameIndex(+e.target.value)}
              style={{width:'100%',marginTop:8}}
            />
          )}
        </div>

        <hr/>
        {/* 1D Profile Settings */}
        <h4>1D Profiles (Relative Position)</h4>
        {/* 方向选择 */}
         <label>
           Direction:&nbsp;
           <select
             value={dir}
             onChange={e => setDir(e.target.value as 'x'|'y')}
           >
             <option value="x">Along x (vary x, y-axis)</option>
             <option value="y">Along y (vary y, x-axis)</option>
           </select>
         </label>
         <br/>

        {positions.map((p,i)=>(
          <div key={i} style={{display:'flex',alignItems:'center',marginBottom:4}}>
            <input
              type="number" step="0.01" min={0} max={1}
              value={p}
              onChange={e=>updatePosition(i,+e.target.value)}
              style={{width:60}}
            />
            <button onClick={()=>removePosition(i)} style={{marginLeft:4}}>×</button>
          </div>
        ))}
        {positions.length<4 && <button onClick={addPosition}>+ Add</button>}
      </aside>

      {/* —— 主视图区：总标题 + 双图 —— */}
    <main
      style={{
        flex: 1,
        display: 'flex',
        flexDirection: 'column',
        background: '#F7F7F7'
      }}
    >
      {/* —— 全局总标题 —— */}
      <div
        style={{
          padding: '12px 0',
          textAlign: 'center',
          background: '#fff',
          borderBottom: '1px solid #eee',
          marginBottom: 8
        }}
      >
        <h2 style={{ margin: 0, fontSize: '1.3rem', color: '#333' }}>
          {runParams.problemType === 'channel'
            ? '2D Channel Flow'
            : '2D Lid-Driven Cavity'}{' '}
          —{' '}
          {runParams.method === 'projection'
            ? 'Projection Method'
            : 'Vorticity–Stream Function'}{' '}
          — Results: {field.toUpperCase()}
        </h2>
      </div>

      {/* —— 图表区域：左右两栏 —— */}
      <div style={{ flex: 1, display: 'flex' }}>
        {/* —— 左边 2D 等高线 —— */}
        <div style={{ flex: 1, borderRight: '1px solid #eee', background: '#fff' }}>
          <ContourPlot
            key={runKey}
            params={runParams}
            field={field}
            data={data}
            frameIndex={frameIndex}
          />
        </div>
        {/* —— 右边 1D 剖面 —— */}
        <div style={{ flex: 1, paddingLeft: 8, background: '#fff' }}>
          <ProfilePlot
            data={data}
            params={runParams}
            field={field}
            frameIndex={frameIndex}
            positions={positions}
            direction={dir}
          />
        </div>
      </div>
    </main>
  </div>
  );
}

export default App;
