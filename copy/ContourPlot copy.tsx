// src/components/ContourPlot.tsx

import React, { useEffect, useState, useRef } from 'react';
import Plot from 'react-plotly.js';
import { SimParams, SimResult, simulate } from '../api';

interface ContourPlotProps {
  params: SimParams        // 包含 method 和 tempBC
  field: 'u'|'v'|'p'|'psi'|'omega'|'T'
  interval?: number
}

const ContourPlot: React.FC<ContourPlotProps> = ({
  params,
  field,
  interval = 500,
}) => {
  const [data, setData] = useState<SimResult | null>(null);
  const [frameIndex, setFrameIndex] = useState(0);
  const timerRef = useRef<number | null>(null);

  // 请求模拟
  useEffect(() => {
    setData(null);
    simulate(params).then(res => {
      setData(res);
      setFrameIndex(0);
    });
    return () => {
      if (timerRef.current !== null) {
        window.clearInterval(timerRef.current);
        timerRef.current = null;
      }
    };
  }, [params]);

  // 启动动画
  useEffect(() => {
    if (!data) return;
    if (timerRef.current !== null) {
      window.clearInterval(timerRef.current);
    }
    timerRef.current = window.setInterval(() => {
      setFrameIndex(idx =>
        (idx + 1) % data.u_frames.length
      );
    }, interval);
    return () => {
      if (timerRef.current !== null) {
        window.clearInterval(timerRef.current);
        timerRef.current = null;
      }
    };
  }, [data, interval]);

  if (!data) {
    return <div>正在模拟，请稍候…</div>;
  }

// 在组件里，拿到 data: SimResult 之后
const zData = (() => {
  switch (field) {
    case 'u':     return data.u_frames[frameIndex]
    case 'v':     return data.v_frames[frameIndex]
    case 'p':     return data.p_frames![frameIndex]
    case 'psi':   return data.psi_frames![frameIndex]
    case 'omega': return data.omega_frames![frameIndex]
    case 'T':     return data.T_frames![frameIndex]
  }
})()

  // 扁平化并计算 z 的最小/最大值，用于等高线级别划分
  const flatZ = zData.reduce<number[]>(
    (acc, row) => acc.concat(row),
    []
  );
  const zMin = Math.min(...flatZ);
  const zMax = Math.max(...flatZ);

  // 想要 50 条等高线
  const nLevels = 12;
  const step = (zMax - zMin) / (nLevels - 1);
  const labelMap: Record<typeof field, string> = {
  u:     'Velocity u (m/s)',
  v:     'Velocity v (m/s)',
  p:     'Pressure (Pa)',
  psi:   'Stream Function ψ',
  omega: 'Vorticity ω',
  T:     'Temperature (K)',
}

  // 计算总模拟时间
  const totalTime = params.dt * params.nt;  // 秒

  return (
    <Plot
      data={[{
        z:          zData,
        x:          data.x,
        y:          data.y,
        type:       'contour',
        colorscale: 'Jet',
        autocontour:false,
        contours: {
          start: zMin,  end: zMax,
          size:  step, 
          coloring: 'heatmap',
          showlines: false
        },
        colorbar: {
          title: { text: labelMap[field] }
        },
        // 如果你想给速度叠加流线（仅 u/v），可以：
        ...(field==='u' || field==='v'
          ? {  
              stream: {  // Plotly 并不内置 stream 参数，这需要你改成两个 trace：一个 heatmap，一个 streamline
                type: 'streamline',
                u: data.u_frames[frameIndex],
                v: data.v_frames[frameIndex],
              }
            }
          : {}
        )
      }]}
      layout={{
        title: { text: `${labelMap[field]} — Time: ${data.t_frames[frameIndex].toFixed(3)} s / Total: ${totalTime.toFixed(3)} s` },
        xaxis: { title: { text:'x'} },
        yaxis: { title: { text:'y'}, scaleanchor:'x' },
        autosize:true
      }}
      useResizeHandler
      style={{ width:'100%', height:'100%' }}
    />

  );
};

export default ContourPlot;
