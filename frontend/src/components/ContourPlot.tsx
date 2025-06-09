// src/components/ContourPlot.tsx

import React from 'react';
import Plot from 'react-plotly.js';
import { SimParams, SimResult, Field } from '../api';

interface ContourPlotProps {
  data: SimResult | null;                 // 从父组件传入
  params: SimParams;                      // 用于 totalTime 计算
  field:      Field;
  frameIndex: number;                     // 当前帧索引
}

const ContourPlot: React.FC<ContourPlotProps> = ({
  data,
  params,
  field,
  frameIndex,
}) => {
  if (!data) {
    return <div>Simulating, please wait...</div>;
  }

  // 1. 取出当前帧的二维场
  let zData: number[][];
  if (field === 'speed') {
    const U2 = data.u_frames[frameIndex];
    const V2 = data.v_frames[frameIndex];
    zData = U2.map((row, i) =>
      row.map((u, j) => Math.hypot(u, V2[i][j]))
    );
  } else {
    switch (field) {
      case 'u':     zData = data.u_frames[frameIndex];    break;
      case 'v':     zData = data.v_frames[frameIndex];    break;
      case 'p':     zData = data.p_frames![frameIndex];   break;
      case 'psi':   zData = data.psi_frames![frameIndex]; break;
      case 'omega': zData = data.omega_frames![frameIndex]; break;
      case 'T':     zData = data.T_frames![frameIndex];   break;
    }
  }

  // 2. 计算等高线范围和步长
  const flatZ = zData.flat();
  const zMin  = Math.min(...flatZ);
  const zMax  = Math.max(...flatZ);
  const nLevels = 20;
  const step    = (zMax - zMin) / (nLevels - 1);

  // 3. 字段标签
  const labelMap: Record<typeof field, string> = {
    speed: 'Speed (m/s)',
    u:     'Velocity u (m/s)',
    v:     'Velocity v (m/s)',
    p:     'Pressure (Pa)',
    psi:   'Stream Function ψ',
    omega: 'Vorticity ω',
    T:     'Temperature (K)',
  };

  const t = data.t_frames[frameIndex];
  const total = params.dt * params.nt;

  return (
    <Plot
      data={[{
        z:          zData,
        x:          data.x,
        y:          data.y,
        type:       'contour',
        colorscale: 'Jet',
        autocontour: false,
        contours: {
          start:      zMin,
          end:        zMax,
          size:       step,
          coloring:   'heatmap',
          showlines:  false,
        },
        line: {
          width: 0.2,
          color: 'white'
        },
        colorbar: {
          title: { text: labelMap[field], font: { size: 14 } },
        },
      }]}
      layout={{
        // —— 全局字体 —— 
        font: {
          family: 'Times New Roman, serif',
          size:   22,
          color:  '#333'
        },
        // —— 主标题 —— 
        title: {
          text: `${labelMap[field]} — t=${t.toFixed(3)}s / total=${total.toFixed(3)}s`,
        },
        // —— 坐标轴 —— （会自动继承 layout.font）
        xaxis: {
          title: {
            text: 'x (m)',
            font: { family:'Times New Roman, serif', size:25 }
          }
          // tickfont 不必写，直接继承全局 font.size=14
        },
        yaxis: {
          title: {
            text: 'y (m)',
            font: { family:'Times New Roman, serif', size:25 }
          },
          scaleanchor: 'x'
        },
    
        // autosize: true,
        // margin: { l:60, r:60, t:80, b:60 }
      }}
      useResizeHandler
      style={{ width:'100%', height:'100%' }}
    />
  )
}

export default ContourPlot;
