// src/components/ProfilePlot.tsx

import React from 'react';
import Plot from 'react-plotly.js';
import { SimParams, SimResult, Field } from '../api';

interface Props {
  data: SimResult | null;
  params: SimParams;
  // field: 'u'|'v'|'p'|'psi'|'omega'|'T';
  field:      Field;
  frameIndex: number;
  positions: number[];   // x/Lx 数组
  direction: 'x' | 'y';     // 新增
}

const ProfilePlot: React.FC<Props> = ({
  data, params, field, frameIndex, positions, direction
}) => {
  if (!data) return <div>Please run a simulation first...</div>;

  // 1) 取当前帧的二维场
  let field2D: number[][];
  if (field === 'speed') {
    const U2 = data.u_frames[frameIndex];
    const V2 = data.v_frames[frameIndex];
    field2D = U2.map((row,i) =>
      row.map((u,j) => Math.hypot(u, V2[i][j]))
    );
  } else {
    field2D = (
      field === 'u' ? data.u_frames[frameIndex] :
      field === 'v' ? data.v_frames[frameIndex] :
      field === 'p' ? data.p_frames![frameIndex] :
      field === 'psi' ? data.psi_frames![frameIndex] :
      field === 'omega' ? data.omega_frames![frameIndex] :
      data.T_frames![frameIndex]
    );
  }


  // 3) 定义样式数组
  const lineStyles = ['solid','dash','dashdot','dot'] as const;
  const markers    = ['circle','square','triangle-up','diamond'] as const;

  // 3) 构造 traces
  const traces = positions.map((p,i) => {
    if (direction === 'x') {
      // x-dir: 剖面在固定 x，纵轴 y
      const ix = Math.round(p * (data.x.length - 1));
      const prof = data.y.map((_, j) => field2D[j][ix]);
      return {
        x: prof,
        y: data.y,
        mode: 'lines+markers' as const,
        name: `x=${(p*params.Lx).toFixed(2)}m`,
        line:   { dash: lineStyles[i], width: 2 },
        marker: { symbol: markers[i], size: 6 },
      };
    } else {
      // y-dir: 剖面在固定 y，横轴 x
      const iy = Math.round(p * (data.y.length - 1));
      const prof = data.x.map((_, j) => field2D[iy][j]);
      return {
        x: data.x,
        y: prof,
        mode: 'lines+markers' as const,
        name: `y=${(p*params.Ly).toFixed(2)}m`,
        line:   { dash: lineStyles[i], width: 2 },
        marker: { symbol: markers[i], size: 6 },
      };
    }
  });

  // 5) 坐标标签
  const physLabel  = field === 'p' ? 'Pressure (Pa)' : 
                     field === 'T' ? 'Temperature (K)' : 
                     field === 'speed' ? 'Speed (m/s)':
                     `Velocity ${field} (m/s)`;

  // 6) 当前时间
  const t = data.t_frames[frameIndex];

  return (
    <Plot
      data={traces}
      layout={{
        // —— 全局字体 —— 
        font: {
          family: 'Times New Roman, serif',
          size:   22,
          color:  '#333'
        },
        title: {
          text: `${field.toUpperCase()} Profiles (${direction}-dir) — t=${data.t_frames[frameIndex].toFixed(3)} s`
        },
        xaxis: { title: { text: direction==='x' ? physLabel : 'x (m)' } },
        yaxis: { title: { text: direction==='x' ? 'y (m)' : physLabel } },

        // —— 这里添加 legend 配置 —— 
        showlegend: true,
        legend: {
          // 把图例放在画布内部 (x:0–1, y:0–1)
          x: 0.75,           // 距左侧 2% 的位置
          y: 0.28,           // 距顶部 2% 的位置
          xanchor: 'left',   // x=0.02 指的是图例左边
          yanchor: 'top',    // y=0.98 指的是图例顶部
          bgcolor: 'rgba(255, 255, 255, 0.7)',  // 半透明背景
          bordercolor: '#ccc',
          borderwidth: 1,
          font: {
            family: 'Times New Roman, serif',
            size:   22,
            color:  '#333'
          }
        },

        autosize: true,
      }}
      useResizeHandler
      style={{ width:'100%', height:'100%' }}
    />
  );
};


export default ProfilePlot;
