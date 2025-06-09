# CFD-Webapp
A web solver for 2D channel flow and 2D lid-driven cavity simulations. There are many parameters for you to customize the input.
# 0. 环境准备

## 0.1 安装 Node.js 与 npm

前往 [https://nodejs.org](https://nodejs.org/) 下载并安装稳定版（node.js 默认配置了npm）

```bash
node -v   # 查看版本，确保 >=14
npm -v    # 查看版本
```

## 0.2 后端依赖

```bash
pip install fastapi uvicorn numpy
```

## 0.3 cfd-web 项目创建

### 0.3.1 前端：React + TypeScript

```bash
cd cfd-web
npx create-react-app frontend --template typescript
cd frontend
npm install react-plotly.js plotly.js axios
npm start
```

### 0.3.2 将核心文件copy进项目

创建文件：

src/App.tsx  (替换其中的内容即可)

src/api.ts      (如果新建这个文件后，App.tsx找不到./api，重启IDE)

src/components/ContourPlot.tsx
src/components/ProfilePlot.tsx     (可以直接将提交的components文件夹复制过去)

### 0.3.2 后端（在提交文件中已包含）

创建文件

backend/app/main.py

backend/app/solver.py

# 1. 启动Web

## 1.1 项目结构

```python
cfd-web/
├── backend/          # FastAPI 服务
│   ├── app/
│   │   ├── **main.py**           # 应用入口
│   │   └── **solver.py**         # projection 方法实现
│   ├── requirements.txt
│   └── .venv/                # Python venv
└── frontend/         # React 应用
    ├── public/
    ├── src/
    │   ├── components/
    │   │   └── **ContourPlot.tsx**
    │   │   └── **ProfilePlot.tsx**
    │   ├── **App.tsx    # 界面相关**
    │   ├── index.tsx
    │   └── api.ts             # 封装 fetch 调用
    ├── package.json
    └── tsconfig.json
```

## 1.2 启动后端服务

`cd backend`

`uvicorn app.main:app --reload --host 0.0.0.0 --port 8000`
(打开浏览器访问 `http://localhost:8000/docs`，可以看到自动生成的 Swagger 文档)

## 1.3 前端

1. 创建 React+TypeScript 项目（首次创建项目时，后续不需要）
    1. 在 `cfd-web` 根目录下运行（首次编译）：
        
        `npx create-react-app frontend --template typescript`
        
        `npm install react-plotly.js plotly.js axios`
        
        (如果报错无法找到模块“react-plotly.js”) `npm install --save-dev @types/react-plotly.js`
        
2. `npm start`
3. 浏览器打开 `http://localhost:3000`，即可看到左侧参数面板，右侧等值线实时更新的动画。
4. **停止服务**：在对应 PowerShell 窗口按 `Ctrl+C` 即可终止后端或前端服务器。    
