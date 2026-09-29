# MSPC — 质谱同位素分布预计算 API

基于 **Python + FastAPI + SQLite** 的质谱预计算服务。客户端提交含元素计数的分子式（可带同位素标记），服务按各元素天然同位素丰度将分布表示为"质量-强度多项式"，逐元素卷积得到完整同位素包络，合并质量差小于 0.0001 Da 的峰，按基峰归一化后返回前 N 个峰及分子单同位素质量。

## 快速启动（单条命令）

```bash
docker compose up --build
```

启动过程会自动完成依赖安装（镜像构建时）与 SQLite 数据库初始化（服务启动时建表），无需任何手工迁移。启动后服务监听 **8000 端口**，API 根地址为 `http://localhost:8000`（容器外同样通过该地址访问）。

不想使用 compose 时，也可以用单条命令构建镜像并直接启动容器（数据写入容器内 `/data`，删除容器即清除；加 `-v mspc-data:/data` 可持久化）：

```bash
docker run --rm -p 8000:8000 $(docker build -q .)
```

服务启动后，API 即位于 `http://localhost:8000`（例如质量轴校准接口 `http://localhost:8000/api/v1/mass-calibrate`、共识谱接口 `http://localhost:8000/api/v1/consensus-spectrum`；Swagger 文档：`http://localhost:8000/docs`）。

一条命令启动时所用的关键配置项（均可在 `docker-compose.yml` 的 `environment` 段覆盖，完整列表见[配置（环境变量）](#配置环境变量)）：

| 配置项 | 默认值 | 说明 |
|---|---|---|
| `MSPC_DB_PATH` | `/data/mspc.db` | SQLite 数据库文件路径（compose 已挂载命名卷 `mspc-data`） |
| `MSPC_MERGE_TOLERANCE` | `0.0001` | 峰合并质量容差（Da），同时用于分子式反查的包络生成 |
| `MSPC_PRUNE_THRESHOLD` | `1e-12` | 卷积过程中按相对强度剪枝的阈值 |
| `MSPC_TOP_N_MAX` / `MSPC_TOP_N_DEFAULT` | `200` / `20` | `top_n` 上限与默认值 |

容器端口固定映射到宿主机 **localhost 的 8000 端口**：`ports: ["8000:8000"]`。

## 访问地址

| 地址 | 说明 |
|---|---|
| http://localhost:8000 | 服务根地址（端口 8000） |
| http://localhost:8000/docs | Swagger 交互式文档 |
| http://localhost:8000/health | 健康检查 |

## API 用法

### `POST /api/v1/envelope`

请求体：

```json
{
  "formula": "C6H12O6",
  "top_n": 5
}
```

- `formula`（必填）：分子式，如 `C6H12O6`；同位素标记用方括号表示，如 `C6H6[13C]2`（2 个纯 ¹³C 原子）、`C2H5[15N]O[18O]`。
- `top_n`（可选，默认 20，最大 200）：返回强度最高的前 N 个峰。

响应：

```json
{
  "formula": "C6H12O6",
  "composition": {"C": 6, "H": 12, "O": 6},
  "labels": {},
  "monoisotopic_mass": 180.063388,
  "base_peak_intensity": 100.0,
  "peak_count": 5,
  "peaks": [
    {"mass": 180.063388, "intensity": 100.0},
    {"mass": 181.066743, "intensity": 6.489437}
  ],
  "cached": false
}
```

- `monoisotopic_mass`：单同位素质量（各元素取天然丰度最高的同位素，标记原子取标记同位素质量）。
- `peaks`：基峰强度归一化为 100，按质量升序返回前 N 个峰。
- `cached`：结果是否来自 SQLite 缓存（相同 `formula` + `top_n` 的重复请求直接命中缓存）。

### `POST /api/v1/ion-envelope`

在中性分子式上施加**加合物**并返回带电离子的 m/z 同位素包络：先按加合物片段增减原子得到离子组成，复用现有分子式解析与多项式卷积生成中性质量包络，再转换为 m/z。

请求体：

```json
{
  "formula": "C6H12O6",
  "adduct": "+H",
  "charge": 1,
  "top_n": 5
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/envelope`（可含 `[13C]` 等纯同位素标记）。
- `adduct`（必填）：由连续的"+/-[系数]分子式片段"组成，每个片段必须以显式 `+` 或 `-` 开头；系数为正整数（缺省为 1），**乘整个片段**。片段本身的语法与分子式相同（元素计数、`[13C]` 标记等），整个字符串必须被完整消费。例如：
  - `+H`：质子化（[M+H]⁺）；
  - `+Na-H`：钠离子加合物减一个氢；
  - `+2H-H2O`：加两个氢并脱去一个水分子（系数 2 作用于整个 `H`，即加 2 个 H）；
  - `-[13C]`：减去一个纯 ¹³C 原子。
  - 应用片段后任一元素（或同位素标记）计数变为负数时返回 422（零允许），`position` 指向导致**首次**负计数的片段首字符（`+`/`-` 与系数之后的第一个字符）。
- `charge`（必填）：带符号整数电荷，取值 **-10～-1 或 1～10**（不允许 0）。
- `top_n`（可选，默认 20，最大 200）：返回强度最高的前 N 个峰。

计算流程：

1. 用 `|charge| × MSPC_MERGE_TOLERANCE` 作为中性质量峰的合并容差（卷积及最终合并均用该放大容差，按默认 0.0001 Da 即电荷 2 时为 0.0002 Da）。
2. 对调整后的离子组成复用同一套卷积/剪枝/归一化（基峰 = 100）。
3. 每个峰按下式从离子质量（Da）转为 m/z（正离子扣除电子质量，负离子加回电子质量）：

   ```
   m/z = (mass - charge × 0.000548579909065 Da) / |charge|
   ```

响应：

```json
{
  "formula": "C6H12O6",
  "adduct": "+H",
  "charge": 1,
  "composition": {"C": 6, "H": 13, "O": 6},
  "labels": {},
  "monoisotopic_mz": 181.070665,
  "base_peak_intensity": 100.0,
  "peak_count": 5,
  "peaks": [
    {"mz": 181.070665, "intensity": 100.0},
    {"mz": 182.07402, "intensity": 6.489437}
  ]
}
```

- `composition` / `labels`：应用加合物并乘系数之后的离子组成（规则同 `/envelope` 的 `labels` 序列化）。
- `monoisotopic_mz`：调整后组成的单同位素中性质量经上式转换得到的单同位素 m/z。
- `peaks`：基峰归一化为 100 的前 N 个 m/z 峰，按 m/z 升序；双键电荷下相邻同位素峰的 m/z 间距约为 1/|charge| Da。

示例：`C6H12O6` + `+2H`、`charge` 2 的单同位素 m/z 为 `(180.063388 + 2×1.007825 - 2×0.0005486) / 2 ≈ 91.038971`；`-H`、`charge` -1 的去质子负离子 m/z 约为 `179.056112`。

服务默认监听 **8000 端口**（见 [访问地址](#访问地址)，容器内外均为 `http://localhost:8000`），直接用 curl 请求示例：

```bash
curl -X POST http://localhost:8000/api/v1/ion-envelope \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C6H12O6","adduct":"+Na-H","charge":1,"top_n":5}'
```

### `POST /api/v1/neutral-profile`

将**中性分子的同位素棒谱**（stick spectrum）与**仪器轮廓**（instrument profile）卷积，模拟连续质量轴上的高斯展宽谱线。先复用 `/api/v1/envelope` 的分子式解析与多项式卷积得到全部同位素棒（不做 `top_n` 截断），再为每根质量为 m 的棒生成以该质量为中心的高斯核：

```
FWHM(m) = m / resolving_power
σ       = FWHM / (2√(2 ln 2))
```

每个核仅在棒质量的 **±4σ** 内采样，并累加到**统一质量网格**（间距 `mass_step`，单位 Da）；落在请求区间外的核被截断，但区间内的贡献原样保留（不做截断后归一化）。全部核累加完成后，将最高强度归一化为 **100**，按质量升序返回整条 profile。

请求体：

```json
{
  "formula": "C6H12O6",
  "resolving_power": 100000,
  "mass_step": 0.001,
  "mass_min": 180.0,
  "mass_max": 190.0
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/envelope`（可含 `[13C]` 等纯同位素标记）。
- `resolving_power`（必填，>0）：仪器分辨率 R = m/FWHM；不同质量的棒按各自 FWHM = m/R 展宽（质量越大峰越宽）。
- `mass_step`（必填，>0）：输出统一质量网格的间距（Da）。
- `mass_min` / `mass_max`（可选，>0）：输出质量区间（Da，闭区间）。**两者均可省略**：
  - 都省略：区间取**所有棒峰各自 ±4σ 支撑窗的并集**边界（最小下界、最大上界）；
  - 只给一个：另一个边界取上述并集边界；
  - 网格必须为非空区间（`mass_min < mass_max`），且总点数 **最多 20000**（点数 = `floor((mass_max − mass_min)/mass_step) + 1`），否则返回 422。

响应：

```json
{
  "formula": "C6H12O6",
  "resolving_power": 100000.0,
  "mass_min": 180.0,
  "mass_max": 190.0,
  "mass_step": 0.001,
  "point_count": 10001,
  "profile": [
    {"mass": 180.0, "intensity": 2.315421},
    {"mass": 180.001, "intensity": 2.453812}
  ]
}
```

- `mass_min` / `mass_max`：实际使用的区间（省略时回填为 ±4σ 并集边界）。
- `point_count`：统一网格点数（`len(profile)`），上限 20000。
- `profile`：`mass` 从 `mass_min` 起按 `mass_step` 等间距递增、最高强度归一化为 100 的连续轮廓；区间外的高斯核被截断、区间内贡献保留，因此区间边缘的峰可能因截断而矮于 100。
- 区间内没有任何网格点落入任一根棒的 ±4σ 支撑窗时（如窗口离谱图太远，或窗口窄于一个步长）返回 422。

```bash
curl -X POST http://localhost:8000/api/v1/neutral-profile \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C6H12O6","resolving_power":100000,"mass_step":0.001}'
```

### `POST /api/v1/profile-centroid`

将客户端提交的**等间距连续轮廓**（mass/intensity 采样点）转换为**质心峰**（centroid）。典型用法是把 `/api/v1/neutral-profile` 的输出直接回本接口，得到可用于富集拟合/离子识别的棒谱。

请求体：

```json
{
  "profile": [
    {"mass": 100.0, "intensity": 0.0},
    {"mass": 100.001, "intensity": 0.5},
    {"mass": 100.002, "intensity": 4.0}
  ],
  "relative_threshold": 1.0,
  "min_separation": 0.005
}
```

- `profile`（必填，**3～20000 个点**）：连续轮廓采样点，字段为 `mass`/`intensity`：
  - `mass` 必须**严格递增**；
  - 网格必须等间距——以第一个间距为参考步长，其后每个间距与参考步长的**相对误差不得超过 1e-6**（即 `|gap − step| / step ≤ 1e-6`），否则按首个违规间距的 0 基索引返回 422；
  - `intensity` 非负（允许 0），且**至少一个为正**；mass/intensity 都必须是有限数。
- `relative_threshold`（必填，**0～100**）：相对阈值（百分比）。轮廓先按基峰归一化（最高点 = 100），归一化强度**低于**该阈值的点在寻峰前丢弃（恰好等于阈值保留）。
- `min_separation`（必填，**>0**，Da）：非极大值抑制（NMS）的最小峰间距。

处理流程：

1. **归一化与阈值过滤**：按基峰归一化到 100，丢弃阈值以下的点（比较含 1e-9 量级的浮点容差，保证阈值 100 时基峰平台不被舍入误差误删）。
2. **平顶局部极大候选**：在阈值以上的点中，把每个"平顶局部极大区间"视为一个候选——一段强度相同、且两侧紧邻点（边界外侧视为无限低）都严格更低的**连续等强游程**。候选**位置取平台中点**（平台两端点质量的平均，可能落在两个网格点之间），候选强度取平台归一化强度。
3. **非极大值抑制**：候选按强度**降序**贪心处理，强度相同时**低质量优先**；若候选与任一已保留候选的距离小于 `min_separation`（恰等于 `min_separation` 允许，含浮点容差）则丢弃。
4. **分配与质心**：将每个阈值以上的点分配给质量最近的保留候选，**等距时归低质量候选**；每个候选用**原始（未归一化）强度**计算：
   - 质心质量 = Σ(m·I) / ΣI（只计分配给它的点）；
   - 峰顶归一化强度 = 平台原始强度 × 100 / 基峰原始强度（平顶取平台高度，而非分配区域内的最大值）；
   - 占总分配强度百分比 = 分配给它的原始强度和 / 全部阈值以上点原始强度和 × 100。
5. 峰按**质量升序**返回。

响应：

```json
{
  "point_count": 27,
  "relative_threshold": 1.0,
  "min_separation": 0.005,
  "candidate_count": 2,
  "peak_count": 2,
  "peaks": [
    {"mass": 100.007, "intensity": 100.0, "total_intensity_percent": 87.356322},
    {"mass": 100.021, "intensity": 25.0, "total_intensity_percent": 12.643678}
  ]
}
```

- `point_count`：输入轮廓点数；`candidate_count` 为 NMS 之前的平顶局部极大候选数；`peak_count` 为 NMS 之后实际返回的峰数。
- `peaks`：每个峰含质心 `mass`（强度加权平均质量）、`intensity`（峰顶相对基峰的归一化强度）与 `total_intensity_percent`（占全部阈值以上分配强度的百分比，所有峰之和为 100），按质量升序。
- 被 NMS 抑制候选附近的阈值以上点不会丢失，而是照常分配给最近的保留候选，因此保留峰的质心质量可能偏离其峰顶位置。

**大强度数值稳定性**：`intensity` 允许取任何非负有限浮点数（最大约 1.8×10³⁰⁸）。阈值筛选与候选峰的归一化强度都先计算 `intensity / base` 比值（必在 [0, 1]）再乘 100，避免 `100 × intensity` 在除以同量级基峰之前先溢出为无穷；各峰的原始强度加权质心与占总强度百分比在累计时即把强度除以基峰（每个权重落在 [0, 1]，累计和只随点数增长），因此既不会在 `m·I` 乘积上溢出，也不会在总强度和上溢出成 `inf/inf = NaN`（质量本身接近浮点上界时再按最大质量等比缩放一次）。因此质量为 100、101、102、强度为 8e307、1e308、8e307 的输入（直接计算 `100 × 1e308`、`102 × 1e308` 及三者强度之和都会溢出）仍返回有限结果：

```json
{
  "point_count": 3,
  "relative_threshold": 5.0,
  "min_separation": 0.5,
  "candidate_count": 1,
  "peak_count": 1,
  "peaks": [
    {"mass": 101.0, "intensity": 100.0, "total_intensity_percent": 100.0}
  ]
}
```

按同一因子整体放大或缩小所有强度不改变质心、归一化强度与百分比拆分；平顶峰判定、等距点归低质量候选的最近候选分配以及质量升序排序均不受缩放影响。

**极大质量平台的数值稳定性**：除强度外，平台两端点的 `mass` 也允许接近浮点上界（最大约 1.8×10³⁰⁸）。平顶候选的中点按 `left/2 + right/2`（先各除以 2 再相加）计算，而非会溢出的 `(left + right)/2`——两个端点都接近浮点上界时直接相加得到无穷大，会令相距超过 `min_separation` 的两个平台中点都变成 `inf`，NMS 以 `abs(inf − inf) = NaN` 把它们误判为"同位置"而抑制掉一个，随后最近候选归属也把全部点并给唯一幸存峰，最终响应只剩一个峰。无溢出中点同时参与**候选排序、最小间隔抑制与最近候选归属**；因此两个极大质量平台各自返回有限质心与各自强度份额（平台质量过大致使 `m·w` 累计溢出时，质心仍走既有的按最大成员质量等比缩放的兜底路径）。普通质量范围的结果与既有网格校验规则（严格递增、1e-6 等间距容差）保持完全一致。

例如质量约 1e308、网格步长 1e306 的两个平顶峰（两端点之和超过浮点最大值），`min_separation=2e306` 时两峰都保留，质心与份额均有限；接口地址为 `http://localhost:8000/api/v1/profile-centroid`：

```bash
curl -X POST http://localhost:8000/api/v1/profile-centroid \
  -H 'Content-Type: application/json' \
  -d '{"relative_threshold":0.0,"min_separation":2e306,
       "profile":[{"mass":1.0e308,"intensity":10.0},{"mass":1.01e308,"intensity":10.0},
                  {"mass":1.02e308,"intensity":0.0},{"mass":1.03e308,"intensity":5.0},
                  {"mass":1.04e308,"intensity":5.0}]}'
```

非法输入返回带字段（必要时带 0 基点索引）的 422：

- 点数不在 3～20000、强度为负或非有限：`error` 为 `invalid_profile`，`field` 为 `profile`（索引问题附 `index`）；
- 质量非严格递增或网格间距相对误差超过 1e-6：`error` 为 `invalid_mass_grid`，`field` 为 `profile`，`index` 指向首个违规间距的右端点；
- 全部强度为 0：`error` 为 `zero_profile`；
- 阈值以上不存在任何局部极大候选：`error` 为 `no_candidates`，`field` 为 `relative_threshold`（正常输入下基峰平台恒为候选，仅在极端情形触发）；
- `relative_threshold` 越界（不在 0～100）、`min_separation` ≤ 0、字段缺失等由请求校验返回标准 422。

```bash
curl -X POST http://localhost:8000/api/v1/profile-centroid \
  -H 'Content-Type: application/json' \
  -d '{"relative_threshold":1.0,"min_separation":0.005,
       "profile":[{"mass":100.0,"intensity":0.0},{"mass":100.001,"intensity":60.0},
                  {"mass":100.002,"intensity":100.0},{"mass":100.003,"intensity":60.0}]}'
```

大强度（约 1e308）输入的示例——返回的质心、归一化强度与百分比均为有限数：

```bash
curl -X POST http://localhost:8000/api/v1/profile-centroid \
  -H 'Content-Type: application/json' \
  -d '{"relative_threshold":5.0,"min_separation":0.5,
       "profile":[{"mass":100.0,"intensity":8e307},
                  {"mass":101.0,"intensity":1e308},
                  {"mass":102.0,"intensity":8e307}]}'
```

### `POST /api/v1/ion-identify`

在一个中性分子式上同时评估 **1～20 个候选加合物/电荷组合**：对每个候选复用 `/api/v1/ion-envelope` 的加合物解析与离子包络计算（含 `|charge| ×` 合并容差缩放与 m/z 换算、基峰归一化为 100），将预测 m/z 峰与观测 m/z 峰在容差内一对一配对，以全部配对（含未匹配峰）的强度均方误差对候选排序，返回全部候选的得分、配对明细与最佳候选。

请求体：

```json
{
  "formula": "C6H12O6",
  "peaks": [
    {"mz": 181.070665, "intensity": 100.0},
    {"mz": 182.07402, "intensity": 6.49}
  ],
  "mass_tolerance": 0.01,
  "candidates": [
    {"adduct": "+Na-H", "charge": 1},
    {"adduct": "+H", "charge": 1},
    {"adduct": "+H", "charge": 2},
    {"adduct": "-99H", "charge": -1}
  ]
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/ion-envelope`（可含 `[13C]` 等纯同位素标记）。
- `peaks`（必填，**1～200 个**）：观测峰，字段为 `mz`/`intensity`；m/z 必须**严格递增**，强度非负且至少一个为正。观测峰在评分前按基峰归一化为 100。
- `mass_tolerance`（必填，>0）：观测峰与预测 m/z 峰配对的质量容差（m/z 单位）。
- `candidates`（必填，**1～20 个**）：候选加合物/电荷组合，每个包含：
  - `adduct`（必填）：加合物字符串，语法与 `/api/v1/ion-envelope` 完全一致（`+/-[系数]片段`，如 `+H`、`+2H-H2O`、`+Na-H`）。
  - `charge`（必填）：带符号电荷，取值 **-10～-1 或 1～10**（不允许 0）。

配对与评分规则（与富集拟合接口同一套实现）：所有容差内的（观测, 预测）候选对按质量差绝对值从小到大贪心占用，每个峰至多参与一对；未匹配的观测峰以预测强度 0 计，未匹配的预测峰以观测强度 0 计；MSE 为全部配对强度差平方的均值。候选预测包络不经 `top_n` 截断（卷积剪枝阈值仍然生效）。

响应：

```json
{
  "formula": "C6H12O6",
  "mass_tolerance": 0.01,
  "observed_peak_count": 2,
  "best_candidate": {"index": 1, "adduct": "+H", "charge": 1, "mse": 0.0003},
  "candidates": [
    {
      "index": 1,
      "valid": true,
      "adduct": "+H",
      "charge": 1,
      "mse": 0.0003,
      "monoisotopic_mz": 181.070665,
      "composition": {"C": 6, "H": 13, "O": 6},
      "labels": {},
      "pair_count": 103,
      "pairs": [
        {"observed_mz": 181.070665, "observed_intensity": 100.0, "predicted_mz": 181.070665, "predicted_intensity": 100.0}
      ]
    },
    {
      "index": 0,
      "valid": true,
      "adduct": "+Na-H",
      "charge": 1,
      "mse": 185.99,
      "monoisotopic_mz": 202.044784,
      "composition": {"C": 6, "H": 11, "O": 6, "Na": 1},
      "labels": {},
      "pair_count": 108,
      "pairs": []
    },
    {
      "index": 3,
      "valid": false,
      "adduct": "-99H",
      "charge": -1,
      "mse": null,
      "monoisotopic_mz": null,
      "composition": null,
      "labels": null,
      "pair_count": null,
      "pairs": null,
      "error": "invalid_adduct",
      "message": "Adduct fragment would drive element 'H' count from 12 to -87; counts must not become negative",
      "position": 3
    }
  ]
}
```

- `candidates`：**全部候选**的结果，顺序为有效候选按 `mse` 升序、MSE 相同时按请求索引（输入顺序）稳定排序；非法候选排在末尾，彼此仍按请求索引升序。
  - 有效候选：`valid: true`，含 `mse`、`monoisotopic_mz`、`composition`、`labels`（含义同 `/api/v1/ion-envelope`）以及 `pairs` 配对明细。`pairs` 按 m/z 升序（未匹配预测峰按预测 m/z 插入）；未匹配观测峰的 `predicted_*` 为 `null`，未匹配预测峰的 `observed_*` 为 `null`。`index` 为该候选在请求中的 **0 基索引**（排序后保持不变）。
  - 非法候选：`valid: false`，评分相关字段为 `null`，并带 `error`（`invalid_adduct`）、`message` 与 `position`——即加合物解析既有的位置错误（含缺少 `+`/`-` 前导符、空白、片段内分子式错误，以及应用片段后计数变负，`position` 语义与 `/api/v1/ion-envelope` 完全一致）。非法候选**不会**导致整个请求 422。
- `best_candidate`：MSE 最小的候选（平局取输入顺序最前者），给出其请求 `index`、`adduct`、`charge`、`mse`；**全部候选均非法时为 `null`**。
- 请求级非法输入（分子式非法、峰列表为空/非递增/全零强度、`mass_tolerance` 非正、候选数不在 1～20、候选电荷为 0 或超出 ±10、峰数超过 200 等）仍按既有规则返回 422。

```bash
curl -X POST http://localhost:8000/api/v1/ion-identify \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C6H12O6","mass_tolerance":0.01,
       "peaks":[{"mz":181.0707,"intensity":100.0},{"mz":182.074,"intensity":6.5}],
       "candidates":[{"adduct":"+H","charge":1},{"adduct":"+Na-H","charge":1}]}'
```

### `POST /api/v1/ion-mixture-fit`

在 ion-identify 的基础上进一步拟合**两种离子的混合谱**：接收 2～10 个候选加合物/电荷，对每个有效候选复用 `/api/v1/ion-envelope` 的加合物解析与离子包络计算（基峰归一化为 100），枚举所有无序候选对与 0～100% 的混合比例，将两个归一化包络按比例叠加后与观测峰配对评分，返回最佳组合、全部有效组合的排名与最佳配对明细。

请求体：

```json
{
  "formula": "C2H6O",
  "peaks": [
    {"mz": 47.049141, "intensity": 100.0},
    {"mz": 48.052614, "intensity": 2.281747},
    {"mz": 68.02326, "intensity": 66.7}
  ],
  "mass_tolerance": 0.01,
  "step_percent": 1.0,
  "candidates": [
    {"adduct": "+H", "charge": 1},
    {"adduct": "+Na-H", "charge": 1},
    {"adduct": "-99H", "charge": -1}
  ]
}
```

- `formula`（必填）：中性分子式，规则同 `/api/v1/ion-identify`（可含 `[13C]` 等纯同位素标记）。
- `peaks`（必填，**1～200 个**）：观测峰，字段为 `mz`/`intensity`；m/z 必须**严格递增**，强度非负且至少一个为正。观测峰在评分前按基峰归一化为 100。
- `mass_tolerance`（必填，>0）：既用于叠加混合谱时融合两个包络间容差内的预测峰，也用于观测峰与混合谱的一对一配对（m/z 单位）。
- `candidates`（必填，**2～10 个**）：候选加合物/电荷组合，`adduct`/`charge` 语法与校验同 `/api/v1/ion-identify`。
- `step_percent`（可选，默认 1.0，范围 0.5～10）：混合比例网格步长（百分点）。第一候选比例的候选点为 0、step、2·step、…、100（100% 端点始终包含）。

拟合规则：

1. 候选对为**无序对**：按请求索引取 `(i, j)`（`i < j`），索引较小者为第一候选；第一候选比例为 r 时第二候选占 100−r。
2. 每个候选包络先各自基峰归一化为 100，再按比例加权叠加；两个包络间 m/z 差小于 `mass_tolerance` 的预测峰融合为一个（质量取强度加权平均、强度求和），混合谱随后**再次**基峰归一化为 100。比例为 0%/100% 时即为另一候选的纯包络。
3. 混合谱与（同样基峰归一化的）观测峰按既有规则一对一配对：容差内候选对按质量差绝对值从小到大贪心占用，每个峰至多参与一对；未匹配的观测峰与预测峰都以强度 0 计入 MSE（全部配对强度差平方的均值）。
4. 全部组合按 MSE 升序排名；**平局依次取较小候选索引对（先比第一索引）、较低第一候选比例**。
5. 非法候选在 `candidates` 中以 `valid: false` 返回其错误（同 ion-identify 的位置感知错误），并**跳过所有包含它的组合**；有效候选少于两个时 `combinations` 为空、`best_combination` 为 `null`。

响应：

```json
{
  "formula": "C2H6O",
  "mass_tolerance": 0.01,
  "step_percent": 1.0,
  "observed_peak_count": 15,
  "valid_candidate_count": 2,
  "combination_count": 21,
  "candidates": [
    {"index": 0, "valid": true, "adduct": "+H", "charge": 1, "monoisotopic_mz": 47.049141, "error": null, "message": null, "position": null},
    {"index": 1, "valid": true, "adduct": "+Na-H", "charge": 1, "monoisotopic_mz": 68.02326, "error": null, "message": null, "position": null},
    {"index": 2, "valid": false, "adduct": "-99H", "charge": -1, "monoisotopic_mz": null,
     "error": "invalid_adduct",
     "message": "Adduct fragment would drive element 'H' count from 6 to -93; counts must not become negative",
     "position": 3}
  ],
  "combinations": [
    {"rank": 1, "first_index": 0, "second_index": 1, "first_adduct": "+H", "first_charge": 1,
     "second_adduct": "+Na-H", "second_charge": 1,
     "first_fraction_percent": 60.0, "second_fraction_percent": 40.0, "mse": 0.0},
    {"rank": 2, "first_index": 0, "second_index": 1, "first_adduct": "+H", "first_charge": 1,
     "second_adduct": "+Na-H", "second_charge": 1,
     "first_fraction_percent": 65.0, "second_fraction_percent": 35.0, "mse": 10.963347}
  ],
  "best_combination": {
    "first_index": 0,
    "second_index": 1,
    "first_adduct": "+H",
    "first_charge": 1,
    "second_adduct": "+Na-H",
    "second_charge": 1,
    "first_fraction_percent": 60.0,
    "second_fraction_percent": 40.0,
    "mse": 0.0,
    "pair_count": 15,
    "pairs": [
      {"observed_mz": 47.049141, "observed_intensity": 100.0, "predicted_mz": 47.049141, "predicted_intensity": 100.0},
      {"observed_mz": 48.052614, "observed_intensity": 2.281747, "predicted_mz": 48.052614, "predicted_intensity": 2.281747}
    ]
  }
}
```

- `candidates`：**全部候选**按请求顺序的结果。有效候选含 `monoisotopic_mz`；非法候选 `valid: false` 并带 `error`（`invalid_adduct`）、`message`、`position`（语义与 `/api/v1/ion-identify` 完全一致），不会导致整个请求 422。
- `combinations`：**全部有效组合的排名**，按 MSE 升序（平局按候选索引对、第一候选比例升序），`rank` 从 1 开始；每项给出候选对的索引/加合物/电荷、两个候选的混合比例（百分比，和为 100）与 `mse`。组合总数为 `C(有效候选数, 2) × 比例网格点数`。
- `best_combination`：排名第一的组合，附带 `pair_count` 与 `pairs` 配对明细（按 m/z 升序；未匹配观测峰的 `predicted_*` 为 `null`，未匹配预测峰的 `observed_*` 为 `null`，两侧强度均为基峰=100 的归一化值）。**有效候选少于两个时为 `null`**。
- 请求级非法输入（分子式非法、峰列表为空/非递增/全零强度、`mass_tolerance` 非正、候选数不在 2～10、候选电荷为 0 或超出 ±10、`step_percent` 超出 0.5～10 等）仍按既有规则返回 422。

```bash
curl -X POST http://localhost:8000/api/v1/ion-mixture-fit \
  -H 'Content-Type: application/json' \
  -d '{"formula":"C2H6O","mass_tolerance":0.01,"step_percent":1.0,
       "peaks":[{"mz":47.0491,"intensity":100.0},{"mz":68.0233,"intensity":66.7}],
       "candidates":[{"adduct":"+H","charge":1},{"adduct":"+Na-H","charge":1}]}'
```

### `POST /api/v1/formula-search`

**分子式反推**：给定目标单同位素质量 `target_mass`、按质量递增的观测同位素峰，以及 **2～6 个不重复元素**各自的计数闭区间（`0 ≤ min_count ≤ max_count ≤ 200`），枚举所有单同位素质量落入 `target_mass ± mass_tolerance` 的元素组成，复用既有同位素卷积、基峰归一化（观测峰与预测峰均归一化为 100）与一对一峰配对（未匹配峰按强度 0 计入 MSE），按 MSE、绝对质量误差、元素符号字典序分子式三级排序后返回前 `max_results`（1～50）个候选。

枚举采用**分支限界**：元素按单同位素质量**从重到轻**依次分配（使界尽早收紧）；在每个树节点上，用"尚未分配元素按各自 min/max 计数还能贡献的最小/最大质量"（后缀质量界）求出当前元素计数的可行整数窗口——

```
⌈(target − tol − partial − suffix_max) / mᵢ⌉ ≤ cᵢ ≤ ⌊(target + tol − partial − suffix_min) / mᵢ⌉
```

窗口外的计数值直接计入剪枝数、不递归；窗口与 `[min_count, max_count]` 取交。全零组成（质量 0）不会成为命中。访问节点数（根节点计 1）**超过 100 万立即中止并返回 422**。

请求体：

```json
{
  "target_mass": 180.063388,
  "observed_peaks": [
    {"mass": 180.063388, "intensity": 100.0},
    {"mass": 181.066743, "intensity": 6.49},
    {"mass": 182.070098, "intensity": 0.28}
  ],
  "elements": [
    {"symbol": "C", "min_count": 0, "max_count": 12},
    {"symbol": "H", "min_count": 0, "max_count": 30},
    {"symbol": "O", "min_count": 0, "max_count": 8}
  ],
  "mass_tolerance": 0.01,
  "max_results": 10
}
```

- `target_mass`（必填，>0）：目标中性分子单同位素质量（Da），须为有限数。
- `observed_peaks`（必填，**1～200 个**）：观测峰，字段 `mass`/`intensity`，质量必须**严格递增**，强度非负且至少一个为正；评分前按基峰归一化为 100。
- `elements`（必填，**2～6 个、元素符号不重复**）：每项含 `symbol`（须为 `/api/v1/elements` 支持的元素）、`min_count`/`max_count`（均为 0～200 的整数且 `min_count ≤ max_count`）。元素在请求中的顺序不影响结果。
- `mass_tolerance`（必填，>0）：分子式单同位素质量落入 `target_mass ± mass_tolerance`（Da，闭区间）才算命中，同时也是观测/预测峰一对一配对的质量容差。
- `max_results`（必填，**1～50**）：返回候选数上限；窗口内全部命中数另由 `total_hits` 给出。

响应：

```json
{
  "target_mass": 180.063388,
  "mass_tolerance": 0.01,
  "max_results": 10,
  "observed_peak_count": 3,
  "visited_nodes": 30,
  "pruned_branches": 704,
  "total_hits": 1,
  "candidates": [
    {
      "rank": 1,
      "formula": "C6H12O6",
      "composition": {"C": 6, "H": 12, "O": 6},
      "monoisotopic_mass": 180.063388,
      "mass_error": 0.0,
      "mse": 0.01562,
      "pair_count": 103,
      "pairs": [
        {"observed_mass": 180.063388, "observed_intensity": 100.0, "predicted_mass": 180.063388, "predicted_intensity": 100.0}
      ]
    }
  ]
}
```

- `visited_nodes`：分支限界树实际访问的节点数（根计 1）；`pruned_branches`：在各节点被后缀质量界剪掉的计数值总数；`total_hits`：容差窗口内的组成总数（不受 `max_results` 截断）。
- `candidates`：按 **(MSE 升序, |质量误差| 升序, 分子式字典序)** 稳定排序后的前 `max_results` 项，`rank` 从 1 开始；三项排序键完全相同的两个组成保留枚举时**先遍历**到者。`formula` 中元素符号按**字符串字典序**排列（计数 1 省略、计数 0 的元素不出现），`composition` 仅含计数 >0 的元素；`mass_error = 预测单同位素质量 − target_mass`。
- 评分随枚举**逐条流式**进行：每到达一个叶节点命中就立即计算包络与配对，并只在一个容量为 `max_results` 的有界堆中保留当前最优候选；被挤出的候选其组成、预测峰谱与配对立即释放、不做累计，因此常驻内存只与 `max_results` 相关，与命中总数无关（全部可达命中仍会被完整遍历，以保证 `total_hits` 与计数准确）。
- `pairs`：与富集拟合接口同一套配对明细（容差内质量差全局从小到大贪心一对一，未匹配峰以对方字段 `null`、强度 0 计入 MSE），按质量升序；`pair_count` 为配对条目数。
- **没有任何组成落入窗口时正常返回 200**：`candidates` 为空列表、`total_hits` 为 0，`visited_nodes`/`pruned_branches` 照常返回。
- **访问节点超过 100 万时返回 422**（建议收窄元素计数范围或质量容差后重试）：

```json
{
  "detail": {
    "error": "search_limit_exceeded",
    "message": "Formula search visited 1000001 nodes, more than the limit of 1000000; narrow the element count ranges or the mass tolerance",
    "visited_nodes": 1000001,
    "pruned_branches": 0,
    "limit": 1000000
  }
}
```

- 其余请求级非法输入（元素不在 2～6 个、元素符号重复或未知、`min_count > max_count`、计数超出 0～200、`mass_tolerance`/`target_mass` 非正或非有限、`max_results` 不在 1～50、观测峰非递增/全零/超过 200 等）由请求校验或既有 `invalid_peaks` 规则返回 422。

```bash
curl -X POST http://localhost:8000/api/v1/formula-search \
  -H 'Content-Type: application/json' \
  -d '{"target_mass":180.063388,"mass_tolerance":0.01,"max_results":10,
       "observed_peaks":[{"mass":180.063388,"intensity":100.0},{"mass":181.066743,"intensity":6.49}],
       "elements":[{"symbol":"C","min_count":0,"max_count":12},
                   {"symbol":"H","min_count":0,"max_count":30},
                   {"symbol":"O","min_count":0,"max_count":8}]}'
```

### `POST /api/v1/mass-calibrate`

**质量轴校准**：接收 **3～200 个校准点**（每点给出观测到的单同位素 m/z `observed_mz`，以及对应参比离子的中性 `formula`、`adduct`、`charge`）、ppm 容差 `tolerance_ppm` 和一批待校正 m/z。每个校准点的**理论单同位素 m/z** 复用 `/api/v1/ion-envelope` 的现有规则：先解析分子式、施加加合物，再按电子质量规则换算

```
理论 m/z = (单同位素质量(加合后组成) − charge × 0.000548579909065 Da) / |charge|
```

随后按"理论 m/z（x）→ 观测 m/z（y）"拟合仿射直线 `y = slope·x + intercept`：

1. **稳健初模**：枚举所有理论 m/z 不相同的有效点对，每对点定义一个斜率与该直线的截距；初模斜率取**所有点对斜率的中位数**，初模截距取**对应点对截距的中位数**。
2. **内点筛选**：计算每个有效点相对初模的残差 ppm = `(观测 − 初模预测) / 理论 × 10⁶`，残差绝对值 **≤ `tolerance_ppm`（闭区间）** 的点记为内点。
3. **最小二乘**：对全部内点做一次普通最小二乘直线拟合，得到最终 `slope`/`intercept`。
4. **逆变换**：对待校正 m/z 返回 `calibrated_mz = (mz − intercept) / slope`。

响应同时返回每个校准点相对**最终模型**的残差（Da 与 ppm）以及该点的内点（`inlier: true`）/外点（`false`）标记（标记由初模筛选决定），**各点严格保持输入顺序**；分子式/加合物解析失败的点以 `valid: false` 附在对应位置（带既有 `error`/`message`），不影响其他点与整体拟合。

请求地址：`http://localhost:8000/api/v1/mass-calibrate`

请求体：

```json
{
  "calibration_points": [
    {"observed_mz": 33.033492, "formula": "CH4O", "adduct": "+H", "charge": 1},
    {"observed_mz": 47.049142, "formula": "C2H6O", "adduct": "+H", "charge": 1},
    {"observed_mz": 181.070756, "formula": "C6H12O6", "adduct": "+H", "charge": 1},
    {"observed_mz": 179.055184, "formula": "C6H12O6", "adduct": "-H", "charge": -1},
    {"observed_mz": 136.061821, "formula": "C5H5N5", "adduct": "+H", "charge": 1}
  ],
  "tolerance_ppm": 5.0,
  "mz": [181.071, 91.039]
}
```

- `calibration_points`（必填，**3～200 个**）：每项含 `observed_mz`（>0 的有限数）、`formula`（规则同 `/api/v1/ion-envelope`，可含 `[13C]` 等标记）、`adduct`（与 `/api/v1/ion-envelope` 完全相同的 `+/-[系数]片段` 语法）、`charge`（**-10～-1 或 1～10**，不允许 0）。
- `tolerance_ppm`（必填，>0）：初模选内点的绝对残差容差（ppm，闭区间）。
- `mz`（必填，**1～200 个**正值）：待校正观测 m/z 列表；也接受别名 `mzs` 或 `mz_values`。结果与输入一一对应、保持顺序。

响应：

```json
{
  "tolerance_ppm": 5.0,
  "slope": 1.0000005,
  "intercept": -0.000123,
  "calibration_point_count": 5,
  "valid_point_count": 5,
  "inlier_count": 5,
  "calibration_points": [
    {
      "index": 0,
      "valid": true,
      "observed_mz": 33.033492,
      "formula": "CH4O",
      "adduct": "+H",
      "charge": 1,
      "theoretical_mz": 33.033491,
      "residual_mz": 0.0000003,
      "residual_ppm": 0.009,
      "inlier": true,
      "error": null,
      "message": null
    }
  ],
  "calibrated_mz": [
    {"mz": 181.071, "calibrated_mz": 181.070882},
    {"mz": 91.039, "calibrated_mz": 91.038942}
  ]
}
```

- `slope` / `intercept`：内点最小二乘得到的最终仿射模型（观测 = slope × 理论 + intercept）。
- `valid_point_count`：成功解析并参与拟合的点数；`inlier_count`：初模容差内的内点数（最小二乘只用这些点）。
- `calibration_points`：**按输入顺序**返回各点；`index` 为 0 基输入索引，`theoretical_mz` 为理论单同位素 m/z，`residual_mz`/`residual_ppm` 为相对**最终模型**的残差（观测 − 模型预测；ppm 以理论 m/z 为分母），`inlier` 为初模筛选标记。解析失败的点 `valid: false`、理论/残差/标记字段为 `null`，并带 `error`（`invalid_formula` 或 `invalid_adduct`）与 `message`。
- `calibrated_mz`：与请求 `mz` 顺序一致，每项回显原 `mz` 并给出逆变换后的 `calibrated_mz`。

以下情形返回 **422**（不返回任何校准结果）：

- 有效点为空（所有点分子式/加合物均解析失败）：`error` 为 `no_valid_points`；
- 所有有效点的理论 m/z 相同、无法形成点对：`error` 为 `no_valid_point_pairs`；
- 内点少于 3 个：`error` 为 `too_few_inliers`（信息中含内点数与 `tolerance_ppm`，可放宽容差或替换异常点后重试）；
- 最小二乘斜率非正（≤ 0，逆变换无定义）：`error` 为 `non_positive_slope`。

```bash
curl -X POST http://localhost:8000/api/v1/mass-calibrate \
  -H 'Content-Type: application/json' \
  -d '{"tolerance_ppm":5.0,
       "calibration_points":[
         {"observed_mz":33.033492,"formula":"CH4O","adduct":"+H","charge":1},
         {"observed_mz":47.049142,"formula":"C2H6O","adduct":"+H","charge":1},
         {"observed_mz":181.070756,"formula":"C6H12O6","adduct":"+H","charge":1}],
       "mz":[181.071]}'
```

合法的大观测 m/z（最高可到双精度浮点上限量级）同样适用：最小二乘累加与逆变换均按二进制指数缩放后计算，中间项经过溢出区间不会产生虚假的无穷/非数，只有数学结果本身超出有限浮点范围才溢出。例如约 1.2×10¹⁵⁴ 量级的参考点：

```bash
curl -X POST http://localhost:8000/api/v1/mass-calibrate \
  -H 'Content-Type: application/json' \
  -d '{"tolerance_ppm":5.0,
       "calibration_points":[
         {"observed_mz":1.1999999999999999e154,"formula":"C1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000","adduct":"+H","charge":1},
         {"observed_mz":2.3999999999999999e154,"formula":"C2000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000","adduct":"+H","charge":1},
         {"observed_mz":3.6e154,"formula":"C3000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000","adduct":"+H","charge":1}],
       "mz":[1.1999999999999999e154]}'
```

返回 200，`slope`/`intercept`、逐点残差与 `calibrated_mz` 均为有限数（该例斜率约为 1、截距约为 0，`calibrated_mz` 约为 1.2e154）。

### `POST /api/v1/isotope-clusters`

**未知质谱同位素簇检测**（不需要分子式）：给定按 m/z 严格递增的质心峰、电荷范围 `min_charge`～`max_charge`（均在 1～10 内）、间距容差 `spacing_tolerance` 与最少成员数 `min_peaks`，在每个候选电荷 *z* 下以固定名义步长

```
step(z) = 1.003354835 / z
```

（¹³C−¹²C 中子质量差按电荷均分的 m/z 间距）枚举同位素链：

1. **候选生成**：每个峰都可作为簇的单同位素（最低 m/z）成员向后生长。相邻成员的名义间隔可以跨 **1 步或 2 步**，且实测间隔与对应整数步长之差的绝对值 **不得超过** `spacing_tolerance`（闭区间）；跨 2 步表示中间恰好缺一位同位素，**每条候选链至多包含一次** 2 步跳跃。成员数达到 `min_peaks` 的链（及其每一条合格前缀）都是候选。
2. **强度得分**：输入强度先按全谱基峰归一化（最高强度 = 100；全零谱按 0 处理）。候选得分 = 成员归一化强度之和 **− 10 × 缺失位数**。
3. **加权区间动态规划**：候选占用的资源是其成员在峰列表中的**闭索引区间** `[首成员, 末成员]`；选出索引跨度互不重叠（前一簇末索引 < 后一簇首索引）且**总得分最高**的候选集合。总分相同时依次取：成员总数更多 → 各簇单同位素 m/z 组成的序列字典序更小 → 电荷序列更小。
4. 返回每个选中簇的电荷、单同位素 m/z、得分、成员在输入峰列表中的 0 基索引、缺失位（相对单同位素成员的 0 基名义同位素位置），以及不属于任何簇的**未分配峰索引**。没有合格候选时正常返回 200，簇列表为空、所有峰均未分配。

请求体：

```json
{
  "peaks": [
    {"mz": 100.0, "intensity": 100.0},
    {"mz": 101.003355, "intensity": 60.0},
    {"mz": 102.00671, "intensity": 10.0},
    {"mz": 103.010065, "intensity": 2.0},
    {"mz": 150.0, "intensity": 1.0}
  ],
  "min_charge": 1,
  "max_charge": 10,
  "spacing_tolerance": 0.005,
  "min_peaks": 3
}
```

- `peaks`（必填，**1～200 个**）：质心峰，字段为 `mz`/`intensity`。m/z 必须**严格递增**（重复或逆序返回带 0 基 `index` 的 422）；强度必须为**非负有限数**（允许 0，允许全零谱；`NaN`/`Infinity` 由请求校验拦截）。
- `min_charge`（可选，默认 1）/ `max_charge`（可选，默认 10）：整数，均在 **1～10**，且 `min_charge ≤ max_charge`。
- `spacing_tolerance`（必填，>0 的有限数，m/z 单位）：相邻成员实测间隔与整数名义步长的最大允许偏差。
- `min_peaks`（必填，**2～200** 的整数）：候选簇至少需要的实测成员数。

响应：

```json
{
  "min_charge": 1,
  "max_charge": 10,
  "spacing_tolerance": 0.005,
  "min_peaks": 3,
  "peak_count": 5,
  "cluster_count": 1,
  "clusters": [
    {
      "charge": 1,
      "monoisotopic_mz": 100.0,
      "score": 172.0,
      "member_indices": [0, 1, 2, 3],
      "missing_positions": []
    }
  ],
  "unassigned_peaks": [4]
}
```

- `clusters`：按首成员索引升序（谱图从左到右）返回各簇。
  - `charge`：该簇的电荷状态；`monoisotopic_mz`：首成员（最低 m/z）的 m/z。
  - `score`：基峰归一化成员强度和减去 10 × 缺失位数。
  - `member_indices`：成员峰在请求峰列表中的 0 基索引，m/z 升序；不同簇的索引区间互不重叠。
  - `missing_positions`：因 2 步跳跃而缺失的名义同位素位置（相对单同位素成员的 0 基偏移），为空或仅一个位置。
- `unassigned_peaks`：未被任何选中簇使用的峰的 0 基索引。
- **无合格候选**（峰数不足、间距全部超差等）：`clusters` 为 `[]`、`cluster_count` 为 0，全部峰索引出现在 `unassigned_peaks` 中。

含缺失位的示例：峰位于名义位置 0、2、3（第一次相邻间隔跨 2 步）时，`member_indices` 为对应三峰、`missing_positions` 为 `[1]`，得分为三峰归一化强度和 − 10。

```bash
curl -X POST http://localhost:8000/api/v1/isotope-clusters \
  -H 'Content-Type: application/json' \
  -d '{"min_charge":1,"max_charge":10,"spacing_tolerance":0.005,"min_peaks":3,
       "peaks":[{"mz":100.0,"intensity":100.0},{"mz":101.003355,"intensity":60.0},
                {"mz":102.00671,"intensity":10.0},{"mz":103.010065,"intensity":2.0}]}'
```

### `POST /api/v1/consensus-spectrum`

**重复测量共识谱**：接收 **2～20 组**同一样品的重复质心峰列表（每组 **1～200 个**峰，m/z 严格递增、强度非负且至少一个为正），先把每组谱各自按基峰归一化（基峰 = 100），再把全部谱的峰合并到同一 m/z 轴上做**确定性完全链接聚类**，最后保留覆盖谱数达到 `min_replicates` 的簇并输出共识峰。

聚类规则：

1. **初始**：每个输入峰是一个单成员簇。
2. **合法合并**：两个簇仅当同时满足以下两条才可合并——
   - 合并后的簇不含同一输入谱的两个峰（每条谱在一个簇中至多贡献一个峰）；
   - 合并后簇的**端点跨度**（最高 − 最低 m/z）不超过 `max(tolerance_da, 端点均值 × tolerance_ppm / 1e6)`，端点均值为合并簇最低、最高 m/z 的算术平均。
3. **择优**：每轮在全部合法合并对中按 **(合并后跨度, 合并后最低 m/z, 成员(谱索引, 峰索引)序列)** 字典序取最优者合并，直到不存在合法合并为止（同一输入下结果完全确定）。
4. **保留与输出**：覆盖谱数 ≥ `min_replicates` 的簇输出为共识峰：
   - `mz`：成员 m/z 按**原始（请求）强度**加权的平均（成员原始强度全为 0 时退化为算术平均）；
   - `intensity`：成员**基峰归一化强度**的中位数——缺失谱（不在簇中的谱）按 0 计入，即对全部 N 条谱取中位数；
   - `members`：成员峰的 (谱索引, 峰索引)，均为 0 基。
5. 共识峰按 `mz` 升序返回。

**大强度数值稳定性**：峰强度允许取任何非负有限浮点数（最大约 1.8×10³⁰⁸）。逐谱基峰归一化先计算 `intensity / base` 比值（必在 [0, 1]）再乘 100，避免 `100 × intensity` 在除以同量级基峰之前先溢出为无穷；共识 m/z 的原始强度加权平均在直接求和会溢出时自动按簇内最大权重等比缩放后再计算，因此两组谱各含 m/z 100.000、100.001 且强度均为 1e308 时仍返回有限的 m/z 与强度，零强度峰（归一化恰好为 0）、强度跨多个数量级的权重以及全零强度簇（回退为算术平均）同样得到有限结果。聚类规则、成员索引与 6 位小数的输出精度均不受影响。

**极大 m/z 合簇边界稳定性**：合并资格的边界 `span ≤ max(tolerance_da, 端点均值 × tolerance_ppm / 1e6)` 以无溢出形式判断——端点均值按 `lo/2 + hi/2` 计算（而非会溢出的 `(lo + hi)/2`），ppm 一侧改写为比值 `span / 端点均值 ≤ tolerance_ppm / 1e6`；候选窗口先乘 ppm 比例，再与全部输入峰的 m/z 全跨度取小，保证窗口与合并判定一致（不漏掉合法合并，也不接纳超限峰，超限与否最终由合并判定把关）。因此 m/z 接近浮点上界时（例如 1e308 与 1.12e308、`tolerance_ppm=110000`：端点和 2.12e308 已超过浮点最大值）容差不再因求和溢出变成无穷大：该对跨度 1.2e307 超过边界 1.166e307，各自保持单峰簇；而 1e308 与 1.1e308（跨度 1e307，在 1.155e307 边界内）仍正常合并。完全链接、同谱峰互斥与确定性合并顺序不变。

请求体：

```json
{
  "spectra": [
    {"peaks": [{"mz": 100.0, "intensity": 50.0}, {"mz": 200.0, "intensity": 100.0}]},
    {"peaks": [{"mz": 100.004, "intensity": 80.0}, {"mz": 200.006, "intensity": 40.0}]},
    {"peaks": [{"mz": 100.002, "intensity": 30.0}]}
  ],
  "tolerance_da": 0.005,
  "tolerance_ppm": 0.0,
  "min_replicates": 2
}
```

- `spectra`（必填，**2～20 组**）：每组为 `{"peaks": [...]}`，峰字段为 `mz`/`intensity`；每组内 m/z 必须**严格递增**（重复或逆序返回带 `spectrum_index`/`peak_index` 的 422），强度为非负有限数且每组至少一个正强度（全零谱无法做基峰归一化，返回 `zero_spectrum` 422）。
- `tolerance_da`（必填，>0 的有限数，Da）/ `tolerance_ppm`（必填，≥0 的有限数）：合并簇端点跨度的上限为两者的 `max()` 组合——低质量端由 Da 容差主导，高质量端由 ppm 容差主导。
- `min_replicates`（必填，**2～20** 的整数）：簇至少覆盖的谱数；超过实际谱数时不报错，只是结果为空。

响应：

```json
{
  "spectrum_count": 3,
  "tolerance_da": 0.005,
  "tolerance_ppm": 0.0,
  "min_replicates": 2,
  "peak_count": 5,
  "cluster_count": 1,
  "clusters": [
    {
      "mz": 100.002375,
      "intensity": 100.0,
      "spectrum_count": 3,
      "members": [
        {"spectrum_index": 0, "peak_index": 0},
        {"spectrum_index": 1, "peak_index": 0},
        {"spectrum_index": 2, "peak_index": 0}
      ]
    }
  ]
}
```

- `spectrum_count`：输入谱数；`peak_count`：全部谱的输入峰总数；`cluster_count`：达到 `min_replicates` 的共识峰数。
- `clusters`：按 `mz` 升序的共识峰。`mz` 为成员原始强度加权平均 m/z；`intensity` 为缺失谱按 0 计入的归一化强度中位数；`spectrum_count` 为该簇覆盖的谱数；`members` 为成员峰的 0 基 `spectrum_index`/`peak_index`。
- 上例中 200 Da 附近的两个峰相距 0.006 Da 超出容差，各自保持单成员簇而被 `min_replicates` 过滤，故只返回 100 Da 处的共识峰。
- **没有任何簇达到 `min_replicates` 时正常返回 200**：`clusters` 为空列表、`cluster_count` 为 0。

```bash
curl -X POST http://localhost:8000/api/v1/consensus-spectrum \
  -H 'Content-Type: application/json' \
  -d '{"tolerance_da":0.005,"tolerance_ppm":5.0,"min_replicates":2,
       "spectra":[
         {"peaks":[{"mz":100.0,"intensity":50.0},{"mz":200.0,"intensity":100.0}]},
         {"peaks":[{"mz":100.004,"intensity":80.0},{"mz":200.006,"intensity":40.0}]},
         {"peaks":[{"mz":100.002,"intensity":30.0}]}]}'
```

### `POST /api/v1/enrichment-fit`

估计分子中单个同位素标记的丰度（富集度）：在 0%～100% 范围内按给定步长网格搜索候选丰度，对每个候选值复用卷积生成归一化包络，与观测峰配对比较，返回均方误差最小的丰度。

请求体：

```json
{
  "formula": "C6",
  "label": "13C",
  "peaks": [
    {"mass": 72.0, "intensity": 5.0},
    {"mass": 73.003355, "intensity": 30.0},
    {"mass": 74.00671, "intensity": 75.0},
    {"mass": 75.010065, "intensity": 100.0},
    {"mass": 76.013419, "intensity": 75.0},
    {"mass": 77.016774, "intensity": 30.0},
    {"mass": 78.020129, "intensity": 5.0}
  ],
  "step_percent": 0.5,
  "mass_tolerance": 0.05
}
```

- `formula`（必填）：分子式，规则同 `/api/v1/envelope`（可含 `[13C]` 等纯同位素标记，标记贡献固定质量位移）。
- `label`（必填）：待估计的单个同位素标签，格式为"质量数+元素符号"，如 `13C`、`15N`、`2H`。
- `peaks`（必填，1～200 个）：观测峰，质量必须**严格递增**，强度非负且至少一个为正。
- `step_percent`（可选，默认 1.0，范围 0.5～10）：丰度搜索步长（百分点）。候选点为 0、step、2·step、…、100（100% 端点始终包含）。
- `mass_tolerance`（必填，>0）：观测峰与预测峰配对的质量容差（Da）。

拟合模型：候选丰度为 p 时，被标记元素的同位素分布改为"标记同位素占 p，其余天然同位素按比例分享 1−p"（p 等于天然丰度时即还原为天然分布），再与其他元素的天然分布卷积得到预测包络；预测包络与观测峰均按各自基峰归一化为 100 后比较。

响应：

```json
{
  "formula": "C6",
  "label": "13C",
  "best_abundance_percent": 50.0,
  "mse": 0.0,
  "search_points": 201,
  "pair_count": 7,
  "pairs": [
    {"observed_mass": 72.0, "observed_intensity": 5.0, "predicted_mass": 72.0, "predicted_intensity": 5.0},
    {"observed_mass": 73.003355, "observed_intensity": 30.0, "predicted_mass": 73.003355, "predicted_intensity": 30.0},
    {"observed_mass": 74.00671, "observed_intensity": 75.0, "predicted_mass": 74.00671, "predicted_intensity": 75.0},
    {"observed_mass": 75.010065, "observed_intensity": 100.0, "predicted_mass": 75.010065, "predicted_intensity": 100.0}
  ]
}
```

- `best_abundance_percent`：最佳丰度（百分比）。MSE 相同时取较低丰度。
- `mse`：最佳丰度下全部配对的强度均方误差。
- `search_points`：搜索点数（0%～100% 按步长划分的候选数）。
- `pairs`：最佳丰度下的逐峰配对，按质量升序。每个观测峰与容差内最近的未占用预测峰配对；未匹配的观测峰 `predicted_*` 为 `null`（预测强度按 0 计），未匹配的预测峰 `observed_*` 为 `null`（观测强度按 0 计）。两侧强度均为基峰=100 的归一化值。

配对与误差细节：所有容差内的（观测, 预测）候选对按质量差从小到大贪心占用，每个峰至多参与一对；MSE 在全部配对（含未匹配峰）上计算。

### `POST /api/v1/dual-enrichment-fit`

联合估计同一分子中**两个不同元素**的目标同位素丰度：对两个标签各自在 0%～100% 范围内按各自步长网格搜索，构成二维候选网格；每个候选点组合两种富集分布、复用卷积生成归一化包络，再用与 `/api/v1/enrichment-fit` 完全相同的一对一峰配对与含未匹配峰的 MSE 评分。

请求体：

```json
{
  "formula": "C2H5NO2",
  "labels": [
    {"label": "13C", "step_percent": 1.0},
    {"label": "15N", "step_percent": 1.0}
  ],
  "peaks": [
    {"mass": 75.032028, "intensity": 42.86},
    {"mass": 76.035383, "intensity": 100.0}
  ],
  "mass_tolerance": 0.05
}
```

- `formula`（必填）：分子式，规则同 `/api/v1/enrichment-fit`。
- `labels`（必填，**恰好两个**）：两个待估计同位素标签，每个包含：
  - `label`（必填）：同位素标签，格式同单标签接口，如 `13C`、`15N`。
  - `step_percent`（可选，默认 1.0，范围 0.5～10）：该标签自己的丰度搜索步长。
- 两个标签必须属于**不同元素**（如 `13C` 与 `12C` 同属 C，属非法组合）。
- `peaks`（必填，1～200 个）与 `mass_tolerance`（必填，>0）：含义与校验规则同单标签接口。

拟合模型与单标签接口一致，只是两个被标记元素的分布同时替换为各自的富集分布；除两个被拟合元素外，其余元素的天然包络与两个候选丰度都无关，只计算一次并在整个二维网格中复用；第一个元素的富集多项式在外层丰度固定时同样只计算一次。

响应：

```json
{
  "formula": "C2H5NO2",
  "best_abundances": [
    {"label": "13C", "best_abundance_percent": 30.0},
    {"label": "15N", "best_abundance_percent": 70.0}
  ],
  "mse": 0.0,
  "search_points": 10201,
  "pair_count": 12,
  "pairs": [
    {"observed_mass": 75.032028, "observed_intensity": 42.86, "predicted_mass": 75.032028, "predicted_intensity": 42.86}
  ]
}
```

- `best_abundances`：两个标签各自的最佳丰度（百分比），顺序与请求中的 `labels` 一致。
- `search_points`：二维搜索点数，即两个标签网格点数之积（如各自 1.0 步长时为 101×101 = 10201；两个网格都始终包含 100% 端点）。
- `mse`、`pair_count`、`pairs`：含义与单标签接口相同，配对规则（容差内最近、一对一、未匹配峰按强度 0 计）也完全一致。
- 平局规则：MSE 相同时，取两项丰度按请求顺序（第一标签优先、第二标签其次）升序的字典序最小解，等价于"第一标签丰度尽量低；仍相同则第二标签丰度尽量低"。

### 错误定位

非法分子式或标记返回 HTTP 422，并指出出错字符的 **0 基位置**：

```json
{
  "detail": {
    "error": "invalid_formula",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 2
  }
}
```

可检测的错误包括：未知元素符号、元素计数为 0、未闭合的 `[`、标记格式错误、未知标记质量数、孤立数字、非法字符等。

`/api/v1/enrichment-fit` 的 422 错误同样定位明确：

```json
{
  "detail": {
    "error": "invalid_label",
    "message": "Element 'C' has no known isotope with mass number 20; check the label mass",
    "position": 0
  }
}
```

`/api/v1/dual-enrichment-fit` 中两个标签属于同一元素时返回可定位的 422（`index` 为标签列表中出错标签的 0 基索引，指向第二个标签）：

```json
{
  "detail": {
    "error": "invalid_labels",
    "message": "The two labels must target different elements, but '12C' targets element 'C' already named by '13C'",
    "index": 1
  }
}
```

- `invalid_formula`：分子式非法（`position` 为分子式中的 0 基字符位置）。
- `invalid_adduct`：`/api/v1/ion-envelope` 的加合物非法——片段缺少 `+`/`-` 前导符、系数为 0、片段为空/含空白、片段内分子式错误（`position` 为加合物字符串中的 0 基字符位置；片段内错误的位置已换算为绝对位置），或应用片段后某元素/标记计数变为负数（`position` 指向导致首次负计数的片段首字符）：

```json
{
  "detail": {
    "error": "invalid_adduct",
    "message": "Adduct fragment would drive element 'H' count from 2 to -1; counts must not become negative",
    "position": 2
  }
}
```

  注意：`/api/v1/ion-identify` 与 `/api/v1/ion-mixture-fit` 对候选加合物**不**返回上述 422——非法候选以 `valid: false` 出现在 200 响应的 `candidates` 中，并携带其请求 `index` 与同一套 `error`/`message`/`position`；仅请求级输入（分子式、峰列表、参数边界）非法时才返回 422。

- `invalid_label`：标签非法——格式错误（应为 `13C` 形式）、未知元素符号、未知质量数（`position` 为标签字符串中的 0 基字符位置）。
- `invalid_labels`：双标签组合非法——两个标签属于同一元素（`index` 为标签列表中出错标签的 0 基索引）。
- `invalid_peaks`：观测峰列表为空、质量非递增或重复（`index` 为出错峰的 0 基索引，空列表时为 `null`）。`/api/v1/isotope-clusters` 复用同一错误结构报告峰 m/z 非严格递增或重复（该接口允许全零强度，不做"至少一个正强度"校验）。`/api/v1/consensus-spectrum` 也以 `invalid_peaks` 报告某组谱的 m/z 非严格递增或重复，但定位字段为 `spectrum_index`（出错谱的 0 基索引）与 `peak_index`（谱内出错峰的 0 基索引）。
- `invalid_mass_step`：`/api/v1/neutral-profile` 的 `mass_step` 非正（Pydantic 字段校验同样会拦截 0/负数并指明 `field`）。
- `invalid_mass_range`：`/api/v1/neutral-profile` 的质量区间为空或反向（`mass_min ≥ mass_max`，含单侧给定时与自动 ±4σ 边界比较的情形）；`field` 为 `mass_min`。
- `grid_too_large`：统一网格点数超过 20000（放宽 `mass_step` 或缩窄质量区间）；`field` 为 `mass_step`。
- `no_profile_points`：区间内没有任何网格点落入任一棒峰的 ±4σ 支撑窗（窗口离谱图太远或窄于一个步长）；`field` 为 `mass_min`。
- `invalid_mass_grid`：`/api/v1/profile-centroid` 的质量网格非法——质量非严格递增，或某间距与首个间距的相对误差超过 1e-6；`field` 为 `profile`，`index` 为首个违规间距右端点的 0 基索引。
- `invalid_profile`：`/api/v1/profile-centroid` 的轮廓点数不在 3～20000、强度为负或存在非有限值（`field` 为 `profile`，索引问题附 `index`）。
- `zero_profile`：`/api/v1/profile-centroid` 的轮廓强度全部为零（至少需要一个正强度点）；`field` 为 `profile`。
- `zero_spectrum`：`/api/v1/consensus-spectrum` 的某组谱强度全部为零（基峰归一化至少需要一个正强度峰）；随附 `spectrum_index`（出错谱的 0 基索引），`peak_index` 为 `null`。
- `no_candidates`：`/api/v1/profile-centroid` 阈值过滤后不存在任何平顶局部极大候选；`field` 为 `relative_threshold`。
- `search_limit_exceeded`：`/api/v1/formula-search` 的分支限界枚举访问节点超过 100 万（随附 `visited_nodes`、`pruned_branches` 与 `limit`）；收窄元素计数范围或质量容差后重试。
- `no_valid_points` / `no_valid_point_pairs` / `too_few_inliers` / `non_positive_slope`：`/api/v1/mass-calibrate` 的拟合级失败——有效校准点为空、有效点理论 m/z 全相同无法成对点对、初模内点少于 3 个、最终最小二乘斜率非正，均返回 422 且 `detail.message` 说明原因（内点不足时含实际内点数与 `tolerance_ppm`）。注意：**单个**校准点的分子式/加合物非法不是请求失败，该点以 `valid: false`、`error` 为 `invalid_formula`/`invalid_adduct` 出现在 200 响应中。
- `too_many_candidates`：`/api/v1/isotope-clusters` 的候选枚举数超过 10 万（通常是 `spacing_tolerance` 过大使大量峰互为后继）；`field` 为 `spacing_tolerance`，应收小容差或增大 `min_peaks` 后重试。该接口其余参数级非法输入（峰数不在 1～200、强度为负或非有限、电荷超出 1～10 或范围倒置、`spacing_tolerance` 非正或非有限、`min_peaks` 超出 2～200）由请求校验返回标准 422 并指明字段位置；无合格候选不是错误，而是 200 响应中 `clusters` 为空。

  ```json
  {
    "detail": {
      "error": "grid_too_large",
      "message": "Mass grid would contain 25001 points; at most 20000 are allowed (widen mass_step or narrow the mass range)",
      "field": "mass_step"
    }
  }
  ```

- 参数越界（`step_percent` 超出 0.5～10、`mass_tolerance` 非正、`labels` 不是恰好两个、峰数超过 200、强度为负、`resolving_power`/`mass_step` 非正、`mass_min`/`mass_max` 为 0 或非有限值等）由请求校验返回标准 422 并指明字段位置。

### 其他端点

- `GET /api/v1/elements`：列出全部支持的元素及其天然同位素质量/丰度。
- `GET /health`：健康检查。

## 配置（环境变量）

均可在 `docker-compose.yml` 的 `environment` 段覆盖：

| 变量 | 默认值 | 说明 |
|---|---|---|
| `MSPC_DB_PATH` | 容器内 `/data/mspc.db`，本地 `./mspc.db` | SQLite 数据库文件路径 |
| `MSPC_MERGE_TOLERANCE` | `0.0001` | 峰合并质量容差（Da） |
| `MSPC_PRUNE_THRESHOLD` | `1e-12` | 卷积过程中按相对强度剪枝的阈值 |
| `MSPC_TOP_N_MAX` | `200` | `top_n` 上限 |
| `MSPC_TOP_N_DEFAULT` | `20` | `top_n` 默认值 |

## 本地开发（不用 Docker）

```bash
pip install -r requirements-dev.txt
uvicorn app.main:app --reload          # 服务: http://127.0.0.1:8000
pytest                                 # 运行测试
```

## 实现说明

- **多项式卷积**：每种元素的同位素模式表示为多项式 `Σ aᵢ·x^mᵢ`（丰度 × 质量的幂），分子包络为各元素多项式（按原子数取幂，使用快速幂）之积；多项式乘法即 (质量, 强度) 峰列表的卷积。纯同位素标记贡献固定质量位移。
- **峰合并**：排序后相邻质量差 < 0.0001 Da 的峰融合为一个，质量取强度加权平均，强度取和。
- **离子包络**：`POST /api/v1/ion-envelope` 将加合物（`app/adduct.py`，形如 `+2H-H2O` 的连续片段，系数乘整个片段）施加到解析后的原子计数上，用 `|charge| × MSPC_MERGE_TOLERANCE` 的中性质量容差卷积，再按 `(mass − charge·mₑ)/|charge|`（`mₑ = 0.000548579909065 Da`）逐峰换算为 m/z。
- **仪器轮廓模拟**：`POST /api/v1/neutral-profile`（`app/profile.py`）复用中性包络卷积得到全部棒峰，对每根质量 m 的棒构造 `FWHM = m/resolving_power`、`σ = FWHM/(2√(2 ln 2))` 的高斯核，只在 ±4σ 内逐点累加到同一张等间距网格（`mass_step`），区间外的核截断、区间内贡献原样保留，再将全图基峰归一化为 100。省略范围时取所有棒 ±4σ 支撑窗的并集；网格点数超过 20000、区间为空/反向、步长非正或区间内无有效采样点均返回带 `field` 的 422。
- **轮廓质心化**：`POST /api/v1/profile-centroid`（`app/centroid.py`）先校验 3～20000 点严格递增、等间距（各间距与首个间距相对误差 ≤ 1e-6）、非负且至少一个正强度；按基峰归一化后丢弃 `relative_threshold` 以下点（比较含微小浮点容差），将连续等强且两侧严格更低的游程识别为平顶局部极大（位置取平台两端质量中点），候选按归一化强度降序（同强度低质量优先）做 `min_separation` 非极大值抑制（恰等间距保留），再用二分最近邻把每个阈值以上点分配给最近保留候选（等距归低质量候选），以**原始强度**做强度加权得到质心质量、峰顶归一化强度（平台高度）与占总分配强度百分比，按质量升序返回；非法网格、全零轮廓或无候选返回带 `field`/`index` 的 422。
- **离子识别**：`POST /api/v1/ion-identify`（`app/identify.py`）对 1～20 个候选（加合物, 电荷）逐个复用同一套加合物解析、离子包络计算与基峰归一化，观测峰同样归一化一次后复用；每个候选用与富集拟合相同的 `pair_and_score` 做容差内质量差全局从小到大的一对一贪心匹配，未匹配峰以零强度计入归一化 MSE。有效候选按 MSE 升序、输入索引作稳定平局键排序；非法候选捕获既有 `AdductError`（携带位置）后以 `valid: false` 附在结果末尾，不影响其他候选评分。
- **离子混合拟合**：`POST /api/v1/ion-mixture-fit`（`app/mixture.py`）对 2～10 个候选中的每个有效候选复用 `app/identify.py` 的基峰归一化离子包络；枚举无序候选对（按请求索引升序，索引较小者为第一候选）与 0～100% 的第一候选比例网格（复用富集拟合的 `abundance_grid`，100% 端点始终包含），将两个归一化包络按比例加权叠加——两包络间 `mass_tolerance` 内的峰按强度加权质量融合、强度求和（复用 `merge_peaks`；比例为 0%/100% 时跳过另一侧，即纯包络）——混合谱再次基峰归一化后与观测峰做同一套 `pair_and_score` 一对一配对，双侧未匹配峰以零强度计入 MSE。全部组合按 (MSE, 候选索引对, 第一候选比例) 升序排名，首名即最佳组合（平局依次取较小索引对、较低首候选比例）；非法候选以 `valid: false` 返回既有位置错误并跳过其所有组合，有效候选少于两个时组合列表为空、最佳组合为 `null`。
- **质量轴校准**：`POST /api/v1/mass-calibrate`（`app/calibration.py`）对 3～200 个校准点复用 `app/parser.py`、`app/adduct.py` 与电子质量规则 `(mass − charge·mₑ)/|charge|` 计算理论单同位素 m/z，解析失败的点按既有 `invalid_formula`/`invalid_adduct` 错误在结果中原位标记并剔除。对全部有效点枚举理论 m/z 不同的点对，斜率取点对斜率中位数、截距取点对截距中位数构成稳健初模（Theil–Sen 式）；以初模残差绝对 ppm ≤ `tolerance_ppm`（闭区间，含 1e-9 ppm 浮点保护）选内点，再对内点做一次普通最小二乘拟合。有效点为空、无不同理论 m/z 的点对、内点少于 3 个或最小二乘斜率非正时抛 `CalibrationError` 转 422。响应按输入顺序返回每点的最终模型残差（Da/ppm）与内点标记，并用逆变换 `(mz − intercept)/slope` 顺序返回待校正 m/z。
- **单同位素质量**：按 IUPAC 惯例取各元素天然丰度最高的同位素质量之和。
- **富集拟合**：对每个候选丰度 p（0%～100% 按步长网格搜索），将被标记元素的同位素分布替换为"标记同位素占 p、其余天然同位素按比例分享 1−p"，复用同一套卷积生成归一化包络（基峰=100，观测峰同样归一化）；观测峰与容差内最近的预测峰贪心配对，未匹配峰强度按 0 计，以全部配对的强度均方误差最小者为最佳丰度，同分取较低丰度。除被拟合元素外其余元素的包络与候选丰度无关，只计算一次。
- **双元素联合拟合**：`POST /api/v1/dual-enrichment-fit` 在两个不同元素标签的丰度上做二维网格搜索，组合两种富集分布后沿用同一套卷积、一对一配对与含未匹配峰的 MSE。与两个候选丰度都无关的"其余元素包络"只卷积一次并在全网格复用；第一元素的富集多项式在其外层丰度固定时只计算一次。MSE 相同按两项丰度（请求顺序）升序字典序取解，`search_points` 为两网格点数之积。
- **分子式反推**：`POST /api/v1/formula-search`（`app/formula_search.py`）对 2～6 个不重复元素及各自 0～200 的计数闭区间做深度优先分支限界枚举，枚举实现为**流式生成器** `iter_compositions`（显式栈、逐叶产出 `(访问序计数, 质量)`，完整遍历结束时经生成器返回值给出节点/剪枝计数；`enumerate_compositions` 仍作为在其之上收集成列表的兼容入口保留）。元素按单同位素质量从重到轻访问；在每个节点用后缀元素最小/最大可贡献质量（`suffix_min`/`suffix_max`）把当前元素的可行计数夹到整数窗口 ⌈…⌉～⌊…⌉ 并与 `[min_count, max_count]` 取交，窗口外的计数值累加为 `pruned_branches` 且不递归（含 1e-10 Da 浮点保护，叶节点再做一次窗口校验）。生产搜索先排空一条**不保留任何命中数据、不计算包络**的廉价探测流，仅用于确定节点/剪枝计数与 100 万节点守卫（超限即在算出任何同位素包络前抛 `SearchLimitExceeded`，与旧链路"先枚举后评分"的中止时机一致）；随后再走评分流，每收到一个命中就立即用 `ParsedFormula` 复用天然丰度卷积、基峰归一化并用 `pair_and_score` 与同样归一化的观测峰一对一配对（未匹配峰按零强度计入 MSE），随后把该候选投入容量恰为 `max_results` 的有界最小堆：堆顶是当前最差候选，更优者替换堆顶并立即释放被挤出候选的组成、预测峰谱与配对，除正在评分的候选外不累计任何淘汰数据，故常驻内存不随命中数增长。堆键为 (-MSE, -|质量误差|, 逆字典序分子式, -遍历序号)，即最终按 (MSE, |质量误差|, 元素符号字典序分子式) 升序取舍、整键相同保留先遍历者；`total_hits` 返回未截断的命中总数。全零组成不命中；访问节点（根计 1）超过 100 万抛出 `SearchLimitExceeded`（携带当时的节点/剪枝计数）由专用异常处理器转 422；无命中返回 200 空候选列表。
- **未知同位素簇检测**：`POST /api/v1/isotope-clusters`（`app/isotope_clusters.py`）不依赖分子式，对每个候选电荷 z 以名义步长 `1.003354835/z` 从每个峰出发做深度优先的链生长：用二分查找在后继窗口 `[尾峰 + d·步长 ± 容差]`（d ∈ {1,2}）内取下一个成员，d=2 记录一个缺失位且全链至多一次；每条成员数 ≥ `min_peaks` 的前缀都作为候选（相同成员序列去重）。强度按全谱基峰归一化（全零谱按 0 处理），候选得分为成员强度和 − 10 × 缺失位数。选择阶段是两趟后缀加权区间动态规划：第一趟以 `(总分, 成员总数)` 标量键求出各后缀的最优值对（取候选则跳到其末成员索引 +1，否则前移一位）；第二趟在达到最优的候选中比较**完整后缀序列**（先单同位素 m/z 元组字典序、再电荷元组字典序）重建唯一最优解，避免单状态增量 DP 在前缀延长后平局关系翻转的问题。枚举候选超过 10 万抛 `IsotopeClusterError`（`too_many_candidates`）转 422。选中簇按首成员索引升序返回，未出现在任何簇中的峰索引汇总为 `unassigned_peaks`；无合格候选时簇列表为空。
- **重复测量共识谱**：`POST /api/v1/consensus-spectrum`（`app/consensus.py`）对 2～20 组质心峰先做逐谱基峰归一化（全零谱抛 `zero_spectrum` 转 422；归一化先求 `intensity/base` 比值再乘 100，使 1e308 量级的合法大强度也不会溢出），再把全部峰按 m/z 排序为单成员簇，做确定性完全链接凝聚聚类：合并对合法当且仅当两簇不共享输入谱且合并簇端点跨度 ≤ `max(tolerance_da, 端点均值 × tolerance_ppm/1e6)`；每轮用最小堆取出按 (合并跨度, 最低 m/z, 成员(谱,峰)索引) 字典序最优的合法对合并（失效堆项惰性删除，新簇只与跨度上界窗口内的簇再评估），直到无合法合并。覆盖谱数 ≥ `min_replicates` 的簇以原始强度加权 m/z（直接求和会溢出时按簇内最大权重缩放后计算，全零权重回退算术平均）、成员索引与缺失谱按 0 计入的归一化强度中位数输出，按 m/z 升序；谱 m/z 非严格递增以带 `spectrum_index`/`peak_index` 的 `invalid_peaks` 转 422。
- **同位素数据**：`app/isotopes.py` 内置约 50 种常见元素的 AME/CIAAW 同位素质量与天然丰度，可按需扩充。
- **持久化**：SQLite 存计算结果缓存（`envelope_cache`）与请求日志（`request_log`），schema 随服务启动自动创建。
