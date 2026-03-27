# Octahedral Bond Length Analyzer

A lightweight Python tool designed to extract, classify, and analyze coordination bond lengths (e.g., in perovskite structures) by automatically identifying **axial** and **equatorial** bonds based on crystal geometry.

## ✨ Key Features
- **Smart Classification**: Distinguishes between axial and equatorial bonds using 3D Periodic Boundary Conditions (PBC) and a user-defined reference vector.
- **Hierarchical Statistics**: Computes the Average, Min, Max, and Count for:
  - Specific individual bonds (e.g., `Pb238-I141(axial)`).
  - Single-center atom summaries (e.g., `Pb238-I(total)`).
  - Overall system averages (e.g., `Overall_Pb-I(axial)`).
- **One-Click Excel Export**: Automatically aggregates and sorts all data into a clean `.xlsx` report for easy plotting and publication.

## 🛠️ Requirements
Python 3.6+ is required. Install the necessary packages via terminal:
```bash
pip install numpy pandas openpyxl

# 钙钛矿/八面体键长智能统计工具 (Octahedral Bond Length Analyzer)

本脚本当初专为分析钙钛矿（如铅卤钙钛矿 $\text{APbI}_3$）及类似具有八面体配位几何的晶体结构而设计。它能够结合初始静态结构（`POSCAR`）与距离数据文件（`.dat`），自动识别、分类并统计配位键长信息，并最终生成直观的 Excel 数据报告。

---

## ✨ 核心功能亮点

1. **精准的键型识别 (轴向 vs 赤道向)**
   - 考虑三维周期性边界条件（PBC / 最小镜像约定）。
   - 用户可自定义轴向参考向量（如 c 轴 `[0,0,1]`），脚本会自动计算夹角并识别**轴向键 (axial)**与**赤道向键 (equatorial)**。
2. **多层级详细统计**
   - **单键追踪**：追踪具体原子编号间的距离演变（例：`Pb238-I141(axial)`），不再掩盖独立原子的行为。
   - **单中心汇总**：统计某个中心原子的平均配位环境（例：`Pb238-I(total)`）。
   - **全局统计**：整个体系同一类键的宏观平均值（例：`Overall_Pb-I(axial)`）。
   - 自动计算所有键的 **平均值、最小值、最大值** 和 **计数**。
3. **一键生成 Excel 报告**
   - 自动将所有统计数据结构化，按照分类、文件名排序后导出为 `.xlsx` 格式，方便后续绘图和发表。

---

## 🛠️ 环境依赖

在运行本脚本之前，请确保你的电脑上安装了 Python 3.6+ 以及以下第三方依赖库：

```bash
pip install numpy pandas openpyxl
```
*(注：`openpyxl` 是 Pandas 导出 Excel 必不可少的底层引擎。)*

---

## 📂 文件目录要求

请将所有相关文件放置在**同一个文件夹**内，并在该文件夹下运行脚本：

```text
当前文件夹/
 ├── Slect_atom_bandv2.py    # 本 Python 脚本
 ├── POSCAR                  # VASP 5+ 格式的初始结构文件（必须包含元素符号行）
 ├── data_01.dat             # 包含键长数据的源文件 (可以是多个)
 ├── data_02.dat             # (同上)
 └── ...
```

> **注意：** `.dat` 文件中的格式必须包含类似 `distance = 3.1415` 的距离数值行，以及随后的成键原子对（如 `Pb238-I141`）以供正则匹配。

---

## ⚙️ 核心参数配置指南

在使用前，请使用文本编辑器（如 VSCode 或 Notepad++）打开 Python 脚本，并根据你的具体体系修改顶部的 **核心配置区**：

```python
# ------------------------ 核心配置区 ------------------------
CENTER_ATOM = "Pb"      # 中心原子（例如钙钛矿中的 Pb 或 Sn）
COORD_ATOM = "I"        # 配位原子（例如 I, Br, Cl, O）

# 数据文件 (.dat) 的键长过滤范围（用于过滤异常数据点）
MIN_DIST = 1.5          
MAX_DIST = 5.0          

# --- 结构几何判定参数 (重点修改) ---
MAX_COORD_DIST = 3.8    # 最大配位截断半径(Å)：用来判定在 POSCAR 中哪些原子属于同一个配位八面体。
AXIAL_VECTOR =[1.0, 0.0, 0.0] # 轴向参考向量（笛卡尔坐标）。比如希望 x 轴方向的键为轴向，则填 [1,0,0]。
TOL_DEG = 30.0          # 角度容忍度：与 AXIAL_VECTOR 夹角在 0±30° 或 180±30° 内的被识别为轴向键。
# -------------------------------------------------------------
```

---

## 🚀 运行脚本

配置好参数并确认文件就绪后，打开终端/命令行，切换到该目录下执行：

```bash
python Slect_atom_bandv2.py
```

执行完毕后：
1. **控制台**会分层级输出每个 `.dat` 文件的详细统计结果。
2. **当前目录下**会生成一个名为 `Bond_Statistics_Result.xlsx` 的 Excel 文件。

---

## 📊 Excel 输出格式说明

生成的 Excel 包含以下字段（已自动整理并排序）：
- **来源文件 (Source File)**：数据来自哪一个 `.dat` 文件。
- **统计分类 (Category)**：分为三大类：
  - `1_Detailed_Specific_Bonds`: 详细具体键（带编号）
  - `2_Single_Center_Summary`: 单中心原子汇总
  - `3_Overall_Statistics`: 全局整体统计
- **键型/编号 (Bond Type)**：键的标识（如 `Pb238-I141(axial)`）。
- **平均值 (Avg / Å)**：该键的平均距离。
- **最小值 / 最大值 (Min / Max / Å)**：该键的距离极值范围。
- **计数 (Count)**：样本出现次数。

---

## ⚠️ 重要科研提示 (关于分子动力学 MD)

如果你的 `.dat` 文件提取自**高温分子动力学 (AIMD / MD) 轨迹**，请务必注意：
本脚本是基于**第一帧（静态 POSCAR）**来建立邻居列表和轴向/赤道身份映射的。
- 如果在长时间 MD 模拟中，晶格发生了**剧烈的八面体翻转**（轴向变成了赤道向），或发生了**离子扩散/键的断裂与重组**（原有的配位 I 跑了，新的 I 补进来），由于脚本使用的是静态模板，它可能会忽略新形成的键，或将翻转后的键继续视为原先的构型。
- **适用场景**：本工具非常适合处理**结构优化结果**、**低温 MD（无明显扩散）** 或 **短时间/小幅度的热振动**的数据统计分析。
