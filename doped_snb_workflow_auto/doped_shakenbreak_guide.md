# `doped` + ShakeNBreak 自动缺陷构型搜索工作流

> **目标**  
> 从一个已经弛豫的无缺陷 `bulk POSCAR` 出发，自动完成：
>
> ```text
> 超胞构建
> → 本征缺陷与可选外源缺陷自动生成
> → 电荷态自动建议与人工覆盖
> → ShakeNBreak 局部扰动结构自动生成
> → 每个缺陷 / 电荷态 / 扰动结构目录自动建立
> → 上传集群
> → 使用你自己的 POTCAR 脚本自动补齐 POTCAR
> → 自动根据 POTCAR ZVAL 与缺陷电荷态写入 NELECT
> → 批量提交 PBE 粗筛
> ```
>
> 本流程不要求在本地 Windows 安装或配置 POTCAR 资源库。  
> 本地只生成结构和目录；`POTCAR` 在集群端通过你自己的脚本补齐。

---

# 1. 文件清单

压缩包中包含：

```text
doped_snb_auto_workflow_bundle/
├── doped_shakenbreak_auto_workflow_guide_zh.md
├── 00_workflow_config.toml.example
├── 01_generate_doped_snb_structures.py
├── 02_finalize_vasp_inputs.py
├── 03_submit_all_snb_slurm.sh
└── templates/
    ├── INCAR_PBE_SNB.template
    ├── KPOINTS_Gamma
    └── job_slurm.sh.template
```

各文件作用：

| 文件 | 运行位置 | 作用 |
|---|---|---|
| `00_workflow_config.toml.example` | 本地编辑 | 控制超胞、缺陷、电荷态、ShakeNBreak 扰动和集群端 POTCAR 处理 |
| `01_generate_doped_snb_structures.py` | 本地 Windows / WSL | 从 bulk 自动生成缺陷与全部扰动结构目录，不需要 POTCAR |
| `02_finalize_vasp_inputs.py` | 集群登录节点 | 批量补齐 `POTCAR`、`NELECT`、模板文件，并生成审计报告 |
| `03_submit_all_snb_slurm.sh` | 集群登录节点 | 检查输入完整性后批量 `sbatch`，避免重复提交 |
| `templates/INCAR_PBE_SNB.template` | 本地或集群 | PBE 粗筛起始模板 |
| `templates/KPOINTS_Gamma` | 本地或集群 | Γ 点粗筛模板 |
| `templates/job_slurm.sh.template` | 本地编辑 | Slurm 模板，需按集群调整 |

---

# 3. 总体逻辑

## 3.1 本地阶段：只做结构

```text
relaxed bulk POSCAR
        ↓
doped.DefectsGenerator
        ↓
自动生成 bulk supercell
自动生成 vacancies
自动生成 substitutions / antisites
自动生成 interstitials
可选生成 extrinsic defects
自动建议 charge states
        ↓
按 TOML 过滤目标缺陷与电荷态
        ↓
ShakeNBreak Distortions.apply_distortions()
        ↓
自动建立全部目录并写入 POSCAR
```

## 3.2 集群阶段：补齐 VASP 输入

```text
02_snb_structures/**/POSCAR
        ↓
复制 INCAR / KPOINTS / job.sh 模板
        ↓
调用你自己的 POTCAR 脚本
        ↓
解析 POTCAR 中 TITEL 与 ZVAL
        ↓
自动计算 NELECT = 中性价电子数 - 缺陷电荷态
        ↓
写入 INCAR
        ↓
输出 audit CSV
        ↓
批量提交
```

---

# 4. 推荐的逐步执行策略

不要第一次就生成全部缺陷、全部电荷态、全部 ±60% 扰动。

## 4.1 第一轮：证明自动化链条可用

```toml
[defects]
include_types = ["vacancy", "substitution"]
generate_interstitials = false

[shakenbreak]
bond_distortions = [-0.2, 0.0, 0.2]
include_dimer = false
```

目标：

```text
目录生成正确
POTCAR 脚本可批量运行
NELECT 正确
模板正确
Slurm 正确
```

## 4.2 第二轮：增加间隙缺陷

```toml
generate_interstitials = true
```

## 4.3 第三轮：扩大扰动范围

```toml
bond_distortions = []
include_dimer = true
```

恢复 ShakeNBreak 默认搜索。

## 4.4 第四轮：根据目标裁剪

如果只研究特定 nonrad 缺陷：

```toml
include_name_regex = ["^I_i", "^v_I"]
include_charges = [-1, 0]
```

---

# 5. 最小命令速查

## Windows 本地

```powershell
Copy-Item .\00_workflow_config.toml.example .\00_workflow_config.toml
notepad .\00_workflow_config.toml

python .\01_generate_doped_snb_structures.py `
  --config .\00_workflow_config.toml
```

## Linux 集群

```bash
python 02_finalize_vasp_inputs.py \
  --config doped_snb_workflow/00_workflow_config.used.toml

cat doped_snb_workflow/00_reports/finalize_vasp_inputs_audit.csv

bash 03_submit_all_snb_slurm.sh \
  doped_snb_workflow/02_snb_structures
```

## ShakeNBreak 后处理

```bash
cd doped_snb_workflow/02_snb_structures
snb-parse
snb-analyse
snb-plot -cb
snb-regenerate
```

---

# 6. 官方文档

- `doped` 缺陷生成教程：  
  <https://doped.readthedocs.io/en/latest/generation_tutorial.html>

- `doped.generation.DefectsGenerator` API：  
  <https://doped.readthedocs.io/en/latest/doped.generation.html>

- `doped` 自定义超胞说明：  
  <https://doped.readthedocs.io/en/latest/advanced_analysis_tutorial.html>

- `doped` GGA 工作流：  
  <https://doped.readthedocs.io/en/latest/GGA_workflow_tutorial.html>

- ShakeNBreak Python API 教程：  
  <https://shakenbreak.readthedocs.io/en/latest/ShakeNBreak_Example_Workflow.html>

- ShakeNBreak `Distortions` API：  
  <https://shakenbreak.readthedocs.io/en/latest/shakenbreak.input.html>

---

# 7. 最终执行顺序

```text
1. 准备 relaxed bulk_POSCAR
2. 复制并编辑 00_workflow_config.toml
3. 先做 vacancy + substitution、小范围扰动、无 dimer 冒烟测试
4. 运行 01_generate_doped_snb_structures.py
5. 检查 00_reports 与代表性 POSCAR
6. 上传工作流目录到集群
7. 编辑 potcar_command 指向你的 POTCAR 生成脚本
8. 运行 02_finalize_vasp_inputs.py
9. 检查 finalize_vasp_inputs_audit.csv
10. 运行 03_submit_all_snb_slurm.sh
11. snb-parse / snb-analyse / snb-plot
12. snb-regenerate 跨电荷态复测
13. 根据需要扩大 interstitials、扰动范围与 dimer 搜索
14. 选择结构不同的候选进入 HSE
15. HSE 端点 → ΔQ → CCD → PES → Wif → 捕获系数
```
