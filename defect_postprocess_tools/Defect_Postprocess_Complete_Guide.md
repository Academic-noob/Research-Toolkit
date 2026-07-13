# FaSnI₃ 缺陷后处理工作流 v2.2.1：完整操作指南

> 适用环境：Windows PowerShell + `sci_env`  
> 项目目录：`C:\Users\123\Desktop\FaSnI3_Fh`  
> 工具目录：`C:\Users\123\Desktop\defect_postprocess_tools`  
> 当前工作流版本：v2.2.1  
> 核心原则：项目根目录保持简洁，所有自动生成结果统一进入 `_defect_postprocess`

---

# 1. 工作流目标

该工作流用于把已经完成的 VASP 缺陷计算目录整理成可供后续 Spinney/KO/FNV 缺陷形成能分析使用的数据库。

主要功能：

1. 扫描缺陷目录与电荷态；
2. 检查中性态 `NELECT`；
3. 检查各电荷态 `OUTCAR` 是否完整；
4. 将体相 `POSCAR` 与缺陷初始 `POSCAR` 进行独立结构审计；
5. 识别实际缺陷事件：
   - 空位；
   - 间隙；
   - 替位；
   - 复合缺陷；
   - FA 等有机分子缺陷；
6. 对普通原子缺陷自动生成 KO/FNV 修正中心；
7. 对 FA 等分子缺陷跳过修正中心坐标的寻找与填写；
8. 生成最终缺陷数据库；
9. 验证数据库是否与当前结构、配置和计算结果一致。

---

# 2. 最终目录结构

## 2.1 工具目录

```text
C:\Users\123\Desktop\defect_postprocess_tools\
├─ defect_config_workflow.py
├─ find_cluster_ko_center_V1.py
├─ defect_components.example.toml
├─ species_groups.example.toml
├─ README.md
└─ 其他说明文件
```

说明：

- `defect_config_workflow.py`：总控脚本；
- `find_cluster_ko_center_V1.py`：普通缺陷 KO/FNV 修正中心查找脚本；
- `*.example.toml`：仅为示例，不是当前项目的实际配置；
- 工具目录不存放项目运行结果。

---

## 2.2 项目目录

```text
C:\Users\123\Desktop\FaSnI3_Fh\
├─ defect_components.toml
├─ species_groups.toml
├─ manifest.txt
├─ spinneyv1_cluster_legacy_V2.py
│
├─ supercell\
│  ├─ POSCAR
│  ├─ OUTCAR
│  └─ vasprun.xml
│
├─ 2int_I_vac_Sn\
├─ int_I_vac_Sn\
├─ int_Sn_2_vac_Sn\
├─ int_Sn_vac_Sn\
├─ vac_Sn\
├─ vac_Sn_2vac_FA\
├─ vac_Sn_2vac_I\
├─ vac_Sn_2vac_I_2\
├─ vac_Sn_vac_FA\
├─ vac_Sn_vac_I\
├─ vac_Sn_vac_I_1\
│
└─ _defect_postprocess\
   ├─ init_report.json
   ├─ scan_report.json
   ├─ check_report.json
   ├─ structure_audit_summary.json
   ├─ generated_defect_db.py
   ├─ generated_defect_db_report.json
   ├─ generated_defect_db_validation.json
   ├─ structure_audit\
   ├─ ko_center_reports\
   ├─ logs\
   └─ backups\
```

规则：

- `defect_components.toml` 和 `species_groups.toml` 放在项目根目录；
- 不再使用 `--workspace`；
- 所有自动输出进入 `_defect_postprocess`；
- 不要把 `logs`、`ko_center_reports`、`structure_audit` 散落到项目根目录。

---

# 3. 每个缺陷目录的基本结构

本项目的缺陷初始结构放在：

```text
<缺陷目录>\charge_state_0\POSCAR
```

电荷态 SCF 输出放在：

```text
<缺陷目录>\charge_state_q\scf\OUTCAR
```

中性态 INCAR 放在：

```text
<缺陷目录>\charge_state_0\scf\INCAR
```

示例：

```text
vac_Sn\
├─ charge_state_-2\
│  └─ scf\
│     └─ OUTCAR
├─ charge_state_-1\
│  └─ scf\
│     └─ OUTCAR
└─ charge_state_0\
   ├─ POSCAR
   └─ scf\
      ├─ INCAR
      └─ OUTCAR
```

因此所有命令都应带：

```powershell
--defect-poscar-name "charge_state_0/POSCAR"
```

---

# 4. 启动 PowerShell 并设置变量

进入项目目录：

```powershell
cd C:\Users\123\Desktop\FaSnI3_Fh
```

确认当前环境：

```powershell
Get-Location
```

应显示：

```text
C:\Users\123\Desktop\FaSnI3_Fh
```

设置变量：

```powershell
$PYTHON = (Get-Command python).Source
$TOOL = "C:\Users\123\Desktop\defect_postprocess_tools\defect_config_workflow.py"
$FINDER = "C:\Users\123\Desktop\defect_postprocess_tools\find_cluster_ko_center_V1.py"
```

检查：

```powershell
$PYTHON
$TOOL
$FINDER
```

检查文件：

```powershell
Test-Path $PYTHON
Test-Path $TOOL
Test-Path $FINDER
```

都应返回：

```text
True
```

检查脚本语法：

```powershell
& $PYTHON -m py_compile $TOOL
& $PYTHON -m py_compile $FINDER
```

无输出表示语法正常。

查看帮助：

```powershell
& $PYTHON $TOOL -h
```

应至少包含：

```text
init
scan
audit
centers
check
build
validate
all
```

---

# 5. 初始化项目

运行：

```powershell
& $PYTHON $TOOL init
```

预期：

```text
initialized=C:\Users\123\Desktop\FaSnI3_Fh\_defect_postprocess
```

该命令会创建：

```text
defect_components.toml
species_groups.toml
_defect_postprocess\
```

注意：

- 若这两个 TOML 已存在，先检查内容，不要盲目覆盖；
- 项目根目录中的 TOML 是实际配置；
- 工具目录中的 example TOML 只是示例。

---

# 6. 扫描缺陷目录

运行：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    scan
```

该步骤会：

- 识别缺陷目录；
- 识别 `charge_state_*`；
- 读取中性态 `NELECT`；
- 检查 OUTCAR；
- 将发现的目录写入 `defect_components.toml`；
- 生成：
  ```text
  _defect_postprocess\scan_report.json
  ```

首次扫描后，通常需要人工填写 `components`。

打开配置：

```powershell
notepad .\defect_components.toml
```

---

# 7. `defect_components.toml` 配置规则

## 7.1 component 语法

支持形式：

```text
vac_Sn
2vac_I
int_I
2int_I
sub_I_Sn
2sub_I_Sn
```

含义：

| 写法 | 含义 |
|---|---|
| `vac_Sn` | 一个 Sn 空位 |
| `2vac_I` | 两个 I 空位 |
| `int_I` | 一个 I 间隙 |
| `2int_I` | 两个 I 间隙 |
| `sub_I_Sn` | I 替位 Sn |
| `2sub_I_Sn` | 两个 I 替位 Sn |

TOML 数组中必须加引号：

```toml
components = ["vac_Sn", "2vac_I"]
```

错误写法：

```toml
components = ["vac_Sn", 2vac_I]
```

---

## 7.2 当前项目建议配置

```toml
# FaSnI3 active defect definitions.
# correction_center_mode:
#   "auto" = 自动寻找并写入 correction_center
#   "skip" = 不寻找、不填写 correction_center，但仍生成完整缺陷条目

["2int_I_vac_Sn"]
enabled = true
components = ["sub_I_Sn", "int_I"]
label = "I_Sn + I_i"
correction_center_mode = "auto"

["int_I_vac_Sn"]
enabled = true
components = ["vac_Sn", "int_I"]
label = "V_Sn + I_i"
correction_center_mode = "auto"

["int_Sn_2_vac_Sn"]
enabled = true
components = ["int_Sn", "vac_Sn"]
label = "Sn_i + V_Sn"
correction_center_mode = "auto"

["int_Sn_vac_Sn"]
enabled = true
components = ["int_Sn", "vac_Sn"]
label = "Sn_i + V_Sn"
correction_center_mode = "auto"

["vac_Sn"]
enabled = true
components = ["vac_Sn"]
label = "V_Sn"
correction_center_mode = "auto"

["vac_Sn_2vac_FA"]
enabled = true
components = ["vac_Sn", "2vac_FA"]
label = "vac_Sn + 2vac_FA"
correction_center_mode = "skip"

["vac_Sn_2vac_I"]
enabled = true
components = ["vac_Sn", "2vac_I"]
label = "vac_Sn + 2vac_I"
correction_center_mode = "auto"

["vac_Sn_2vac_I_2"]
enabled = true
components = ["vac_Sn", "2vac_I"]
label = "vac_Sn + 2vac_I_2"
correction_center_mode = "auto"

["vac_Sn_vac_FA"]
enabled = true
components = ["vac_Sn", "vac_FA"]
label = "vac_Sn + vac_FA"
correction_center_mode = "skip"

["vac_Sn_vac_I"]
enabled = true
components = ["vac_Sn", "vac_I"]
label = "vac_Sn + vac_I"
correction_center_mode = "auto"

["vac_Sn_vac_I_1"]
enabled = true
components = ["vac_Sn", "vac_I"]
label = "vac_Sn + vac_I_1"
correction_center_mode = "auto"
```

---

# 8. `correction_center_mode = "skip"` 的准确含义

这是当前流程最重要的规则之一。

## 8.1 `skip` 不代表禁用缺陷

错误理解：

```text
skip = 不生成该缺陷
```

正确理解：

```text
skip = 该缺陷仍然完整生成，只是不寻找、不填写 correction_center
```

例如：

```toml
["vac_Sn_vac_FA"]
enabled = true
components = ["vac_Sn", "vac_FA"]
label = "vac_Sn + vac_FA"
correction_center_mode = "skip"
```

最终应生成：

```python
'vac_Sn_vac_FA': {
    'type': 'cluster_legacy',
    'label': 'vac_Sn + vac_FA',
    'components': ['vac_Sn', 'vac_FA'],
    'charges': [0, -1, -2, 1, 2],
    'conc': {
        'defect': 2,
        'electron': 615,
        'hole': 615
    }
},
```

注意：

- 条目仍然存在；
- `charges` 仍然存在；
- `conc` 仍然存在；
- 只是不包含：
  ```python
  'correction_center': ...
  ```

---

## 8.2 普通缺陷仍保留修正中心

例如：

```python
'vac_Sn_vac_I_1': {
    'type': 'cluster_legacy',
    'label': 'vac_Sn + vac_I_1',
    'components': ['vac_Sn', 'vac_I'],
    'charges': [0, -1, -2, 1, 2],
    'correction_center': {
        'kind': 'fractional_coordinates',
        'coords': [
            0.4999811865487445,
            0.6222131944325223,
            0.6666666860128694
        ]
    },
    'conc': {
        'defect': 2,
        'electron': 627,
        'hole': 627
    }
},
```

---

## 8.3 `enabled = false` 与 `skip` 的区别

| 配置 | 是否审计 | 是否生成条目 | 是否写 correction_center |
|---|---:|---:|---:|
| `enabled = true`, `auto` | 是 | 是 | 是 |
| `enabled = true`, `skip` | 是 | 是 | 否 |
| `enabled = false` | 否 | 否 | 否 |

FA 缺陷应使用：

```toml
enabled = true
correction_center_mode = "skip"
```

而不是：

```toml
enabled = false
```

---

# 9. `species_groups.toml` 配置

内容：

```toml
[FA]
composition = { C = 1, H = 5, N = 2 }
anchor_species = "C"
radius = 2.2
```

含义：

| 字段 | 含义 |
|---|---|
| `composition` | 一个 FA 分子的组成 |
| `anchor_species` | 分子分组时的锚点元素 |
| `radius` | 分子成员识别半径，单位 Å |

FA 对应：

```text
C₁H₅N₂
```

---

# 10. Windows BOM 问题

这是本项目实际遇到过的重要问题。

## 10.1 症状

```powershell
Get-Content .\species_groups.toml
```

看起来正常：

```toml
[FA]
composition = { C = 1, H = 5, N = 2 }
anchor_species = "C"
radius = 2.2
```

但 Python 读取报错：

```text
tomllib.TOMLDecodeError: Invalid statement (at line 1, column 1)
```

原因通常是文件开头存在 UTF-8 BOM：

```text
EF BB BF
```

---

## 10.2 去除 BOM

运行：

```powershell
& $PYTHON -c "from pathlib import Path; p=Path('species_groups.toml'); s=p.read_text(encoding='utf-8-sig'); p.write_text(s, encoding='utf-8', newline='\n'); print('BOM removed:', p)"
```

也建议对 `defect_components.toml` 同样执行：

```powershell
& $PYTHON -c "from pathlib import Path; p=Path('defect_components.toml'); s=p.read_text(encoding='utf-8-sig'); p.write_text(s, encoding='utf-8', newline='\n'); print('BOM removed:', p)"
```

---

## 10.3 验证 TOML

```powershell
& $PYTHON -c "import tomllib, pprint; d=tomllib.load(open('species_groups.toml','rb')); pprint.pp(d)"
```

正确输出：

```text
{'FA': {'composition': {'C': 1, 'H': 5, 'N': 2},
        'anchor_species': 'C',
        'radius': 2.2}}
```

检查缺陷配置：

```powershell
& $PYTHON -c "import tomllib, pprint; d=tomllib.load(open('defect_components.toml','rb')); pprint.pp(d)"
```

---

# 11. 独立结构审计 `audit`

## 11.1 审计做什么

比较：

```text
supercell\POSCAR
```

与：

```text
<缺陷目录>\charge_state_0\POSCAR
```

检查：

- 晶格是否一致；
- 原子数变化；
- 空位；
- 间隙；
- 替位；
- FA 分子缺失；
- 缺陷数量；
- 原子映射；
- 最小原子间距；
- 人工 components 与实际结构是否一致。

---

## 11.2 运行全部审计

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit
```

目标：

```text
selected=11 pass=11 fail=0 unresolved=0 disabled=0
```

报告位置：

```text
_defect_postprocess\structure_audit_summary.json
_defect_postprocess\structure_audit\<缺陷名>.json
```

---

## 11.3 只审计一个缺陷

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit `
    --only "vac_Sn_vac_FA"
```

---

## 11.4 查看汇总

```powershell
$r = Get-Content `
    ".\_defect_postprocess\structure_audit_summary.json" `
    -Raw | ConvertFrom-Json

$r.defects | ForEach-Object {
    [PSCustomObject]@{
        Defect   = $_.defect_name
        Status   = $_.status
        Expected = ($_.expected_components_expanded -join ", ")
        Detected = ($_.detected_components -join ", ")
    }
} | Format-Table -AutoSize
```

---

# 12. FA 审计的正确结果

## 12.1 `vac_Sn_vac_FA`

体相与缺陷组成差：

```text
C  : -1
H  : -5
N  : -2
Sn : -1
```

应识别为：

```text
vac_FA + vac_Sn
```

报告应出现：

```text
Status        : PASS
ExpectedDelta : {"C":-1,"H":-5,"N":-2,"Sn":-1}
ActualDelta   : {"C":-1,"H":-5,"N":-2,"Sn":-1}
Expected      : vac_FA, vac_Sn
Detected      : vac_FA, vac_Sn
```

---

## 12.2 `vac_Sn_2vac_FA`

体相与缺陷组成差：

```text
C  : -2
H  : -10
N  : -4
Sn : -1
```

应识别为：

```text
2vac_FA + vac_Sn
```

---

# 13. 修改 `species_groups.toml` 后必须重新审计全部缺陷

这是当前实际流程中最容易忽略的步骤。

`species_groups.toml` 被纳入审计元数据哈希。只要该文件发生变化：

- 即使只修复 BOM；
- 即使内容肉眼看起来一样；
- 即使普通原子缺陷不含 FA；

原有结构审计报告都可能被判定为 stale。

因此，修改或修复 `species_groups.toml` 后，必须运行：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit
```

不能只重新审计一个 FA 缺陷后就直接运行 `centers`。

否则会出现：

```text
selected=11 generated=0 current=0 skipped=1 unresolved=10 failed=0
```

含义：

- 只有刚重新审计的一个 FA 缺陷有效；
- 其余 10 个审计报告已过期；
- 所以它们被标记为 unresolved。

正确做法：

```text
修改 species_groups.toml
→ audit 全部
→ centers
→ check
→ build
→ validate
```

---

# 14. 生成 KO/FNV 修正中心 `centers`

## 14.1 运行

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers
```

理想结果：

```text
selected=11 generated=0 current=9 skipped=2 unresolved=0 failed=0
```

或者首次生成时：

```text
selected=11 generated=9 current=0 skipped=2 unresolved=0 failed=0
```

关键指标：

```text
skipped=2
unresolved=0
failed=0
```

两个 skipped 对应：

```text
vac_Sn_vac_FA
vac_Sn_2vac_FA
```

---

## 14.2 `centers` 的状态含义

| 状态 | 含义 |
|---|---|
| `generated` | 本次新生成中心 |
| `current` | 已有中心仍有效 |
| `skipped` | 配置为 `correction_center_mode = "skip"` |
| `unresolved` | 审计、配置或输入未满足 |
| `failed` | Finder 执行失败 |

---

## 14.3 强制重算中心

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers `
    --force
```

适用于：

- 修改了缺陷 POSCAR；
- 修改了体相 POSCAR；
- 修改了 components；
- 修改了 Finder；
- 怀疑旧中心报告不可信。

---

## 14.4 只处理一个缺陷

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers `
    --only "vac_Sn" `
    --force
```

---

# 15. 查看 KO 中心报告

例如：

```powershell
Get-Content `
    ".\_defect_postprocess\ko_center_reports\vac_Sn.json"
```

权威字段：

```json
"recommended_correction_center": {
  "kind": "fractional_coordinates",
  "coords": [
    0.4999449662,
    0.5000049806,
    0.6666667391
  ]
}
```

最终应使用：

```text
recommended_correction_center
```

不要误用：

```text
diagnostic_relaxed_mixed_centroid
```

---

# 16. 检查项目 `check`

运行：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    check
```

目标：

```text
discovered=11 generated=11 unresolved=0 disabled=0
```

报告：

```text
_defect_postprocess\check_report.json
```

---

## 16.1 查看 unresolved

```powershell
$r = Get-Content `
    ".\_defect_postprocess\check_report.json" `
    -Raw | ConvertFrom-Json

$r.defects |
Where-Object { $_.status -eq "unresolved" } |
ForEach-Object {
    [PSCustomObject]@{
        Defect     = $_.defect_name
        CenterMode = $_.correction_center_mode
        Reasons    = ($_.unresolved_reasons -join " | ")
        Warnings   = ($_.warnings -join " | ")
    }
} | Format-List
```

---

# 17. 检查 NELECT

```powershell
Get-ChildItem . -Directory |
Where-Object {
    $_.Name -notin @("supercell", "_defect_postprocess")
} |
ForEach-Object {
    $incar = Join-Path $_.FullName "charge_state_0\scf\INCAR"
    $line = if (Test-Path $incar) {
        Select-String -Path $incar -Pattern '^\s*NELECT\s*=' |
            ForEach-Object { $_.Line }
    }

    [PSCustomObject]@{
        Defect = $_.Name
        INCAR  = Test-Path $incar
        NELECT = ($line -join "; ")
    }
} | Format-Table -AutoSize
```

每个启用缺陷都必须有明确的：

```text
NELECT = ...
```

---

# 18. 检查 OUTCAR

```powershell
Get-ChildItem . -Directory |
Where-Object {
    $_.Name -notin @("supercell", "_defect_postprocess")
} |
ForEach-Object {
    Get-ChildItem $_.FullName -Directory -Filter "charge_state_*" |
    ForEach-Object {
        $outcar = Join-Path $_.FullName "scf\OUTCAR"
        $text = if (Test-Path $outcar) {
            Get-Content $outcar -Raw
        } else {
            ""
        }

        [PSCustomObject]@{
            Defect        = Split-Path $_.Parent.FullName -Leaf
            Charge        = $_.Name
            OUTCAR        = Test-Path $outcar
            TOTEN         = $text -match "free\s+energy\s+TOTEN"
            NormalEnd     = $text -match "General timing and accounting"
            CorePotential = $text -match "potential at core"
        }
    }
} | Format-Table -AutoSize
```

正常要求：

- `OUTCAR = True`
- `TOTEN = True`
- `NormalEnd = True`
- `CorePotential = True`

---

# 19. 生成数据库 `build`

只有在：

```text
unresolved=0
```

后运行：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    build
```

成功后生成：

```text
_defect_postprocess\generated_defect_db.py
```

若出现：

```text
ERROR: unresolved=...; generated_defect_db.py was not updated.
```

不要继续反复运行 build，应先查看 `check_report.json`。

---

# 20. 验证数据库 `validate`

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    validate
```

成功：

```text
generated_defect_db.py is current and valid.
```

验证报告：

```text
_defect_postprocess\generated_defect_db_validation.json
```

---

# 21. 推荐完整运行顺序

## 21.1 首次运行

```powershell
cd C:\Users\123\Desktop\FaSnI3_Fh

$PYTHON = (Get-Command python).Source
$TOOL = "C:\Users\123\Desktop\defect_postprocess_tools\defect_config_workflow.py"

& $PYTHON -m py_compile $TOOL

& $PYTHON $TOOL init

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    scan

# 人工编辑 defect_components.toml
# 人工编辑 species_groups.toml
# 去除 BOM 并验证 TOML

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    check

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    build

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    validate
```

---

## 21.2 修改 `species_groups.toml` 后

必须：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    check

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    build

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    validate
```

不能只审计一个缺陷。

---

## 21.3 修改一个缺陷 POSCAR 后

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit `
    --only "缺陷名"

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers `
    --only "缺陷名" `
    --force

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    check

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    build

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    validate
```

---

# 22. 常见问题

## 22.1 `tomllib.TOMLDecodeError` 第 1 行第 1 列

原因：

- UTF-8 BOM。

处理：

```powershell
& $PYTHON -c "from pathlib import Path; p=Path('species_groups.toml'); s=p.read_text(encoding='utf-8-sig'); p.write_text(s, encoding='utf-8', newline='\n')"
```

---

## 22.2 FA 被识别成 `vac_C + 5vac_H + 2vac_N`

原因：

- `species_groups.toml` 未成功读取；
- 常见原因是 BOM；
- 或 FA 配置缺失。

处理：

1. 验证 TOML；
2. 去除 BOM；
3. 重新 audit 全部缺陷。

---

## 22.3 `centers` 显示 skipped=1 unresolved=10

原因：

- 只重新审计了一个 FA 缺陷；
- 修改 `species_groups.toml` 后，其余审计报告 stale。

处理：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit
```

重新审计全部。

---

## 22.4 `Cannot parse component 'vac_FA'`

原因：

- Finder 不支持分子 component。

正确处理：

```toml
correction_center_mode = "skip"
```

不是把 FA 展开成八个独立原子缺陷。

---

## 22.5 `Expected exactly one unique correction center; found 0`

原因：

- KO JSON 字段未被识别；
- 或报告缺失；
- 或报告过期。

当前 v2.2.1 应识别：

```text
recommended_correction_center
```

---

## 22.6 `build` 提示 unresolved

先运行：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    check
```

再查看 unresolved 原因。

---

# 23. 最终验收标准

正式进入 Spinney 之前，应满足：

```text
audit:
pass=11 fail=0 unresolved=0

centers:
current/generated=9
skipped=2
unresolved=0
failed=0

check:
generated=11 unresolved=0 disabled=0

build:
成功生成 generated_defect_db.py

validate:
generated_defect_db.py is current and valid.
```

---

# 24. 最终数据库预期

普通原子缺陷：

```python
'vac_Sn_vac_I_1': {
    'type': 'cluster_legacy',
    'label': 'vac_Sn + vac_I_1',
    'components': ['vac_Sn', 'vac_I'],
    'charges': [0, -1, -2, 1, 2],
    'correction_center': {
        'kind': 'fractional_coordinates',
        'coords': [0.4999811865, 0.6222131944, 0.6666666860]
    },
    'conc': {
        'defect': 2,
        'electron': 627,
        'hole': 627
    }
},
```

FA 分子缺陷：

```python
'vac_Sn_vac_FA': {
    'type': 'cluster_legacy',
    'label': 'vac_Sn + vac_FA',
    'components': ['vac_Sn', 'vac_FA'],
    'charges': [0, -1, -2, 1, 2],
    'conc': {
        'defect': 2,
        'electron': 615,
        'hole': 615
    }
},
```

两者区别仅在于：

```text
FA 缺陷没有 correction_center
```

而不是：

```text
FA 缺陷不生成
```

---

# 25. 当前最简后续操作

你目前已经修复 `species_groups.toml`，下一步应直接：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    audit
```

然后依次：

```powershell
& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    centers

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    check

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    build

& $PYTHON $TOOL `
    --defect-poscar-name "charge_state_0/POSCAR" `
    validate
```
