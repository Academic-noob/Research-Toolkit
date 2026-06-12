# doped 严格坐标锁缺陷结构生成：详细使用说明

## 1. 定位工作目录与文件

本说明适用于以下 Windows 工作目录：

```text
C:\Users\123\Desktop\doped_snb_work
```

打开 PowerShell，进入工作目录：

```powershell
Set-Location C:\Users\123\Desktop\doped_snb_work
```

先查看根目录：

```powershell
Get-ChildItem
```

根目录应至少包含：

```text
doped_snb_work
│
├── bulk_POSCAR
├── 00_workflow_config.strict.toml
│
└── scripts
    ├── 01_generate_doped_snb_structures_strict.py
    ├── strict_poscar_lock.py
    └── 04_verify_exact_defect_text.py
```

其中：

| 文件 | 作用 |
|---|---|
| `bulk_POSCAR` | 你提供的最终超胞，是唯一权威母版。 |
| `00_workflow_config.strict.toml` | 当前工作流的配置文件。 |
| `scripts\01_generate_doped_snb_structures_strict.py` | 缺陷结构生成入口。 |
| `scripts\strict_poscar_lock.py` | 严格坐标锁模块。 |
| `scripts\04_verify_exact_defect_text.py` | 独立文本审计脚本。 |

检查文件是否存在：

```powershell
Test-Path .\bulk_POSCAR
Test-Path .\00_workflow_config.strict.toml
Test-Path .\scripts\01_generate_doped_snb_structures_strict.py
Test-Path .\scripts\strict_poscar_lock.py
Test-Path .\scripts\04_verify_exact_defect_text.py
```

五条命令均应返回：

```text
True
```

---

## 2. 本流程的严格约束

当前流程只讨论未扰动缺陷结构生成，即 Stage 1。

`bulk_POSCAR` 是唯一权威结构。生成缺陷时：

- 不允许重新扩胞；
- 不允许改变比例因子；
- 不允许改变三行晶格矢量；
- 不允许改变坐标模式；
- 不允许重新格式化未受影响的坐标；
- 不允许改变未受影响坐标的最后一位小数；
- 不允许用 pymatgen 或 ASE 重新写出未扰动缺陷 POSCAR。

允许发生的变化仅限于描述缺陷所必需的内容。

### 2.1 空位 vacancy

例如：

```text
v_S
```

只允许：

1. 对应元素数量减 1；
2. 删除目标原子的原始坐标文本行。

其余内容必须保持不变。

### 2.2 替位 substitution

例如：

```text
Zn_Cu
```

只允许：

1. Cu 数量减 1；
2. Zn 数量加 1；
3. 被替换位点的元素归属改变。

对应坐标文本必须原样保留。

由于 POSCAR 按元素分组排列坐标，替位后该坐标行可能移动到新元素分组中，但坐标字符串本身不能变化。

### 2.3 间隙 interstitial

当前建议关闭间隙缺陷生成。

如果以后启用，只允许：

1. 对应元素数量加 1；
2. 新增一条间隙原子坐标行。

原有坐标文本必须全部保留。

---

## 3. 为什么能够保证最后一位小数不变化

普通结构处理通常会经历：

```text
POSCAR 文本
  ↓
解析为浮点数
  ↓
Structure 或 Atoms 对象
  ↓
重新格式化并写出
```

即使物理结构不变，重新写出时最后一位小数仍可能变化。

当前严格流程不这样做。

它将原始 `bulk_POSCAR` 逐行保存为字符串，只使用 doped 内部结构识别缺陷位点。定位完成后，最终输出通过原始文本行手术生成：

```text
原始 bulk_POSCAR 文本
  ↓
识别需要删除、替换或新增的位点
  ↓
只修改必要文本
  ↓
其余文本原样复制
```

内部浮点容差仅用于判断“doped 内部位点对应母版中的哪一个原子”，不会用于写出最终坐标。

---

## 4. 检查 bulk_POSCAR

### 4.1 确保使用最终超胞

`bulk_POSCAR` 必须是你已经确认的目标超胞。

不要在严格工作流中再次扩胞。

### 4.2 删除坐标区后的非结构尾部行

静态 POSCAR 应包含：

```text
标题
比例因子
三行晶格矢量
元素名称
元素数量
可选：Selective dynamics
Direct 或 Cartesian
N 行原子坐标
```

坐标区之后不要保留：

- MD 速度；
- 全零速度块；
- predictor-corrector 数据；
- 其他辅助尾部信息。

如果存在尾部全零行，可以手动删除。

删除时不要修改：

- 晶格行；
- 元素行；
- 计数行；
- 原子坐标行；
- 坐标小数位；
- 空格格式。

---

## 5. 配置文件定位与设置

打开：

```text
C:\Users\123\Desktop\doped_snb_work\00_workflow_config.strict.toml
```

### 5.1 路径设置

确认：

```toml
[paths]
bulk = "bulk_POSCAR"
workflow_root = "doped_snb_workflow"
```

`bulk = "bulk_POSCAR"` 表示读取：

```text
C:\Users\123\Desktop\doped_snb_work\bulk_POSCAR
```

### 5.2 禁止自动扩胞

确认：

```toml
[supercell]
generate_supercell = false
```

该参数必须为 `false`。

含义是：doped 只能使用你提供的超胞，不允许重新选择其他超胞。

### 5.3 严格坐标锁

确认：

```toml
[strict_lock]
enabled = true
mapping_tolerance_angstrom = 1.0e-4
lattice_transform_tolerance = 1.0e-6
```

参数含义：

| 参数 | 作用 |
|---|---|
| `enabled` | 开启严格文本锁。 |
| `mapping_tolerance_angstrom` | 判断 doped 内部位点与原始 POSCAR 位点是否对应。 |
| `lattice_transform_tolerance` | 判断 doped 内部晶格是否只存在可接受的浮点噪声。 |

这些容差只用于内部识别，不会修改最终 POSCAR。

### 5.4 关闭 ShakeNBreak

当前只生成未扰动缺陷结构，因此设置：

```toml
[shakenbreak]
enabled = false
```

不要在本阶段讨论：

- bond distortion；
- Dimer；
- rattle；
- SnB 扰动目录。

### 5.5 关闭无意义的 extrinsic 设置

如果只生成本征缺陷，设置：

```toml
[defects]
extrinsic = []
```

例如，Cu2ZnSnS4 中本身已经存在 S。生成 `v_S` 不需要：

```toml
extrinsic = ["S"]
```

### 5.6 关闭间隙缺陷

当前建议设置：

```toml
[defects]
generate_interstitials = false
```

---

## 6. 首次冒烟测试：只生成中性 S 空位

首次运行不要直接生成全部缺陷。

先设置：

```toml
[defects]
include_types = ["vacancy"]
generate_interstitials = false
extrinsic = []
include_name_regex = ["^v_S$"]
include_charges = [0]
```

这样只生成：

```text
v_S_0
```

目标是验证：

```text
S 数量减 1
删除一个 S 坐标行
其他文本完全一致
```

---

## 7. 运行缺陷结构生成

在 PowerShell 中执行：

```powershell
Set-Location C:\Users\123\Desktop\doped_snb_work

python .\scripts\01_generate_doped_snb_structures_strict.py `
    --config .\00_workflow_config.strict.toml
```

单行写法：

```powershell
python .\scripts\01_generate_doped_snb_structures_strict.py --config .\00_workflow_config.strict.toml
```

---

## 8. 输出目录定位

成功运行后，根目录下会新增：

```text
doped_snb_workflow
```

目录结构：

```text
doped_snb_workflow
│
├── 00_reports
│   ├── doped_selected_defect_inventory.csv
│   ├── strict_coordinate_lock_audit.csv
│   └── generation_summary.json
│
└── 01_doped_defects
    ├── bulk_supercell_POSCAR
    ├── bulk_supercell_POSCAR_doped_internal
    ├── defect_gen_all.json.gz
    ├── selected_defect_entries.json.gz
    │
    ├── by_charge
    │   └── <缺陷名称_电荷态>
    │       ├── POSCAR
    │       └── DEFECT_ENTRY_META.json
    │
    └── representative_defects
        └── <缺陷名称>
            ├── defect_POSCAR
            └── DEFECT_SPECIES_META.json
```

### 8.1 权威 pristine bulk

使用：

```text
doped_snb_workflow\01_doped_defects\bulk_supercell_POSCAR
```

该文件必须与原始：

```text
bulk_POSCAR
```

逐字节完全相同。

### 8.2 doped 内部结构

文件：

```text
doped_snb_workflow\01_doped_defects\bulk_supercell_POSCAR_doped_internal
```

仅用于调试。

不要将它作为 VASP pristine bulk 输入。

### 8.3 实际缺陷结构

按电荷态划分的结构位于：

```text
doped_snb_workflow\01_doped_defects\by_charge
```

例如：

```text
doped_snb_workflow\01_doped_defects\by_charge\v_S_0\POSCAR
```

---

## 9. 运行独立文本审计

生成完成后，必须运行独立审计：

```powershell
python .\scripts\04_verify_exact_defect_text.py `
    --bulk .\bulk_POSCAR `
    --workflow-root .\doped_snb_workflow
```

单行写法：

```powershell
python .\scripts\04_verify_exact_defect_text.py --bulk .\bulk_POSCAR --workflow-root .\doped_snb_workflow
```

通过时应输出：

```text
========================================================================================
Strict Stage-1 POSCAR text audit
========================================================================================
authoritative_bulk        = C:\Users\123\Desktop\doped_snb_work\bulk_POSCAR
exported_bulk             = C:\Users\123\Desktop\doped_snb_work\doped_snb_workflow\01_doped_defects\bulk_supercell_POSCAR
bulk_byte_identical       = True
defect_directories        = ...
failed_defect_directories = 0
audit                     = C:\Users\123\Desktop\doped_snb_work\doped_snb_workflow\00_reports\strict_stage1_text_audit.csv
========================================================================================
PASS: Every Stage-1 defect POSCAR satisfies strict text preservation.
```

只有出现：

```text
PASS
```

才允许继续使用这些缺陷结构。

---

## 10. 审计脚本检查内容

审计脚本不调用：

- doped；
- pymatgen；
- ASE；
- ShakeNBreak。

它只比较文本。

### 10.1 pristine bulk 检查

比较：

```text
bulk_POSCAR
```

与：

```text
doped_snb_workflow\01_doped_defects\bulk_supercell_POSCAR
```

是否逐字节完全一致。

成功标志：

```text
bulk_byte_identical = True
```

### 10.2 vacancy 检查

对空位，检查：

```text
原子数变化 = -1
恰好删除 1 条原始坐标文本
没有新增坐标文本
晶格文本完全一致
```

### 10.3 substitution 检查

对替位，检查：

```text
原子数变化 = 0
全部坐标字符串集合完全一致
晶格文本完全一致
```

### 10.4 interstitial 检查

对间隙，检查：

```text
原子数变化 = +1
全部原始坐标字符串仍然存在
只新增 1 条坐标文本
晶格文本完全一致
```

---

## 11. 二进制核对 pristine bulk

除自动审计外，还可以手动执行：

```powershell
fc.exe /b `
    .\bulk_POSCAR `
    .\doped_snb_workflow\01_doped_defects\bulk_supercell_POSCAR
```

通过时显示：

```text
FC: no differences encountered
```

这说明两个文件逐字节完全一致，包括：

- 标题；
- 空格；
- 换行；
- 晶格；
- 坐标；
- 最后一位小数。

---

## 12. 查看审计报告

文本审计报告位于：

```text
doped_snb_workflow\00_reports\strict_stage1_text_audit.csv
```

PowerShell 中查看：

```powershell
Import-Csv `
    .\doped_snb_workflow\00_reports\strict_stage1_text_audit.csv |
    Format-Table
```

严格锁内部映射报告位于：

```text
doped_snb_workflow\00_reports\strict_coordinate_lock_audit.csv
```

查看：

```powershell
Import-Csv `
    .\doped_snb_workflow\00_reports\strict_coordinate_lock_audit.csv |
    Format-Table
```

重点关注：

| 字段 | 含义 |
|---|---|
| `entry_name` | 缺陷名称与电荷态。 |
| `defect_type` | vacancy、substitution 或 interstitial。 |
| `operation` | 删除、替换或新增操作。 |
| `source_site_index` | 原始 `bulk_POSCAR` 中对应位点索引。 |
| `bulk_mapping_max_distance_angstrom` | doped 内部位点映射回原始超胞时的最大距离偏差。 |
| `lattice_transform_max_residual` | doped 内部晶格与原始晶格之间的浮点残差。 |
| `stage1_export_mode` | 应为 `LITERAL_TEXT_SURGERY`。 |

---

## 13. 从冒烟测试扩展到全部空位

确认 `v_S_0` 通过后，修改：

```toml
[defects]
include_types = ["vacancy"]
generate_interstitials = false
extrinsic = []
include_name_regex = []
include_charges = []
```

重新运行：

```powershell
python .\scripts\01_generate_doped_snb_structures_strict.py `
    --config .\00_workflow_config.strict.toml
```

然后重新审计：

```powershell
python .\scripts\04_verify_exact_defect_text.py `
    --bulk .\bulk_POSCAR `
    --workflow-root .\doped_snb_workflow
```

只有全部通过，才进入下一步。

---

## 14. 扩展到替位缺陷

全部空位通过后，再设置：

```toml
[defects]
include_types = ["vacancy", "substitution"]
generate_interstitials = false
extrinsic = []
include_name_regex = []
include_charges = []
```

重新生成并重新审计。

替位结构必须满足：

```text
原子总数不变
坐标字符串集合不变
仅元素归属和计数发生必要变化
```

---

## 15. 重新运行时的 overwrite

配置文件中通常包含：

```toml
[workflow]
overwrite = true
```

启用后，脚本会删除旧的：

```text
doped_snb_workflow
```

并重新生成。

如果旧目录中已经有需要保留的结果，先备份：

```powershell
Rename-Item `
    .\doped_snb_workflow `
    .\doped_snb_workflow_backup
```

然后再运行。

---

## 16. 常见错误与处理

### 16.1 找不到脚本

错误：

```text
[Errno 2] No such file or directory
```

确认使用：

```powershell
python .\scripts\01_generate_doped_snb_structures_strict.py `
    --config .\00_workflow_config.strict.toml
```

不要遗漏：

```text
scripts\
```

### 16.2 POSCAR 含尾部数据

错误可能包含：

```text
Authoritative POSCAR contains non-empty trailing rows
```

处理：

1. 打开 `bulk_POSCAR`；
2. 确认元素数量总和；
3. 从 `Direct` 或 `Cartesian` 后数出 N 行原子坐标；
4. 删除坐标区之后的全零速度行或辅助数据；
5. 不要修改有效坐标区。

### 16.3 晶格浮点残差过严

错误可能包含：

```text
Maximum transform residual
allowed
```

例如：

```text
Maximum transform residual = 7.134507e-07
allowed = 1.000000e-08
```

这表示内部识别阈值过严。

设置：

```toml
[strict_lock]
lattice_transform_tolerance = 1.0e-6
```

该参数只用于内部识别，不会修改最终输出。

### 16.4 extrinsic 警告

如果出现：

```text
Specified 'extrinsic' elements ['S'] are present in the host structure
```

将：

```toml
extrinsic = ["S"]
```

改为：

```toml
extrinsic = []
```

### 16.5 审计失败

如果输出：

```text
ERROR: At least one Stage-1 defect POSCAR violated strict text preservation.
```

查看：

```text
doped_snb_workflow\00_reports\strict_stage1_text_audit.csv
```

重点看：

```text
status
errors
removed_coordinate_rows
added_coordinate_rows
lattice_text_identical
```

不要跳过审计直接提交计算。

---

## 17. 每次运行后的最低检查清单

每次重新生成后，必须依次执行：

```powershell
python .\scripts\01_generate_doped_snb_structures_strict.py `
    --config .\00_workflow_config.strict.toml
```

```powershell
python .\scripts\04_verify_exact_defect_text.py `
    --bulk .\bulk_POSCAR `
    --workflow-root .\doped_snb_workflow
```

```powershell
fc.exe /b `
    .\bulk_POSCAR `
    .\doped_snb_workflow\01_doped_defects\bulk_supercell_POSCAR
```

最终必须同时满足：

```text
bulk_byte_identical = True
failed_defect_directories = 0
PASS: Every Stage-1 defect POSCAR satisfies strict text preservation.
FC: no differences encountered
```

---

## 18. 当前阶段的完成标准

Stage 1 完成的判断标准：

```text
1. bulk_POSCAR 是唯一权威母版
2. generate_supercell = false
3. shakenbreak.enabled = false
4. exported bulk 与原始 bulk 逐字节一致
5. vacancy 只删除必要坐标行
6. substitution 只改变元素归属
7. 未受影响坐标文本完全不变
8. 独立审计全部通过
```

满足以上条件后，才进入后续局域扰动结构生成讨论。
