# Research-Toolkit
# neb_submit,sh
A robust Bash script designed to automate the tedious process of NEB (Nudged Elastic Band) pre-processing.

Smart Logic: Supports both shared and independent initial states (IS).

Automation: Automatically links POSCARs, generates POTCARs, and updates Slurm job.sh headers.

Efficiency: Batch submits multiple defect migration paths with a single execution.
# 🚀 VASP NEB 前处理自动化提交流程脚本使用说明书

## 1. 脚本简介

本脚本专为 VASP 计算中的 **NEB（微动弹性带）过渡态搜索前处理** 编写。
其核心旨在解决 NEB 计算前期**极度繁琐且容易出错**的建文件夹、改名、配对、改 Slurm 任务名及作业提交工作。

**核心特性：**
* **高度自动化**：一键生成标准化的 `ini`（初始态）和 `fin`（末态）目录树。
* **双模式切换**：支持“多末态共用单初始态”与“独立初始态一对一配对”两种物理场景。
* **防呆防覆盖**：内置目录检测机制，已存在或正在计算的任务会被自动跳过，绝不覆盖计算数据。
* **动态任务追踪**：自动修改 `job.sh` 中的 `#SBATCH --job-name`，让您在 `squeue` 中清晰分辨每一个缺陷任务。

---

## 2. 运行前环境与文件准备

在执行本脚本前，请务必确保**当前工作目录**下存在以下文件，且命名绝对规范：

1. **公共参数文件**：
   * `INCAR`、`KPOINTS`、`job.sh`（Slurm 作业模板）。
2. **结构文件（POSCAR）**：
   * 必须严格以 `POSCAR_` 作为前缀。
   * 初始态必须以 `_IS` 结尾（如：`POSCAR_Vac_I_1_39_IS`）。
   * 末态必须以 `_FS` 结尾（如：`POSCAR_Vac_I_1_39_4_FS`）。
3. **依赖环境**：
   * 脚本依赖预设的 POTCAR 生成工具（`/home/student/bin/pbe5.4.4_pot/potcar_5.4.4.py`），请确保当前环境已加载 Python，且该路径有效。

---

## 3. 核心计算模式详解 (💡 重点)

NEB 计算中，初始态（IS）的选取对应着不同的物理场景。通过修改脚本中的 `PAIRED_INI` 变量，可以在以下两种模式间无缝切换。

### 模式 A：共享初始态模式 (`PAIRED_INI=false`)
* **物理场景**：体系的初始态是完全相同的（例如：一个完美的超胞，或者同一个基态缺陷结构），但原子可能向周边多个不同位点发生跃迁，产生多个不同的末态。
* **工作逻辑**：脚本会在根目录下**只创建一个共享的 `ini` 文件夹**，并提交一次计算。随后，为每一个指定的末态单独创建文件夹和 `fin` 子目录。
* **所需变量**：必须正确填写 `MASTER_IS`。
* **生成的目录拓扑**：

```text
neb/
├── ini/                    <-- (提取 MASTER_IS 的结构，全批次共用)
├── Vac_I_1_39_4/
│   └── fin/                <-- (缺陷 4 的末态任务)
└── Vac_I_1_39_22/
    └── fin/                <-- (缺陷 22 的末态任务)
```
	
### 模式 B：独立初始态配对模式 (PAIRED_INI=true)
* **物理场景**：每一个末态，都有一个与之在构型、原子排序上严格对应的专属初始态。必须做到“一对一”计算。
* **工作逻辑**：脚本会为指定的每一个缺陷创建独立的主文件夹。在主文件夹内，同时生成独立的 ini 和 fin 目录。脚本会根据末态名字，自动去当前目录下寻找同名的 _IS 结构文件并塞入 ini 中。
* **所需变量**：MASTER_IS 变量失效，无需理会。但要求同组的 IS 和 FS 命名严格对应（除后缀外完全一致）。
* **生成的目录拓扑**：
```text
neb/
├── Vac_I_1_39_4/
│   ├── ini/                <-- (自动寻找并读取 POSCAR_Vac_I_1_39_4_IS)
│   └── fin/                <-- (读取 POSCAR_Vac_I_1_39_4_FS)
└── Vac_I_1_39_22/
    ├── ini/                <-- (自动寻找并读取 POSCAR_Vac_I_1_39_22_IS)
    └── fin/                <-- (读取 POSCAR_Vac_I_1_39_22_FS)
```

## 4. 变量修改指南 (每次计算必看)
打开脚本，您只需修改前 30 行中的以下几个变量：

### 1. 作业名称前缀：显示在超算队列中，建议每次换体系都修改
taskname='NEB_i_v' 

### 2. 模式切换：false 为共享模式，true 为独立配对模式
PAIRED_INI=false  

### 3. 输出总目录：建议每次换体系修改名字，如 "neb_Vac_O"
NEB_DIR="neb"     

### 4. 共用初始态名称（仅当 PAIRED_INI=false 时有效！）
规则：去掉 "POSCAR_" 和 "_IS"
MASTER_IS="Vac_I_1_39" 

### 5. 待计算的末态列表清单
规则：去掉 "POSCAR_" 前缀，必须保留 "_FS" 后缀，用双引号包裹
DEFECT_LIST=(
"Vac_I_1_39_4_FS"
"Vac_I_1_39_22_FS"
"Vac_I_1_39_29_FS"
)

## 5. ⚠️ 深度排错与注意事项 (极重要)

### 🔴 1. 增量追加任务与防覆盖
如果您昨天提交了 4 个缺陷，今天想追加计算第 5 个缺陷：
* **正确做法**：直接把第 5 个缺陷的名字加入 `DEFECT_LIST` 数组中，然后再次运行脚本。脚本会打印 `Skip xxx: already exists.` 跳过前 4 个，**安全且精准地只生成并提交第 5 个**。
* **切勿担心**：脚本绝不会覆盖已存在文件夹中的 `CHGCAR`、`WAVECAR` 等重要数据。

### 🔴 2. 任务报错后的干预策略
如果超算因为断电、步数耗尽导致某个缺陷的 `fin` 或 `ini` 计算失败，或者您需要**单独调控某个任务的参数**：
* **方案 A（推倒重来）**：在根目录下直接删除该缺陷对应的整个文件夹（例如 `rm -rf neb/Vac_I_1_39_4`），修改好参数后，重新执行本脚本。
* **方案 B（精准微调）**：进入报错的 `fin` 或 `ini` 文件夹，手动修改 `INCAR`，然后直接在该目录下输入 `sbatch job.sh`。**请勿在父级目录再次运行本自动化脚本**。

### 🔴 3. 数组书写规范陷阱
在修改 `DEFECT_LIST` 时，请务必注意 Bash 的语法：
* 每行一个文件名，**必须带双引号** `""`。
* 严禁在文件名内部或双引号外部包含不必要的空格，不要添加逗号 `,`。

### 🔴 4. Windows 换行符引发的血案 (`\r` 报错)
如果您习惯在 Windows 系统中编辑此脚本，然后再上传到超算服务器运行，极大可能会遇到如下报错：
> `bash: syntax error near unexpected token` 或 `\r: command not found`

**终极解决方案**：在超算终端中执行一次以下命令，剔除 Windows 换行符即可：
```bash
sed -i 's/\r$//' 您的脚本文件名.sh
```
