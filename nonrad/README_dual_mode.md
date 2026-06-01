# nonrad workflow：命令行与脚本内部配置双模式

## 1. 修改原则

这组脚本保留原有的命令行传参方式，同时增加脚本内部填写参数的方式。

- `cli` 模式：保持原有用法。适合批处理、保留命令历史和重复计算。
- `internal` 模式：在 Python 脚本顶部的 `INTERNAL_CONFIG` 中填写路径和参数，然后直接运行脚本。适合学习阶段和单个体系反复测试。
- 显式写在命令行中的参数优先级最高。即使启用了 `internal` 模式，仍可临时覆盖某一个参数。

默认仍为：

```python
CONFIG_MODE = "cli"
```

需要内部填写时，修改为：

```python
CONFIG_MODE = "internal"
```

也可以不修改文件，单次运行时增加：

```bash
--mode internal
```

## 2. `00_setup_nonrad_case.sh`

该脚本原本已经包含部分内部变量，但依赖当前工作目录。新版增加了 `PROJECT_ROOT`：

```bash
PROJECT_ROOT="/home/student/works/ljd/nonrad/Int_I"
```

填写后，可以从任意目录运行：

```bash
bash /path/to/scripts/00_setup_nonrad_case.sh relax
bash /path/to/scripts/00_setup_nonrad_case.sh ccd
bash /path/to/scripts/00_setup_nonrad_case.sh static
bash /path/to/scripts/00_setup_nonrad_case.sh wswq
```

仍然保留原来的相对路径方式：不填写 `PROJECT_ROOT`，进入缺陷根目录后运行即可。

## 3. `01_generate_ccd_structures.py`

### 原有命令行方式

```bash
python scripts/01_generate_ccd_structures.py \
  --ground-contcar nonrad_case/01_relax/ground/CONTCAR \
  --excited-contcar nonrad_case/01_relax/excited/CONTCAR \
  --outdir nonrad_case/02_ccd_structures \
  --npoints 9 \
  --qmin -0.5 \
  --qmax 0.5
```

### 内部填写方式

修改脚本顶部：

```python
CONFIG_MODE = "internal"
INTERNAL_CONFIG = {
    "ground_contcar": "/absolute/path/to/nonrad_case/01_relax/ground/CONTCAR",
    "excited_contcar": "/absolute/path/to/nonrad_case/01_relax/excited/CONTCAR",
    "outdir": "/absolute/path/to/nonrad_case/02_ccd_structures",
    "npoints": 9,
    "qmin": -0.5,
    "qmax": 0.5,
}
```

然后运行：

```bash
python scripts/01_generate_ccd_structures.py
```

## 4. `02_extract_pes.py`

### 原有命令行方式

```bash
python scripts/02_extract_pes.py \
  --case-root nonrad_case \
  --dE 1.234
```

### 内部填写方式

```python
CONFIG_MODE = "internal"
INTERNAL_CONFIG = {
    "case_root": "/absolute/path/to/nonrad_case",
    "dE": 1.234,
    "ref_index": 4,
}
```

然后运行：

```bash
python scripts/02_extract_pes.py
```

## 5. `03_get_Wif_from_WSWQ.py`

### 原有命令行方式

```bash
python scripts/03_get_Wif_from_WSWQ.py \
  --case-root nonrad_case \
  --branch ground \
  --ref-index 4 \
  --def-index 192 \
  --bulk-index 189 190 191 \
  --spin 1 \
  --kpoint 1
```

### 内部填写方式

```python
CONFIG_MODE = "internal"
INTERNAL_CONFIG = {
    "case_root": "/absolute/path/to/nonrad_case",
    "branch": "ground",
    "ref_index": 4,
    "def_index": 192,
    "bulk_index": [189, 190, 191],
    "spin": 1,
    "kpoint": 1,
    "include_zero": False,
}
```

然后运行：

```bash
python scripts/03_get_Wif_from_WSWQ.py
```

## 6. `04_calculate_capture.py`

### 原有命令行方式

```bash
python scripts/04_calculate_capture.py \
  --case-root nonrad_case \
  --dQ 1.6858758484 \
  --dE 1.234 \
  --wi 0.02764 \
  --wf 0.02900 \
  --Wif 0.050 \
  --volume 1234.5 \
  --g 1.0
```

### 内部填写方式

将相同参数写入脚本顶部的 `INTERNAL_CONFIG` 后直接运行：

```bash
python scripts/04_calculate_capture.py
```

## 7. `06_scaling_corrections_vasp.py`

### 原有命令行方式

```bash
python scripts/06_scaling_corrections_vasp.py \
  --case-root nonrad_case \
  --Z -1 \
  --m-eff 0.18 \
  --eps-static 8.9 \
  --wavecar /absolute/path/to/WAVECAR \
  --band-index 189 \
  --def-index 192 \
  --tag full_scaling
```

### 内部填写方式

填写脚本顶部 `INTERNAL_CONFIG` 中所需字段，再运行：

```bash
python scripts/06_scaling_corrections_vasp.py
```

Sommerfeld 修正与 charged-supercell scaling 是独立可选的：只填写某一组参数时，另一部分会跳过。

## 8. 实际建议

不要把所有体系都永久改成内部硬编码。内部模式适合学习和单次复现；正式批量计算仍应优先使用命令行模式，因为命令本身就是可追溯记录。
