import os
import re

# ================= 可配置参数 =================

PREFIX = "POSCAR_"        # 文件前缀筛选

DELIMITERS = r"[ _ ]"     # 分隔符（用于字段提取）

MODE = "full"            # full  = 完整文件名
                          # field = 提取字段

FIELD_INDEX = 1           # MODE="field" 时生效

OUTPUT_STYLE = "one_line"
# one_line     = 一行，空格分隔
# multi_line   = 每行一个
# quoted_comma = 加双引号，并用逗号+空格分隔

OUTPUT_FILE = "file_list.txt"

# ==============================================

files = sorted(os.listdir("."))

results = []

for f in files:

    if not os.path.isfile(f):
        continue

    if not f.startswith(PREFIX):
        continue

    name = os.path.splitext(f)[0]

    if MODE == "full":
        results.append(name)

    elif MODE == "field":
        parts = re.split(DELIMITERS, name)

        if len(parts) > FIELD_INDEX:
            results.append(parts[FIELD_INDEX])

# ========= 写入文件 =========

with open(OUTPUT_FILE, "w", encoding="utf-8") as f:

    if OUTPUT_STYLE == "one_line":
        f.write(" ".join(results))

    elif OUTPUT_STYLE == "multi_line":
        for r in results:
            f.write(r + "\n")

    elif OUTPUT_STYLE == "quoted_comma":
        f.write(", ".join(f'"{r}"' for r in results))

    else:
        raise ValueError(
            f"未知 OUTPUT_STYLE: {OUTPUT_STYLE}。"
            "可选值：one_line、multi_line、quoted_comma"
        )

print(f"已生成 {OUTPUT_FILE}")
