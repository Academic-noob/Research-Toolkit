#!/bin/bash

taskname='NEB_i_v' #这是你的化学物质的标签。这个脚本它会自动修改你的job文件，逻辑是任务名称+缺陷文件名称+始模态
BASE_DIR=$(pwd)
#如果某个缺陷文件夹里面fin、ini出现问题或者是需要单独修改INCAR/KPOINTS参数，要么把该缺陷文件夹删除重新使用该脚本提交，要么进入fin或者是ini手动提交。

PAIRED_INI=false  # false: 所有缺陷共用ini；true: 每个缺陷独立ini。
#两种模式的目录结构
# neb/
# ├── ini/ (Master IS)
# ├── Vac_I_1_39_4/fin/
# └── Vac_I_1_39_22/fin/
# neb/
# ├── Vac_I_1_39_4/
# │   ├── ini/
# │   └── fin/
# └── Vac_I_1_39_22/
    # ├── ini/
    # └── fin/
# POSCAR 文件名列表（去掉前缀 "POSCAR_"）；补充后继续计算记得删除正在进行计算的文件名；每个缺陷独立ini的话，同一个缺陷的文件名称只能是后缀不一样。
DEFECT_LIST=(
"Vac_I_1_39_4_FS"
"Vac_I_1_39_22_FS"
"Vac_I_1_39_29_FS"
"Vac_I_1_39_38_FS"
)

NEB_DIR="neb" # 存在则不会新建，对于不同的共用IS，最好分开放置文件夹，不然分不清。
mkdir -p "$NEB_DIR" || exit 1
cd "$NEB_DIR" || exit 1

if [ "$PAIRED_INI" = false ]; then
    MASTER_IS="Vac_I_1_39" # 对于共用时记得修改IS名称，目的是以共用IS名称作为分类标签，防止覆盖。
    INI_DIR="./ini"

    if [ ! -d "$INI_DIR" ]; then
        mkdir -p "$INI_DIR"
        cp "$BASE_DIR/POSCAR_${MASTER_IS}_IS" "$INI_DIR/POSCAR" || exit 1

        (
            cd "$INI_DIR" || exit 1
            python /home/student/bin/pbe5.4.4_pot/potcar_5.4.4.py &&
            cp "$BASE_DIR/KPOINTS" ./KPOINTS &&
            cp "$BASE_DIR/INCAR" ./INCAR &&
            cp "$BASE_DIR/job.sh" ./job.sh &&
            sed -i "2c #SBATCH --job-name=${taskname}_ini_${MASTER_IS}" job.sh &&
            sbatch job.sh
        )
    fi
fi

for DEFECT in "${DEFECT_LIST[@]}"; do
    DEFECT_SHORT=${DEFECT%_*}  
    DEFECT_DIR="$DEFECT_SHORT"

    if [ ! -d "$DEFECT_DIR" ]; then
        FS_DIR="$DEFECT_DIR/fin"
        mkdir -p "$FS_DIR" || exit 1
        cp "$BASE_DIR/POSCAR_$DEFECT" "$FS_DIR/POSCAR" || exit 1

        (
            cd "$FS_DIR" || exit 1
            python /home/student/bin/pbe5.4.4_pot/potcar_5.4.4.py &&
            cp "$BASE_DIR/KPOINTS" ./KPOINTS &&
            cp "$BASE_DIR/INCAR" ./INCAR &&
            cp "$BASE_DIR/job.sh" ./job.sh &&
            sed -i "2c #SBATCH --job-name=${taskname}_${DEFECT_SHORT}_fin" job.sh &&
            sbatch job.sh
        )

        if [ "$PAIRED_INI" = true ]; then
            INI_DIR="$DEFECT_DIR/ini"
            IS_FILE="$BASE_DIR/POSCAR_${DEFECT_SHORT}_IS"
            if [ -f "$IS_FILE" ]; then
                mkdir -p "$INI_DIR"
                cp "$IS_FILE" "$INI_DIR/POSCAR" || exit 1

                (
                    cd "$INI_DIR" || exit 1
                    python /home/student/bin/pbe5.4.4_pot/potcar_5.4.4.py &&
                    cp "$BASE_DIR/KPOINTS" ./KPOINTS &&
                    cp "$BASE_DIR/INCAR" ./INCAR &&
                    cp "$BASE_DIR/job.sh" ./job.sh &&
                    sed -i "2c #SBATCH --job-name=${taskname}_ini_${DEFECT_SHORT}" job.sh &&
                    sbatch job.sh
                )
            else
                echo "Warning: IS file $IS_FILE not found! Skipping ini for $DEFECT_SHORT."
            fi
        fi
    else
        echo "Skip $DEFECT_DIR: already exists."
    fi
done