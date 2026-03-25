# -*- coding: utf-8 -*-

import os
import glob
import re
import numpy as np

# ------------------------ 核心配置区 ------------------------
CENTER_ATOM = "Pb"      # 中心原子
COORD_ATOM = "I"        # 配位原子

# 数据文件 (.dat) 的键长过滤范围
MIN_DIST = 1.5          
MAX_DIST = 5.0          

# --- 结构几何判定参数 (重点修改) ---
MAX_COORD_DIST = 3.8    # 最大配位键长(Å)：用于识别哪些原子属于同一个八面体，请根据实际成键情况调整
AXIAL_VECTOR =[1.0, 0.0, 0.0] # 指定轴向键的参考向量 (笛卡尔坐标系)。[0,0,1] 代表平行于 c 轴的键为轴向键
TOL_DEG = 30.0          # 角度容忍度：与 AXIAL_VECTOR 夹角在 0±30° 或 180±30° 内的被识别为轴向键
# -------------------------------------------------------------

def read_poscar(poscar_path):
    """读取 POSCAR，返回晶格基矢、元素列表和分数坐标"""
    with open(poscar_path, 'r') as f:
        lines =[l.strip() for l in f if l.strip()]
    
    # 缩放系数与晶格基矢
    scale = float(lines[1])
    cell = np.array([[float(x) for x in lines[2].split()],
        [float(x) for x in lines[3].split()],
        [float(x) for x in lines[4].split()]
    ]) * scale

    # 解析元素 (兼容 VASP 5+)
    if lines[5].split()[0].isalpha():
        element_names = lines[5].split()
        element_counts = list(map(int, lines[6].split()))
        coord_start = 8
        if lines[7].lower()[0] == 's': # Selective dynamics
            coord_start = 9
            is_direct = lines[8].lower()[0] == 'd'
        else:
            is_direct = lines[7].lower()[0] == 'd'
    else:
        raise ValueError("仅支持包含元素符号的 VASP 5+ POSCAR 格式。")

    elements =[]
    for name, count in zip(element_names, element_counts):
        elements.extend([name] * count)
    
    coords =[]
    for line in lines[coord_start:coord_start+len(elements)]:
        x, y, z = map(float, line.split()[:3])
        coords.append([x, y, z])
    coords = np.array(coords)

    # 如果是 Cartesian 坐标，统一转化为 Direct (分数坐标) 方便计算周期性边界
    if not is_direct:
        inv_cell = np.linalg.inv(cell)
        coords = np.dot(coords, inv_cell)

    return cell, elements, coords

def build_center_coord_map(cell, elements, coords, center_atom, coord_atom):
    """
    考虑周期性边界条件(PBC)，找出每个中心原子的配位原子，并区分轴向和赤道向
    """
    center_map = {}
    
    # 归一化参考轴向量
    axis_vec = np.array(AXIAL_VECTOR, dtype=float)
    axis_vec = axis_vec / np.linalg.norm(axis_vec)

    center_indices =[i for i, e in enumerate(elements) if e == center_atom]
    coord_indices =[i for i, e in enumerate(elements) if e == coord_atom]

    for c_idx in center_indices:
        c_coord = coords[c_idx]
        axial = []
        equatorial =[]
        
        for i in coord_indices:
            # 计算分数坐标差
            diff_frac = coords[i] - c_coord
            # 应用周期性边界条件 (Minimum Image Convention)
            diff_frac = diff_frac - np.round(diff_frac)
            
            # 转换为真实的笛卡尔空间向量
            diff_cart = np.dot(diff_frac, cell)
            dist = np.linalg.norm(diff_cart)
            
            # 如果距离大于配位截断半径，说明不是配位八面体内的键，直接跳过
            if dist > MAX_COORD_DIST or dist < 0.1:
                continue
                
            # 计算与指定轴向向量的夹角
            v_dir = diff_cart / dist
            cos_theta = np.clip(np.dot(v_dir, axis_vec), -1.0, 1.0)
            angle_deg = np.degrees(np.arccos(cos_theta))
            
            # 判断是否在轴向容忍角度内 (接近 0 度 或 接近 180 度)
            if angle_deg <= TOL_DEG or angle_deg >= (180 - TOL_DEG):
                axial.append(i + 1)       # VESTA 等软件通常从 1 开始编号
            else:
                equatorial.append(i + 1)
                
        center_map[c_idx + 1] = {'axial': axial, 'equatorial': equatorial}
    return center_map

def analyze_dat_file(dat_file, center_map, min_dist, max_dist):
    bond_stats = {}
    with open(dat_file, 'r') as f:
        current_distance = 0.0
        for line in f:
            line = line.strip()
            if not line:
                continue
            match = re.search(r'distance\s*=\s*(\d+\.\d+)', line)
            if match:
                current_distance = float(match.group(1))
                continue
            
            # 全局键长过滤
            if current_distance < min_dist or current_distance > max_dist:
                continue
                
            # 分析原子对
            parts = line.split()
            for part in parts:
                if '-' not in part:
                    continue
                atoms = re.findall(r'([A-Za-z]+)(\d+)', part)
                if len(atoms) != 2:
                    continue
                elem1, num1 = atoms[0]
                elem2, num2 = atoms[1]
                num1, num2 = int(num1), int(num2)

                # 遍历字典判断键型
                for center, local_bonds in center_map.items():
                    orientation = None
                    # 中心在左，配位在右
                    if center == num1 and elem2 == COORD_ATOM:
                        if num2 in local_bonds['axial']:
                            orientation = 'axial'
                        elif num2 in local_bonds['equatorial']:
                            orientation = 'equatorial'
                    # 中心在右，配位在左
                    elif center == num2 and elem1 == COORD_ATOM:
                        if num1 in local_bonds['axial']:
                            orientation = 'axial'
                        elif num1 in local_bonds['equatorial']:
                            orientation = 'equatorial'
                    
                    if orientation:
                        keys_to_update =[
                            f"{CENTER_ATOM}{center}-{COORD_ATOM}({orientation})",
                            f"{CENTER_ATOM}{center}-{COORD_ATOM}(total)",
                            f"Overall_{CENTER_ATOM}-{COORD_ATOM}({orientation})",
                            f"Overall_{CENTER_ATOM}-{COORD_ATOM}(total)"
                        ]

                        for bond_key in keys_to_update:
                            if bond_key not in bond_stats:
                                bond_stats[bond_key] = {'sum': 0.0, 'count': 0, 'min': float('inf'), 'max': float('-inf')}
                            stats = bond_stats[bond_key]
                            stats['sum'] += current_distance
                            stats['count'] += 1
                            stats['min'] = min(stats['min'], current_distance)
                            stats['max'] = max(stats['max'], current_distance)
    return bond_stats

def print_stats(bond_stats, dat_file):
    print(f"\n[{os.path.basename(dat_file)}] 统计结果：")
    if not bond_stats:
        print("没有找到匹配的键数据，请检查配位判定参数(MAX_COORD_DIST)或键长范围(MIN_DIST/MAX_DIST)。")
        return
    
    print(f"{'键型 / Bond Type':<30} | {'平均值(Å)':<10} | {'最小值(Å)':<10} | {'最大值(Å)':<10} | {'数量':<5}")
    print("=" * 80)
    
    individual_keys = sorted([k for k in bond_stats.keys() if not k.startswith("Overall_")])
    overall_keys = sorted([k for k in bond_stats.keys() if k.startswith("Overall_")])

    for k in individual_keys:
        v = bond_stats[k]
        avg = v['sum'] / v['count']
        print(f"{k:<30} | {avg:<10.5f} | {v['min']:<10.5f} | {v['max']:<10.5f} | {v['count']:<5}")

    print("-" * 80)
    print("【全局整体统计 / Overall Statistics】")
    for k in overall_keys:
        v = bond_stats[k]
        avg = v['sum'] / v['count']
        print(f"{k:<30} | {avg:<10.5f} | {v['min']:<10.5f} | {v['max']:<10.5f} | {v['count']:<5}")
    print("=" * 80)

def main():
    current_dir = os.path.dirname(os.path.abspath(__file__))
    poscar_path = os.path.join(current_dir, 'POSCAR')
    
    if not os.path.exists(poscar_path):
        print("错误：未在当前目录下找到 POSCAR 文件。")
        return

    # 解析真实的晶格信息
    cell, elements, coords = read_poscar(poscar_path)
    
    # 建立映射（包含了严谨的 PBC 和角度判定）
    center_map = build_center_coord_map(cell, elements, coords, CENTER_ATOM, COORD_ATOM)

    dat_files = glob.glob(os.path.join(current_dir, '*.dat'))
    if not dat_files:
        print("当前文件夹下没有找到任何 .dat 文件")
        return
    
    for f in dat_files:
        stats = analyze_dat_file(f, center_map, MIN_DIST, MAX_DIST)
        print_stats(stats, f)

if __name__ == '__main__':
    main()