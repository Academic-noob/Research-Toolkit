F_MAX_CONVERGENCE = 0.02    # 最大力收敛阈值 (eV/Å)
DRIFT_CONVERGENCE = 0.005   # 总漂移收敛阈值 (eV/Å)
ENERGY_CONVERGENCE = None    # 能量阈值，可不显示
PRESSURE_CONVERGENCE = None  # 压力阈值，可不显示
NELM_CONVERGENCE = None      # SCF步阈值，可不显示
# ============================================================

import sys
import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import re
import os

class VaspDoctor:
    def __init__(self, outcar_file="OUTCAR", oszicar_file="OSZICAR", poscar_file="POSCAR",
                 F_max_thresh=F_MAX_CONVERGENCE, Drift_thresh=DRIFT_CONVERGENCE):
        self.outcar_path = outcar_file
        self.oszicar_path = oszicar_file
        self.poscar_path = poscar_file

        self.steps = []
        self.nelm_list = []
        self.atom_names = []
        self.atom_counts = []
        self.full_atom_list = []
        self.poscar_lattice = []

        self.F_max_thresh = F_max_thresh
        self.Drift_thresh = Drift_thresh
        self.Energy_thresh = ENERGY_CONVERGENCE
        self.Pressure_thresh = PRESSURE_CONVERGENCE
        self.NELM_thresh = NELM_CONVERGENCE

    # ---------------- 1. 解析 POSCAR ----------------
    def read_poscar(self):
        if not os.path.exists(self.poscar_path):
            print(f" 提示: 找不到 {self.poscar_path}，原子名称将显示为 '?'")
            return
        with open(self.poscar_path,'r') as f:
            lines = f.readlines()
        try:
            self.poscar_lattice = [list(map(float, lines[i].split())) for i in range(2,5)]
            line5_parts = lines[5].strip().split()
            line6_parts = lines[6].strip().split()
            if line5_parts[0].isdigit():
                self.atom_counts = list(map(int,line5_parts))
                self.atom_names = ["X"]*len(self.atom_counts)
            elif line6_parts[0].isdigit():
                self.atom_names = line5_parts
                self.atom_counts = list(map(int,line6_parts))
            self.full_atom_list = []
            for name,count in zip(self.atom_names,self.atom_counts):
                self.full_atom_list.extend([name]*count)
            print(f" POSCAR 信息: 共 {sum(self.atom_counts)} 个原子 {self.atom_names}")
        except Exception as e:
            print(f"️ POSCAR 解析警告: {e}")

    # ---------------- 2. 解析 OSZICAR ----------------
    def parse_oszicar(self):
        if not os.path.exists(self.oszicar_path):
            print(" 错误: 找不到 OSZICAR，无法绘制电子步图")
            return
        with open(self.oszicar_path,'r') as f:
            lines = f.readlines()
        nelm_counts=[]
        current_nelm=0
        for line in lines:
            if re.search(r'(DAV:|RMM:|Broyden:|CG:)',line):
                current_nelm+=1
            elif "F=" in line:
                if current_nelm==0: current_nelm=1
                nelm_counts.append(current_nelm)
                current_nelm=0
        self.nelm_list=nelm_counts
        print(f" OSZICAR 解析: {len(self.nelm_list)} 个离子步的电子步数据")

    # ---------------- 3. 解析 OUTCAR (核心修复) ----------------
    def parse_outcar(self):
        if not os.path.exists(self.outcar_path):
            print(" 错误: 找不到 OUTCAR")
            return

        # 1. 严格匹配正则: 匹配 "free energy" + 任意空格 + "TOTEN" + 任意空格 + "=" + 任意空格 + 数值 + 任意空格 + "eV"
        # 即使 VASP 输出格式中的空格数量变化，\s+ 也能确保匹配，且不会匹配到 "energy without entropy"
        re_energy = re.compile(r"free\s+energy\s+TOTEN\s+=\s+(-?\d+\.\d+)\s+eV")
        
        re_drift  = re.compile(r"total drift:\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)")
        re_pos_head = re.compile(r"^\s*POSITION\s+TOTAL-FORCE")
        re_press = re.compile(r"external pressure\s+=\s+(-?\d+\.\d+)\s+kB")

        with open(self.outcar_path,'r',errors='replace') as f:
            lines=f.readlines()

        current_step={}
        pos_buffer=[]
        in_pos_block=False

        for line in lines:
            if re_pos_head.search(line):
                in_pos_block=True
                pos_buffer=[]
                continue
            if in_pos_block:
                if "---" in line: continue
                parts=line.split()
                if len(parts)==6:
                    try:
                        coords=[float(x) for x in parts[0:3]]
                        forces=[float(x) for x in parts[3:6]]
                        pos_buffer.append({'c':coords,'f':forces})
                    except:
                        in_pos_block=False
                else:
                    in_pos_block=False
                    if pos_buffer:
                        f_norms=[np.linalg.norm(a['f']) for a in pos_buffer]
                        max_f=max(f_norms)
                        max_idx=np.argmax(f_norms)
                        
                        current_step['atoms']=pos_buffer
                        current_step['max_force']=max_f
                        current_step['max_atom_idx']=max_idx
                        pos_buffer=[]

            # --- 解析能量 (Energy) ---
            m_en=re_energy.search(line)
            if m_en:
                
                if 'max_force' in current_step:
                     self.steps.append(current_step)
                     current_step={} 
                
                current_step['energy']=float(m_en.group(1))

            # --- 解析 Drift ---
            m_dr=re_drift.search(line)
            if m_dr:
                drift_vec=np.array([float(x) for x in m_dr.groups()])
                current_step['drift']=np.linalg.norm(drift_vec)
                current_step['drift_vec']=drift_vec

            # --- 解析 Pressure ---
            m_pr=re_press.search(line)
            if m_pr:
                current_step['pressure']=float(m_pr.group(1))

        # 循环结束后，确保最后一步也被保存
        if current_step and 'max_force' in current_step and 'energy' in current_step:
            self.steps.append(current_step)

        print(f" OUTCAR 解析: {len(self.steps)} 个有效离子步数据")

    def diagnose_and_suggest(self):
        if not self.steps: return
        n_steps = min(len(self.steps), len(self.nelm_list)) if self.nelm_list else len(self.steps)

        print("\n"+"="*40)
        print(f"🔍 综合收敛诊断 (基于 {n_steps} 步数据)")
        print("="*40)

        all_max_atoms=[s['max_atom_idx'] for s in self.steps[:n_steps]]
        c=Counter(all_max_atoms)
        for atom_idx,count in c.most_common(3):
            max_f_history=max([s['max_force'] for s in self.steps if s['max_atom_idx']==atom_idx])
            name=self.full_atom_list[atom_idx] if self.full_atom_list and atom_idx<len(self.full_atom_list) else "?"
            print(f"  原子 #{atom_idx+1} ({name})")
            print(f"     -> 成为瓶颈次数: {count}/{n_steps} 步")
            print(f"     -> 历史最大受力: {max_f_history:.4f} eV/Å")
            print(f"     -> 建议: 检查该原子周围是否过度紧密或键长异常")
            print("-"*30)

        last_step=self.steps[n_steps-1]
        print(f"\n>>> 2. 最后一步状态 (Step {n_steps})")
        print(f"  当前能量: {last_step['energy']:.6f} eV")
        print(f"  当前最大力: {last_step['max_force']:.4f} eV/Å (原子 #{last_step['max_atom_idx']+1})")
        if 'drift_vec' in last_step:
            drift_vec=last_step['drift_vec']
            if last_step['drift']>self.Drift_thresh:
                max_dir=np.argmax(np.abs(drift_vec))
                print(f"  Drift 过大: {last_step['drift']:.4f} eV/Å, 最大方向: {['X','Y','Z'][max_dir]} 分量: {drift_vec[max_dir]:.4f}")

        if n_steps<=len(self.nelm_list):
            print(f"  当前电子步数: {self.nelm_list[n_steps-1]}")

        f_max=last_step['max_force']
        drift=last_step.get('drift',0)
        if f_max>1.0:
            print("力极大 (>1.0)。初始结构极差。")
        elif f_max>self.F_max_thresh:
            print(f"力较大 (> {self.F_max_thresh})。结构仍在调整。")
        elif drift>self.Drift_thresh:
            print(f"Drift ({drift:.3f}) 偏高 (> {self.Drift_thresh})")
        else:
            print("结构接近收敛或已收敛。")

    # ---------------- 5. 绘图 ----------------
    def plot_all(self):
        if not self.steps: return

        n_out=len(self.steps)
        n_osz=len(self.nelm_list)
        n_min=min(n_out,n_osz) if n_osz>0 else n_out
        nelms=self.nelm_list[:n_min] if n_osz>0 else [0]*n_min
        x_axis=list(range(1,n_min+1))
        energies=[s['energy'] for s in self.steps[:n_min]]
        forces=[s['max_force'] for s in self.steps[:n_min]]
        drifts=[s.get('drift',0) for s in self.steps[:n_min]]
        pressures=[s.get('pressure',0) for s in self.steps[:n_min]]

        fig_sum, axes=plt.subplots(3,2,figsize=(14,12),sharex=True)
        fig_sum.suptitle(f'VASP Optimization Summary (Total Steps: {n_min})',fontsize=16)

        # 能量
        ax=axes[0,0]; ax.plot(x_axis,energies,'b-o',ms=3)
        ax.set_ylabel('Energy (eV)'); ax.set_title('Total Energy'); ax.grid(True)
        if self.Energy_thresh is not None: ax.axhline(self.Energy_thresh,color='red',linestyle='--',label=f'Threshold={self.Energy_thresh}'); ax.legend()

        # 力
        ax=axes[0,1]; ax.plot(x_axis,forces,'r-x',ms=3); ax.set_yscale('log')
        ax.set_ylabel('Force (eV/Å)'); ax.set_title('Max Force (Log)'); ax.grid(True)
        ax.axhline(self.F_max_thresh,color='green',linestyle='--',label=f'F_thresh={self.F_max_thresh}'); ax.legend()

        # 漂移
        ax=axes[1,0]; ax.plot(x_axis,drifts,'g--',ms=2); ax.set_ylabel('Drift (eV/Å)'); ax.set_title('Total Drift'); ax.grid(True)
        ax.axhline(self.Drift_thresh,color='orange',linestyle='--',label=f'Drift_thresh={self.Drift_thresh}'); ax.legend()

        # NELM
        ax=axes[1,1]; ax.bar(x_axis,nelms,color='purple',alpha=0.6); ax.set_ylabel('SCF Steps'); ax.set_title('Electronic Steps'); ax.grid(True,axis='y')
        if self.NELM_thresh is not None: ax.axhline(self.NELM_thresh,color='red',linestyle='--',label=f'Threshold={self.NELM_thresh}'); ax.legend()

        # Pressure
        ax=axes[2,0]; ax.plot(x_axis,pressures,'k-s',ms=3); ax.set_ylabel('Pressure (kB)'); ax.set_title('External Pressure'); ax.set_xlabel('Ionic Step'); ax.grid(True)
        if self.Pressure_thresh is not None: ax.axhline(self.Pressure_thresh,color='purple',linestyle='--',label=f'Threshold={self.Pressure_thresh}'); ax.legend()

        # 信息框
        ax=axes[2,1]; ax.axis('off')
        last_f=forces[-1]; last_d=drifts[-1]; last_e=energies[-1]; last_nelm=nelms[-1]
        y_text=0.5 if n_min<20 else 0.3
        info_text=(
            f"Step {n_min}:\n\n"
            f"E = {last_e:.4f} eV\n"
            f"F_max = {last_f:.4f} eV/Å\n"
            f"Drift = {last_d:.4f}\n"
            f"NELM = {last_nelm}\n"
        )
        ax.text(0.1,y_text,info_text,fontsize=12,family='monospace')
        plt.tight_layout(rect=[0,0.03,1,0.95]); plt.savefig('plot_summary.png',dpi=150); plt.close(fig_sum)

        # 分图
        def save_fig(x,y,ylabel,title,filename,kind='line',ylog=False,thresh=None):
            plt.figure(figsize=(8,6))
            if kind=='line': plt.plot(x,y,'o-',ms=4,color='tab:blue')
            elif kind=='bar': plt.bar(x,y,alpha=0.7,color='tab:purple')
            plt.title(title); plt.xlabel('Ionic Step'); plt.ylabel(ylabel)
            if ylog: plt.yscale('log')
            if thresh is not None: plt.axhline(thresh,color='red',linestyle='--',label=f'Threshold={thresh}'); plt.legend()
            plt.grid(True); plt.tight_layout(); plt.savefig(filename,dpi=150); plt.close()

        save_fig(x_axis,energies,'Energy (eV)','Total Energy','plot_energy.png',thresh=self.Energy_thresh)
        save_fig(x_axis,forces,'Max Force (eV/Å)','Max Force (Log Scale)','plot_force.png',ylog=True,thresh=self.F_max_thresh)
        save_fig(x_axis,drifts,'Drift (eV/Å)','Total Drift','plot_drift.png',thresh=self.Drift_thresh)
        save_fig(x_axis,nelms,'SCF Steps','Electronic Steps per Ionic Step','plot_nelm.png',kind='bar',thresh=self.NELM_thresh)
        save_fig(x_axis,pressures,'Pressure (kB)','External Pressure','plot_pressure.png',thresh=self.Pressure_thresh)

        print(f"\n 绘图完成！已保存 6 张图片")

    # ---------------- 6. 保存最新结构 ----------------
    def save_poscar(self):
        if not self.steps: return
        last_step=self.steps[-1]
        with open("POSCAR_latest.vasp","w") as f:
            f.write("Generated by VaspDoctor\n1.0\n")
            if self.poscar_lattice:
                for v in self.poscar_lattice: f.write(f" {v[0]:12.8f} {v[1]:12.8f} {v[2]:12.8f}\n")
            else:
                f.write(" 10.0 0.0 0.0\n 0.0 10.0 0.0\n 0.0 0.0 10.0\n")
            if self.atom_names: f.write(" "+" ".join(self.atom_names)+"\n")
            if self.atom_counts: f.write(" "+" ".join(map(str,self.atom_counts))+"\n")
            f.write("Cartesian\n")
            for a in last_step['atoms']:
                c=a['c']
                f.write(f" {c[0]:12.8f} {c[1]:12.8f} {c[2]:12.8f}\n")
        print(" 结构已提取: POSCAR_latest.vasp")

# ------------------- 主程序 -------------------
if __name__=="__main__":
    f_out=sys.argv[1] if len(sys.argv)>1 else "OUTCAR"
    f_osz=sys.argv[2] if len(sys.argv)>2 else "OSZICAR"
    f_pos=sys.argv[3] if len(sys.argv)>3 else "POSCAR"

    doc=VaspDoctor(f_out,f_osz,f_pos)
    doc.read_poscar()
    doc.parse_oszicar()
    doc.parse_outcar()
    doc.diagnose_and_suggest()
    doc.save_poscar()
    doc.plot_all()