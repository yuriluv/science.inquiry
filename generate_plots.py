#!/usr/bin/env python3
"""C-C 검증 실험 — 결과 그래프 생성 스크립트"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import json
import math
import os

# 한글 폰트 설정
plt.rcParams['font.family'] = ['Noto Sans CJK KR', 'DejaVu Sans', 'sans-serif']
plt.rcParams['axes.unicode_minus'] = False

# 데이터 로드
with open(os.path.join(os.path.dirname(__file__), 'plot_data.json'), 'r') as f:
    data = json.load(f)

# 물리 상수
R_J = 8.314
T_TRIPLE_K = 273.16
P_TRIPLE_ATM = 0.00603
DH_VAP_LIT = 40.7
DH_VAP_LIT_J = 40700

def P_water_lit(T_c):
    T = T_c + 273.15
    return math.exp(13.722 - 5120.6 / T)

# 문헌값 C-C 선
T_range = np.linspace(0, 100, 200)
T_K = T_range + 273.15
P_lit = np.array([P_water_lit(T) for T in T_range])
lit_x = 1.0 / T_K
lit_y = np.log(P_lit)
lit_slope = -DH_VAP_LIT_J / R_J
lit_intercept = math.log(P_TRIPLE_ATM) - lit_slope / T_TRIPLE_K

# ===== Figure 1: C-C Plot (전체 Case 비교) =====
fig, axes = plt.subplots(1, 2, figsize=(16, 7))

# Left: C-C plot (ln P vs 1/T)
ax = axes[0]
ax.plot(lit_x, lit_y, 'k-', linewidth=2, label=f'Literature (ΔH={DH_VAP_LIT} kJ/mol)', zorder=5)

# Case A 회귀선
a = data['case_a']
x_fit = np.linspace(1/373.15, 1/273.16, 200)
y_fit_a = a['slope'] * x_fit + a['intercept']
ax.plot(x_fit, y_fit_a, 'r--', linewidth=1.5, label=f'A: ΔH={a["dH"]:.1f} kJ/mol, P\u2083={a["P_triple"]:.4f}')

# Case C 회귀선
c = data['case_c']
y_fit_c = c['slope'] * x_fit + c['intercept']
ax.plot(x_fit, y_fit_c, 'b-.', linewidth=1.5, label=f'C: ΔH={c["dH"]:.1f} kJ/mol, P\u2083={c["P_triple"]:.4f}')

# Case D 회귀선
d = data['case_d']
y_fit_d = d['slope'] * x_fit + d['intercept']
ax.plot(x_fit, y_fit_d, 'g-', linewidth=2, label=f'D: ΔH={d["dH"]:.1f} kJ/mol, P\u2083={d["P_triple"]:.4f}')

# 데이터 포인트
ax.scatter(data['heat_x'], data['heat_y'], c='red', s=50, zorder=10, label='Heating data (Case A)')
ax.scatter(data['cool_x'], data['cool_y'], c='blue', s=50, marker='s', zorder=10, label='Cooling data')

# 앵커 포인트
ax.scatter([data['anchor_x']], [data['anchor_y']], c='green', s=100, marker='*', zorder=11, 
           label=f'Cooling extrapolation anchor')

# 삼중점 표시
ax.scatter([1/273.16], [math.log(P_TRIPLE_ATM)], c='black', s=80, marker='D', zorder=11,
           label=f'Triple point (0.00603 atm)')

ax.set_xlabel('1/T (K$^{-1}$)', fontsize=13)
ax.set_ylabel('ln(P$_{water}$ / atm)', fontsize=13)
ax.set_title('Clausius-Clapeyron: ln(P) vs 1/T', fontsize=14, fontweight='bold')
ax.legend(fontsize=9, loc='upper right')
ax.grid(True, alpha=0.3)
ax.invert_xaxis()

# Right: Trade-off scatter
ax2 = axes[1]
cases_data = [
    ('A: Heat only', 10.5, 1.5, 0.978, 'red', 'o'),
    ('B: +Cool 20°C', 4.6, 0.2, 0.534, 'blue', 's'),
    ('C: All (cost)', 2.6, 0.7, 0.213, 'purple', '^'),
    ('C: All (R²)', 5.7, 2.6, 0.772, 'cyan', 'D'),
    ('D: w=1', 10.4, 0.1, 0.990, 'lime', 'p'),
    ('D: w=3', 9.6, 0.1, 0.994, 'green', 'P'),
    ('D: w=5', 8.5, 0.2, 0.994, 'darkgreen', '*'),
]

for name, dh_err, tri_err, r2, color, marker in cases_data:
    size = 150 if 'D' in name else 100
    ax2.scatter(dh_err, tri_err, c=color, s=size, marker=marker, 
                label=f'{name} (R²={r2:.3f})', zorder=5, edgecolors='black', linewidth=0.5)

ax2.axvline(x=5, color='gray', linestyle=':', alpha=0.5, label='5% threshold')
ax2.axhline(y=5, color='gray', linestyle=':', alpha=0.5)
ax2.axvspan(0, 5, alpha=0.1, color='green')
ax2.axhspan(0, 5, alpha=0.1, color='green')

# 교집합 강조
ax2.fill_between([0, 5], [0, 0], [5, 5], alpha=0.15, color='green', label='Target zone (< 5% both)')

ax2.set_xlabel('ΔH$_{vap}$ Error (%)', fontsize=13)
ax2.set_ylabel('Triple Point Error (%)', fontsize=13)
ax2.set_title('Trade-off: ΔH vs Triple Point Accuracy', fontsize=14, fontweight='bold')
ax2.legend(fontsize=8, loc='upper right')
ax2.grid(True, alpha=0.3)
ax2.set_xlim(0, 12)
ax2.set_ylim(0, 4)

plt.tight_layout()
plt.savefig('/data/data/com.termux/files/home/science.inquiry/results/cc_comparison.png', dpi=200, bbox_inches='tight')
print("그래프 1 저장: cc_comparison.png")

# ===== Figure 2: Case D 상세 (앵커 효과) =====
fig2, ax3 = plt.subplots(figsize=(10, 7))

# 문헌값
ax3.plot(1.0/T_K, np.log(P_lit), 'k-', linewidth=2.5, label=f'Literature (ΔH={DH_VAP_LIT} kJ/mol)', zorder=1)

# Case A 회귀선
y_a = a['slope'] * x_fit + a['intercept']
ax3.plot(x_fit, y_a, 'r--', linewidth=1.8, alpha=0.7, label=f'Case A: ΔH={a["dH"]:.1f}, R²=0.978')

# Case D 회귀선
y_d = d['slope'] * x_fit + d['intercept']
ax3.plot(x_fit, y_d, 'g-', linewidth=2.2, label=f'Case D (w=3): ΔH={d["dH"]:.1f}, R²=0.994')

# 데이터 포인트
ax3.scatter(data['heat_x'], data['heat_y'], c='red', s=60, zorder=10, label='Heating data', edgecolors='darkred')
ax3.scatter(data['cool_x'], data['cool_y'], c='royalblue', s=60, marker='s', zorder=10, label='Cooling data (not used in D)', edgecolors='navy')

# 앵커
ax3.scatter([data['anchor_x']], [data['anchor_y']], c='limegreen', s=200, marker='*', zorder=11, 
           edgecolors='darkgreen', linewidth=1.5, label='Extrapolation anchor (0°C)')

# 삼중점
ax3.scatter([1/273.16], [math.log(P_TRIPLE_ATM)], c='black', s=100, marker='D', zorder=11,
           label='Triple point (literature)')

# 삼중점 추정치 (Case A)
ax3.scatter([1/273.16], [a['slope']/273.16 + a['intercept']], c='red', s=80, marker='v', zorder=11,
           label=f'Case A triple pt est. ({a["P_triple"]:.4f} atm)')

# 삼중점 추정치 (Case D)  
ax3.scatter([1/273.16], [d['slope']/273.16 + d['intercept']], c='green', s=80, marker='^', zorder=11,
           label=f'Case D triple pt est. ({d["P_triple"]:.4f} atm)')

# 화살표 (앵커가 개선하는 방향)
ax3.annotate('Anchor pulls\nregression toward\ntriple point', 
             xy=(1/273.16 + 0.00005, math.log(P_TRIPLE_ATM) + 0.15),
             fontsize=10, color='darkgreen', fontweight='bold',
             arrowprops=dict(arrowstyle='->', color='darkgreen', lw=1.5),
             xytext=(1/310, math.log(0.03)))

ax3.set_xlabel('1/T (K$^{-1}$)', fontsize=14)
ax3.set_ylabel('ln(P$_{water}$ / atm)', fontsize=14)
ax3.set_title('Effect of Extrapolation Anchor on Triple Point Estimation', fontsize=14, fontweight='bold')
ax3.legend(fontsize=9, loc='upper right')
ax3.grid(True, alpha=0.3)
ax3.invert_xaxis()

plt.tight_layout()
plt.savefig('/data/data/com.termux/files/home/science.inquiry/results/cc_anchor_detail.png', dpi=200, bbox_inches='tight')
print("그래프 2 저장: cc_anchor_detail.png")

# ===== Figure 3: 민감도 분석 (δT vs n_air landscape) =====
fig3, axes3 = plt.subplots(1, 2, figsize=(16, 7))

# Cost landscape (Case A)
from cc_analysis import compute_cc, P_GAS_HEAT, HEAT_DATA as HD, MENISCUS as MENS

dT_vals = np.arange(0.5, 8.1, 0.2)
n_vals = np.arange(8e-5, 1.4e-4, 1e-6)
cost_grid = np.zeros((len(n_vals), len(dT_vals)))
dH_grid = np.zeros_like(cost_grid)
tri_grid = np.zeros_like(cost_grid)

for i, n in enumerate(n_vals):
    for j, dT in enumerate(dT_vals):
        r = compute_cc(dT, n, case="A")
        if r and r['valid']:
            cost_grid[i, j] = r['cost_dH'] + r['cost_triple']
            dH_grid[i, j] = r['cost_dH'] * 100
            tri_grid[i, j] = r['cost_triple'] * 100
        else:
            cost_grid[i, j] = np.nan
            dH_grid[i, j] = np.nan
            tri_grid[i, j] = np.nan

# Left: ΔH error heatmap
ax_l = axes3[0]
im = ax_l.pcolormesh(dT_vals, n_vals * 1e4, dH_grid, cmap='RdYlGn_r', vmin=0, vmax=25)
ax_l.contour(dT_vals, n_vals * 1e4, dH_grid, levels=[5, 10, 15, 20], colors='black', linewidths=0.8)
plt.colorbar(im, ax=ax_l, label='ΔH$_{vap}$ error (%)')
ax_l.set_xlabel('δT (°C)', fontsize=12)
ax_l.set_ylabel('n$_{air}$ (×10$^{-4}$ mol)', fontsize=12)
ax_l.set_title('ΔH$_{vap}$ Error Landscape (Case A)', fontsize=13, fontweight='bold')
ax_l.scatter([1.6], [1.06], c='white', s=100, marker='*', zorder=10, edgecolors='black', linewidth=1)

# Right: Triple point error heatmap
ax_r = axes3[1]
im2 = ax_r.pcolormesh(dT_vals, n_vals * 1e4, tri_grid, cmap='RdYlGn_r', vmin=0, vmax=25)
ax_r.contour(dT_vals, n_vals * 1e4, tri_grid, levels=[1, 2, 5, 10, 20], colors='black', linewidths=0.8)
plt.colorbar(im2, ax=ax_r, label='Triple point error (%)')
ax_r.set_xlabel('δT (°C)', fontsize=12)
ax_r.set_ylabel('n$_{air}$ (×10$^{-4}$ mol)', fontsize=12)
ax_r.set_title('Triple Point Error Landscape (Case A)', fontsize=13, fontweight='bold')
ax_r.scatter([1.6], [1.06], c='white', s=100, marker='*', zorder=10, edgecolors='black', linewidth=1)

plt.tight_layout()
plt.savefig('/data/data/com.termux/files/home/science.inquiry/results/cc_sensitivity.png', dpi=200, bbox_inches='tight')
print("그래프 3 저장: cc_sensitivity.png")

# ===== Figure 4: 최종 비교 바 차트 =====
fig4, (ax_b1, ax_b2) = plt.subplots(1, 2, figsize=(14, 6))

names = ['A: Heat', 'B: +20°C', 'C: cost', 'C: R²', 'D: w=1', 'D: w=3', 'D: w=5']
dh_errors = [10.5, 4.6, 2.6, 5.7, 10.4, 9.6, 8.5]
tri_errors = [1.5, 0.2, 0.7, 2.6, 0.1, 0.1, 0.2]
r2_vals = [0.978, 0.534, 0.213, 0.772, 0.990, 0.994, 0.994]
colors = ['#e74c3c', '#3498db', '#9b59b6', '#1abc9c', '#2ecc71', '#27ae60', '#1e8449']

x = np.arange(len(names))
width = 0.35

bars1 = ax_b1.bar(x - width/2, dh_errors, width, label='ΔH$_{vap}$ error (%)', color=colors, alpha=0.8, edgecolor='black')
bars2 = ax_b1.bar(x + width/2, tri_errors, width, label='Triple point error (%)', color=colors, alpha=0.5, 
                   edgecolor='black', hatch='//')
ax_b1.axhline(y=5, color='red', linestyle='--', alpha=0.7, linewidth=1.5, label='5% target')
ax_b1.set_ylabel('Error (%)', fontsize=12)
ax_b1.set_title('Error by Case', fontsize=13, fontweight='bold')
ax_b1.set_xticks(x)
ax_b1.set_xticklabels(names, fontsize=9, rotation=15)
ax_b1.legend(fontsize=9)
ax_b1.grid(True, alpha=0.3, axis='y')

# R² comparison
bars3 = ax_b2.bar(x, r2_vals, 0.6, color=colors, edgecolor='black')
ax_b2.axhline(y=0.95, color='orange', linestyle='--', alpha=0.7, linewidth=1.5, label='R²=0.95')
for i, v in enumerate(r2_vals):
    ax_b2.text(i, v + 0.02, f'{v:.3f}', ha='center', fontsize=9, fontweight='bold')
ax_b2.set_ylabel('R²', fontsize=12)
ax_b2.set_title('R² by Case', fontsize=13, fontweight='bold')
ax_b2.set_xticks(x)
ax_b2.set_xticklabels(names, fontsize=9, rotation=15)
ax_b2.set_ylim(0, 1.1)
ax_b2.legend(fontsize=9)
ax_b2.grid(True, alpha=0.3, axis='y')

plt.tight_layout()
plt.savefig('/data/data/com.termux/files/home/science.inquiry/results/cc_bar_comparison.png', dpi=200, bbox_inches='tight')
print("그래프 4 저장: cc_bar_comparison.png")

print("\n모든 그래프 생성 완료!")