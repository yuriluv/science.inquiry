#!/usr/bin/env python3
"""보고서 표를 이미지로 변환"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# 폰트 설정
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 11

def save_table(fig, path, title_text=None):
    if title_text:
        fig.text(0.5, 0.97, title_text, ha='center', va='top', fontsize=13, fontweight='bold')
    fig.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"저장: {path}")

# ============================================================
# 표 1: 가열 데이터 (12포인트)
# ============================================================
heat_data = [
    ['2.8', '21.0', '292.55', '2.6', '0.99759'],
    ['3.0', '32.0', '303.55', '2.8', '0.99759'],
    ['3.1', '37.5', '309.05', '2.9', '0.99759'],
    ['3.2', '43.0', '314.55', '3.0', '0.99759'],
    ['3.3', '47.5', '318.55', '3.1', '0.99759'],
    ['3.4', '52.0', '323.55', '3.2', '0.99759'],
    ['3.5', '54.0', '325.55', '3.3', '0.99759'],
    ['3.6', '56.0', '327.55', '3.4', '0.99759'],
    ['3.7', '58.0', '329.55', '3.5', '0.99759'],
    ['3.8', '60.0', '331.55', '3.6', '0.99759'],
    ['3.9', '61.0', '332.55', '3.7', '0.99759'],
    ['4.0', '62.0', '333.55', '3.8', '0.99759'],
]
col_labels = ['V_raw\n(mL)', 'T_water\n(°C)', 'T_gas\n(K)', 'V_corr\n(mL)', 'P_gas\n(atm)']

fig, ax = plt.subplots(figsize=(8, 5))
ax.axis('off')
ax.set_title('[Table 1] Heating Data (12 points)', fontsize=14, fontweight='bold', pad=20)
ax.text(0.5, 0.95, 'T_gas = T_water - δT (δT: correction param)\nV_corr = V_raw - 0.2 mL (meniscus)\nP_gas = P_ATM - ρgh,  V_exp = 2.0 mL (heating)',
        transform=ax.transAxes, ha='center', va='top', fontsize=9, style='italic', color='gray')

table = ax.table(cellText=heat_data, colLabels=col_labels, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.0, 1.4)
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')
    elif key[0] % 2 == 0:
        cell.set_facecolor('#ecf0f1')
    else:
        cell.set_facecolor('#ffffff')
    cell.set_edgecolor('#bdc3c7')

save_table(fig, '/data/data/com.termux/files/home/science.inquiry/results/table1_heat_data.png')

# ============================================================
# 표 2: 냉각 데이터 (3포인트)
# ============================================================
cool_data = [
    ['3.0', '20.0', '2.8', '0.99744'],
    ['2.9', '14.0', '2.7', '0.99744'],
    ['2.8', '13.0', '2.6', '0.99744'],
]
col_labels2 = ['V_raw\n(mL)', 'T_water\n(°C)', 'V_corr\n(mL)', 'P_gas\n(atm)']

fig, ax = plt.subplots(figsize=(7, 3.5))
ax.axis('off')
ax.set_title('[Table 2] Cooling Data (3 points, deduplicated)', fontsize=14, fontweight='bold', pad=20)
ax.text(0.5, 0.88, 'V_exp = 2.2 mL (cooling)\nPoints 3,4 (13°C, 2.8mL) are duplicates — 1 used',
        transform=ax.transAxes, ha='center', va='top', fontsize=9, style='italic', color='gray')

table = ax.table(cellText=cool_data, colLabels=col_labels2, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(11)
table.scale(1.1, 1.5)
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')
    elif key[0] % 2 == 0:
        cell.set_facecolor('#ecf0f1')
    else:
        cell.set_facecolor('#ffffff')
    cell.set_edgecolor('#bdc3c7')

save_table(fig, '/data/data/com.termux/files/home/science.inquiry/results/table2_cool_data.png')

# ============================================================
# 표 3: Case별 종합 비교
# ============================================================
case_data = [
    ['A\n(heating only)', '1.6', '1.06×10⁻⁴', '—', '44.98', '0.00612', '0.978', '10.5%', '1.5%'],
    ['B\n(+cool 20°C)', '0.0', '1.07×10⁻⁴', '1.01×10⁻⁴', '42.58', '0.00604', '0.534', '4.6%', '0.2%'],
    ['C\n(full, cost)', '1.5', '1.08×10⁻⁴', '1.04×10⁻⁴', '41.76', '0.00607', '0.213', '2.6%', '0.7%'],
    ['C\n(full, R²)', '0.8', '1.07×10⁻⁴', '1.10×10⁻⁴', '43.03', '0.00619', '0.772', '5.7%', '2.6%'],
    ['D\n(anchor, w=1)', '1.4', '1.06×10⁻⁴', 'anchor', '44.95', '0.00604', '0.990', '10.4%', '0.1%'],
    ['D\n(anchor, w=3)', '1.1', '1.06×10⁻⁴', 'anchor', '44.59', '0.00604', '0.994', '9.6%', '0.1%'],
    ['D\n(anchor, w=5)', '0.7', '1.06×10⁻⁴', 'anchor', '44.17', '0.00602', '0.994', '8.5%', '0.2%'],
]
col_labels3 = ['Case', 'δT\n(°C)', 'n_air_H\n(mol)', 'n_air_C\n(mol)', 'ΔH_vap\n(kJ/mol)', 'P_triple\n(atm)', 'R²', 'ΔH\nerror', 'triple\nerror']

fig, ax = plt.subplots(figsize=(12, 5.5))
ax.axis('off')
ax.set_title('[Table 3] Comparative Analysis: All Cases', fontsize=14, fontweight='bold', pad=20)
ax.text(0.5, 0.95, 'Literature values: ΔH_vap = 40.7 kJ/mol,  P_triple = 0.00603 atm',
        transform=ax.transAxes, ha='center', va='top', fontsize=10, style='italic', color='#2c3e50')

table = ax.table(cellText=case_data, colLabels=col_labels3, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(9.5)
table.scale(1.0, 1.8)

# Highlight best Case D (w=5) row
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')
    elif key[0] == 7:  # Case D w=5 row
        cell.set_facecolor('#d5f5e3')
        cell.set_text_props(fontweight='bold')
    elif key[0] % 2 == 0:
        cell.set_facecolor('#ecf0f1')
    else:
        cell.set_facecolor('#ffffff')
    cell.set_edgecolor('#bdc3c7')

save_table(fig, '/data/data/com.termux/files/home/science.inquiry/results/table3_case_comparison.png')

# ============================================================
# 표 4: 보정 파라미터 요약
# ============================================================
param_data = [
    ['δT', 'Thermal offset\n(T_gas = T_water - δT)', 'Independent scan\n(optimal: 0.7~1.6°C)'],
    ['n_air_H', 'Air moles (heating)', 'Independent scan\n(optimal: ~1.06×10⁻⁴)'],
    ['n_air_C', 'Air moles (cooling)', 'Independent scan\n(separate from heating)'],
    ['V_exp (H)', 'Hydrostatic height\n(heating)', 'Fixed: 2.0 mL'],
    ['V_exp (C)', 'Hydrostatic height\n(cooling)', 'Fixed: 2.2 mL'],
    ['Meniscus', 'Meniscus correction', 'Fixed: 0.2 mL'],
    ['d_inner', 'Cylinder inner diameter', 'Inverse calc: ≈1.01 cm'],
    ['P_ATM', 'Atmospheric pressure', '1.000 atm (assumed)'],
]
col_labels4 = ['Parameter', 'Role', 'Treatment']

fig, ax = plt.subplots(figsize=(11, 6))
ax.axis('off')
ax.set_title('[Table 4] Correction Parameters', fontsize=14, fontweight='bold', pad=15)

table = ax.table(cellText=param_data, colLabels=col_labels4, loc='center', cellLoc='left')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.0, 1.8)
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')
    elif key[0] % 2 == 0:
        cell.set_facecolor('#ecf0f1')
    else:
        cell.set_facecolor('#ffffff')
    cell.set_edgecolor('#bdc3c7')
    # First column bold
    if key[0] > 0 and key[1] == 0:
        cell.set_text_props(fontweight='bold')

# Set column widths
col_widths = [0.15, 0.35, 0.50]
for key, cell in table.get_celld().items():
    cell.set_width(col_widths[key[1]])

save_table(fig, '/data/data/com.termux/files/home/science.inquiry/results/table4_parameters.png')

# ============================================================
# 표 5: 내경 추정 교차 검증
# ============================================================
diam_data = [
    ['120', '10.30', '0.8333', 'Smallest base assumption'],
    ['122', '10.22', '0.8197', ''],
    ['124', '10.13', '0.8065', ''],
    ['125', '10.09', '0.8000', '← Adopted (d ≈ 1.01 cm)'],
    ['126', '10.05', '0.7937', ''],
    ['128', '9.97', '0.7812', ''],
    ['130', '9.90', '0.7692', 'Largest base assumption'],
]
col_labels5 = ['Marked height\nh_marked (mm)', 'd_inner\n(mm)', 'A\n(cm²)', 'Note']

fig, ax = plt.subplots(figsize=(10, 4.5))
ax.axis('off')
ax.set_title('[Table 5] Cylinder Inner Diameter Cross-Validation', fontsize=14, fontweight='bold', pad=15)
ax.text(0.5, 0.93, 'd_inner = 2 × √(V / (π × h_marked)),  V = 10.0 mL,  d_adopted ≈ 1.01 cm',
        transform=ax.transAxes, ha='center', va='top', fontsize=9, style='italic', color='gray')

table = ax.table(cellText=diam_data, colLabels=col_labels5, loc='center', cellLoc='center')
table.auto_set_font_size(False)
table.set_fontsize(10)
table.scale(1.1, 1.6)
for key, cell in table.get_celld().items():
    if key[0] == 0:
        cell.set_facecolor('#2c3e50')
        cell.set_text_props(color='white', fontweight='bold')
    elif key[0] == 4:  # Adopted row
        cell.set_facecolor('#d5f5e3')
        cell.set_text_props(fontweight='bold')
    elif key[0] % 2 == 0:
        cell.set_facecolor('#ecf0f1')
    else:
        cell.set_facecolor('#ffffff')
    cell.set_edgecolor('#bdc3c7')

save_table(fig, '/data/data/com.termux/files/home/science.inquiry/results/table5_diameter.png')

print("\n모든 표 이미지 생성 완료!")