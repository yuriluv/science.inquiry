#!/usr/bin/env python3
"""
Clausius-Clapeyron 검증 실험 — 파이프라인 v2
독립 scanning: (δT, n_air_heat, n_air_cool) 3차원 최적화
물리적 제약: P_water > 0 (모든 포인트)
종료 조건: ΔH_vap ≈ 40.7 kJ/mol + 삼중점 P ≈ 0.00603 atm
"""

import math
import json
import sys

# ===== 상수 =====
R_atm = 0.08206   # L·atm/(mol·K)
R_J = 8.314        # J/(mol·K)
P_ATM = 1.000      # atm
rho = 998.0        # kg/m³ (물 20°C)
g = 9.80665        # m/s²
MENISCUS = 0.2     # mL

# 메스실린더 (10mL, 외경반지름 17mm, 높이 152mm)
d_inner_cm = 1.01  # cm (추정)
A_cm2 = math.pi * (d_inner_cm / 2) ** 2

# V_exp 고정값
V_EXP_HEAT = 2.0  # mL
V_EXP_COOL = 2.2  # mL

# 수압 보정 (고정)
h_heat_m = V_EXP_HEAT / A_cm2 * 1e-2
h_cool_m = V_EXP_COOL / A_cm2 * 1e-2
P_gas_heat = P_ATM - rho * g * h_heat_m / 101325
P_gas_cool = P_ATM - rho * g * h_cool_m / 101325

# 실험 데이터: (V_raw mL, T_water °C)
heat_data = [
    (2.8, 21.0), (3.0, 32.0), (3.1, 37.5), (3.2, 43.0),
    (3.3, 47.5), (3.4, 52.0), (3.5, 54.0), (3.6, 56.0),
    (3.7, 58.0), (3.8, 60.0), (3.9, 61.0), (4.0, 62.0),
]
# 냉각: 포인트 3,4 중복 제거 → 3개만
cool_data = [(3.0, 20.0), (2.9, 14.0), (2.8, 13.0)]

# 삼중점
T_TRIPLE_K = 273.16
P_TRIPLE_ATM = 0.00603  # 문헌값
DH_VAP_LIT = 40.7       # kJ/mol 문헌값


def compute_cc(dT, n_air_h, n_air_c, use_cool=True):
    """
    (δT, n_air_heat, n_air_cool) → C-C 회귀 결과
    P_water < 0인 포인트는 제외
    반환: dict or None (포인트 부족 시)
    """
    points_x = []
    points_y = []
    detail_heat = []
    detail_cool = []
    neg_count = 0

    for V_raw, T_w in heat_data:
        T_gas_K = T_w + 273.15 - dT
        V_corr = (V_raw - MENISCUS) / 1000.0
        P_air = n_air_h * R_atm * T_gas_K / V_corr
        P_water = P_gas_heat - P_air
        P_lit = math.exp(13.722 - 5120.6 / T_gas_K) if T_gas_K > 0 else 0
        detail_heat.append({
            'T_w': T_w, 'T_gas_K': T_gas_K, 'V_raw': V_raw,
            'V_corr': V_corr * 1000, 'P_gas': P_gas_heat,
            'P_air': P_air, 'P_water': P_water, 'P_lit': P_lit
        })
        if P_water <= 0:
            neg_count += 1
            continue
        points_x.append(1.0 / T_gas_K)
        points_y.append(math.log(P_water))

    if use_cool:
        for V_raw, T_w in cool_data:
            T_gas_K = T_w + 273.15 - dT
            V_corr = (V_raw - MENISCUS) / 1000.0
            P_air = n_air_c * R_atm * T_gas_K / V_corr
            P_water = P_gas_cool - P_air
            P_lit = math.exp(13.722 - 5120.6 / T_gas_K) if T_gas_K > 0 else 0
            detail_cool.append({
                'T_w': T_w, 'T_gas_K': T_gas_K, 'V_raw': V_raw,
                'V_corr': V_corr * 1000, 'P_gas': P_gas_cool,
                'P_air': P_air, 'P_water': P_water, 'P_lit': P_lit
            })
            if P_water <= 0:
                neg_count += 1
                continue
            points_x.append(1.0 / T_gas_K)
            points_y.append(math.log(P_water))

    n = len(points_x)
    if n < 4:
        return None

    # 선형 회귀
    sx = sum(points_x)
    sy = sum(points_y)
    sxy = sum(x * y for x, y in zip(points_x, points_y))
    sx2 = sum(x ** 2 for x in points_x)

    slope = (n * sxy - sx * sy) / (n * sx2 - sx ** 2)
    intercept = (sy - slope * sx) / n

    y_mean = sy / n
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(points_x, points_y))
    ss_tot = sum((y - y_mean) ** 2 for y in points_y)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    dH = -slope * R_J / 1000  # kJ/mol
    P_triple = math.exp(slope / T_TRIPLE_K + intercept)

    # 물리적 제약: P_water > 0
    if neg_count > 0:
        valid = False
    else:
        valid = True

    return {
        'dT': dT, 'n_air_h': n_air_h, 'n_air_c': n_air_c,
        'slope': slope, 'intercept': intercept, 'R2': r2,
        'dH': dH, 'P_triple': P_triple,
        'n_points': n, 'neg_count': neg_count, 'valid': valid,
        'cost_dH': abs(dH - DH_VAP_LIT) / DH_VAP_LIT,
        'cost_triple': abs(P_triple - P_TRIPLE_ATM) / P_TRIPLE_ATM,
        'detail_heat': detail_heat, 'detail_cool': detail_cool,
    }


def scan(dT_range, n_air_range, use_cool=True, label=""):
    """2차원 scanning: δT × n_air"""
    results = []
    best = None
    best_cost = float('inf')

    dT_values = [dT_range[0] + i * dT_range[2] for i in range(int((dT_range[1] - dT_range[0]) / dT_range[2]) + 1)]
    n_values = [n_air_range[0] + i * n_air_range[2] for i in range(int((n_air_range[1] - n_air_range[0]) / n_air_range[2]) + 1)]

    for dT in dT_values:
        for n_h in n_values:
            n_c = n_h  # 같은 스케일에서 시작
            # 냉각 n_air도 독립 탐색
            for n_c_i in [n_h] if not use_cool else [n_h * 0.85, n_h, n_h * 1.15]:
                r = compute_cc(dT, n_h, n_c_i, use_cool=use_cool)
                if r is None:
                    continue
                cost = r['cost_dH'] + r['cost_triple']
                # P_water < 0 이면 패널티
                if not r['valid']:
                    cost += 10.0  # 큰 패널티
                r['cost'] = cost
                results.append(r)
                if cost < best_cost:
                    best_cost = cost
                    best = r

    if best:
        print(f"\n{'='*70}")
        print(f"SCAN RESULT [{label}]")
        print(f"{'='*70}")
        print(f"최적 δT = {best['dT']:.2f} °C")
        print(f"n_air_가열 = {best['n_air_h']:.8f} mol")
        print(f"n_air_냉각 = {best['n_air_c']:.8f} mol")
        print(f"ΔH_vap = {best['dH']:.2f} kJ/mol (문헌 {DH_VAP_LIT}, 오차 {best['cost_dH']*100:.1f}%)")
        print(f"삼중점 = {best['P_triple']:.5f} atm (문헌 {P_TRIPLE_ATM}, 오차 {best['cost_triple']*100:.1f}%)")
        print(f"R² = {best['R2']:.6f}")
        print(f"포인트 수 = {best['n_points']}, P_water<0 = {best['neg_count']}")
        print(f"cost = {best['cost']:.4f}")

        # 상세 데이터
        print(f"\n가열 데이터 상세:")
        print(f"{'T_w':>6} {'T_gas':>8} {'V_raw':>6} {'P_gas':>10} {'P_air':>10} {'P_water':>10} {'P_lit':>10} {'비율':>8}")
        for d in best['detail_heat']:
            ratio = d['P_water'] / d['P_lit'] if d['P_lit'] > 0 else float('inf')
            mark = " ***NEG" if d['P_water'] <= 0 else ""
            print(f"{d['T_w']:6.1f} {d['T_gas_K']:8.2f} {d['V_raw']:6.1f} "
                  f"{d['P_gas']:10.6f} {d['P_air']:10.6f} {d['P_water']:10.6f} {d['P_lit']:10.6f} {ratio:8.3f}{mark}")

        if use_cool:
            print(f"\n냉각 데이터 상세:")
            for d in best['detail_cool']:
                ratio = d['P_water'] / d['P_lit'] if d['P_lit'] > 0 else float('inf')
                mark = " ***NEG" if d['P_water'] <= 0 else ""
                print(f"{d['T_w']:6.1f} {d['T_gas_K']:8.2f} {d['V_raw']:6.1f} "
                      f"{d['P_gas']:10.6f} {d['P_air']:10.6f} {d['P_water']:10.6f} {d['P_lit']:10.6f} {ratio:8.3f}{mark}")

    return best, results


if __name__ == "__main__":
    print("C-C 검증 파이프라인 v2")
    print(f"메스실린더 내경 = {d_inner_cm} cm, A = {A_cm2:.4f} cm²")
    print(f"P_gas_가열 = {P_gas_heat:.6f} atm")
    print(f"P_gas_냉각 = {P_gas_cool:.6f} atm")
    print()

    # Phase 1: 넓은 범위 coarse scan (가열만)
    print("Phase 1: Coarse scan (가열 데이터만)")
    print("δT: 0.0~15.0 step 0.5, n_air: 8e-5~1.5e-4 step 1e-5")
    best1, _ = scan(
        dT_range=(0.0, 15.0, 0.5),
        n_air_range=(8e-5, 1.5e-4, 1e-5),
        use_cool=False,
        label="가열 only (coarse)"
    )

    # Phase 2: Fine scan around best (가열만)
    if best1:
        dT_center = best1['dT']
        n_center = best1['n_air_h']
        print(f"\n\nPhase 2: Fine scan (가열 데이터만)")
        print(f"δT: {dT_center-1:.1f}~{dT_center+1:.1f} step 0.1, n_air: {n_center-2e-5:.2e}~{n_center+2e-5:.2e} step 2e-6")
        best2, _ = scan(
            dT_range=(max(0, dT_center - 1), dT_center + 1, 0.1),
            n_air_range=(max(1e-6, n_center - 2e-5), n_center + 2e-5, 2e-6),
            use_cool=False,
            label="가열 only (fine)"
        )

        # Phase 3: 가열+냉각 (가열 최적값 근처에서 냉각 n_air 탐색)
        if best2:
            print(f"\n\nPhase 3: 가열+냉각 scanning")
            best3, results3 = scan(
                dT_range=(max(0, best2['dT'] - 2), best2['dT'] + 2, 0.2),
                n_air_range=(max(1e-6, best2['n_air_h'] - 1e-5), best2['n_air_h'] + 1e-5, 1e-6),
                use_cool=True,
                label="가열+냉각"
            )
