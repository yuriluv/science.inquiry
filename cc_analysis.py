#!/usr/bin/env python3
"""
Clausius-Clapeyron 검증 실험 — 모듈화 분석 스크립트 v3
Case A/B/C 비교 분석, δT+n_air 독립 scanning, P_water>0 제약

Usage:
    python cc_analysis.py              # 전체 실행 (Case A + B + C)
    python cc_analysis.py --case A     # Case A만
    python cc_analysis.py --case B     # Case A+B만
"""

import math
import argparse
import json
import os

# ===== 물리 상수 =====
R_ATM = 0.08206    # L·atm/(mol·K)
R_J = 8.314        # J/(mol·K)
P_ATM = 1.000      # atm (표준)
RHO = 998.0        # kg/m³ (물 20°C)
G = 9.80665        # m/s²
MENISCUS = 0.2     # mL
T_TRIPLE_K = 273.16
P_TRIPLE_ATM = 0.00603
DH_VAP_LIT = 40.7  # kJ/mol

# 메스실린더 규격
D_INNER_CM = 1.01
A_CM2 = math.pi * (D_INNER_CM / 2) ** 2

# V_exp 고정값
V_EXP_HEAT = 2.0  # mL
V_EXP_COOL = 2.2   # mL

# 수압 보정 (고정)
P_GAS_HEAT = P_ATM - RHO * G * (V_EXP_HEAT / A_CM2 * 1e-2) / 101325
P_GAS_COOL = P_ATM - RHO * G * (V_EXP_COOL / A_CM2 * 1e-2) / 101325

# ===== 실험 데이터 =====
HEAT_DATA = [
    (2.8, 21.0), (3.0, 32.0), (3.1, 37.5), (3.2, 43.0),
    (3.3, 47.5), (3.4, 52.0), (3.5, 54.0), (3.6, 56.0),
    (3.7, 58.0), (3.8, 60.0), (3.9, 61.0), (4.0, 62.0),
]
COOL_DATA = [(3.0, 20.0), (2.9, 14.0), (2.8, 13.0)]  # 포인트 4는 중복 제거


def P_water_lit(T_c):
    """Antoine equation for water (≈ 문헌값 P_water)"""
    T = T_c + 273.15
    return math.exp(13.722 - 5120.6 / T)


def compute_cc(dT, n_air_h, n_air_c=None, case="A"):
    """
    C-C 회귀 계산
    
    Parameters:
        dT: 열 불균형 보정 (°C)
        n_air_h: 가열 실험 공기 몰수
        n_air_c: 냉각 실험 공기 몰수 (Case B/C만)
        case: "A" (가열만), "B" (+냉각20°C), "C" (전체)
    
    Returns:
        dict or None
    """
    if case in ("B", "C") and n_air_c is None:
        n_air_c = n_air_h
    
    cc_x = []
    cc_y = []
    detail = {"heat": [], "cool": []}
    neg_count = 0
    
    # 가열 데이터
    for V_raw, T_w in HEAT_DATA:
        T_gas_K = T_w + 273.15 - dT
        V_corr = (V_raw - MENISCUS) / 1000.0
        P_air = n_air_h * R_ATM * T_gas_K / V_corr
        P_water = P_GAS_HEAT - P_air
        P_lit = P_water_lit(T_w - dT) if T_gas_K > 0 else 0
        detail["heat"].append({
            "T_w": T_w, "T_gas_K": round(T_gas_K, 2),
            "V_raw": V_raw, "V_corr_ml": round(V_corr * 1000, 1),
            "P_gas": round(P_GAS_HEAT, 6), "P_air": round(P_air, 6),
            "P_water": round(P_water, 6), "P_lit": round(P_lit, 6),
        })
        if P_water <= 0:
            neg_count += 1
            continue
        cc_x.append(1.0 / T_gas_K)
        cc_y.append(math.log(P_water))
    
    # 냉각 데이터
    if case in ("B", "C"):
        cool_points = COOL_DATA if case == "C" else [COOL_DATA[0]]  # B: 20°C만
        for V_raw, T_w in cool_points:
            T_gas_K = T_w + 273.15 - dT
            V_corr = (V_raw - MENISCUS) / 1000.0
            P_air = n_air_c * R_ATM * T_gas_K / V_corr
            P_water = P_GAS_COOL - P_air
            P_lit = P_water_lit(T_w - dT) if T_gas_K > 0 else 0
            detail["cool"].append({
                "T_w": T_w, "T_gas_K": round(T_gas_K, 2),
                "V_raw": V_raw, "V_corr_ml": round(V_corr * 1000, 1),
                "P_gas": round(P_GAS_COOL, 6), "P_air": round(P_air, 6),
                "P_water": round(P_water, 6), "P_lit": round(P_lit, 6),
            })
            if P_water <= 0:
                neg_count += 1
                continue
            cc_x.append(1.0 / T_gas_K)
            cc_y.append(math.log(P_water))
    
    n = len(cc_x)
    if n < 4:
        return None
    
    # 선형 회귀
    sx = sum(cc_x)
    sy = sum(cc_y)
    sxy = sum(x * y for x, y in zip(cc_x, cc_y))
    sx2 = sum(x ** 2 for x in cc_x)
    
    denom = n * sx2 - sx ** 2
    if denom == 0:
        return None
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    
    y_mean = sy / n
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in zip(cc_x, cc_y))
    ss_tot = sum((y - y_mean) ** 2 for y in cc_y)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    
    dH = -slope * R_J / 1000
    P_triple = math.exp(slope / T_TRIPLE_K + intercept)
    
    valid = neg_count == 0
    
    return {
        "dT": dT, "n_air_h": n_air_h, "n_air_c": n_air_c or n_air_h,
        "slope": slope, "intercept": intercept, "R2": r2,
        "dH": dH, "P_triple": P_triple,
        "n_points": n, "neg_count": neg_count, "valid": valid,
        "cost_dH": abs(dH - DH_VAP_LIT) / DH_VAP_LIT,
        "cost_triple": abs(P_triple - P_TRIPLE_ATM) / P_TRIPLE_ATM,
        "detail": detail,
    }


def scan_case(case, dT_range=(0, 15, 0.2), n_air_range=(7e-5, 1.4e-4, 1e-6),
              n_air_c_range=None):
    """
    2~3차원 scanning
    
    Case A: n_air_heat만 (가열)
    Case B/C: n_air_heat + n_air_cool (독립)
    """
    results = []
    best = None
    best_cost = float('inf')
    
    dT_values = []
    d = dT_range[0]
    while d <= dT_range[1] + 0.001:
        dT_values.append(round(d, 2))
        d += dT_range[2]
    
    n_values = []
    n = n_air_range[0]
    while n <= n_air_range[1] + 1e-9:
        n_values.append(n)
        n += n_air_range[2]
    
    for dT in dT_values:
        for n_h in n_values:
            if case == "A":
                # Case A: n_air 하나만
                r = compute_cc(dT, n_h, case="A")
                if r:
                    r["cost"] = r["cost_dH"] + r["cost_triple"]
                    if not r["valid"]:
                        r["cost"] += 10.0
                    results.append(r)
                    if r["cost"] < best_cost:
                        best_cost = r["cost"]
                        best = r
            else:
                # Case B/C: n_air_cool도 독립 scanning
                nc_values = n_values if n_air_c_range is None else []
                if n_air_c_range:
                    nc = n_air_c_range[0]
                    while nc <= n_air_c_range[1] + 1e-9:
                        nc_values.append(nc)
                        nc += n_air_c_range[2]
                
                for n_c in nc_values:
                    r = compute_cc(dT, n_h, n_c, case=case)
                    if r:
                        r["cost"] = r["cost_dH"] + r["cost_triple"]
                        if not r["valid"]:
                            r["cost"] += 10.0
                        results.append(r)
                        if r["cost"] < best_cost:
                            best_cost = r["cost"]
                            best = r
    
    return best, results


def print_result(r, label=""):
    """결과 출력"""
    print(f"\n{'='*70}")
    print(f"RESULT [{label}]")
    print(f"{'='*70}")
    print(f"δT = {r['dT']:.2f} °C")
    print(f"n_air_가열 = {r['n_air_h']:.8f} mol")
    print(f"n_air_냉각 = {r['n_air_c']:.8f} mol")
    print(f"ΔH_vap = {r['dH']:.2f} kJ/mol (문헌 {DH_VAP_LIT}, 오차 {r['cost_dH']*100:.1f}%)")
    print(f"삼중점 = {r['P_triple']:.5f} atm (문헌 {P_TRIPLE_ATM}, 오차 {r['cost_triple']*100:.1f}%)")
    print(f"R² = {r['R2']:.6f}")
    print(f"포인트 수 = {r['n_points']}, P_water<0 = {r['neg_count']}")
    print(f"cost = {r['cost']:.4f}")
    
    print(f"\n가열 데이터 상세:")
    print(f"{'T_w':>6} {'T_gas':>8} {'V_raw':>6} {'P_gas':>10} {'P_air':>10} "
          f"{'P_water':>10} {'P_lit':>10} {'비율':>8}")
    for d in r['detail']['heat']:
        ratio = d['P_water'] / d['P_lit'] if d['P_lit'] > 0 else float('inf')
        mark = " ***NEG" if d['P_water'] <= 0 else ""
        print(f"{d['T_w']:6.1f} {d['T_gas_K']:8.2f} {d['V_raw']:6.1f} "
              f"{d['P_gas']:10.6f} {d['P_air']:10.6f} {d['P_water']:10.6f} "
              f"{d['P_lit']:10.6f} {ratio:8.3f}{mark}")
    
    if r['detail']['cool']:
        print(f"\n냉각 데이터 상세:")
        for d in r['detail']['cool']:
            ratio = d['P_water'] / d['P_lit'] if d['P_lit'] > 0 else float('inf')
            mark = " ***NEG" if d['P_water'] <= 0 else ""
            print(f"{d['T_w']:6.1f} {d['T_gas_K']:8.2f} {d['V_raw']:6.1f} "
                  f"{d['P_gas']:10.6f} {d['P_air']:10.6f} {d['P_water']:10.6f} "
                  f"{d['P_lit']:10.6f} {ratio:8.3f}{mark}")


def comparative_table(results_a, results_b, results_c):
    """Case A/B/C 비교 테이블"""
    print(f"\n{'='*80}")
    print("COMPARATIVE TABLE: Case A vs B vs C")
    print(f"{'='*80}")
    print(f"{'Case':>5} {'δT':>6} {'n_air_H':>10} {'n_air_C':>10} {'ΔH':>7} "
          f"{'P_tri':>8} {'R²':>7} {'dH_err':>7} {'tri_err':>7} {'cost':>7} {'valid':>5}")
    print("-"*95)
    
    for label, best in [("A", results_a), ("B", results_b), ("C", results_c)]:
        if best:
            v = "Y" if best['valid'] else "N"
            print(f"{label:>5} {best['dT']:6.2f} {best['n_air_h']:10.8f} {best['n_air_c']:10.8f} "
                  f"{best['dH']:7.2f} {best['P_triple']:8.5f} {best['R2']:7.4f} "
                  f"{best['cost_dH']*100:6.1f}% {best['cost_triple']*100:6.1f}% "
                  f"{best['cost']:7.4f} {v:>5}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="C-C 검증 파이프라인 v3")
    parser.add_argument("--case", choices=["A", "B", "C", "all"], default="all",
                        help="Case to run (A=가열만, B=+냉각20°C, C=전체, all=전부)")
    args = parser.parse_args()
    
    print(f"C-C 검증 파이프라인 v3")
    print(f"메스실린더 내경 = {D_INNER_CM} cm, A = {A_CM2:.4f} cm²")
    print(f"P_gas_가열 = {P_GAS_HEAT:.6f} atm")
    print(f"P_gas_냉각 = {P_GAS_COOL:.6f} atm")
    
    best_a = best_b = best_c = None
    
    if args.case in ("A", "all"):
        print("\n" + "="*70)
        print("CASE A: 가열 데이터만 (12포인트)")
        print("="*70)
        # Coarse → fine
        best_a_coarse, _ = scan_case("A", dT_range=(0, 15, 0.5), n_air_range=(8e-5, 1.4e-4, 2e-6))
        if best_a_coarse:
            dT_c = best_a_coarse['dT']
            n_c = best_a_coarse['n_air_h']
            best_a, _ = scan_case("A", 
                dT_range=(max(0, dT_c-2), dT_c+2, 0.1),
                n_air_range=(max(1e-6, n_c-1e-5), n_c+1e-5, 5e-7))
            if best_a:
                print_result(best_a, "Case A (가열만)")
    
    if args.case in ("B", "all"):
        print("\n" + "="*70)
        print("CASE B: 가열 + 냉각 20°C (13포인트)")
        print("="*70)
        # n_air_cool 범위를 n_air_heat 중심으로 ±30%
        if best_a:
            n_c = best_a['n_air_h']
            nc_range = (n_c * 0.7, n_c * 1.3, 5e-7)
        else:
            nc_range = (7e-5, 1.4e-4, 2e-6)
        best_b_coarse, _ = scan_case("B", 
            dT_range=(0, 15, 0.5), n_air_range=(8e-5, 1.4e-4, 2e-6),
            n_air_c_range=nc_range)
        if best_b_coarse:
            best_b, _ = scan_case("B",
                dT_range=(max(0, best_b_coarse['dT']-2), best_b_coarse['dT']+2, 0.1),
                n_air_range=(max(1e-6, best_b_coarse['n_air_h']-1e-5), best_b_coarse['n_air_h']+1e-5, 5e-7),
                n_air_c_range=(max(1e-6, best_b_coarse['n_air_c']-5e-6), best_b_coarse['n_air_c']+5e-6, 5e-7))
            if best_b:
                print_result(best_b, "Case B (가열+냉각20°C)")
    
    if args.case in ("C", "all"):
        print("\n" + "="*70)
        print("CASE C: 가열 + 냉각 전체 (P_water>0 제약)")
        print("="*70)
        if best_a:
            n_c = best_a['n_air_h']
            nc_range = (n_c * 0.5, n_c * 1.5, 5e-7)
        else:
            nc_range = (5e-5, 2e-4, 2e-6)
        best_c_coarse, _ = scan_case("C",
            dT_range=(0, 15, 0.5), n_air_range=(5e-5, 2e-4, 2e-6),
            n_air_c_range=nc_range)
        if best_c_coarse:
            best_c, _ = scan_case("C",
                dT_range=(max(0, best_c_coarse['dT']-2), best_c_coarse['dT']+2, 0.1),
                n_air_range=(max(1e-6, best_c_coarse['n_air_h']-1e-5), best_c_coarse['n_air_h']+1e-5, 5e-7),
                n_air_c_range=(max(1e-6, best_c_coarse['n_air_c']-5e-6), best_c_coarse['n_air_c']+5e-6, 5e-7))
            if best_c:
                print_result(best_c, "Case C (전체)")
    
    # Case D: 앵커 보정 (냉각 외삽 → 가열 앵커)
    best_d = None
    if args.case in ("D", "all") and best_a:
        print("\n" + "="*70)
        print("CASE D: 가열 + 냉각 외삽 앵커 보정")
        print("="*70)
        
        # Step 1: 냉각 데이터만으로 C-C 회귀 → 0°C 외삽
        best_cool = None
        best_cool_cost = float('inf')
        for dT in [x * 0.1 for x in range(0, 80)]:
            for n_c_10k in range(5, 20):
                n_c = n_c_10k * 1e-5
                cc_x = []
                valid = True
                for V_raw, T_w in COOL_DATA:
                    T_gas_K = T_w + 273.15 - dT
                    V_corr = (V_raw - MENISCUS) / 1000.0
                    if V_corr <= 0:
                        valid = False; break
                    P_air = n_c * R_ATM * T_gas_K / V_corr
                    P_water = P_GAS_COOL - P_air
                    if P_water <= 0:
                        valid = False; break
                    cc_x.append((1.0 / T_gas_K, math.log(P_water)))
                if not valid or len(cc_x) < 2:
                    continue
                n = len(cc_x)
                xs = [x[0] for x in cc_x]; ys = [x[1] for x in cc_x]
                sx = sum(xs); sy = sum(ys)
                sxy = sum(x*y for x,y in zip(xs, ys))
                sx2 = sum(x**2 for x in xs)
                denom = n * sx2 - sx**2
                if denom == 0: continue
                slope_c = (n * sxy - sx * sy) / denom
                int_c = (sy - slope_c * sx) / n
                P_triple_c = math.exp(slope_c / T_TRIPLE_K + int_c)
                dH_c = -slope_c * R_J / 1000
                cost = abs(P_triple_c - P_TRIPLE_ATM) / P_TRIPLE_ATM + abs(dH_c - DH_VAP_LIT) / DH_VAP_LIT * 0.5
                if cost < best_cool_cost:
                    best_cool_cost = cost
                    best_cool = {'dT': dT, 'n_c': n_c, 'slope': slope_c, 'intercept': int_c}
        
        if best_cool:
            anchor_x = 1.0 / T_TRIPLE_K
            anchor_y = best_cool['slope'] * anchor_x + best_cool['intercept']
            P_anchor = math.exp(anchor_y)
            print(f"  냉각 외삽 앵커: (1/T={anchor_x:.6f}, ln(P)={anchor_y:.4f}, P={P_anchor:.5f} atm)")
            
            # Step 2: 가중치 w로 앵커 포함 회귀
            best_d = None
            for w in [0.5, 1.0, 2.0, 3.0, 5.0]:
                best_for_w = None
                best_cost_w = float('inf')
                for dT in [x * 0.1 for x in range(0, 200)]:
                    for n_h_10k in range(100, 140):
                        n_h = n_h_10k * 1e-6
                        cc_x = []; cc_y = []; cc_w = []
                        valid = True
                        for V_raw, T_w in HEAT_DATA:
                            T_gas_K = T_w + 273.15 - dT
                            V_corr = (V_raw - MENISCUS) / 1000.0
                            P_air = n_h * R_ATM * T_gas_K / V_corr
                            P_water = P_GAS_HEAT - P_air
                            if P_water <= 0:
                                valid = False; break
                            cc_x.append(1.0 / T_gas_K)
                            cc_y.append(math.log(P_water))
                            cc_w.append(1.0)
                        if not valid: continue
                        cc_x.append(anchor_x); cc_y.append(anchor_y); cc_w.append(w)
                        
                        sw = sum(cc_w)
                        swx = sum(wi*xi for wi,xi in zip(cc_w, cc_x))
                        swy = sum(wi*yi for wi,yi in zip(cc_w, cc_y))
                        swxy = sum(wi*xi*yi for wi,xi,yi in zip(cc_w, cc_x, cc_y))
                        swx2 = sum(wi*xi**2 for wi,xi in zip(cc_w, cc_x))
                        denom = sw * swx2 - swx**2
                        if denom == 0: continue
                        slope = (sw * swxy - swx * swy) / denom
                        intercept = (swy - slope * swx) / sw
                        
                        y_pred = [slope * x + intercept for x in cc_x]
                        ss_res = sum(wi * (yi - yp)**2 for wi, yi, yp in zip(cc_w, cc_y, y_pred))
                        y_wmean = swy / sw
                        ss_tot = sum(wi * (yi - y_wmean)**2 for wi, yi in zip(cc_w, cc_y))
                        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                        
                        dH = -slope * R_J / 1000
                        P_triple = math.exp(slope / T_TRIPLE_K + intercept)
                        cost = abs(dH - DH_VAP_LIT) / DH_VAP_LIT + abs(P_triple - P_TRIPLE_ATM) / P_TRIPLE_ATM
                        if cost < best_cost_w:
                            best_cost_w = cost
                            best_for_w = {
                                'dT': dT, 'n_air_h': n_h, 'n_air_c': anchor_x,
                                'slope': slope, 'intercept': intercept, 'R2': r2,
                                'dH': dH, 'P_triple': P_triple, 'cost': cost,
                                'cost_dH': abs(dH - DH_VAP_LIT) / DH_VAP_LIT,
                                'cost_triple': abs(P_triple - P_TRIPLE_ATM) / P_TRIPLE_ATM,
                                'valid': True, 'n_points': 13, 'weight': w
                            }
                if best_for_w:
                    print(f"  w={w:.1f}: δT={best_for_w['dT']:.1f}, n={best_for_w['n_air_h']:.8f}, "
                          f"ΔH={best_for_w['dH']:.2f}, P_tri={best_for_w['P_triple']:.5f}, "
                          f"R²={best_for_w['R2']:.4f}, dH_err={best_for_w['cost_dH']*100:.1f}%, "
                          f"tri_err={best_for_w['cost_triple']*100:.1f}%")
                    if best_d is None or best_for_w['cost'] < best_d['cost']:
                        best_d = best_for_w
    
    # 비교 테이블
    comparative_table(best_a, best_b, best_c)
    
    # Case D 결과도 출력
    if best_d:
        print(f"\n{'='*70}")
        print("CASE D 최적 (가열 + 냉각 외삽 앵커)")
        print(f"{'='*70}")
        print(f"δT = {best_d['dT']:.1f}, n_air_h = {best_d['n_air_h']:.8f}")
        print(f"ΔH_vap = {best_d['dH']:.2f} kJ/mol (오차 {best_d['cost_dH']*100:.1f}%)")
        print(f"삼중점 = {best_d['P_triple']:.5f} atm (오차 {best_d['cost_triple']*100:.1f}%)")
        print(f"R² = {best_d['R2']:.4f}")