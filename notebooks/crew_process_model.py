"""
Coupled Wastewater Treatment and Carbon Mineralization Model
Inspired by CREW's process intensification technology.

Models:
1. ASM1 (Activated Sludge Model No. 1, Henze et al. 1987) - biological treatment
2. CO2 mineral sequestration - carbon capture via alkalinity enhancement

Same mathematical framework as Dong et al. (2025) AWL reactor model.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import warnings
warnings.filterwarnings('ignore')

plt.rcParams.update({'figure.dpi': 120, 'font.size': 11})


# ============================================================
# ASM1 KINETIC PARAMETERS (Henze et al. 1987, T=20 C)
# ============================================================

ASM1_PARAMS = {
    'mu_H':          6.0,    # Max specific growth rate, heterotrophs [1/d]
    'K_S':           20.0,   # Half-saturation, substrate [g COD/m3]
    'K_OH':          0.2,    # Half-saturation, O2 (heterotrophs) [g O2/m3]
    'K_NO':          0.5,    # Half-saturation, NO3 [g N/m3]
    'b_H':           0.62,   # Decay rate, heterotrophs [1/d]
    'eta_g':         0.8,    # Anoxic growth correction factor
    'Y_H':           0.67,   # Yield, heterotrophs [g COD biomass / g COD substrate]
    'mu_A':          0.8,    # Max specific growth rate, autotrophs [1/d]
    'K_NH':          1.0,    # Half-saturation, NH4 [g N/m3]
    'K_OA':          0.4,    # Half-saturation, O2 (autotrophs) [g O2/m3]
    'b_A':           0.05,   # Decay rate, autotrophs [1/d]
    'Y_A':           0.24,   # Yield, autotrophs [g COD / g N oxidized]
    'i_XB':          0.086,  # N content of biomass [g N / g COD]
    'f_P':           0.08,   # Fraction of biomass as endogenous residue
    'CO2_per_COD':   1.375,  # g CO2 produced per g COD oxidized
}

MINERAL_PARAMS = {
    'k_min':           0.05,   # Mineralization rate constant [1/d per kg/m3 mineral]
    'KH_CO2':          3.4e-2, # Henry constant CO2 [mol/L/atm] at 20C
    'pCO2_atm':        420e-6, # Atmospheric pCO2 [atm]
    'mineral_loading': 10.0,   # kg mineral / m3 reactor
    'MW_CO2':          44.01,  # g/mol
}

REACTOR_PARAMS = {
    'HRT':     0.5,    # Hydraulic retention time [d] = 12 hours
    'SRT':     10.0,   # Sludge retention time [d]
    'Q_O2':    200.0,  # Aeration rate [g O2/m3/d]
    'SS_in':   200.0,  # Influent COD [g COD/m3]
    'SO_in':   0.0,    # Influent DO [g O2/m3]
    'SNH_in':  40.0,   # Influent ammonium [g N/m3]
    'SNO_in':  0.0,    # Influent nitrate [g N/m3]
    'SALK_in': 7.0,    # Influent alkalinity [mol HCO3/m3]
    'SCO2_in': 0.5,    # Influent dissolved CO2 [g CO2/m3]
    'XBH_in':  0.0,    # Influent heterotroph biomass
    'XBA_in':  0.0,    # Influent autotroph biomass
}


# ============================================================
# KINETIC FUNCTIONS
# ============================================================

def monod(S, K):
    return max(S, 0.0) / (K + max(S, 0.0))

def switching(S, K):
    return K / (K + max(S, 0.0))

def co2_eq_g_m3(mineral_params):
    """Equilibrium dissolved CO2 with atmosphere [g CO2/m3]"""
    CO2_eq_mol_L = mineral_params['KH_CO2'] * mineral_params['pCO2_atm']
    return CO2_eq_mol_L * mineral_params['MW_CO2'] * 1000.0

def co2_mineralization_rate(S_CO2, mineral_loading, mp):
    """Rate of CO2 sequestration by mineral dissolution [g CO2/m3/d]"""
    CO2_eq = co2_eq_g_m3(mp)
    excess = max(S_CO2 - CO2_eq, 0.0)
    return mp['k_min'] * mineral_loading * excess

def asm1_rates(state, p):
    """Calculate ASM1 process rates. Returns (rho1, rho2, rho3, rho4, rho5)."""
    S_S, S_O, S_NH, S_NO, X_BH, X_BA = [max(v, 0.0) for v in state[:6]]

    rho1 = p['mu_H'] * monod(S_S, p['K_S']) * monod(S_O, p['K_OH']) * X_BH
    rho2 = (p['mu_H'] * monod(S_S, p['K_S']) * switching(S_O, p['K_OH'])
            * monod(S_NO, p['K_NO']) * p['eta_g'] * X_BH)
    rho3 = p['mu_A'] * monod(S_NH, p['K_NH']) * monod(S_O, p['K_OA']) * X_BA
    rho4 = p['b_H'] * X_BH
    rho5 = p['b_A'] * X_BA
    return rho1, rho2, rho3, rho4, rho5


# ============================================================
# COUPLED ODE SYSTEM
# ============================================================
# State: [S_S, S_O, S_NH, S_NO, X_BH, X_BA, S_ALK, S_CO2]
# Units: g/m3 (except S_ALK in mol HCO3/m3)

def coupled_odes(t, y, rp, p, mp):
    S_S, S_O, S_NH, S_NO, X_BH, X_BA, S_ALK, S_CO2 = y

    inv_HRT = 1.0 / rp['HRT']
    inv_SRT = 1.0 / rp['SRT']

    rho1, rho2, rho3, rho4, rho5 = asm1_rates(y, p)

    # Dissolved CO2 at atmospheric equilibrium
    CO2_eq = co2_eq_g_m3(mp)
    k_La_CO2 = 5.0  # d-1, CO2 gas transfer coefficient

    # Mineralization rate
    r_min = co2_mineralization_rate(max(S_CO2, 0.0), mp['mineral_loading'], mp)

    # Alkalinity gain from mineralization: 2 mol HCO3 per mol CO2 mineralized
    r_alk_mineral = r_min * 2.0 / mp['MW_CO2']  # mol HCO3/m3/d

    # Mass balances
    dSS  = inv_HRT * (rp['SS_in']  - S_S)  - (1.0/p['Y_H']) * (rho1 + rho2)
    dSO  = inv_HRT * (rp['SO_in']  - S_O)  + rp['Q_O2'] \
           - ((1.0 - p['Y_H']) / p['Y_H']) * rho1 \
           - (4.57 - p['Y_A']) / p['Y_A'] * rho3
    dSNH = inv_HRT * (rp['SNH_in'] - S_NH) \
           - p['i_XB'] * (rho1 + rho2) \
           - (p['i_XB'] + 1.0/p['Y_A']) * rho3
    dSNO = inv_HRT * (rp['SNO_in'] - S_NO) \
           + (1.0/p['Y_A']) * rho3 \
           - ((1.0 - p['Y_H']) / (2.86 * p['Y_H'])) * rho2
    dXBH = inv_HRT * (rp['XBH_in'] - X_BH) - inv_SRT * X_BH + rho1 + rho2 - rho4
    dXBA = inv_HRT * (rp['XBA_in'] - X_BA) - inv_SRT * X_BA + rho3 - rho5
    dSALK = inv_HRT * (rp['SALK_in'] - S_ALK) \
            - 0.07 * rho3 / 14.0 \
            + 0.07 * rho2 / 14.0 \
            + r_alk_mineral
    dSCO2 = inv_HRT * (rp['SCO2_in'] - S_CO2) \
            + p['CO2_per_COD'] * (1.0 - p['Y_H']) * rho1 \
            + 0.15 * rho3 \
            - r_min \
            - k_La_CO2 * (max(S_CO2, 0.0) - CO2_eq)

    return [dSS, dSO, dSNH, dSNO, dXBH, dXBA, dSALK, dSCO2]


# ============================================================
# SIMULATION
# ============================================================

def run_simulation(reactor_params, asm_params, mineral_params,
                   y0=None, t_end=30.0, n_points=3000):
    """Run coupled wastewater + carbon sequestration simulation."""
    rp = reactor_params.copy()
    rp.setdefault('XBH_in', 0.0)
    rp.setdefault('XBA_in', 0.0)

    if y0 is None:
        y0 = [rp['SS_in'], 2.0, rp['SNH_in'], 0.0,
              100.0, 10.0, rp['SALK_in'], 0.5]

    sol = solve_ivp(
        coupled_odes,
        (0, t_end), y0,
        args=(rp, asm_params, mineral_params),
        t_eval=np.linspace(0, t_end, n_points),
        method='LSODA', rtol=1e-6, atol=1e-9
    )
    return sol


# ============================================================
# MAIN DEMONSTRATION
# ============================================================

if __name__ == '__main__':
    print("Running coupled ASM1 + carbon mineralization model...")

    sol = run_simulation(REACTOR_PARAMS, ASM1_PARAMS, MINERAL_PARAMS)
    t = sol.t
    S_S, S_O, S_NH, S_NO, X_BH, X_BA, S_ALK, S_CO2 = sol.y

    # Compute derived quantities
    CO2_eq_val = co2_eq_g_m3(MINERAL_PARAMS)

    rho_arr = np.array([asm1_rates(sol.y[:, i], ASM1_PARAMS) for i in range(len(t))])
    CO2_prod = ASM1_PARAMS['CO2_per_COD'] * (1 - ASM1_PARAMS['Y_H']) * rho_arr[:, 0]
    CO2_min  = np.array([co2_mineralization_rate(max(S_CO2[i], 0.0),
                         MINERAL_PARAMS['mineral_loading'], MINERAL_PARAMS)
                         for i in range(len(t))])
    seq_eff = np.where(CO2_prod > 0, CO2_min / CO2_prod, 0.0)

    # --- Print steady state ---
    print(f"\nSteady-state results (day {t[-1]:.0f}):")
    print(f"  Effluent COD:        {S_S[-1]:.1f} g/m3  (removal: {(1-S_S[-1]/REACTOR_PARAMS['SS_in'])*100:.1f}%)")
    print(f"  Effluent NH4+:       {S_NH[-1]:.2f} g N/m3  (removal: {(1-S_NH[-1]/REACTOR_PARAMS['SNH_in'])*100:.1f}%)")
    print(f"  Effluent NO3-:       {S_NO[-1]:.1f} g N/m3")
    print(f"  Dissolved O2:        {S_O[-1]:.2f} g O2/m3")
    print(f"  Heterotroph biomass: {X_BH[-1]:.0f} g COD/m3")
    print(f"  Nitrifier biomass:   {X_BA[-1]:.1f} g COD/m3")
    print(f"  Alkalinity:          {S_ALK[-1]:.2f} mol HCO3/m3  (gain: {S_ALK[-1]-REACTOR_PARAMS['SALK_in']:.2f})")
    print(f"  Dissolved CO2:       {S_CO2[-1]:.2f} g/m3  (equilibrium: {CO2_eq_val:.2f})")
    print(f"  CO2 sequestration:   {seq_eff[-1]*100:.1f}%  of biogenic CO2")

    # --- Plot ---
    fig, axes = plt.subplots(3, 3, figsize=(14, 10))
    fig.suptitle('Coupled ASM1 + Carbon Mineralization Model\n'
                 'Municipal Wastewater Treatment with CO2 Sequestration',
                 fontsize=13, fontweight='bold')

    panels = [
        (axes[0,0], S_S,  'COD (g/m3)',       'Substrate Removal',   'b',      REACTOR_PARAMS['SS_in']),
        (axes[0,1], S_NH, 'NH4+ (g N/m3)',     'Ammonium',            'orange', REACTOR_PARAMS['SNH_in']),
        (axes[0,2], S_NO, 'NO3- (g N/m3)',     'Nitrate',             'r',      None),
        (axes[1,0], S_O,  'DO (g O2/m3)',      'Dissolved Oxygen',    'cyan',   None),
        (axes[1,1], X_BH, 'X_BH (g COD/m3)',  'Heterotrophs',        'g',      None),
        (axes[1,2], X_BA, 'X_BA (g COD/m3)',  'Nitrifiers',          'm',      None),
        (axes[2,0], S_ALK,'Alk (mol HCO3/m3)','Alkalinity',          'purple', REACTOR_PARAMS['SALK_in']),
        (axes[2,1], S_CO2,'CO2 (g/m3)',        'Dissolved CO2',       'brown',  CO2_eq_val),
    ]

    for ax, data, ylabel, title, color, ref in panels:
        ax.plot(t, data, color=color, lw=2)
        if ref is not None:
            ax.axhline(ref, color=color, ls='--', alpha=0.4, label='Reference')
            ax.legend(fontsize=8)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xlabel('Time (days)')
        ax.grid(alpha=0.3)

    # Sequestration panel
    ax = axes[2,2]
    ax.plot(t, CO2_prod, 'r-', lw=2, label='CO2 produced')
    ax.plot(t, CO2_min,  'g-', lw=2, label='CO2 mineralized')
    ax.set_ylabel('g CO2/m3/d')
    ax.set_title('Carbon Sequestration')
    ax.legend(fontsize=9)
    ax.set_xlabel('Time (days)')
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig('/mnt/user-data/outputs/crew_model_results.png', dpi=150, bbox_inches='tight')
    print("\nPlot saved.")
    plt.show()
