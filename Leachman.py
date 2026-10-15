"""
CALCULADOR / FLASH DE PROPRIEDADES PARA FLUIDO PURO
Modelos: Peng-Robinson (multifluido) e Leachman (H2 normal)

Pares de entrada aceitos
------------------------
    TP   : T [°C]       e P [bar]
    TRHO : T [°C]       e rho [kg/m3]
    PRHO : P [bar]      e rho [kg/m3]
    PH   : P [bar]      e H [kJ/kg]
    TH   : T [°C]       e H [kJ/kg]
    PS   : P [bar]      e S [kJ/(kg K)]
    TS   : T [°C]       e S [kJ/(kg K)]

Propriedades calculadas
-----------------------
    T       [°C]
    P       [bar]
    H       [kJ/kg]
    S       [kJ/(kg K)]
    Cp      [kJ/(kg K)]
    mu_JT   [K/bar]
    rho     [kg/m3]
    fase
    título x, quando houver duas fases

Regra para a região de saturação
--------------------------------
Se o par de entrada cair na região líquido-vapor, a PRIMEIRA variável
é mantida como variável especificada:

    PH, PS, PRHO -> P é fixada e calcula-se Tsat(P)
    TH, TS, TRHO -> T é fixada e calcula-se Psat(T)

A segunda variável (H, S ou rho) determina o título x quando o estado
cai dentro do envelope líquido-vapor. Para os pares com densidade,
o título é obtido pela regra da alavanca aplicada ao volume específico:

    v = (1-x) vL + x vV

Para TP exatamente sobre a curva de saturação, T e P não definem o título.
O programa sinaliza "saturação subespecificada" e fornece os limites
saturados, sem inventar uma qualidade.

IMPORTANTE SOBRE Cp E muJT NA REGIÃO BIFÁSICA
---------------------------------------------
Para 0 < x < 1, Cp=(dH/dT)_P não é uma propriedade finita de fase
homogênea; em equilíbrio bifásico a temperatura permanece Tsat enquanto
muda o título. O programa reporta Cp = inf.

O coeficiente Joule-Thomson na região bifásica é calculado numericamente
como (dT/dP)_H, realizando dois flashes PH próximos ao estado especificado.

REFERÊNCIAS DE H E S
--------------------
LEACHMAN:
    Usa a referência nativa da formulação de Leachman implementada no
    código (mesma base do script P-h fornecido anteriormente).

PENG-ROBINSON:
    Para cada fluido, aplica-se um deslocamento de referência de modo que
    H = 0 kJ/kg e S = 0 kJ/(kg K) em 25 °C e 1 bar.
    Portanto, valores de H e S usados como entrada com PR devem obedecer
    a essa mesma referência.

Dependências
------------
    numpy
    scipy

Observação
----------
Leachman está disponível apenas para H2 normal neste arquivo.
Peng-Robinson pode ser usado com qualquer fluido puro cadastrado no banco.
"""

import numpy as np
from scipy.optimize import brentq
#from functools import lru_cache
from pathlib import Path


# =============================================================================
# 1. CONFIGURAÇÃO DO USUÁRIO
# =============================================================================

MODELO = "LEACHMAN"       # "PR" ou "LEACHMAN"
FLUIDO = "H2"             # para LEACHMAN somente "H2" disponível.
PAR_ENTRADA = "TP"        # TP, TRHO, PRHO, PH, TH, PS ou TS

# Exemplos de entrada:
# TP   -> VALOR1 = T [°C] ; VALOR2 = P [bar]
# TRHO -> VALOR1 = T [°C] ; VALOR2 = rho [kg/m3]
# PRHO -> VALOR1 = P [bar]; VALOR2 = rho [kg/m3]
# PH   -> VALOR1 = P [bar]; VALOR2 = H [kJ/kg]
# TH   -> VALOR1 = T [°C] ; VALOR2 = H [kJ/kg]
# PS   -> VALOR1 = P [bar]; VALOR2 = S [kJ/(kg K)]
# TS   -> VALOR1 = T [°C] ; VALOR2 = S [kJ/(kg K)]

VALOR1 = 25
VALOR2 = 20

# Limites numéricos do flash
P_MIN_BAR = 0.01
P_MAX_BAR = 2000.0
T_MAX_K_PR = 2000.0

# Passos numéricos para Cp e muJT
DT_NUM = 0.02             # K
DP_JT_BAR = 0.01          # bar (ajustado automaticamente em alta pressão)

# Tolerância para identificar TP sobre a saturação
SAT_REL_TOL = 2.0e-5


# =============================================================================
# 2. BANCO DE FLUIDOS PARA PENG-ROBINSON
# =============================================================================

# Formato:
# [nome, Tc(K), Pc(bar), omega, M(kg/mol), Cp_ig(J/mol/K)]
#
# Cp_ig é assumido constante nesta versão didática. Para faixas amplas de
# temperatura, substitua por uma correlação Cp_ig(T).

FLUIDOS_PR = [
    ["H2",    33.18,   13.00, -0.219, 0.00201588, 28.84],
    ["CH4",  190.564,  45.99,  0.011, 0.01604246, 35.69],
    ["CO2",  304.18,   73.80,  0.224, 0.04400950, 37.14],
    ["C2H6", 305.32,   48.70,  0.099, 0.03006904, 52.49],
    ["C3H8", 369.80,   42.50,  0.152, 0.04409562, 73.60],
]

BANCO_PR = {
    item[0].upper(): {
        "nome": item[0],
        "Tc": item[1],
        "Pc_bar": item[2],
        "Pc": item[2] * 1.0e5,
        "omega": item[3],
        "M": item[4],
        "Cp_ig": item[5],
    }
    for item in FLUIDOS_PR
}


# =============================================================================
# 3. CONSTANTES E CONVERSÕES
# =============================================================================

R = 8.314462618            # J mol^-1 K^-1
T_REF_PR = 298.15          # K
P_REF_PR = 1.0e5           # Pa


def C_para_K(T_C):
    return T_C + 273.15


def K_para_C(T_K):
    return T_K - 273.15


# =============================================================================
# 4. UTILIDADES NUMÉRICAS
# =============================================================================


def _finite(x):
    return x is not None and np.isfinite(x)


def root_scan(func, lo, hi, n=240, log=False):
    """Busca uma mudança de sinal e resolve a raiz por Brent."""
    if not (_finite(lo) and _finite(hi)) or hi <= lo:
        raise ValueError("Intervalo inválido para busca de raiz.")

    if log and lo > 0:
        xs = np.geomspace(lo, hi, n)
    else:
        xs = np.linspace(lo, hi, n)

    x_prev = None
    f_prev = None

    for x in xs:
        try:
            fx = func(x)
        except Exception:
            continue

        if not np.isfinite(fx):
            continue

        if abs(fx) < 1.0e-10:
            return x

        if f_prev is not None and fx * f_prev < 0.0:
            return brentq(func, x_prev, x, maxiter=200)

        x_prev = x
        f_prev = fx

    raise ValueError(
        f"Não foi possível localizar uma raiz no intervalo [{lo:g}, {hi:g}]."
    )


# =============================================================================
# 5. MODELO PENG-ROBINSON
# =============================================================================


def pr_parametros(T, props):
    Tc = props["Tc"]
    Pc = props["Pc"]
    omega = props["omega"]

    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    sqrt_Tr = np.sqrt(T / Tc)
    g = 1.0 + kappa * (1.0 - sqrt_Tr)
    alpha = g**2

    ac = 0.45724 * R**2 * Tc**2 / Pc
    a = ac * alpha
    b = 0.07780 * R * Tc / Pc

    dalpha_dT = -kappa * g / (Tc * sqrt_Tr)
    da_dT = ac * dalpha_dT

    return a, b, da_dT


def pr_roots(T, P, props):
    a, b, _ = pr_parametros(T, props)

    A = a * P / (R**2 * T**2)
    B = b * P / (R * T)

    coef = [
        1.0,
        -(1.0 - B),
        A - 3.0 * B**2 - 2.0 * B,
        -(A * B - B**2 - B**3),
    ]

    roots = np.roots(coef)
    roots = roots[np.abs(roots.imag) < 1.0e-8].real
    roots = roots[roots > B + 1.0e-12]
    roots.sort()

    return roots, A, B


def pr_lnphi(Z, A, B):
    sqrt2 = np.sqrt(2.0)
    if B < 1.0e-14:
        return 0.0

    return (
        Z - 1.0
        - np.log(Z - B)
        - A / (2.0 * sqrt2 * B)
        * np.log(
            (Z + (1.0 + sqrt2) * B)
            / (Z + (1.0 - sqrt2) * B)
        )
    )


def pr_psat(T, props):
    """Pressão de saturação [Pa] pela igualdade de fugacidades."""
    Tc = props["Tc"]
    Pc = props["Pc"]

    if T >= Tc or T <= 0.0:
        return None

    P_scan = np.geomspace(max(10.0, 1.0e-8 * Pc), 0.999999 * Pc, 600)
    prev = None
    prevP = None

    for P in P_scan:
        roots, A, B = pr_roots(T, P, props)
        if len(roots) == 3:
            F = pr_lnphi(roots[0], A, B) - pr_lnphi(roots[-1], A, B)
            if prev is not None and F * prev < 0.0:
                return brentq(
                    lambda p: (
                        lambda rr_ab: pr_lnphi(rr_ab[0][0], rr_ab[1], rr_ab[2])
                        - pr_lnphi(rr_ab[0][-1], rr_ab[1], rr_ab[2])
                    )(pr_roots(T, p, props)),
                    prevP,
                    P,
                    maxiter=120,
                )
            prev = F
            prevP = P
        else:
            prev = None
            prevP = None

    return None


def pr_tsat(P, props):
    """Temperatura de saturação [K] para P<Pc."""
    if P >= props["Pc"] or P <= 0.0:
        return None

    Tc = props["Tc"]
    Tmin = max(1.0, 0.25 * Tc)
    Tmax = 0.99999 * Tc

    def f(T):
        ps = pr_psat(T, props)
        if ps is None:
            return np.nan
        return ps - P

    try:
        return root_scan(f, Tmin, Tmax, n=160)
    except Exception:
        return None


def pr_root_for_phase(T, P, props, phase):
    roots, A, B = pr_roots(T, P, props)
    if len(roots) == 0:
        raise ValueError("Nenhuma raiz física da PR encontrada.")

    if phase == "liquido":
        Z = roots[0]
    else:
        Z = roots[-1]

    return Z, A, B


def pr_departures(T, P, props, phase):
    """H^R [J/mol], S^R [J/mol/K], V [m3/mol]."""
    a, b, da_dT = pr_parametros(T, props)
    Z, _, B = pr_root_for_phase(T, P, props, phase)
    sqrt2 = np.sqrt(2.0)

    log_term = np.log(
        (Z + (1.0 + sqrt2) * B)
        / (Z + (1.0 - sqrt2) * B)
    )

    Hres = (
        R * T * (Z - 1.0)
        + (T * da_dT - a) / (2.0 * sqrt2 * b) * log_term
    )

    Sres = (
        R * np.log(Z - B)
        + da_dT / (2.0 * sqrt2 * b) * log_term
    )

    V = Z * R * T / P

    return Hres, Sres, V


def pr_phase_stable(T, P, props):
    if T >= props["Tc"]:
        if P >= props["Pc"]:
            return "condicao_supercritica"
        return "gas"

    ps = pr_psat(T, props)
    if ps is None:
        # fallback pela raiz mais externa
        return "vapor"

    if abs(P - ps) / ps <= SAT_REL_TOL:
        return "saturacao"

    return "vapor" if P < ps else "liquido"


def pr_raw_hs(T, P, props, phase):
    Cp = props["Cp_ig"]
    Hres, Sres, V = pr_departures(T, P, props, phase)

    h_ig = Cp * (T - T_REF_PR)
    s_ig = Cp * np.log(T / T_REF_PR) - R * np.log(P / P_REF_PR)

    return h_ig + Hres, s_ig + Sres, V


def pr_reference_offsets(props):
    phase = pr_phase_stable(T_REF_PR, P_REF_PR, props)
    if phase == "saturacao":
        phase = "vapor"
    h0, s0, _ = pr_raw_hs(T_REF_PR, P_REF_PR, props, phase)
    return h0, s0


def pr_hs_phase(T, P, props, phase):
    h_raw, s_raw, V = pr_raw_hs(T, P, props, phase)
    h0, s0 = pr_reference_offsets(props)

    h = (h_raw - h0) / props["M"] / 1000.0
    s = (s_raw - s0) / props["M"] / 1000.0
    rho_mass = props["M"] / V

    return h, s, V, rho_mass


def pr_sat_state_T(T, props):
    ps = pr_psat(T, props)
    if ps is None:
        raise ValueError("PR: saturação não disponível para esta temperatura.")

    hL, sL, VL, rhoL = pr_hs_phase(T, ps, props, "liquido")
    hV, sV, VV, rhoV = pr_hs_phase(T, ps, props, "vapor")

    return {
        "T_K": T, "P_Pa": ps,
        "hL": hL, "hV": hV, "sL": sL, "sV": sV,
        "rhoL": rhoL, "rhoV": rhoV,
    }


def pr_sat_state_P(P, props):
    T = pr_tsat(P, props)
    if T is None:
        raise ValueError("PR: saturação não disponível para esta pressão.")
    return pr_sat_state_T(T, props)


def pr_basic_TP(T, P, props, force_phase=None):
    phase = force_phase or pr_phase_stable(T, P, props)
    if phase == "saturacao":
        return {"ambiguous_sat": True, "sat": pr_sat_state_T(T, props)}

    h, s, V, rho = pr_hs_phase(T, P, props, phase)
    return {
        "T_K": T,
        "P_Pa": P,
        "H": h,
        "S": s,
        "V_molar": V,
        "rho": rho,
        "phase": phase,
        "x": None,
        "ambiguous_sat": False,
    }


def pr_cp_mu(T, P, props, phase):
    """Cp [kJ/(kg K)] e muJT [K/bar], numericamente."""
    dT = max(DT_NUM, 1.0e-5 * T)

    hp = pr_hs_phase(T + dT, P, props, phase)[0]
    hm = pr_hs_phase(T - dT, P, props, phase)[0]
    Cp = (hp - hm) / (2.0 * dT)

    V = pr_hs_phase(T, P, props, phase)[2]
    Vp = pr_hs_phase(T + dT, P, props, phase)[2]
    Vm = pr_hs_phase(T - dT, P, props, phase)[2]
    dVdT = (Vp - Vm) / (2.0 * dT)

    # Cp: kJ/kg/K -> J/mol/K
    Cp_molar = Cp * 1000.0 * props["M"]
    mu_Pa = (T * dVdT - V) / Cp_molar
    mu_bar = mu_Pa * 1.0e5

    return Cp, mu_bar


# =============================================================================
# 6. MODELO LEACHMAN PARA H2 NORMAL
# =============================================================================

# Constantes da formulação utilizada no script P-h anexado
R_L = 8.314472
M_L = 0.00201588
TC_L = 33.145
RHOC_L = 15508.0
PC_L = 1_296_400.0
TTRIPLE_L = 13.957
PTRIPLE_L = 7357.828141607672
TMAX_L = 1000.0

# Parte ideal
A1_L = -1.4579856475
A2_L = 1.888076782
ALOG_L = 1.5
NPE_L = np.array([1.616, -0.4117, -0.792, 0.758, 1.217])
VPE_L = np.array([531.0, 751.0, 1989.0, 2484.0, 6859.0])

# Residual potência
D_L = np.array([1, 4, 1, 1, 2, 2, 3, 1, 3], dtype=float)
L_L = np.array([0, 0, 0, 0, 0, 0, 0, 1, 1], dtype=float)
N_L = np.array([-6.93643, 0.01, 2.1101, 4.52059, 0.732564,
                -1.34086, 0.130985, -0.777414, 0.351944])
TEXP_L = np.array([0.6844, 1.0, 0.989, 0.489, 0.803,
                   1.1444, 1.409, 1.754, 1.311])

# Residual gaussianos
BETA_L = np.array([0.171, 0.2245, 0.1304, 0.2785, 0.3967])
DG_L = np.array([2, 1, 3, 1, 1], dtype=float)
EPS_L = np.array([1.506, 0.156, 1.736, 0.67, 1.662])
ETA_L = np.array([1.685, 0.489, 0.103, 2.506, 1.607])
GAMMA_L = np.array([0.7164, 1.3444, 1.4517, 0.7204, 1.5445])
NG_L = np.array([-0.0211716, 0.0226312, 0.032187, -0.0231752, 0.0557346])
TG_L = np.array([4.187, 5.646, 0.791, 7.249, 2.986])

# Ancillaries de saturação
PSN_L = np.array([-0.5553870929135406, -7.0852972599030295,
                  0.006868538454095541, 0.8295514415153659,
                  -0.08773773693832483, 3.754274235303921])
PST_L = np.array([0.882, 1.072, 1.355, 2.61, 8.191, 1.239])

RHOLN_L = np.array([1.9546584452513824, -5.715534926385584,
                    714.2366728890577, -1666.5447521212081,
                    2698.8598951804593, -12755.576112629806])
RHOLT_L = np.array([0.407, 3.401, 7.095, 8.095, 10.911, 17.353])

RHOVN_L = np.array([-3.036644018774332, -773.9211864728516,
                    1530.809493601263, -3087.995404817007,
                    3472.491545729249, -15649.585518005188])
RHOVT_L = np.array([0.473, 6.618, 7.083, 9.094, 10.577, 19.579])


def l_alpha0(tau, delta):
    x = VPE_L * tau / TC_L
    return (
        np.log(delta) + A1_L + A2_L * tau + ALOG_L * np.log(tau)
        + np.sum(NPE_L * np.log(1.0 - np.exp(-x)))
    )


def l_alpha0_tau(tau):
    x = VPE_L * tau / TC_L
    return (
        A2_L + ALOG_L / tau
        + np.sum(
            NPE_L * (VPE_L / TC_L) * np.exp(-x) / (1.0 - np.exp(-x))
        )
    )


def l_alphar_derivatives(tau, delta):
    expo_power = np.where(L_L == 0.0, 0.0, -delta**L_L)
    terms = N_L * delta**D_L * tau**TEXP_L * np.exp(expo_power)

    ar = np.sum(terms)
    ar_delta = np.sum(
        terms * (
            D_L / delta
            - np.where(L_L == 0.0, 0.0, L_L * delta**(L_L - 1.0))
        )
    )
    ar_tau = np.sum(terms * TEXP_L / tau)

    terms_g = (
        NG_L * delta**DG_L * tau**TG_L
        * np.exp(-ETA_L * (delta - EPS_L)**2 - BETA_L * (tau - GAMMA_L)**2)
    )

    ar += np.sum(terms_g)
    ar_delta += np.sum(terms_g * (DG_L / delta - 2.0 * ETA_L * (delta - EPS_L)))
    ar_tau += np.sum(terms_g * (TG_L / tau - 2.0 * BETA_L * (tau - GAMMA_L)))

    return ar, ar_delta, ar_tau


def l_properties_Trho(T, rho):
    tau = TC_L / T
    delta = rho / RHOC_L

    a0 = l_alpha0(tau, delta)
    a0_tau = l_alpha0_tau(tau)
    ar, ar_delta, ar_tau = l_alphar_derivatives(tau, delta)

    P = rho * R_L * T * (1.0 + delta * ar_delta)
    h_molar = R_L * T * (
        1.0 + tau * (a0_tau + ar_tau) + delta * ar_delta
    )
    s_molar = R_L * (
        tau * (a0_tau + ar_tau) - a0 - ar
    )

    return P, h_molar, s_molar


def l_sat_ancillary(T):
    if not (TTRIPLE_L <= T < TC_L):
        raise ValueError("Leachman: temperatura fora da região líquido-vapor.")

    theta = 1.0 - T / TC_L
    Psat = PC_L * np.exp(TC_L / T * np.sum(PSN_L * theta**PST_L))
    rhoL = RHOC_L * (1.0 + np.sum(RHOLN_L * theta**RHOLT_L))
    rhoV = RHOC_L * np.exp(TC_L / T * np.sum(RHOVN_L * theta**RHOVT_L))

    return Psat, rhoL, rhoV


def l_sat_state_T(T):
    Psat, rhoL, rhoV = l_sat_ancillary(T)
    _, hL, sL = l_properties_Trho(T, rhoL)
    _, hV, sV = l_properties_Trho(T, rhoV)

    return {
        "T_K": T,
        "P_Pa": Psat,
        "hL": hL / M_L / 1000.0,
        "hV": hV / M_L / 1000.0,
        "sL": sL / M_L / 1000.0,
        "sV": sV / M_L / 1000.0,
        "rhoL": rhoL * M_L,
        "rhoV": rhoV * M_L,
    }


def l_tsat(P):
    if P < PTRIPLE_L or P >= PC_L:
        return None
    try:
        return brentq(lambda T: l_sat_ancillary(T)[0] - P,
                      TTRIPLE_L, TC_L * (1.0 - 1.0e-10), maxiter=150)
    except Exception:
        return None


def l_sat_state_P(P):
    T = l_tsat(P)
    if T is None:
        raise ValueError("Leachman: saturação não disponível para esta pressão.")
    return l_sat_state_T(T)


def l_solve_density_phase(T, P, phase=None):
    if T <= 0.0 or T > TMAX_L:
        raise ValueError("Leachman: temperatura fora da faixa implementada.")

    def f(rho):
        return l_properties_Trho(T, rho)[0] - P

    if T < TC_L and T >= TTRIPLE_L:
        Psat, rhoLsat, rhoVsat = l_sat_ancillary(T)

        if phase is None:
            phase = "vapor" if P < Psat else "liquido"

        if phase == "vapor":
            lo = max(1.0e-12, P / (R_L * T) * 1.0e-6)
            hi = max(rhoVsat * 1.000001, lo * 10.0)
            return brentq(f, lo, hi, maxiter=200)

        if phase == "liquido":
            lo = rhoLsat * 0.999999
            hi = 120000.0
            return brentq(f, lo, hi, maxiter=200)

    # supercrítico / vapor fora da região de saturação
    lo = max(1.0e-12, P / (R_L * T) * 1.0e-7)
    hi = 120000.0
    return brentq(f, lo, hi, maxiter=250)


def l_phase_stable(T, P):
    if T >= TC_L:
        if P >= PC_L:
            return "condicao_supercritica"
        return "gas"
    if T < TTRIPLE_L:
        raise ValueError("Leachman: região abaixo do ponto triplo não tratada como fluido estável.")

    ps = l_sat_ancillary(T)[0]
    if abs(P - ps) / ps <= SAT_REL_TOL:
        return "saturacao"
    return "vapor" if P < ps else "liquido"


def l_hs_phase(T, P, phase):
    rho_molar = l_solve_density_phase(T, P, phase)
    _, h, s = l_properties_Trho(T, rho_molar)
    V_molar = 1.0 / rho_molar
    rho_mass = rho_molar * M_L

    return h / M_L / 1000.0, s / M_L / 1000.0, V_molar, rho_mass


def l_basic_TP(T, P, force_phase=None):
    phase = force_phase or l_phase_stable(T, P)
    if phase == "saturacao":
        return {"ambiguous_sat": True, "sat": l_sat_state_T(T)}

    h, s, V, rho = l_hs_phase(T, P, phase)
    return {
        "T_K": T,
        "P_Pa": P,
        "H": h,
        "S": s,
        "V_molar": V,
        "rho": rho,
        "phase": phase,
        "x": None,
        "ambiguous_sat": False,
    }


def l_cp_mu(T, P, phase):
    dT = max(DT_NUM, 1.0e-5 * T)

    hp = l_hs_phase(T + dT, P, phase)[0]
    hm = l_hs_phase(T - dT, P, phase)[0]
    Cp = (hp - hm) / (2.0 * dT)

    V = l_hs_phase(T, P, phase)[2]
    Vp = l_hs_phase(T + dT, P, phase)[2]
    Vm = l_hs_phase(T - dT, P, phase)[2]
    dVdT = (Vp - Vm) / (2.0 * dT)

    Cp_molar = Cp * 1000.0 * M_L
    mu_Pa = (T * dVdT - V) / Cp_molar
    mu_bar = mu_Pa * 1.0e5

    return Cp, mu_bar


# =============================================================================
# 7. INTERFACE COMUM ENTRE OS MODELOS
# =============================================================================


def get_model_data(modelo, fluido):
    modelo = modelo.upper()
    fluido = fluido.upper()

    if modelo == "PR":
        if fluido not in BANCO_PR:
            raise ValueError(f"Fluido '{fluido}' não cadastrado para PR.")
        return modelo, BANCO_PR[fluido]

    if modelo == "LEACHMAN":
        if fluido != "H2":
            raise ValueError("Leachman nesta versão está disponível apenas para H2 normal.")
        return modelo, None

    raise ValueError("MODELO deve ser 'PR' ou 'LEACHMAN'.")


def model_limits(modelo, props):
    if modelo == "LEACHMAN":
        return TTRIPLE_L + 1.0e-5, TMAX_L, P_MIN_BAR * 1e5, P_MAX_BAR * 1e5

    Tmin = max(1.0, 0.25 * props["Tc"])
    return Tmin, T_MAX_K_PR, P_MIN_BAR * 1e5, P_MAX_BAR * 1e5


def model_Tc_Pc(modelo, props):
    if modelo == "LEACHMAN":
        return TC_L, PC_L
    return props["Tc"], props["Pc"]


def sat_T(modelo, T, props):
    if modelo == "LEACHMAN":
        return l_sat_state_T(T)
    return pr_sat_state_T(T, props)


def sat_P(modelo, P, props):
    if modelo == "LEACHMAN":
        return l_sat_state_P(P)
    return pr_sat_state_P(P, props)


def basic_TP(modelo, T, P, props, force_phase=None):
    if modelo == "LEACHMAN":
        return l_basic_TP(T, P, force_phase)
    return pr_basic_TP(T, P, props, force_phase)


def cp_mu_TP(modelo, T, P, props, phase):
    if modelo == "LEACHMAN":
        return l_cp_mu(T, P, phase)
    return pr_cp_mu(T, P, props, phase)


def pressure_Trho(modelo, T, rho_mass, props):
    """
    Calcula P [Pa] a partir de T [K] e rho [kg/m3].

    Para PR, usa diretamente a forma explícita P(T,V).
    Para Leachman, usa a equação fundamental de Helmholtz.
    """
    if rho_mass <= 0.0:
        raise ValueError("A densidade deve ser positiva.")

    if modelo == "LEACHMAN":
        rho_molar = rho_mass / M_L
        return l_properties_Trho(T, rho_molar)[0]

    # Peng-Robinson
    M = props["M"]
    V = M / rho_mass
    a, b, _ = pr_parametros(T, props)

    if V <= b:
        raise ValueError(
            "PR: densidade especificada implica V <= b; estado não físico."
        )

    denom = V * (V + b) + b * (V - b)

    return (
        R * T / (V - b)
        - a / denom
    )


def quality_from_density(rho, rhoL, rhoV):
    """
    Calcula o título x a partir da densidade mássica usando
    a regra da alavanca no volume específico:

        v = (1-x) vL + x vV
    """
    if rho <= 0.0 or rhoL <= 0.0 or rhoV <= 0.0:
        raise ValueError("Densidades devem ser positivas.")

    v = 1.0 / rho
    vL = 1.0 / rhoL
    vV = 1.0 / rhoV

    if abs(vV - vL) < 1.0e-30:
        raise ValueError("Não foi possível calcular x: vV ≈ vL.")

    return (v - vL) / (vV - vL)


# =============================================================================
# 8. CONSTRUÇÃO DE ESTADO BIFÁSICO
# =============================================================================


def two_phase_from_sat(sat, x, specified=None):
    x = float(np.clip(x, 0.0, 1.0))

    H = sat["hL"] + x * (sat["hV"] - sat["hL"])
    S = sat["sL"] + x * (sat["sV"] - sat["sL"])
    rho = 1.0 / (
        (1.0 - x) / sat["rhoL"] + x / sat["rhoV"]
    )

    if x <= 1.0e-10:
        phase = "liquido_saturado"
    elif x >= 1.0 - 1.0e-10:
        phase = "vapor_saturado"
    else:
        phase = "bifasico"

    return {
        "T_K": sat["T_K"],
        "P_Pa": sat["P_Pa"],
        "H": H,
        "S": S,
        "rho": rho,
        "phase": phase,
        "x": x,
        "ambiguous_sat": False,
        "specified": specified,
        "sat": sat,
    }


# =============================================================================
# 9. FLASH BÁSICO (SEM Cp E muJT)
# =============================================================================


def flash_core(modelo, fluido, par, valor1, valor2):
    modelo, props = get_model_data(modelo, fluido)
    par = par.upper().replace("-", "").replace("_", "")

    if par not in {"TP", "TRHO", "PRHO", "PH", "TH", "PS", "TS"}:
        raise ValueError("Par de entrada deve ser TP, TRHO, PRHO, PH, TH, PS ou TS.")

    Tmin, Tmax, Pmin, Pmax = model_limits(modelo, props)
    Tc, Pc = model_Tc_Pc(modelo, props)

    # -------------------------------------------------------------------------
    # TP
    # -------------------------------------------------------------------------
    if par == "TP":
        T = C_para_K(valor1)
        P = valor2 * 1.0e5

        if not (Tmin <= T <= Tmax):
            raise ValueError("Temperatura fora dos limites numéricos do modelo.")
        if not (Pmin <= P <= Pmax):
            raise ValueError("Pressão fora dos limites numéricos do modelo.")

        st = basic_TP(modelo, T, P, props)
        st["specified"] = "T"
        return st

    # -------------------------------------------------------------------------
    # TRHO: temperatura é especificada; calcula P
    # -------------------------------------------------------------------------
    if par == "TRHO":
        T = C_para_K(valor1)
        rho_target = valor2

        if not (Tmin <= T <= Tmax):
            raise ValueError("Temperatura fora dos limites numéricos do modelo.")
        if rho_target <= 0.0:
            raise ValueError("A densidade deve ser positiva.")

        sat = None
        force_phase = None

        if T < Tc:
            try:
                sat = sat_T(modelo, T, props)
            except Exception:
                sat = None

        if sat is not None:
            rhoL = sat["rhoL"]
            rhoV = sat["rhoV"]
            tol_rho = 1.0e-9 * max(1.0, abs(rhoL), abs(rhoV))

            # Dentro do envelope líquido-vapor: T é mantida fixa
            if rhoV - tol_rho <= rho_target <= rhoL + tol_rho:
                x = quality_from_density(rho_target, rhoL, rhoV)
                st = two_phase_from_sat(sat, x, specified="T")
                st["rho"] = rho_target
                return st

            force_phase = "vapor" if rho_target < rhoV else "liquido"

        P = pressure_Trho(modelo, T, rho_target, props)

        if not np.isfinite(P) or P <= 0.0:
            raise ValueError(
                "O par T-rho fornecido resultou em pressão não física."
            )
        if not (Pmin <= P <= Pmax):
            raise ValueError(
                f"Pressão calculada ({P/1e5:.6g} bar) fora dos limites numéricos."
            )

        st = basic_TP(
            modelo,
            T,
            P,
            props,
            force_phase=force_phase,
        )
        st["specified"] = "T"
        # Preserva exatamente a densidade especificada
        st["rho"] = rho_target
        return st

    # -------------------------------------------------------------------------
    # PRHO: pressão é especificada; calcula T
    # -------------------------------------------------------------------------
    if par == "PRHO":
        P = valor1 * 1.0e5
        rho_target = valor2

        if not (Pmin <= P <= Pmax):
            raise ValueError("Pressão fora dos limites numéricos do modelo.")
        if rho_target <= 0.0:
            raise ValueError("A densidade deve ser positiva.")

        sat = None

        if P < Pc:
            try:
                sat = sat_P(modelo, P, props)
            except Exception:
                sat = None

        if sat is not None:
            rhoL = sat["rhoL"]
            rhoV = sat["rhoV"]
            tol_rho = 1.0e-9 * max(1.0, abs(rhoL), abs(rhoV))

            # Dentro do envelope: P é mantida fixa e T = Tsat(P)
            if rhoV - tol_rho <= rho_target <= rhoL + tol_rho:
                x = quality_from_density(rho_target, rhoL, rhoV)
                st = two_phase_from_sat(sat, x, specified="P")
                st["rho"] = rho_target
                return st

            if rho_target < rhoV:
                # vapor superaquecido: T > Tsat
                Tlo = sat["T_K"] * (1.0 + 2.0e-7)
                Thi = Tmax
                force_phase = "vapor"
            else:
                # líquido comprimido: T < Tsat
                Tlo = Tmin
                Thi = sat["T_K"] * (1.0 - 2.0e-7)
                force_phase = "liquido"

            def f_rho_T(T):
                return pressure_Trho(modelo, T, rho_target, props) - P

            Tsol = root_scan(f_rho_T, Tlo, Thi, n=360)
            st = basic_TP(
                modelo,
                Tsol,
                P,
                props,
                force_phase=force_phase,
            )
            st["specified"] = "P"
            st["rho"] = rho_target
            return st

        # Sem saturação possível: resolve diretamente P(T,rho)=P_especificada
        def f_rho_T(T):
            return pressure_Trho(modelo, T, rho_target, props) - P

        Tsol = root_scan(f_rho_T, Tmin, Tmax, n=500)
        st = basic_TP(modelo, Tsol, P, props)

        if st.get("ambiguous_sat", False):
            raise ValueError(
                "Solução P-rho caiu numericamente sobre a saturação; "
                "use a densidade para definir o título dentro do envelope."
            )

        st["specified"] = "P"
        st["rho"] = rho_target
        return st

    # -------------------------------------------------------------------------
    # PH ou PS: pressão é a PRIMEIRA variável e permanece especificada
    # -------------------------------------------------------------------------
    if par in {"PH", "PS"}:
        P = valor1 * 1.0e5
        target = valor2
        key = "H" if par == "PH" else "S"

        if not (Pmin <= P <= Pmax):
            raise ValueError("Pressão fora dos limites numéricos do modelo.")

        # Tenta primeiro a região de saturação a esta pressão
        sat = None
        if P < Pc:
            try:
                sat = sat_P(modelo, P, props)
            except Exception:
                sat = None

        if sat is not None:
            low = sat["hL"] if key == "H" else sat["sL"]
            high = sat["hV"] if key == "H" else sat["sV"]

            tol = 1.0e-9 * max(1.0, abs(low), abs(high))
            if low - tol <= target <= high + tol:
                x = (target - low) / (high - low)
                st = two_phase_from_sat(sat, x, specified="P")
                # preserva exatamente a variável de entrada
                st[key] = target
                return st

            # Abaixo do líquido saturado: busca T<Tsat
            if target < low:
                Thi = sat["T_K"] * (1.0 - 2.0e-7)
                Tlo = Tmin
                force_phase = "liquido"
            else:
                # Acima do vapor saturado: busca T>Tsat
                Tlo = sat["T_K"] * (1.0 + 2.0e-7)
                Thi = Tmax
                force_phase = "vapor"

            def f(T):
                st = basic_TP(modelo, T, P, props, force_phase=force_phase)
                return st[key] - target

            Tsol = root_scan(f, Tlo, Thi, n=280)
            st = basic_TP(modelo, Tsol, P, props, force_phase=force_phase)
            st["specified"] = "P"
            st[key] = target
            return st

        # Sem saturação possível: estado monofásico
        def f(T):
            st = basic_TP(modelo, T, P, props)
            if st.get("ambiguous_sat", False):
                return np.nan
            return st[key] - target

        Tsol = root_scan(f, Tmin, Tmax, n=320)
        st = basic_TP(modelo, Tsol, P, props)
        st["specified"] = "P"
        st[key] = target
        return st

    # -------------------------------------------------------------------------
    # TH ou TS: temperatura é a PRIMEIRA variável e permanece especificada
    # -------------------------------------------------------------------------
    if par in {"TH", "TS"}:
        T = C_para_K(valor1)
        target = valor2
        key = "H" if par == "TH" else "S"

        if not (Tmin <= T <= Tmax):
            raise ValueError("Temperatura fora dos limites numéricos do modelo.")

        sat = None
        if T < Tc:
            try:
                sat = sat_T(modelo, T, props)
            except Exception:
                sat = None

        if sat is not None:
            low = sat["hL"] if key == "H" else sat["sL"]
            high = sat["hV"] if key == "H" else sat["sV"]

            tol = 1.0e-9 * max(1.0, abs(low), abs(high))
            if low - tol <= target <= high + tol:
                x = (target - low) / (high - low)
                st = two_phase_from_sat(sat, x, specified="T")
                st[key] = target
                return st

            if target < low:
                # líquido comprimido: P > Psat
                Plo = sat["P_Pa"] * (1.0 + 2.0e-7)
                Phi = Pmax
                force_phase = "liquido"
            else:
                # vapor: P < Psat
                Plo = Pmin
                Phi = sat["P_Pa"] * (1.0 - 2.0e-7)
                force_phase = "vapor"

            def f(P):
                st = basic_TP(modelo, T, P, props, force_phase=force_phase)
                return st[key] - target

            Psol = root_scan(f, Plo, Phi, n=320, log=True)
            st = basic_TP(modelo, T, Psol, props, force_phase=force_phase)
            st["specified"] = "T"
            st[key] = target
            return st

        # Monofásico: busca P
        def f(P):
            st = basic_TP(modelo, T, P, props)
            if st.get("ambiguous_sat", False):
                return np.nan
            return st[key] - target

        Psol = root_scan(f, Pmin, Pmax, n=360, log=True)
        st = basic_TP(modelo, T, Psol, props)
        st["specified"] = "T"
        st[key] = target
        return st

    raise RuntimeError("Fluxo de flash não tratado.")


# =============================================================================
# 10. Cp E JOULE-THOMSON
# =============================================================================


def muJT_two_phase(modelo, fluido, state):
    """muJT=(dT/dP)_H [K/bar] por diferença central de flashes PH."""
    P0_bar = state["P_Pa"] / 1.0e5
    H0 = state["H"]

    dP = max(DP_JT_BAR, 1.0e-4 * max(1.0, P0_bar))
    Pminus = max(P_MIN_BAR, P0_bar - dP)
    Pplus = min(P_MAX_BAR, P0_bar + dP)

    if Pplus <= Pminus:
        return np.nan

    try:
        sm = flash_core(modelo, fluido, "PH", Pminus, H0)
        sp = flash_core(modelo, fluido, "PH", Pplus, H0)
        Tm = sm["T_K"]
        Tp = sp["T_K"]
        return (Tp - Tm) / (Pplus - Pminus)
    except Exception:
        return np.nan


def augment_properties(modelo, fluido, state):
    if state.get("ambiguous_sat", False):
        state["Cp"] = np.nan
        state["muJT"] = np.nan
        return state

    phase = state["phase"]

    if phase in {"bifasico", "liquido_saturado", "vapor_saturado"}:
        state["Cp"] = np.inf if phase == "bifasico" else np.nan
        state["muJT"] = muJT_two_phase(modelo, fluido, state)
        return state

    modelo_u, props = get_model_data(modelo, fluido)
    Cp, mu = cp_mu_TP(
        modelo_u,
        state["T_K"],
        state["P_Pa"],
        props,
        phase,
    )
    state["Cp"] = Cp
    state["muJT"] = mu
    return state


# =============================================================================
# 11. FUNÇÃO PÚBLICA PRINCIPAL
# =============================================================================


def calcular_estado(modelo, fluido, par, valor1, valor2):
    """
    Executa o flash e retorna um dicionário com todas as propriedades.
    """
    st = flash_core(modelo, fluido, par, valor1, valor2)
    st = augment_properties(modelo, fluido, st)
    return st


# =============================================================================
# 12. IMPRESSÃO DOS RESULTADOS
# =============================================================================


def imprimir_estado(st, modelo, fluido, par):
    print("\n" + "=" * 74)
    print("FLASH DE PROPRIEDADES - FLUIDO PURO")
    print("=" * 74)
    print(f"Modelo               : {modelo.upper()}")
    print(f"Fluido               : {fluido}")
    print(f"Par de entrada       : {par.upper()}")

    if modelo.upper() == "PR":
        print("Referência PR         : H=S=0 em 25 °C e 1 bar")
    else:
        print("Referência H,S        : nativa da EOS de Leachman")

    if st.get("ambiguous_sat", False):
        sat = st["sat"]
        print("\nESTADO SOBRE A CURVA DE SATURAÇÃO COM TP")
        print("T e P, sozinhos, não determinam a qualidade x.")
        print(f"T                    : {K_para_C(sat['T_K']):.6f} °C")
        print(f"P                    : {sat['P_Pa']/1e5:.6f} bar")
        print(f"H líquido saturado   : {sat['hL']:.6f} kJ/kg")
        print(f"H vapor saturado     : {sat['hV']:.6f} kJ/kg")
        print(f"S líquido saturado   : {sat['sL']:.6f} kJ/(kg K)")
        print(f"S vapor saturado     : {sat['sV']:.6f} kJ/(kg K)")
        print("Cp e muJT             : não definidos sem especificar x/fase")
        print("=" * 74)
        return

    print("-" * 74)
    print(f"Fase                 : {st['phase']}")
    if st.get("x") is not None:
        print(f"Título, x            : {st['x']:.8f}")

    print(f"T                    : {K_para_C(st['T_K']):.6f} °C")
    print(f"P                    : {st['P_Pa']/1e5:.6f} bar")
    print(f"H                    : {st['H']:.6f} kJ/kg")
    print(f"S                    : {st['S']:.6f} kJ/(kg K)")

    if np.isinf(st["Cp"]):
        print("Cp                   : infinito / não finito na região bifásica")
    elif np.isnan(st["Cp"]):
        print("Cp                   : não reportado neste estado de saturação")
    else:
        print(f"Cp                   : {st['Cp']:.6f} kJ/(kg K)")

    if np.isnan(st["muJT"]):
        print("mu_JT                : não disponível")
    else:
        print(f"mu_JT                : {st['muJT']:.8f} K/bar")

    if "rho" in st and _finite(st["rho"]):
        print(f"densidade            : {st['rho']:.6f} kg/m3")

    print("=" * 74)


# =============================================================================
# 13. RELATÓRIO EM ARQUIVO .TXT
# =============================================================================

def salvar_relatorio_txt(
    st,
    modelo,
    fluido,
    par,
    valor1,
    valor2,
    nome_arquivo="relatorio_propriedades.txt",
):
    """
    Grava em arquivo .txt os dados de entrada e todas as propriedades
    termodinâmicas calculadas.

    Unidades:
        T      : °C
        P      : bar
        H      : kJ/kg
        S      : kJ/(kg K)
        Cp     : kJ/(kg K)
        muJT   : K/bar
        rho    : kg/m3
    """

    par = par.upper().replace("-", "").replace("_", "")
    modelo = modelo.upper()

    nomes_entrada = {
        "TP": ("Temperatura", "Pressão"),
        "TRHO": ("Temperatura", "Densidade"),
        "PRHO": ("Pressão", "Densidade"),
        "PH": ("Pressão", "Entalpia específica"),
        "TH": ("Temperatura", "Entalpia específica"),
        "PS": ("Pressão", "Entropia específica"),
        "TS": ("Temperatura", "Entropia específica"),
    }

    unidades_entrada = {
        "TP": ("°C", "bar"),
        "TRHO": ("°C", "kg/m3"),
        "PRHO": ("bar", "kg/m3"),
        "PH": ("bar", "kJ/kg"),
        "TH": ("°C", "kJ/kg"),
        "PS": ("bar", "kJ/(kg K)"),
        "TS": ("°C", "kJ/(kg K)"),
    }

    nome1, nome2 = nomes_entrada[par]
    un1, un2 = unidades_entrada[par]

    linhas = []

    linhas.append("=" * 78)
    linhas.append("RELATÓRIO DE PROPRIEDADES TERMODINÂMICAS - FLUIDO PURO")
    linhas.append("=" * 78)
    linhas.append("")
    linhas.append(f"Modelo termodinâmico : {modelo}")
    linhas.append(f"Fluido               : {fluido}")
    linhas.append(f"Par de entrada       : {par}")

    if modelo == "PR":
        linhas.append("Referência de H e S   : H = S = 0 em 25 °C e 1 bar")
    else:
        linhas.append("Referência de H e S   : nativa da EOS de Leachman")

    linhas.append("")
    linhas.append("VARIÁVEIS ESPECIFICADAS")
    linhas.append("-" * 78)
    linhas.append(f"{nome1:<26s}: {valor1:.10g} {un1}")
    linhas.append(f"{nome2:<26s}: {valor2:.10g} {un2}")
    linhas.append("")

    # -----------------------------------------------------------------
    # Caso especial: TP exatamente sobre a curva de saturação
    # -----------------------------------------------------------------
    if st.get("ambiguous_sat", False):

        sat = st["sat"]

        linhas.append("ESTADO TERMODINÂMICO")
        linhas.append("-" * 78)
        linhas.append("Fase / condição          : saturação subespecificada")
        linhas.append(
            "Observação              : T e P na curva de saturação não "
            "determinam a qualidade x."
        )
        linhas.append("")
        linhas.append("PROPRIEDADES DE SATURAÇÃO")
        linhas.append("-" * 78)
        linhas.append(
            f"{'Temperatura':<26s}: {K_para_C(sat['T_K']):.8f} °C"
        )
        linhas.append(
            f"{'Pressão':<26s}: {sat['P_Pa']/1.0e5:.8f} bar"
        )
        linhas.append(
            f"{'H líquido saturado':<26s}: {sat['hL']:.8f} kJ/kg"
        )
        linhas.append(
            f"{'H vapor saturado':<26s}: {sat['hV']:.8f} kJ/kg"
        )
        linhas.append(
            f"{'S líquido saturado':<26s}: {sat['sL']:.8f} kJ/(kg K)"
        )
        linhas.append(
            f"{'S vapor saturado':<26s}: {sat['sV']:.8f} kJ/(kg K)"
        )
        linhas.append(
            f"{'rho líquido saturado':<26s}: {sat['rhoL']:.8f} kg/m3"
        )
        linhas.append(
            f"{'rho vapor saturado':<26s}: {sat['rhoV']:.8f} kg/m3"
        )
        linhas.append(f"{'Cp':<26s}: não definido sem especificar x/fase")
        linhas.append(f"{'mu_JT':<26s}: não definido sem especificar x/fase")

        linhas.append("")
        linhas.append("=" * 78)
        linhas.append("FIM DO RELATÓRIO")
        linhas.append("=" * 78)

        Path(nome_arquivo).write_text(
            "\n".join(linhas),
            encoding="utf-8",
        )

        print(f"Relatório salvo em: {nome_arquivo}")
        return

    # -----------------------------------------------------------------
    # Estado completamente especificado
    # -----------------------------------------------------------------
    linhas.append("RESULTADOS CALCULADOS")
    linhas.append("-" * 78)

    linhas.append(f"{'Fase':<26s}: {st['phase']}")

    if st.get("x") is not None:
        linhas.append(f"{'Título, x':<26s}: {st['x']:.10f}")

    linhas.append(
        f"{'Temperatura':<26s}: {K_para_C(st['T_K']):.8f} °C"
    )
    linhas.append(
        f"{'Pressão':<26s}: {st['P_Pa']/1.0e5:.8f} bar"
    )
    linhas.append(
        f"{'Entalpia específica, H':<26s}: {st['H']:.8f} kJ/kg"
    )
    linhas.append(
        f"{'Entropia específica, S':<26s}: {st['S']:.8f} kJ/(kg K)"
    )

    # Cp
    if np.isinf(st["Cp"]):
        linhas.append(
            f"{'Cp':<26s}: infinito / não finito na região bifásica"
        )
    elif np.isnan(st["Cp"]):
        linhas.append(
            f"{'Cp':<26s}: não disponível"
        )
    else:
        linhas.append(
            f"{'Cp':<26s}: {st['Cp']:.8f} kJ/(kg K)"
        )

    # Joule-Thomson
    if np.isnan(st["muJT"]):
        linhas.append(
            f"{'Coef. Joule-Thomson':<26s}: não disponível"
        )
    else:
        linhas.append(
            f"{'Coef. Joule-Thomson':<26s}: "
            f"{st['muJT']:.10f} K/bar"
        )

    # Densidade
    if "rho" in st and _finite(st["rho"]):
        linhas.append(
            f"{'Densidade':<26s}: {st['rho']:.8f} kg/m3"
        )

    # Informações adicionais de saturação
    if "sat" in st and isinstance(st["sat"], dict):
        sat = st["sat"]

        linhas.append("")
        linhas.append("PROPRIEDADES DE SATURAÇÃO ASSOCIADAS")
        linhas.append("-" * 78)

        if "T_K" in sat:
            linhas.append(
                f"{'T_sat':<26s}: {K_para_C(sat['T_K']):.8f} °C"
            )

        if "P_Pa" in sat:
            linhas.append(
                f"{'P_sat':<26s}: {sat['P_Pa']/1.0e5:.8f} bar"
            )

        if "hL" in sat:
            linhas.append(
                f"{'H líquido saturado':<26s}: {sat['hL']:.8f} kJ/kg"
            )

        if "hV" in sat:
            linhas.append(
                f"{'H vapor saturado':<26s}: {sat['hV']:.8f} kJ/kg"
            )

        if "sL" in sat:
            linhas.append(
                f"{'S líquido saturado':<26s}: {sat['sL']:.8f} kJ/(kg K)"
            )

        if "sV" in sat:
            linhas.append(
                f"{'S vapor saturado':<26s}: {sat['sV']:.8f} kJ/(kg K)"
            )

        if "rhoL" in sat:
            linhas.append(
                f"{'rho líquido saturado':<26s}: {sat['rhoL']:.8f} kg/m3"
            )

        if "rhoV" in sat:
            linhas.append(
                f"{'rho vapor saturado':<26s}: {sat['rhoV']:.8f} kg/m3"
            )

    linhas.append("")
    linhas.append("OBSERVAÇÕES")
    linhas.append("-" * 78)

    if st["phase"] == "bifasico":
        linhas.append(
            "Estado bifásico: o título x é determinado pela propriedade "
            "secundária especificada (H, S ou densidade)."
        )
        linhas.append(
            "Cp não é tratado como uma propriedade finita de fase homogênea."
        )

    if modelo == "PR":
        linhas.append(
            "Os valores absolutos de H e S dependem do estado de referência "
            "definido no modelo PR."
        )

    if modelo == "LEACHMAN":
        linhas.append(
            "Os valores de H e S seguem a referência termodinâmica da "
            "formulação de Leachman implementada."
        )

    linhas.append("")
    linhas.append("=" * 78)
    linhas.append("FIM DO RELATÓRIO")
    linhas.append("=" * 78)

    Path(nome_arquivo).write_text(
        "\n".join(linhas),
        encoding="utf-8",
    )

    print(f"Relatório salvo em: {nome_arquivo}")




# =============================================================================
# 14. EXECUÇÃO
# =============================================================================

if __name__ == "__main__":

    estado = calcular_estado(
        MODELO,
        FLUIDO,
        PAR_ENTRADA,
        VALOR1,
        VALOR2,
    )

    imprimir_estado(
        estado,
        MODELO,
        FLUIDO,
        PAR_ENTRADA,
    )

    # ----------------------------------------------------------
    # Nome automático do relatório
    # ----------------------------------------------------------

    nome_relatorio = (
        f"report_{MODELO}_{FLUIDO}_{PAR_ENTRADA}.txt"
    )

    salvar_relatorio_txt(
        st=estado,
        modelo=MODELO,
        fluido=FLUIDO,
        par=PAR_ENTRADA,
        valor1=VALOR1,
        valor2=VALOR2,
        nome_arquivo=nome_relatorio,
    )

