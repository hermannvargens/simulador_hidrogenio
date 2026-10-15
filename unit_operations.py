import numpy as np

def simular_compressor(stream_in, P_out, eta_s, thermo):
    st_s = thermo.flash_PS(P_out, stream_in["S"])
    dh_s = st_s["H"] - stream_in["H"]
    dh_real = dh_s / eta_s
    h_out = stream_in["H"] + dh_real
    st_out = thermo.flash_PH(P_out, h_out)
    st_out["m_dot"] = stream_in["m_dot"]
    W_dot = stream_in["m_dot"] * dh_real
    S_gen = stream_in["m_dot"] * (st_out["S"] - stream_in["S"])
    E_d = thermo.T0 * S_gen
    metrics = {
        "Tipo": "Compressor",
        "W_dot [kW]": W_dot,
        "Q_dot [kW]": 0.0,
        "S_gen [kW/K]": S_gen,
        "E_destruida [kW]": E_d,
        "T_out [°C]": st_out["T_C"]
    }
    return st_out, metrics

def simular_aftercooler(stream_in, T_out_target, thermo, cp_w=4.18, Tw_in=20.0, Tw_out=30.0):
    P_out = stream_in["P_bar"]
    st_out = thermo.flash_TP(T_out_target, P_out)
    st_out["m_dot"] = stream_in["m_dot"]
    Q_H2 = stream_in["m_dot"] * (st_out["H"] - stream_in["H"])
    m_water = abs(Q_H2) / (cp_w * (Tw_out - Tw_in))
    Tw_in_K, Tw_out_K = Tw_in + 273.15, Tw_out + 273.15
    ds_w = cp_w * np.log(Tw_out_K / Tw_in_K)
    S_gen = (stream_in["m_dot"] * (st_out["S"] - stream_in["S"])) + (m_water * ds_w)
    E_d = thermo.T0 * S_gen
    metrics = {
        "Tipo": "Aftercooler",
        "W_dot [kW]": 0.0,
        "Q_dot [kW]": Q_H2,
        "Vazao_Agua [kg/s]": m_water,
        "S_gen [kW/K]": S_gen,
        "E_destruida [kW]": E_d,
        "T_out [°C]": st_out["T_C"]
    }
    return st_out, metrics

def simular_valvula(stream_in, P_out, thermo):
    st_out = thermo.flash_PH(P_out, stream_in["H"])
    st_out["m_dot"] = stream_in["m_dot"]
    S_gen = stream_in["m_dot"] * (st_out["S"] - stream_in["S"])
    E_d = thermo.T0 * S_gen
    metrics = {
        "Tipo": "Válvula",
        "W_dot [kW]": 0.0,
        "Q_dot [kW]": 0.0,
        "S_gen [kW/K]": S_gen,
        "E_destruida [kW]": E_d,
        "T_out [°C]": st_out["T_C"]
    }
    return st_out, metrics

def simular_expansor(stream_in, P_out, eta_s, thermo):
    st_s = thermo.flash_PS(P_out, stream_in["S"])
    dh_s = stream_in["H"] - st_s["H"]
    dh_real = eta_s * dh_s
    h_out = stream_in["H"] - dh_real
    st_out = thermo.flash_PH(P_out, h_out)
    st_out["m_dot"] = stream_in["m_dot"]
    W_dot = stream_in["m_dot"] * dh_real
    S_gen = stream_in["m_dot"] * (st_out["S"] - stream_in["S"])
    E_d = thermo.T0 * S_gen
    metrics = {
        "Tipo": "Expansor",
        "W_dot [kW]": -W_dot,
        "W_recuperado [kW]": W_dot,
        "Q_dot [kW]": 0.0,
        "S_gen [kW/K]": S_gen,
        "E_destruida [kW]": E_d,
        "T_out [°C]": st_out["T_C"]
    }
    return st_out, metrics
