import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from thermo_engine import ThermoEngine
import unit_operations as uo
import Leachman

st.set_page_config(page_title="Simulador de Condicionamento de H2", layout="wide")
st.title("⚡ Simulador de Condicionamento e Expansão de $H_2$")

# ==========================================
# 1. BARRA LATERAL: ENTRADAS E FLUXOGRAMA
# ==========================================
st.sidebar.header("1. Modelo e Condições de Entrada")
modelo_opt = st.sidebar.radio("Modelo de Cálculo Padrão", ["REAL (Leachman)", "IDEAL (Gás Ideal)"])
model_key = "REAL" if "REAL" in modelo_opt else "IDEAL"

st.sidebar.subheader("Corrente de Entrada (Estado 1)")
m_dot_kgh = st.sidebar.number_input("Vazão mássica [kg/h]", value=500.0, step=50.0)
m_dot = m_dot_kgh / 3600.0
T1 = st.sidebar.number_input("Temperatura inicial [°C]", value=25.0)
P1 = st.sidebar.number_input("Pressão inicial [bar]", value=20.0)

st.sidebar.header("2. Estrutura do Fluxograma")
equipamentos_possiveis = ["Compressor", "Aftercooler", "Válvula de Controle", "Expansor"]
flowsheet_sequence = st.sidebar.multiselect(
    "Selecione e ordene os equipamentos:",
    options=equipamentos_possiveis,
    default=["Compressor", "Aftercooler", "Válvula de Controle"]
)

params = {}
with st.sidebar.expander("Parâmetros dos Equipamentos", expanded=True):
    if "Compressor" in flowsheet_sequence:
        st.markdown("**Compressor**")
        params["comp_P2"] = st.number_input("P descarga [bar]", value=100.0, key="c_p2")
        params["comp_eta"] = st.number_input("Eficiência Isentrópica Comp.", value=0.80, key="c_eta")
    if "Aftercooler" in flowsheet_sequence:
        st.markdown("**Aftercooler**")
        params["ac_T3"] = st.number_input("T saída [°C]", value=35.0, key="ac_t3")
        params["ac_Twin"] = st.number_input("T água in [°C]", value=20.0, key="ac_twin")
        params["ac_Twout"] = st.number_input("T água out [°C]", value=30.0, key="ac_twout")
    if "Válvula de Controle" in flowsheet_sequence:
        st.markdown("**Válvula de Controle**")
        params["valv_P4"] = st.number_input("P saída (Válvula) [bar]", value=60.0, key="v_p4")
    if "Expansor" in flowsheet_sequence:
        st.markdown("**Expansor**")
        params["exp_P4"] = st.number_input("P descarga (Expansor) [bar]", value=60.0, key="e_p4")
        params["exp_eta"] = st.number_input("Eficiência Isentrópica Exp.", value=0.80, key="e_eta")

# ==========================================
# 2. MOTOR DE SIMULAÇÃO CENTRALIZADO
# ==========================================
def rodar_simulacao(mod_key):
    """Roda todo o fluxograma para um determinado modelo (IDEAL ou REAL) e retorna os dados."""
    thermo = ThermoEngine(model_type=mod_key)
    st_in = thermo.flash_TP(T1, P1)
    st_in["m_dot"] = m_dot
    
    current_stream = st_in
    stream_results = [{"Corrente": "Entrada (1)", "T [°C]": current_stream["T_C"], 
                       "P [bar]": current_stream["P_bar"], "H [kJ/kg]": current_stream["H"], 
                       "S [kJ/(kg K)]": current_stream["S"]}]
    equipment_results = []

    for idx, equip in enumerate(flowsheet_sequence, start=1):
        if equip == "Compressor":
            current_stream, metrics = uo.simular_compressor(current_stream, params["comp_P2"], params["comp_eta"], thermo)
        elif equip == "Aftercooler":
            current_stream, metrics = uo.simular_aftercooler(current_stream, params["ac_T3"], thermo, 
                                                            Tw_in=params["ac_Twin"], Tw_out=params["ac_Twout"])
        elif equip == "Válvula de Controle":
            current_stream, metrics = uo.simular_valvula(current_stream, params["valv_P4"], thermo)
        elif equip == "Expansor":
            current_stream, metrics = uo.simular_expansor(current_stream, params["exp_P4"], params["exp_eta"], thermo)

        metrics["Etapa"] = f"{idx}. {equip}"
        equipment_results.append(metrics)
        stream_results.append({
            "Corrente": f"Saída {equip} ({idx+1})",
            "T [°C]": current_stream["T_C"],
            "P [bar]": current_stream["P_bar"],
            "H [kJ/kg]": current_stream["H"],
            "S [kJ/(kg K)]": current_stream["S"]
        })
    return stream_results, equipment_results

# ==========================================
# 3. LÓGICA DE EXIBIÇÃO E INTERFACE
# ==========================================
col1, col2 = st.columns(2)
run_simple = col1.button(f"Executar Simulação ({model_key})", use_container_width=True)
run_compare = col2.button("🔍 Modo Comparativo (Ideal vs Real)", type="primary", use_container_width=True)

if run_simple or run_compare:
    try:
        if run_compare:
            str_ideal, eq_ideal = rodar_simulacao("IDEAL")
            str_real, eq_real = rodar_simulacao("REAL")
            stream_results, equipment_results = str_real, eq_real 
        else:
            stream_results, equipment_results = rodar_simulacao(model_key)
            
        st.success("Simulação concluída com sucesso!")
        
        # --- BLOCO: TABELA DE COMPARAÇÃO DIRETA ---
        if run_compare:
            st.header("Tabela de Comparação (Gás Ideal vs. Gás Real)")
            
            eq_i_dict = {m["Tipo"]: m for m in eq_ideal}
            eq_r_dict = {m["Tipo"]: m for m in eq_real}
            comp_data = []

            if "Compressor" in eq_r_dict:
                comp_data.extend([
                    {"Grandeza": "T2s Descarga Isentrópica (°C)", "Gás ideal": eq_i_dict["Compressor"]["st_s"]["T_C"], "Leachman": eq_r_dict["Compressor"]["st_s"]["T_C"]},
                    {"Grandeza": "T2 Descarga Real (°C)", "Gás ideal": eq_i_dict["Compressor"]["T_out [°C]"], "Leachman": eq_r_dict["Compressor"]["T_out [°C]"]},
                    {"Grandeza": "Potência do Compressor W_c (kW)", "Gás ideal": eq_i_dict["Compressor"]["W_dot [kW]"], "Leachman": eq_r_dict["Compressor"]["W_dot [kW]"]},
                    {"Grandeza": "Geração Entropia S_gen,c (kW/K)", "Gás ideal": eq_i_dict["Compressor"]["S_gen [kW/K]"], "Leachman": eq_r_dict["Compressor"]["S_gen [kW/K]"]}
                ])
            if "Aftercooler" in eq_r_dict:
                comp_data.extend([
                    {"Grandeza": "Carga Térmica Q_HX (kW)", "Gás ideal": eq_i_dict["Aftercooler"]["Q_dot [kW]"], "Leachman": eq_r_dict["Aftercooler"]["Q_dot [kW]"]},
                    {"Grandeza": "Vazão de Água m_w (kg/s)", "Gás ideal": eq_i_dict["Aftercooler"]["Vazao_Agua [kg/s]"], "Leachman": eq_r_dict["Aftercooler"]["Vazao_Agua [kg/s]"]},
                    {"Grandeza": "Geração Entropia S_gen,HX (kW/K)", "Gás ideal": eq_i_dict["Aftercooler"]["S_gen [kW/K]"], "Leachman": eq_r_dict["Aftercooler"]["S_gen [kW/K]"]}
                ])
            if "Válvula" in eq_r_dict:
                JT_ideal = eq_i_dict["Válvula"]["T_out [°C]"] - eq_i_dict["Válvula"]["T_in [°C]"]
                JT_real = eq_r_dict["Válvula"]["T_out [°C]"] - eq_r_dict["Válvula"]["T_in [°C]"]
                comp_data.extend([
                    {"Grandeza": "T4 Saída Válvula (°C)", "Gás ideal": eq_i_dict["Válvula"]["T_out [°C]"], "Leachman": eq_r_dict["Válvula"]["T_out [°C]"]},
                    {"Grandeza": "ΔT (Efeito Joule-Thomson) (°C)", "Gás ideal": JT_ideal, "Leachman": JT_real},
                    {"Grandeza": "Geração Entropia S_gen,V (kW/K)", "Gás ideal": eq_i_dict["Válvula"]["S_gen [kW/K]"], "Leachman": eq_r_dict["Válvula"]["S_gen [kW/K]"]}
                ])
            if "Expansor" in eq_r_dict:
                comp_data.extend([
                    {"Grandeza": "Potência Recuperada W_t (kW)", "Gás ideal": eq_i_dict["Expansor"]["W_recuperado [kW]"], "Leachman": eq_r_dict["Expansor"]["W_recuperado [kW]"]},
                    {"Grandeza": "Geração Entropia S_gen,t (kW/K)", "Gás ideal": eq_i_dict["Expansor"]["S_gen [kW/K]"], "Leachman": eq_r_dict["Expansor"]["S_gen [kW/K]"]}
                ])
            
            # Totais
            E_d_tot_i = sum(m["E_destruida [kW]"] for m in eq_ideal)
            E_d_tot_r = sum(m["E_destruida [kW]"] for m in eq_real)
            comp_data.append({"Grandeza": "Exergia Destruída Total E_D (kW)", "Gás ideal": E_d_tot_i, "Leachman": E_d_tot_r})

            for row in comp_data:
                r_val = row["Leachman"]
                i_val = row["Gás ideal"]
                row["Diferença (%)"] = (abs(i_val - r_val) / abs(r_val) * 100) if r_val != 0 else 0.0

            df_comp = pd.DataFrame(comp_data)
            st.dataframe(df_comp.style.format({"Gás ideal": "{:.4f}", "Leachman": "{:.4f}", "Diferença (%)": "{:.2f}%"}))
            st.divider()

        # --- BLOCO: RESULTADOS DA CORRENTE ---
        col_res1, col_res2 = st.columns([1, 1])
        with col_res1:
            st.subheader(f"Correntes do Processo ({model_key})")
            st.dataframe(pd.DataFrame(stream_results).style.format({"T [°C]": "{:.2f}", "P [bar]": "{:.2f}", 
                                                                     "H [kJ/kg]": "{:.2f}", "S [kJ/(kg K)]": "{:.4f}"}))
            
            st.subheader("Desempenho dos Equipamentos")
            df_
