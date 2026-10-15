import streamlit as st
import pandas as pd
from thermo_engine import ThermoEngine
import unit_operations as uo

st.set_page_config(page_title="Simulador de Condicionamento de H2", layout="wide")
st.title("⚡ Simulador de Condicionamento e Expansão de $H_2$")

st.sidebar.header("1. Modelo e Condições de Entrada")
modelo_opt = st.sidebar.radio("Modelo Termodinâmico", ["REAL (Leachman)", "IDEAL (Gás Ideal)"])
model_key = "REAL" if "REAL" in modelo_opt else "IDEAL"

st.sidebar.subheader("Corrente de Entrada (Estado 1)")
m_dot_kgh = st.sidebar.number_input("Vazão mássica [kg/h]", value=500.0, step=50.0)
m_dot = m_dot_kgh / 3600.0
T1 = st.sidebar.number_input("Temperatura inicial [°C]", value=25.0)
P1 = st.sidebar.number_input("Pressão inicial [bar]", value=20.0)

try:
    thermo = ThermoEngine(model_type=model_key)
    st1 = thermo.flash_TP(T1, P1)
    st1["m_dot"] = m_dot
except Exception as e:
    st.error(f"Erro ao inicializar o estado: {e}. Verifique o modelo Leachman.")
    st.stop()

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
        params["comp_eta"] = st.number_input("Eficiência Isentrópica", value=0.80, key="c_eta")
    if "Aftercooler" in flowsheet_sequence:
        st.markdown("**Aftercooler**")
        params["ac_T3"] = st.number_input("T saída [°C]", value=35.0, key="ac_t3")
        params["ac_Twin"] = st.number_input("T água in [°C]", value=20.0, key="ac_twin")
        params["ac_Twout"] = st.number_input("T água out [°C]", value=30.0, key="ac_twout")
    if "Válvula de Controle" in flowsheet_sequence:
        st.markdown("**Válvula**")
        params["valv_P4"] = st.number_input("P saída [bar]", value=60.0, key="v_p4")
    if "Expansor" in flowsheet_sequence:
        st.markdown("**Expansor**")
        params["exp_P4"] = st.number_input("P descarga [bar]", value=60.0, key="e_p4")
        params["exp_eta"] = st.number_input("Eficiência Isentrópica", value=0.80, key="e_eta")

if st.button("Executar Simulação", type="primary"):
    current_stream = st1
    stream_results = [{"Corrente": "Entrada (Estado 1)", "T [°C]": current_stream["T_C"], 
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
            "Corrente": f"Saída {equip} (Estado {idx+1})",
            "T [°C]": current_stream["T_C"],
            "P [bar]": current_stream["P_bar"],
            "H [kJ/kg]": current_stream["H"],
            "S [kJ/(kg K)]": current_stream["S"]
        })

    st.subheader("1. Tabela de Correntes")
    st.dataframe(pd.DataFrame(stream_results).style.format({"T [°C]": "{:.2f}", "P [bar]": "{:.2f}", 
                                                             "H [kJ/kg]": "{:.2f}", "S [kJ/(kg K)]": "{:.4f}"}))

    st.subheader("2. Desempenho dos Equipamentos e Análise Exergética")
    df_eq = pd.DataFrame(equipment_results)
    st.dataframe(df_eq.style.format(precision=2))

    if "E_destruida [kW]" in df_eq.columns:
        st.subheader("3. Destruição de Exergia por Bloco ($\dot{E}_D$)")
        st.bar_chart(df_eq.set_index("Etapa")["E_destruida [kW]"])
