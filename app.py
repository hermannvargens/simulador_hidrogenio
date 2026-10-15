# app.py
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
        # Roda o selecionado se for modo simples, ou obrigatoriamente os dois se for comparação
        if run_compare:
            str_ideal, eq_ideal = rodar_simulacao("IDEAL")
            str_real, eq_real = rodar_simulacao("REAL")
            # Para exibição dos gráficos, usaremos o real como base
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

            # Calcula % de diferença baseada no Leachman
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
            df_eq = pd.DataFrame(equipment_results).drop(columns=["st_s", "T_in [°C]"], errors="ignore")
            st.dataframe(df_eq.style.format(precision=2))

        # --- BLOCO: DIAGRAMA P-h (PLOTLY) ---
        with col_res2:
            st.subheader("Diagrama P-h Interativo")
            fig = go.Figure()

            # Plota o domo de saturação do H2 usando as rotinas do Leachman
            try:
                # Intervalo do ponto triplo até quase o ponto crítico
                T_sat_range = np.linspace(Leachman.TTRIPLE_L + 0.1, Leachman.TC_L - 0.2, 80)
                h_L, h_V, p_sat = [], [], []
                for T in T_sat_range:
                    sat_state = Leachman.l_sat_state_T(T)
                    h_L.append(sat_state['hL'])
                    h_V.append(sat_state['hV'])
                    p_sat.append(sat_state['P_Pa'] / 1e5)  # bar

                fig.add_trace(go.Scatter(x=h_L, y=p_sat, mode='lines', line=dict(color='gray', dash='dash'), name='Líquido Saturado'))
                fig.add_trace(go.Scatter(x=h_V, y=p_sat, mode='lines', line=dict(color='gray', dash='dash'), name='Vapor Saturado'))
            except Exception as e:
                st.warning(f"Não foi possível desenhar o envelope de saturação: {e}")

            # Extrai os caminhos para traçar a evolução do processo
            H_vals = [s['H'] for s in stream_results]
            P_vals = [s['P_bar'] for s in stream_results]
            text_labels = [f"Est. {i+1}" for i in range(len(stream_results))]
            
            # Linha principal do processo
            fig.add_trace(go.Scatter(
                x=H_vals, y=P_vals, mode='lines+markers+text',
                text=text_labels, textposition="top center",
                name='Caminho Real', line=dict(color='red', width=3),
                marker=dict(size=10, symbol='circle', color='darkred')
            ))

            # Para processos isentrópicos, adicionar as linhas tracejadas (1->2s ou 3->4s)
            for m in equipment_results:
                if "st_s" in m:
                    idx = int(m["Etapa"][0]) # Pega o índice inicial baseado na string (Ex: "1. Compressor")
                    H_in = stream_results[idx-1]["H"]
                    P_in = stream_results[idx-1]["P_bar"]
                    H_out_s = m["st_s"]["H"]
                    P_out_s = m["st_s"]["P_bar"]
                    
                    fig.add_trace(go.Scatter(
                        x=[H_in, H_out_s], y=[P_in, P_out_s],
                        mode='lines+markers', name=f'Processo Isentrópico ({m["Tipo"]})',
                        line=dict(color='blue', width=2, dash='dot'),
                        marker=dict(size=6, color='blue')
                    ))

            fig.update_layout(
                xaxis_title="Entalpia (kJ/kg)",
                yaxis_title="Pressão (bar)",
                yaxis_type="log", # Eixo y logarítmico padrão para diagramas P-h
                margin=dict(l=20, r=20, t=30, b=20),
                legend=dict(x=0.02, y=0.98),
                height=500
            )
            
            st.plotly_chart(fig, use_container_width=True)
            
            st.caption("No processo de estrangulamento da Válvula (Est. 3 $\\rightarrow$ 4), observe o comportamento vertical da linha (Entalpia constante). O aquecimento residual (Efeito Joule-Thomson) ocorre porque, para fluídos reais nestas condições, a isoterma tem curvatura no diagrama P-h.")

    except Exception as e:
        st.error(f"Erro durante a execução: {str(e)}")
