import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from streamlit_flow import streamlit_flow
from streamlit_flow.elements import StreamlitFlowNode, StreamlitFlowEdge
from thermo_engine import ThermoEngine
import unit_operations as uo
import Leachman

st.set_page_config(page_title="Simulador H2 - PFD Flow", layout="wide")
st.title("⚡ Simulador Interativo: Condicionamento de $H_2$")

# ==========================================
# 1. BARRA LATERAL: ENTRADAS E FLUXOGRAMA
# ==========================================
st.sidebar.header("1. Condições de Alimentação")
modelo_opt = st.sidebar.radio("Equação de Estado", ["REAL (Leachman)", "IDEAL (Gás Ideal)"])
model_key = "REAL" if "REAL" in modelo_opt else "IDEAL"

m_dot_kgh = st.sidebar.number_input("Vazão mássica [kg/h]", value=500.0, step=50.0)
m_dot = m_dot_kgh / 3600.0
T1 = st.sidebar.number_input("Temperatura inicial [°C]", value=25.0)
P1 = st.sidebar.number_input("Pressão inicial [bar]", value=20.0)

st.sidebar.header("2. Fluxograma de Processo")
equipamentos_possiveis = ["Compressor", "Aftercooler", "Válvula de Controle", "Expansor"]
flowsheet_sequence = st.sidebar.multiselect(
    "Construa sua topologia (ordem de fluxo):",
    options=equipamentos_possiveis,
    default=["Compressor", "Aftercooler", "Válvula de Controle"]
)

params = {}
with st.sidebar.expander("Especificações dos Equipamentos", expanded=True):
    if "Compressor" in flowsheet_sequence:
        st.markdown("**Compressor (K-101)**")
        params["comp_P2"] = st.number_input("P descarga [bar]", value=100.0, key="c_p2")
        params["comp_eta"] = st.number_input("Eficiência Isentrópica", value=0.80, key="c_eta")
    if "Aftercooler" in flowsheet_sequence:
        st.markdown("**Aftercooler (E-101)**")
        params["ac_T3"] = st.number_input("T saída [°C]", value=35.0, key="ac_t3")
        params["ac_Twin"] = st.number_input("T água in [°C]", value=20.0, key="ac_twin")
        params["ac_Twout"] = st.number_input("T água out [°C]", value=30.0, key="ac_twout")
    if "Válvula de Controle" in flowsheet_sequence:
        st.markdown("**Válvula JT (V-101)**")
        params["valv_P4"] = st.number_input("P saída [bar]", value=60.0, key="v_p4")
    if "Expansor" in flowsheet_sequence:
        st.markdown("**Turbina/Expansor (EX-101)**")
        params["exp_P4"] = st.number_input("P descarga [bar]", value=60.0, key="e_p4")
        params["exp_eta"] = st.number_input("Eficiência Isentrópica", value=0.80, key="e_eta")

# ==========================================
# 2. MOTOR DE CÁLCULO
# ==========================================
def rodar_simulacao(mod_key):
    thermo = ThermoEngine(model_type=mod_key)
    st_in = thermo.flash_TP(T1, P1)
    st_in["m_dot"] = m_dot
    
    current_stream = st_in
    stream_results = [{"Corrente": "Alimentação (S1)", "T [°C]": current_stream["T_C"], 
                       "P [bar]": current_stream["P_bar"], "H [kJ/kg]": current_stream["H"], 
                       "S [kJ/(kg K)]": current_stream["S"]}]
    equipment_results = []

    for idx, equip in enumerate(flowsheet_sequence, start=1):
        if equip == "Compressor":
            current_stream, metrics = uo.simular_compressor(current_stream, params["comp_P2"], params["comp_eta"], thermo)
        elif equip == "Aftercooler":
            current_stream, metrics = uo.simular_aftercooler(current_stream, params["ac_T3"], thermo, Tw_in=params["ac_Twin"], Tw_out=params["ac_Twout"])
        elif equip == "Válvula de Controle":
            current_stream, metrics = uo.simular_valvula(current_stream, params["valv_P4"], thermo)
        elif equip == "Expansor":
            current_stream, metrics = uo.simular_expansor(current_stream, params["exp_P4"], params["exp_eta"], thermo)

        metrics["Etapa"] = f"{idx}. {equip}"
        equipment_results.append(metrics)
        stream_results.append({
            "Corrente": f"Saída {equip} (S{idx+1})",
            "T [°C]": current_stream["T_C"],
            "P [bar]": current_stream["P_bar"],
            "H [kJ/kg]": current_stream["H"],
            "S [kJ/(kg K)]": current_stream["S"]
        })
    return stream_results, equipment_results

# ==========================================
# 3. INTERFACE E RENDERIZAÇÃO
# ==========================================
col_btn1, col_btn2 = st.columns(2)
executar = col_btn1.button(f"Executar Simulação ({model_key})", use_container_width=True, type="primary")
comparar = col_btn2.button("🔍 Comparar Ideal vs Real", use_container_width=True)

if executar or comparar:
    with st.spinner("Processando balanços termodinâmicos..."):
        if comparar:
            str_ideal, eq_ideal = rodar_simulacao("IDEAL")
            str_real, eq_real = rodar_simulacao("REAL")
            stream_results, equipment_results = str_real, eq_real
        else:
            stream_results, equipment_results = rodar_simulacao(model_key)
    
    # ------------------------------------------
    # RENDERIZAÇÃO DO PFD INTERATIVO
    # ------------------------------------------
    st.subheader("Process Flow Diagram (PFD)")
    nodes = []
    edges = []
    
    # Nó de Alimentação
    nodes.append(StreamlitFlowNode(
        id='feed',
        pos=(50, 150),
        data={'label': f"🟡 Alimentação\n{P1} bar | {T1}°C"},
        node_type='input',
        source_position='right',
        target_position='left'
    ))

    prev_node = 'feed'
    x_pos = 50
    
    # Construção iterativa de Nós e Arestas
    for idx, (equip, metric) in enumerate(zip(flowsheet_sequence, equipment_results)):
        node_id = f"block_{idx}"
        x_pos += 300
        
        # Formata o texto do nó
        node_text = f"⚙️ {equip}\n"
        if "W_dot [kW]" in metric and metric["W_dot [kW]"] != 0:
            node_text += f"W = {metric['W_dot [kW]']:.1f} kW\n"
        if "Q_dot [kW]" in metric and metric["Q_dot [kW]"] != 0:
            node_text += f"Q = {metric['Q_dot [kW]']:.1f} kW\n"
        node_text += f"Exergia Dest. = {metric['E_destruida [kW]']:.1f} kW"

        nodes.append(StreamlitFlowNode(
            id=node_id,
            pos=(x_pos, 150),
            data={'label': node_text},
            node_type='default',
            source_position='right',
            target_position='left'
        ))

        # Aresta (Corrente) conectando o bloco anterior a este
        st_out = stream_results[idx]
        edges.append(StreamlitFlowEdge(
            id=f"edge_{idx}",
            source=prev_node,
            target=node_id,
            animated=True,
            label=f"S{idx+1} ({st_out['P [bar]']:.0f} bar, {st_out['T [°C]']:.1f} °C)"
        ))
        prev_node = node_id

    # Nó final (Produto)
    x_pos += 300
    st_final = stream_results[-1]
    nodes.append(StreamlitFlowNode(
        id='product',
        pos=(x_pos, 150),
        data={'label': f"📦 Produto H2\n{st_final['P [bar]']:.0f} bar | {st_final['T [°C]']:.1f}°C"},
        node_type='output',
        source_position='right',
        target_position='left'
    ))
    edges.append(StreamlitFlowEdge(
        id="edge_final",
        source=prev_node,
        target='product',
        animated=True,
        label=f"S{len(stream_results)}"
    ))

    # Desenha o PFD na tela
    streamlit_flow("pfd_canvas", nodes=nodes, edges=edges, fit_view=True, height=400)
    
    st.divider()

    # ------------------------------------------
    # TABELA DE COMPARAÇÃO (IDEAL VS REAL)
    # ------------------------------------------
    if comparar:
        st.subheader("Comparativo de Desempenho Termodinâmico")
        eq_i_dict = {m["Tipo"]: m for m in eq_ideal}
        eq_r_dict = {m["Tipo"]: m for m in eq_real}
        comp_data = []

        if "Compressor" in eq_r_dict:
            comp_data.extend([
                {"Métrica": "T2s Descarga Isentrópica (°C)", "Gás Ideal": eq_i_dict["Compressor"]["st_s"]["T_C"], "Gás Real (Leachman)": eq_r_dict["Compressor"]["st_s"]["T_C"]},
                {"Métrica": "T2 Descarga Real (°C)", "Gás Ideal": eq_i_dict["Compressor"]["T_out [°C]"], "Gás Real (Leachman)": eq_r_dict["Compressor"]["T_out [°C]"]},
                {"Métrica": "Potência do Compressor W_c (kW)", "Gás Ideal": eq_i_dict["Compressor"]["W_dot [kW]"], "Gás Real (Leachman)": eq_r_dict["Compressor"]["W_dot [kW]"]}
            ])
        if "Aftercooler" in eq_r_dict:
            comp_data.extend([
                {"Métrica": "Carga Térmica Q_HX (kW)", "Gás Ideal": eq_i_dict["Aftercooler"]["Q_dot [kW]"], "Gás Real (Leachman)": eq_r_dict["Aftercooler"]["Q_dot [kW]"]},
                {"Métrica": "Vazão de Água M_w (kg/s)", "Gás Ideal": eq_i_dict["Aftercooler"]["Vazao_Agua [kg/s]"], "Gás Real (Leachman)": eq_r_dict["Aftercooler"]["Vazao_Agua [kg/s]"]}
            ])
        if "Válvula" in eq_r_dict:
            JT_ideal = eq_i_dict["Válvula"]["T_out [°C]"] - eq_i_dict["Válvula"]["T_in [°C]"]
            JT_real = eq_r_dict["Válvula"]["T_out [°C]"] - eq_r_dict["Válvula"]["T_in [°C]"]
            comp_data.extend([
                {"Métrica": "ΔT de Inversão (Joule-Thomson) (°C)", "Gás Ideal": JT_ideal, "Gás Real (Leachman)": JT_real}
            ])
        if "Expansor" in eq_r_dict:
            comp_data.extend([
                {"Métrica": "Potência Recuperada W_t (kW)", "Gás Ideal": eq_i_dict["Expansor"]["W_recuperado [kW]"], "Gás Real (Leachman)": eq_r_dict["Expansor"]["W_recuperado [kW]"]}
            ])
        
        comp_data.append({"Métrica": "Exergia Destruída Total (kW)", "Gás Ideal": sum(m["E_destruida [kW]"] for m in eq_ideal), "Gás Real (Leachman)": sum(m["E_destruida [kW]"] for m in eq_real)})

        for row in comp_data:
            r_val = row["Gás Real (Leachman)"]
            i_val = row["Gás Ideal"]
            row["Desvio (%)"] = (abs(i_val - r_val) / abs(r_val) * 100) if r_val != 0 else 0.0

        st.dataframe(pd.DataFrame(comp_data).style.format({"Gás Ideal": "{:.3f}", "Gás Real (Leachman)": "{:.3f}", "Desvio (%)": "{:.2f}%"}))

    # ------------------------------------------
    # DIAGRAMA P-h E TABELAS
    # ------------------------------------------
    col_tab, col_graf = st.columns([1, 1.2])
    
    with col_tab:
        st.subheader("Correntes Termodinâmicas")
        st.dataframe(pd.DataFrame(stream_results).style.format({"T [°C]": "{:.2f}", "P [bar]": "{:.2f}", "H [kJ/kg]": "{:.2f}", "S [kJ/(kg K)]": "{:.4f}"}))
        st.subheader("Balanço de Exergia")
        df_eq = pd.DataFrame(equipment_results).drop(columns=["st_s", "T_in [°C]"], errors="ignore")
        st.dataframe(df_eq.style.format(precision=2))

    with col_graf:
        st.subheader("Diagrama P-h do Ciclo")
        fig = go.Figure()

        try:
            T_sat = np.linspace(Leachman.TTRIPLE_L + 0.1, Leachman.TC_L - 0.2, 80)
            h_L, h_V, p_sat = zip(*[(Leachman.l_sat_state_T(t)['hL'], Leachman.l_sat_state_T(t)['hV'], Leachman.l_sat_state_T(t)['P_Pa']/1e5) for t in T_sat])
            fig.add_trace(go.Scatter(x=h_L, y=p_sat, mode='lines', line=dict(color='gray', dash='dash'), name='Líq. Saturado'))
            fig.add_trace(go.Scatter(x=h_V, y=p_sat, mode='lines', line=dict(color='gray', dash='dash'), name='Vap. Saturado'))
        except:
            pass

        H_vals = [s['H [kJ/kg]'] for s in stream_results]
        P_vals = [s['P [bar]'] for s in stream_results]
        fig.add_trace(go.Scatter(x=H_vals, y=P_vals, mode='lines+markers+text', text=[f"S{i+1}" for i in range(len(stream_results))], textposition="top center", name='Evolução Real', line=dict(color='red', width=3), marker=dict(size=10)))

        for m in equipment_results:
            if "st_s" in m:
                idx = int(m["Etapa"][0])
                fig.add_trace(go.Scatter(x=[stream_results[idx-1]["H [kJ/kg]"], m["st_s"]["H"]], y=[stream_results[idx-1]["P [bar]"], m["st_s"]["P_bar"]], mode='lines+markers', name=f'Isentrópica ({m["Tipo"]})', line=dict(color='blue', width=2, dash='dot'), marker=dict(size=6)))

        fig.update_layout(xaxis_title="Entalpia (kJ/kg)", yaxis_title="Pressão (bar)", yaxis_type="log", margin=dict(l=20, r=20, t=30, b=20), height=500)
        st.plotly_chart(fig, use_container_width=True)
