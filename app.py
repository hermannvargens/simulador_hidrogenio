import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from streamlit_flow import streamlit_flow
from streamlit_flow.elements import StreamlitFlowNode, StreamlitFlowEdge
from streamlit_flow.state import StreamlitFlowState
from thermo_engine import ThermoEngine
import unit_operations as uo
import Leachman

st.set_page_config(page_title="Simulador H2 - PFD Flow", layout="wide")
st.title("Simulador Interativo: Condicionamento de $H_2$")

# ==========================================
# 0. INICIALIZAÇÃO DE MEMÓRIA E ESTADOS
# ==========================================
if "simulacao_realizada" not in st.session_state:
    st.session_state.simulacao_realizada = False
if "modo_comparativo" not in st.session_state:
    st.session_state.modo_comparativo = False
if "sequence" not in st.session_state:
    # Sequência inicial padrão sugerida no problema
    st.session_state.sequence = ["Compressor", "Aftercooler", "Válvula de Controle"]

# ==========================================
# 1. BARRA LATERAL: ENTRADAS GERAIS
# ==========================================
st.sidebar.header("1. Condições de Alimentação")
modelo_opt = st.sidebar.radio("Equação de Estado", ["REAL (Leachman)", "IDEAL (Gás Ideal)"])
model_key = "REAL" if "REAL" in modelo_opt else "IDEAL"

m_dot_kgh = st.sidebar.number_input("Vazão mássica [kg/h]", value=500.0, step=50.0)
m_dot = m_dot_kgh / 3600.0
T1 = st.sidebar.number_input("Temperatura inicial [°C]", value=25.0)
P1 = st.sidebar.number_input("Pressão inicial [bar]", value=20.0)

# ==========================================
# 2. BARRA LATERAL: CONSTRUTOR DE FLUXOGRAMA
# ==========================================
st.sidebar.header("2. Estrutura do Processo")

equipamentos_possiveis = ["Compressor", "Aftercooler", "Válvula de Controle", "Expansor"]

# Controles para ADICIONAR novo equipamento
add_col1, add_col2 = st.sidebar.columns([3, 2])
novo_eq = add_col1.selectbox("Novo equipamento", equipamentos_possiveis, label_visibility="collapsed")
if add_col2.button("➕ Adicionar"):
    st.session_state.sequence.append(novo_eq)
    st.session_state.simulacao_realizada = False # Reseta a tela se o processo mudou
    st.rerun()

st.sidebar.markdown("---")
st.sidebar.subheader("Especificações dos Equipamentos")

params = []
to_remove = None

# Renderiza as caixas de configuração dinamicamente e permite REMOVER equipamentos
for i, eq in enumerate(st.session_state.sequence):
    with st.sidebar.expander(f"{i+1}. {eq}", expanded=True):
        p = {}
        if eq == "Compressor":
            p["P_out"] = st.number_input("P descarga [bar]", value=100.0, key=f"pout_{i}")
            p["eta"] = st.number_input("Eficiência Isentrópica", value=0.80, key=f"eta_{i}")
        elif eq == "Aftercooler":
            p["T_out"] = st.number_input("T saída do H2 [°C]", value=35.0, key=f"tout_{i}")
            p["Tw_in"] = st.number_input("T água in [°C]", value=20.0, key=f"twin_{i}")
            p["Tw_out"] = st.number_input("T água out [°C]", value=30.0, key=f"twout_{i}")
        elif eq == "Válvula de Controle":
            p["P_out"] = st.number_input("P saída [bar]", value=60.0, key=f"pout_{i}")
        elif eq == "Expansor":
            p["P_out"] = st.number_input("P descarga [bar]", value=60.0, key=f"pout_{i}")
            p["eta"] = st.number_input("Eficiência Isentrópica", value=0.80, key=f"eta_{i}")
        
        params.append(p)
        
        # Botão de remover equipamento específico
        if st.button("❌ Remover etapa", key=f"rem_{i}"):
            to_remove = i

# Remove o item da lista e recarrega a barra lateral corretamente
if to_remove is not None:
    st.session_state.sequence.pop(to_remove)
    st.session_state.simulacao_realizada = False
    st.rerun()

# ==========================================
# 3. MOTOR DE CÁLCULO DINÂMICO
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

    for idx, (equip, p) in enumerate(zip(st.session_state.sequence, params), start=1):
        if equip == "Compressor":
            current_stream, metrics = uo.simular_compressor(current_stream, p["P_out"], p["eta"], thermo)
        elif equip == "Aftercooler":
            current_stream, metrics = uo.simular_aftercooler(current_stream, p["T_out"], thermo, Tw_in=p["Tw_in"], Tw_out=p["Tw_out"])
        elif equip == "Válvula de Controle":
            current_stream, metrics = uo.simular_valvula(current_stream, p["P_out"], thermo)
        elif equip == "Expansor":
            current_stream, metrics = uo.simular_expansor(current_stream, p["P_out"], p["eta"], thermo)

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
# 4. CONTROLES DE EXECUÇÃO
# ==========================================
col_btn1, col_btn2 = st.columns(2)
executar_click = col_btn1.button(f"Executar Simulação ({model_key})", type="primary", width="stretch")
comparar_click = col_btn2.button("🔍 Comparar Ideal vs Real", width="stretch")

if executar_click or comparar_click:
    if len(st.session_state.sequence) == 0:
        st.warning("Adicione pelo menos um equipamento no fluxograma antes de simular.")
    else:
        with st.spinner("Processando balanços termodinâmicos e gerando PFD..."):
            if comparar_click:
                str_ideal, eq_ideal = rodar_simulacao("IDEAL")
                str_real, eq_real = rodar_simulacao("REAL")
                
                st.session_state.eq_ideal = eq_ideal
                st.session_state.eq_real = eq_real
                st.session_state.stream_results = str_real
                st.session_state.equipment_results = eq_real
                st.session_state.modo_comparativo = True
            else:
                stream_results, equipment_results = rodar_simulacao(model_key)
                st.session_state.stream_results = stream_results
                st.session_state.equipment_results = equipment_results
                st.session_state.modo_comparativo = False
            
            # GERA O GRAFO DO PFD INTERATIVO
            nodes = []
            edges = []
            
            # Formatação em múltiplas linhas com HTML (<br>) garantindo a quebra
            feed_text = (
                f"🟡 <b>Alimentação</b><br>"
                f"T = {P1} °C<br>"
                f"P = {T1} bar"
            )
            nodes.append(StreamlitFlowNode(
                id='feed',
                pos=(50, 150),
                data={'content': feed_text},
                node_type='input',
                source_position='right',
                target_position='left'
            ))

            prev_node = 'feed'
            x_pos = 50
            
            for idx, (equip, metric) in enumerate(zip(st.session_state.sequence, st.session_state.equipment_results), start=1):
                node_id = f"block_{idx}"
                x_pos += 300
                st_out = st.session_state.stream_results[idx]
                
                # Monta os dados separados por linha
                node_text = f"⚙️ <b>{idx}. {equip}</b><br>"
                node_text += f"T = {st_out['T [°C]']:.1f} °C<br>"
                node_text += f"P = {st_out['P [bar]']:.1f} bar<br>"
                
                if "W_dot [kW]" in metric and metric["W_dot [kW]"] != 0:
                    node_text += f"Trabalho = {metric['W_dot [kW]']:.1f} kW<br>"
                if "Q_dot [kW]" in metric and metric["Q_dot [kW]"] != 0:
                    node_text += f"Calor = {metric['Q_dot [kW]']:.1f} kW<br>"
                    
                node_text += f"Exergia = {metric['E_destruida [kW]']:.1f} kW"

                nodes.append(StreamlitFlowNode(
                    id=node_id,
                    pos=(x_pos, 150),
                    data={'content': node_text},
                    node_type='default',
                    source_position='right',
                    target_position='left'
                ))

                edges.append(StreamlitFlowEdge(
                    id=f"edge_{idx}",
                    source=prev_node,
                    target=node_id,
                    animated=True,
                    label=f"S{idx}"
                ))
                prev_node = node_id

            x_pos += 300
            st_final = st.session_state.stream_results[-1]
            
            prod_text = (
                f"📦 <b>Produto H2</b><br>"
                f"T = {st_final['T [°C]']:.1f} °C<br>"
                f"P = {st_final['P [bar]']:.1f} bar"
            )
            nodes.append(StreamlitFlowNode(
                id='product',
                pos=(x_pos, 150),
                data={'content': prod_text},
                node_type='output',
                source_position='right',
                target_position='left'
            ))
            
            edges.append(StreamlitFlowEdge(
                id="edge_final",
                source=prev_node,
                target='product',
                animated=True,
                label=f"S{len(st.session_state.stream_results)}"
            ))

            st.session_state.flow_state = StreamlitFlowState(nodes, edges)
            st.session_state.simulacao_realizada = True

# ==========================================
# 5. RENDERIZAÇÃO DA INTERFACE 
# ==========================================
if st.session_state.simulacao_realizada and len(st.session_state.sequence) > 0:
    
    st.subheader("Process Flow Diagram (PFD)")
    streamlit_flow("pfd_canvas", state=st.session_state.flow_state, fit_view=True, height=450)
    st.divider()

    # TABELA DE COMPARAÇÃO (IDEAL VS REAL)
    if st.session_state.modo_comparativo:
        st.subheader("Comparativo de Desempenho Termodinâmico")
        
        # Como agora pode haver equipamentos repetidos, comparamos pelas Etapas iterativas
        eq_i_dict = {m["Etapa"]: m for m in st.session_state.eq_ideal}
        eq_r_dict = {m["Etapa"]: m for m in st.session_state.eq_real}
        comp_data = []

        for etapa_nome, m_real in eq_r_dict.items():
            m_ideal = eq_i_dict[etapa_nome]
            equip_tipo = m_real["Tipo"]
            
            if equip_tipo == "Compressor":
                comp_data.extend([
                    {"Equipamento": etapa_nome, "Métrica": "T2s Descarga Isentrópica (°C)", "Gás Ideal": m_ideal["st_s"]["T_C"], "Gás Real (Leachman)": m_real["st_s"]["T_C"]},
                    {"Equipamento": etapa_nome, "Métrica": "T2 Descarga Real (°C)", "Gás Ideal": m_ideal["T_out [°C]"], "Gás Real (Leachman)": m_real["T_out [°C]"]},
                    {"Equipamento": etapa_nome, "Métrica": "Potência do Compressor W_c (kW)", "Gás Ideal": m_ideal["W_dot [kW]"], "Gás Real (Leachman)": m_real["W_dot [kW]"]}
                ])
            elif equip_tipo == "Aftercooler":
                comp_data.extend([
                    {"Equipamento": etapa_nome, "Métrica": "Carga Térmica Q_HX (kW)", "Gás Ideal": m_ideal["Q_dot [kW]"], "Gás Real (Leachman)": m_real["Q_dot [kW]"]},
                    {"Equipamento": etapa_nome, "Métrica": "Vazão de Água M_w (kg/s)", "Gás Ideal": m_ideal["Vazao_Agua [kg/s]"], "Gás Real (Leachman)": m_real["Vazao_Agua [kg/s]"]}
                ])
            elif equip_tipo == "Válvula":
                JT_ideal = m_ideal["T_out [°C]"] - m_ideal["T_in [°C]"]
                JT_real = m_real["T_out [°C]"] - m_real["T_in [°C]"]
                comp_data.extend([
                    {"Equipamento": etapa_nome, "Métrica": "ΔT (Efeito Joule-Thomson) (°C)", "Gás Ideal": JT_ideal, "Gás Real (Leachman)": JT_real}
                ])
            elif equip_tipo == "Expansor":
                comp_data.extend([
                    {"Equipamento": etapa_nome, "Métrica": "Potência Recuperada W_t (kW)", "Gás Ideal": m_ideal["W_recuperado [kW]"], "Gás Real (Leachman)": m_real["W_recuperado [kW]"]}
                ])
        
        comp_data.append({"Equipamento": "Sistema Total", "Métrica": "Exergia Destruída Total (kW)", "Gás Ideal": sum(m["E_destruida [kW]"] for m in st.session_state.eq_ideal), "Gás Real (Leachman)": sum(m["E_destruida [kW]"] for m in st.session_state.eq_real)})

        for row in comp_data:
            r_val = row["Gás Real (Leachman)"]
            i_val = row["Gás Ideal"]
            row["Desvio (%)"] = (abs(i_val - r_val) / abs(r_val) * 100) if r_val != 0 else 0.0

        st.dataframe(pd.DataFrame(comp_data).style.format({"Gás Ideal": "{:.3f}", "Gás Real (Leachman)": "{:.3f}", "Desvio (%)": "{:.2f}%"}))

    # DIAGRAMA P-h E TABELAS TÉCNICAS
    col_tab, col_graf = st.columns([1, 1.2])
    
    with col_tab:
        st.subheader("Correntes Termodinâmicas")
        st.dataframe(pd.DataFrame(st.session_state.stream_results).style.format({"T [°C]": "{:.2f}", "P [bar]": "{:.2f}", "H [kJ/kg]": "{:.2f}", "S [kJ/(kg K)]": "{:.4f}"}))
        st.subheader("Balanço de Exergia")
        df_eq = pd.DataFrame(st.session_state.equipment_results).drop(columns=["st_s", "T_in [°C]"], errors="ignore")
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

        H_vals = [s['H [kJ/kg]'] for s in st.session_state.stream_results]
        P_vals = [s['P [bar]'] for s in st.session_state.stream_results]
        fig.add_trace(go.Scatter(
            x=H_vals, y=P_vals, mode='lines+markers+text', 
            text=[f"S{i+1}" for i in range(len(st.session_state.stream_results))], 
            textposition="top center", name='Evolução Real', 
            line=dict(color='red', width=3), marker=dict(size=10)
        ))

        for m in st.session_state.equipment_results:
            if "st_s" in m:
                # Divide pela string do ponto para pegar sempre o índice numérico correto, caso ultrapasse 9 equipamentos
                idx = int(m["Etapa"].split(".")[0])
                fig.add_trace(go.Scatter(
                    x=[st.session_state.stream_results[idx-1]["H [kJ/kg]"], m["st_s"]["H"]], 
                    y=[st.session_state.stream_results[idx-1]["P [bar]"], m["st_s"]["P_bar"]], 
                    mode='lines+markers', name=f'Isentrópica ({m["Tipo"]}-{idx})', 
                    line=dict(color='blue', width=2, dash='dot'), marker=dict(size=6)
                ))

        fig.update_layout(xaxis_title="Entalpia (kJ/kg)", yaxis_title="Pressão (bar)", yaxis_type="log", margin=dict(l=20, r=20, t=30, b=20), height=500)
        st.plotly_chart(fig, use_container_width=True)
