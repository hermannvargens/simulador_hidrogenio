# Simulador Interativo de Condicionamento e Expansão de H₂

LINK DE ACESSO [https://pvhwbeder33zgds2rpugc2.streamlit.app/#diagrama-p-h-interativo]

Aplicativo web desenvolvido em Python e Streamlit para modelagem termodinâmica e análise exergética de sistemas de condicionamento de hidrogênio puro. A ferramenta permite a construção dinâmica de fluxogramas de processo e resolve balanços de massa, energia e entropia operando sob diferentes equações de estado.

## Como o Simulador Funciona

O motor de cálculo termodinâmico (`thermo_engine.py`) opera como uma interface unificada que abstrai a complexidade do cálculo de propriedades. Ele recebe parâmetros operacionais (Temperatura, Pressão, Eficiência Isentrópica) e executa *flashes* termodinâmicos em série para cada equipamento definido no fluxograma.

O usuário pode alternar entre duas abordagens termodinâmicas:

* **Gás Ideal:** Utiliza capacidades caloríficas constantes ($c_p = 14,30 \text{ kJ/kg.K}$ e $R = 4,124 \text{ kJ/kg.K}$) para cálculos analíticos de entalpia e entropia.
* **Gás Real (Leachman et al., 2009):** Integra a equação de estado fundamental para o hidrogênio normal, capturando desvios de idealidade como o efeito Joule-Thomson e fornecendo propriedades exatas (entalpia e entropia residuais) de forma numérica.



## Funcionalidades Principais

* **Construtor Dinâmico de Topologia:** Permite adicionar, remover e reordenar equipamentos livremente através de uma barra lateral interativa. O usuário pode simular processos multiestágios inserindo sequências customizadas (ex: Compressor $\to$ Aftercooler $\to$ Compressor $\to$ Aftercooler).
* **Equipamentos Suportados:**
* **Compressor:** Calcula potência consumida, temperatura de descarga real baseada na eficiência isentrópica e geração de entropia.
* **Aftercooler:** Calcula carga térmica de resfriamento e a vazão de utilidade (água) necessária por balanço de energia adiabático.
* **Válvula de Controle:** Simula a expansão isoentálpica, capturando aquecimento ou resfriamento intrínseco (inversão térmica).
* **Expansor (Turbina):** Calcula a potência mecânica recuperada e o estado real de descarga.


* **Modo Comparativo:** Uma função de auditoria que executa a mesma topologia de processo em *background* utilizando simultaneamente os modelos Ideal e Real, gerando um relatório direto de divergências.

## Resultados e Gráficos Obtidos

Ao executar a simulação, o sistema gera instantaneamente os seguintes artefatos de engenharia:

### 1. Process Flow Diagram (PFD) Interativo

Um diagrama de blocos renderizado na tela (`streamlit-flow`). Os nós representam os equipamentos e exibem diretamente os resultados consolidados de Trabalho ($\dot{W}$), Calor ($\dot{Q}$) e Exergia Destruída ($\dot{E}_D$). As arestas (correntes) mostram a Pressão e Temperatura de cada estágio.

### 2. Diagrama P-h (Pressão-Entalpia)

Um gráfico interativo construído com Plotly que mapeia termodinamicamente o ciclo.

* Plota o envelope de saturação do hidrogênio (curvas de líquido e vapor saturado).


* Traça o caminho real de evolução da corrente no processo.
* Sobrepõe as trajetórias isentrópicas teóricas (em linhas tracejadas) para compressores e expansores, permitindo a visualização direta do afastamento da reversibilidade.
* Evidencia graficamente a linha vertical de entalpia constante durante o estrangulamento na válvula e o respectivo efeito sobre a temperatura.

### 3. Tabelas de Balanço Termodinâmico

* **Correntes do Processo:** Matriz contendo $T$, $P$, $h$ e $s$ mapeados da alimentação (S1) até o produto final.
* **Balanço de Exergia:** Relatório por bloco detalhando potência, transferência de calor, taxa de geração de entropia ($\dot{S}_{gen}$) e exergia destruída ($\dot{E}_D$) com base no estado morto de referência ($25^\circ\text{C}$ e $1\text{ bar}$).

### 4. Tabela de Desvios (Modo Comparativo)

Compara os impactos de projeto caso o hidrogênio fosse assumido como gás ideal versus gás real, calculando a diferença percentual de variáveis críticas como:

* Consumo de energia no eixo do compressor (+% em kW).
* Temperatura de descarga e vazão de água de resfriamento.
* Efeito Joule-Thomson explícito na válvula de controle ($\Delta T$ térmico vs. $\Delta T = 0$ do gás ideal).

## Estrutura do Repositório

* `app.py`: Interface principal do usuário (Streamlit) e renderização do PFD e gráficos.
* `thermo_engine.py`: Motor unificado de *flashes* termodinâmicos (Ideal/Real).
* `unit_operations.py`: Módulo contendo a matemática dos equipamentos e balanços de primeira e segunda lei da termodinâmica.
* `Leachman.py`: Core numérico e sub-rotinas da Equação de Estado de Leachman (fornecido/pré-existente).



## Instalação e Execução Local

1. Instale as bibliotecas necessárias presentes no ambiente:

```bash
pip install -r requirements.txt

```

2. Certifique-se de que o arquivo numérico base `Leachman.py` está no diretório raiz do projeto.


3. Inicie o servidor local do Streamlit:

```bash
streamlit run app.py

```
