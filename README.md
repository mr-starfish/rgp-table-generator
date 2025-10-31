# Gerador de Tabelas `.rgp` para ANSYS CFX (CoolProp)

Este projeto fornece um **gerador robusto de tabelas termodinâmicas** no formato `.rgp`, usado pelo **ANSYS CFX**.\
Ele utiliza a biblioteca [CoolProp](http://www.coolprop.org/) para cálculos termodinâmicos e adiciona várias **melhorias de usabilidade e confiabilidade**:

- Barras de progresso (`tqdm`) durante a geração e reescrita.
- Warnings coloridos no terminal.
- Correção automática de “buracos” (NaN) e valores constantes suspeitos.
- Sistema de retomada automática (continua de onde parou sem recalcular tudo).
- Relatório consolidado de warnings ao final.
- Geração de **dois arquivos**:
  - `<NOME>.rgp` → original (cálculo direto)
  - `<NOME>_CORRIGIDO.rgp` → versão pós-processada (valores corrigidos e suavizados)

---

## 📦 Instalação

Clone este repositório e instale as dependências:

```bash
git clone https://github.com/seu-usuario/rgp-table-generator.git
cd rgp-table-generator
python -m venv venv
source venv/bin/activate  # Linux/Mac
venv\\Scripts\\activate     # Windows

pip install -r requirements.txt
```

Dependências principais:

- `coolprop`
- `tqdm`
- `colorama`

---

## 🚀 Como usar

Execute o script de forma interativa:

```bash
python rgp_table_generator.py
```

O programa guiará o usuário passo a passo:

1. Escolha entre **fluido puro** ou **mistura**.
   - Fluido puro → basta informar o nome (ex.: `CO2`).
   - Mistura → informe o **primeiro fluido**, depois o **segundo fluido** e em seguida as **frações molares** (ex.: `0.7` e `0.3`).
2. Escolha se deseja usar **faixas automáticas sugeridas** de **Temperatura (K)** e **Pressão (Pa)** ou inserir manualmente.
3. Defina o **tamanho da malha** (número de pontos em temperatura e pressão). Os valores devem ser múltiplos de 5.
4. Informe o **nome da tabela** (máx. 8 caracteres) e o **nome do arquivo de saída**.

Ao final da execução, serão gerados dois arquivos:

- `CO2.rgp`
- `CO2_CORRIGIDO.rgp`

---

## 📂 Estrutura do arquivo `.rgp`

Cada arquivo contém:

- Cabeçalho com propriedades críticas e limites.
- 9 tabelas de propriedades (densidade, entalpia, energia interna, Cp, Cv, viscosidade, condutividade etc.).
- Blocos de saturação (`Tsat` e `φsat`).

---

## ⚙️ Pós-processamento e correção

Após gerar o `.rgp` base, o script executa o **PASSO 2/2**:

- Varre os coletores (`tables`, `sat_t`, `sat_phi`).
- Aplica correções com:
  - `fix_grid`: normaliza grades 2D (P×T).
  - `fix_vector`: corrige vetores 1D (Tsat e φ).
- Remove **zeros e NaN**, preenchendo lacunas por interpolação e suavização.
- Mantém a monotonicidade sempre que possível.

### Parâmetros internos (personalizáveis)

- `zero_tol=1e-9` → valores abaixo disso são tratados como zero.
- `min_const_run=2` → mínimo de pontos iguais para ser considerado platô.
- `eps_equal=1e-12` → tolerância para igualdade numérica.

Exemplo de uso direto:

```python
fix_grid(valores, n_t=25, zero_tol=1e-8, min_const_run=3)
fix_vector(valores, zero_tol=1e-10)
```

---

## 🔄 Retomada automática

Se o cálculo for interrompido, ao reiniciar:

- O script detecta a última tabela concluída.
- Continua a partir dela.
- Evita recalcular o que já foi feito.

---

## 📊 Exemplo prático

### Fluido puro (CO₂)

```bash
python rgp_table_generator.py
```

Entradas:

```
1   # Fluido puro
CO2 # Nome do fluido
N   # Usar faixas automáticas
25  # Pontos em T
25  # Pontos em P
CO2 # Nome da tabela
CO2 # Nome do arquivo de saída
```

Saídas:

```
CO2.rgp
CO2_CORRIGIDO.rgp
```

### Mistura (Metano + CO₂)

```
2   # Mistura
METHANE  # Primeiro fluido
CO2      # Segundo fluido
0.4      # Fração molar do primeiro fluido
0.6      # Fração molar do segundo fluido
N        # Usar faixas automáticas
25       # Pontos em T
25       # Pontos em P
METHCO2  # Nome da tabela
METHCO2  # Nome do arquivo
```

---

## 🛠️ Personalização avançada

O código expõe funções auxiliares (`fix_grid`, `fix_vector`, coletores e parsers) que permitem ajustes manuais e integração com outros fluxos de trabalho.

---

## 📑 Notas avançadas e exemplos

### WarningCounter

- O programa pode gerar muitos avisos por tabela durante o cálculo.
- Esses **warnings foram silenciados** no loop principal (print comentado).
- Apenas um **resumo final** é exibido com a contagem total.
- É possível reativar os prints detalhados para depuração.

### Limitações práticas

- **Fluido puro**: mesmo em malhas grandes (ex.: 800×1000 pontos), leva cerca de ~10 minutos em um PC comum.
- **Misturas**: tempo cresce muito mais. Exemplo real: mistura 97% CO₂ + 3% CH₄ → ~36 horas.
- **Tabela 6** (derivadas parciais) é sempre a mais lenta.

### Cobertura do CoolProp

- Nem todos os fluidos/misturas são suportados. Apenas os casos simples disponíveis no CoolProp (QPROP).
- Se o cálculo falhar, significa que o CoolProp não possui suporte.
- Próximo da região crítica, o código tenta suavizar e ajustar automaticamente.

### Retomada e reescrita

- Sempre são gerados **dois arquivos**:
  - `<NOME>.rgp` (original, cálculo direto).
  - `<NOME>_CORRIGIDO.rgp` (pós-processado: zeros removidos, buracos interpolados).
- O arquivo original **não é sobrescrito**.
- O **PASSO 2/2** (reescrita) sempre é executado para garantir padronização.

### Compatibilidade ANSYS CFX

- A estrutura do `.rgp` segue fielmente as exigências do CFX.
- Caso contrário, o software rejeita o arquivo.
- Restrições importantes: nome ≤8 caracteres e número de pontos múltiplo de 5.

### Exemplo de saída (`CO2`, 10×15)

Executando:

```bash
python rgp_table_generator.py
```

Com entradas pequenas (10×15), o script gera:

```
CO2_10_15.rgp
CO2_10_15_CORRIGIDO.rgp
```

Trecho inicial do `.rgp`:

```text
$$$$HEADER 
$$$CO2 
    1 
$$PARAM 
    26 
DESCRIPTION 
CO2 [HelmholtzEOSBackend, CoolProp 6.8.0] 
NAME 
CO2 
INDEX 
CO2 
MODEL 
    3 
UNITS 
    1 
PMIN_SUPERHEAT 
    1.0000000E+006 
PMAX_SUPERHEAT 
    1.0100000E+008 
TMIN_SUPERHEAT 
    2.4000000E+002 
TMAX_SUPERHEAT 
    1.0400000E+003 
...
```

Isso dá ao usuário uma visão clara da estrutura antes de rodar malhas maiores.

---

## 📜 Licença

MIT — Uso livre, com aviso de copyright mantido. Sem garantia de funcionamento ou responsabilidade por danos.

