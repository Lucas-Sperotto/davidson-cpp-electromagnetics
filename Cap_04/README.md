# Capítulo 4 - Método dos Momentos

Este diretório reúne a tradução para C++ dos códigos MATLAB do Capítulo 4 do livro de Davidson. O foco do capítulo é o método dos momentos em três contextos complementares:

- espalhamento TM por um cilindro PEC em 2D;
- distribuição de carga em um fio retilíneo no regime estático;
- corrente em um dipolo fino com modelos de excitação por franja magnética e delta gap.

## Status

Estado atual do capítulo:

- `MoM_2D_TM.m` -> traduzido em `MoM_2D_TM.cpp`;
- `MoM_TM_solver.m` -> traduzido em `src/mom_tm_solver.cpp/.hpp`;
- `cyl_TM_echo_width.m` -> traduzido em `src/cyl_tm_echo_width.cpp/.hpp`;
- `static_mom.m` -> traduzido em `static_mom.cpp`;
- `thin_dipole.m` -> traduzido em `thin_dipole.cpp`.

O núcleo numérico do capítulo já está funcional. As principais adaptações desta versão em C++ são:

- substituição da interação por `input(...)` por argumentos de linha de comando quando isso ajuda a manter o fluxo reproduzível;
- gravação padronizada de resultados em `Cap_04/out/`;
- geração de CSVs e figuras em vez de depender apenas dos `plot(...)` interativos do MATLAB.

## Teoria

O capítulo introduz o método dos momentos como uma forma de transformar equações integrais em sistemas lineares.

No caso de `MoM_2D_TM`, o contorno circular PEC é discretizado em segmentos. O campo elétrico incidente em `E_z` é amostrado nos centros desses segmentos, e a corrente superficial desconhecida é expandida em funções base constantes por trecho. O resultado é um sistema denso cuja matriz pode ser montada por:

- avaliação pontual simples dos termos fora da diagonal;
- quadratura de baixa ordem para aproximar a integral do núcleo da EFIE.

Em `static_mom`, o mesmo espírito é usado em regime eletrostático: o fio é particionado em segmentos, a carga linear é tomada constante em cada trecho e resolve-se um sistema para a densidade de carga.

Em `thin_dipole`, o fio fino é modelado com funções base senoidais por partes, o que melhora a representação da corrente ao longo da antena. O capítulo compara dois modelos de alimentação:

- `magnetic frill`;
- `delta gap`.

Além da corrente, o programa calcula a impedância de entrada e o coeficiente de reflexão em um sistema de 75 ohms.

## Estrutura

```text
Cap_04/
├── CMakeLists.txt
├── MoM_2D_TM.cpp
├── static_mom.cpp
├── thin_dipole.cpp
├── include/
│   ├── cyl_tm_echo_width.hpp
│   └── mom_tm_solver.hpp
├── src/
│   ├── cyl_tm_echo_width.cpp
│   └── mom_tm_solver.cpp
├── scripts/
│   └── plot_cap04.py
├── out/
│   ├── mom_tm_current.csv
│   ├── mom_tm_current_summary.csv
│   ├── mom_tm_rcs.csv
│   ├── mom_tm_rcs_summary.csv
│   ├── static_mom_charge.csv
│   ├── static_mom_summary.csv
│   ├── thin_dipole_current.csv
│   └── thin_dipole_summary.csv
└── README.md
```

## Arquivos MATLAB de Referência

Arquivos usados como base local para a tradução:

- `original_matlab/Chapter 4/MoM_2D_TM.m`
- `original_matlab/Chapter 4/MoM_TM_solver.m`
- `original_matlab/Chapter 4/cyl_TM_echo_width.m`
- `original_matlab/Chapter 4/static_mom.m`
- `original_matlab/Chapter 4/thin_dipole.m`

## Compilação

```bash
cd Cap_04
cmake -S . -B build
cmake --build build -j$(nproc)
```

Executáveis gerados:

```bash
./build/MoM_2D_TM
./build/static_mom
./build/thin_dipole
```

## Execução

Espalhamento TM com análise de corrente em uma frequência:

```bash
./build/MoM_2D_TM --analysis-type current
```

Também são aceitos os códigos legados `1` e `2` para refletir o script MATLAB:

```bash
./build/MoM_2D_TM --analysis-type 1
./build/MoM_2D_TM --analysis-type 2
```

Sweep de eco TM em frequência:

```bash
./build/MoM_2D_TM --analysis-type rcs
```

Distribuição de carga eletrostática:

```bash
./build/static_mom
```

Com parâmetros explícitos:

```bash
./build/static_mom --length 1.0 --radius 0.001 --segments 5 --voltage 1.0
```

Dipolo fino:

```bash
./build/thin_dipole --n-seg 80 --length 0.48 --radius 0.005
```

Ajuda rápida:

```bash
./build/MoM_2D_TM --help
./build/static_mom --help
./build/thin_dipole --help
```

## Saídas

Os arquivos são gerados em `Cap_04/out/`:

- `mom_tm_current.csv` - corrente superficial normalizada para o caso de uma frequência;
- `mom_tm_current_summary.csv` - parâmetros da rodada e números de condição;
- `mom_tm_rcs.csv` - radar cross section/echo width no sweep em frequência;
- `mom_tm_rcs_summary.csv` - metadados da rodada RCS;
- `static_mom_charge.csv` - distribuição de carga ao longo do fio;
- `static_mom_summary.csv` - parâmetros da simulação estática;
- `thin_dipole_current.csv` - correntes do dipolo para `delta gap` e `magnetic frill`;
- `thin_dipole_summary.csv` - impedâncias de entrada e coeficientes de reflexão.

## Visualização

```bash
cd Cap_04
python3 scripts/plot_cap04.py
```

Figuras geradas:

- `out/mom_tm_current.png`
- `out/mom_tm_rcs_pi_a.png`
- `out/mom_tm_rcs_lambda.png`
- `out/static_mom_charge.png`
- `out/thin_dipole_current.png`

## Resultados de Referência

Na validação inicial desta tradução:

- `MoM_2D_TM --analysis-type current` gerou números de condição próximos de `7.40476` e `7.72303` para as montagens single-point e quadrature;
- `thin_dipole` com `--n-seg 80 --length 0.48 --radius 0.005` produziu `Zin_mag_frill ≈ 82.5003 + j22.1967` e `Zin_delta_gap ≈ 85.9360 + j16.1577`;
- `static_mom` com os defaults do script original gerou a distribuição de carga correspondente ao caso `L = 1 m`, `a = 0.001 m` e `N = 5`.

## Limitações e Observações

- A interface do MATLAB baseada em prompts foi convertida para CLI para facilitar automação e repetibilidade.
- Os gráficos são gerados a partir dos CSVs por script Python, em vez de `plot(...)` direto no solver.
- O objetivo aqui continua sendo fidelidade algorítmica; por isso, a formulação permanece densa e explícita, sem tentar “otimizar” o método original além do necessário para a tradução.

## Ligação com o Repositório

Este é o módulo correspondente ao **Capítulo 4**. Consulte o [README principal](../README.md) e a [matriz de tradução](../TRANSLATION_MATRIX.md) para acompanhar o restante do projeto.
