# Capítulo 7 - Integrais de Sommerfeld e Dipolo Impresso

Este diretório reúne a tradução para C++ dos códigos MATLAB do Capítulo 7 do livro de Davidson. O capítulo estuda o potencial escalar em microfita aterrada por meio de integrais de Sommerfeld e usa esse mesmo potencial na formulação de uma equação integral mista para um dipolo impresso fino.

## Status

Estado atual do capítulo:

- `scalar_pot.m` -> traduzido em `scalar_pot.cpp`;
- `V_pot_eps.m` -> traduzido em `V_pot_eps.cpp`;
- `V_pot_height.m` -> traduzido em `V_pot_height.cpp`;
- `MoM_Som.m` -> traduzido em `MoM_Som.cpp`;
- `D_TE.m`, `D_TM.m`, `F.m`, `F_reg1.m`, `F_reg2.m`, `F_reg3.m`, `F_static.m`, `V_int.m` e `root_D_TM.m` -> traduzidos em arquivos homônimos em `src/`.

O núcleo numérico do capítulo já está funcional. As principais adaptações desta versão em C++ são:

- substituição de `input(...)` por CLI no `MoM_Som`;
- gravação padronizada de CSVs em `Cap_07/out/`;
- geração de gráficos por script Python, em vez de `plot(...)` interativo;
- uso de formas assintóticas estáveis para `coth`, `sech` e `csch` quando a integração empurra `lambda` para valores altos.

Os quatro programas principais do capítulo já foram confrontados diretamente com o MATLAB original.

## Teoria

O capítulo parte da representação espectral do potencial escalar sobre um substrato dielétrico aterrado. A variável espectral `lambda` separa o comportamento em três regiões:

- região 1: `0 <= lambda <= k_0`, tratada com `lambda = k_0 cos(t)`;
- região 2: `k_0 <= lambda <= k_0 sqrt(eps_r')`, onde a singularidade do modo TM é extraída explicitamente;
- região 3: `lambda >= k_0 sqrt(eps_r')`, onde o termo estático é subtraído para acelerar a convergência.

As funções auxiliares do capítulo implementam exatamente essa decomposição:

- `D_TE` e `D_TM` são os denominadores associados às ondas de superfície TE e TM;
- `F` é o integrando original do potencial;
- `F_reg1`, `F_reg2` e `F_reg3` representam as formas especiais usadas em cada região;
- `F_static` é o termo estático extraído na região 3;
- `V_int` recompõe o potencial escalar total;
- `root_D_TM` localiza o polo TM por bisseção no eixo real.

No `MoM_Som`, o potencial escalar é usado dentro da Mixed Potential Integral Equation. A corrente no dipolo impresso é expandida em funções base cossenoidais de domínio inteiro, e a impedância de entrada é calculada ao longo de uma banda de frequências.

## Estrutura

```text
Cap_07/
├── CMakeLists.txt
├── MoM_Som.cpp
├── V_pot_eps.cpp
├── V_pot_height.cpp
├── scalar_pot.cpp
├── include/
│   ├── chapter7_common.hpp
│   └── chapter7_functions.hpp
├── src/
│   ├── D_TE.cpp
│   ├── D_TM.cpp
│   ├── F.cpp
│   ├── F_reg1.cpp
│   ├── F_reg2.cpp
│   ├── F_reg3.cpp
│   ├── F_static.cpp
│   ├── V_int.cpp
│   ├── chapter7_common.cpp
│   └── root_D_TM.cpp
├── scripts/
│   └── plot_cap07.py
├── out/
│   ├── scalar_pot_*.csv
│   ├── v_pot_eps.csv
│   ├── v_pot_height_*.csv
│   ├── mom_som_*.csv
│   └── *.png
└── README.md
```

## Arquivos MATLAB de Referência

Arquivos usados como base local para a tradução:

- `original_matlab/Chapter 7/scalar_pot.m`
- `original_matlab/Chapter 7/V_pot_eps.m`
- `original_matlab/Chapter 7/V_pot_height.m`
- `original_matlab/Chapter 7/MoM_Som.m`
- `original_matlab/Chapter 7/D_TE.m`
- `original_matlab/Chapter 7/D_TM.m`
- `original_matlab/Chapter 7/F.m`
- `original_matlab/Chapter 7/F_reg1.m`
- `original_matlab/Chapter 7/F_reg2.m`
- `original_matlab/Chapter 7/F_reg3.m`
- `original_matlab/Chapter 7/F_static.m`
- `original_matlab/Chapter 7/V_int.m`
- `original_matlab/Chapter 7/root_D_TM.m`

## Compilação

```bash
cd Cap_07
cmake -S . -B build
cmake --build build -j$(nproc)
```

Executáveis gerados:

```bash
./build/scalar_pot
./build/V_pot_eps
./build/V_pot_height
./build/MoM_Som
```

## Execução

Estudo do integrando e das três regiões da integral:

```bash
./build/scalar_pot
```

Potencial escalar para vários valores de permissividade relativa:

```bash
./build/V_pot_eps
```

Potencial escalar para várias alturas normalizadas:

```bash
./build/V_pot_height
```

Dipolo impresso com MoM:

```bash
./build/MoM_Som --max-mode-index 4 --num-int 12
```

Exemplo de smoke test curto:

```bash
./build/MoM_Som --max-mode-index 4 --num-int 12 --freq-start-scale 0.9 --freq-stop-scale 1.0 --freq-step-scale 0.05
```

Ajuda rápida:

```bash
./build/MoM_Som --help
```

## Saídas

Os arquivos são gerados em `Cap_07/out/`:

- `scalar_pot_dtm.csv`, `scalar_pot_dte.csv`, `scalar_pot_integrand.csv` e arquivos auxiliares `scalar_pot_region*.csv` para as figuras da seção de potencial escalar;
- `scalar_pot_summary.csv` com os parâmetros do caso base;
- `v_pot_eps.csv` e `v_pot_eps_summary.csv`;
- `v_pot_height_mag.csv`, `v_pot_height_phase.csv` e `v_pot_height_summary.csv`;
- `mom_som_impedance.csv`, `mom_som_currents.csv`, `mom_som_vpot_table.csv` e `mom_som_summary.csv`.

## Visualização

```bash
cd Cap_07
python3 scripts/plot_cap07.py
```

Figuras geradas:

- `out/scalar_pot_dtm.png`
- `out/scalar_pot_dte.png`
- `out/scalar_pot_integrand.png`
- `out/scalar_pot_integrand_zoom.png`
- `out/scalar_pot_region1.png`
- `out/scalar_pot_region2_lambda.png`
- `out/scalar_pot_region2_t.png`
- `out/scalar_pot_region3.png`
- `out/v_pot_eps.png`
- `out/v_pot_height_mag.png`
- `out/v_pot_height_phase.png`
- `out/mom_som_impedance.png`
- `out/mom_som_currents.png`

## Resultados de Referência e Validação MATLAB

Na validação recente desta tradução:

- `scalar_pot` localizou `lambda_p ≈ 270.076` para o caso base do script;
- `V_pot_eps` e `V_pot_height` rodaram com os defaults completos sem gerar `NaN` nos CSVs;
- `MoM_Som --max-mode-index 4 --num-int 12 --freq-start-scale 0.9 --freq-stop-scale 1.0 --freq-step-scale 0.05` gerou impedâncias aproximadas de `51.73 + j232.93 ohm`, `72.17 + j287.74 ohm` e `101.43 + j350.47 ohm` para `9.0`, `9.5` e `10.0 GHz`.
- `scalar_pot` coincidiu com o MATLAB com erros típicos entre `1e-12` e `1e-8`, com a única ressalva intencional de o primeiro ponto exportado da `region3` começar ligeiramente acima da singularidade;
- `V_pot_eps` coincidiu com o MATLAB com erro máximo de aproximadamente `3.09e-2` em `|V|`;
- `V_pot_height` coincidiu com o MATLAB com erro máximo de aproximadamente `1.05e-2` em `|V|` e `0.313` grau na fase;
- `MoM_Som` coincidiu com o MATLAB com erro máximo de aproximadamente `0.0036` na parte real de `Zin` e `0.0093` na parte imaginária.

## Limitações e Observações

- A própria pasta original do livro já registra limitação para tangente de perdas não nula; essa observação continua valendo aqui.
- No `scalar_pot_region3.csv`, a amostragem começa ligeiramente acima da singularidade TE/TM compartilhada, para evitar `NaN` no primeiro ponto exportado.
- O `MoM_Som` continua sendo a parte mais cara do capítulo, porque o potencial escalar é pré-computado por integração numérica em vários pontos antes da montagem da matriz.

## Ligação com o Repositório

Este é o módulo correspondente ao **Capítulo 7**. Consulte o [README principal](../README.md) e a [matriz de tradução](../TRANSLATION_MATRIX.md) para acompanhar o restante do projeto.
