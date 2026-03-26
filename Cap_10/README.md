# Capítulo 10 - FEM 2D, Modos TE e Microstrip

Este diretório reúne a tradução para C++ dos códigos MATLAB do Capítulo 10 do livro de Davidson ligados a:

- solução eletrostática 2D para uma microstrip encapsulada;
- autovalores e autovetores TE com elementos de Whitney;
- autovalores TE com elementos LTQN.

## Status

Estado atual do capítulo:

- todos os arquivos principais do MATLAB local já possuem equivalente em C++;
- `Static2D`, `Eigen2D` e `Eigen2D_LTQN` compilam e rodam;
- os defaults do `Static2D` foram realinhados com o MATLAB original;
- o capítulo agora gera saídas e resumos numéricos mais explícitos em `out/`;
- a documentação ainda pode crescer em teoria e validação formal C++ versus MATLAB.

## Arquivos MATLAB de Referência

Arquivos usados como base local:

- `Static2D.m`
- `Eigen2D.m`
- `Eigen2D_LTQN.m`
- `sandt.m`
- `sandt_LTQN.m`
- `whitney.m`
- `LTQN.m`
- `curl_LTQN.m`
- `plot_field.m`
- e os auxiliares de malha, topologia, quadratura e contorno do capítulo.

## Estrutura

Principais executáveis:

- `Static2D.cpp` - modo quasi-TEM em microstrip encapsulada;
- `Eigen2D.cpp` - modos TE em guia retangular com elementos de Whitney;
- `Eigen2D_LTQN.cpp` - modos TE com elementos LTQN.

Principais auxiliares:

- malha: `trimesh.cpp`, `edgemake.cpp`;
- geometria/topologia: `simplex2D.cpp`, `find_local_dofs.cpp`, `renumber_dof.cpp`, `renumber_dof_LTQN.cpp`;
- elementos: `whitney.cpp`, `LTQN.cpp`, `curl_LTQN.cpp`;
- matrizes locais: `sandt.cpp`, `s_nodal.cpp`, `sandt_LTQN.cpp`;
- contorno: `free_nodes.cpp`, `free_nodes_mstrip.cpp`, `free_dof.cpp`, `prescr_nodes_mstrip.cpp`;
- pós-processamento: `plot_field.cpp`, `TEeig_err.cpp`.

## Compilação

```bash
cd Cap_10
cmake -S . -B build
cmake --build build -j$(nproc)
```

## Execução

### Static2D

Executa o problema eletrostático da microstrip. Os defaults agora seguem o script MATLAB original:

```bash
./build/Static2D
```

Exemplo com parâmetros explícitos:

```bash
./build/Static2D --a 2.5 --b 1.25 --h 0.5 --w 0.5 --eps-r-sub 2.55 --x-mesh 20 --y-mesh 20
```

### Eigen2D

Executa o problema modal TE com elementos de Whitney:

```bash
./build/Eigen2D
```

Exemplo:

```bash
./build/Eigen2D --x-mesh 8 --y-mesh 4 --num-physical-modes 6 --num-spurious-modes 6
```

### Eigen2D_LTQN

Executa o problema modal TE com elementos LTQN:

```bash
./build/Eigen2D_LTQN
```

Exemplo:

```bash
./build/Eigen2D_LTQN --x-mesh 4 --y-mesh 2
```

Use `--help` em qualquer um dos executáveis para ver as opções disponíveis.

## Saídas

Arquivos gerados em `Cap_10/out/`:

### Static2D

- `field_map.csv` - mapa de potencial e campo elétrico;
- `strip_nodes.csv` - nós prescritos do condutor central;
- `geom_params.csv` - geometria, malha, material e capacitância calculada;
- `static2d_summary.csv` - resumo da rodada.

### Eigen2D

- `field_kc_*.csv` - campos dos modos físicos exportados;
- `spur_field_kc_*.csv` - campos dos modos espúrios exportados;
- `cpp_S_mat.csv`, `cpp_T_mat.csv` - matrizes globais montadas;
- `cpp_eigvals.csv`, `cpp_eigvecs.csv` - autovalores e autovetores do problema generalizado;
- `eigdata_whitney_modes.csv` - tabela consolidada de modos físicos e espúrios;
- `eigen2d_summary.csv` - resumo da rodada.

### Eigen2D_LTQN

- `eigdata_LTQN_modes.csv` - tabela consolidada dos modos LTQN;
- `eigdata_LTQN_relerr.csv` - erro relativo dos primeiros modos úteis;
- `eigen2d_ltqn_summary.csv` - resumo da rodada.

### Figuras

- `microstrip_field.png` - figura do problema estático;
- `img/*.png` - campos vetoriais gerados a partir dos CSVs de modos.

## Scripts de Visualização

Todos os scripts operam diretamente sobre `Cap_10/out/`.

Plot da microstrip:

```bash
cd Cap_10
python3 scripts/plot_field_microstrip.py
```

Plot de um campo específico:

```bash
python3 scripts/plot_field.py field_kc_137.224379.csv --save out/img/field_kc_137.224379.png
```

Plot de todos os campos exportados:

```bash
python3 scripts/plot_all_fields.py
```

## Observações Técnicas

- O capítulo ainda usa estruturas globais para malha e conectividade, o que manteve a tradução próxima do estilo MATLAB.
- O `Static2D` voltou a usar, por padrão, `eps_r_sub = 2.55` e malha `20 x 20`, em linha com o script original.
- O solver LTQN segue funcional e produz valores físicos coerentes, mas ainda merece uma validação documentada e reprodutível contra o MATLAB.
- O `Eigen2D` agora exporta explicitamente modos físicos e espúrios, aproximando melhor o comportamento do script original.

## Próximos Passos Naturais

- documentar a teoria do capítulo em maior profundidade;
- registrar uma validação numérica C++ versus MATLAB com os mesmos parâmetros;
- reduzir gradualmente a dependência de estado global sem perder a legibilidade da tradução.
