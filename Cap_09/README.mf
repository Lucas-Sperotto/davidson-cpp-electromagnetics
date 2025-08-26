# FEM 1D - Linha de Transmissão

Este projeto implementa a **análise 1D por Elementos Finitos (FEM)** de uma linha de transmissão, baseada no livro *Davidson - Computational Electromagnetics for RF and Microwave Engineering*.  
O código calcula o perfil de tensão ao longo da linha e compara com a solução analítica, permitindo avaliar a convergência do método.

## 📂 Estrutura do projeto
- `CMakeLists.txt` → configuração de build com CMake.
- `FEM_1D_solver.cpp/.h` → solver FEM 1D (montagem e solução do sistema).
- `FEM_pp.cpp/.h` → pós-processamento (interpolação FEM linear).
- `FEM_1D_1st.cpp` → programa principal que executa a simulação.
- `plot_fem.py` → script Python para gerar gráficos de convergência e perfis.

## ⚙️ Compilação
É necessário ter **Eigen3** e **CMake** instalados.

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
````

O executável gerado será `fem_tl`.

## ▶️ Execução

O programa pode ser executado com parâmetros opcionais:

```bash
./fem_tl [N_elem_init] [num_meshes] [f] [L] [C] [V_in]
```

* `N_elem_init` → número inicial de elementos (default = 2)
* `num_meshes` → número de refinamentos sucessivos (default = 1)
* `f` → frequência \[Hz] (default = 1.0)
* `L` → indutância por unidade de comprimento (default = 1.0)
* `C` → capacitância por unidade de comprimento (default = 1.0)
* `V_in` → tensão aplicada no final da linha (default = 1.0)

Exemplo:

```bash
./fem_tl 4 3 1.0 1.0 1.0 1.0
```

## 📊 Saídas

Após a execução, arquivos `.csv` serão gerados em `build/`:

* `fem_profile_stage_#.csv` → perfil FEM vs analítico por estágio.
* `fem_nodes_stage_#.csv` → valores nodais.
* `fem_convergence.csv` → curva de convergência (erro RMS).

## 📈 Visualização

Use o script Python para gerar os gráficos:

```bash
cd build
python3 ../plot_fem.py
```

Arquivos gerados:

* `convergence.png` → curva log-log de convergência.
* `fem_profile_stage_#.png` → comparação FEM vs analítico.

## 📖 Referência

* Davidson, D. B. *Computational Electromagnetics for RF and Microwave Engineering*. Cambridge University Press, 2nd ed., 2010.