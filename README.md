# Davidson C++ Electromagnetics

Repositório de tradução dos códigos MATLAB do livro de D. B. Davidson para C++, com três objetivos principais:

- preservar a fidelidade numérica e algorítmica dos scripts originais;
- tornar cada capítulo mais didático, com teoria e instruções de uso em Markdown;
- deixar os resultados reproduzíveis com build, execução e pós-processamento claros.

Referência principal:

**D. B. Davidson, _Computational Electromagnetics for RF and Microwave Engineering_, 2nd Edition, Cambridge University Press.**

Link da editora:
https://www.cambridge.org/br/universitypress/subjects/engineering/rf-and-microwave-engineering/computational-electromagnetics-rf-and-microwave-engineering-2nd-edition?format=HB&isbn=9780521518918

## Estado Atual

| Capítulo | Tema | Originais disponíveis | Status da tradução | Documentação |
| --- | --- | --- | --- | --- |
| `Cap_02` | FDTD 1D | Sim | Traduzido e funcional | Mais completa |
| `Cap_03` | FDTD 2D/3D | Sim | Demos principais e scripts auxiliares principais traduzidos | Precisa validacao fina |
| `Cap_04` | MoM 2D / eletrostática / dipolo fino | Sim | Traduzido e funcional | Boa base, pode aprofundar |
| `Cap_09` | FEM 1D | Sim | Traduzido e funcional | Boa, mas pode padronizar |
| `Cap_10` | FEM 2D / modos / microstrip | Sim | Conjunto principal traduzido e funcional | Precisa aprofundar |

## Estrutura

- [`Cap_02/`](Cap_02/) - capítulo 2 em C++ com scripts Python de visualização.
- [`Cap_03/`](Cap_03/) - capítulo 3 em C++ com demos 2D, PML e 3D.
- [`Cap_04/`](Cap_04/) - capítulo 4 em C++ com espalhamento TM, MoM estático e dipolo fino.
- [`Cap_09/`](Cap_09/) - capítulo 9 em C++ com solver FEM 1D.
- [`Cap_10/`](Cap_10/) - capítulo 10 em C++ com elementos de Whitney e LTQN.
- `original_matlab/` - base local de comparação com os scripts originais do livro.

Observação: `original_matlab/` é usado apenas como referência local para validação e continuidade da tradução.

## Fluxo de Uso

Cada capítulo possui build, executáveis, scripts e saídas próprios.

Exemplo geral:

```bash
cd Cap_03
cmake -S . -B build
cmake --build build -j$(nproc)
./build/fdtd_2d_demo
```

Os detalhes de cada capítulo estão nos respectivos `README.md`.

## Dependências

Pacotes C++/build:

```bash
sudo apt install build-essential gfortran cmake pkg-config libeigen3-dev libfftw3-dev
```

Pacotes Python para visualização:

```bash
sudo apt install python3 python3-pip python3-numpy python3-matplotlib python3-pil ffmpeg
```

## Documentação e Roadmap

- [`Cap_02/README.md`](Cap_02/README.md) - referência mais completa hoje para o padrão didático desejado.
- [`Cap_03/README.md`](Cap_03/README.md) - estado atual do capítulo 3 e limitações conhecidas.
- [`Cap_04/README.md`](Cap_04/README.md) - tradução do método dos momentos do capítulo 4.
- [`Cap_09/README.md`](Cap_09/README.md) - solver FEM 1D.
- [`Cap_10/README.md`](Cap_10/README.md) - estrutura e status do capítulo 10.
- [`TRANSLATION_MATRIX.md`](TRANSLATION_MATRIX.md) - mapa atual entre arquivos MATLAB e equivalentes C++.
- [`ROADMAP.md`](ROADMAP.md) - plano de adequação, padronização e continuidade.

## Disclaimer

Este projeto é não oficial e não possui afiliação com o autor ou com a editora. As traduções são mantidas para fins educacionais, de estudo e de pesquisa.
