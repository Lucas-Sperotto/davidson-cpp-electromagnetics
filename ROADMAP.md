# Roadmap

Este arquivo consolida o plano atual para transformar o repositório em uma base:

- fiel aos códigos originais em MATLAB;
- didática para estudo por capítulo;
- consistente em build, execução e pós-processamento.

## Fase 1 - Adequação da Base

Objetivo: padronizar a infraestrutura mínima do repositório.

- unificar a convenção de diretórios por capítulo: `src/`, `include/`, `scripts/`, `out/`, `README.md`, `CMakeLists.txt`;
- alinhar nomes de executáveis e arquivos principais com os demos originais;
- padronizar a forma de escrever saídas em `out/`;
- remover padrões frágeis, como `#include` de arquivos `.cpp` em vez de cabeçalhos;
- manter os arquivos MATLAB originais apenas como referência local, fora do fluxo normal de versionamento.

## Fase 2 - Padronização da Documentação

Objetivo: deixar cada capítulo autoexplicativo.

Template desejado por capítulo:

1. visão geral do problema físico;
2. teoria e equações principais;
3. mapa dos arquivos MATLAB originais;
4. mapa dos arquivos C++ equivalentes;
5. instruções de compilação;
6. instruções de execução;
7. arquivos de saída gerados;
8. scripts de visualização;
9. resultados esperados;
10. limitações conhecidas e diferenças intencionais.

## Fase 3 - Validação por Fidelidade

Objetivo: formalizar a equivalência entre MATLAB e C++.

- criar uma matriz 1:1 entre arquivos originais e traduzidos;
- marcar cada item como `traduzido`, `traduzido com desvios intencionais` ou `faltando`;
- validar parâmetros, malhas, resultados numéricos e figuras;
- registrar diferenças de implementação que sejam necessárias em C++ mas não alterem a física.

## Fase 4 - Continuação Técnica

Prioridade atual:

- `Cap_03`: completar a parte de validação do capítulo, especialmente os auxiliares do 2D/PML e os trechos ainda marcados como não testados;
- `Cap_10`: expandir a documentação didática e consolidar a validação numérica dos solvers principais;
- depois disso, abrir a próxima frente de tradução em capítulos ainda ausentes do repositório.

## Critério de Qualidade

Uma tradução só deve ser considerada pronta quando atender aos três critérios abaixo:

- compila e roda sem intervenção manual inesperada;
- reproduz o comportamento físico esperado do script MATLAB de origem;
- possui documentação suficiente para alguém entender, executar e interpretar os resultados.
