# Regressão Birnbaum-Saunders Generalizado via GAMLSS

Códigos desenvolvidos para a dissertação de mestrado com o título "Novo Modelo de Regressão Birnbaum-Saunders Gerenalizado: Uma Abordagem via GAMLSS". Neste repositório você encontra a implementação computacional do novo modelo de regressão Birnbaum-Saunders Generalizado (BSG) proposto via GAMLSS, a partir da biblioteca do R de mesmo nome. 

## Códigos Principais

* **BSG_GAMLSS**: implementação computacional da distribuição BSG na arquitetura da biblioteca **gamlss** do R.
* **SIMULATION_SCENARIO**: script R que resulta o procedimento de simulação considerado na dissertação.
* **SAMPLE**: script R que resulta nas matrizes designs e variável resposta.
* **BOXPLOT_ESTIMATES**: script R que resulta os _boxplots_ das estimativas dos coeficientes de regressão.

## Simulação

Os diretórios **CENÁRIO X** correspondem aos cenários de simulação de considerados, em que X identifica o cenário e pode assumir valores de 1 a 6. Além disso, os cenários 1 e 2 apresentam três subdivisões cada. Em cada diretório contém os seguintes arquivos: **BSG_GAMLSS**, **SIMULATION_SCENARIO**, **SAMPLE** e **BOXPLOT_ESTIMATES**. Além disso, nestes diretórios também estão os resultados das 5000 réplicas para cada tamanho amostral considerado, bem como uma arquivo contendo um resumo do procecimento de simulação.


## Aplicação

Os códigos utilizados para executar as aplicações a dados reais estão no diretório [APLICAÇÃO](https://github.com/braga0m/bsg_gamlss/tree/main/APLICA%C3%87%C3%83O). Neste diretório, estão os gráficos utilizados em cada uma das aplicações,  estão contidos nas pastas
"CHEESE" E "DENTAL". O script R **mods_dissertacao_final** é o código utilizado na modelagem dos dados, enquanto **descritivas** é a análise descritiva.

##  Autor

* **[Matheus Braga Milhomem](https://github.com/braga0m)**

## Colaboração

* **[Terezinha Kessia de Assis Ribeiro](https://github.com/terezinharibeiro/)** como orientadora da dissertação de mestrado.
