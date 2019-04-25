# DATA CLUSTERING
# Disciplina de Exploração e Mineração de Dados UFAL -2018.2

# INTRODUÇÃO
Problemas de otimização combinatória são difíceis de lidar e resolver, para encontrar as
melhores soluções entre as alternativas para um determinado problema. Todos os anos, suas
aplicações na ciência e na indústria estão se tornando cada vez maiores e complexas,
impulsionando necessidade de métodos mais poderosos e especializados. Tradicionalmente,
existem duas abordagens principais usadas para resolver esses problemas: métodos exatos e
heurísticos.

Neste trabalho, o enfoque é analisar a paralelização em Graphics Processing Unit (GPU)
da meta-heurística GRASP para o Problemas de Clusterização Automática (PCA).
A típica análise de clusterização consiste de quatro etapas que são o recurso seleção, que
é um problema de otimização, onde o problema é pesquisar através do espaço de subconjuntos
de recursos para identificar o ideal ou próximo do ideal em relação a um medida de
desempenho.

(I) Feature Extraction: Onde algumas transformações são usados para gerar recursos úteis
e inovadores dos originais.

(II) The Clustering Algorithm Design: Que geralmente é combinado com a seleção de uma
correspondente medida de proximidade.

(III) Cluster Validation: Os índices externos, índices internos e índices relativos são usados
para análise de validade dos cluster.

(IV) Results Interpretation: São selecionados os especialistas nos campos mais relevantes, e depois  interpretam a partição
de dados a fim de garantir a confiabilidade da extração que foi extraída na base de conhecimento.

De acordo com Satoru Ochi (2015), o Problema de Clusterização Automática (PCA) é apresentado
da seguinte forma: dado um conjunto X de n objetos X = {x1, x2, x3, ... , Xn}, onde cada objeto Xi é uma tupla (xi1, xi2, ..., Xip), e cada coordenada xij está relacionada a um atributo do objeto (cada objeto pode ser um ponto p no espaço. O objetivo é encontrar uma partição X = {C1,C2, ... ,Ck} de clusters, onde K não é conhecido, de modo que a similaridade entre os objetos do mesmo cluster são maximizados e, a similaridade entre objetos de diferentes clusters é minimizada, sob as seguintes condições adicionais:

 Formulação matemáticado problema segue abaixo:
 
![](https://github.com/uncisal/Dataclustering/blob/master/clusters.PNG)

Onde a formulação matemática constitui de $C = \{C_1, C_2,..., Ck} é uma partição particular de clusters e C é um conjunto de todas as partições possíveis do conjunto X, com k = 2,…, n-1 e S(Xi) é um valor de Silhouete de cada ponto Xi  pertence a X.
A função, denominada Índice de Silhouete, foi proposta por (Kaufman e Rousseeuw, 2009). Por fim retorna os valores dentro do intervalo [-1,1] e não é necessário definir o valor de k,
porque o valor mais apropriado de k é alcançado maximizando essa função.

  Foi utilizados dois algoritmos bases para resolução do Problema de Clusterização  automática (K-means e GRASP)
abaixo segue uma breve explicação dos dois algoritmos.

  Em mineração de dados, agrupamento k-means é um método de Clustering que objetiva particionar n observações dentre k grupos onde cada observação pertence ao grupo mais próximo da média. Isso resulta em uma divisão do espaço de dados em um Diagrama de Voronoi.
   O problema é computacionalmente difícil (NP-difícil), no entanto, existem algoritmos heurísticos eficientes que são comumente empregados e convergem rapidamente para um local ótimo. Estes são geralmente semelhantes ao algoritmo de maximização da expectativa para misturas de distribuições gaussianas através de uma abordagem de refinamento iterativo utilizado por ambos os algoritmos. Além disso, ambos usam os centros de clusters para modelar dados, no entanto, a clusterização k-means tende a encontrar clusters de extensão espacial comparáveis enquanto o mecanismo de maximização da expectativa permite ter diferentes formas. 
   
   A meta-heurística GRASP (do inglês Greedy Randomized Adaptative Search Procedure) é um procedimento iterativo, em que cada iteração do GRASP constitui de duas fases: uma fase de construção e uma fase de busca local (FEO e RESENDE, 1995). A melhor solução global é mantida como resultado.
   
   Na Figura 1 é apresentado um pseudocódigo do algoritmo GRASP. Na linha 1, é inicializada a instância do problema. A partir das instâncias, é feita a construção das soluções e a busca local, onde serão atualizados as melhores soluções encontradas dentro das instâncias. Por fim, será retornado o melhor valor da solução que foi encontrada.
   
   ![](https://github.com/uncisal/Dataclustering/blob/master/grasp1.jpg)
   
   Figura 1: Um pseudo-código genérico GRASP
   
Por sua vez, na Figura 2 é feita a fase de construção de novas soluções. Nesta etapa, é criada uma lista restrita de candidatos e atribui os melhores valores das soluções da lista restrita candidatos a uma função gulosa.
   
   ![](https://github.com/uncisal/Dataclustering/blob/master/grasp2.jpg)
   
   Figura 2: Fase de construção do Pseudo-código GRASP
   
   A eficácia destes métodos depende da sua capacidade para se adaptar a uma realização particular, evitar o aprisionamento em ótimos locais e explorar a estrutura básica do problema, como por exemplo uma rede ou um ordenamento natural entre seus componentes. Além disso, procedimentos como: reinicialização, randomização controlada, estruturas de dados eficientes, e pré-processamento também são benéficos.
# METODOLOGIA
As ferramentas que foram utilizadas nesse projeto foi: 

Para fase inicial da pesquisa foi  utilizado o ambiente "RStudio" que é um software livre de ambiente de desenvolvimento integrado para R, uma linguagem de programação para gráficos e cálculos estatísticos.

Foi utilizada a linguagem de programação C++ (gcc - 4.8.4) para implementação da
solução proposta, além da biblioteca CPLEX 12.5, que é um dos pacotes de software de
otimização linear mista mais utilizados na literatura.

Compute Unified Device Architecture (CUDA), que é uma plataforma de programação paralela projetada pela NVIDIA para suas GPUs, podendo ser usado para aumentar o poder computacional em aplicações de propósitos gerais. O CUDA fornece um pequeno conjunto de extensões para linguagens de programação
padrão ( A linguagem utilizada foi C++).



Foi utilizada a base de dados Íris composta por um conjunto de dados que contém 3 classes de 50 instâncias cada, onde cada classe se refere a um tipo de planta da íris. Os tipos são (Iris Setosa) (Iris Versicolour) (Iris Virginica). Para mais informações acesse: http://archive.ics.uci.edu/ml/datasets/iris  e http://wilkelab.org/classes/SDS348/data_sets/biopsy.csv

# RESULTADOS ATRAVEZ DA LINGUAGEM R E C++

Primeiro será descrito os resultados obtidos na linguagem R

Representação da base íris.

<img src=https://github.com/uncisal/Dataclustering/blob/master/representacao%20da%20base.PNG width="700" height="500">

Fase separação da base íris.

<img src=https://github.com/uncisal/Dataclustering/blob/master/fase%20separacao.PNG width="700" height="500">

Centroides da base Íris

<img src=https://github.com/uncisal/Dataclustering/blob/master/centroides%20da%20%C3%ADris.PNG width="700" height="500">

Plots finais da base Íris segue na ordem: Correlação da base, números de clusters, média, mediana, quartil e amplitude, serparação das espécies.

<img src=https://github.com/uncisal/Dataclustering/blob/master/plots%20finais.PNG width="700" height="500">

Resultados obtidos com C++  usando a plataforma CUDA em GPU 

Interação 01

<img src=https://github.com/uncisal/Dataclustering/blob/master/intera%C3%A7%C3%A3o%2001.PNG width="400" height="250">

Interação 02

<img src=https://github.com/uncisal/Dataclustering/blob/master/intera%C3%A7%C3%A3o%2002.PNG width="400" height="250">

Resultados das interações:

<img src=https://github.com/uncisal/Dataclustering/blob/master/resultados%20%20em%20GPU.PNG width="450" height="300">


# CONCLUSÃO

Este trabalho propôs uma visão de como é feita a clusterização de dados utlizando o Dataset ÍRIS, como podemos observar a base pode ser bem explorada pois é bem organizada e seus dados são bem consistes e organizados. Devido a isso, podemos entender o comportamento da base através dos algoritmos que foram utilizados tanto o (K-means e GRASP) que foram bem sucetives a base. Por fim os resultados que foram dados pela meta-heurística GRASP mostra que o ganho de desempenho é bastante interessante, pois a velocidade em GPU mostrou-se bastante eficaz para fazer a clusterização dos dodos. 


# REFERÊNCIAS

Dib Cruz, Marcelo & Ochi, Luiz. (2015). A Hybrid Method Using Evolutionary and a Linear Integer Model to Solve the Automatic Clustering Problem. Learning and Nonlinear Models. 13. 73-91. 10.21528/LNLM-vol13-no2-art2. 

FEO.T.A., RESENDE, M. G. C. Greedy Randomized Adaptive Search Procedures. Journal of Global Optimization, v. 6, p. 109–134, 1995. 

Nogueira, B., Tavares, E., Araujo, J., & Callou, G. (2019). Accelerating continuous GRASP with a GPU. Journal of Supercomputing.

The wilke Lab: Computational Evolutionary Biology disponível em: https://wilkelab.org/

UCI:  Center for Machine Learning and Intelligent Systems. Machine Leaning Repository Disponível em:https://archive.ics.uci.edu/ml/datasets/iris


