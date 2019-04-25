library(lubridate)
library(grid)
library(gridExtra)
library(scales)
library(rattle)
library(randomForest)
library(caret)
library(Hmisc)
library(e1071)
library(corrplot)
library(dplyr)
library(miscTools)
library(ggplot2)
theme_set(theme_bw(base_size=12)) # set default ggplot2 theme
library(dplyr)


#Importação da base íris
biopsy <- read.csv("http://wilkelab.org/classes/SDS348/data_sets/biopsy.csv")
head(biopsy)

iris %>% select(-Species) %>% # remove as colunas das espécies
  kmeans(centers=3) ->        # k-means clustering com 3 centros
  km                          # armazena os resultados em `km`

now display the results from the analysis

K-means clustering with 3 clusters of sizes 50, 38, 62

Cluster means:
  Sepal.Length Sepal.Width Petal.Length Petal.Width
1     5.006000    3.428000     1.462000    0.246000
2     6.850000    3.073684     5.742105    2.071053
3     5.901613    2.748387     4.393548    1.433871

Clustering vector:
  [1]   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
  [36]  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
  [71]  3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 2 2 2
  [106] 2 3 2 2 2 2 2 2 3 3 2 2 2 2 3 2 3 2 3 2 2 3 3 2 2 2 2 2 3 2 2 2 2 3 2
  [141] 2 2 3 2 2 2 3 2 2 3

Within cluster sum of squares by cluster:
  [1] 15.15100 23.87947 39.82097
(between_SS / total_SS =  88.4 %)
Available components:
  
  [1] "cluster"      "centers"      "totss"        "withinss"    
  [5] "tot.withinss" "betweenss"    "size"         "iter"        
  [9] "ifault"

km$centers

Sepal.Length Sepal.Width Petal.Length Petal.Width
1     5.006000    3.428000     1.462000    0.246000
2     6.850000    3.073684     5.742105    2.071053
3     5.901613    2.748387     4.393548    1.433871

km$cluster

[1]   1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
[36]  1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3
[71]  3 3 3 3 3 3 3 2 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 2 3 2 2 2
[106] 2 3 2 2 2 2 2 2 3 3 2 2 2 2 3 2 3 2 3 2 2 3 3 2 2 2 2 2 3 2 2 2 2 3 2
[141] 2 2 3 2 2 2 3 2 2 3

iris_clustered <- data.frame(iris, cluster=factor(km$cluster))
ggplot(iris_clustered, aes(x=Petal.Width, y=Sepal.Width, color=cluster, shape=Species)) + geom_point()

iris %>% select(-Species) %>% kmeans(centers=3) -> km
iris_clustered <- data.frame(iris, cluster=factor(km$cluster))
ggplot(iris_clustered, aes(x=Petal.Width, y=Sepal.Width, color=cluster, shape=Species)) + geom_point()

iris %>% select(-Species) %>% kmeans(centers=3, nstart=10) -> km
iris_clustered <- data.frame(iris, cluster=factor(km$cluster))
ggplot(iris_clustered, aes(x=Petal.Width, y=Sepal.Width, color=cluster, shape=Species)) + geom_point()

iris_numeric <- select(iris, -Species) # make a data table with only the numeric measurements from iris
wss <- (nrow(iris_numeric)-1)*sum(apply(iris_numeric,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(iris_numeric,
                                     nstart=10,
                                     centers=i)$withinss)
wss_data <- data.frame( centers=1:15, wss)
ggplot(wss_data, aes(x=centers, y=wss)) + geom_point() + geom_line() +
  xlab("Number of Clusters") + ylab("Within groups sum of squares")

clump_thickness uniform_cell_size uniform_cell_shape marg_adhesion
1               5                 1                  1             1
2               5                 4                  4             5
3               3                 1                  1             1
4               6                 8                  8             1
5               4                 1                  1             3
6               8                10                 10             8
epithelial_cell_size bare_nuclei bland_chromatin normal_nucleoli mitoses
1                    2           1               3               1       1
2                    7          10               3               2       1
3                    2           2               3               1       1
4                    3           4               3               7       1
5                    2           1               3               1       1
6                    7          10               9               7       1
outcome
1    benign
2    benign
3    benign
4    benign
5    benign
6 malignant
biopsy %>% select(-outcome) %>% scale() %>% prcomp() -> pca

biopsy %>% select(-outcome) %>% kmeans(centers=2, nstart=10) -> km
cluster_data <- data.frame(pca$x, cluster=factor(km$cluster), outcome=biopsy$outcome)
ggplot(cluster_data, aes(x=PC1, y=PC2, color=cluster, shape=outcome)) + geom_point()

biopsy %>% select(-outcome) %>% kmeans(centers=3, nstart=10) -> km
cluster_data <- data.frame(pca$x, cluster=factor(km$cluster), outcome=biopsy$outcome)
ggplot(cluster_data, aes(x=PC1, y=PC2, color=cluster, shape=outcome)) + geom_point()

biopsy %>% select(-outcome) %>% kmeans(centers=4, nstart=10) -> km
cluster_data <- data.frame(pca$x, cluster=factor(km$cluster), outcome=biopsy$outcome)
ggplot(cluster_data, aes(x=PC1, y=PC2, color=cluster, shape=outcome)) + geom_point()

biopsy_numeric <- select(biopsy, -outcome) # faz uma tabela de dados com apenas as medidas numéricas de biopsy
wss <- (nrow(biopsy_numeric)-1)*sum(apply(biopsy_numeric,2,var))
for (i in 2:15) wss[i] <- sum(kmeans(biopsy_numeric,
                                     nstart=10,
                                     centers=i)$withinss)
wss_data <- data.frame( centers=1:15, wss)
ggplot(wss_data, aes(x=centers, y=wss)) + geom_point() + geom_line() +
  xlab("Número de Clusters") + ylab("soma dos quadrados dentro dos grupos")

iris %>% select(-Species) %>% kmeans(centers=3, nstart=10) -> km
iris_clustered <- data.frame(iris, cluster=factor(km$cluster))
ggplot(iris_clustered, aes(x=Petal.Width, y=Sepal.Width, color=cluster, shape=Species)) + geom_point()

centroids <- data.frame(km$centers)
centroids

Sepal.Length Sepal.Width Petal.Length Petal.Width
1     5.901613    2.748387     4.393548    1.433871
2     5.006000    3.428000     1.462000    0.246000
3     6.850000    3.073684     5.742105    2.071053

then we need to add a column that indicates the cluster number. Again we use
factor() to express that the cluster numbers are discrete categories.
centroids <- data.frame(centroids, cluster=factor(1:3))
now we can plot

ggplot(iris_clustered, aes(x=Petal.Width, y=Sepal.Width, color=cluster)) + 
  geom_point(aes(shape=Species)) + # pontos individuais do quadro de dados  `iris_clustered`
  geom_point(data=centroids, aes(fill=cluster), shape=21, color="black", size=4, stroke=1) # centroids


ggplot(iris, aes(x=Species, fill=Species)) +geom_bar()

ggplot(iris, aes(x=Sepal.Length, y=Petal.Length)) +geom_point(aes(color=Species)) +geom_smooth(aes(group=Species))

ggplot(iris, aes(x=Sepal.Length, y=Petal.Length)) + stat_smooth(aes(group=Species), geom="ribbon", fill='gray70') +geom_point(aes(color=Species))

ggplot(iris, aes(x=Sepal.Length, y=Petal.Length)) + geom_point()+ scale_x_log10()

ggplot(iris, aes(x=Species, fill=Species))+geom_bar() + scale_fill_grey()

ggplot(iris, aes(x=Species, fill=Species))+geom_bar()+ scale_fill_brewer()


par(mfrow = c(1,2), oma = c(4,1,1,1))

plot(iris$Sepal.Length, iris$Sepal.Width)

plot(iris$Petal.Length, iris$Petal.Width)

par(mfrow = c(1,2), oma = c(4,1,1,1))

plot(iris$Sepal.Length, iris$Sepal.Width, xlab = "Comprimento da Sépala", ylab = "Largura da Sépala")

plot(iris$Petal.Length, iris$Petal.Width, xlab = "Comprimento da Pétala", ylab = "Largura da Pétala")






