require(gdata);
ATV = read.xls ("ATV_Vectors.xlsx", sheet = 1, header = TRUE);
set.seed(9850);
gp <- runif(nrow(ATV));
ATV_Ordered <- ATV[order(gp),] //Sort on order of gp
options(scipen = 999)

//Normalization not needed here an optional so not applied following:

normalize <- function(x) {
	return((x-min(x)) / (max(x) - min(x)))
}
ATV_Ordered <- as.data.frame(lapply(ATV_Ordered[,c(3:212)],normalize))

////////////
ATV_Ordered_dataonly <- ATV_Ordered[,c(3:212)]
ATV_Ordered_train <- ATV_Ordered_dataonly[1:2000,]
ATV_Ordered_test <- ATV_Ordered_dataonly[2001:8430,]
ATV_Ordered_train_target <- ATV_Ordered[1: 2000, 1]
ATV_Ordered_test_target <- ATV_Ordered[2001:8430, 1]

require(class)
m1 <- knn(train=ATV_Ordered_train, test=ATV_Ordered_test, cl= ATV_Ordered_train_target, k=5)
table(ATV_Ordered_test_target, m1)

//K Means Clustering

cluster <- kmeans(ATV_Ordered[,3:212], 10, nstart=20)
table(cluster$cluster, ATV_Ordered$Identifiers)


//Selecting optimum number of clusters using elbow method

set.seed(123)
k.max <- 155
data <- ATV_Ordered[,3:212]
wss <- sapply(145:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(145:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


tree <- read.csv(file="Selected_Mutants_ATV.csv",header=TRUE,sep=",") 

tree <- as.data.frame(tree)

colors <- colors[as.numeric(tree$Kvalue)]

colors <- c("#999999", "#E69F00", "#56B4E9", '#ff0000')

png(file="ANN_ATV_2.png",width=400,height=400,res=72)

scatterplot3d(tree[,1:2], pch = 16, color=colors,  main="ATV KNN CLUSTERING",
              xlab = "Resistance (Fold)",
               ylab = "Number of Mutations"
               )

legend("bottom", legend = levels(tree$Kvalue),
       col =  c("#999999", "#E69F00", "#56B4E9", '#ff0000'), 
       pch = c(16, 17, 18,19), 
       inset = -0.25, xpd = TRUE, horiz = TRUE)
       
       
//K-Means Clustering updated start from below

ATV <- read.csv(file="ATV_Vectors.csv",header=TRUE,sep=",") 
DRV <- read.csv(file="DRV_Vectors.csv",header=TRUE,sep=",") 
FPV <- read.csv(file="FPV_Vectors.csv",header=TRUE,sep=",") 
IDV <- read.csv(file="IDV_Vectors.csv",header=TRUE,sep=",") 
LPV <- read.csv(file="LPV_Vectors.csv",header=TRUE,sep=",") 
NFV <- read.csv(file="NFV_Vectors.csv",header=TRUE,sep=",") 
SQV <- read.csv(file="SQV_Vectors.csv",header=TRUE,sep=",") 
TPV <- read.csv(file="TPV_Vectors.csv",header=TRUE,sep=",") 

//Selecting optimum number of clusters 

//Method 1: Using elbow method

set.seed(123)
k.max <- 15
data <- TPV[,4:213]
wss <- sapply(1:k.max, 
              function(k){kmeans(data, k, nstart=50,iter.max = 15 )$tot.withinss})
plot(1:k.max, wss,type="b", pch = 19, frame = FALSE, xlab="Number of clusters K", ylab="Total within-clusters sum of squares")
abline(v = 3, lty = 2)


//Method 2: Using clusGap

library(cluster)
ATV_opt <- clusGap(ATV[,4:213], kmeans, 160, B = 100, verbose = interactive())
DRV_opt <- clusGap(DRV[,4:213], kmeans, 160, B = 100, verbose = interactive())
FPV_opt <- clusGap(FPV[,4:213], kmeans, 160, B = 100, verbose = interactive())
IDV_opt <- clusGap(IDV[,4:213], kmeans, 160, B = 100, verbose = interactive())
LPV_opt <- clusGap(LPV[,4:213], kmeans, 160, B = 100, verbose = interactive())
NFV_opt <- clusGap(NFV[,4:213], kmeans, 160, B = 100, verbose = interactive())
SQV_opt <- clusGap(SQV[,4:213], kmeans, 160, B = 100, verbose = interactive())
TPV_opt <- clusGap(TPV[,4:213], kmeans, 160, B = 100, verbose = interactive())

print(ATV_opt, method = "firstmax")
plot(ATV_opt, frame = FALSE, xlab = "Number of clusters k")
abline(v = 150, lty = 2)


//Methods 3 and 4 Didnt Used Memory inefficient: 

d_clust <- Mclust(as.matrix(ATV[,4:213]), G=1:50, modelNames = mclust.options("emModelNames"))

nb <- NbClust(as.matrix(ATV[,4:213]), diss="NULL", distance = "euclidean", min.nc=2, max.nc=15, method = "kmeans", index = "alllong", alphaBeale = 0.1)

     
//Clustering based on optimum clusters

ATV_cluster <- kmeans(ATV[,4:213], 3, nstart=20)
DRV_cluster <- kmeans(DRV[,4:213], 3, nstart=20)
FPV_cluster <- kmeans(FPV[,4:213], 4, nstart=20)
IDV_cluster <- kmeans(IDV[,4:213], 3, nstart=20)
LPV_cluster <- kmeans(LPV[,4:213], 4, nstart=20)
NFV_cluster <- kmeans(NFV[,4:213], 4, nstart=20)
SQV_cluster <- kmeans(SQV[,4:213], 4, nstart=20)
TPV_cluster <- kmeans(TPV[,4:213], 3, nstart=20)

//Tables to get accuracy predictions of clusters

table(ATV_cluster$cluster, ATV$class)
table(DRV_cluster$cluster, DRV$class)
table(FPV_cluster$cluster, FPV$class)
table(IDV_cluster$cluster, IDV$class)
table(LPV_cluster$cluster, LPV$class)
table(NFV_cluster$cluster, NFV$class)
table(SQV_cluster$cluster, SQV$class)
table(TPV_cluster$cluster, TPV$class)

//Getting all clustering vectors to export

ATV_Clustering_vector <- data.frame(ATV, ATV_cluster$cluster) 
DRV_Clustering_vector <- data.frame(DRV, DRV_cluster$cluster)
FPV_Clustering_vector <- data.frame(FPV, FPV_cluster$cluster)
IDV_Clustering_vector <- data.frame(IDV, IDV_cluster$cluster)
LPV_Clustering_vector <- data.frame(LPV, LPV_cluster$cluster)
NFV_Clustering_vector <- data.frame(NFV, NFV_cluster$cluster) 
SQV_Clustering_vector <- data.frame(SQV, SQV_cluster$cluster)
TPV_Clustering_vector <- data.frame(TPV, TPV_cluster$cluster)      

write.csv(ATV_Clustering_vector, file = "ATV_Clustering_vector.csv")
write.csv(DRV_Clustering_vector, file = "DRV_Clustering_vector.csv")
write.csv(FPV_Clustering_vector, file = "FPV_Clustering_vector.csv")
write.csv(IDV_Clustering_vector, file = "IDV_Clustering_vector.csv")
write.csv(LPV_Clustering_vector, file = "LPV_Clustering_vector.csv")
write.csv(NFV_Clustering_vector, file = "NFV_Clustering_vector.csv")
write.csv(SQV_Clustering_vector, file = "SQV_Clustering_vector.csv")
write.csv(TPV_Clustering_vector, file = "TPV_Clustering_vector.csv")

//Hierarchial Clustering

ATV_dist = dist(ATV[,4:213])
ATV.hclust = hclust(ATV_dist)
groups.6 = cutree(ATV.hclust,6)
#plot(ATV.hclust)
table(groups.6)
ATV_Hierarchial_Clustering <- data.frame(ATV, groups.6)


DRV_dist = dist(DRV[,4:213])
DRV.hclust = hclust(DRV_dist)
groups.6 = cutree(DRV.hclust,6)
#plot(DRV.hclust)
table(groups.6)
DRV_Hierarchial_Clustering <- data.frame(DRV, groups.6)


FPV_dist = dist(FPV[,4:213])
FPV.hclust = hclust(FPV_dist)
groups.6 = cutree(FPV.hclust,6)
#plot(FPV.hclust)
table(groups.6)
FPV_Hierarchial_Clustering <- data.frame(FPV, groups.6)


IDV_dist = dist(IDV[,4:213])
IDV.hclust = hclust(IDV_dist)
groups.6 = cutree(IDV.hclust,6)
#plot(IDV.hclust)
table(groups.6)
IDV_Hierarchial_Clustering <- data.frame(IDV, groups.6)


LPV_dist = dist(LPV[,4:213])
LPV.hclust = hclust(LPV_dist)
groups.6 = cutree(LPV.hclust,6)
#plot(LPV.hclust)
table(groups.6)
LPV_Hierarchial_Clustering <- data.frame(LPV, groups.6)


NFV_dist = dist(NFV[,4:213])
NFV.hclust = hclust(NFV_dist)
groups.6 = cutree(NFV.hclust,6)
#plot(NFV.hclust)
table(groups.6)
NFV_Hierarchial_Clustering <- data.frame(NFV, groups.6)


SQV_dist = dist(SQV[,4:213])
SQV.hclust = hclust(SQV_dist)
groups.6 = cutree(SQV.hclust,6)
#plot(SQV.hclust)
table(groups.6)
SQV_Hierarchial_Clustering <- data.frame(SQV, groups.6)


TPV_dist = dist(TPV[,4:213])
TPV.hclust = hclust(TPV_dist)
groups.6 = cutree(TPV.hclust,6)
#plot(TPV.hclust)
table(groups.6)
TPV_Hierarchial_Clustering <- data.frame(TPV, groups.6)

write.csv(ATV_Hierarchial_Clustering, file = "ATV_Hierarchial_Clustering.csv")
write.csv(DRV_Hierarchial_Clustering, file = "DRV_Hierarchial_Clustering.csv")
write.csv(FPV_Hierarchial_Clustering, file = "FPV_Hierarchial_Clustering.csv")
write.csv(IDV_Hierarchial_Clustering, file = "IDV_Hierarchial_Clustering.csv")
write.csv(LPV_Hierarchial_Clustering, file = "LPV_Hierarchial_Clustering.csv")
write.csv(NFV_Hierarchial_Clustering, file = "NFV_Hierarchial_Clustering.csv")
write.csv(SQV_Hierarchial_Clustering, file = "SQV_Hierarchial_Clustering.csv")
write.csv(TPV_Hierarchial_Clustering, file = "TPV_Hierarchial_Clustering.csv")


//AGNES: Agglomerative Nesting TRY DIVISIVE INSTEAD!!

ATV_dv <- diana(ATV_dist, diss = TRUE, stand = FALSE)
ATV.hclust = hclust(ATV_dv)


library(cluster)
ATV_dist = dist(ATV[,4:213])
ATV.agnes = agnes(ATV_dist)
table(cutree(ATV.agnes,6))

DRV_dist = dist(DRV[,4:213])
DRV.agnes = agnes(DRV_dist)
table(cutree(DRV.agnes,6))

FPV_dist = dist(FPV[,4:213])
FPV.agnes = agnes(FPV_dist)
table(cutree(FPV.agnes,6))

IDV_dist = dist(IDV[,4:213])
IDV.agnes = agnes(IDV_dist)
table(cutree(IDV.agnes,6))

LPV_dist = dist(LPV[,4:213])
LPV.agnes = agnes(LPV_dist)
table(cutree(LPV.agnes,6))

NFV_dist = dist(NFV[,4:213])
NFV.agnes = agnes(NFV_dist)
table(cutree(NFV.agnes,6))

SQV_dist = dist(SQV[,4:213])
SQV.agnes = agnes(SQV_dist)
table(cutree(SQV.agnes,6))

TPV_dist = dist(TPV[,4:213])
TPV.agnes = agnes(TPV_dist)
table(cutree(TPV.agnes,6))



























      
       