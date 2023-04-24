#Loading the relevant libraries
library(MASS)
library(ggplot2)
library(PlackettLuce)
library(PLMIX)
library(BayesMallows)


#Removing names from the data
fruittable <- as.data.frame(Ranking_Table)


#Basic statistics; number of people asked
nrow(fruittable)

#Plot displaying the proportions of Genders from either location
genloc <- ggplot(fruittable, aes(x = Location, fill = Gender)) + 
  ylab("Frequency")
genloc + geom_bar(position = "fill") + 
  scale_x_discrete(name ="Location", labels=c("Exeter","Oxford"))

#Assessing the frequency of people of a given gender from the locations
lwgenloc <- function(array, loc, gen){
  return(length(which(array$Location==loc & array$Gender==gen)))
}

lwgenloc(fruittable, "E", "M")
lwgenloc(fruittable, "E", "F")
lwgenloc(fruittable, "O", "M")
lwgenloc(fruittable, "O", "F")

#Removing gender identity and location from the different data tables
fruits <- fruittable[-c(1,2)]

#Producing a plot of the rank totals
sums <- as.data.frame(colSums(fruits))
sums$key <- colnames(fruits)
sums$key <- factor(sums$key, levels=unique(sums$key))
ggplot(sums, aes(x=key, y=colSums(fruits), fill = key)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Fruits") + ylab("Rank Totals") + theme(legend.position = "none")

#Producing a mosaic plot of the data
mosfruit <- t(cbind(table(fruits$Apples), table(fruits$Oranges), 
                    table(fruits$Bananas), table(fruits$Grapes), 
                    table(fruits$Raspberries), table(fruits$Strawberries),
                    table(fruits$Pineapples), table(fruits$Blueberries),
                    table(fruits$Cherries), table(fruits$Peaches)))

rownames(mosfruit) <- colnames(fruits)
mosaicplot(mosfruit, main = "", xlab = "Fruit", ylab = "Rank Given", 
           color = TRUE, cex.axis = 0.9, las = 2)


#Creating an explicit ordering of each assessors rankings
mrank <- apply(as.matrix.noquote(fruits),2,as.numeric)
ranks <- as.rankings(mrank)
ranks



#####
#FOLLOWING USED FOR COMPARISON ONLY
#Fitting a Plackett-Luce Model using MLE via the PlackettLuce package
fruitPL <- PlackettLuce(ranks, npseudo = 0)
summary(fruitPL)

#Getting the fitted v_i values for each object
coef(fruitPL, log = FALSE)

#Calculation and plot of the quasi variances of the logs of the parameters
qv <- qvcalc(fruitPL)
plot(qv, ylab = "Worth (log)", main = NULL)

#Creating a function that determines the probability of seeing a ranking given
#a set of Plackett-Luce parameters
PL_prob <- function(rankvec, par) {
  if (length(rankvec) != length(par)) {
    print("Please input a vector with the correct number of items")
  } else {
    par <- unname(par)
    pars <- par
    tot <- 1
    for (x in 1:length(rankvec)) {
      tot <- (par[which(rankvec == x)]/sum(pars)) * tot
      pars[which(rankvec == x)] <- 0
    }
    return(tot)
  }
}

#Given that there are 10! possible rankings, these probabilities are bound to
#be very low, though comparison of magnitudes may be helpful
PL_prob(c(8,7,6,2,3,1,4,10,9,5), coef(fruitPL, log = FALSE))
PL_prob(c(3,4,5,9,8,10,7,1,2,6), coef(fruitPL, log = FALSE))
PL_prob(c(1,2,3,4,5,6,7,8,9,10), coef(fruitPL, log = FALSE))





#Using the PLMIX package to get object parameters via MAP estimation and Gibbs
#Sampling

#Ordering the items according to their ranks
RR <- as.top_ordering(ranks)
asnum <- nrow(RR)

#Function that returns a Mallows Model or Mallows Mixture Model with mixnum 
#cluster and Gamma prior
nPLMix <- function(ord, mixnum, nassess, shape, rate, alpha) {
  if (shape == 0 & rate == 0) {
    print("Uniform priors chosen")
    return(mapPLMIX(ord, K=ncol(ord), G=mixnum, plot_objective = TRUE))
  } else {
    return(mapPLMIX(ord, K=ncol(ord), G=mixnum, 
                    init = list(p = NULL, omega = rep(1/mixnum, mixnum)),
                    hyper = list(shape0 = matrix(shape, nrow = mixnum, 
                                                 ncol = nassess), 
                                 rate0 = rep(rate, mixnum),
                                 alpha0 = rep(alpha, mixnum)), 
                    plot_objective = TRUE))
  }
}

#Non-Mixture PL with flat prior
MFlat <- nPLMix(RR, 1, asnum, 0, 0, 1)
summary(MFlat)
MFlat$P

#Non-Mixture PL with Gamma(2, 1) priors
MPL <- nPLMix(RR, 1, asnum, 2, 1, 1) 
summary(MPL)
MPL$P

#With 2 mixture components, initially equally weighted, with Gamma(2,1) priors
MPL2 <- nPLMix(RR, 2, asnum, 2, 1, 1) 
summary(MPL2)
MPL2$P
MPL2$class_map

#With 3 mixture components
MPL3 <- nPLMix(RR, 3, asnum, 2, 1, 1) 
MPL3$P
MPL3$class_map

#With 4 mixture components
MPL4 <- nPLMix(RR, 4, asnum, 2, 1, 1) 
MPL4$P
MPL4$class_map

# With 5 mixture components
MPL5<- nPLMix(RR, 5, asnum, 2, 1, 1) 
MPL5$P
MPL5$class_map

#With 6 mixture components
MPL6 <- nPLMix(RR, 6, asnum, 2, 1, 1) 
MPL6$P
MPL6$class_map


#Using Gibbs sampling with previously found MAP estimates with vague priors 
#to give model parameter estimates; if no MAP estimates are supplied, assume
#randomly uniform component membership
GibbsEst <- function(ord, mixnum, MAPmod = NULL){
  if (is.null(MAPmod) == TRUE) {
    return(gibbsPLMIX(ord, K=ncol(ord), G=mixnum, init = list(p = NULL, z=NULL),
                      n_iter = 10^5))
  } else {
    return(gibbsPLMIX(ord, K=ncol(ord), G=mixnum, init = list(p = MAPmod$P, 
                           z=binary_group_ind(MAPmod$class_map, G=mixnum)), 
                      n_iter = 10^5))
  }
}

#Getting the respective Gibbs Sample estimates using MAP estimates and uniform 
#prior membership
GPL <- GibbsEst(RR, 1, MPL)
GPLU <- GibbsEst(RR, 1)
matrix(colMeans(GPL$P),ncol=ncol(RR))

GPL2 <- GibbsEst(RR, 2, MPL2)
GPL2U <- GibbsEst(RR, 2)
matrix(colMeans(GPL2$P),ncol=ncol(RR))

GPL3 <- GibbsEst(RR, 3, MPL3)
GPL3U <- GibbsEst(RR, 3)
matrix(colMeans(GPL3$P),ncol=ncol(RR))

GPL4 <- GibbsEst(RR, 4, MPL4)
GPL4U <- GibbsEst(RR, 4)
matrix(colMeans(GPL4$P),ncol=ncol(RR))

GPL5<- GibbsEst(RR, 5, MPL5)
GPL5U <- GibbsEst(RR, 5)
matrix(colMeans(GPL5$P),ncol=ncol(RR))

GPL6<- GibbsEst(RR, 6, MPL6)
GPL6U <- GibbsEst(RR, 6)
matrix(colMeans(GPL6$P),ncol=ncol(RR))




#Generating various values concerning model selection criteria, comparing the
#6 models produced above using the posterior mode (MAP) and the Gibbs Sampling 
#processes that used the MAP estimate as an initialisation
modsel <- selectPLMIX(pi_inv=RR, seq_G=1:6, parallel=FALSE, 
            MAPestP=list(MPL$P, MPL2$P, MPL3$P, MPL4$P, MPL5$P, MPL6$P), 
            MAPestW=list(MPL$W, MPL2$W, MPL3$W, MPL4$W, MPL5$W, MPL6$W),
            deviance=list(GPL$deviance, GPL2$deviance, 
                          GPL3$deviance, GPL4$deviance, 
                          GPL5$deviance, GPL6$deviance))

#Doing the same for where initial component membership and PL parameters were 
#uniformly random, with posterior summary statistic being the mean
modselMCMC <- selectPLMIX(pi_inv=RR, seq_G=1:6, parallel=FALSE, 
                       MAPestP=NULL, MAPestW=NULL, post_est = "mean",
                       MCMCsampleP =list(GPLU$P, GPL2U$P, GPL3U$P, 
                                    GPL4U$P, GPL5U$P, GPL6U$P), 
                       MCMCsampleW =list(GPLU$W, GPL2U$W, GPL3U$W, 
                                         GPL4U$W, GPL5U$W, GPL6U$W),
                       deviance=list(GPLU$deviance, GPL2U$deviance, 
                                     GPL3U$deviance, GPL4U$deviance, 
                                     GPL5U$deviance, GPL6U$deviance))

#As before but with summary statistic being the median
modselMCMCmed <- selectPLMIX(pi_inv=RR, seq_G=1:6, parallel=FALSE, 
                          MAPestP=NULL, MAPestW=NULL, post_est = "median",
                          MCMCsampleP =list(GPLU$P, GPL2U$P, GPL3U$P, 
                                            GPL4U$P, GPL5U$P, GPL6U$P), 
                          MCMCsampleW =list(GPLU$W, GPL2U$W, GPL3U$W, 
                                            GPL4U$W, GPL5U$W, GPL6U$W),
                          deviance=list(GPLU$deviance, GPL2U$deviance, 
                                        GPL3U$deviance, GPL4U$deviance, 
                                        GPL5U$deviance, GPL6U$deviance))

#Looking at the criterion explicitly
modsel$criteria
modselMCMC$criteria
modselMCMCmed$criteria

#Most criterion seem to favour 2 mixture components
classes <- as.matrix(MPL2$class_map)
rwclass <- cbind(classes, fruittable)

mg1 <- rwclass[which(rwclass$classes == 1),] 
mg2 <- rwclass[which(rwclass$classes == 2),]

#Proportion of males in the overall data set compared to the different 
#components
(lwgenloc(fruittable,"O","M") + lwgenloc(fruittable,"E","M"))/nrow(fruittable)
(lwgenloc(mg1,"O","M") + lwgenloc(mg1,"E","M"))/nrow(mg1)
(lwgenloc(mg2,"O","M") + lwgenloc(mg2,"E","M"))/nrow(mg2)

#Same with the proportion of data set taken from Oxford
(lwgenloc(fruittable,"O","M") + lwgenloc(fruittable,"O","F"))/nrow(fruittable)
(lwgenloc(mg1,"O","M") + lwgenloc(mg1,"O","F"))/nrow(mg1)
(lwgenloc(mg2,"O","M") + lwgenloc(mg2,"O","F"))/nrow(mg2)

#Resolving the label switching phenomena using the Pivotal Reordering Algorithm
labswitch2 <- label_switchPLMIX(pi_inv = RR, seq_G = 2, MCMCsampleP = list(GPL2U$P), 
                               MCMCsampleW =list(GPL2U$W), MAPestP=list(MPL2$P), 
                               MAPestW=list(MPL2$W), parallel = FALSE)

#Doing so yields PL parameter estimates
LSpars <- round(apply(labswitch2$final_sampleP$G_2, 2, rowMeans), 3)
colnames(LSpars) <- c("Apl", "Orn", "Ban", "Grp", "Rsp", "Str", 
                      "Pin", "Blu", "Chr", "Pch")
LSpars




#####
#KENDALL DISTANCE
#Computing Mallows Mixture models with Kendall Distance for between 1 and 16 
#mixture components
fruitMalkenmix <- compute_mallows_mixtures(n_clusters = c(1:16), mrank, 
                                           metric = "kendall", nmc = 10^5L, 
                                           include_wcd=TRUE)

#Elbow plot to determine optimal number of components; suggests 2
plot_elbow(fruitMalkenmix, burnin = 500)

#Creating a model with 2 mixture components
Kenranks <- compute_mallows(mrank, metric="kendall", n_clusters = 2L, 
                            save_clus = TRUE, save_aug = TRUE, nmc = 10^5)
Kenranks$burnin <- 5000


#Computing statistics of and plotting the posterior density of the shape
#parameter
compute_posterior_intervals(Kenranks, parameter = "alpha")
plot(Kenranks)

# Compute the CP consensus for each component
KenCP <- as.data.frame(compute_consensus(Kenranks, type = "CP"))
KenCP$cumprob <- NULL
KC<-stats::reshape(KenCP, direction = "wide", idvar = "ranking", 
               timevar = "cluster", varying = list(as.character(unique(KenCP$cluster))))

# Compute the MAP consensus for each component
KenMAP <- compute_consensus(Kenranks, type = "MAP")
KenMAP$probability <- NULL
KM<-stats::reshape(KenMAP, direction = "wide", idvar = "map_ranking",
               timevar = "cluster", varying = list(as.character(unique(KenMAP$cluster))))

#Combine the two to compare
Kdf <- data.frame(KC,KM)
colnames(Kdf) <- c("Ranking","C1 (CP)","C2 (CP)","","C1 (MAP)","C2 (MAP)")
Kdf[c(1,2,5,3,6)] 

#Checking which assessors belong to which component
Kclasses <- assign_cluster(Kenranks, soft = FALSE, expand = FALSE)$map_cluster
MMClass <- cbind(Kclasses, fruittable)

K1 <- MMClass[which(MMClass$Kclasses == "Cluster 1"),] 
K2 <- MMClass[which(MMClass$Kclasses == "Cluster 2"),]

#Proportion of males in the overall data set compared to the different 
#components
(lwgenloc(fruittable,"O","M") + lwgenloc(fruittable,"E","M"))/nrow(fruittable)
(lwgenloc(K1,"O","M") + lwgenloc(K1,"E","M"))/nrow(K1)
(lwgenloc(K2,"O","M") + lwgenloc(K2,"E","M"))/nrow(K2)

#Same with the proportion of data set taken from Oxford
(lwgenloc(fruittable,"O","M") + lwgenloc(fruittable,"O","F"))/nrow(fruittable)
(lwgenloc(K1,"O","M") + lwgenloc(K1,"O","F"))/nrow(K1)
(lwgenloc(K2,"O","M") + lwgenloc(K2,"O","F"))/nrow(K2)

#Plotting the posterior probabilities of the Peaches consensus
plot(Kenranks, parameter = "rho", items = "Peaches")

#Checking which assessors are in Cluster 1 of both PL and Mallows
intersect(rownames(mg1[-1]),rownames(K1[-1]))





#CAYLEY DISTANCE
fruitMalcaymix <- compute_mallows_mixtures(n_clusters = c(1:16), mrank, 
                                           metric = "cayley", nmc = 10^5L, 
                                           include_wcd=TRUE)

#Elbow plot suggests non-mixture scenario
plot_elbow(fruitMalcaymix, burnin = 500)

#Creating a non-mixture model
Cayranks <- compute_mallows(mrank, metric = "cayley", n_clusters = 1L, 
                            save_clus = TRUE, save_aug = TRUE, nmc = 10^5)
Cayranks$burnin <- 5000

#Computing statistics of and plotting the posterior density of the shape
#parameter
compute_posterior_intervals(Cayranks, parameter = "alpha")
plot(Cayranks)


plot(Cayranks, parameter = "rho", items = "Apples")
plot(Cayranks, parameter = "rho", items = "Raspberries")

# Compute the CP consensus for each component
CayCP <- as.data.frame(compute_consensus(Cayranks, type = "CP"))
CayCP$cumprob <- NULL

# Compute the MAP consensus for each component
CayMAP <- compute_consensus(Cayranks, type = "MAP")
CayMAP$probability <- NULL

#Combine the two to compare
Cdf <- data.frame(CayCP,CayMAP)
colnames(Cdf) <- c("Ranking","CP","MAP","")
Cdf[c(1,2,3)] 






#HAMMING DISTANCE
fruitMalhammix <- compute_mallows_mixtures(n_clusters = c(1:16), mrank, 
                                           metric = "hamming", nmc = 10^5L, 
                                           include_wcd=TRUE)

#Elbow plot again suggests non-mixture model
plot_elbow(fruitMalhammix, burnin = 500)

#Creating a non-mixture model
Hamranks <- compute_mallows(mrank, metric = "hamming", n_clusters = 1L, 
                            save_clus = TRUE, save_aug = TRUE, nmc = 10^5)
Hamranks$burnin <- 5000

#Computing statistics of and plotting the posterior density of the shape
#parameter
compute_posterior_intervals(Hamranks, parameter = "alpha")
plot(Hamranks)

plot(Hamranks, parameter = "rho", items = "Raspberries")
plot(Hamranks, parameter = "rho", items = "Pineapples")

# Compute the CP consensus for each component
HamCP <- as.data.frame(compute_consensus(Hamranks, type = "CP"))
HamCP$cumprob <- NULL

# Compute the MAP consensus for each component
HamMAP <- compute_consensus(Hamranks, type = "MAP")
HamMAP$probability <- NULL

#Combine the two to compare
Hdf <- data.frame(HamCP,HamMAP)
colnames(Hdf) <- c("Ranking","CP","MAP","")
Hdf[c(1,2,3)] 
