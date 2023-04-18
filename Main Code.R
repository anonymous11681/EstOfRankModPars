#Loading the relevant libraries
library(MASS)
library(ggplot2)
library(PlackettLuce)
library(PLMIX)
library(BayesMallows)



#####
#DATA ANALYSIS ON PERSONALLY COLLECTED FRUIT DATA
#Removing names from the data
fruittable <- as.data.frame(Ranking_Table[-c(1)])


#Basic statistics; number of people asked
nrow(fruittable)

#Plot displaying the proportions of Genders from either location
genloc <- ggplot(fruittable, aes(x = Location, fill = Gender)) + 
  ylab("Frequency")
genloc + geom_bar(position = "fill") + scale_x_discrete(name ="Location", 
                                                        labels=c("Exeter","Oxford"))

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

#Rank totals across the data
colSums(fruits)

#Producing a plot of the rank totals
sums <- as.data.frame(colSums(fruits))
sums$key <- colnames(fruits)
sums$key <- factor(sums$key, levels=unique(sums$key))
ggplot(sums, aes(x=key, y=colSums(fruits), fill = key)) + 
  geom_bar(stat = "identity") + theme(axis.text.x = element_text(angle = 90)) + 
  xlab("Fruits") + ylab("Rank Totals") + theme(legend.position = "none")

#Creating an explicit ordering of each assessors rankings
mrank <- apply(as.matrix.noquote(fruits),2,as.numeric)
ranks <- as.rankings(mrank)
ranks

#####
#Fitting a Plackett-Luce Model using MLE via the PlackettLuce package
#prior <- list(mu = rep(0, ncol(ranks)), Sigma = diag(rep(9, ncol(ranks))))
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
#be very low
PL_prob(c(8,7,6,2,3,1,4,10,9,5), coef(fruitPL, log = FALSE))
PL_prob(c(3,4,5,9,8,10,7,1,2,6), coef(fruitPL, log = FALSE))
PL_prob(c(1,2,3,4,5,6,7,8,9,10), coef(fruitPL, log = FALSE))

#Probabilities of seeing the rankings seen in the fruits dataset
probdata <- c()
for (x in 1:nrow(mrank)) {
  probdata <- append(probdata, PL_prob(mrank[x,], coef(fruitPL, log = FALSE)))
}

#Which such data points had the most and least likely rankings according to
#the data?
which.max(probdata)
max(probdata)
which.min(probdata)
min(probdata)


#Using the PLMIX package to get object parameters via MAP estimation and Gibbs
#Sampling

#Ordering the items according to their ranks
RR <- as.top_ordering(ranks)
asnum <- nrow(RR)

#Function that returns a set of MAP estimates given noninformative priors
MAPFlat <- function(ordering, nmix, omega = c(1), sumry = TRUE, plot = TRUE) {
  if (sum(omega) != 1 | length(omega) != nmix) {
    return("Please choose a valid set of weights")
  } else {
    MPL <- mapPLMIX(ordering, K=ncol(ordering), G=nmix, plot_objective = plot)
    if (sumry == TRUE) {
      print(summary(MPL))
    }
    return(MPL)
  }
}

#With a single mixture component (Non-mixture PL); corresponds to the same
#method using PlackettLuce package
MFlat <- MAPFlat(RR, 1, omega = 1, sumry = FALSE)
summary(MFlat)
bicPLMIX(MFlat$max_objective, RR, G=1)$bic
MFlat$P

#Function that returns a Mallows Mixture Model with mixnum cluster and Gamma prior
nPLMix <- function(ord, mixnum, nassess, shape, rate, alpha) {
  return(mapPLMIX(ord, K=ncol(ord), G=mixnum, 
                  init = list(p = NULL, omega = rep(1/mixnum, mixnum)),
                  hyper = list(shape0 = matrix(shape, nrow = mixnum, 
                                               ncol = nassess), 
                               rate0 = rep(rate, mixnum),
                               alpha0 = rep(alpha, mixnum)), 
                  plot_objective = TRUE))
}

#Non-Mixture PL with Gamma(3, 0.5) priors
MPL <- nPLMix(RR, 1, asnum, 3, 0.5, 1) 
summary(MPL)
MPL$P


#With 2 mixture components, initially equally weighted, with Gamma(3,0.5) priors
MPL2 <- nPLMix(RR, 2, asnum, 3, 0.5, 1) 
summary(MPL2)
MPL2$P
MPL2$class_map

#With 3 mixture components
MPL3 <- nPLMix(RR, 3, asnum, 3, 0.5, 1) 
MPL3$P
MPL3$class_map

#With 4 mixture components
MPL4 <- nPLMix(RR, 4, asnum, 3, 0.5, 1) 
MPL4$P
MPL4$class_map


# With 5 mixture components
MPL5<- nPLMix(RR, 5, asnum, 3, 0.5, 1) 
MPL5$P
MPL5$class_map


#With 6 mixture components
MPL6 <- nPLMix(RR, 6, asnum, 3, 0.5, 1) 
MPL6$P
MPL6$class_map

#With 6 mixture components
MPL7 <- nPLMix(RR, 7, asnum, 3, 0.5, 1) 
MPL7$P
MPL7$class_map

#Using Gibbs sampling with previously found MAP estimates with vague priors 
#to give model parameter estimates; if no MAP estimates are supplied, assume
#randomly uniform component membership
GibbsEst <- function(ord, mixnum, MAPmod = NULL){
  if (is.null(MAPmod) == TRUE) {
    return(gibbsPLMIX(ord, K=ncol(ord), G=mixnum, init = list(p = NULL, z=NULL)))
  } else {
    return(gibbsPLMIX(ord, K=ncol(ord), G=mixnum, init = list(p = MAPmod$P, 
                           z=binary_group_ind(MAPmod$class_map, G=mixnum))))
  }
}

#Getting the respective Gibbs Sample estimates using MAP estimates and uniform 
#prior membership

GPL <- GibbsEst(RR, 1, MPL)
GPLU <- GibbsEst(RR, 1)
matrix(colMeans(GPL$P),ncol=ncol(RR))

#####
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

GPL7<- GibbsEst(RR, 7, MPL7)
GPL7U <- GibbsEst(RR, 7)
matrix(colMeans(GPL7$P),ncol=ncol(RR))

#####
#Generating various values concerning model selection criteria, comparing the
#6 models produced above
modsel <- selectPLMIX(pi_inv=RR, seq_G=1:7, parallel=TRUE, 
            MAPestP=list(MPL$P, MPL2$P, MPL3$P, MPL4$P, MPL5$P, MPL6$P, MPL7$P), 
            MAPestW=list(MPL$W, MPL2$W, MPL3$W, MPL4$W, MPL5$P, MPL6$W, MPL7$P),
            deviance=list(GPL$deviance, GPL2$deviance, 
                          GPL3$deviance, GPL4$deviance, 
                          GPL5$deviance, GPL6$deviance, GPL7$deviance))

#Doing the same for where initial component membership was uniformly random
modselU <- selectPLMIX(pi_inv=RR, seq_G=1:7, parallel=TRUE, 
                       MAPestP=NULL, MAPestW=NULL,  
                       MCMCsampleP =list(GPLU$P, GPL2U$P, GPL3U$P, 
                                    GPL4U$P, GPL5U$P, GPL6U$P), 
                       MCMCsampleW =list(GPLU$W, GPL2U$W, GPL3U$W, 
                                         GPL4U$W, GPL5U$W, GPL6U$W), 
                       post_est = "mean",
                       deviance=list(GPLU$deviance, GPL2U$deviance, 
                                     GPL3U$deviance, GPL4U$deviance, 
                                     GPL5U$deviance, GPL6U$deviance))

#Model with 2 components seems to fit the most appropriately, though some 
#criterion heavily suggest 5 components (e.g. DIC1)
modsel$criteria
modselU$criteria

postpred <- ppcheckPLMIX(pi_inv = RR, seq_G=1:7, parallel=TRUE, 
                         MCMCsampleP =list(GPLU$P, GPL2U$P, GPL3U$P, 
                                           GPL4U$P, GPL5U$P, GPL6U$P, GPL7U$P), 
                         MCMCsampleW =list(GPLU$W, GPL2U$W, GPL3U$W, 
                                           GPL4U$W, GPL5U$W, GPL6U$W, GPL7U$W))
postpred$post_pred_pvalue


classes <- as.matrix(MPL5$class_map)
rwclass <- cbind(classes, fruittable)

mg1 <- rwclass[which(rwclass$classes == 1),] 
mg2 <- rwclass[which(rwclass$classes == 2),]
mg3 <- rwclass[which(rwclass$classes == 3),]
mg4 <- rwclass[which(rwclass$classes == 4),]
mg5 <- rwclass[which(rwclass$classes == 5),]

#Proportion of males in the overall dataset compared to the different 
#components
(lwgenloc(fruittable,"O","M") + lwgenloc(fruittable,"E","M"))/nrow(fruittable)
(lwgenloc(mg1,"O","M") + lwgenloc(mg1,"E","M"))/nrow(mg1)
(lwgenloc(mg2,"O","M") + lwgenloc(mg2,"E","M"))/nrow(mg2)
(lwgenloc(mg3,"O","M") + lwgenloc(mg3,"E","M"))/nrow(mg3)
(lwgenloc(mg4,"O","M") + lwgenloc(mg4,"E","M"))/nrow(mg4)
(lwgenloc(mg5,"O","M") + lwgenloc(mg5,"E","M"))/nrow(mg5)

#Same with the proportion of dataset taken from Oxford
(lwgenloc(fruittable,"O","M") + lwgenloc(fruittable,"O","F"))/nrow(fruittable)
(lwgenloc(mg1,"O","M") + lwgenloc(mg1,"O","F"))/nrow(mg1)
(lwgenloc(mg2,"O","M") + lwgenloc(mg2,"O","F"))/nrow(mg2)
(lwgenloc(mg3,"O","M") + lwgenloc(mg3,"O","F"))/nrow(mg3)
(lwgenloc(mg4,"O","M") + lwgenloc(mg4,"O","F"))/nrow(mg4)
(lwgenloc(mg5,"O","M") + lwgenloc(mg5,"O","F"))/nrow(mg5)



#####
#Fitting a Mallows Model to the data with Kendall Distance
fruitMalken <- compute_mallows(mrank, metric = "kendall", nmc = 5000L)
assess_convergence(fruitMalken)

fruitMalkenmix <- compute_mallows_mixtures(n_clusters = c(1:16), mrank, 
                                           metric = "kendall", nmc = 5000L, 
                                           include_wcd=TRUE)
plot_elbow(fruitMalkenmix, burnin = 500)


compute_posterior_intervals(fruitMalken, parameter = "alpha", burnin = 0)
plot(fruitMalken, parameter = "alpha", burnin = 0)

compute_posterior_intervals(fruitMalken, parameter = "rho", burnin = 0)
compute_consensus(fruitMalken, type = "CP", burnin = 0)
compute_consensus(fruitMalken, type = "MAP", burnin = 0)

#Doing the same with Cayley Distance
fruitMalcay <- compute_mallows(mrank, metric = "cayley", nmc = 10000L)
assess_convergence(fruitMalcay)

fruitMalcaymix <- compute_mallows_mixtures(n_clusters = c(1:16), mrank, 
                                           metric = "cayley", nmc = 5000L, 
                                           include_wcd=TRUE)
plot_elbow(fruitMalcaymix, burnin = 500)

compute_posterior_intervals(fruitMalcay, parameter = "alpha", burnin = 0)
plot(fruitMalcay, parameter = "alpha", burnin = 0)

compute_posterior_intervals(fruitMalcay, parameter = "rho", burnin = 0)
compute_consensus(fruitMalcay, type = "CP", burnin = 0)
compute_consensus(fruitMalcay, type = "MAP", burnin = 0)


#And again with Hamming Distance
fruitMalham <- compute_mallows(mrank, metric = "hamming", nmc = 10000L)
assess_convergence(fruitMalham)

fruitMalhammix <- compute_mallows_mixtures(n_clusters = c(1:16), mrank, 
                                           metric = "hamming", nmc = 5000L, 
                                           include_wcd=TRUE)
plot_elbow(fruitMalhammix, burnin = 500)

compute_posterior_intervals(fruitMalham, parameter = "alpha", burnin = 0)
plot(fruitMalham, parameter = "alpha", burnin = 0)

compute_posterior_intervals(fruitMalham, parameter = "rho", burnin = 0)
compute_consensus(fruitMalham, type = "CP", burnin = 0)
compute_consensus(fruitMalham, type = "MAP", burnin = 0)


