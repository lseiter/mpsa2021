library(igraph)
library(RSQLite)
library(jsonlite)
library(R6)

VotingSimulation <- 
  R6Class ("VotingSimulation",
           public = list(
             #setup
             configuration="tmp",
             DISTANCENEIGHBORS=10,  # average degree =  2*distance for small world graphs
             PROBREWIRE=0.0,  #0 for BA graph
             graphType="BA_min2", 
             N_VECTOR = c(1000),   #size of graph
             ITERATIONS=100,   
             
             adjMatrix = NULL,
             invDiagDegMatrix = NULL,
             graph=NULL,
             degreeVector = NULL,
             MAXSTEP=100,
             seeds = c( 0.18, 0.20, 0.22),
             tbs = c( 0.6, 0.8, 1.0),
             SEEDPERCENT_VECTOR=NULL,
             TBPERCENT_VECTOR =NULL,
             PERCENTDECAY_VECTOR = c(0.0),
             
             K=3,  #num candidates
             
             #simulation state
             n=1,
             step = 1,
             percentDecay = 0,
             winner = 0,
             currentFreq = NULL,
             prevPrefVector = NULL,
             
             #simulation data structures
             thresholdMatrix = NULL,
             thresholdMatrixTN = NULL,
             preferenceVectorT1 = NULL,
             preferenceVector = NULL,
             preferenceMatrixT1 = NULL,
             preferenceMatrix = NULL,
             numNeighborsMatrix = NULL,
             seedIndexVector = NULL,
             seedCountVector = NULL,
             
             #shock
             timeShock=3,
             lengthShock=NULL,
             sizeShock=NULL,
             percentShocked=1.0,
             preShockMatrix = NULL,
             SHOCK_VECTOR = c( -0.20, -0.10, 0.0, 0.10, 0.20),
             SHOCKLENGTH_VECTOR = c(3,4),
             restored = TRUE,
             
             #database
             db = NULL,
             
             
             initialize = function() {
               self$db <- dbConnect(SQLite(), dbname="tmp8.sqlite")
             },
             
             
             newGraph = function() {
               #self$graph<- sample_smallworld(dim=1, size=self$n, nei=self$DISTANCENEIGHBORS, p=self$PROBREWIRE)
               self$graph<- sample_pa( n=self$n, m=self$DISTANCENEIGHBORS,directed=FALSE)
               
               #ensure a minimum degree 2
               degreeVector <- degree(self$graph)
               for (i in 1:self$n)  {
                 if (degreeVector[i]==0) {
                   self$graph<-add_edges(self$graph,c(i, sample(1:vcount(self$graph),1))) #add random edge
                   self$graph<-add_edges(self$graph,c(i, sample(1:vcount(self$graph),1))) #add random edge
                 }
                 else if (degreeVector[i]==1) {
                   self$graph<-add_edges(self$graph,c(i, sample(1:vcount(self$graph),1))) #add random edge
                 }
                 
               }
               self$degreeVector = degree(self$graph)  #recompute
               self$adjMatrix= t(as.matrix(as_adjacency_matrix(self$graph)))  #transpose
               self$invDiagDegMatrix=diag(sapply(self$degreeVector, function(i) (1/i)),self$n,self$n)
             },
             
             newPreferenceVectorT1 = function() {
               self$preferenceVectorT1 = rep((self$K+1), self$n)  #initialize to all undecided value K+1
               offset=0
               for (candidate in 1:self$K) {
                 seedsThisCandidate <- self$seedCountVector[candidate]
                 for (i in 1:seedsThisCandidate) {
                   #seedIndexVector is random ordering of vertices
                   r <- self$seedIndexVector[offset+i]  #get next voter row using random ordering of vertices
                   self$preferenceVectorT1[r]=candidate    
                 }
                 offset <- offset + seedsThisCandidate
               }
             },
             
             newPreferenceMatrixT1 = function() {
               self$preferenceMatrixT1 = matrix( 0, nrow = self$n, ncol = self$K)  #all 0, no Preference
               for (r in (1:self$n)) {
                 if (self$preferenceVectorT1[r]<= self$K) {   #skip if undecided k+1
                   self$preferenceMatrixT1[r,self$preferenceVectorT1[r]] = 1   #seed
                 }
               }
             },
             
             newThresholdMatrixTN = function() {
               self$thresholdMatrixTN = matrix( -1, nrow = self$n, ncol = self$K )  #init all to -1
               offset=0
               for (candidate in 1:self$K) {
                 nSeeds = self$seedCountVector[candidate]
                 nTrueBelievers = ceiling(nSeeds * self$TBPERCENT_VECTOR[candidate])
                 if (nSeeds > 0)  {     #ugh, loop iterates once when  nSeeds 0!
                   for (i in 1:nSeeds) {
                     r = self$seedIndexVector[offset+i]
                     if (i<=nTrueBelievers) {
                       #true believer, threshold 0 for seed candidate,  threshold above 100% for other candidates
                       for (c in 1:self$K) {
                         self$thresholdMatrixTN[r,c] = ifelse(c==candidate,0,self$degreeVector[r]+1) #Nov 15
                       }
                     }
                     else {
                       #adherent,  seed candidate threshold=min random/degree, other candidates unique random/degree
                       #UPDATE - pick k random numbers between 2..degree.  possible duplicate values if degree-1<k
                       if (self$degreeVector[r] < self$K) {
                         randomDegreeVector = rep(2,self$K)
                       }
                       else {
                         randomDegreeVector = sample(2:self$degreeVector[r], self$K, replace= (self$degreeVector[r]-1<self$K ))
                       }
                       index = which.min(randomDegreeVector)
                       #seed candidate should get min threshold, so swap with min
                       tmp=randomDegreeVector[index]
                       randomDegreeVector[index] = randomDegreeVector[candidate]
                       randomDegreeVector[candidate] = tmp
                       
                       for (c in 1:self$K) {
                         self$thresholdMatrixTN[r,c] = randomDegreeVector[c]   #Nov 15
                      }
                     }
                   }
                 }
                 offset = offset + nSeeds
               }
               #for non-seed voters,  threshold = random(2..degree)
               for (r in 1:self$n) {
                 if (self$degreeVector[r] < self$K) {
                   randomDegreeVector = rep(2,self$K)
                 }
                 else {
                   randomDegreeVector = sample(2:self$degreeVector[r], self$K, replace= TRUE )
                 }
                 if (self$thresholdMatrixTN[r,1] == -1)  {
                   for (c in 1:self$K) {
                     self$thresholdMatrixTN[r,c]  = randomDegreeVector[c]  #Nov 15
                   }
                 }
               }
               
             },
             
             newPreferenceMatrixTn = function() {
               self$preferenceMatrix[]=0    #reset values of existing matrix 
               self$preferenceVector[]=(self$K+1)  #init to all undecided
               #assign based on neighbor preferences and node threshold
               for (r in (1:self$n))   {
                 min<-self$degreeVector[r]+1  #initialize to max threshold exceeding degree
                 for (c in (1:self$K))  {
                   if ( self$numNeighborsMatrix[r,c]>=self$thresholdMatrix[r,c] )  {
                       #threshold met, test for new min threshold
                     if (self$thresholdMatrix[r,c]<min)  {
                       min<-self$thresholdMatrix[r,c]
                     }
                   }
                 }
                 if (min<self$degreeVector[r]+1)  {
                   #could be several candidates, select random
                   min_indices <- c()
                   neighb <-c()
                   for (c in (1:self$K)) {
                     if (self$thresholdMatrix[r,c] == min & self$numNeighborsMatrix[r,c] >= min) {
                       min_indices <- c(min_indices,c)
                       neighb <- c(neighb,self$numNeighborsMatrix[r,c])
                     }
                   }
                   #sample - if only one item in vector sample generates random number from 1..item so need to test length
                   if (length(min_indices) <= 1) {
                     preference<-min_indices[1]
                   } 
                   else {
                     
                     #version 7, pick candidate with highest neighbor preference
                     #max_neighbor = which.max(neighb)  #index of max
                     #preference = min_indices[max_neighbor]
                     
                     
                     #version 6, bandwagon, pick lowest candidate affected by shock
                     if (self$sizeShock == 0) {
                       preference<-sample(min_indices,1) #random choice when no shock
                     }
                     else {
                     # version 8, if positive shock, assign preference to candidate 1 
                       if (self$sizeShock>0) {
                        if (1 %in% min_indices) {
                         preference = 1  
                        }
                        else {
                          preference<-3 #positive shock would impact candidate 2 negatively, so pick 3
                        }
                       }
                       else {
                         #negative shock, assign preference to next candidate
                         if (1 %in% min_indices) {
                           preference<-min_indices[2]
                         }
                         else {
                          preference<-2 #bandwagon, negative shock would impact candidate 2 positively,, so pick 2
                         }
                       }
                     }
                  
                     
                   }
                   self$preferenceMatrix[r,preference] = 1
                   self$preferenceVector[r]=preference  
                 }
               }
             },
             
             addShock = function() {
               thresholdAdjustSign = -1*sign(self$sizeShock)  #adjust in opposite direction of shock
               
               for (r in 1:self$n) {
                 #positive shock should reduce threshold for candiate 1, negative shock should increase threshold
                 #don't shock true believers (candidate 1, 2, or 3). 
                 
                 if ( self$thresholdMatrix[r,1] > 0 & self$thresholdMatrix[r,2] > 0 & self$thresholdMatrix[r,3] > 0)
                 { 
                   degreeShock = ceiling(self$degreeVector[r] * abs(self$sizeShock))
                   
                   #self$thresholdMatrix[r,1] =  self$thresholdMatrix[r,1] + (ceiling(self$thresholdMatrix[r,1] * abs(self$sizeShock)) * thresholdAdjustSign)
                   self$thresholdMatrix[r,1] =  self$thresholdMatrix[r,1] + (degreeShock * thresholdAdjustSign)
                   #maintain minimum of 2 and max of degree
                   if (self$thresholdMatrix[r,1] < 2)
                   {
                     self$thresholdMatrix[r,1] = 2
                
                   }
                   if (self$thresholdMatrix[r,1] > self$degreeVector[r]) 
                   {
                     self$thresholdMatrix[r,1] = self$degreeVector[r]
                   
                   }
                   
                   #self$thresholdMatrix[r,2] =  self$thresholdMatrix[r,2] - (ceiling(self$thresholdMatrix[r,2] * abs(self$sizeShock)) * thresholdAdjustSign)
                  
                   #UPDATE TMP8, DON"T ADJUST CANDIDATE 2
                   #self$thresholdMatrix[r,2] =  self$thresholdMatrix[r,2] - (degreeShock * thresholdAdjustSign)
                   #if (self$thresholdMatrix[r,2]< 2)
                   #{
                   #  self$thresholdMatrix[r,2] = 2
                   #  #cat('adjust to min 2, ')
                   #}
                   
                   #if (self$thresholdMatrix[r,2] > self$degreeVector[r]) 
                   #{
                   #  self$thresholdMatrix[r,2] = self$degreeVector[r]
                     #cat('adjust to max degree, ')
                   #}
                   #cat ("after thresholds", self$thresholdMatrix[r,],"\n")
                 }
                 #else {
                 #  cat("tb, do not shock", self$thresholdMatrix[r,],"\n")
                 #}
               }
             },
             
             preserveThresholds = function() {
               self$restored = FALSE
               y=as.vector(self$thresholdMatrix)
               self$preShockMatrix = matrix(y,nrow=self$n,ncol = self$K)
               
             },
             
             restoreThresholds = function() {
               self$restored = TRUE
               y=as.vector(self$preShockMatrix)
               self$thresholdMatrix =  matrix(y,nrow=self$n,ncol = self$K)
               
             },
             
            
             
             #weight edges for decay
             randomWeightVector = function ()  {
               weightVector <- rep(1, self$n)  #all voters have default weight 1
               m=ceiling(self$percentDecay * self$n)
               randomSample <- sample(1:self$n, m, replace=FALSE )  #select m voters
               for (i in 1:m) {
                 r <- randomSample[i]  #voter row
                 weightVector[r]<-runif(1, 0.5, 0.9)    #set weight to random value between 0.5 and 0.9
               }
               weightVector
             },
             
             updatenumNeighborsMatrix = function() {
               if (self$percentDecay > 0.0)  {
                 weightVector<- self$randomWeightVector()
                 weightDiagMatrix<- diag(weightVector) 
                 weightedPreferenceMatrix<- weightDiagMatrix  %*% self$preferenceMatrix
                 self$numNeighborsMatrix = self$invDiagDegMatrix %*% self$adjMatrix %*% weightedPreferenceMatrix
               }
               else {
                 self$numNeighborsMatrix = self$adjMatrix %*% self$preferenceMatrix
               }
             },
             
             freqTable = function(cVector) {
               kplus1=self$K + 1
               n=self$n
               tmp=table(factor(cVector,levels=c(1:kplus1)))
               #round to 4 digits, proportion of n
               sapply(tmp,function(i) round(1.0*i/n, 4))
             },
             
             checkWinner = function() { 
               self$currentFreq<- self$freqTable(self$preferenceVector)
               self$winner=as.vector(which.max(self$currentFreq))
             },
             
             saveResult = function()  {
               result=ifelse(self$step ==self$MAXSTEP,"maxsteps","converge")
               ranking<-order(- self$currentFreq)
               winnerPercent<- self$currentFreq[ranking[1]]
               winnerLead <- winnerPercent - self$currentFreq[ranking[2]]
               
               seedJSON <-toJSON(self$SEEDPERCENT_VECTOR)
               tbJSON <- toJSON(self$TBPERCENT_VECTOR)
               degreeDistributionJSON <- toJSON(degree_distribution(self$graph))
               insertStmt <- sprintf("INSERT INTO session (configuration,iterations,k,type,DISTANCENEIGHBORS,probRewire,n,step,result,finalPercentages,seed,tb,decay,timeShock,lengthShock,sizeShock,percentShocked,winner,winnerPercent,winnerLead,meanDistance,clusterCoeff,meanDegree, degreeDistribution) 
                                     VALUES ('%s',%d,%d,'%s',%d,%f,%d,%d,'%s','%s','%s','%s',%f,%d,%d,%f,%f,%d,%f,%f,%f,%f, %f,'%s')",
                                     self$configuration, self$ITERATIONS, self$K,self$graphType,self$DISTANCENEIGHBORS,self$PROBREWIRE, 
                                     self$n,self$step,result,toString(self$currentFreq,sep=','), seedJSON, tbJSON, self$percentDecay, 
                                     self$timeShock, self$lengthShock, self$sizeShock, self$percentShocked,
                                     self$winner,winnerPercent,winnerLead,
                                     mean_distance(self$graph),transitivity(self$graph), mean(degree(self$graph)),degreeDistributionJSON)
               dbSendQuery(conn = self$db,insertStmt)
             },
             
             run = function()  {
               for (n in self$N_VECTOR) {
                 self$n=n   #used to be 100, 500, 1000.  Now just 1000
                 for (i in 1:self$ITERATIONS)  {
                   
                   self$newGraph()
                   for (s in self$seeds)  {
                     self$SEEDPERCENT_VECTOR <- c(s, 0.15, 0.15)  #candidate 2 and 3 both 10% UPDATED 3/31/21 to 15%
                     
                     #dependent on n and seed % vector.  compute outside tb loop so all tb use same seed assignment
                     self$seedCountVector=sapply(self$SEEDPERCENT_VECTOR, function(i) (i*self$n))
                     
                     self$seedIndexVector = sample(1:self$n, sum(self$seedCountVector), replace=FALSE )
                     # compute so same vector used for each tb and decay
                     self$newPreferenceVectorT1()   #based on seed vectors 
                     self$newPreferenceMatrixT1()   #based on preferenceVectorT1
                     
                     for (t in self$tbs)  {
                       self$TBPERCENT_VECTOR <- c(t, 1.0, 1.0)
                       cat("seed " , self$SEEDPERCENT_VECTOR, " tb ",self$TBPERCENT_VECTOR, n,i,"\n")
                       self$newThresholdMatrixTN()  #based on seed and tb.  create thresholdMatrixTN
                       for (percentDecay in self$PERCENTDECAY_VECTOR)  {
                         #cat("decay",percentDecay,"\n")
                         for (sizeShock in self$SHOCK_VECTOR) {
                           #cat("sizeShock",sizeShock,"\n")
                           for (lengthShock in self$SHOCKLENGTH_VECTOR) {
                             #cat("lengthShock",lengthShock,"\n")
                             self$lengthShock = lengthShock
                             self$sizeShock = sizeShock
                             self$percentDecay=percentDecay
                             self$step=1
                             
                             #reset so each decay/shock starts with same initial Preferences
                             self$preferenceVector = sapply(self$preferenceVectorT1 , function(z) z)
                             
                             y=as.vector(self$preferenceMatrixT1)
                             self$preferenceMatrix = matrix(y,nrow=self$n,ncol = self$K)  
                             
                             self$updatenumNeighborsMatrix()  #updates numNeighborsMatrix based on preferenceMatrix
                             
                             #reset thresholdmatrix to thresholdMatrixTN so each shock starts with same configuration
                             y=as.vector(self$thresholdMatrixTN)
                             self$thresholdMatrix = matrix(y,nrow=self$n,ncol = self$K)
                             
                             #compute next iteration of simulation #4/2 just loop 10 times to avoid history
                             #while (!self$votingComplete && self$step<self$MAXSTEP) {  
                             for (numsteps in 1:10) {
                               self$step=self$step+1
                               if (self$sizeShock != 0.0) {
                                 if (self$step == self$timeShock){
                                   self$preserveThresholds()
                                   self$addShock()
                                 }
                                 else if (self$step == self$timeShock + self$lengthShock) {
                                   self$restoreThresholds()
                                 }
                               }
                               self$newPreferenceMatrixTn()  #updates preferenceVecor and preferenceMatrix based on numNeighborsMatrix and thresholdMatrix
                               self$updatenumNeighborsMatrix()  #updates numNeighborsMatrix based on preferenceMatrix
                             } 
                             
                             #current configuration done (looped 10 times)
                             self$checkWinner() #compute based on preferenceVector
                             self$saveResult()  #store result in database
                           }
                         }
                         
                       }
                     }
                   }
                 }
               }
               
             }
           )     
           
  )

euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))


main <- function() {
  mySimulation <- VotingSimulation$new()
  mySimulation$run()
}

#run the simulation
main()
