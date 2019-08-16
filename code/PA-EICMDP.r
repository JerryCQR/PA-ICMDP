###PA-EICMDP algorithm###
###QiRong Cai,2019.7.14###

#The code uses OV as an example. It can be applied to other cancer type data by modifying the input data.
#Read mutation matrix--SNV and expression matrix--GE
SNV_data<-read.csv('F:/我的学习资料/协作通路github上传文件/data/OV/NEW_OV_SNV.csv')
GE_data<-read.csv('F:/我的学习资料/协作通路github上传文件/data/OV/NEW_OV_GE.csv')

rownames(SNV_data)<-SNV_data[,1]
SNV_data<-SNV_data[,-1]
rownames(GE_data)<-GE_data[,1]
GE_data<-GE_data[,-1]

#Data preprocessing, extracting common parts of SNV and GE
snv_colsum<-colSums(SNV_data) 
SNVdata<-SNV_data[,snv_colsum>4]
GEdata<-GE_data[,which(colnames(GE_data)%in%colnames(SNVdata))]

#Extract TP53 from A and B
TP53_snv <- as.matrix(SNVdata[,which(colnames(SNVdata)=='TP53')])
TP53_ge <- as.matrix(GEdata[,which(colnames(GEdata)=='TP53')])
SNVdata <- SNVdata[,-which(colnames(SNVdata)=='TP53')]
GEdata <- GEdata[,-which(colnames(GEdata)=='TP53')]

#Pvalue--Testing the correlation between each gene by t-test
p_value <- vector(length = ncol(GEdata))
every_Pvalue <- matrix(0,ncol(SNVdata),ncol(GEdata),dimnames = list(colnames(SNVdata),colnames(GEdata)))
for(i in seq(ncol(SNVdata))){
  n1<-GEdata[which(SNVdata[,i]==1),]
  n2<-GEdata[which(SNVdata[,i]==0),]
  
  for(j in seq(ncol(GEdata))){
    p_value[j]<- t.test(n1[,j],n2[,j])$p.value
  }
  every_Pvalue[i,] <- p_value
}
Pvalue <- every_Pvalue


#Calculate the Pearson coefficient between k gene
cor_sum <- function(E,x){
  k <- length(x)
  result <- 0
  if(k>1){
    for(i in 1:(k-1)){
      for(j in (i+1):k){
        temp <- cor(E[,x[i]],E[,x[j]])
        result <- result+abs(temp)
      }
    }
  }
  else{
    result <- 1
  }
  return(result)
}

# Objective function of two co-ocurrence submatrices
fitness_double_mat <- function(A,E,P,x2){
  M <- TP53_snv #first geneset M
  temp_N <- matrix(0,1,n)
  temp_N[x2] <- 1
  index_N <- which(temp_N==1)
  N <- as.matrix(A[,index_N]) #second geneset N
  Nweight <- t(as.matrix(P[,index_N]))
  m1 <- TP53_snv
  n1 <- as.matrix(rowSums(N))
  c <- length(intersect(which(m1>0),which(n1>0))) #c(M,N)
  un <- length(union(which(m1>0),which(n1>0)))#b(M,N)
  d <- un-c #b(M,N)
  w_M <- 0#w(M)
  w_N <- length(which(n1>1))#w(N)
  R_N <- cor_sum(E,x2) 
  
  H <- 10*((k-sum(Nweight))/k)*c-d- w_M - w_N+5*R_N
  L <- list(H,c,un)
  return(L)
}


#Roulette wheel selection
select <- function(pop){
  I <- order(pop[,k+1],decreasing = T)
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+3)
  pop_next[1,] <- pop[1,]
  if(sum(pop[,k+1])<0){
    p=pop[,k+1]/(-sum(pop[,k+1]))
    for(i in 2:popsize){
      random_data=(-runif(1))
      p_cumsum <- cumsum(p)
      temp <- which(p_cumsum<=random_data)
      index1 <- temp[1]
      pop_next[i,] <- pop[index1,]
    }
  }
  else{
    p=pop[,k+1]/sum(pop[,k+1])
    for(i in 2:popsize){
      random_data=runif(1)
      p_cumsum <- cumsum(p)
      temp <- which(p_cumsum>=random_data)
      index1 <- temp[1]
      pop_next[i,] <- pop[index1,]
    }
  }
  return(pop_next)
}

#A greedy based recombination operator is presented to generate new offsprings
mutation <- function(pop){
  chouqu <- sapply(1:popsize,function(x){
    sample(k,1)
  }) 
  for(j in 1:popsize){
    Kgene <- pop[j,1:k]
    f1<- pop[j,k+1]
    temp<- chouqu[j]
    ind_index<- pop[j,k+3]
    beixuan <- ind[ind_index,]
    beixuan <- beixuan[!beixuan%in%Kgene] 
    for(i in 1:(n-k)){
      Kgene[temp] <- beixuan[i]
      f2 <- fitness_double_mat(SNVdata,GEdata,MinP,Kgene)[[1]]
      if(f2>f1){
        pop[j,k+1] <- f2
        pop[j,1:k] <- Kgene
      }
    }
  }
  return(pop)
}

#significance test
significance_A_E <- function(SNVdata,GEdata,MinP,subset_M){
  m <- ncol(SNVdata)
  w <- matrix(0,1,1000)
  l <- length(subset_M)
  for(j in 1:1000){
    index <- round(runif(l,min = 1,max = m))
    w[j] <- fitness_double_mat(SNVdata,GEdata,MinP,index)[[1]]

  }
  p <- sum(w>=fitness_double_mat(SNVdata,GEdata,MinP,subset_M)[[1]])/1000
  return(p)
}



#Begin,Data initialization
n <- ncol(SNVdata)
geneName <- colnames(SNVdata)
popsize <- 60 #population size
iteration <- 200 # maximum evolution generation
k <-4  #Number of genes 
obj_value <- matrix(0,iteration,k+1)#Optimal value for each iteration
#1:k  records the collaborative genes,
#k+1  records the weight of the selected gene sets，
#k+2  records the co-occurrence significance level of the two gene sets
#k+3  record the location of the individual in the parent


pop <- matrix(0, popsize,k+3)#population
ind <- matrix(0,popsize,n)#individual

#Gene weight Wt(g)
MinP <- sapply(1:nrow(Pvalue),function(x){
  mean(sort(Pvalue[x,])[1:(n/5)])})
MinP <- matrix(MinP, nrow = 1, ncol = n,dimnames = list('1',rownames(Pvalue)))

#Initial
total <- 1:n
  for(i in 1:popsize){
    B <- sample(total,n) 
    ind[i,] <- B
    for(j in 1:floor(n/k)){
      fit <- fitness_double_mat(SNVdata,GEdata,MinP,B[1:k])[[1]]
      if(pop[i,k+1]<fit){
        pop[i,1:k] <- B[1:k]
        pop[i,k+1] <-fit 
      }
      B <- B[k+1:n]
    }
    pop[i,k+3] <- i 
  }

#iteration
R <- 0
i <- 1
while(i<=iteration & R<10){
  fit_vector <-as.numeric(pop[,k+1]) 
  temp_maxweight <- max(fit_vector)
  pop<- select(pop)
  pop <- mutation(pop)
  I <- order(pop[,k+1],decreasing = T)
  pop <- pop[I,]
  obj_value[i,] <- pop[1,1:(k+1)]
  i <- i+1
  maxweight <- max(pop[,k+1])
  if(maxweight==temp_maxweight)
    R=R+1
  else
    R=0
}
I <- order(pop[,k+1],decreasing = T)
pop <- pop[I,]

maxpop <- pop[!duplicated(pop[,1:k]),,drop=F]

# significance test
m <- nrow(maxpop)
for(i in 1:m){
  maxpop[i,k+2] <- significance_A_E(SNVdata,GEdata,MinP,maxpop[i,1:k])
}

maxpop[,1:k] <- apply(maxpop[,1:k,drop=F],2,function(x) {geneName[x]})
