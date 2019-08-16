###PA-ICMDP algorithm###
###QiRong Cai,2019.7.14###

#The code uses GBM as an example. It can be applied to other cancer type data by modifying the input data.
#Read mutation matrix--SNV and expression matrix--GE
SNV_data<-read.csv('F:/data/GBM/NEW_GBM_SNV.csv')
GE_data<-read.csv('F:/data/GBM/NEW_GBM_GE.csv')

rownames(SNV_data)<-SNV_data[,1]
SNV_data<-SNV_data[,-1]
rownames(GE_data)<-GE_data[,1]
GE_data<-GE_data[,-1]

#Data preprocessing, extracting common parts of SNV and GE
snv_colsum<-colSums(SNV_data) 
SNVdata<-SNV_data[,snv_colsum>3]
GEdata<-GE_data[,which(colnames(GE_data)%in%colnames(SNVdata))]

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

#Objective function of a submatrix
fitness_single_mat <- function(A,E,P,x){
  temp <- matrix(0,1,n)
  temp[x] <- 1
  y <- temp
  index <- which(y==1)
  a <- as.matrix(A[,index])
  b <- t(as.matrix(P[,index]))
  a_indexsum <- rowSums(a)
  a_colsum <- colSums(a)
  f <- (sum(a_indexsum>0)*2)-(sum(a_indexsum))
  coverage <- sum(a_indexsum>0)
  L <- list(f,coverage)
  return(L)
}

# Objective function of two co-ocurrence submatrices
fitness_double_mat <- function(A,E,P,x1,x2){
  temp_M <- matrix(0,1,n)
  temp_M[x1] <- 1
  index_M <- which(temp_M==1)
  M <- as.matrix(A[,index_M]) #first geneset M
  Mweight <- t(as.matrix(P[,index_M]))
  temp_N <- matrix(0,1,n)
  temp_N[x2] <- 1
  index_N <- which(temp_N==1)
  N <- as.matrix(A[,index_N]) #second geneset N
  Nweight <- t(as.matrix(P[,index_N]))
  
  m1 <- as.matrix(rowSums(M))
  n1 <- as.matrix(rowSums(N))
  c <- length(intersect(which(m1>0),which(n1>0))) #c(M,N)
  un <- length(union(which(m1>0),which(n1>0))) #b(M,N)
  d <- un-c #b(M,N)
  
  w_M <- length(which(m1>1))#w(M)
  w_N <- length(which(n1>1))#w(N)
  
  x_MN <- c(x1,x2)
  R_MN <- cor_sum(E,x_MN) 
  H <- alp*((k-sum(Mweight)-sum(Nweight))/k)*c-d- w_M - w_N+bet*R_MN
  L <- list(H,c,un)
  return(L)
}

#Roulette wheel selection
select <- function(pop){
  I <- order(pop[,k+1],decreasing = T)
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+9)
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
mutation <- function(pop,ite){
  if(ite<=0.7*iteration){
    chouqu <- sapply(1:popsize,function(x){
      sample(k,1)
    })
    chouqu_beixuan <- sapply(1:popsize,function(x){
      sample(n-k,1)
    })
    for(j in 1:popsize){
      Kgene <- pop[j,1:k]
      M_num <- pop[j,k+2]
      N_num <- pop[j,k+3]
      if(M_num+N_num<k){
        g <- sample(n,1)
        while(g %in% Kgene ){
          g <- sample(n,1)
        }
        pop[j,k] <- g
      }
      
      Kgene <- pop[j,1:k]
      temp<- chouqu[j]
      temp_beixuan<- chouqu_beixuan[j]
      ind_index<- pop[j,k+7]
      beixuan <- ind[ind_index,]
      beixuan <- beixuan[!beixuan%in%Kgene]
      Kgene[temp] <- beixuan[temp_beixuan]
      
      if(j==1){
        M <- Kgene[1:M_num]
        N <- Kgene[(M_num+1):k]
        temp <- fitness_double_mat(SNVdata,GEdata,MinP,M,N)[[1]]
        if(pop[j,k+1]<temp){
          pop[j,k+1] <- temp
          pop[j,1:k] <- Kgene
          pop[j,k+2] <- length(M) 
          pop[j,k+3] <- length(N) 
          pop[j,k+4] <- fitness_double_mat(SNVdata,GEdata,MinP,M,N)[[2]]
          pop[j,k+5] <- fitness_double_mat(SNVdata,GEdata,MinP,M,N)[[3]]
          pop[j,k+8] <- fitness_single_mat(SNVdata,GEdata,MinP,M)[[2]]
          pop[j,k+9] <- fitness_single_mat(SNVdata,GEdata,MinP,N)[[2]]
          
        }
      }
      else{
        pop[j,1:k] <- Kgene
        M <- pop[j,1:M_num]  #1:k1
        N <- pop[j,(M_num+1):k]# k1+1:k
        pop[j,k+1] <- fitness_double_mat(SNVdata,GEdata,MinP,M,N)[[1]]
        pop[j,k+2] <- length(M) 
        pop[j,k+3] <- length(N) 
        pop[j,k+4] <- fitness_double_mat(SNVdata,GEdata,MinP,M,N)[[2]]
        pop[j,k+5] <- fitness_double_mat(SNVdata,GEdata,MinP,M,N)[[3]]
        pop[j,k+8] <- fitness_single_mat(SNVdata,GEdata,MinP,M)[[2]]
        pop[j,k+9] <- fitness_single_mat(SNVdata,GEdata,MinP,N)[[2]]
      }
      
    }
  }
  
  else{
    chouqu <- sapply(1:popsize,function(x){
      sample(k,1)
    }) 
    for(j in 1:popsize){
      Kgene <- pop[j,1:k]
      H1<- pop[j,k+1]
      M_num <- pop[j,k+2]
      N_num <- pop[j,k+3]
      
      temp<- chouqu[j]
      ind_index<- pop[j,k+7]
      beixuan <- ind[ind_index,]
      beixuan <- beixuan[!beixuan%in%Kgene] 
      beixuan <- as.matrix(append(beixuan,Kgene[temp]))
      M <-Kgene[1:M_num]  
      N <-Kgene[(M_num+1):k] 
      
      if(temp>M_num){
        M1 <- M
        N1 <- N
        N_temp <- N1
        M_temp <- M1
        H_N <- H1
        H_M <- H1
        for(i in 1:(n-k+1)){
          N1[temp-M_num]=beixuan[i]
          H_tempN <- fitness_double_mat(SNVdata,GEdata,MinP,M1,N1)[[1]]
          if(H_tempN>=H_N){
            N_temp<- N1 
            H_N <- H_tempN 
          }
        }
        

        N_sub<- N1[-(temp-M_num)]
        M_add <- M
        for(i in 1:(n-k+1)){
          M_add[M_num+1] <- beixuan[i]
          H_tempM <- fitness_double_mat(SNVdata,GEdata,MinP,M_add,N_sub)[[1]]
          if(H_tempM>H_M){
            H_M<- H_tempM
            M_temp<- M_add
          }
        }

        if(H_N>H1 || H_M>H1){
          if(H_N>=H_M){
            pop[j,1:length(M1)] <- M1
            pop[j,(length(M1)+1):(length(M1)+length(N_temp))] <- N_temp
            pop[j,k+1] <- H_N 
            pop[j,k+2] <- length(M1) 
            pop[j,k+3] <- length(N_temp) 
            pop[j,k+4] <- fitness_double_mat(SNVdata,GEdata,MinP,M1,N_temp)[[2]]
            pop[j,k+5] <- fitness_double_mat(SNVdata,GEdata,MinP,M1,N_temp)[[3]]
            pop[j,k+8] <- fitness_single_mat(SNVdata,GEdata,MinP,M1)[[2]]
            pop[j,k+9] <- fitness_single_mat(SNVdata,GEdata,MinP,N_temp)[[2]]
            
          }
          else{
            pop[j,1:length(M_temp)] <- M_temp
            pop[j,(length(M_temp)+1):(length(M_temp)+length(N_sub))] <- N_sub 
            pop[j,k+1] <- H_M 
            pop[j,k+2] <- length(M_temp) 
            pop[j,k+3] <- length(N_sub) 
            pop[j,k+4] <- fitness_double_mat(SNVdata,GEdata,MinP,M_temp,N_sub)[[2]]
            pop[j,k+5] <- fitness_double_mat(SNVdata,GEdata,MinP,M_temp,N_sub)[[3]]
            pop[j,k+8] <- fitness_single_mat(SNVdata,GEdata,MinP,M_temp)[[2]]
            pop[j,k+9] <- fitness_single_mat(SNVdata,GEdata,MinP,N_sub)[[2]]
            
          }
        }
      }
      else{
        M1 <- M
        N1 <- N
        N_temp <- N1
        M_temp <- M1
        H_N <- H1
        H_M <- H1
        for(i in 1:(n-k+1)){
          M1[temp] <- beixuan[i]
          H_tempM <- fitness_double_mat(SNVdata,GEdata,MinP,M1,N1)[[1]]
          if(H_tempM>H_M){
            H_M<- H_tempM
            M_temp<- M1 
          }
        }
        
        M_sub<- M1[-temp]
        N_add <- N
        for(i in 1:(n-k+1)){
          N_add[N_num+1] <- beixuan[i]
          H_tempN <- fitness_double_mat(SNVdata,GEdata,MinP,M_sub,N_add)[[1]]
          if(H_tempN>H_N){
            H_N<- H_tempN
            N_temp<- N_add
          }
        }
        
        if(H_N>H1 || H_M>H1){
          if(H_M>H_N){
            pop[j,1:length(M_temp)] <- M_temp
            pop[j,(length(M_temp)+1):(length(M_temp)+length(N1))] <- N1
            pop[j,k+1] <- H_M 
            pop[j,k+2] <- length(M_temp) 
            pop[j,k+3] <- length(N1) 
            pop[j,k+4] <- fitness_double_mat(SNVdata,GEdata,MinP,M_temp,N1)[[2]]
            pop[j,k+5] <- fitness_double_mat(SNVdata,GEdata,MinP,M_temp,N1)[[3]]
            pop[j,k+8] <- fitness_single_mat(SNVdata,GEdata,MinP,M_temp)[[2]]
            pop[j,k+9] <- fitness_single_mat(SNVdata,GEdata,MinP,N1)[[2]]
          }
          else{
            pop[j,1:length(M_sub)] <- M_sub
            pop[j,(length(M_sub)+1):(length(M_sub)+length(N_temp))] <- N_temp 
            pop[j,k+1] <- H_N 
            pop[j,k+2] <- length(M_sub) 
            pop[j,k+3] <- length(N_temp) 
            pop[j,k+4] <- fitness_double_mat(SNVdata,GEdata,MinP,M_sub,N_temp)[[2]]
            pop[j,k+5] <- fitness_double_mat(SNVdata,GEdata,MinP,M_sub,N_temp)[[3]]
            pop[j,k+8] <- fitness_single_mat(SNVdata,GEdata,MinP,M_sub)[[2]]
            pop[j,k+9] <- fitness_single_mat(SNVdata,GEdata,MinP,N_temp)[[2]]
          }
        }
        
      }
    }
  }
  return(pop)
}


#significance test
significance_A_E <- function(SNVdata,GEdata,MinP,x1,x2){
  m <- ncol(SNVdata)
  w <- matrix(0,1,1000)
  l <- length(x1)+length(x2)
  for(j in 1:1000){
    index <- round(runif(l,min = 1,max = m))
    l1 <- index[1:length(x1)]
    l2 <- index[(length(x1)+1):l]
    w[j] <- fitness_double_mat(SNVdata,GEdata,MinP,l1,l2)[[1]]

  }
  p <- sum(w>=fitness_double_mat(SNVdata,GEdata,MinP,x1,x2)[[1]])/1000
  return(p)
}

#Begin,Data initialization
n <- ncol(SNVdata)
geneName <- colnames(SNVdata)
popsize <- 60 #population size
iteration <- 200 # maximum evolution generation
k <-4  #Number of genes 
alp <- 10 #alpha
bet <-5 #beta
obj_value <- matrix(0,iteration,k+1) #Optimal value for each iteration
#1:k1  records the first selected genes
#(k1+1):k  records the second selected genes
#k+1  records the weight of the selected gene sets£¬
#k+2  records the number of genes in the first gene set
#k+3  records the number of genes in the second gene set
#k+4  records the common coverage of the two gene sets
#k+5  records the union coverage of the two gene sets
#k+6  records the co-occurrence significance level of the two gene sets
#k+7  record the location of the individual in the parent
#k+8  records the Objective function of first gene set
#k+9  records the Objective function of second gene set

pop <- matrix(0, popsize,k+9)#population
ind <- matrix(0,popsize,n)#individual

#Gene weight Wt(g)
MinP <- sapply(1:nrow(Pvalue),function(x){
  mean(sort(Pvalue[x,])[1:(n/5)])})#Indicates that only the first n/5 genes are considered
MinP <- matrix(MinP, nrow = 1, ncol = n,dimnames = list('1',rownames(Pvalue)))

#Initial
total <- 1:n
for(i in 1:popsize){
  B <- sample(total,n)
  ind[i,] <- B
  fit <- matrix(0,floor(n/(floor(k/2))),(floor(k/2)+1))
  for(j in 1:floor(n/(floor(k/2)))){
    fit[j,1:(floor(k/2))] <- B[1:floor(k/2)]
    fit[j,(floor(k/2)+1)]<- fitness_single_mat(SNVdata,GEdata,MinP,B[1:floor(k/2)])[[2]]
    B <- B[floor(k/2)+1:n]
  }
  I <- order(fit[,(floor(k/2)+1)],decreasing = T)
  fit <- fit[I,]
  
  pop[i,1:floor(k/2)] <- fit[1,1:(floor(k/2))] #first geneset M
  pop[i,k+8] <- fit[1,(floor(k/2)+1)]
  pop[i,k+2] <- floor(k/2)# the number of gene in first geneset
  pop[i,(floor(k/2)+1):(floor(k/2)*2)] <- fit[2,1:(floor(k/2))]#second geneset N
  pop[i,k+9] <- fit[2,(floor(k/2)+1)]
  pop[i,k+3] <- floor(k/2)# the number of gene in second geneset
  
  pop[i,k+7] <- i 
  pop[i,k+1] <- fitness_double_mat(SNVdata,GEdata,MinP,fit[1,1:(floor(k/2))],fit[2,1:(floor(k/2))])[[1]]
  
  
  
}

#iteration
R <- 0
i <- 1
while(i<=iteration){
  fit_vector <-as.numeric(pop[,k+1]) 
  temp_maxweight <- max(fit_vector)
  pop<- select(pop)
  pop <- mutation(pop,i)
  I <- order(pop[,k+1],decreasing = T)
  pop <- pop[I,]
  obj_value[i,] <- pop[1,1:(k+1)]
  i <- i+1
  maxweight <- max(pop[,k+1])
  
  if((i>0.5*iteration)&(maxweight==temp_maxweight))
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
  maxpop[i,k+6] <- significance_A_E(SNVdata,GEdata,MinP,maxpop[i,1:maxpop[i,k+2]],maxpop[i,(maxpop[i,k+2]+1):k])
}
# pop[,1:k] <- apply(pop[,1:k], 2, function(x) {geneName[x]})
maxpop[,1:k] <- apply(maxpop[,1:k,drop=F],2,function(x) {geneName[x]})

#Êä³ö½á¹û
# write.csv(maxpop,'E:/R_CODE/cancerData/maxpop.csv')

