library(entropy)
library(varrank)
library(sfsmisc)
library(scales)


interaction_matrix <- function(size,  frac,  mu, e=1){
  
  S <- length(size)
  probs.M <- size/sum(size)              
  pred_sp <- which(size %in% sample(size, floor(S*frac), prob=probs.M))  
  m_j <- (1/4+2/3) # consumer scaling Rall et al 2012
  m_i <- 2/3       # resource scaling Rall et al 2012
  f <- function(a0){
    IntePP0 <- (matrix(rep(size,S),ncol=S, byrow=F)^m_i *  matrix(rep(size,S),ncol=S, byrow=T)^m_j )
    diag(IntePP0) <- NA
    IntePP <- scales::rescale(IntePP0)
    IntePP[lower.tri(IntePP)] <- -IntePP[lower.tri(IntePP)]
    diag(IntePP) <- 1
    
    InteC0 <- (matrix(rep(size,S),ncol=S, byrow=T)/matrix(rep(size,S),ncol=S, byrow=F))^(3/4)
    diag(InteC0) <- NA
    InteC <- scales::rescale(InteC0)
    diag(InteC) <- 1
    
    A0 <- InteC 
    for(i in pred_sp){
      A0[i, 1:(i-1) ] <- e * IntePP[i, 1:(i-1) ] + InteC[i, 1:(i-1) ]
      A0[1:(i-1), i ] <- e * IntePP[1:(i-1), i ] + InteC[1:(i-1), i ]
    }
    A <- a0 * A0
    diag(A) <- 1
    out <- ((sum(A) - S)/((S-1)*S) - mu)^2  
    
    return(out)
  }
  
  a0 <- optimize(f,interval = c(0,100),maximum = F,tol = 10^-14)$minimum
  IntePP0 <- (matrix(rep(size,S),ncol=S, byrow=F)^m_i *  matrix(rep(size,S),ncol=S, byrow=T)^m_j )
  diag(IntePP0) <- NA
  IntePP <- scales::rescale(IntePP0)
  IntePP[lower.tri(IntePP)] <- -IntePP[lower.tri(IntePP)]
  diag(IntePP) <- 1
  
  InteC0 <- (matrix(rep(size,S),ncol=S, byrow=T)/matrix(rep(size,S),ncol=S, byrow=F))^(3/4)
  diag(InteC0) <- NA
  InteC <- scales::rescale(InteC0)
  diag(InteC) <- 1
  
  A0 <- InteC
  for(i in pred_sp){
    A0[i, 1:(i-1) ] <- e * IntePP[i, 1:(i-1) ] + InteC[i, 1:(i-1) ]
    A0[1:(i-1), i ] <- e * IntePP[1:(i-1), i ] + InteC[1:(i-1), i ]
  }
  A <- a0*A0
  diag(A) <- 1
  
  return(list(A=A, pred=pred_sp))
  
}

check.integer <- function(x) {
  x == round(x)
}

mut.info <- function(x, y, discretization.method=NULL){
  x <- na.omit(x)
  y <- na.omit(y)
  
  if(all(x == 0) | all(y==0)) {
    MI.norm <- NA
  } else {
    if( all(check.integer(x)) & all(check.integer(y)) ){
      H1 <- entropy.plugin(table(x), unit="log2")
      H2 <- entropy.plugin(table(y), unit="log2")
      H12 <- entropy.plugin(table(x,y), unit="log2")
      
    } else {
      x_disc <- table(discretization(data.df = data.frame(x), discretization.method=discretization.method, freq = F))
      y_disc <- table(discretization(data.df = data.frame(y), discretization.method=discretization.method, freq = F))
      xy_disc <- table(discretization(data.df = data.frame(x,y), discretization.method=discretization.method, freq = F))
      H1 <- entropy.plugin(x_disc, unit="log2")
      H2 <- entropy.plugin(y_disc, unit="log2")
      H12 <- entropy.plugin(xy_disc, unit="log2")
      
    }
    MI <-  H1+H2-H12
    MI.max <-  min(H1, H2)
    MI.norm <-  MI / MI.max
  }
  
  return(MI.norm)
}

GL_func <- function(df0, time){
  
  df0 <- cbind(time, df0)
  df <- df0[,-1]
  
  t <- length(sort(time))-1
  ST <- G <- G.ab <- L <- L.ab <- P <- G.num <- L.num <- P.num <-  NULL
  for(q in 1:t){
    gl <- ifelse(df[1+q,]>0, 1, 0) - ifelse(df[q,]>0,1,0)
    if(all(is.na(gl))){
      ST[[q]] <- NA
      G[[q]] <- NA
      L[[q]] <- NA
      P[[q]] <- NA
      G.ab[[q]] <- NA
      L.ab[[q]] <- NA
      G.num[[q]] <- NA
      L.num[[q]] <- NA
      P.num[[q]] <- NA
    } else {
      gl_ab <- df[1+q,] - df[q,]
      gained <- length(which(gl==1))
      gained_ab <- sum(gl_ab[which(gl_ab>0)])
      lost <- length(which(gl==-1))
      lost_ab <- -sum(gl_ab[which(gl_ab<0)])
      persist <- length(which((ifelse(df[1+q,]>0,1,0) + ifelse(df[q,]>0,1,0) )==2 ))
      total <- length(which((ifelse(df[1+q,]>0,1,0) + ifelse(df[q,]>0,1,0) ) > 0))
      total_ab <- sum(df[1+q,] + df[q,]) 
      ST[[q]] <- (gained+lost) / total
      G[[q]] <- gained/total
      L[[q]] <- lost/total
      P[[q]] <- persist/total
      G.ab[[q]] <- gained_ab/total_ab
      L.ab[[q]] <- lost_ab/total_ab
      G.num[[q]] <- gained
      L.num[[q]] <- lost
      P.num[[q]] <- persist
    }
  }
  
  out <- data.frame(time1=time[-length(time)], time2=time[-1], turnover=unlist(ST), 
                    gains=unlist(G), losses=unlist(L), persist=unlist(P),
                    gains_ab=unlist(G.ab), losses_ab=unlist(L.ab),
                    gains_num = unlist(G.num), losses_num = unlist(L.num),
                    persist_num = unlist(P.num) )
  
  return(out)
}

evenness <- function(ab){ 
  N <-  ab 
  NN <- N[which(N>0)]
  p <- NN / sum( NN, na.rm = T)   
  S <- length(NN)
  ev <- - sum(p*log(p) ) / log(S)
  ent <- - sum(p*log(p) )
  return(list(J=ev, H=ent, S=length(NN) ))
}


shdist <- function(g1,g2){
  m1 <- wgtMatrix(g1@graph, transpose = FALSE)
  m1[m1 != 0] <- 1
  m1 <- m1[sort(rownames(m1)), sort(rownames(m1)) ]
  m2 <- wgtMatrix(g2@graph, transpose = FALSE)
  m2[m2 != 0] <- 1
  m2 <- m2[sort(rownames(m2)), sort(rownames(m2)) ]
  shd <- 0
  s1 <- m1 + t(m1)
  s2 <- m2 + t(m2)
  s1[s1 == 2] <- 1
  s2[s2 == 2] <- 1
  ds <- s1 - s2
  ind <- which(ds > 0)
  m1[ind] <- 0
  shd <- shd + length(ind)/2
  ind <- which(ds < 0)
  m1[ind] <- m2[ind]
  shd <- shd + length(ind)/2
  d <- abs(m1 - m2)
  shd + sum((d + t(d)) > 0)/2
}




