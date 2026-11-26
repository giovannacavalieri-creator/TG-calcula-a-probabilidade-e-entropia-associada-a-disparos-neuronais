### Esse código corresponde ao cálculo das probabilidades preditivas e entropias para o cenário 2 da monografia ####
## T = 12000 - Cenário 2  ##

# beta = -3 (phi = 0.04742587), w: equilíbrio sináptico homogêneo forte 

library(glmnet)

# Configurações de entrada 

T_final <- 12000 #tamanho da amostra 
n = 5 #numero de neuronios da rede
w = matrix(c(0, 0, 3, 3, 3,
             0, 0, 3, 3, 3,
             3, 3, 0, 3,-12,
             3, 3, 3, 0, -12,
             3, 3, -12,-12, 0), ncol=n, byrow = T)  #w_0 verdadeira

beta <- -3 #taxa basal
H <- 4 #horizonte de predição
B <- 3000 #simulações MC
R_repl <- 110 #réplicas

t <- T_final/0.8 #0.8 da fração faz burn-in: geramos uma amostra maior para depois tirar 20% dela


# Funções necessárias para gerar X ----------------------------------------

#Função que seleciona todas as linhas da matriz "m" menos a do "indice" indicado 
matriz_pesos <- function(m,indice){
  matriz_pesos = matrix(c(m[!(row(m) %in% c(indice))]), nrow=nrow(m)-1)
  return(matriz_pesos)
}

#Função que seleciona a linha "indice" da matriz "m" indicada
matriz_x <- function(m,indice){
  matriz_x = c(m[(row(m) %in% c(indice))])
  return(matriz_x)
}

#Função taxa de disparo logísstica
phi <- function(x) {
  plogis(x)
}


#Função verifica indice do ultimo disparo que ocorreu na sequência de disparos x, que á a linha do neurônio especificada
ultimo <- function(x,tempo_inicial,tempo_final){
  u = tail(which(x[tempo_inicial:tempo_final] == 1), n = 1)
  return(u)
}

#Função que calcula o potencial de membrana do neurônio "indice_neurônio" no tempo t com taxa basal beta
V_novo <- function(X,indice_neuronio,tempo){
  potencial = if(ultimo(matriz_x(X,indice_neuronio),1,tempo) == tempo) beta  
  else if(ultimo(matriz_x(X,indice_neuronio),1,tempo) != (tempo-1)) beta + sum((1/(2^(tempo - ultimo(matriz_x(X,indice_neuronio),1,tempo))))*rowSums(matriz_pesos(X,indice_neuronio)[,(ultimo(matriz_x(X,indice_neuronio),1,tempo)+1):tempo])%*%matriz_pesos(w,indice_neuronio)[,indice_neuronio])
  else beta + sum((1/(2^(tempo - ultimo(matriz_x(X,indice_neuronio),1,tempo))))*matriz_pesos(X,indice_neuronio)[,(ultimo(matriz_x(X,indice_neuronio),1,tempo)+1):tempo]%*%matriz_pesos(w,indice_neuronio)[,indice_neuronio])
  return(potencial)
}


##Função que gera a cadeia X
gera_cadeia <- function(X){
  for (j in 5:(t+4)){
    for (i in 1:n){
      X[i,j] = ifelse(runif(1,0,1) <= phi(V_novo(X,i,j-1)), 1, 0)
    }
  }
  return(X)
} #sem considerar alfa


# Estimação ---------------------------------------------------------------

#Função que calcula o potencial de membrana do neurônio "indice_neurônio" no tempo t com taxa basal beta
#com a matriz de pesos como argumento, precisamos para calcular as coisas que dependem de w_est
V_geral <- function(X, indice_neuronio, tempo, W){
  potencial = if(ultimo(matriz_x(X,indice_neuronio),1,tempo) == tempo) beta  
  else if(ultimo(matriz_x(X,indice_neuronio),1,tempo) != (tempo-1)) beta + sum((1/(2^(tempo - ultimo(matriz_x(X,indice_neuronio),1,tempo))))*rowSums(matriz_pesos(X,indice_neuronio)[,(ultimo(matriz_x(X,indice_neuronio),1,tempo)+1):tempo])%*%matriz_pesos(W,indice_neuronio)[,indice_neuronio])
  else beta + sum((1/(2^(tempo - ultimo(matriz_x(X,indice_neuronio),1,tempo))))*matriz_pesos(X,indice_neuronio)[,(ultimo(matriz_x(X,indice_neuronio),1,tempo)+1):tempo]%*%matriz_pesos(W,indice_neuronio)[,indice_neuronio])
  return(potencial)
}


#Função que gera a cadeia recebendo W como parâmetro com H passos a mais no X_hist (historico observado da rede)
simula_futuro <- function(X_hist, W, H){
  Tn <- ncol(X_hist)
  X <- cbind(X_hist, matrix(0, nrow = n, ncol = H))
  for (j in (Tn + 1):(Tn + H)){
    for (i in 1:n){
      X[i, j] <- ifelse(runif(1, 0, 1) <= phi(V_geral(X, i, j - 1, W)), 1, 0)
    }
  }
  return(X)
}

#Funçaõ que calcula a probabilidade de disparo p_{i,h} para todos os neuronios i=1..n e horizontes, dado o histórico e W.

prob_disparo_todos <- function(X_hist, H, W, B){
  Tn <- ncol(X_hist)
  p <- matrix(NA_real_, nrow = n, ncol = H)
  
  for (i in 1:n) p[i, 1] <- phi(V_geral(X_hist, i, Tn, W))   # h = 1 (calculado direto do modelo)
  
  if (H > 1){
    soma <- matrix(0, nrow = n, ncol = H - 1)
    for (b in 1:B){
      X_fut <- simula_futuro(X_hist, W, H)
      soma <- soma + X_fut[, (Tn + 2):(Tn + H), drop = FALSE]
    }
    p[, 2:H] <- soma / B
  }
  return(p)   
}


#Replicas


# erros de cada replica: phat_{i,h}^(r) - p0_{i,h}^(r) 
erros <- array(NA_real_, dim = c(n, H, R_repl))
P0_todos   <- array(NA_real_, dim = c(n, H, R_repl))
Phat_todos <- array(NA_real_, dim = c(n, H, R_repl))

W_est_todos <- array(NA_real_, dim = c(n, n, R_repl))
X_todos     <- vector("list", R_repl)
conv_todos  <- matrix(NA,       nrow = n, ncol = R_repl)
iter_todos  <- matrix(NA_real_, nrow = n, ncol = R_repl)
ndisp_todos <- matrix(NA_real_, nrow = n, ncol = R_repl)

for (r in 1:R_repl){
  set.seed(100+r)
  
  X = matrix(c(1,0,0,1,
               0,1,1,0,
               0,0,1,0,
               0,1,0,1,    
               1,0,0,0), nrow=n, byrow=T) 
  
  X = cbind(X, matrix(0,nrow=n,ncol=t))
  X <- gera_cadeia(X)
  
  for (j in 1: n){
    z = matrix(0,nrow=t, ncol = n)
    for ( i in 1:t){
      X_aux = X[,1:(i+3)]
      u = ultimo(X_aux[j,],1,(i+3)) 
      z[i,]  = if (u==(i+3)) matrix(0,ncol=n,nrow=1)  else if (u==(i+2))  (1/(2^((i+3) - ultimo(matriz_x( X_aux,j),1,(i+3)))))*t(X_aux[,(ultimo(X_aux[j,],1,(i+3))+1):(i+3)]) else (1/(2^((i+3) - ultimo(matriz_x( X_aux,j),1,(i+3)))))*t(rowSums(X_aux[,(ultimo(X_aux[j,],1,(i+3))+1):(i+3)]))
    }
    z = z[-(1:(0.2*t)),]
    nome = paste0("matriz", j)
    assign(nome, z)
  }
  
  for (i in 1:n){
    resp <- X[i,(5 + 0.2*t):(t+4)]
    nome_resp <- paste0('y',i)
    assign(nome_resp,resp)
  }
  
  w_est <- matrix(0,nrow=n,ncol=n)
  for (i in 1:n){
    y_i <- get(paste0('y',i))
    mod <- glm(as.factor(y_i) ~ 0 + ., data = as.data.frame(get(paste0('matriz',i))),
               family = binomial(link = "logit"), offset = rep(beta, length(y_i)))
    w_est[,i] <- coef(mod)
    conv_todos[i, r]  <- mod$converged
    iter_todos[i, r]  <- mod$iter
    ndisp_todos[i, r] <- sum(y_i)
  }
  
  # Histórico no instante final da amostra, usamos t+4 por que começamos considerando X com 4 colunas
  F_T <- X[, 1:(t+4)]  
  
  W_est_todos[, , r] <- w_est
  X_todos[[r]]       <- X
  
  P0   <- prob_disparo_todos(F_T, H, w,     B)   # probabilidade oraculo
  Phat <- prob_disparo_todos(F_T, H, w_est, B)   # probabilidade plug-in
  
  erros[, , r] <- Phat - P0
  P0_todos[, , r]   <- P0      # oraculo da replica r
  Phat_todos[, , r] <- Phat    # plug-in da replica r
  
}


# Calculando os erros das probabilidades preditivas 
MAE_P  <- apply(abs(erros),  c(1, 2), mean)
RMSE_P <- sqrt(apply(erros^2, c(1, 2), mean))

# Probabilidades médias 
P0_medio <- apply(P0_todos, c(1, 2), mean)
Phat_medio <- apply(Phat_todos, c(1, 2), mean)

############# Cálculo da Entropia ##############

#convencao: 0*log(0) = 0 pros cados em que p=1 ou p=0

h2 <- function(p){
  ifelse(p <= 0 | p >= 1, 0, -p*log2(p) - (1-p)*log2(1-p))
}

#
U_0   <- h2(P0_todos)    # entropia oráculo por réplica 
U_hat <- h2(Phat_todos)  # entropia plug-in por réplica


erro_U <- U_hat - U_0      # erro da entropia por réplica

# Erros e Víes da entropia 
Vies_U <- apply(erro_U,      c(1, 2), mean)
MAE_U  <- apply(abs(erro_U), c(1, 2), mean)
RMSE_U <- sqrt(apply(erro_U^2, c(1, 2), mean))

# Entropia média e desvios-padrão 
U_0_medio   <- apply(U_0,   c(1, 2), mean)
U_hat_medio <- apply(U_hat, c(1, 2), mean)

U_0_dp    <- apply(U_0, c(1, 2), sd) 
Uhat_dp    <- apply(U_hat, c(1, 2), sd) 
