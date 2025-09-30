#matriz de métricas de diversidade_microbacias nível 7
metricas<- read.csv("metricas_characiformes.csv", header=T)
#matriz de métricas de diversidade_microbacias nível 4
metricas_lev4<- read.csv("metricas_lev4.csv", header = T)
#matriz de métricas de diversidade_microbacias nível 5
metricas_lev5<- read.csv("metricas_lev5", header= T)


#baixe o pacote
devtools::install_github("GabrielNakamura/Dimensionality_package", force = TRUE)
#carregue o pacote
library(Dimensionality)


##Calculando a dimensionalidade para a bacia

EE_metric <- Dimensionality::EvennessEigen(matrix.M = as.matrix(metricas),
                                           scale = TRUE,
                                           method = "standardize",
                                           evenness = "Camargo")
IVs_diversity <- Dimensionality::ImportanceVal(matrix.M = as.matrix(metricas),
                                               scale= TRUE,
                                               method = "max",
                                               stopRule = TRUE)



######################################
##### Esse código está adaptado para que a função funcione com a presença
##### microbacias para as quais não há informações de métricas. 


### análises microbacias lev 4


metrics_com<- metricas_lev4[,-c(1,9,10,11)]
rownames(metrics_com)<- NULL
metric_matrix<- metricas_lev4[,-c(1,9,10)]
vetor_com<- metric_matrix$HYBAS_ID

metrics_com <- as.data.frame(lapply(metrics_com, function(col) as.numeric(as.character(col))))
metric_matrix <- as.data.frame(lapply(metric_matrix, function(col) as.numeric(as.character(col))))



dimensionality_level4 <- matrix(
  data = unlist(lapply(unique(vetor_com), function(x) {
    sub_mat <- metrics_com[which(metric_matrix$HYBAS_ID == x), ]
    
    # Converte para matriz
    mat_x <- as.matrix(sub_mat)
    
    # Remove colunas com variância zero (evita divisão por zero)
    mat_x <- mat_x[, apply(mat_x, 2, function(col) sd(col, na.rm = TRUE) != 0), drop = FALSE]
    
    # Remove linhas com NA, NaN ou Inf
    mat_x <- mat_x[complete.cases(mat_x) & apply(mat_x, 1, function(row) all(is.finite(row))), , drop = FALSE]
    
    # Se matriz final está vazia, retorna NA
    if (nrow(mat_x) == 0 || ncol(mat_x) == 0) {
      return(rep(NA, ncol(metrics_com) + 1))
    }
    
    # Calcula as métricas
    EE_metric_com <- Dimensionality::EvennessEigen(
      matrix.M = mat_x,
      scale = TRUE,
      method = "standardize",
      evenness = "Camargo"
    )
    
    IVs_diversity_com <- Dimensionality::ImportanceVal(
      matrix.M = mat_x,
      scale = TRUE,
      method = "max",
      stopRule = TRUE
    )
    
    if (is.matrix(IVs_diversity_com$IV.obs_stopRule)) {
      IVs_total_com <- colSums(IVs_diversity_com$IV.obs_stopRule)
    } else {
      IVs_total_com <- IVs_diversity_com$IV.obs_stopRule
    }
    
    dimensionality_obs <- c(EE_metric_com, IVs_total_com)
    
    # Completa com NA se tiver menos colunas do que o total esperado
    length_diff <- (ncol(metrics_com)) - length(IVs_total_com)
    if (length_diff > 0) {
      dimensionality_obs <- c(dimensionality_obs, rep(NA, length_diff))
    }
    
    return(dimensionality_obs)
  })),
  nrow = length(unique(vetor_com)),
  ncol = ncol(metrics_com) + 1,
  byrow = TRUE,
  dimnames = list(unique(vetor_com), c("EE", colnames(metrics_com)))
)


##################################################
#análises microbacias lev_5


metrics_com2<- metricas_lev5[,-c(1,9,10,11)]
rownames(metrics_com)<- NULL
metric_matrix2<- metricas_lev5[,-c(1,9,10)]
vetor_com2<- metric_matrix2$HYBAS_ID

dimensionality_level5 <- matrix(
  data = unlist(lapply(unique(vetor_com2), function(x) {
    sub_mat <- metrics_com2[which(metric_matrix2$HYBAS_ID == x), ]
    
    # Converte para matriz
    mat_x <- as.matrix(sub_mat)
    
    # Remove colunas com variância zero (evita divisão por zero)
    mat_x <- mat_x[, apply(mat_x, 2, function(col) sd(col, na.rm = TRUE) != 0), drop = FALSE]
    
    # Remove linhas com NA, NaN ou Inf
    mat_x <- mat_x[complete.cases(mat_x) & apply(mat_x, 1, function(row) all(is.finite(row))), , drop = FALSE]
    
    # Se matriz final está vazia, retorna NA
    if (nrow(mat_x) == 0 || ncol(mat_x) == 0) {
      return(rep(NA, ncol(metrics_com2) + 1))
    }
    
    # Calcula as métricas
    EE_metric_com <- Dimensionality::EvennessEigen(
      matrix.M = mat_x,
      scale = TRUE,
      method = "standardize",
      evenness = "Camargo"
    )
    
    IVs_diversity_com <- Dimensionality::ImportanceVal(
      matrix.M = mat_x,
      scale = TRUE,
      method = "max",
      stopRule = TRUE
    )
    
    if (is.matrix(IVs_diversity_com$IV.obs_stopRule)) {
      IVs_total_com <- colSums(IVs_diversity_com$IV.obs_stopRule)
    } else {
      IVs_total_com <- IVs_diversity_com$IV.obs_stopRule
    }
    
    dimensionality_obs <- c(EE_metric_com, IVs_total_com)
    
    # Completa com NA se tiver menos colunas do que o total esperado
    length_diff <- (ncol(metrics_com2)) - length(IVs_total_com)
    if (length_diff > 0) {
      dimensionality_obs <- c(dimensionality_obs, rep(NA, length_diff))
    }
    
    return(dimensionality_obs)
  })),
  nrow = length(unique(vetor_com2)),
  ncol = ncol(metrics_com2) + 1,
  byrow = TRUE,
  dimnames = list(unique(vetor_com2), c("EE", colnames(metrics_com2)))
)
