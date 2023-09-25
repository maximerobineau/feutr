###################################################################
########################## Fonctions ##############################
###################################################################

# Les fonctions ont toutes le préfixe "feutr" (pour "Forecast efficiency and unbiasedness tests in R").
# La notation est typiquement celle utilisée pour des méthodes plutôt que des fonctions, en vue de la transition 
# éventuelle du code vers cette structure plus universelle pour des paquetages


# Fonction pour mettre les codes de significativité de R dans les tableaux et les sommaires avec plusieurs h
feutr.signif <- function(arg){ 
  if (is.na(arg)==TRUE) {
    return(" ")
  } else if(arg < 0.001){
    return("***")
  } else if(arg < 0.01){
    return("**")
  } else if(arg < 0.05){
    return("*")
  } else if(arg < 0.1){
    return(".")
  } else {
    return(" ")
  }
}

# Fonction pour que les "diagonales" appropriées soient des colonnes 
# Ainsi, la diagonale principale devient la colonne des prévision à h=0 périodes à l'avance, 
# La diagonale en-dessous devient la colonne des prévision à h=1 périodes à l'avance, etc. 
feutr.ForecastData <- function(var, hmax, avh = 2){
  ForecastData <- data.frame(1:nrow(var))
  for (i in 1:(ncol(var)-(avh+1))){
    ForecastData[i, 1] <- var[i,i+avh+1]
  }
  for(k in 0:hmax){
    for (i in 1:(ncol(var)-1)){
      ForecastData[i+k, k+2] <- var[i+k,i+1]
    }
  }
  ForecastData <- head(ForecastData, ncol(var)-(avh+1))
  colnames(ForecastData) <- c("Actual", 0:hmax)
  ForecastData
}

# Fonction de calcul des erreurs prévisionnelles
feutr.ErrorData <- function(var, hmax, avh = 2){
  ForecastData <- feutr.ForecastData(var = var, hmax = hmax, avh = avh)
  ErrorData <- data.frame(ForecastData[,1]) 
  for (i in 0:hmax) {
    ErrorData[,i+2] <- ForecastData[,1]-ForecastData[,i+2]
  }
  colnames(ErrorData) <- c("Actual", 0:hmax)
  ErrorData
} 

# Fonction de calcul des révisions
feutr.RevisionData <- function(var, hmax, avh = 2){
  ForecastData <- feutr.ForecastData(var = var, hmax = hmax, avh = avh)
  RevisionData <- data.frame(ForecastData[,1])
  for (i in 0:(hmax-1)) {
    RevisionData[,i+2] <- ForecastData[,i+2]-ForecastData[,i+3]
  }
  colnames(RevisionData) <- c("Actual", 0:(hmax-1))
  RevisionData
}

# Fonction du test de Mincer et Zarnowitz (1969)
feutr.mz <- function(var, h, avh = 2, binary = NULL){
  if (is.null(binary)==TRUE) {
    if (length(h)==1) {
      mz <- lm(feutr.ForecastData(var = var, hmax = h, avh = avh)[,1]~feutr.ForecastData(var = var, hmax = h, avh = avh)[,(h+2)]);mz
    } else {
      mz <- data.frame(1:6)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,1]~feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,(h[k]+2)]))
        mz[1,k] <- round(results$coefficients[1], digits = 4)
        mz[2,k] <- round(results$coefficients[3], digits = 4)
        mz[3,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[1,k])-0)/as.numeric(mz[2,k])), results$df[2], lower.tail = FALSE))
        mz[4,k] <- round(results$coefficients[2], digits = 4)
        mz[5,k] <- round(results$coefficients[4], digits = 4)
        mz[6,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[4,k])-1)/as.numeric(mz[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(mz) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(mz) <- h
      mz
    }
  } else {
    if (length(h)==1) {
      mz <- lm(feutr.ForecastData(var = var, hmax = h, avh = avh)[,1]~feutr.ForecastData(var = var, hmax = h, avh = avh)[,(h+2)]+binary+I(binary*feutr.ForecastData(var = var, hmax = h, avh = avh)[,(h+2)]));mz
    } else {
      mz <- data.frame(1:12)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,1]~feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,(h[k]+2)]+binary+I(binary*feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,(h[k]+2)])))
        mz[1,k] <- round(results$coefficients[1], digits = 4)
        mz[2,k] <- round(results$coefficients[5], digits = 4)
        mz[3,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[1,k])-0)/as.numeric(mz[2,k])), results$df[2], lower.tail = FALSE))
        mz[4,k] <- round(results$coefficients[2], digits = 4)
        mz[5,k] <- round(results$coefficients[6], digits = 4)
        mz[6,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[4,k])-1)/as.numeric(mz[5,k])), results$df[2], lower.tail = FALSE))
        mz[7,k] <- round(results$coefficients[3], digits = 4)
        mz[8,k] <- round(results$coefficients[7], digits = 4)
        mz[9,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[7,k])-0)/as.numeric(mz[8,k])), results$df[2], lower.tail = FALSE))
        mz[10,k] <- round(results$coefficients[4], digits = 4)
        mz[11,k] <- round(results$coefficients[8], digits = 4)
        mz[12,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[10,k])-1)/as.numeric(mz[11,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(mz) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency", "Bias (binary)", "S.E. of bias (binary)", "Signif. of bias (binary)", "Efficiency (binary)", "S.E. of efficiency (binary)", "Signif. of efficiency (binary)")
      colnames(mz) <- h
      mz
    }
  }   
}

# Fonction du test de Holden et Peel (1990)
feutr.hp <- function(var, h, avh = 2, binary = NULL){
  if (is.null(binary)==TRUE) {
    if (length(h)==1) {
      mz <- lm(feutr.ForecastData(var = var, hmax = h, avh = avh)[,1]-feutr.ForecastData(var = var, hmax = h, avh = avh)[,(h+2)]~1);mz
    } else {
      mz <- data.frame(1:3)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,1]-feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,(h[k]+2)]~1))
        mz[1,k] <- round(results$coefficients[1], digits = 4)
        mz[2,k] <- round(results$coefficients[2], digits = 4)
        mz[3,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[1,k]))/as.numeric(mz[2,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(mz) <- c("Bias", "S.E. of bias", "Signif. of bias")
      colnames(mz) <- h
      mz
    }
  } else {
    if (length(h)==1) {
      mz <- lm(feutr.ForecastData(var = var, hmax = h, avh = avh)[,1]-feutr.ForecastData(var = var, hmax = h, avh = avh)[,(h+2)]~binary);mz
    } else {
      mz <- data.frame(1:6)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,1]-feutr.ForecastData(var = var, hmax = h[k], avh = avh)[,(h[k]+2)]~binary))
        mz[1,k] <- round(results$coefficients[1], digits = 4)
        mz[2,k] <- round(results$coefficients[3], digits = 4)
        mz[3,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[1,k])-0)/as.numeric(mz[2,k])), results$df[2], lower.tail = FALSE))
        mz[4,k] <- round(results$coefficients[2], digits = 4)
        mz[5,k] <- round(results$coefficients[4], digits = 4)
        mz[6,k] <- feutr.signif(2 * pt(abs((as.numeric(mz[4,k])-1)/as.numeric(mz[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(mz) <- c("Bias", "S.E. of bias", "Signif. of bias", "Bias (binary)", "S.E. of bias (binary)", "Signif. of bias (binary)")
      colnames(mz) <- h
      mz
    }
  }
}

# Fonction du test de Nordhaus modifié 
feutr.mn <- function(var, h, avh = 2, binary = NULL){
  if (is.null(binary)==TRUE) {
    if (length(h)==1) {
      mn <- lm(feutr.ErrorData(var = var, hmax = h, avh = avh)[,h+2]~feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]);mn
    } else {
      mn <- data.frame(1:6)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ErrorData(var = var, hmax = h[k], avh = avh)[,h[k]+2]~feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]))
        mn[1,k] <- round(results$coefficients[1], digits = 4)
        mn[2,k] <- round(results$coefficients[3], digits = 4)
        mn[3,k] <- feutr.signif(2 * pt(abs(as.numeric(mn[1,k])/as.numeric(mn[2,k])), results$df[2], lower.tail = FALSE))
        mn[4,k] <- round(results$coefficients[2], digits = 4)
        mn[5,k] <- round(results$coefficients[4], digits = 4)
        mn[6,k] <- feutr.signif(2 * pt(abs(as.numeric(mn[4,k])/as.numeric(mn[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(mn) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(mn) <- h
      mn
    }
  } else {
    if (length(h)==1) {
      mn <- lm(feutr.ErrorData(var = var, hmax = h, avh = avh)[,h+2]~feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]+binary+I(binary*feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]));mn
    } else {
      mn <- data.frame(1:12)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ErrorData(var = var, hmax = h[k], avh = avh)[,h[k]+2]~feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]+binary+I(binary*feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2])))
        mn[1,k] <- round(results$coefficients[1], digits = 4)
        mn[2,k] <- round(results$coefficients[5], digits = 4)
        mn[3,k] <- feutr.signif(2 * pt(abs(as.numeric(mn[1,k])/as.numeric(mn[2,k])), results$df[2], lower.tail = FALSE))
        mn[4,k] <- round(results$coefficients[2], digits = 4)
        mn[5,k] <- round(results$coefficients[6], digits = 4)
        mn[6,k] <- feutr.signif(2 * pt(abs(as.numeric(mn[4,k])/as.numeric(mn[5,k])), results$df[2], lower.tail = FALSE))
        mn[7,k] <- round(results$coefficients[3], digits = 4)
        mn[8,k] <- round(results$coefficients[7], digits = 4)
        mn[9,k] <- feutr.signif(2 * pt(abs(as.numeric(mn[7,k])/as.numeric(mn[8,k])), results$df[2], lower.tail = FALSE))
        mn[10,k] <- round(results$coefficients[4], digits = 4)
        mn[11,k] <- round(results$coefficients[8], digits = 4)
        mn[12,k] <- feutr.signif(2 * pt(abs(as.numeric(mn[10,k])/as.numeric(mn[11,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(mn) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency", "Bias (binary)", "S.E. of bias (binary)", "Signif. of bias (binary)", "Efficiency (binary)", "S.E. of efficiency (binary)", "Signif. of efficiency (binary)")
      colnames(mn) <- h
      mn
    }
  }
}

# Fonction du test de Nordhaus modifié avec la nouvelle spécification
feutr.rm.set <- function(var, h, avh = 2, slope.null = 1, binary = NULL){
  if (is.null(binary)==TRUE) {
    if (length(h)==1) {
      rm <- lm(feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3]~feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]);rm
    } else {
      rm <- data.frame(1:6)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3]~feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]))
        rm[1,k] <- round(results$coefficients[1], digits = 4)
        rm[2,k] <- round(results$coefficients[3], digits = 4)
        rm[3,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[1,k])-0)/as.numeric(rm[2,k])), results$df[2], lower.tail = FALSE))
        rm[4,k] <- round(results$coefficients[2], digits = 4)
        rm[5,k] <- round(results$coefficients[4], digits = 4)
        rm[6,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[4,k])-slope.null)/as.numeric(rm[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(rm) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(rm) <- h
      rm
    }
  } else {
    if (length(h)==1) {
      rm <- lm(feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3]~feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]+binary+I(binary*feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]));rm
    } else {
      rm <- data.frame(1:12)
      for (k in 1:length(h)) {
        results <- summary(lm(feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3]~feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]+binary+I(binary*feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2])))
        rm[1,k] <- round(results$coefficients[1], digits = 4)
        rm[2,k] <- round(results$coefficients[5], digits = 4)
        rm[3,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[1,k])-0)/as.numeric(rm[2,k])), results$df[2], lower.tail = FALSE))
        rm[4,k] <- round(results$coefficients[2], digits = 4)
        rm[5,k] <- round(results$coefficients[6], digits = 4)
        rm[6,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[4,k])-slope.null)/as.numeric(rm[5,k])), results$df[2], lower.tail = FALSE))
        rm[7,k] <- round(results$coefficients[3], digits = 4)
        rm[8,k] <- round(results$coefficients[7], digits = 4)
        rm[9,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[7,k])-0)/as.numeric(rm[8,k])), results$df[2], lower.tail = FALSE))
        rm[10,k] <- round(results$coefficients[4], digits = 4)
        rm[11,k] <- round(results$coefficients[8], digits = 4)
        rm[12,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[10,k])-slope.null)/as.numeric(rm[11,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(rm) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency", "Bias (binary)", "S.E. of bias (binary)", "Signif. of bias (binary)", "Efficiency (binary)", "S.E. of efficiency (binary)", "Signif. of efficiency (binary)")
      colnames(rm) <- h
      rm
    }
  }
}

# Fonction du nouveau test de rigidité pour le sous-ensemble des révisions dans le bon sens
feutr.rm.ImprovedSubset <- function(var, h, avh = 2, binary = NULL){
  if (is.null(binary)==TRUE) {
    if (length(h)==1) {
      ImprovedSubset <- feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]/feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3]>0
      rm <- lm(subset(feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3],ImprovedSubset)~subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2],ImprovedSubset));rm
    } else {
      rm <- data.frame(1:6)
      for (k in 1:length(h)) {  
        ImprovedSubset <- feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]/feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3]>0
        results <- summary(lm(subset(feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3],ImprovedSubset)~subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2],ImprovedSubset)))
        rm[1,k] <- round(results$coefficients[1], digits = 4)
        rm[2,k] <- round(results$coefficients[3], digits = 4)
        rm[3,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[1,k])-0)/as.numeric(rm[2,k])), results$df[2], lower.tail = FALSE))
        rm[4,k] <- round(results$coefficients[2], digits = 4)
        rm[5,k] <- round(results$coefficients[4], digits = 4)
        rm[6,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[4,k])-1)/as.numeric(rm[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(rm) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(rm) <- h
      rm
    }
  } else {
    if (length(h)==1) {
      ImprovedSubset <- feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]/feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3]>0
      rm <- lm(subset(feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3],ImprovedSubset)~subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2],ImprovedSubset)+subset(binary,ImprovedSubset)+I(subset(binary,ImprovedSubset)*subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2],ImprovedSubset)))
    } else {
      rm <- data.frame(1:12)
      for (k in 1:length(h)) {
        ImprovedSubset <- feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]/feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3]>0
        results <- summary(lm(subset(feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3],ImprovedSubset)~subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2],ImprovedSubset)+subset(binary,ImprovedSubset)+I(subset(binary,ImprovedSubset)*subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2],ImprovedSubset))))
        rm[1,k] <- round(results$coefficients[1], digits = 4)
        rm[2,k] <- round(results$coefficients[5], digits = 4)
        rm[3,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[1,k])-0)/as.numeric(rm[2,k])), results$df[2], lower.tail = FALSE))
        rm[4,k] <- round(results$coefficients[2], digits = 4)
        rm[5,k] <- round(results$coefficients[6], digits = 4)
        rm[6,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[4,k])-1)/as.numeric(rm[5,k])), results$df[2], lower.tail = FALSE))
        rm[7,k] <- round(results$coefficients[3], digits = 4)
        rm[8,k] <- round(results$coefficients[7], digits = 4)
        rm[9,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[7,k])-0)/as.numeric(rm[8,k])), results$df[2], lower.tail = FALSE))
        rm[10,k] <- round(results$coefficients[4], digits = 4)
        rm[11,k] <- round(results$coefficients[8], digits = 4)
        rm[12,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[10,k])-1)/as.numeric(rm[11,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(rm) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency", "Bias (binary)", "S.E. of bias (binary)", "Signif. of bias (binary)", "Efficiency (binary)", "S.E. of efficiency (binary)", "Signif. of efficiency (binary)")
      colnames(rm) <- h
      rm
    }
  }
}

# Fonction du nouveau test de rigidité pour le sous-ensemble des révisions dans le mauvais sens
feutr.rm.WorsenedSubset <- function(var, h, avh = 2, binary = NULL){
  if (is.null(binary)==TRUE) {
    if (length(h)==1) {
      WorsenedSubset <- feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]/feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3]<0
      rm <- lm(subset(feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3],WorsenedSubset)~subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2],WorsenedSubset));rm
    } else {
      rm <- data.frame(1:6)
      for (k in 1:length(h)) {  
        WorsenedSubset <- feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]/feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3]<0
        results <- summary(lm(subset(feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3],WorsenedSubset)~subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2],WorsenedSubset)))
        rm[1,k] <- round(results$coefficients[1], digits = 4)
        rm[2,k] <- round(results$coefficients[3], digits = 4)
        rm[3,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[1,k])-0)/as.numeric(rm[2,k])), results$df[2], lower.tail = FALSE))
        rm[4,k] <- round(results$coefficients[2], digits = 4)
        rm[5,k] <- round(results$coefficients[4], digits = 4)
        rm[6,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[4,k])+1)/as.numeric(rm[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(rm) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(rm) <- h
      rm
    }
  } else {
    if (length(h)==1) {
      WorsenedSubset <- feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]/feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3]<0
      rm <- lm(subset(feutr.ErrorData(var = var, hmax = h+1, avh = avh)[,h+3],WorsenedSubset)~subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2],WorsenedSubset)+subset(binary,WorsenedSubset)+I(subset(binary,WorsenedSubset)*subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2],WorsenedSubset)))
    } else {
      rm <- data.frame(1:12)
      for (k in 1:length(h)) {
        WorsenedSubset <- feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]/feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3]<0
        results <- summary(lm(subset(feutr.ErrorData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+3],WorsenedSubset)~subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2],WorsenedSubset)+subset(binary,WorsenedSubset)+I(subset(binary,WorsenedSubset)*subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2],WorsenedSubset))))
        rm[1,k] <- round(results$coefficients[1], digits = 4)
        rm[2,k] <- round(results$coefficients[5], digits = 4)
        rm[3,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[1,k])-0)/as.numeric(rm[2,k])), results$df[2], lower.tail = FALSE))
        rm[4,k] <- round(results$coefficients[2], digits = 4)
        rm[5,k] <- round(results$coefficients[6], digits = 4)
        rm[6,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[4,k])+1)/as.numeric(rm[5,k])), results$df[2], lower.tail = FALSE))
        rm[7,k] <- round(results$coefficients[3], digits = 4)
        rm[8,k] <- round(results$coefficients[7], digits = 4)
        rm[9,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[7,k])-0)/as.numeric(rm[8,k])), results$df[2], lower.tail = FALSE))
        rm[10,k] <- round(results$coefficients[4], digits = 4)
        rm[11,k] <- round(results$coefficients[8], digits = 4)
        rm[12,k] <- feutr.signif(2 * pt(abs((as.numeric(rm[10,k])+1)/as.numeric(rm[11,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(rm) <- c("Bias", "S.E. of bias", "Signif. of bias", "Efficiency", "S.E. of efficiency", "Signif. of efficiency", "Bias (binary)", "S.E. of bias (binary)", "Signif. of bias (binary)", "Efficiency (binary)", "S.E. of efficiency (binary)", "Signif. of efficiency (binary)")
      colnames(rm) <- h
      rm
    }
  }
}

# Fonction "mère" appelant les trois nouveaux tests présentés
feutr.rm <- function(var, h, avh = 2, slope.null = 1, subset = FALSE, ImprovedSubset = TRUE, binary = NULL){
  if (subset == FALSE) {
    rm <- feutr.rm.set(var = var, h = h, slope.null = slope.null, avh = avh, binary = binary);rm
  } else {
    if (ImprovedSubset == TRUE) {
      rm <- feutr.rm.ImprovedSubset(var = var, h = h, avh = avh, binary = binary);rm
    } else {
      rm <- feutr.rm.WorsenedSubset(var = var, h = h, avh = avh, binary = binary);rm
    }
  }
}

# Fonction du test d'Ashiya 
feutr.ash <- function(var, h, avh = 2, positive = TRUE){
  if (positive==TRUE) {
    if (length(h)==1) {
      PositiveSubset <- feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]>0
      ash <- lm(subset(feutr.ErrorData(var = var, hmax = h, avh = avh)[,h+2], PositiveSubset)~subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2], PositiveSubset));ash
    } else {
      ash <- data.frame(1:6)
      for (k in 1:length(h)) {
        PositiveSubset <- feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]>0
        results <- summary(lm(subset(feutr.ErrorData(var = var, hmax = h[k], avh = avh)[,h[k]+2], PositiveSubset)~subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2], PositiveSubset)))
        ash[1,k] <- round(results$coefficients[1], digits = 4)
        ash[2,k] <- round(results$coefficients[3], digits = 4)
        ash[3,k] <- feutr.signif(2 * pt(abs(as.numeric(ash[1,k])/as.numeric(ash[2,k])), results$df[2], lower.tail = FALSE))
        ash[4,k] <- round(results$coefficients[2], digits = 4)
        ash[5,k] <- round(results$coefficients[4], digits = 4)
        ash[6,k] <- feutr.signif(2 * pt(abs(as.numeric(ash[4,k])/as.numeric(ash[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(ash) <- c("Sentiment", "S.E. of sentiment", "Signif. of sentiment", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(ash) <- h
      ash
    }
  } else {
    if (length(h)==1) {
      NegativeSubset <- feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2]<0
      ash <- lm(subset(feutr.ErrorData(var = var, hmax = h, avh = avh)[,h+2], NegativeSubset)~subset(feutr.RevisionData(var = var, hmax = h+1, avh = avh)[,h+2], NegativeSubset));ash
    } else {
      ash <- data.frame(1:6)
      for (k in 1:length(h)) {
        NegativeSubset <- feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2]<0
        results <- summary(lm(subset(feutr.ErrorData(var = var, hmax = h[k], avh = avh)[,h[k]+2], NegativeSubset)~subset(feutr.RevisionData(var = var, hmax = h[k]+1, avh = avh)[,h[k]+2], NegativeSubset)))
        ash[1,k] <- round(results$coefficients[1], digits = 4)
        ash[2,k] <- round(results$coefficients[3], digits = 4)
        ash[3,k] <- feutr.signif(2 * pt(abs(as.numeric(ash[1,k])/as.numeric(ash[2,k])), results$df[2], lower.tail = FALSE))
        ash[4,k] <- round(results$coefficients[2], digits = 4)
        ash[5,k] <- round(results$coefficients[4], digits = 4)
        ash[6,k] <- feutr.signif(2 * pt(abs(as.numeric(ash[4,k])/as.numeric(ash[5,k])), results$df[2], lower.tail = FALSE))
      }
      rownames(ash) <- c("Sentiment", "S.E. of sentiment", "Signif. of sentiment", "Efficiency", "S.E. of efficiency", "Signif. of efficiency")
      colnames(ash) <- h
      ash
    }
  }
}

# Fonction "mère" appelant tous les tests précédents
feutr <- function(var, h, avh = 2, test = "rm", slope.null = 1, subset = FALSE, ImprovedSubset = TRUE, binary = NULL, positive = TRUE){
  if (test == "rm") {
    ans <- feutr.rm(var = var, h = h, avh = avh, slope.null = slope.null, subset = subset, ImprovedSubset = ImprovedSubset, binary = binary);ans
  } else if (test == "mz") {
    ans <- feutr.mz(var = var, h = h, avh = avh, binary = binary);ans
  } else if (test == "hp") {
    ans <- feutr.hp(var = var, h = h, avh = avh, binary = binary);ans
  } else if (test == "mn") {
    ans <- feutr.mn(var = var, h = h, avh = avh, binary = binary);ans
  } else if (test == "ash") {
    ans <- feutr.ash(var = var, h = h, avh = avh, positive = positive);ans
  }
}

# Nouvelle fonction sommaire pour tester contre les bonnes hypothèses nulles et nommer les coefficients
feutr.summary <- function(test, slope.null=1){
  if (is.data.frame(test)==TRUE) {
    test
  } else if (as.character(test$call[2]) == "feutr.ForecastData(var = var, hmax = h, avh = avh)[, 1] - feutr.ForecastData(var = var, hmax = h, avh = avh)[, (h + 2)] ~ 1") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias")
    ans
  } else if (as.character(test$call[2]) == "feutr.ForecastData(var = var, hmax = h, avh = avh)[, 1] - feutr.ForecastData(var = var, hmax = h, avh = avh)[, (h + 2)] ~ binary") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Bias (binary)")
    ans
  } else if (as.character(test$call[2]) == "feutr.ForecastData(var = var, hmax = h, avh = avh)[, 1] ~ feutr.ForecastData(var = var, hmax = h, avh = avh)[, (h + 2)]") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency")
    ans$coefficients[6] <- (ans$coefficients[2]-1)/ans$coefficients[4]
    ans$coefficients[8] <- 2 * pt(abs(ans$coefficients[6]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "feutr.ForecastData(var = var, hmax = h, avh = avh)[, 1] ~ feutr.ForecastData(var = var, hmax = h, avh = avh)[, (h + 2)] + binary + I(binary * feutr.ForecastData(var = var, hmax = h, avh = avh)[, (h + 2)])") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency","Bias (binary)","Efficiency (binary)")
    ans$coefficients[10] <- (ans$coefficients[2]-1)/ans$coefficients[6]
    ans$coefficients[12] <- (ans$coefficients[4]-1)/ans$coefficients[8]
    ans$coefficients[14] <- 2 * pt(abs(ans$coefficients[10]), ans$df[2], lower.tail = FALSE)
    ans$coefficients[16] <- 2 * pt(abs(ans$coefficients[12]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "feutr.ErrorData(var = var, hmax = h, avh = avh)[, h + 2] ~ feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2]") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency")
    ans
  } else if (as.character(test$call[2]) == "feutr.ErrorData(var = var, hmax = h, avh = avh)[, h + 2] ~ feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2] + binary + I(binary * feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2])") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency","Bias (binary)","Efficiency (binary)")
    ans
  } else if (as.character(test$call[2]) == "feutr.ErrorData(var = var, hmax = h + 1, avh = avh)[, h + 3] ~ feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2]") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency")
    ans$coefficients[6] <- (ans$coefficients[2]-slope.null)/ans$coefficients[4]
    ans$coefficients[8] <- 2 * pt(abs(ans$coefficients[6]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "feutr.ErrorData(var = var, hmax = h + 1, avh = avh)[, h + 3] ~ feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2] + binary + I(binary * feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2])") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency","Bias (binary)","Efficiency (binary)")
    ans$coefficients[10] <- (ans$coefficients[2]-slope.null)/ans$coefficients[6]
    ans$coefficients[12] <- (ans$coefficients[4]-slope.null)/ans$coefficients[8]
    ans$coefficients[14] <- 2 * pt(abs(ans$coefficients[10]), ans$df[2], lower.tail = FALSE)
    ans$coefficients[16] <- 2 * pt(abs(ans$coefficients[12]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "subset(feutr.ErrorData(var = var, hmax = h + 1, avh = avh)[, h + 3], ImprovedSubset) ~ subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], ImprovedSubset)") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency")
    ans$coefficients[6] <- (ans$coefficients[2]-1)/ans$coefficients[4]
    ans$coefficients[8] <- 2 * pt(abs(ans$coefficients[6]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "subset(feutr.ErrorData(var = var, hmax = h + 1, avh = avh)[, h + 3], ImprovedSubset) ~ subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], ImprovedSubset) + subset(binary, ImprovedSubset) + I(subset(binary, ImprovedSubset) * subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], ImprovedSubset))") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency","Bias (binary)","Efficiency (binary)")
    ans$coefficients[10] <- (ans$coefficients[2]-1)/ans$coefficients[6]
    ans$coefficients[12] <- (ans$coefficients[4]-1)/ans$coefficients[8]
    ans$coefficients[14] <- 2 * pt(abs(ans$coefficients[10]), ans$df[2], lower.tail = FALSE)
    ans$coefficients[16] <- 2 * pt(abs(ans$coefficients[12]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "subset(feutr.ErrorData(var = var, hmax = h + 1, avh = avh)[, h + 3], WorsenedSubset) ~ subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], WorsenedSubset)") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency")
    ans$coefficients[6] <- (ans$coefficients[2]+1)/ans$coefficients[4]
    ans$coefficients[8] <- 2 * pt(abs(ans$coefficients[6]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "subset(feutr.ErrorData(var = var, hmax = h + 1, avh = avh)[, h + 3], WorsenedSubset) ~ subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], WorsenedSubset) + subset(binary, WorsenedSubset) + I(subset(binary, WorsenedSubset) * subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], WorsenedSubset))") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Bias","Efficiency","Bias (binary)","Efficiency (binary)")
    ans$coefficients[10] <- (ans$coefficients[2]+1)/ans$coefficients[6]
    ans$coefficients[12] <- (ans$coefficients[4]+1)/ans$coefficients[8]
    ans$coefficients[14] <- 2 * pt(abs(ans$coefficients[10]), ans$df[2], lower.tail = FALSE)
    ans$coefficients[16] <- 2 * pt(abs(ans$coefficients[12]), ans$df[2], lower.tail = FALSE)
    ans
  } else if (as.character(test$call[2]) == "subset(feutr.ErrorData(var = var, hmax = h, avh = avh)[, h + 2], PositiveSubset) ~ subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], PositiveSubset)") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Sentiment","Efficiency")
    ans
  } else if (as.character(test$call[2]) == "subset(feutr.ErrorData(var = var, hmax = h, avh = avh)[, h + 2], NegativeSubset) ~ subset(feutr.RevisionData(var = var, hmax = h + 1, avh = avh)[, h + 2], NegativeSubset)") {
    ans <- summary(test)
    rownames(ans$coefficients) <- c("Sentiment","Efficiency")
    ans
  } else {
    warning("Forecast efficiency and unbiasedness test not recongnized. Using summary() from {base}.")
    summary(test)
  }
}