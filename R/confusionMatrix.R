
# modified from PresenceAbsence::cmx
cmx<- function (DATA, threshold = 0.5, which.model = 1, na.rm = FALSE)
{
  if (is.logical(na.rm) == FALSE) {
    stop("'na.rm' must be of logical type")
  }
  if (sum(is.na(DATA)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(DATA), 1, sum)
      warning(length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
      DATA <- DATA[NA.rows == 0, ]
    }
    else {
      return(NA)
    }
  }
  if (length(which.model) != 1) {
    stop("this function will only work for a single model, 'which.model' must be of length one")
  }
  if (which.model < 1 || round(which.model) != which.model) {
    stop("'which.model' must be a positive integer")
  }
  if (which.model + 2 > ncol(DATA)) {
    stop("'which.model' must not be greater than number of models in 'DATA'")
  }
  if (length(threshold) != 1) {
    stop("'threshold' must be a single number between zero and one")
  }
  if (max(threshold) > 1) {
    stop("'threshold' must be a single number between zero and one")
  }
  if (min(threshold) < 0) {
    stop("'threshold' must be a single number between zero and one")
  }
  # OBS.ind <- DATA[, 2] > 0
  OBS.ind<- DATA[, 2]
  if (threshold == 0) {
    PRED.ind = DATA[, which.model + 2] >= threshold
  }
  else {
    #PRED.ind = DATA[, which.model + 2] > threshold
    PRED.ind = ifelse(DATA[, which.model + 2] > threshold, 1, 0)
  }
  C = data.frame(predicted=PRED.ind,
                 observed=OBS.ind)
  #C = as.table(matrix(C, nrow = 2))
  #dimnames(C) = list(predicted = c(1, 0), observed = c(1, 0))
  #storage.mode(C) = "double"
  return(C)
}

# modiefied from mecofun::evalSDM: https://gitup.uni-potsdam.de/macroecology/mecofun/-/blob/master/R/evalSDM.R
confusionMatrix <- function(observation, predictions, thresh=NULL, thresh.method='MaxSens+Spec', req.sens=0.85, req.spec = 0.85, FPC=1, FNC=1, weigths=rep(1, length(observation))){

  thresh.dat <- data.frame(ID=seq_len(length(observation)),
                           obs = observation,
                           pred = predictions)

  if (is.null(thresh)) {
    thresh.mat <- PresenceAbsence::optimal.thresholds(DATA= thresh.dat, req.sens=req.sens, req.spec = req.spec, FPC=FPC, FNC=FNC)
    thresh <- thresh.mat[thresh.mat$Method==thresh.method,2]
  }

  cmx.opt <- cmx(DATA= thresh.dat, threshold=thresh)
  return(cmx.opt)
}


