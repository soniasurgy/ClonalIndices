norm_eblup100 <- function(dframe, blups, meanvec, criteria){
  eblup_norm <- dframe
  for (col in blups) {
    if (!is.null(meanvec)){
      eblup_norm[[col]] <- eblup_norm[[col]] / meanvec[[col]] * 100
    }
    if (!is.null(criteria)){
      eblup_norm[[col]] <- eblup_norm[[col]] * criteria[[col]]
    }
  }
  return(eblup_norm)
}

norm_eblup <- function(dframe, cols, meanvec, criteria){
  eblup_norm <- dframe
  for (col in cols) {
    if (!is.null(meanvec)){
      eblup_norm[[col]] <- eblup_norm[[col]] / meanvec[[col]]
    }
    if (!is.null(criteria)){
      eblup_norm[[col]] <- eblup_norm[[col]] * criteria[[col]]
    }
  }
  return(eblup_norm)
}

all_zero <- function(ref, data){
  auxlength <- length(data[,-1])
  if (is.null(ref)){
    auxconst <- names(data[,-1])
  }else{
    numcol <- which(names(data) == ref)
    auxconst <- names(data[,-numcol])
  }

  relbase <- c( rep(">=", auxlength))
  rhsbase <- c( rep(0, auxlength))
  ctr <- data.frame(a = auxconst, b = relbase, c = rhsbase)

  return(ctr)
}
