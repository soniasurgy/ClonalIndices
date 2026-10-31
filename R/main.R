#' Indices Calculation
#'
#' This function implements the following genotypic selection indices: Balanced selection index (BSI), Multiplicative Balanced Selection Index (MBSI), Adapted Smith Selection Index (SSI), Adapted Elston Selection Index (ESI) and Adapted Mulamba & Mock Selection Index (MSI).
#'
#' Two indices, BSI and MBSI, were newly developed for clonal selection, while SSI, ESI and MSI are adaptations of existing selection indices for use with genotypic data.
#'
#' For more information about these indices see 'README' file.
#'
#' @details
#' The BSI and MBSI are original indices developed for the selection of superior clones based on genotypic values. The SSI, ESI and MSI are adaptations of previously established selection indices to accommodate genotypic data obtained from mixed models.
#'
#'The function allows the calculation of individual indices, as well as the comparison of clone rankings across all implemented indices.
#'
#' @param blups A vector with the names of the columns in the data corresponding to the predictors of the target traits. Needed for 'indextype' 'bsi', 'mbsi', 'ssi' and 'msi'.
#' @param preds A vector with the names of the columns in the data corresponding to the predicted values of the target traits. Needed for 'indextype' 'esi'.
#' @param ref Name of the reference column used to identify individual clones (e.g., clone ID). Defaults to the first column.
#' @param clnum An integer specifying the number of top clones to be seen in the output. If omitted, equal to the total number of clones.
#' @param meanvec A named numeric vector of trait means; if omitted, data are assumed to be already normalized by the mean. Not used for 'esi'.
#' @param criteria A named numeric vector indicating the selection criterion for each trait: 1 for traits to be increased, -1 for traits to be decreased. If omitted, all traits are assumed to be selected for increase.
#' @param indextype A string indicating the selection index to be calculated. 'bsi' for Balanced selection index, 'mbsi' for Multiplicative Balanced Selection Index, 'ssi' for the Adapted Smith selection index, 'msi' for the Adapted Mulamba & Mock selection index, or 'esi' for the Adapted Elston selection index . If omitted, 'bsi' is calculated.
#' @param data A data frame comprising the input data consisting of the EBLUPs of genotypic effects and/or the predicted values obtained by the fitting of the mixed models, which serve as the basis for the selection procedure.
#'
#' @return
#' A list with the following objects:
#'  - `ranking`  with the top 'clnum' clones in the index mentioned in 'indextype' and the corresponding index value
#'  - `rankall`  with the top 'clnum' clones in all indices and the corresponding index value
#'  - `rankglobal` with index value and rank position for each clone in all indices
#'
#' @references
#' Smith, H., 1936. A discriminant function for plant selection. Ann. Eugen, Volume 7, pp. 240-250.
#'
#' Elston, R., 1963. A weight-free index for the purpose of ranking or selection with respect to several traits at a time. Biometrics, Volume 19, pp. 85-97.
#'
#' Mulamba, N. N. and Mock, J., 1978. Improvement of yield potential of the ETO blanco maize (Zea mays L.) population by breeding for plant traits. Egyptian Journal of Genetics and Cytology, Volume 7, pp. 40-57.
#'
#' Smith, H., 1936. A discriminant function for plant selection. Ann. Eugen, Volume 7, pp. 240-250.
#'
#'
#'
#' @export
#'
#' @examples
#' myeblups <- c("ydblup","tssblup", "tablup", "bwblup", "phblup" )
#' mypreds <- c("ydpred","tsspred", "tapred", "bwpred", "phpred" )
#' mymeanvec <- c(ydblup = 4.841, tssblup = 23.16, tablup = 5.589, bwblup = 1.342, phblup = 3.682)
#' mycriteria <- c(ydblup = 1, tssblup = 1, tablup = 1, bwblup = -1, phblup = -1, ydpred = 1,
#'    tsspred = 1, tapred = 1, bwpred = -1, phpred = -1)
#' exemp <- clonind(
#'    blups = myeblups,
#'    preds = mypreds,
#'    clnum = 25,
#'    meanvec = mymeanvec,
#'    criteria = mycriteria,
#'    indextype = "esi",
#'    data = Arinto
#'    )
#' exemp
#' exemp$rankall
#' exemp$rankglobal
clonind <- function(blups = NULL, preds = NULL, ref = NULL, clnum = NULL, meanvec = NULL, criteria = NULL, indextype = c("bsi", "mbsi", "ssi", "esi", "msi"), data){


indextype <- match.arg(indextype)


  # Input checks
  if (length(blups) < 1 && indextype != "esi") stop("There must be at least one trait in selection (argument 'blups' for this 'indextype' ")
  if (length(preds) < 1 && indextype == "esi") stop("There must be at least one trait in selection (argument 'preds' for this 'indextype' ")
  if ((!is.null(blups) && !is.null(preds)) && (length(blups) != length(preds))) stop("The number of traits indicated in 'blups' must be the same than in 'preds' (if none is null) ")
  if (!is.null(ref) && length(ref) != 1) stop("Only one reference column is allowed (ref).")
  if (is.null(ref)) ref <- names(data)[1]

  if (!all(blups %in% names(data))) stop("Some columns in 'blups' do not exist in `data`.")
  if (!all(preds %in% names(data))) stop("Some columns in 'preds' do not exist in `data`.")
  if (!all(ref %in% names(data))) stop("The reference column (ref) does not exist in `data`.")

  if (is.null(clnum)) clnum <- nrow(data)
  if (clnum < 1 || clnum > nrow(data)) stop("'clnum' must be >= 0 and <= than the total of genotypes.")
  clnum <- as.integer(clnum)

  # Trait columns
  cols <- unique(c(blups, preds))

  # inicialing indexcol, indexrow

  indexcol <- length(cols) + 1
  indexrow <- nrow(data)


  # Validate vector lengths
  if (!is.null(meanvec) && length(meanvec) != length(blups)) {
    stop("Length of 'meanvec' must contain values for all traits in `blups`.")
  }

  if (!is.null(criteria) && length(criteria) != (length(blups)+length(preds))) {
    stop("Length of 'criteria' must contain values for all traits.")
  }


  # Set defaults
  if (is.null(criteria)) criteria <- stats::setNames(rep(1, length(cols)), cols)
  if (is.null(meanvec)) meanvec <- stats::setNames(rep(1, length(blups)), blups)

  # Normalize trait
  if (!is.null(blups)){
    auxeblups <- norm_eblup100(data[, blups, drop = FALSE], blups, meanvec, criteria)
    index <- data.frame(data[, ref, drop = FALSE], auxeblups)
  }

  if (!is.null(preds)){
    indexesi <- data.frame(data[, ref, drop = FALSE], data[, preds, drop = FALSE])
  }


  # Calculate indices  (index matrix)
  if (!is.null(blups)){
  # BSI
    indexbsi <- index
    indexbsi$numer <- rowSums(indexbsi[,c(blups)])
    indexbsi$denom <- rowSums(abs(indexbsi[,c(blups)]))
    indexbsi$BSI <- indexbsi$numer/indexbsi$denom
    indexbsi$BSIrank <- rank(-indexbsi$BSI, ties.method = "first")
    indexordbsi <- indexbsi[order(indexbsi$BSI, indexbsi$numer, decreasing = TRUE),, drop = FALSE]

    # Object all
    indexall <- data.frame(Position = c(1:clnum))
    indexall$Ref.BSI <- indexordbsi[1:clnum, 1]
    indexall$BSI <- indexordbsi$BSI[1:clnum]

    # object global

    indexglobal <- data.frame(Genotype = data[,1])
    indexglobal$BSI = indexbsi$BSI
    indexglobal$BSIrank = indexbsi$BSIrank


    # MBSI
    indexmbsi <- index
    mbsicols <- ncol(indexmbsi)
    indexmbsi$numer <- rep(1, indexrow)
    indexmbsi$denom <- rep(1, indexrow)

    for (i in 1:indexrow){
      for (j in 2:mbsicols){

        valor <- indexmbsi[[j]][i]

        if (valor >= 0){
          indexmbsi$numer[i] <- indexmbsi$numer[i] * (abs(valor)+1)
        }else{
          indexmbsi$denom[i] <- indexmbsi$denom[i] * (abs(valor)+1)
        }
      }
    }
    indexmbsi$MBSI <- indexmbsi$numer/indexmbsi$denom
    indexmbsi$MBSIrank <- rank(-indexmbsi$MBSI, ties.method = "first")
    indexordmbsi <- indexmbsi[order(indexmbsi$MBSI, decreasing = TRUE),, drop = FALSE]

    # Object all

    indexall$Ref.MBSI <- indexordmbsi[1:clnum, 1]
    indexall$MBSI <- indexordmbsi$MBSI[1:clnum]

    # object global

    indexglobal$MBSI = indexmbsi$MBSI
    indexglobal$MBSIrank = indexmbsi$MBSIrank



  # SSI
    indexssi <- index
    indexssi$SSI <- rowSums(indexssi[,c(blups)])
    indexssi$SSIrank <- rank(-indexssi$SSI, ties.method = "first")
    indexordssi <- indexssi[order(indexssi$SSI,  decreasing = TRUE),]

    # Object all

    indexall$Ref.SSI <- indexordssi[1:clnum, 1]
    indexall$SSI <- indexordssi$SSI[1:clnum]

    # object global

    indexglobal$SSI <- indexssi$SSI
    indexglobal$SSIrank <- indexssi$SSIrank


  # MSI
    indexmsi <- index

    for (i in blups){
      nomecol <- paste("r",i, sep = "_")
      indexmsi[nomecol] <- rank(indexmsi[,i])
    }

    indexmsi$MSI <- rowSums(indexmsi[,c((length(blups)+2):(2*length(blups)+1))])
    indexmsi$MSIrank <- rank(-indexmsi$MSI, ties.method = "first")
    indexordmsi <- indexmsi[order(indexmsi$MSI,  decreasing = TRUE),]

    indice = toupper(indextype)

    # Object all

    indexall$Ref.MSI <- indexordmsi[1:clnum, 1]
    indexall$MSI <- indexordmsi$MSI[1:clnum]

    # object global

    indexglobal$MSI <- indexmsi$MSI
    indexglobal$MSIrank <- indexmsi$MSIrank
  }

  if (!is.null(preds)){
      # ESI

      minval <- NULL
      for (i in preds){
        if (criteria[i]==(-1)){
          minval <- c(minval, max(indexesi[,i]))
        }else{
          minval <- c(minval, min(indexesi[,i]))
        }
      }

      termos <- sweep(indexesi[, -1], 2, minval, "-")
      for (i in preds){
        if (criteria[i]==(-1)){
          termos[, i]= termos[, i]*(-1)
        }
      }
      indexesi <- data.frame(data[, ref, drop = FALSE], termos)
      indexesi$ESI <- apply(termos, 1, prod)
      indexesi$ESIrank <- rank(-indexesi$ESI, ties.method = "first")
      indexordesi <- indexesi[order(indexesi$ESI,  decreasing = TRUE),]

      # Object all

      if (is.null(blups)){
        indexall <- data.frame(Position = c(1:clnum))
      }
      indexall$Ref.ESI <- indexordesi[1:clnum, 1]
      indexall$ESI <- indexordesi$ESI[1:clnum]

      # object global

      if (is.null(blups)){
        indexglobal <- data.frame(data[,1])
      }
      indexglobal$ESI <- indexesi$ESI
      indexglobal$ESIrank <- indexesi$ESIrank

    }

    # Resultados individuais
    if (indextype == "bsi"){
      indexsi <- data.frame(Ref = indexall$Ref.BSI, BSI = indexall$BSI)
    }

    if (indextype == "mbsi"){
      indexsi <- data.frame(Ref = indexall$Ref.MBSI, MBSI = indexall$MBSI)
    }

    if (indextype == "ssi"){
      indexsi <- data.frame(Ref = indexall$Ref.SSI, SSI = indexall$SSI)
    }

    if (indextype == "esi"){
      indexsi <- data.frame(Ref = indexall$Ref.ESI, ESI = indexall$ESI)
    }

    if (indextype == "msi"){
      indexsi <- data.frame(Ref = indexall$Ref.MSI, MSI = indexall$MSI)
    }

  # results

  result <- list(
    ranking = indexsi,
    rankall = indexall[,-1],
    rankglobal = indexglobal
  )

  attr(result, "indice") <- indice

  #global <- list()

  class(result) <- "index"
  #result$call <- match.call()
  return(result)

}


#"FIM"
