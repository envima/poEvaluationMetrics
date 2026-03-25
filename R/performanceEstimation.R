#' @title SDM performance evaluation
#'
#' @description
#' Evaluates the performance of species distribution models (SDMs).
#' The evaluation supports three data types: presence-absence (PA), presence-artificial-absence (PAA), and presence-background (PBG).
#' Metrics include predictive accuracy (AUC, COR), predictive error measures (MAE, BIAS), across cross-validation replicates.
#'
#' @param prediction A `terra::SpatRaster` object with the prediction map.
#' @param presence **Required.** An `sf` object of known presence locations.
#' @param absence Optional. An `sf` object of known absence locations, or `FALSE` to skip PA metrics.
#' @param background Optional. An `sf` object of background points, or logical `TRUE` to auto-generate, or `FALSE` to skip PBG.
#' @param aa Optional. An `sf` object of artificial absence points, or logical `TRUE` to derive via AOA, or `FALSE` to skip PAA.
#' @param environmentalVariables A `terra::SpatRaster` of covariates when `background` or `aa` is `TRUE`, or `sf` points are not supplied for those; may be `NA` only if both are `FALSE` or only fixed `sf` points are used.
#' @param noPointsTesting Integer. Number of background or artificial absence points to generate when sampling.
#' @param replicates Integer. Number of replicate draws for PBG/PAA; column-wise
#'   means across replicates are returned for those scenarios.
#'
#' @return A `data.frame` with one row per computed scenario. The first column is
#'   `scenario` (`"PA"`, `"PAA"`, or `"PBG"`). Remaining columns match
#'   [calculateMetrics()] / [indexCalculation()].
#'
#' @details
#' Metrics combine discrimination (e.g. AUC, COR), probability calibration errors
#' (MAE, BIAS), and threshold-based skill from [mecofun::evalSDM()] and helper
#' indices ([Fbp()], [sedi()], [orss()], [omission()]). PAA uses CAST `aoa()` to
#' mask environmentally dissimilar areas before sampling artificial absences.
#'
#' @examples
#' \dontrun{
#'   result <- performanceEstimation(
#'     prediction = prediction_raster,
#'     presence = presence_points,
#'     background = TRUE,
#'     environmentalVariables = env_rasters,
#'     replicates = 50
#'   )
#'   result[result$scenario == "PBG", ]
#' }
#'
#' @export
performanceEstimation <- function(
    prediction,
    presence,
    absence = FALSE,
    background = TRUE,
    aa = TRUE,
    environmentalVariables = NA,
    noPointsTesting = NA,
    replicates = 100
) {

  # -------------------------------------------------------------------
  # Input validation
  # -------------------------------------------------------------------
  if (!inherits(prediction, "SpatRaster")) stop("'prediction' must be an spatRaster object.")
  if (!inherits(presence, "sf")) stop("'presence' must be an sf object.")
  if (!(inherits(absence, "sf") || isFALSE(absence))) stop("'absence' must be an sf object or FALSE")
  if (!(inherits(background, "sf") || is.logical(background))) stop("'background' must be an sf object or a logical (TRUE/FALSE)")
  if (!(inherits(aa, "sf") || is.logical(aa))) stop("'aa' must be an sf object or a logical (TRUE/FALSE)")
  if (!(inherits(environmentalVariables, "SpatRaster") || is.na(environmentalVariables))) stop("'environmentalVariables' must be either a terra::SpatRaster object or NA.")
  if (!(is.numeric(noPointsTesting) || is.na(noPointsTesting))) stop("'noPointsTesting' must be a numeric value or NA")

  if (inherits(absence, "sf") && nrow(absence) < 1) {
    absence <- FALSE
    message("Number of observations in absence is < 1. Set absence to FALSE. No metric on presence-absence data calculated.")
  }

  if (isFALSE(absence)[1] && isFALSE(background)[1] && isFALSE(aa)[1]) {
    stop("At least one of absence, background, or artificial absence (aa) must be provided.")
  }

  if ((isTRUE(background) || isTRUE(aa)) && !inherits(environmentalVariables, "SpatRaster")) {
    stop("Environmental variables must be provided to generate background or artificial absence data.")
  }

  if (is.na(noPointsTesting)) {
    noPointsTesting <- nrow(presence)
  }

  # -------------------------------------------------------------------
  # 1. Presence-Background (PBG)
  # -------------------------------------------------------------------
  if (isTRUE(background)) {
    message(paste("Calculating metrics on presence-background with", replicates, "replicates."))
    indexPBG <- do.call("rbind", lapply(1:replicates, function(i) {
      bg <- generateBackgroundPoints(environmentalVariables, noPointsTesting)
      calculateMetrics(prediction, presence, bg)
    }))
    indexPBG <- as.data.frame(lapply(indexPBG, mean, na.rm = TRUE))
  } else {
    indexPBG <- NA
  }

  # -------------------------------------------------------------------
  # 2. Presence-Artificial-Absence (PAA)
  # -------------------------------------------------------------------
  if (is.logical(aa) && isTRUE(aa)) {
    message(paste("Calculating metrics on presence-artificial-absence with", replicates, "replicates."))

    extr <- terra::extract(environmentalVariables, presence, ID = FALSE)
    aoa_result <- suppressMessages(CAST::aoa(newdata = environmentalVariables, train = extr, variables = "all", verbose = FALSE))
    aa_mask <- aoa_result$AOA
    aa_mask[aa_mask > 0] <- NA

    indexPAA <- do.call("rbind", lapply(1:replicates, function(i) {
      aa <- generateAAPoints(aa_mask, noPointsTesting)
      calculateMetrics(prediction, presence, aa)
    }))
    indexPAA <- as.data.frame(lapply(indexPAA, mean, na.rm = TRUE))
  } else {
    indexPAA <- NA
  }

  # -------------------------------------------------------------------
  # 3. Presence-Absence (PA)
  # -------------------------------------------------------------------
  if (!is.logical(absence) || !isFALSE(absence)) {
    indexPA <- calculateMetrics(prediction, presence, absence)
  } else {
    indexPA <- NA
  }

  # -------------------------------------------------------------------
  # Combine as one data.frame; first column identifies scenario
  # -------------------------------------------------------------------
  append_scenario <- function(df, name) {
    if (identical(df, NA)) return(NULL)
    cbind(scenario = name, df, row.names = NULL, stringsAsFactors = FALSE)
  }

  parts <- list(
    append_scenario(indexPA, "PA"),
    append_scenario(indexPAA, "PAA"),
    append_scenario(indexPBG, "PBG")
  )
  parts <- parts[!vapply(parts, is.null, logical(1))]

  if (!length(parts)) {
    stop("No evaluation scenarios produced results.")
  }

  out <- do.call(rbind, parts)
  rownames(out) <- NULL

  gc()
  return(out)
}
