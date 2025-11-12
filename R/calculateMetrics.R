#' @title Calculate Performance Metrics for a single replicate/model/map
#'
#' @description
#' Calculates a suite of performance metrics for a single replicate of SDM evaluation.
#' Supports presence-background (PBG), presence-artificial-absence (PAA), and presence-absence (PA) datasets.
#'
#' @param prediction A \code{terra::SpatRaster} object containing model predictions.
#' @param presence An \code{sf} object of known presence locations.
#' @param absence_or_bg_sf An \code{sf} object of absence, background, or artificial absence points.
#'
#' @return A \code{data.frame} with calculated metrics for the replicate.
#'

calculateMetrics <- function(prediction, presence, absence_or_bg_sf) {

  # Combine predicted values for presence and absence/background
  inputDF <- stats::na.omit(rbind(
    data.frame(predicted = terra::extract(prediction, presence)[[2]], observed = 1),
    data.frame(predicted = terra::extract(prediction, absence_or_bg_sf)[[2]], observed = 0)
  ))

  # Compute all metrics using indexCalculation
  metrics <- indexCalculation(inputDF, prediction = prediction)
  return(metrics)
}
