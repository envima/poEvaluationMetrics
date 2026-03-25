#' Smoothed Boyce-index style correlations (SBI)
#'
#' Six Spearman correlations between a suitability gradient and GAM-smoothed
#' presence–background responses (multiple basis types), following Liu et al.
#' (methods in Ecography; see references).
#'
#' @param prd1 Numeric vector of predicted suitability at presence sites.
#' @param prd0 Numeric vector of predicted suitability at random (e.g. background) sites.
#' @param ktry Positive integer: basis dimension passed to [mgcv::gam()] smooths
#'   (capped by the number of unique prediction values).
#'
#' @return Numeric vector of length six:
#'   `SBI_tp`, `SBI_cr`, `SBI_bs`, `SBI_ps`, `SBI_ad`, `SBI_m` (mean curve).
#'   `NA` may appear if a smoother fails (e.g. `SBI_cr`).
#'
#' @references
#' Liu, C., Newell, G., White, M., & Machunter, J. Improving the estimation of the
#' Boyce index using statistical smoothing methods for evaluating species
#' distribution models with presence-only data. *Ecography*.
#' \doi{10.1111/ecog.07218}
#'
#' @author Canran Liu, Graeme Newell, Matt White, Josephine Machunter (original implementation).
#'
#' @keywords internal
#' @noRd
sfbi <- function(prd1, prd0, ktry = 10) {

  p <- c(prd1, prd0)
  n1 <- length(prd1)
  n0 <- length(prd0)
  prd <- seq(min(p), max(p), length=n0)
  oc <- c(rep(1, n1), rep(0, n0))

  prd_cr <- tryCatch({
    md_cr = mgcv::gam(oc ~ s(p,bs="cr",k=min(ktry,length(unique(p)))), family=stats::binomial)
    prd_cr = predict(md_cr,newdata=data.frame(p=prd),type='response')
  }, error = function(e) {
    # Handle mgcvError specifically
    return(NA)
  })

  md_tp = mgcv::gam(oc ~ s(p,bs="tp",k=min(ktry,length(unique(p)))), family=stats::binomial)
  prd_tp = stats::predict(md_tp,newdata=data.frame(p=prd),type='response')
  md_bs = mgcv::gam(oc ~ s(p,bs="bs",k=min(ktry,length(unique(p)))), family=stats::binomial)
  prd_bs = stats::predict(md_bs,newdata=data.frame(p=prd),type='response')
  md_ps = mgcv::gam(oc ~ s(p,bs="ps",k=min(ktry,length(unique(p)))), family=stats::binomial)
  prd_ps = stats::predict(md_ps,newdata=data.frame(p=prd),type='response')
  md_ad = mgcv::gam(oc ~ s(p, bs = "ad",k=min(ktry,length(unique(p)))), family=stats::binomial)
  prd_ad = stats::predict(md_ad,newdata=data.frame(p=prd),type='response')
  prd_m = (prd_tp + prd_cr + prd_bs + prd_ps + prd_ad)/5
  SBI_tp <- stats::cor(prd,prd_tp,method="spearman")
  SBI_bs <- stats::cor(prd,prd_bs,method="spearman")
  SBI_ps <- stats::cor(prd,prd_ps,method="spearman")
  SBI_ad <- stats::cor(prd,prd_ad,method="spearman")

  SBI_cr <- tryCatch({

    SBI_cr <- stats::cor(prd,prd_cr,method="spearman")
    #return(SBI_cr)
  }
  , error = function(e) {return(NA)}
  )



  SBI_m <- stats::cor(prd,prd_m,method="spearman")

  return(c(SBI_tp, SBI_cr, SBI_bs, SBI_ps, SBI_ad, SBI_m))


}
