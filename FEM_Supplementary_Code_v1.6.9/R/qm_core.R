# Shared by nested validation and final production mapping (SI S3.2).
fem_fit_qm <- function(oob, observed, response, seed=42L, n_knots=1000L) {
  if (!response %in% c("Rich_tree", "Shannon_wiener")) stop("Unknown QM response")
  if (length(oob) != length(observed)) stop("QM observation/prediction length mismatch")
  ok <- is.finite(oob) & is.finite(observed)
  x <- as.numeric(oob[ok]); y <- as.numeric(observed[ok])
  if (length(y) < 10L || any(y < 0)) stop("QM requires >=10 finite nonnegative observations")
  if (response == "Rich_tree") {
    # A fixed local RNG must not change later fold assignment or model RNG.
    had_seed <- exists(".Random.seed", envir=.GlobalEnv, inherits=FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir=.GlobalEnv)
    on.exit(if (had_seed) assign(".Random.seed",old_seed,envir=.GlobalEnv)
            else if (exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
              rm(".Random.seed",envir=.GlobalEnv), add=TRUE)
    set.seed(seed); y <- pmax(0, y + runif(length(y), -.5, .5))
  }
  probs <- seq(0,1,length.out=n_knots)
  xq <- as.numeric(quantile(x,probs,names=FALSE))
  yq <- as.numeric(quantile(y,probs,names=FALSE))
  tab <- aggregate(yq,list(x=xq),mean); names(tab) <- c("x","y")
  if (nrow(tab)<2L) stop("Degenerate OOB distribution cannot define QM")
  list(x_knots=tab$x,y_knots=cummax(tab$y),ymin=0,
       ymax=if(response=="Rich_tree")150 else 4.5,response=response,
       n_knots=n_knots,seed=seed,continuity_correction=response=="Rich_tree",
       method="empirical_QM_shared_v1.6.1")
}
fem_apply_qm <- function(q,p) {
  ans <- pmin(q$ymax,pmax(q$ymin,stats::approx(q$x_knots,q$y_knots,
                   xout=p,rule=2,ties="ordered")$y))
  ans[!is.finite(p)] <- NA_real_
  ans
}
