# Gold-standard validation (DeBruine & Barr, 2021): simulate one large data set from
# known ground-truth parameters, fit the maximal mixed model, and confirm it recovers the
# fixed effect, the by-subject and by-item random-effect SDs and correlations, and the
# residual SD that were specified in the design.

args <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
source(file.path(here, "..", "R", "core.R"))
source(file.path(here, "..", "R", "simulate.R"))
suppressPackageStartupMessages(library(lmerTest))

spec <- load_spec(file.path(here, "..", "..", "..", "spec", "examples", "crossed_mixed_rt.json"))
spec$units$subject$n <- 80   # enlarge so recovery is tight
spec$units$item$n <- 60
shift <- spec$response$shift

d <- simulate_design(spec)
d$cond_num <- ifelse(d$condition == "unrelated", 0.5, -0.5)   # the specified contrast
d$y <- log(d$RT - shift)                                      # shifted-lognormal model scale

fit <- suppressWarnings(suppressMessages(
  lmer(y ~ cond_num + (1 + cond_num | subject) + (1 + cond_num | item), data = d)))

fe <- fixef(fit); vc <- as.data.frame(VarCorr(fit))
getsd  <- function(grp, v1) {
  if (is.na(v1)) return(vc$sdcor[vc$grp == grp])
  vc$sdcor[vc$grp == grp & vc$var1 == v1 & is.na(vc$var2)]
}
getcor <- function(grp) vc$sdcor[vc$grp == grp & !is.na(vc$var2)]

cat(sprintf("N = %d rows (80 subjects x 60 items x 2 conditions)\n\n", nrow(d)))
cat("Parameter                     specified   recovered\n")
cat("-----------------------------------------------------\n")
cat(sprintf("fixed intercept              %8.3f   %8.3f\n", spec$fixed$intercept, fe[["(Intercept)"]]))
cat(sprintf("fixed effect (condition)     %8.3f   %8.3f\n", spec$fixed$coefficients$cond, fe[["cond_num"]]))
cat(sprintf("by-subject intercept SD      %8.3f   %8.3f\n", spec$random$subject$intercept_sd, getsd("subject", "(Intercept)")))
cat(sprintf("by-subject slope SD          %8.3f   %8.3f\n", spec$random$subject$slopes$cond, getsd("subject", "cond_num")))
cat(sprintf("by-subject corr(int,slope)   %8.3f   %8.3f\n", spec$random$subject$correlations[["intercept,cond"]], getcor("subject")))
cat(sprintf("by-item intercept SD         %8.3f   %8.3f\n", spec$random$item$intercept_sd, getsd("item", "(Intercept)")))
cat(sprintf("by-item slope SD             %8.3f   %8.3f\n", spec$random$item$slopes$cond, getsd("item", "cond_num")))
cat(sprintf("by-item corr(int,slope)      %8.3f   %8.3f\n", spec$random$item$correlations[["intercept,cond"]], getcor("item")))
cat(sprintf("residual SD                  %8.3f   %8.3f\n", spec$response$sigma, getsd("Residual", NA)))
