#' @title Risk Assessment Plot.
#' @description  The function raplots plots the Sensitivity and 1-Specificity curves against the calculated risk for the null (reference) and alternative (new) models, thus graphically displaying the IDIs for those with and without the events.  These plots can aid interpretation of the NRI and IDI metrics.
#' @param x1 Calculated probabilities (eg through a logistic regression model) of the reference (null) model.  Must be between 0 & 1
#' @param x2 Calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1.
#' @return A list is returned with the following items: 
#' \item{Null.p.sens}{The horizontal axis coordinates (calculated probablities) for the sensitivity of the reference (null) model}
#'  \item{Null.sens}{The vertical axis coordinates (sensitivities) of the reference (null) model}
#'   \item{Null.p.1spec}{The horizontal axis coordinates (calculated probablities) for 1-specificity of the reference (null) model}
#'    \item{Null.1spec}{The vertical axis coordinates (1-specificities) of the reference (null) model}
#' \item{Alt.p.sens}{The horizontal axis coordinates (calculated probablities) for the sensitivity of the new (alternative) model}
#'  \item{Alt.sens}{The vertical axis coordinates (sensitivities) of the new (alternative) model}
#'   \item{Alt.p.1spec}{The horizontal axis coordinates (calculated probablities) for 1-specificity of the new (alternative) model}
#'    \item{Alt.1spec}{The vertical axis coordinates (1-specificities) of the new (alternative) model}
#' @export
#' @examples
#'\dontrun{
#'data(data_risk)
#'outcome<-data_risk$outcome 
#'null_mod<-data_risk$ref
#'alt_mod<-data_risk$new
#'output<-raplot(null_mod,alt_mod,outcome)  
#'}
#' @section Further reading : The Risk Assessment Plot in this form was described by Pickering, J. W., & Endre, Z. H. (2012). New Metrics for Assessing Diagnostic Potential of Candidate Biomarkers. Clinical Journal of the American Society of Nephrology, 7, 1355–1364. doi:10.2215/CJN.09590911
#' 
raplot <- function(x1, x2, y) {
  
  require(pROC)
  
  roc.model1 <- roc(y, x1)
  roc.model2 <- roc(y, x2)
  
  sens.model1 <- roc.model1$sensitivities
  spec.model1 <- 1 - roc.model1$specificities
  n.model1 <- length(sens.model1)
  thresh.model1 <- roc.model1$thresholds
  thresh.model1 <- thresh.model1[c(-1,-n.model1)]
  sens.model1 <- sens.model1[c(-1,-n.model1)]
  spec.model1 <- spec.model1[c(-1,-n.model1)]
  
  sens.model2 <- roc.model2$sensitivities
  spec.model2 <- 1 - roc.model2$specificities
  n.model2 <- length(sens.model2)
  thresh.model2 <- roc.model2$thresholds
  thresh.model2[1]<-0
  thresh.model2[length(thresh.model2)]<-1
  thresh.model2 <- thresh.model2[c(-1,-n.model2)]
  sens.model2 <- sens.model2[c(-1,-n.model2)]
  spec.model2 <- spec.model2[c(-1,-n.model2)]
  
  n.model1 <- length(sens.model1)
  n.model2 <- length(sens.model2)
  
  # actual plotting
  
  plot(thresh.model1, sens.model1, xlim = c(0, 1), ylim = c(0, 1), type = "n", 
       lty = 2, lwd = 2, xlab = "Risk of Event", ylab = "", col = "black")
  box(lwd = 2)
  
  polygon(x = c(thresh.model1, thresh.model2[n.model2:1]), 
          y = c(sens.model1, sens.model2[n.model2:1]), border = NA, col = gray(0.8))
  polygon(x = c(thresh.model1, thresh.model2[n.model2:1]), 
          y = c(spec.model1, spec.model2[n.model2:1]), border = NA, col = gray(0.8))
  
  lines(thresh.model1, sens.model1, type = "l", lty = 2, lwd = 2, col = "black")
  lines(thresh.model2, sens.model2, type = "l", lty = 1, lwd = 2, col = "black")
  
  lines(thresh.model1, spec.model1, type = "l", lty = 2, lwd = 2, col = "red")
  lines(thresh.model2, spec.model2, type = "l", lty = 1, lwd = 2, col = "red")
  
  text(x = -0.15, y = 0.4, labels = "Sensitivity, ", col = "black", xpd = TRUE, srt = 90)
  text(x = -0.15, y = 0.4 + 0.175, labels = "1-Specificity", col = "red", xpd = TRUE, srt = 90)
  legend("topright", c("Event: New model", "Event: Baseline model", 
                       "No Event: New model", "No Event: Baseline model"), 
         col = c("black", "black", "red", "red"), 
         lty = c(1,2, 1, 2), lwd = 2, bty = "n")
  
  return(list("Null.p.sens"=thresh.model1,
              "Null.sens"=sens.model1,
              "Null.p.1spec"=thresh.model1,
              "Null.1spec"=spec.model1,
              "Alt.p.sens"=thresh.model2,
              "Alt.sens"=sens.model2,
              "Alt.p.1spec"=thresh.model2,
              "Alt.1spec"=spec.model2))
  
}


#' @title Reclassification metrics with calculated risk as inputs.
#' @description The function statistics.raplot calculates the reclassification metrics. Used by CI.raplot.
#' @param x1 Calculated probabilities (eg through a logistic regression model) of the reference (null) model.  Must be between 0 & 1
#' @param x2 Calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1. 
#' @param s1 The savings or benefit when am event is reclassified to a higher group by the new model
#' @param s2 The benefit when a non-event is reclassified to a lower even
#' @param  t The risk threshold(s) for groups. eg t<-c(0,0.1,1) is a two group model with a threshold of 0.1 & t<-c(0,0.1,0.3,1)  is a three group model with thresholds at 0.1 and 0.3.
#' @return A matrix of metrics for use within CI.raplot
#' @export
#' @section Further reading :  A useful reference for the weighted NRI is: Van Calster, B., Vickers, A. J., Pencina, M. J., Baker, S. G., Timmerman, D., & Steyerberg, E. W. (2013). Evaluation of markers and risk prediction models: overview of relationships between NRI and decision-analytic measures. Medical Decision Making : an International Journal of the Society for Medical Decision Making, 33(4), 490–501. doi:10.1177/0272989X12470757
statistics.raplot <- function(x1, x2, y,s1,s2,t) {   
  require(Hmisc)
  require(pROC)
  require(pracma)
  require(caret)
  
  s <- is.na(x1 + x2 + y)  #Remove rows with missing data
  if (any(s)) {
    s <- !s
    x1 <- x1[s]
    x2 <- x2[s]
    y <- y[s]
  }
  n <- length(y)
  y <- as.numeric(y)
  u <- sort(unique(y))
  if (length(u) != 2 || u[1] != 0 || u[2] != 1)
    stop("y must have two values: 0 and 1")
  r <- range(x1, x2)
  if (r[1] < 0 || r[2] > 1)
    stop("x1 and x2 must be in [0,1]")
  if(length(x1)!=length(x2))
    stop("Reference (Null) and New (Alt) model vectors must be the same length")
  incidence<-sum(y)/n
  if (missing(t)) {t<-c(0, incidence,1) }
  a <- y == 1
  b <- y == 0
  na <- sum(a)
  nb <- sum(b)
  d <- x2 - x1
  pev<-na/n
  pne<-nb/n
  # NRI 
  n.thresh<-length(t)-1
  risk.class.x1.ev<-cut2(x1[a],t)
  risk.class.x2.ev<-cut2(x2[a],t)
  thresh<-c()
  for (i in 1:(length(t)-1)){
    ifelse(i==(length(t)-1),thresh[i]<-paste("[",toString(t[i]),",",toString(t[i+1]),"]"), 
           thresh[i]<-paste("[",toString(t[i]),",",toString(t[i+1]),")"))
  }
  levels(risk.class.x1.ev)<-thresh
  levels(risk.class.x2.ev)<-thresh
  cM.ev<-confusionMatrix(risk.class.x2.ev,risk.class.x1.ev)
  pup.ev<-0  # P(up|event)
  pdown.ev<-0 # P(down|event)
  for (i in 1:(n.thresh-1)){pup.ev<-pup.ev+sum(cM.ev$table[(i+1):n.thresh,i])}
  for (i in 2:n.thresh){pdown.ev<-pdown.ev+sum(cM.ev$table[1:(i-1),i])}
  npup.ev<-pup.ev
  pup.ev<-pup.ev/na 
  ndown.ev<-pup.ev
  pdown.ev<-pdown.ev/na
  risk.class.x1.ne<-cut2(x1[b],t)
  risk.class.x2.ne<-cut2(x2[b],t)  
  levels(risk.class.x1.ne)<-thresh
  levels(risk.class.x2.ne)<-thresh
  cM.ne<-confusionMatrix(risk.class.x2.ne,risk.class.x1.ne)
  pup.ne<-0 # P(up|nonevent)
  pdown.ne<-0 # P(down|nonevent)
  for (i in 1:(n.thresh-1)){pup.ne<-pup.ev+sum(cM.ne$table[(i+1):n.thresh,i])}
  for (i in 2:n.thresh){pdown.ne<-pdown.ne+sum(cM.ne$table[1:(i-1),i])} 
  ndown.ne<-pdown.ne
  pdown.ne<-pdown.ne/nb
  npup.ne<-pup.ne
  pup.ne<-pup.ne/nb
  nri <- pup.ev - pdown.ev - (pup.ne - pdown.ne)
  se.nri <- sqrt((pup.ev + pdown.ev)/na + (pup.ne + pdown.ne)/nb)
  z.nri <- nri/se.nri
  nri.ev <- pup.ev - pdown.ev
  se.nri.ev <- sqrt((pup.ev + pdown.ev)/na)
  z.nri.ev <- nri.ev/se.nri.ev
  nri.ne <- pdown.ne - pup.ne
  se.nri.ne <- sqrt((pdown.ne + pup.ne)/nb)
  z.nri.ne <- nri.ne/se.nri.ne 
  # weighted NRI
  pup<-(npup.ev+npup.ne)/n
  pdown<-(ndown.ev+ndown.ne)/n
  pevent.up<- pup.ev * pev / pup        # P(event|up)= P(up|event).P(event)/P(up)
  pevent.dn<- pdown.ev * pev / pdown      # P(event|dn)= P(dn|event).P(event)/P(dn)
  pnonevent.up<- pup.ne * pne / pup      
  pnonevent.dn<- pdown.ne * pne / pdown
  wnri<-s1*(pevent.up*pup - pevent.dn * pdown)+s2*(pnonevent.dn*pdown-pnonevent.up*pup)
  se.wnri<-NA
  z.wnri<-NA
  # Category Free NRI calculations
  cfpup.ev <- mean(d[a] > 0)
  cfpup.ne <- mean(d[b] > 0)
  cfpdown.ev <- mean(d[a] < 0)
  cfpdown.ne <- mean(d[b] < 0)
  cfnri <- cfpup.ev - cfpdown.ev - (cfpup.ne - cfpdown.ne)
  se.cfnri <- sqrt((cfpup.ev + cfpdown.ev)/na + (cfpup.ne + cfpdown.ne)/nb)
  z.cfnri <- cfnri/se.cfnri
  cfnri.ev <- cfpup.ev - cfpdown.ev
  se.cfnri.ev <- sqrt((cfpup.ev + cfpdown.ev)/na)
  z.cfnri.ev <- cfnri.ev/se.cfnri.ev
  cfnri.ne <- cfpdown.ne - cfpup.ne
  se.cfnri.ne <- sqrt((cfpdown.ne + cfpup.ne)/nb)
  z.cfnri.ne <- cfnri.ne/se.cfnri.ne
  # IDI calculations
  improveSens <- sum(d[a])/na
  improveSpec <- -sum(d[b])/nb
  idi.ev <- mean(improveSens)
  idi.ne <- mean(improveSpec)
  idi <- idi.ev + idi.ne
  relidi <- 100*((sum(x2[a])/na - sum(x2[b])/nb)/(sum(x1[a])/na-sum(x1[b])/nb)-1) #relative IDI expressed as a percentage
  se.relidi <-NA
  z.relidi <- NA
  var.ev <- var(d[a])/na
  se.idi.ev <- sqrt(var.ev)
  z.idi.ev <- idi.ev/se.idi.ev
  var.ne <- var(d[b])/nb
  se.idi.ne <- sqrt(var.ne)
  z.idi.ne <- idi.ne/se.idi.ne
  se.idi <- sqrt(var.ev + var.ne)
  z.idi <- idi/se.idi
  # AUC calculations
  roc.x1 <- roc(y, x1)
  auc.x1 <- auc(roc.x1)
  ci.auc.x1 <- ci.auc(roc.x1)
  se.auc.x1 <- (ci.auc.x1[3] - auc.x1)/qnorm(0.975)
  roc.x2 <- roc(y, x2)
  auc.x2 <- auc(roc.x2)
  ci.auc.x2 <- ci.auc(roc.x2)
  se.auc.x2 <- (ci.auc.x2[3] - auc.x2)/qnorm(0.975)
  roc.test.x1.x2 <- roc.test(roc.x1, roc.x2)  #Uses the default Delong method
  sens.x1 <- roc.x1$sensitivities
  spec.x1 <- 1 - roc.x1$specificities
  n.x1 <- length(sens.x1)
  x1 <- roc.x1$thresholds
  x1 <- x1[c(-1,-n.x1)]
  sens.x1 <- sens.x1[c(-1,-n.x1)]
  spec.x1 <- spec.x1[c(-1,-n.x1)]
  sens.x2 <- roc.x2$sensitivities
  spec.x2 <- 1 - roc.x2$specificities
  n.x2 <- length(sens.x2)
  x2 <- roc.x2$thresholds
  x2 <- x2[c(-1,-n.x2)]
  sens.x2 <- sens.x2[c(-1,-n.x2)]
  spec.x2 <- spec.x2[c(-1,-n.x2)]
  # Integrated sensitivity & 1-specificity calculations: Note a 1 and 0 are added to the beginning and end of the sens and spec, and a 0 & 1 to the risks, so that it is the area from 0 to 1 and not just partial.
  is.x1 <- trapz(x = c(0,x1,1), y = c(1,sens.x1,0))  # area under curves (relates to integrated sens, 1-spec)
  is.x2 <- trapz(x = c(0,x2,1), y = c(1,sens.x2,0))
  ip.x1 <- trapz(x = c(0,x1,1), y = c(1,spec.x1,0))
  ip.x2 <- trapz(x = c(0,x2,1), y = c(1,spec.x2,0))
  
  # Output
  output <- c(n, na, nb, pup.ev, pup.ne, pdown.ev, pdown.ne, nri, se.nri, z.nri,
              nri.ev, se.nri.ev, z.nri.ev, nri.ne, se.nri.ne, z.nri.ne, 
              cfpup.ev, cfpup.ne, cfpdown.ev, cfpdown.ne, cfnri, se.cfnri, z.cfnri,
              cfnri.ev, se.cfnri.ev, z.cfnri.ev, cfnri.ne, se.cfnri.ne, z.cfnri.ne, 
              improveSens, improveSpec, 
              s1,s2,wnri,se.wnri,z.wnri,
              idi.ev, se.idi.ev, z.idi.ev, idi.ne, 
              se.idi.ne, z.idi.ne, idi, se.idi, z.idi,relidi, se.relidi, z.relidi, is.x1, NA, is.x2, NA, 
              ip.x1, NA, ip.x2, NA, auc.x1, se.auc.x1, auc.x2, se.auc.x2, 
              roc.test.x1.x2$p.value,incidence)
  names(output) <- c("n", "na", "nb", "pup.ev", "pup.ne", "pdown.ev", "pdown.ne", 
                     "nri", "se.nri", "z.nri", "nri.ev", "se.nri.ev", "z.nri.ev",
                     "nri.ne", "se.nri.ne", "z.nri.ne",
                     "cfpup.ev", "cfpup.ne", "cfpdown.ev", "cfpdown.ne", 
                     "cfnri", "se.cfnri", "z.cfnri", "cfnri.ev", "se.cfnri.ev", "z.cfnri.ev",
                     "cfnri.ne", "se.cfnri.ne", "z.cfnri.ne", "improveSens", "improveSpec",
                     "s1","s2","wnri","se.wnri","z.wnri",
                     "idi.ev", "se.idi.ev", "z.idi.ev", "idi.ne", "se.idi.ne", 
                     "z.idi.ne", "idi", "se.idi", "z.idi","relidi", "se.relidi", "z.relidi", 
                     "is.x1", "se.is.x1",
                     "is.x2", "se.is.x2", "ip.x1", "se.ip.x1", "ip.x2", "se.ip.x2", 
                     "auc.x1", "se.auc.x1", "auc.x2", "se.auc.x2", 
                     "roc.test.x1.x2.pvalue","incidence")
  return(output)
}

#' @title Reclassification metrics and confidence intervals with calculated risk as inputs .
#' @description The CI.raplot function produces summary metrics for risk assessment. Outputs the NRI, IDI, weighted NRI and category Free NRI all for those with events and those without events.  Also the AUCs of the two models and the comparison (DeLong) between AUCs. Output includes confidence intervals. Uses statistics.raplot.  Displayed graphically by raplot.
#' @param x1 Calculated probabilities (eg through a logistic regression model) of the reference (null) model.  Must be between 0 & 1
#' @param x2 Calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1. 
#' @param s1 The savings or benefit when am event is reclassified to a higher group by the new model
#' @param s2 The benefit when a non-event is reclassified to a lower group
#' @param  t The risk threshold(s) for groups. eg t<-c(0,0.1,1) is a two group model with a threshold of 0.1 & t<-c(0,0.1,0.3,1)  is a three group model with thresholds at 0.1 and 0.3.
#' @param cis The confidence interval method to use "asymptotic" or "boot".  Normally "boot" is used
#' @param conf.level The confidence interval expressed as a fraction of 1 (ie 0.95 is the 95\% confidence interval )
#' @param n.boot  The number of "bootstraps" to use. Performance slows down with more bootstraps. For trialling result, use a low number (eg 2), for accuracy use a large number (eg 2000)
#' @return A matrix is returned with the following items: 
#' \item{Total (n)}{Total number of subjects}
#'  \item{Events (n)}{Number of subjects with the event (outcome) of intrest}
#'   \item{Non-events (n)}{Number of subjects without the event (outcome) of intrest}
#'    \item{cfNRI events}{The category free (continuous) NRI (Net Reclassification Improvement) with confidence interval for those with the event. While expressed as a percentage this is between -200 and +200}
#' \item{cfNRI non-events}{The category free (continuous) NRI with confidence interval for those without the event. While expressed as a percentage this is between -200 and +200}
#' \item{cfNRI}{The sum of cfNRI events and cfNRI non-events}
#' \item{NRI events}{The NRI with confidence interval for those with the event.}
#' \item{NRI non-events}{The NRI with confidence interval for those without the event.}
#' \item{NRI}{The sum of NRI events and NRI non-events. This is between -200 and +200}
#' \item{wNRI}{The weighted NRI}
#' \item{IDI events}{The IDI (Integrated Discrimination Improvement) with confidence interval for those with the event. Expressed as a fraction}
#' \item{IDI non-events}{The IDI with confidence interval for those without the event. Expressed as a fraction}
#' \item{IDI}{The sum of IDI events and IDI non-events}
#' \item{relIDI}{The relative IDI}
#' \item{IS(null model)}{The Integrated Sensitivity (area under the sensitivity-calculated risk curve) for the reference (null) model}
#' \item{IS(alt model)}{The Integrated Sensitivity for the reference (alt) model. Note, the IDI events should be the difference between IS(alt model) and IS(null model)}
#' \item{IP(null model)}{The Integrated 1-Specificity (area under the 1-specificity-calculated risk curve) for the reference (null) model}
#' \item{IP(alt model)}{The Integrated S1-Specificity for the reference (alt) model. Note, the IDI non-events should be the difference between IP(alt model) and IP(null model)}
#' \item{AUC(null model)}{The Area Under the Receiver Operator Characteristic Curve for the reference (null) model}
#' \item{AUC(alt model)}{The Area Under the Receiver Operator Characteristic Curve for the reference (alt) model}
#' \item{difference (p)}{P value for the difference in AUCs (DeLong method)}
#' \item{incidence}{The incidence of the event}
#' @export 
#'@examples
#'\dontrun{
#'data(data_risk)
#'y<-data_risk$outcome 
#'x1<-data_risk$ref
#'x2<-data_risk$new
#'t<-c(0,0.19,1) 
#'T=0.1 #e.g.
#'s1=1/T 
#'s2=1/(1-T)
#'output<-CI.raplot(x1, x2, y,s1,s2,t,cis = c("boot"), conf.level = 0.95, n.boot = 4, dp = 4) 
#'}
#' @section Further reading :  
#' A useful reference for the weighted NRI is: Van Calster, B., Vickers, A. J., Pencina, M. J., Baker, S. G., Timmerman, D., & Steyerberg, E. W. (2013). Evaluation of markers and risk prediction models: overview of relationships between NRI and decision-analytic measures. Medical Decision Making : an International Journal of the Society for Medical Decision Making, 33(4), 490–501. doi:10.1177/0272989X12470757
#' Reference for the NRI and IDI: Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added predictive ability of a new marker: From area under the ROC curve to reclassification and beyond. Statistics in Medicine, 27(2), 157–172. doi:10.1002/sim.2929
#' Reference for the continuous (category free) NRI: Pencina, M. J., D'Agostino, R. B., & Steyerberg, E. W. (2011). Extensions of net reclassification improvement calculations to measure usefulness of new biomarkers. Statistics in Medicine, 30(1), 11–21. doi:10.1002/sim.4085
#' Reference for the method of comparing AUCs: DeLong, E., DeLong, D., & Clarke-Pearson, D. (1988). Comparing the areas under 2 or more correlated receiver operating characteristic curves - a nonparametric approach. Biometrics, 44(3), 837–845.
CI.raplot <- function(x1, x2, y, s1,s2, t, cis = c("asymptotic", "boot"), conf.level = 0.95, n.boot = 2000, dp = 4) {
  if(length(x1)!=length(x2))
    stop("Reference (Null) and New (Alt) model vectors must be the same length")
  if (missing(s1)){
    s1<-0
    s2<-0
  }
  ifelse(missing(t),
         results <- statistics.raplot(x1, x2,y,s1,s2),  
         results <- statistics.raplot(x1, x2, y, s1,s2,t)
  )
  
  if (cis == "boot") {
    
    results.boot <- matrix(NA, n.boot, length(results))
    
    colnames(results.boot) <- names(results)
    
    for (i in 1:n.boot) {
      boot.index <- sample(length(y), replace = TRUE)
      risk.model1.boot <- x1[boot.index]
      risk.model2.boot <- x2[boot.index]
      cc.status.boot <- y[boot.index]           
      results.boot[i, ] <- statistics.raplot(x1 = risk.model1.boot, 
                                             x2 = risk.model2.boot, 
                                             y = cc.status.boot,
                                             s1,
                                             s2,
                                             t)
    }
    
    results.se.boot <- apply(results.boot, 2, sd) 
    results[grep("se", names(results))] <- results.se.boot[grep("se", names(results)) - 1]
  }
  
  # calculate cis and return 
  z <- abs(qnorm((1 - conf.level)/2))
  
  results.matrix <- matrix(NA, 27, 2)
  
  results.matrix[1, ] <- c("Total (n)", results["n"])
  results.matrix[2, ] <- c("Events (n)", results["na"])
  results.matrix[3, ] <- c("Non-events (n)", results["nb"])
  results.matrix[4, ] <- c("cfNRI and summary statistics","-------------------------")
  results.matrix[5, ] <- c("cfNRI events (%)", 
                           paste(round(100*results["cfnri.ev"], dp-2), " (", 
                                 round(100*results["cfnri.ev"] - z * 100*results["se.cfnri.ev"], dp-2),
                                 " to ", round(100*results["cfnri.ev"] + 
                                                 z * 100*results["se.cfnri.ev"], dp-2), ")", sep = ""))
  results.matrix[6, ] <- c("cfNRI non-events (%)", 
                           paste(round(100*results["cfnri.ne"], dp-2), " (",
                                 round(100*results["cfnri.ne"] - z * 100*results["se.cfnri.ne"], dp-2),
                                 " to ", round(100*results["cfnri.ne"] +  z * 100*results["se.cfnri.ne"], 
                                               dp-2), ")", sep = "")) 
  results.matrix[7, ] <- c("cfNRI (dimensionless)", 
                           paste(round(100*results["cfnri"], dp-2), " (", 
                                 round(100*results["cfnri"] - z * 100*results["se.cfnri"], dp-2), 
                                 " to ", round(100*results["cfnri"] + z * 100*results["se.cfnri"], 
                                               dp-2), ")", sep = ""))
  results.matrix[8, ] <- c("NRI and summary statistics","-------------------------")
  results.matrix[9, ] <- c("NRI events (%)", 
                           paste(round(100*results["nri.ev"], dp-2), " (", 
                                 round(100*results["nri.ev"] - z * 100*results["se.nri.ev"], dp-2),
                                 " to ", round(100*results["nri.ev"] + 
                                                 z * 100*results["se.nri.ev"], dp-2), ")", sep = ""))
  results.matrix[10, ] <- c("NRI non-events (%)", 
                            paste(round(100*results["nri.ne"], dp-2), " (",
                                  round(100*results["nri.ne"] - z * 100*results["se.nri.ne"], dp-2),
                                  " to ", round(100*results["nri.ne"] +  z * 100*results["se.nri.ne"],  dp-2), ")", sep = "")) 
  results.matrix[11, ] <- c("NRI (dimensionless)", 
                            paste(round(100*results["nri"], dp-2), " (", 
                                  round(100*results["nri"] - z * 100*results["se.nri"], dp-2), 
                                  " to ", round(100*results["nri"] + z * 100*results["se.nri"], 
                                                dp-2), ")", sep = ""))
  results.matrix[12, ] <- c("Weighted NRI and summary statistics","-------------------------")
  results.matrix[13, ] <- c("wNRI (dimensionless)", 
                            paste(round(100*results["wnri"], dp-2), " (", 
                                  round(100*results["wnri"] - z * 100*results["se.wnri"], dp-2), 
                                  " to ", round(100*results["wnri"] + z * 100*results["se.wnri"], 
                                                dp-2), ")", sep = ""))
  results.matrix[14, ] <- c("IDI and summary statistics","-------------------------")
  results.matrix[15, ] <- c("IDI events", 
                            paste(round(results["idi.ev"], dp), " (", 
                                  round(results["idi.ev"] - z * results["se.idi.ev"], dp), 
                                  " to ", round(results["idi.ev"] + z * results["se.idi.ev"], 
                                                dp), ")", sep = ""))
  results.matrix[16, ] <- c("IDI non-events", 
                            paste(round(results["idi.ne"], dp), " (", 
                                  round(results["idi.ne"] - z * results["se.idi.ne"], dp), 
                                  " to ", round(results["idi.ne"] + z * results["se.idi.ne"], 
                                                dp), ")", sep = ""))
  results.matrix[17, ] <- c("IDI", 
                            paste(round(results["idi"], dp), " (", 
                                  round(results["idi"] - z * results["se.idi"], dp), 
                                  " to ", round(results["idi"] + z * results["se.idi"], 
                                                dp), ")", sep = ""))
  results.matrix[18, ] <- c("Relative IDI (%)", 
                            paste(round(results["relidi"], dp-2), " (", 
                                  round(results["relidi"] - z * results["se.relidi"], dp-2), 
                                  " to ", round(results["relidi"] + z * results["se.relidi"], 
                                                dp-2), ")", sep = ""))
  results.matrix[19, ] <- c("IS (null model)", 
                            paste(round(results["is.x1"], dp), " (", 
                                  round(results["is.x1"] - z * results["se.is.x1"], dp), 
                                  " to ", round(results["is.x1"] + z * results["se.is.x1"], 
                                                dp), ")", sep = ""))
  results.matrix[20, ] <- c("IS (alt model)", 
                            paste(round(results["is.x2"], dp), " (", 
                                  round(results["is.x2"] - z * results["se.is.x2"], dp), 
                                  " to ", round(results["is.x2"] + z * results["se.is.x2"], 
                                                dp), ")", sep = ""))
  results.matrix[21, ] <- c("IP (null model)", 
                            paste(round(results["ip.x1"], dp), " (", 
                                  round(results["ip.x1"] - z * results["se.ip.x1"], dp), 
                                  " to ", round(results["ip.x1"] + z *  results["se.ip.x1"], 
                                                dp), ")", sep = ""))
  results.matrix[22, ] <- c("IP (alt model)", 
                            paste(round(results["ip.x2"], dp), " (", 
                                  round(results["ip.x2"] - z * results["se.ip.x2"], dp), 
                                  " to ", round(results["ip.x2"] + z * results["se.ip.x2"], 
                                                dp), ")", sep = ""))
  results.matrix[23, ] <- c("AUC","-------------------------")
  results.matrix[24, ] <- c("AUC (null model)", 
                            paste(round(results["auc.x1"], dp), " (", 
                                  round(results["auc.x1"] - z * results["se.auc.x1"], dp), 
                                  " to ", round(results["auc.x1"] + z * results["se.auc.x1"], 
                                                dp), ")", sep = ""))
  results.matrix[25, ] <- c("AUC (alt model)", 
                            paste(round(results["auc.x2"], dp), " (", 
                                  round(results["auc.x2"] - z * results["se.auc.x2"], dp), 
                                  " to ", round(results["auc.x2"] +  z * results["se.auc.x2"], 
                                                dp), ")", sep = ""))
  results.matrix[26, ] <- c("difference (P)", round(results["roc.test.x1.x2.pvalue"], dp))
  results.matrix[27, ] <- c("Incidence", round(results["incidence"], dp))
  
  return(results.matrix)
}


#' @title Reclassification metrics with classes (ordinals) as inputs.
#' @description  The function statistics.classNRI calculates the NRI metrics for reclassification of data already in classes. For use by CI.classNRI.
#' @param c1 Risk class of Reference model (ordinal)
#' @param c2 Risk class of Reference model (new)
#' @param y Binary of outcome of interest. Must be 0 or 1. 
#' @param s1 The savings or benefit when am event is reclassified to a higher group by the new model
#' @param s2 The benefit when a non-event is reclassified to a lower group
#' @return A matrix of metrics for use within CI.classNRI
#' @export
statistics.classNRI <- function(c1, c2, y,s1,s2) {    
  
  s <- is.na(c1 + c2 + y)  #Remove rows with missing data
  if (any(s)) {
    s <- !s
    c1 <- c1[s]
    c2 <- c2[s]
    y <- y[s]
  }
  n <- length(y)
  y <- as.numeric(y)
  u <- sort(unique(y))
  if (length(u) != 2 || u[1] != 0 || u[2] != 1)
    stop("y must have two values: 0 and 1")
  if(length(c1)!=length(c2))
    stop("Reference (Null) and New (Alt) model vectors must be the same length")
  incidence<-sum(y)/n
  
  a <- y == 1
  b <- y == 0
  na <- sum(a)
  nb <- sum(b)
  d <- c2 - c1
  pev<-na/n
  pne<-nb/n
  # NRI 
  nup.ev<-sum(d[a] > 0)
  pup.ev <- mean(d[a] > 0)
  nup.ne<-sum(d[b] > 0)
  pup.ne <- mean(d[b] > 0)
  ndown.ev <- sum(d[a] < 0)
  pdown.ev <- mean(d[a] < 0)
  ndown.ne <- sum(d[b] < 0)
  pdown.ne <- mean(d[b] < 0)
  nri <- pup.ev - pdown.ev - (pup.ne - pdown.ne)
  se.nri <- sqrt((pup.ev + pdown.ev)/na + (pup.ne + pdown.ne)/nb)
  z.nri <- nri/se.nri
  nri.ev <- pup.ev - pdown.ev
  se.nri.ev <- sqrt((pup.ev + pdown.ev)/na)
  z.nri.ev <- nri.ev/se.nri.ev
  nri.ne <- pdown.ne - pup.ne
  se.nri.ne <- sqrt((pdown.ne + pup.ne)/nb)
  z.nri.ne <- nri.ne/se.nri.ne
  
  
  # weighted NRI
  pup<-(nup.ev+nup.ne)/n
  pdown<-(ndown.ev+ndown.ne)/n
  pevent.up<- pup.ev * pev / pup        # P(event|up)= P(up|event).P(event)/P(up)
  pevent.dn<- pdown.ev * pev / pdown      # P(event|dn)= P(dn|event).P(event)/P(dn)
  pnonevent.up<- pup.ne * pne / pup      
  pnonevent.dn<- pdown.ne * pne / pdown
  wnri<-s1*(pevent.up*pup - pevent.dn * pdown)+s2*(pnonevent.dn*pdown-pnonevent.up*pup)
  se.wnri<-NA
  z.wnri<-NA
  
  
  # Output
  output <- c(n, na, nb, pup.ev, pup.ne, pdown.ev, pdown.ne, nri, se.nri, z.nri,
              nri.ev, se.nri.ev, z.nri.ev, nri.ne, se.nri.ne, z.nri.ne,  
              s1,s2,wnri,se.wnri,z.wnri,
              incidence)
  names(output) <- c("n", "na", "nb", "pup.ev", "pup.ne", "pdown.ev", "pdown.ne", 
                     "nri", "se.nri", "z.nri", "nri.ev", "se.nri.ev", "z.nri.ev",
                     "nri.ne", "se.nri.ne", "z.nri.ne",
                     "s1","s2","wnri","se.wnri","z.wnri",
                     "incidence")
  return(output)
}


#' @title Reclassification metrics  and confidence intervals with classes (ordinals) as inputs
#' @description The function CI.classNRI calculates the NRI statistics for reclassification of data already in classes with confidence intervals.  Uses statistics.classNRI.
#' @param c1 Risk class of Reference model (ordinal)
#' @param c2 Risk class of Reference model (new)
#' @param y Binary of outcome of interest. Must be 0 or 1. 
#' @param s1 The savings or benefit when am event is reclassified to a higher group by the new model
#' @param s2 The benefit when a non-event is reclassified to a lower group
#' @param cis The confidence interval method to use "asymptotic" or "boot".  Normally "boot" is used
#' @param conf.level The confidence interval expressed as a fraction of 1 (ie 0.95 is the 95\% confidence interval )
#' @param n.boot  The number of "bootstraps" to use. Performance slows down with more bootstraps. For trialling result, use a low number (eg 2), for accuracy use a large number (eg 2000)
#' @return A matrix is returned with the following items: 
#' \item{Total (n)}{Total number of subjects}
#'  \item{Events (n)}{Number of subjects with the event (outcome) of intrest}
#'   \item{Non-events (n)}{Number of subjects without the event (outcome) of intrest}
#'    \item{NRI events}{The NRI with confidence interval for those with the event.}
#' \item{NRI non-events}{The NRI with confidence interval for those without the event.}
#' \item{NRI}{The sum of NRI events and NRI non-events. This is between -200 and +200}
#' \item{wNRI}{The weighted NRI}
#' @export
#' @examples
#'\dontrun{
#'data(data_class)
#'y<-data_class$outcome 
#'c1<-data_class$ref_class
#'c2<-data_class$new_class
#'T=0.1 # e.g.
#'s1=1/T 
#'s2=1/(1-T)
#'output<-CI.classNRI(c1, c2, y,s1,s2,cis = c("boot"), conf.level = 0.95, n.boot = 4, dp = 4)  
#'}
CI.classNRI <- function(c1, c2, y, s1,s2, cis = c("asymptotic", "boot"), conf.level = 0.95, n.boot = 2000, dp = 4) {
  if (missing(s1)){
    s1<-0
    s2<-0
  }
  if(length(c1)!=length(c2))
    stop("Reference (Null) and New (Alt) model vectors must be the same length") 
  results <- statistics.classNRI(c1, c2,y,s1,s2)  
  
  if (cis == "boot") {
    
    results.boot <- matrix(NA, n.boot, length(results))
    
    colnames(results.boot) <- names(results)
    
    for (i in 1:n.boot) {
      boot.index <- sample(length(y), replace = TRUE)
      risk.model1.boot <- c1[boot.index]
      risk.model2.boot <- c2[boot.index]
      cc.status.boot <- y[boot.index]           
      results.boot[i, ] <- statistics.classNRI(c1 = risk.model1.boot, 
                                               c2 = risk.model2.boot, 
                                               y = cc.status.boot,
                                               s1,
                                               s2)
    }
    
    results.se.boot <- apply(results.boot, 2, sd)
    
    results[grep("se", names(results))] <- results.se.boot[grep("se", names(results)) - 1]
    
  }
  
  # calculate cis and return
  
  z <- abs(qnorm((1 - conf.level)/2))
  
  results.matrix <- matrix(NA, 10, 2)
  
  results.matrix[1, ] <- c("Total (n)", results["n"])
  results.matrix[2, ] <- c("Events (n)", results["na"])
  results.matrix[3, ] <- c("Non-events (n)", results["nb"])
  results.matrix[4, ] <- c("NRI and summary statistics","-------------------------")
  results.matrix[5, ] <- c("NRI events (%)", 
                           paste(round(100*results["nri.ev"], dp-2), " (", 
                                 round(100*results["nri.ev"] - z * 100*results["se.nri.ev"], dp-2),
                                 " to ", round(100*results["nri.ev"] + 
                                                 z * 100*results["se.nri.ev"], dp-2), ")", sep = ""))
  results.matrix[6, ] <- c("NRI non-events (%)", 
                           paste(round(100*results["nri.ne"], dp-2), " (",
                                 round(100*results["nri.ne"] - z * 100*results["se.nri.ne"], dp)-2,
                                 " to ", round(100*results["nri.ne"] +  z * 100*results["se.nri.ne"], 
                                               dp-2), ")", sep = "")) 
  results.matrix[7, ] <- c("NRI (dimensionless)", 
                           paste(round(100*results["nri"], dp-2), " (", 
                                 round(100*results["nri"] - z * 100*results["se.nri"], dp-2), 
                                 " to ", round(100*results["nri"] + z * 100*results["se.nri"], 
                                               dp-2), ")", sep = ""))
  
  results.matrix[8, ] <- c("Weighted NRI and summary statistics","-------------------------")
  results.matrix[9, ] <- c("wNRI (dimensionless)", 
                           paste(round(100*results["wnri"], dp-2), " (", 
                                 round(100*results["wnri"] - z * 100*results["se.wnri"], dp-2), 
                                 " to ", round(100*results["wnri"] + z * 100*results["se.wnri"], 
                                               dp-2), ")", sep = ""))
  results.matrix[10, ] <- c("Incidence", round(results["incidence"], dp))
  
  return(results.matrix)
}
