#' The Risk Assessment Plot
#'
#' The function ggrap() plots the Sensitivity and 1-Specificity curves against the calculated risk for the baseline (reference) and newmodels, thus graphically displaying the IDIs for those with and without the events.  These plots can aid interpretation of the NRI and IDI metrics.
#' 
#' @param x1 Either a logistic regression fitted using glm (base package) or lrm (rms package) or alculated probabilities (eg through a logistic regression model) of the baseline model.  Must be between 0 & 1
#' @param x2 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1 (if fitted models are provided this is extracted from the fit which for an rms fit must have x = TRUE, y = TRUE). 
#' @return a ggplot 
#' @references The Risk Assessment Plot in this form was described by Pickering, J. W., & Endre, Z. H. (2012). New Metrics for Assessing Diagnostic Potential of Candidate Biomarkers. Clinical Journal of the American Society of Nephrology, 7, 1355–1364. doi:10.2215/CJN.09590911
#' @import ggplot2
#' @import tidyr 
#' @import dplyr
#' @export
ggrap <- function(x1, x2, y) {
  
  if (class(x1)[1] == "glm") {
    y = x1$y
    x1 = stats::predict(x1, type = "response")
    data_type = "glm"
  }
  if (class(x2)[1] == "glm") {
    x2 = stats::predict(x2, type = "response")
  }
  if (class(x1)[1] == "lrm") {
    x1 = stats::predict(x1, type = "fitted")
    if (length(x1$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    y = x1$as.numeric(as.character(x1$y))
    data_type = "lrm"
  }
  if (class(x2)[1] == "lrm") {
    if (length(x2$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    x2 = stats::predict(x2, type = "fitted")
  }
  
  if (class(x1)[1] != "glm"  & class(x1)[1] != "lrm" ) {
    data_type = "User supplied"
  }
  if (length(x1) != length(x2))
    stop("Reference (baseline) and New (Alt) model vectors must be the same length")
  
  df <- data.frame(Baseline = x1, New = x2, Event = y)
  df <- df %>% 
    filter(!is.na(Baseline)) %>% 
    filter(!is.na(New)) %>% 
    filter(!is.na(Event))
  df$ID = seq(1,nrow(df),1)
  
  df_long <- df %>% 
    tidyr::pivot_longer(cols = c(Baseline,New), values_to = "Probabilities", names_to = "Model") %>% 
    group_by(Model) %>% 
    arrange(Model, Probabilities) %>% 
    group_by(Model,Probabilities) %>% 
    summarise(n_event = sum(Event == 1),
              n_nonevent = sum(Event == 0),
              IDs = as.character(ID)) %>% 
    ungroup() %>% 
    group_by(Model) %>% 
    mutate(n_event_cumulative = cumsum(n_event)) %>% 
    mutate(n_nonevent_cumulative = cumsum(n_nonevent)) %>% 
    mutate(Sensitivity = (max(n_event_cumulative) - n_event_cumulative)/max(n_event_cumulative)) %>% 
    mutate(`1-Specificity` = (max(n_nonevent_cumulative) - n_nonevent_cumulative)/max(n_nonevent_cumulative)) %>% 
    ungroup
  
  df_long_av <- df_long %>% 
    group_by(Model, Probabilities) %>% 
    summarise(Sensitivity = mean(Sensitivity),
              `1-Specificity` = mean(`1-Specificity`),
              ID_list = list(IDs)) %>% 
    ungroup()
  
  df_g <- df_long_av %>% 
    select(Probabilities, Sensitivity, `1-Specificity`, Model, ID_list) %>%  
    tidyr::pivot_longer(cols = c(Sensitivity, `1-Specificity`), values_to = "Metric", names_to = "Metric name")
  
  g <- ggplot(df_g, aes(x = Probabilities, y = Metric, colour = `Metric name`, linetype = Model )) +
    geom_line(size = 0.75) +
    scale_linetype_manual(values = c("dotted", "solid")) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1))    
  
  return(g)
}

#' The Decision curve
#' 
#' ggdecision plots decsion curves to assess the net benefit at different thresholds
#' 
#' @param x1 Either a logistic regression fitted using glm (base package) or lrm (rms package) or alculated probabilities (eg through a logistic regression model) of the baseline model.  Must be between 0 & 1
#' @param x2 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1 (if fitted models are provided this is extracted from the fit which for an rms fit must have x = TRUE, y = TRUE). 
#' @return a ggplot 
#' @references Vickers AJ, van Calster B, Steyerberg EW. A simple, step-by-step guide to interpreting decision curve analysis. Diagn Progn Res 2019;3(1):18. 2. Zhang Z, Rousson V, Lee W-C, et al. Decision curve analysis: a technical note. Ann Transl Med 2018;6(15):308–308. 
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @export
#'
ggdecision <- function(x1, x2, y) {

  if (class(x1)[1] == "glm") {
    y = x1$y
    x1 = stats::predict(x1, type = "response")
    data_type = "glm"
  }
  if (class(x2)[1] == "glm") {
    x2 = stats::predict(x2, type = "response")
  }
  if (class(x1)[1] == "lrm") {
    x1 = stats::predict(x1, type = "fitted")
    if (length(x1$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    y = x1$as.numeric(as.character(x1$y))
    data_type = "lrm"
  }
  if (class(x2)[1] == "lrm") {
    if (length(x2$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    x2 = stats::predict(x2, type = "fitted")
  }
  
  if (class(x1)[1] != "glm"  & class(x1)[1] != "lrm" ) {
    data_type = "User supplied"
  }
  if (length(x1) != length(x2))
    stop("Reference (baseline) and New (Alt) model vectors must be the same length")
  
  df <- data.frame(Baseline = x1, New = x2, Event = y)
  df <- df %>% 
    filter(!is.na(Baseline)) %>% 
    filter(!is.na(New)) %>% 
    filter(!is.na(Event))
  df$ID = seq(1,nrow(df),1)
  
  n <- nrow(df)
  n_event = sum(df$Event)
  n_nonevent = n - n_event
  incidence <- 100 * n_event/n
  
  
  # Net Benefit for Decision Curves
  benefit.1 <- df %>% 
    mutate(Prediction = Baseline) %>% 
    group_by(Prediction) %>% 
    summarise(n_ev = sum(Event == 1),
              n_nonev = sum(Event == 0)) %>% 
    ungroup() %>% 
    mutate(TP = n_event - cumsum(n_ev),
           FN = n_event - TP, 
           TN = cumsum(n_nonev),
           FP = n_nonevent - TN,
           prevalence = (TP + FN)/n,
           sens = TP/(TP + FN),
           spec = TN/(TN + FP),
           `1-spec` = 1 - spec,
           npv = TN/(TN + FN),
           ppv = TP/(TP + FP),
           treated = TP/n - (FP/n) * Prediction/(1 - Prediction),
           untreated = TN/n - (FN/n) * (1 - Prediction)/Prediction,
           overall = treated + untreated,
           all = prevalence - (1 - prevalence)  * Prediction/(1 - Prediction),
           Model = "Baseline"
    )
  
  benefit.2 <- df %>% 
    mutate(Prediction = New) %>% 
    group_by(Prediction) %>% 
    summarise(n_ev = sum(Event == 1),
              n_nonev = sum(Event == 0)) %>% 
    ungroup() %>% 
    mutate(TP = n_event - cumsum(n_ev),
           FN = n_event - TP, 
           TN = cumsum(n_nonev),
           FP = n_nonevent - TN,
           prevalence = (TP + FN)/n,
           sens = TP/(TP + FN),
           spec = TN/(TN + FP),
           `1-spec` = 1 - spec,
           npv = TN/(TN + FN),
           ppv = TP/(TP + FP),
           treated = TP/n - (FP/n) * Prediction/(1 - Prediction),
           untreated = TN/n - (FN/n) * (1 - Prediction)/Prediction,
           overall = treated + untreated,
           all = prevalence - (1 - prevalence)  * Prediction/(1 - Prediction),
           Model = "New"
    )
  
  benefit = bind_rows(benefit.1, benefit.2)
  
  benefit_all_none <- benefit %>% 
    select(Prediction, all) %>% 
    mutate(none = 0)  %>% 
    tidyr::pivot_longer(cols = c("all", "none"), names_to = "Extreme models", values_to = "extremes") %>% 
    arrange(`Extreme models`)
  
  g <- ggplot() + 
    geom_point(data = benefit, aes(x = Prediction, y = treated, colour = Model), alpha = 0.3, size = 0.5) + 
    geom_smooth(data = benefit, aes(x = Prediction, y = treated, colour = Model), se = FALSE) + 
    geom_line(data = benefit_all_none, aes(x = Prediction, y = extremes, linetype = `Extreme models`)) +
    ylab("Net benefit") + xlab("Prediction threshold") +
    coord_cartesian(ylim = c(1.5 * min(benefit$treated), 1.05 * incidence/100)) +
    NULL
  return(g)
}

#' The Calibration plot
#' 
#' ggcalibrate plots the stats::predicted events against the actual event rate
#' 
#' @param x1 Either a logistic regression fitted using glm (base package) or lrm (rms package) or alculated probabilities (eg through a logistic regression model) of the baseline model.  Must be between 0 & 1
#' @param x2 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1 (if fitted models are provided this is extracted from the fit which for an rms fit must have x = TRUE, y = TRUE). 
#' @param models One of three strings:
#' \itemize{
#' \item{\strong{both}} {plot(s) will contain both the Baseline and New models;} 
#' \item{\strong{x1}} {plot(s) will contain only the Baseline model;} 
#' \item{\strong{x2}} {plot(s) will contain only the New model.}
#' }
#' @param n_cut An integer indicating either the number of intervals of the same width, the number of intervals of the same number of subjects, or the width (as a percentage) of the intervals.  
#' @param cut_type One of three strings:
#' \itemize{
#' \item{\strong{interval}} {applies the cut_interval() function to get n_cut intervals of approximately equal width;}
#' \item{\strong{number}} {applies the cut_number() function to get n_cut intervals of approximately equal number of subjects;} 
#' \item{\strong{width}} {applies the cut_width() function to get ~ 100/n_cut intervals of n_cut width.}
#' }
#' @param include_margin TRUE for including producing a bar plot of the counts of in each of the intervals. Default is FALSE.  Note if the output is saved to my_graphs then using the library gridExtra the function grid.arrange(graphs$g, graphs$g_marg , nrow = 2, heights = c(2,1)) will produce a plot with both the calibration plot and the marginal plot.  
#' @return a list of one or two ggplots
#' @examples
#'\dontrun{
#'data(data_risk)
#'y<-data_risk$outcome 
#'x1<-data_risk$baseline
#'x2<-data_risk$new
#'#e.g.
#'output <- ggcalibrate(x1, x2, y, models = "both", n_cut = 10, cut_type = "number") 
#'}
#' @import forcats
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom pracma trapz
#' @export
ggcalibrate <- function(x1, x2, y, models = c("both","x1","x2"), n_cut = 5, cut_type = c("interval","number","width"), include_margin = FALSE) {
  
  if (class(x1)[1] == "glm") {
    y = x1$y
    x1 = stats::predict(x1, type = "response")
    data_type = "glm"
  }
  if (class(x2)[1] == "glm") {
    x2 = stats::predict(x2, type = "response")
  }
  if (class(x1)[1] == "lrm") {
    x1 = stats::predict(x1, type = "fitted")
    if (length(x1$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    y = x1$as.numeric(as.character(x1$y))
    data_type = "lrm"
  }
  if (class(x2)[1] == "lrm") {
    if (length(x2$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    x2 = stats::predict(x2, type = "fitted")
  }
  
  if (class(x1)[1] != "glm"  & class(x1)[1] != "lrm" ) {
    data_type = "User supplied"
  }
  if (length(x1) != length(x2))
    stop("Reference (baseline) and New (Alt) model vectors must be the same length")
  
  df <- data.frame(Baseline = x1, New = x2, Event = y)
  df <- df %>% 
    filter(!is.na(Baseline)) %>% 
    filter(!is.na(New)) %>% 
    filter(!is.na(Event))
  df$ID = seq(1,nrow(df),1)
  
  n <- nrow(df)
  n_event = sum(df$Event)
  n_nonevent = n - n_event
  incidence <- 100 * n_event/n
  
  
  # Calibration
  
  df_long <- df %>% 
    tidyr::pivot_longer(cols = c(Baseline,New), values_to = "prediction", names_to = "Model") %>% 
    select(Model,Event, prediction) %>% # or whatever your outcome and predictions are
    mutate(prediction = 100 * as.numeric(prediction)) %>%  # If prediciton is in the 0-1 range I prefer it in the 0-100 range
    filter(!is.na(prediction)) 
  
  if(models == "x1"){df_long <- df_long %>% filter(Model == "Baseline")}  
  if(models == "x2"){df_long <- df_long %>% filter(Model == "New")}  
  
  if(cut_type == "interval"){df_long <- df_long %>% group_by(Model) %>% mutate(pred_interval = cut_interval(prediction, n = n_cut))}
  if(cut_type == "number"){df_long <- df_long %>% group_by(Model) %>% mutate(pred_interval = cut_number(prediction, n = n_cut))}
  if(cut_type == "width"){df_long <- df_long %>% group_by(Model) %>% mutate(pred_interval = cut_width(prediction, width = n_cut))}
  
  # prepare file to use in plot - sometimes only the CIs are made for the actual events and not the prediction, but that is just laziness
  gdf <- df_long %>%
    group_by(Model, pred_interval) %>%
    summarise(n = n(),
              mn_pred = mean(prediction, na.rm = TRUE),
              lci_pred = quantile(prediction, c(0.025), na.rm = TRUE),
              uci_pred = quantile(prediction, c(0.975), na.rm = TRUE),
              n_event = sum(Event == "1")) %>%
    mutate(prop_event = 100 * n_event/n) %>%
    mutate(lci_prop_event = prop_event - 1.96 * sqrt(prop_event * (100 - prop_event)/n)) %>%
    mutate(uci_prop_event = prop_event + 1.96 * sqrt(prop_event * (100 - prop_event)/n))
  
  g <- ggplot() +
    geom_abline(slope = 1, intercept = 0, colour = "grey50", linetype = "dashed") 
  
  if(models == "both"){
    g <- g +
      geom_segment(data = gdf, aes(x = mn_pred, xend = mn_pred, y = lci_prop_event, yend = uci_prop_event, colour = Model))  +
      geom_segment(data = gdf, aes(x = lci_pred, xend = uci_pred, y = prop_event, yend = prop_event, colour = Model)) +
      geom_point(data = gdf, aes(x = mn_pred, y = prop_event, colour = Model))
    
  }
  if(models == "x1" | models == "x2"){
    g <- g +
      geom_segment(data = gdf, aes(x = mn_pred, xend = mn_pred, y = lci_prop_event, yend = uci_prop_event))  +
      geom_segment(data = gdf, aes(x = lci_pred, xend = uci_pred, y = prop_event, yend = prop_event)) +
      geom_point(data = gdf, aes(x = mn_pred, y = prop_event))
  }
  
  g <- g +
    scale_x_continuous(breaks = seq(0,100,10)) +
    scale_y_continuous(breaks = seq(0,100,10)) +
    xlab("Predicted percentage") +
    ylab("Actual percentage") +
    coord_cartesian(xlim = c(0,100), ylim = c(0,100)) +
    NULL
  
  ifelse(include_margin == TRUE,
         {
           gdf <- gdf %>% 
             group_by(Model) %>% 
             arrange(mn_pred) %>% 
             mutate(width_ci = uci_pred - lci_pred)
           width = mean(gdf$width_ci, na.rm = TRUE)
           gdf$pred_interval <- fct_inorder(gdf$pred_interval)
           
           
           if(models == "both"){ 
             g_marg <- ggplot(gdf) +
               geom_col(aes(x = mn_pred, y = n, fill = Model), width = width * 2/3, alpha = 0.6)  
           }
           
           if(models == "x1" | models == "x2"){
             g_marg <- ggplot(gdf) +
               geom_col(aes(x = mn_pred, y = n), width = width * 2/3, alpha = 0.6)  
           }
           
           g_marg <-   g_marg +   
             scale_x_continuous(breaks = seq(0,100,10)) +
             xlab("Predicted percentage") +
             ylab("Count") +
             coord_cartesian(xlim = c(0,100)) +
             NULL
           
           output <- list(g = g, g_marg = g_marg)
         },
         {
           output <- list(g = g)
         }
  )  
  
  return(output)
}

#' Statistical metrics
#' 
#' The function statistics.raplot calculates the reclassification metrics. Used by CI.raplot.
#' 
#' @param x1 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the baseline model.  Must be between 0 & 1
#' @param x2 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1 (if fitted models are provided this is extracted from the fit which for an rms fit must have x = TRUE, y = TRUE). 
#' @param  t The risk threshold(s) for groups. eg t<-c(0,0.1,1) is a two group scenario with a threshold of 0.1 & t<-c(0,0.1,0.3,1)  is a three group scenario with thresholds at 0.1 and 0.3. Nb. If no t is provided it defaults to a single threshold at the prevalence of the cohort.  
#' @return A matrix of metrics for use within CI.raplot
#' @import pROC
#' @import dplyr
#' @import tidyr
#' @export
#' 
statistics.raplot <- function(x1, x2, y,  t = NULL) {   
 
  df <- data.frame(x1 = x1, x2 = x2, event = y) # y = 1 for the event, 0 for not the event.  
  
  #Remove rows with missing data
  df <- df %>%  
    filter(!is.na(x1)) %>% 
    filter(!is.na(x2)) %>% 
    filter(!is.na(event))
  
  n <- nrow(df)
  
  df <- df %>%
    mutate(event == as.numeric(event)) # just in case
  
  u <- sort(unique(df$event))
  if (length(u) != 2 || u[1] != 0 || u[2] != 1)
    stop("Event must have two values: 0 and 1")
  
  r <- range(df$x1, df$x2)
  if (r[1] < 0 || r[2] > 1)
    stop("x1 and x2 must be in [0,1]")
  
  n_event = sum(df$event)
  incidence <- 100*n_event/n
  
  if (is.null(t)) {t <- c(0, incidence/100,1) }
  
  # ROC
  roc.1 <- roc(df$event, df$x1)
  roc.2 <- roc(df$event, df$x2)
  auc.x1 <- as.numeric(auc(roc.1))
  auc.x2 <- as.numeric(auc(roc.2))
  auc.difference = auc.x2 - auc.x1
  
  df <- df %>% 
    mutate(risk.class.x1 = Hmisc::cut2(x1,t)) %>% #risk groups based on thresholds
    mutate(risk.class.x2 = Hmisc::cut2(x2,t)) %>% 
    mutate(difference = x2 - x1)  %>%   # difference in probabilities
    mutate(mse_x1 = (event - x1)^2) %>% # mean squared errors of the baseline model
    mutate(mse_x2 = (event - x2)^2)  
  
  n_event = sum(df$event)
  n_non_event = nrow(df) - n_event
  n = n_event + n_non_event
  prevalence = n_event/n
  
  #NRI
  df_NRI <- df %>% 
    group_by(event) %>% 
    mutate(up = ifelse(as.numeric(risk.class.x2) > as.numeric(risk.class.x1),1,0))  %>% 
    mutate(down = ifelse(as.numeric(risk.class.x2) < as.numeric(risk.class.x1),1,0)) %>% 
    dplyr::summarise(
      n = n(),
      n_up = sum(up),
      n_down = sum(down),
      nri = (n_up - n_down)/n # for the non-event take the negative of this (below)
    ) %>% 
    ungroup() %>% 
    mutate(nri = ifelse(event == 0, -nri, nri))
  
  #IDI
  df_IDI <- df %>% 
    group_by(event) %>%  # only report for those with and without the disease seperately
    dplyr::summarise(IDI = mean(difference, na.rm = TRUE)) %>% # the IDI is the mean change in Prediction
    ungroup() %>% 
    mutate(IDI = ifelse(event == 0, -IDI, IDI)) # for the non-event convert this so that a positive mean change means the move in the right direction
  
  #Brier
  df_Brier <- df %>% 
    dplyr::summarise(
      brier_baseline = mean(mse_x1),
      brier_new = mean(mse_x2),
      brier_skill = 100 * (brier_baseline - brier_new)/brier_baseline
    )
  
  df_model.1 = data.frame(sens = roc.1$sensitivities, spec = roc.1$specificities, Prediction = roc.1$thresholds)
  df_model.1 <- df_model.1 %>% mutate(Model = "Baseline")
  
  df_model.2 = data.frame(sens = roc.2$sensitivities, spec = roc.2$specificities, Prediction = roc.2$thresholds)
  df_model.2 <- df_model.2 %>% mutate(Model = "New")
  
  df_model = bind_rows(df_model.1, df_model.2)
  df_model <- df_model %>% 
    mutate(Prediction = ifelse(sens == 1 & spec == 0, 0, ifelse(sens == 0 & spec == 1, 1, Prediction) ))
  
  df_model_Istats <- df_model %>% 
    group_by(Model) %>% 
    arrange(Prediction) %>% 
    dplyr::summarise(
      IS = pracma::trapz(x = c(Prediction), y = c(sens)),  
      IP = pracma::trapz(x = c(Prediction), y = 1 - spec) )
  
  # Output
  output <- list(n = n, 
                 n_event = n_event, 
                 n_non_event = n_non_event, 
                 Prevalence = prevalence,
                 NRI_up_event = as.integer(df_NRI[df_NRI$event == 1,]$n_up),
                 NRI_up_nonevent = as.integer(df_NRI[df_NRI$event == 0,]$n_up),
                 NRI_down_event = as.integer(df_NRI[df_NRI$event == 1,]$n_down),
                 NRI_down_nonevent = as.integer(df_NRI[df_NRI$event == 0,]$n_down),
                 NRI_event = df_NRI[df_NRI$event == 1,]$nri,
                 NRI_nonevent = df_NRI[df_NRI$event == 0,]$nri,
                 IDI_event = df_IDI[df_IDI$event == 1,]$IDI,
                 IDI_nonevent = df_IDI[df_IDI$event == 0,]$IDI,
                 IP_baseline = df_model_Istats[df_model_Istats$Model == "Baseline",]$IP,
                 IS_baseline = df_model_Istats[df_model_Istats$Model == "Baseline",]$IS,
                 IP_new = df_model_Istats[df_model_Istats$Model == "New",]$IP, 
                 IS_new = df_model_Istats[df_model_Istats$Model == "New",]$IS,
                 Brier_baseline = df_Brier$brier_baseline,
                 Brier_new = df_Brier$brier_new,
                 Brier_skill = df_Brier$brier_skill,
                 AUC_baseline = auc.x1,
                 AUC_new = auc.x2,
                 AUC_difference = auc.difference
  )
  return(output)
}


#' Extract confidence interval
#' 
#' Extract a confidence in interval from the bootstrapped results. Used by CI.raplot
#' 
#' @param results.boot The matrix of n.boot metrics from within CI.raplot
#' @param conf.level The confidence interval expressed between 0 & 1 (eg 95\%CI is conf.level = 0.95)
#' @param n.boot The number of bootstrapped samples
#' @param dp the number of decimal places to report the point estimate and confidence interval
#' @return A two column matrix with the metric name and statistic with a confidence interval
#' @import tidyr
#' @import dplyr
#' @export
extractCI <- function(results.boot, conf.level, n.boot, dp){
  
  n_vars = length(results.boot[[1]])
  results.df <- data.frame(matrix(nrow = n.boot, ncol = n_vars))   
  
  for (i in 1:n.boot) {
    temp_df <- as.data.frame(results.boot[[i]])
    results.df[i,] <- temp_df[1,]
  } 
  
  results.matrix_est <- results.df %>% 
    summarise(across(where(is.numeric), ~round(quantile(.x,c(0.5)),dp)))
  results.matrix_lower_CI <- results.df %>% 
    summarise(across(where(is.numeric), ~round(quantile(.x,c((1 - conf.level)/2)),dp)))
  results.matrix_upper_CI <- results.df %>% 
    summarise(across(where(is.numeric), ~round(quantile(.x,c(1 - (1-conf.level)/2)),dp)))
  
  temp <- bind_rows(results.matrix_est, results.matrix_lower_CI , results.matrix_upper_CI)
  results.matrix <- t(temp)
  colnames(results.matrix) <- c("V1", "V2", "V3")

  results.matrix <- as_tibble(results.matrix)

  results.matrix$metric <- names(results.boot[[1]]) 
  
  results.matrix <- results.matrix %>% 
    mutate(statistics = paste0( V1, " (CI: ",V2, " to ", V3,")" )) %>% 
    select(metric, statistics)
  
  return(results.matrix)
}                

#' Extract NRI confidence intervals
#' 
#' Extract a confidence in interval from the bootstrapped results. Used by CI.NRI
#' 
#' @param results.boot The matrix of n.boot metrics from within CI.NRI
#' @param conf.level The confidence interval expressed between 0 & 1 (eg 95\%CI is conf.level = 0.95)
#' @param n.boot The number of bootstrapped samples
#' @param dp the number of decimal places to report the point estimate and confidence interval
#' @return A two column matrix with the metric name and statistic with a confidence interval
#' @import tidyr
#' @import dplyr
#' @export
extract_NRI_CI <- function(results.boot, conf.level, n.boot, dp){
  
  n_vars = sum(unlist(lapply( results.boot[[1]], class)) == "numeric" | unlist(lapply( results.boot[[1]], class)) == "integer" ) #because we don't calculate CIs for the confusion matrices
  results.df <- data.frame(matrix(nrow = n.boot, ncol = n_vars))   
  
  for (i in 1:n.boot) {
    temp <- results.boot[[i]]
    temp <- unlist(temp[1:n_vars])
    results.df[i,] <- temp 
  } 
  
  results.matrix_est <- results.df %>% 
    summarise(across(where(is.numeric), ~round(quantile(.x,c(0.5)),dp)))
  results.matrix_lower_CI <- results.df %>% 
    summarise(across(where(is.numeric), ~round(quantile(.x,c((1-conf.level)/2)),dp)))
  results.matrix_upper_CI <- results.df %>% 
    summarise(across(where(is.numeric), ~round(quantile(.x,c(1 - (1-conf.level)/2)),dp)))
  
  temp <- bind_rows(results.matrix_est, results.matrix_lower_CI , results.matrix_upper_CI)
  results.matrix <- t(temp)
  colnames(results.matrix) <- c("V1", "V2", "V3")
  
  results.matrix <- as_tibble(results.matrix)
  results.matrix$metric <- names(results.boot[[1]])[1:10]

  results.matrix <- results.matrix %>% 
    mutate(statistics = paste0( V1, " (CI: ",V2, " to ", V3,")" )) %>% 
    select(metric, statistics)
  
  return(results.matrix)
} 


#' Statistical metrics with confidence intervals
#' 
#' The CI.raplot function produces summary metrics for risk assessment. Outputs the NRI, IDI, weighted NRI and category Free NRI all for those with events and those without events.  Also the AUCs of the two models and the comparison (DeLong) between AUCs. Output includes confidence intervals. Uses statistics.raplot.  Displayed graphically by raplot.
#' 
#' @param x1 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the baseline model.  Must be between 0 & 1
#' @param x2 Either a logistic regression fitted using glm (base package) or lrm (rms package) or calculated probabilities (eg through a logistic regression model) of the new (alternative) model.   Must be between 0 & 1
#' @param y Binary of outcome of interest. Must be 0 or 1 (if fitted models are provided this is extracted from the fit which for an rms fit must have x = TRUE, y = TRUE). 
#' @param t The risk threshold(s) for groups. eg t<-c(0,0.1,1) is a two group model with a threshold of 0.1 & t<-c(0,0.1,0.3,1)  is a three group model with thresholds at 0.1 and 0.3.
#' @param conf.level The confidence interval expressed as a fraction of 1 (ie 0.95 is the 95\% confidence interval )
#' @param n.boot  The number of "bootstraps" to use. Performance slows down with more bootstraps. For trialling result, use a low number (eg 5), for accuracy use a large number (eg 2000)
#' @param dp The number of decimal places to display
#' @return A list with four items: 
#' \itemize{
#' \item{\strong{1. meta_data}} {Some overall meta data - Confidence Interval, number of bootstraps, thresholds, input type}
#' \item{\strong{2. Metrics}} {Point estimates of the statistical metrics (see list below)}
#' \item{\strong{3. Each_bootstrap_metrics}} {Point estimates of the statistical metrics for each bootstrapped sample (see list below)}
#' \item{\strong{4. Summary Metrics}} {Point estimates with confidence intervals of the statistical metrics. See following list:)}
#' \itemize{
#' \item{\strong{Total (n)}} {Total number of subjects}
#' \item{\strong{Events (n)}} {Number of subjects with the event (outcome) of intrest}
#' \item{\strong{Non-events (n)}} {Number of subjects without the event (outcome) of intrest}
#' \item{\strong{NRI events}} {The NRI with confidence interval for those with the event.}
#' \item{\strong{NRI non-events}} {The NRI with confidence interval for those without the event.}
#' \item{\strong{IDI events}} {The IDI (Integrated Discrimination Improvement) with confidence interval for those with the event. Expressed as a fraction}
#' \item{\strong{IDI non-events}} {The IDI with confidence interval for those without the event. Expressed as a fraction}
#' \item{\strong{IS(baseline model)}} {The Integrated Sensitivity (area under the sensitivity-calculated risk curve) for the baseline model}
#' \item{\strong{IS(new model)}} {The Integrated Sensitivity for the reference (alt) model. Note, the IDI events should be the difference between IS(new model) and IS(baseline model)}
#' \item{\strong{IP(baseline model)}} {The Integrated 1-Specificity (area under the 1-specificity-calculated risk curve) for the baseline model}
#' \item{\strong{IP(new model)}} {The Integrated S1-Specificity for the reference (alt) model. Note, the IDI non-events should be the difference between IP(new model) and IP(baseline model)}
#' \item{\strong{AUC(baseline model)}} {The Area Under the Receiver Operator Characteristic Curve for the baseline model}
#' \item{\strong{AUC(new model)}} {The Area Under the Receiver Operator Characteristic Curve for the new (alt) model}
#' \item{\strong{AUC difference}} {The difference in the AUCs betwen the reference and new model with a confidence interval}
#' \item{\strong{difference (p)}} {P value for the difference in AUCs (DeLong method)}
#' \item{\strong{Brier(baseline model)}} {The Brier score for the baseline model}
#' \item{\strong{Brier(new model)}} {The Brier score for the alternate model}
#' \item{\strong{Brier skill}} {The percent improvement of the alternatve over the baseline model based on the relative change in Brier score}
#' \item{\strong{incidence}} {The incidence of the event}
#' }
#' }
#' @export 
#' @examples
#'\dontrun{
#'data(data_risk)
#'y<-data_risk$outcome 
#'x1<-data_risk$baseline
#'x2<-data_risk$new
#'t<-c(0,0.19,1) 
#'#e.g.
#'output<-CI.raplot(x1, x2, y, t, conf.level = 0.95, n.boot = 5, dp = 2) 
#'}
#' @references  Pencina, M. J., D'Agostino, R. B., & Vasan, R. S. (2008). Evaluating the added stats::predictive ability of a new marker: From area under the ROC curve to reclassification and beyond. Statistics in Medicine, 27(2), 157–172. doi:10.1002/sim.2929
CI.raplot <- function(x1, x2, y = NULL,  t = NULL,  conf.level = 0.95, n.boot = 2000, dp = 3) {
  
  if (class(x1)[1] == "glm") {
    y = x1$y # must come first
    x1 = stats::predict(x1, type = "response")
    data_type = "glm"
  }
  if (class(x2)[1] == "glm") {
    x2 = stats::predict(x2, type = "response")
  }
  if (class(x1)[1] == "lrm") {
    x1 = stats::predict(x1, type = "fitted")
    if (length(x1$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    y = x1$as.numeric(as.character(x1$y))
    data_type = "lrm"
  }
  if (class(x2)[1] == "lrm") {
    if (length(x2$y) == 0 ) {
      stop("Fitted models with the rms package must be made using y = TRUE")
    }
    x2 = stats::predict(x2, type = "fitted")
  }
  if (class(x1)[1] != "glm"  & class(x1)[1] != "lrm" ) {
    data_type = "User supplied"
  }
  if (length(x1) != length(x2))
    stop("Reference (baseline) and New (Alt) model vectors must be the same length")
  
  ifelse(is.null(t),
         results <- statistics.raplot(x1, x2, y, t = NULL),  
         results <- statistics.raplot(x1, x2, y,  t)
  )
  
  # Get results for each bootstrapped sample
  results.boot <- list()
  
  for (i in 1:n.boot) {
    boot.index <- sample(length(y), replace = TRUE)
    risk.model1.boot <- x1[boot.index]
    risk.model2.boot <- x2[boot.index]
    cc.status.boot <- y[boot.index]           
    results.boot[[i]] <- statistics.raplot(x1 = risk.model1.boot, 
                                           x2 = risk.model2.boot, 
                                           y = cc.status.boot,
                                           t)
  }
  
  results.matrix <- extractCI(results.boot = results.boot, conf.level = conf.level, n.boot = n.boot, dp = dp)
  
  thresh = ifelse(is.null(t), "baseline",t)
  
  meta_data <- data.frame(Thresholds = thresh, `Confidence interval` = 100 * conf.level, `Number of bootstraps` = n.boot, `Input data type` = data_type, `# decimal places` = dp)
  return(list(meta_data = meta_data, Metrics = results, Each_bootstrap_metrics = results.boot, Summary_metrics = results.matrix))
}


#' List risk assessment metrics
#' 
#' Display the summary metrics
#' 
#' @param l List returned from CI.raplot
#' @return A tibble
summary.rap = function(l) {
  summary_metrics <- l[["Summary_metrics"]]
  return(summary_metrics)
}

#' List meta data
#' 
#' Display the meta data
#' 
#' @param l List returned from CI.raplot
#' @return A tibble
meta.rap = function(l) {
  meta_data <- l[["meta_data"]]
  return(meta_data)
}

#' Reclassification metrics with classes (ordinals) as inputs
#' 
#' The function statistics.classNRI calculates the NRI metrics for reclassification of data already in classes. For use by CI.classNRI.
#' 
#' @param c1 Risk class of Reference model (ordinal factor).
#' @param c2 Risk class of New model (ordinal factor)
#' @param y Binary of outcome of interest. Must be 0 or 1. 
#' @param s1 The savings or benefit when an event is reclassified to a higher group by the new model. i.e instead of counting as 1 an event classified to a higher group, it is counted as s1.
#' @param s2 The benefit when a non-event is reclassified to a lower group.  i.e instead of counting as 1 an event classified to a lower group, it is counted as s2.  
#' @return A matrix of metrics for use within CI.classNRI
#' @examples
#'\dontrun{
#'data(data_class)
#'y <- data_class$outcome 
#'c1 <- data_class$base_class
#'c2 <- data_class$new_class
#'#e.g.
#'output<-statistics.classNRI(c1, c2, y) 
#'}
#' @export
statistics.classNRI <- function(c1, c2, y,s1 = NULL, s2 = NULL) {    
  
  if (!is.factor(c1)) {c1 <- as.factor(c1) }
  if (!is.factor(c2)) {c1 <- as.factor(c2) }
  
  df <- data.frame(c1 = c1, c2 = c2, event = y)
  
  #Remove rows with missing data
  df <- df %>%  
    mutate(event = as.numeric(event)) %>% 
    filter(!is.na(c1)) %>% 
    filter(!is.na(c2)) %>% 
    filter(!is.na(event))
  
  u <- sort(unique(df$event))
  if (length(u) != 2 || u[1] != 0 || u[2] != 1)
    stop("Outcome/Event must have two values: 0 and 1")

  
  n <- nrow(df)
  n_event = sum(df$event)
  incidence <- 100 * n_event/n
  
  #NRI
  df_NRI <- df %>% 
    group_by(event) %>% 
    mutate(up = ifelse(as.numeric(c2) > as.numeric(c1),1,0))  %>% 
    mutate(down = ifelse(as.numeric(c2) < as.numeric(c1),1,0)) %>% 
    summarise(
      n = n(),
      n_up = sum(up),
      n_down = sum(down),
      nri = (n_up - n_down)/n # for the non-event take the negative of this (below)
    ) %>% 
    ungroup() %>% 
    mutate(nri = ifelse(event == 0, -nri, nri))
  
  
  df_wNRI = NULL
  if (!is.null(s1) & !is.null(s2)) {
    df_wNRI <- df %>% 
      group_by(event) %>% 
      mutate(up = ifelse(as.numeric(c2) > as.numeric(c1),s1,0))  %>% 
      mutate(down = ifelse(as.numeric(c2) < as.numeric(c1),s2,0)) %>% 
      summarise(
        n = n(),
        n_up = sum(up),
        n_down = sum(down),
        nri = (n_up - n_down)/n # for the non-event take the negative of this (below)
      ) %>% 
      ungroup() %>% 
      mutate(nri = ifelse(event == 0, -nri, nri))
  }
  
  df_event <- df %>% filter(event == 1)
  confusion.matrix_event <- table(df_event$c1, df_event$c2, dnn = c("Baseline", "New"))
  
  df_nonevent <- df %>% filter(event == 0)
  confusion.matrix_nonevent <- table(df_nonevent$c1, df_nonevent$c2, dnn = c("Baseline", "New"))
  
  
  # Output
  output <- list(n = n, 
                 n_event = n_event, 
                 n_non_event = n - n_event, 
                 Prevalence = n_event/n,
                 NRI_up_event = as.integer(df_NRI[df_NRI$event == 1,]$n_up),
                 NRI_up_nonevent = as.integer(df_NRI[df_NRI$event == 0,]$n_up),
                 NRI_down_event = as.integer(df_NRI[df_NRI$event == 1,]$n_down),
                 NRI_down_nonevent = as.integer(df_NRI[df_NRI$event == 0,]$n_down),
                 NRI_event = df_NRI[df_NRI$event == 1,]$nri,
                 NRI_nonevent = df_NRI[df_NRI$event == 0,]$nri,
                 wNRI_event = df_wNRI[df_wNRI$event == 1,]$nri,
                 wNRI_nonevent = df_wNRI[df_wNRI$event == 0,]$nri,
                 confusion.matrix_event = confusion.matrix_event,
                 confusion.matrix_nonevent = confusion.matrix_nonevent
  )
  return(output)
}

#' Statistical metrics and confidence intervals for classes
#' 
#' The function CI.classNRI calculates the NRI statistics for reclassification of data already in classes with confidence intervals.  Uses statistics.classNRI.
#' 
#' @param c1 Risk classes of the baseline model (ordinal)
#' @param c2 Risk classes of new model
#' @param y Binary of outcome of interest. Must be 0 or 1. 
#' @param s1 The savings or benefit when am event is reclassified to a higher group by the new model (positive numeric)
#' @param s2 The benefit when a non-event is reclassified to a lower group (positive numeric)
#' @param conf.level The confidence interval expressed as a fraction of 1 (ie 0.95 is the 95\% confidence interval )
#' @param n.boot  The number of "bootstraps" to use. Performance slows down with more bootstraps. For trialling result, use a low number (eg 2), for accuracy use a large number (eg 2000)
#' @param dp The number of decimal places to display
#' @return A list with four items: 
#' \itemize{
#' \item{\strong{1. meta_data}} {Some overall meta data - Confidence Interval, number of bootstraps, s1, s2}
#' \item{\strong{2. Metrics}} {Point estimates of the statistical metrics (see list below)}
#' \item{\strong{3. Each_bootstrap_metrics}} {Point estimates of the statistical metrics for each bootstrapped sample (see list below)}
#' \item{\strong{4. Summary Metrics}} {Point estimates with confidence intervals of the statistical metrics. See following list:)}
#' \itemize{
#'  \item{\strong{Total (n)}} {Total number of subjects}
#' \item{\strong{Events (n)}} {Number of subjects with the event (outcome) of interest}
#' \item{\strong{Non-events (n)}} {Number of subjects without the event (outcome) of interest}
#' \item{\strong{Prevalence}} {The prevalence of the event}
#' \item{\strong{NRI events}} {The NRI with confidence interval for those with the event.}
#' \item{\strong{NRI non-events}} {The NRI with confidence interval for those without the event.}
#' \item{\strong{wNRI-events}} {The weighted NRI for those with the event}
#' \item{\strong{wNRI-nonevents}} {The weighted NRI for those without the event }
#' \item{\strong{Confusion matrix events}} {A confusion matric (table) showing the relationship of the numbers classified into each class with the baseline model compared to those in each class with the new model}
#' \item{\strong{Confusion matrix non-events}} {A confusion matric (table) showing the relationship of the numbers classified into each class with the baseline model compared to those in each class with the new model}
#' }
#' }
#' @return A matrix of metrics
#' @export
CI.classNRI <- function(c1, c2, y, s1 = NULL, s2 = NULL,  conf.level = 0.95, n.boot = 1000, dp = 3) {
  
  if (!is.factor(c1)) {c1 <- as.factor(c1) }
  if (!is.factor(c2)) {c1 <- as.factor(c2) }
  
  results <- statistics.classNRI(c1, c2,y,s1,s2)  
  
  # Get results for each bootstrapped sample
  results.boot <- list()
  
  for (i in 1:n.boot) {
    boot.index <- sample(length(y), replace = TRUE)
    risk.model1.boot <- c1[boot.index]
    risk.model2.boot <- c2[boot.index]
    cc.status.boot <- y[boot.index]           
    results.boot[[i]] <- statistics.classNRI(c1 = risk.model1.boot, 
                                             c2 = risk.model2.boot, 
                                             y = cc.status.boot,
                                             s1,
                                             s2)
  }
  
  results.matrix <- extract_NRI_CI(results.boot = results.boot, conf.level = conf.level, n.boot = n.boot, dp = dp) #Note, no CIs on the confusion matrices
  ifelse(is.null(s1),
         meta_data <- data.frame(`Confidence interval` = 100 * conf.level, `Number of bootstraps` = n.boot, `# decimal places` = dp),
         meta_data <- data.frame(`Confidence interval` = 100 * conf.level, `Number of bootstraps` = n.boot, `Weight s1` = s1, `Weight s2` = s2,`# decimal places` = dp))
  
  
  return(list(meta_data = meta_data, Metrics = results, Each_bootstrap_metrics = results.boot, Summary_metrics = results.matrix))
}