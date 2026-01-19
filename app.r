## app.R --------------------------------------------------------------
library(shiny)
library(bslib)
library(MASS)
library(plotly)
library(ggplot2)
library(shinyWidgets)
library(DT)
library(brglm2)
library(readxl)
library(rlang)
library(glmmTMB)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(knitr)

brglmControl(maxit = 5000)

## ---------- Shared helper functions for Ct app ----------------------

shapeWide <- function(df) {
  df %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    filter(if_any(starts_with("Ct"), ~ !is.na(.))) %>%
    pivot_longer(
      cols      = starts_with("Ct"),
      names_to  = "replicate",
      values_to = "Ct"
    ) %>%
    mutate(
      replicate = as.integer(sub("Ct", "", replicate))
    )
}

add_glmmTMB_preds <- function(mod,
                              level        = 0.95,
                              newdata      = NULL,
                              type         = "response",
                              prefix       = "pred_",
                              include_PI   = TRUE,
                              re.form      = NA,
                              allow.new.levels = FALSE,
                              decimal      = 2) {
  if (is.null(newdata)) {
    newdata <- stats::model.frame(mod)
  }
  
  z <- qnorm(1 - (1 - level) / 2)
  
  pop <- predict(
    mod,
    newdata = newdata,
    type    = type,
    re.form = re.form,
    se.fit  = TRUE
  )
  
  sigma  <- sigma(mod)
  sigma2 <- sigma^2
  
  newdata[[paste0(prefix, "mean")]] <- round(pop$fit, decimal)
  newdata[[paste0(prefix, "lwr")]]  <- round(pop$fit - z * pop$se.fit, decimal)
  newdata[[paste0(prefix, "upr")]]  <- round(pop$fit + z * pop$se.fit, decimal)
  
  if (include_PI) {
    newdata[[paste0(prefix, "pi_lwr")]] <-
      round(pop$fit - z * sqrt(pop$se.fit^2 + sigma2), decimal)
    
    newdata[[paste0(prefix, "pi_upr")]] <-
      round(pop$fit + z * sqrt(pop$se.fit^2 + sigma2), decimal)
  }
  
  newdata
}

plot_fits_by_exp_and_overall <- function(mod,
                                         data       = NULL,
                                         level      = 0.95,
                                         show_ci    = TRUE,
                                         show_points = TRUE,
                                         decimal    = 2) {
  if (is.null(data)) data <- stats::model.frame(mod)
  
  data <- data %>%
    mutate(exp = factor(exp)) %>%
    arrange(exp, dilution)
  
  d_cond <- add_glmmTMB_preds(
    mod,
    level   = level,
    newdata = data,
    prefix  = "cond_",
    include_PI = FALSE,
    re.form = NULL,
    decimal = decimal
  )
  
  d_pop <- add_glmmTMB_preds(
    mod,
    level   = level,
    newdata = data,
    prefix  = "pop_",
    include_PI = FALSE,
    re.form = NA,
    decimal = decimal
  )
  
  pop_df <- d_pop %>%
    distinct(dilution, .keep_all = TRUE) %>%
    arrange(dilution)
  
  p <- ggplot(d_cond, aes(x = dilution, y = Ct)) +
    { if (show_points) geom_point(aes(colour = exp), alpha = 0.6) } +
    geom_line(aes(y = cond_mean, colour = exp, group = exp), linewidth = 0.9) +
    geom_line(
      data        = pop_df,
      aes(x = dilution, y = pop_mean),
      inherit.aes = FALSE,
      linewidth   = 1.2,
      colour      = "white"
    ) +
    labs(x = "Dilution", y = "Ct", colour = "Experiment") +
    theme_minimal(base_size = 13) +
    theme(
      plot.background  = element_rect(fill = "#111827", colour = NA),
      panel.background = element_rect(fill = "#111827", colour = NA),
      legend.background = element_rect(fill = "#111827", colour = NA),
      legend.key       = element_rect(fill = "#111827", colour = NA),
      text             = element_text(colour = "#E5E7EB"),
      axis.text        = element_text(colour = "#D1D5DB"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "#374151"),
      panel.grid.major.y = element_line(color = "#374151")
    )
  
  if (show_ci) {
    p <- p +
      geom_ribbon(
        data        = pop_df,
        aes(x = dilution, ymin = pop_lwr, ymax = pop_upr),
        inherit.aes = FALSE,
        alpha       = 0.2,
        fill        = "#60A5FA"
      )
  }
  
  list(p = p, d_cond = d_cond, d_pop = d_pop)
}

get_sd_summary <- function(mod1, mod2) {
  if (!inherits(mod1, "glmmTMB")) stop("mod1 is not a glmmTMB model.")
  if (!inherits(mod2, "glmmTMB")) stop("mod2 is not a glmmTMB model.")
  
  sm1 <- summary(mod1)
  if (!"disp" %in% names(sm1$coefficients)) {
    stop(
      "mod1 has no dispersion component (no `disp` in summary(mod1)$coefficients). ",
      "Check that you fitted mod1 with `dispformula = ~ 0 + dilu`."
    )
  }
  
  disp_coef       <- sm1$coefficients$disp[, 1]
  var_by_dilution <- exp(disp_coef)
  
  dil_labels <- names(var_by_dilution)
  sd_by_dil  <- sqrt(as.numeric(var_by_dilution))
  
  df_dil <- data.frame(
    label = dil_labels,
    sd    = sd_by_dil,
    row.names  = NULL,
    check.names = FALSE
  )
  
  vc2 <- glmmTMB::VarCorr(mod2)
  if (!"cond" %in% names(vc2)) {
    stop("VarCorr(mod2) has no 'cond' component; unexpected VarCorr structure.")
  }
  cond_vc <- vc2$cond
  
  if (!"exp" %in% names(cond_vc)) {
    stop("VarCorr(mod2)$cond has no group 'exp'; check random-effects structure.")
  }
  
  re_exp           <- cond_vc$exp
  between_sd_mixed <- as.numeric(attr(re_exp, "stddev"))[1L]
  within_sd_mixed  <- sigma(mod2)
  total_sd_mixed   <- sqrt(between_sd_mixed^2 + within_sd_mixed^2)
  
  df_comp <- data.frame(
    label = c("mixed_intra", "mixed_inter", "mixed_total"),
    sd    = c(within_sd_mixed, between_sd_mixed, total_sd_mixed),
    row.names  = NULL,
    check.names = FALSE
  )
  
  mf <- model.frame(mod2)
  
  lm_mod  <- lm(Ct ~ dilution + exp, data = mf)
  aov_tab <- anova(lm_mod)
  
  MS_between <- aov_tab["exp",       "Mean Sq"]
  MS_within  <- aov_tab["Residuals", "Mean Sq"]
  
  n_rep_by_exp <- table(mf$exp)
  k_vec <- as.numeric(n_rep_by_exp)
  
  if (length(unique(k_vec)) != 1L) {
    warning("Design is not perfectly balanced; ANOVA variance components are approximate.")
  }
  k <- mean(k_vec)
  
  sigma2_within_anova  <- MS_within
  sigma2_between_anova <- (MS_between - MS_within) / k
  sigma2_between_anova <- max(sigma2_between_anova, 0)
  
  sd_within_anova  <- sqrt(sigma2_within_anova)
  sd_between_anova <- sqrt(sigma2_between_anova)
  sd_total_anova   <- sqrt(sigma2_within_anova + sigma2_between_anova)
  
  df_anova <- data.frame(
    label = c("ANOVA_between", "ANOVA_within", "ANOVA_total"),
    sd    = c(sd_between_anova, sd_within_anova, sd_total_anova),
    row.names  = NULL,
    check.names = FALSE
  )
  
  rbind(df_dil, df_comp, df_anova)
}

## ---------- Diagnostic 2×2 helper + module -----------------------

diag_prop_ci <- function(x, n, conf.level = 0.95) {
  if (is.na(x) || is.na(n) || n <= 0) {
    return(c(Estimate = NA_real_, Lower = NA_real_, Upper = NA_real_))
  }
  bt <- binom.test(round(x), round(n), conf.level = conf.level)
  c(
    Estimate = unname(bt$estimate),
    Lower    = bt$conf.int[1],
    Upper    = bt$conf.int[2]
  )
}

diag2x2UI <- function(id) {
  ns <- NS(id)
  
  tabsetPanel(
    type = "pills",
    id   = ns("diag_tabs"),
    
    tabPanel(
      "Calculator",
      h3("Diagnostic Test 2×2 Calculator"),
      fluidRow(
        column(
          width = 5,
          tags$h4("Enter 2×2 table"),
          tags$table(
            class = "table table-bordered table-sm diag-input-table",
            tags$thead(
              tags$tr(
                tags$th(""),
                tags$th(colspan = 2, "Condition"),
                tags$th("Totals")
              ),
              tags$tr(
                tags$th(""),
                tags$th("Absent"),
                tags$th("Present"),
                tags$th("")
              )
            ),
            tags$tbody(
              tags$tr(
                tags$th("Test Positive"),
                tags$td(
                  numericInput(
                    ns("fp"), NULL,
                    value = 0, min = 0, step = 1
                  )
                ),
                tags$td(
                  numericInput(
                    ns("tp"), NULL,
                    value = 0, min = 0, step = 1
                  )
                ),
                tags$td(textOutput(ns("row_pos_tot")))
              ),
              tags$tr(
                tags$th("Test Negative"),
                tags$td(
                  numericInput(
                    ns("tn"), NULL,
                    value = 0, min = 0, step = 1
                  )
                ),
                tags$td(
                  numericInput(
                    ns("fn"), NULL,
                    value = 0, min = 0, step = 1
                  )
                ),
                tags$td(textOutput(ns("row_neg_tot")))
              ),
              tags$tr(
                tags$th("Totals"),
                tags$td(textOutput(ns("col_abs_tot"))),
                tags$td(textOutput(ns("col_pres_tot"))),
                tags$td(textOutput(ns("grand_tot")))
              )
            )
          ),
          div(
            style = "margin-top: 10px; display:flex; gap:8px; flex-wrap:wrap;",
            actionButton(ns("calc"), "Calculate", class = "btn btn-primary"),
            actionButton(ns("reset"), "Reset", class = "btn btn-secondary"),
            downloadButton(ns("download_docx"), "Download DOCX", class = "btn btn-success")
          ),
          div(
            style = "margin-top: 8px; max-width: 220px;",
            numericInput(
              ns("digits"),
              label = "Digits to display:",
              value = 6,
              min   = 1,
              max   = 10
            )
          ),
          tags$small(
            class = "text-muted",
            "Cells: TP = Test positive & Condition present; ",
            "TN = Test negative & Condition absent; ",
            "FP = Test positive & Condition absent; ",
            "FN = Test negative & Condition present."
          )
        ),
        
        column(
          width = 7,
          h4("Summary"),
          br(),
          h5("Prevalence, sensitivity, specificity"),
          tableOutput(ns("tbl_main")),
          br(),
          h5("Overall probability that the test result will be:"),
          tableOutput(ns("tbl_test_result")),
          br(),
          h5("For a positive test result, probability that it is:"),
          tableOutput(ns("tbl_pos_result")),
          br(),
          h5("For a negative test result, probability that it is:"),
          tableOutput(ns("tbl_neg_result")),
          br(),
          h5("Likelihood Ratios"),
          tableOutput(ns("tbl_lr"))
        )
      )
    ),
    
    tabPanel(
      "About",
      br(),
      fluidRow(
        column(
          10,
          h3("About the Diagnostic 2×2 Calculator"),
          p(
            "This tab reproduces the core logic of a classical 2×2 diagnostic test calculator. ",
            "You enter counts of true positives (TP), false positives (FP), true negatives (TN), ",
            "and false negatives (FN). The app then derives standard diagnostic performance measures ",
            "with exact binomial confidence intervals."
          ),
          h4("Quantities computed"),
          tags$ul(
            tags$li(
              strong("Prevalence: "),
              "proportion of subjects truly having the condition ",
              "( (TP + FN) / N )."
            ),
            tags$li(
              strong("Sensitivity: "),
              "probability the test is positive among those with the condition ",
              "( TP / (TP + FN) )."
            ),
            tags$li(
              strong("Specificity: "),
              "probability the test is negative among those without the condition ",
              "( TN / (TN + FP) )."
            ),
            tags$li(
              strong("Positive predictive value (PPV): "),
              "probability the condition is present given a positive result ",
              "( TP / (TP + FP) )."
            ),
            tags$li(
              strong("Negative predictive value (NPV): "),
              "probability the condition is absent given a negative result ",
              "( TN / (TN + FN) )."
            ),
            tags$li(
              strong("Likelihood ratios: "),
              "we report both conventional likelihood ratios (based only on sensitivity and specificity) ",
              "and prevalence-weighted odds ratios, which incorporate the underlying prevalence."
            )
          ),
          p(
            "Proportions and confidence intervals are based on the exact binomial test ",
            "(via ",
            code("binom.test"),
            " in R), mirroring the style of established clinical calculators."
          ),
          h4("Reference / original calculator"),
          p(
            "This module is conceptually aligned with the well-known VassarStats 2×2 diagnostic test calculator:"
          ),
          tags$p(
            tags$a(
              href   = "http://www.vassarstats.net/clin1.html",
              target = "_blank",
              "VassarStats: Diagnostic Test Calculator (2×2 table)"
            )
          ),
          p(
            "The Shiny implementation here is not an exact clone, but it computes the same core quantities ",
            "and presents them in a similar structured layout for routine assay validation and method comparison."
          )
        )
      )
    )
  )
}

diag2x2Server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    totals <- reactive({
      tp <- as.numeric(input$tp); fp <- as.numeric(input$fp)
      fn <- as.numeric(input$fn); tn <- as.numeric(input$tn)
      list(
        tp = tp, fp = fp, fn = fn, tn = tn,
        row_pos  = tp + fp,
        row_neg  = fn + tn,
        col_abs  = fp + tn,
        col_pres = tp + fn,
        grand    = tp + fp + fn + tn
      )
    })
    
    output$row_pos_tot  <- renderText(totals()$row_pos)
    output$row_neg_tot  <- renderText(totals()$row_neg)
    output$col_abs_tot  <- renderText(totals()$col_abs)
    output$col_pres_tot <- renderText(totals()$col_pres)
    output$grand_tot    <- renderText(totals()$grand)
    
    observeEvent(input$reset, {
      updateNumericInput(session, "fp", value = 0)
      updateNumericInput(session, "tp", value = 23)
      updateNumericInput(session, "tn", value = 92)
      updateNumericInput(session, "fn", value = 1)
    })
    
    stats_list <- eventReactive(input$calc, ignoreNULL = FALSE, {
      tt  <- totals()
      tp  <- tt$tp; fp <- tt$fp; fn <- tt$fn; tn <- tt$tn
      N   <- tt$grand
      
      if (N <= 0) return(NULL)
      
      prev   <- diag_prop_ci(tp + fn, N)
      sens   <- diag_prop_ci(tp, tp + fn)
      spec   <- diag_prop_ci(tn, tn + fp)
      
      p_pos  <- diag_prop_ci(tp + fp, N)
      p_neg  <- diag_prop_ci(tn + fn, N)
      
      ppv    <- if (tp + fp > 0) diag_prop_ci(tp, tp + fp) else rep(NA, 3)
      fpv    <- if (tp + fp > 0) diag_prop_ci(fp, tp + fp) else rep(NA, 3)
      npv    <- if (tn + fn > 0) diag_prop_ci(tn, tn + fn) else rep(NA, 3)
      fnv    <- if (tn + fn > 0) diag_prop_ci(fn, tn + fn) else rep(NA, 3)
      
      Se       <- sens["Estimate"]
      Sp       <- spec["Estimate"]
      prev_est <- prev["Estimate"]
      
      lr_pos_c_est <- if (!is.na(Se) && !is.na(Sp) && Sp < 1) Se / (1 - Sp) else Inf
      lr_neg_c_est <- if (!is.na(Se) && !is.na(Sp) && Sp > 0) (1 - Se) / Sp else NA_real_
      
      SeL <- sens["Lower"]; SeU <- sens["Upper"]
      SpL <- spec["Lower"]; SpU <- spec["Upper"]
      
      lr_pos_c_ci <- c(NA, NA)
      lr_neg_c_ci <- c(NA, NA)
      
      if (!is.na(SeL) && !is.na(SeU) && !is.na(SpL) && !is.na(SpU) && SpU < 1) {
        lr_pos_c_lo <- SeL / (1 - SpU)
        lr_pos_c_hi <- SeU / (1 - SpL)
        lr_pos_c_ci <- c(lr_pos_c_lo, lr_pos_c_hi)
      }
      
      if (!is.na(SeL) && !is.na(SeU) && !is.na(SpL) && !is.na(SpU) && SpL > 0) {
        lr_neg_c_lo <- (1 - SeU) / SpU
        lr_neg_c_hi <- (1 - SeL) / SpL
        lr_neg_c_ci <- c(lr_neg_c_lo, lr_neg_c_hi)
      }
      
      odds_prev <- if (!is.na(prev_est) && prev_est > 0 && prev_est < 1) {
        prev_est / (1 - prev_est)
      } else NA_real_
      
      lr_pos_w_est <- if (is.finite(lr_pos_c_est) && !is.na(odds_prev)) {
        lr_pos_c_est * odds_prev
      } else NA_real_
      lr_neg_w_est <- if (!is.na(lr_neg_c_est) && !is.na(odds_prev)) {
        lr_neg_c_est * odds_prev
      } else NA_real_
      
      lr_pos_w_ci <- if (!any(is.na(c(lr_pos_c_ci, odds_prev)))) lr_pos_c_ci * odds_prev else c(NA, NA)
      lr_neg_w_ci <- if (!any(is.na(c(lr_neg_c_ci, odds_prev)))) lr_neg_c_ci * odds_prev else c(NA, NA)
      
      list(
        N     = N,
        prev  = prev,
        sens  = sens,
        spec  = spec,
        p_pos = p_pos,
        p_neg = p_neg,
        ppv   = ppv,
        fpv   = fpv,
        npv   = npv,
        fnv   = fnv,
        lr_pos_c = c(Estimate = lr_pos_c_est, Lower = lr_pos_c_ci[1], Upper = lr_pos_c_ci[2]),
        lr_neg_c = c(Estimate = lr_neg_c_est, Lower = lr_neg_c_ci[1], Upper = lr_neg_c_ci[2]),
        lr_pos_w = c(Estimate = lr_pos_w_est, Lower = lr_pos_w_ci[1], Upper = lr_pos_w_ci[2]),
        lr_neg_w = c(Estimate = lr_neg_w_est, Lower = lr_neg_w_ci[1], Upper = lr_neg_w_ci[2])
      )
    })
    
    # format as character to avoid renderTable rounding to 2 digits
    fmt_tbl <- function(x, digits = 6) {
      if (is.null(x)) return(NULL)
      df <- as.data.frame(x)
      rn <- rownames(df)
      df[] <- lapply(df, function(col) {
        if (is.numeric(col)) {
          formatC(col, digits = digits, format = "g")
        } else {
          as.character(col)
        }
      })
      rownames(df) <- rn
      df
    }
    
    output$tbl_main <- renderTable({
      s <- stats_list(); if (is.null(s)) return(NULL)
      digs <- if (!is.null(input$digits)) input$digits else 6
      df <- rbind(
        Prevalence  = s$prev,
        Sensitivity = s$sens,
        Specificity = s$spec
      )
      fmt_tbl(df, digits = digs)
    }, rownames = TRUE)
    
    output$tbl_test_result <- renderTable({
      s <- stats_list(); if (is.null(s)) return(NULL)
      digs <- if (!is.null(input$digits)) input$digits else 6
      df <- rbind(
        `Positive` = s$p_pos,
        `Negative` = s$p_neg
      )
      fmt_tbl(df, digits = digs)
    }, rownames = TRUE)
    
    output$tbl_pos_result <- renderTable({
      s <- stats_list(); if (is.null(s)) return(NULL)
      digs <- if (!is.null(input$digits)) input$digits else 6
      df <- rbind(
        `True Positive (PPV)` = s$ppv,
        `False Positive`      = s$fpv
      )
      fmt_tbl(df, digits = digs)
    }, rownames = TRUE)
    
    output$tbl_neg_result <- renderTable({
      s <- stats_list(); if (is.null(s)) return(NULL)
      digs <- if (!is.null(input$digits)) input$digits else 6
      df <- rbind(
        `True Negative (NPV)` = s$npv,
        `False Negative`      = s$fnv
      )
      fmt_tbl(df, digits = digs)
    }, rownames = TRUE)
    
    output$tbl_lr <- renderTable({
      s <- stats_list(); if (is.null(s)) return(NULL)
      digs <- if (!is.null(input$digits)) input$digits else 6
      df <- rbind(
        `Positive [C]` = s$lr_pos_c,
        `Negative [C]` = s$lr_neg_c,
        `Positive [W]` = s$lr_pos_w,
        `Negative [W]` = s$lr_neg_w
      )
      fmt_tbl(df, digits = digs)
    }, rownames = TRUE)
    
    
    output$download_docx <- downloadHandler(
      filename = function() {
        paste0("diagnostic2x2_summary_", Sys.Date(), ".docx")
      },
      content = function(file) {
        s <- stats_list()
        if (is.null(s)) {
          stop("No results to export yet — click Calculate first.")
        }
        
        if (!requireNamespace("rmarkdown", quietly = TRUE) || !requireNamespace("knitr", quietly = TRUE)) {
          stop("To export DOCX on Posit Cloud, please install packages: rmarkdown and knitr.")
        }
        
        digs <- if (!is.null(input$digits)) input$digits else 6
        
        prep_df <- function(df) {
          df2 <- fmt_tbl(df, digits = digs)
          df2 <- cbind(Measure = rownames(df2), df2)
          rownames(df2) <- NULL
          df2
        }
        
        df_main <- prep_df(rbind(
          Prevalence  = s$prev,
          Sensitivity = s$sens,
          Specificity = s$spec
        ))
        
        df_test_result <- prep_df(rbind(
          Positive = s$p_pos,
          Negative = s$p_neg
        ))
        
        df_pos_result <- prep_df(rbind(
          `True Positive (PPV)` = s$ppv,
          `False Positive`      = s$fpv
        ))
        
        df_neg_result <- prep_df(rbind(
          `True Negative (NPV)` = s$npv,
          `False Negative`      = s$fnv
        ))
        
        df_lr <- prep_df(rbind(
          `Positive [C]` = s$lr_pos_c,
          `Negative [C]` = s$lr_neg_c,
          `Positive [W]` = s$lr_pos_w,
          `Negative [W]` = s$lr_neg_w
        ))
        
        tt <- totals()
        
        # Render a tiny Rmd to Word using Pandoc (more Posit-Cloud friendly than officer)
        rmd <- tempfile(fileext = ".Rmd")
        out_dir <- tempdir()
        out_name <- paste0("diagnostic2x2_summary_", Sys.Date(), ".docx")
        
        writeLines(c(
          "---",
          "title: 'Diagnostic 2×2 summary'",
          "output: word_document",
          "---",
          "",
          sprintf("TP=%s, FP=%s, TN=%s, FN=%s (N=%s)", tt$tp, tt$fp, tt$tn, tt$fn, tt$grand),
          "",
          sprintf("Digits shown: %s", digs),
          "",
          "## Prevalence, sensitivity, specificity",
          "```{r, echo=FALSE}",
          "knitr::kable(df_main)",
          "```",
          "",
          "## Overall probability that the test result will be",
          "```{r, echo=FALSE}",
          "knitr::kable(df_test_result)",
          "```",
          "",
          "## For a positive test result, probability that it is",
          "```{r, echo=FALSE}",
          "knitr::kable(df_pos_result)",
          "```",
          "",
          "## For a negative test result, probability that it is",
          "```{r, echo=FALSE}",
          "knitr::kable(df_neg_result)",
          "```",
          "",
          "## Likelihood ratios",
          "```{r, echo=FALSE}",
          "knitr::kable(df_lr)",
          "```"
        ), con = rmd)
        
        rmarkdown::render(
          input        = rmd,
          output_file  = out_name,
          output_dir   = out_dir,
          quiet        = TRUE,
          envir        = list2env(list(
            df_main        = df_main,
            df_test_result = df_test_result,
            df_pos_result  = df_pos_result,
            df_neg_result  = df_neg_result,
            df_lr          = df_lr
          ), parent = globalenv())
        )
        
        file.copy(file.path(out_dir, out_name), file, overwrite = TRUE)
      }
    )
  })
}

## ---------- Ct regression module -----------------------

ctAppUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        3,
        div(
          class = "ct-panel",
          wellPanel(
            h4(
              span(icon("vial", class = "hero-icon")),
              "Ct regression / variance"
            ),
            fileInput(ns("file_ct"),
                      label  = h5("Upload Ct data (xlsx):"),
                      accept = c(".xlsx")),
            tags$small(
              "Wide format with dilution, experiment, and Ct1/Ct2/... columns."
            ),
            hr(),
            uiOutput(ns("namesY")),
            uiOutput(ns("namesD")),
            uiOutput(ns("expLevels")),
            numericInput(
              ns("decimal"),
              h5("Decimal places to display:"),
              value = 2,
              min   = 0,
              max   = 10
            ),
            hr(),
            tags$p(
              strong("Sample file: "),
              tags$a(href = "testWide.xlsx", "testWide.xlsx", target = "_blank")
            ),
            tags$small(
              style = "color:#9CA3AF;",
              "Tip: each experiment should have Ct1, Ct2, ... columns at the same dilutions."
            )
          )
        )
      ),
      column(
        9,
        tabsetPanel(
          type = "pills",
          id   = ns("ct_tabs"),
          
          tabPanel(
            title = tagList(icon("table"), "Data & predictions"),
            br(),
            DTOutput(ns("table"))
          ),
          
          tabPanel(
            title = tagList(icon("chart-bar"), "SD / variance summary"),
            br(),
            DTOutput(ns("table1"))
          ),
          
          tabPanel(
            title = tagList(icon("line-chart"), "Ct vs dilution plot"),
            br(),
            plotlyOutput(ns("plotly"), height = "600px")
          ),
          
          tabPanel(
            title = tagList(icon("circle-info"), "About (Ct regression)"),
            br(),
            fluidRow(
              column(
                12,
                h3("About this app: Ct regression and variance decomposition"),
                p(
                  "This tab analyses real-time PCR data using Ct values with a regression/mixed-model framework. ",
                  "Ct is treated as a continuous outcome and modelled as a function of dilution, ",
                  "with experiments (runs/plates) included as a random effect."
                ),
                h4("Statistical model"),
                p(
                  "The main model is a Gaussian mixed model of the form ",
                  code("Ct ~ dilution + (1 | exp)"),
                  ", fitted with ",
                  code("glmmTMB"),
                  "."
                ),
                tags$ul(
                  tags$li(
                    strong("Fixed effect of dilution: "),
                    "captures how the mean Ct changes with dilution across all experiments."
                  ),
                  tags$li(
                    strong("Random intercept for exp: "),
                    "allows each experiment (run/plate) to have its own baseline Ct, ",
                    "capturing between-experiment variability."
                  )
                ),
                p(
                  "In one version of the model, a dispersion formula ",
                  code("dispformula = ~ 0 + dilu"),
                  " is used so that the residual variance can differ by dilution. ",
                  "In the simpler model without a dispersion formula, a single residual variance is assumed across all dilutions."
                ),
                
                h4("Interpretation of the SDs shown in the app"),
                p("The app summarises several standard deviations (SD) that correspond to different sources of variation."),
                
                h5("By-dilution SDs from the dispersion model"),
                tags$ul(
                  tags$li(
                    strong("dilu1, dilu2, …: "),
                    "these SDs come from the model with ",
                    code("dispformula = ~ 0 + dilu"),
                    ". ",
                    "Each value represents the residual SD of Ct for that specific dilution, ",
                    "after adjusting for the fixed effect of dilution and the experiment random effect. ",
                    "They describe how variable Ct is at each dilution level."
                  )
                ),
                
                h5("Mixed-model SDs (from glmmTMB without dispersion model)"),
                tags$ul(
                  tags$li(
                    strong("mixed_intra: "),
                    "the within-experiment (residual) SD from the homogeneous-variance mixed model. ",
                    "It reflects repeatability: how much individual Ct values vary around the fitted mean within the same experiment."
                  ),
                  tags$li(
                    strong("mixed_inter: "),
                    "the between-experiment SD (random intercept SD for ",
                    code("exp"),
                    "). ",
                    "This describes how much experiment-level baselines differ from run to run."
                  ),
                  tags$li(
                    strong("mixed_total: "),
                    "the combined SD for a new Ct measurement when you do not know which experiment it will come from. ",
                    "It is computed as ",
                    HTML("&#8730;(mixed_intra<sup>2</sup> + mixed_inter<sup>2</sup>).")
                  )
                ),
                
                h5("ANOVA-based SDs (fixed-effects view)"),
                p(
                  "For comparison, an ordinary linear model ",
                  code("lm(Ct ~ dilution + exp)"),
                  " is also fitted and an ANOVA table is used to derive variance components, ",
                  "treating ",
                  code("exp"),
                  " as if it were a random factor in a classical ANOVA framework."
                ),
                tags$ul(
                  tags$li(
                    strong("ANOVA_within: "),
                    "square root of the residual mean square from the ANOVA. ",
                    "This is the within-experiment SD in the fixed-effects view and is comparable to ",
                    code("mixed_intra"),
                    "."
                  ),
                  tags$li(
                    strong("ANOVA_between: "),
                    "an ANOVA-based estimate of the between-experiment SD, ",
                    "derived from the difference between the experiment mean square and the residual mean square ",
                    "scaled by the number of replicates per experiment (using the expected mean squares approach). ",
                    "Conceptually comparable to ",
                    code("mixed_inter"),
                    ", though estimated in a different way."
                  ),
                  tags$li(
                    strong("ANOVA_total: "),
                    "the total SD combining ANOVA_between and ANOVA_within, ",
                    "computed as ",
                    HTML("&#8730;(ANOVA_between<sup>2</sup> + ANOVA_within<sup>2</sup>). "),
                    "This is the ANOVA analogue of ",
                    code("mixed_total"),
                    "."
                  )
                ),
                
                h4("How to interpret these in practice"),
                tags$ul(
                  tags$li(
                    strong("By-dilution SDs (dilu1–diluK): "),
                    "tell you how variable Ct is at each dilution level."
                  ),
                  tags$li(
                    strong("mixed_intra / ANOVA_within: "),
                    "repeatability: variation of Ct within the same experiment/run."
                  ),
                  tags$li(
                    strong("mixed_inter / ANOVA_between: "),
                    "reproducibility across experiments: how much runs/plates differ in their baseline Ct."
                  ),
                  tags$li(
                    strong("mixed_total / ANOVA_total: "),
                    "overall measurement uncertainty for a single new Ct value in routine use, ",
                    "combining both within- and between-experiment variation."
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

ctAppServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    getNames <- reactive({
      inFile1 <- input$file_ct
      if (is.null(inFile1)) return(NULL)
      theD    <- readxl::read_xlsx(inFile1$datapath)
      theNames <- names(theD)
      list(theNames = theNames, N = nrow(theD))
    })
    
    Data <- reactive({
      inFile1 <- input$file_ct
      if (is.null(inFile1)) return(NULL)
      req(input$D, input$Y, input$decimal)
      
      theD <- readxl::read_xlsx(inFile1$datapath)
      theD <- shapeWide(theD)
      
      myD <- theD %>%
        dplyr::select(Ct, dilution = all_of(input$D), exp = all_of(input$Y)) %>%
        dplyr::mutate(
          exp  = as.factor(exp),
          dilu = as.factor(dilution)
        )
      
      if (!is.null(input$exp_keep) && length(input$exp_keep) > 0) {
        myD <- myD %>% filter(exp %in% input$exp_keep)
      }
      
      mod1 <- glmmTMB::glmmTMB(
        Ct ~ dilution + (1 | exp),
        dispformula = ~ 0 + dilu,
        data        = myD,
        family      = gaussian()
      )
      
      mod2 <- glmmTMB::glmmTMB(
        Ct ~ dilution + (1 | exp),
        data        = myD,
        family      = gaussian()
      )
      
      var_df <- get_sd_summary(mod1, mod2) %>%
        mutate(sd = round(sd, input$decimal))
      
      r <- plot_fits_by_exp_and_overall(
        mod2,
        myD,
        level    = 0.95,
        show_ci  = TRUE,
        decimal  = input$decimal
      )
      
      list(dd = r$d_pop, var = var_df, plot = r$p)
    })
    
    output$namesY <- renderUI({
      theNames <- getNames()$theNames
      if (is.null(theNames)) return(NULL)
      theL <- as.list(theNames); names(theL) <- theNames
      selectInput(
        session$ns("Y"),
        label   = h5("Select column for experiments"),
        choices = theL,
        selected = theL[1]
      )
    })
    
    output$namesD <- renderUI({
      theNames <- getNames()$theNames
      if (is.null(theNames)) return(NULL)
      theL <- as.list(theNames); names(theL) <- theNames
      selectInput(
        session$ns("D"),
        label   = h5("Select column for dilution (log10)"),
        choices = theL,
        selected = theL[2]
      )
    })
    
    output$expLevels <- renderUI({
      req(input$file_ct, input$Y)
      theD   <- readxl::read_xlsx(input$file_ct$datapath)
      exp_vec <- theD[[input$Y]]
      choices <- sort(unique(as.character(exp_vec)))
      choices <- choices[!is.na(choices) & choices != ""]
      
      shinyWidgets::pickerInput(
        inputId  = session$ns("exp_keep"),
        label    = "Select experiment(s) to include",
        choices  = choices,
        selected = choices,
        multiple = TRUE,
        options  = list(`actions-box` = TRUE, `live-search` = TRUE)
      )
    })
    
    output$table <- DT::renderDataTable({
      if (is.null(input$file_ct)) {
        return(
          datatable(
            data.frame(Message = "Upload a Ct data file (xlsx) to get started."),
            options = list(dom = "t"),
            rownames = FALSE
          )
        )
      }
      temp1 <- Data()$dd
      if (is.null(temp1)) temp1 <- data.frame()
      if (nrow(temp1) == 0) {
        datatable(data.frame(No_Data = NULL),
                  options = list(dom = "t"),
                  rownames = FALSE)
      } else {
        datatable(
          temp1,
          rownames   = FALSE,
          extensions = "Buttons",
          options    = list(
            dom        = "Bfrtip",
            buttons    = c("copy", "csv", "excel"),
            pageLength = 25,
            lengthMenu = c(10, 25, 50, 100)
          )
        )
      }
    })
    
    output$table1 <- DT::renderDataTable({
      if (is.null(input$file_ct)) {
        return(
          datatable(
            data.frame(Message = "Upload a Ct data file (xlsx) to see SD / variance summary."),
            options = list(dom = "t"),
            rownames = FALSE
          )
        )
      }
      temp1 <- Data()$var
      if (is.null(temp1)) temp1 <- data.frame()
      if (nrow(temp1) == 0) {
        datatable(data.frame(No_Data = NULL),
                  options = list(dom = "t"),
                  rownames = FALSE)
      } else {
        datatable(
          temp1,
          rownames   = FALSE,
          extensions = "Buttons",
          options    = list(
            dom        = "Bfrtip",
            buttons    = c("copy", "csv", "excel"),
            pageLength = 15,
            lengthMenu = c(10, 15, 30, 50)
          )
        )
      }
    })
    
    output$plotly <- renderPlotly({
      if (is.null(input$file_ct)) return(NULL)
      req(Data()$plot)
      ggplotly(Data()$plot)
    })
  })
}

## ---------- Probit LOD module ------------------------

lodAppUI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(
        3,
        div(
          class = "lod-panel",
          wellPanel(
            h4(
              span(icon("bullseye", class = "hero-icon")),
              "Probit LOD analysis"
            ),
            fileInput(
              ns("file1"),
              label  = h5("Upload data in csv format:"),
              accept = c(".csv")
            ),
            uiOutput(ns("namesN")),
            uiOutput(ns("namesY")),
            checkboxInput(ns("manD"), "Use dilution column in the data", FALSE),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == false", ns("manD")),
              numericInput(
                ns("dilution"),
                h5("Initial dilution in \\( \\log_{10} \\):"),
                -1,
                min = -10,
                max = 0
              )
            ),
            
            numericInput(
              ns("copy"),
              h5("Copy number for initial dilution:"),
              500,
              min = 0,
              max = 1e6
            ),
            numericInput(
              ns("Ncopy"),
              h5("Number of decimal places for copy number results:"),
              4,
              min = 0,
              max = 10
            ),
            
            conditionalPanel(
              condition = sprintf("input['%s'] == true", ns("manD")),
              uiOutput(ns("namesD"))
            ),
            
            selectInput(
              ns("p"),
              h5("Percentile of interest:"),
              c(
                seq(0.01, 0.10, by = 0.01),
                seq(0.15, 0.85, by = 0.05),
                seq(0.90, 0.99, by = 0.01)
              ),
              selected = 0.95
            ),
            hr(),
            tags$p(
              strong("Sample file: "),
              tags$a(href = "Book1.csv", "Book1.csv", target = "_blank")
            ),
            tags$small(
              style = "color:#9CA3AF;",
              "Assumes grouped data: one row per dilution with total N and positives Y."
            )
          )
        )
      ),
      column(
        9,
        tabsetPanel(
          type = "pills",
          id   = ns("lod_tabs"),
          
          tabPanel(
            title = tagList(icon("table"), "Input data"),
            br(),
            DTOutput(ns("table"))
          ),
          
          tabPanel(
            title = tagList(icon("dna"), "LOD: copy number"),
            br(),
            DTOutput(ns("tableRRC"))
          ),
          
          tabPanel(
            title = tagList(icon("flask"), "LOD: dilution"),
            br(),
            DTOutput(ns("tableRR"))
          ),
          
          tabPanel(
            title = tagList(icon("chart-line"), "Percentile plot"),
            br(),
            plotOutput(ns("plot"), height = "600px")
          ),
          
          tabPanel(
            title = tagList(icon("circle-info"), "About (Probit LOD)"),
            br(),
            fluidRow(
              column(
                12,
                h3("About this app: Probit analysis for LOD estimation"),
                p(
                  "This tab performs a probit regression analysis to estimate the limit of detection (LOD) ",
                  "from real-time PCR or similar binary outcome experiments. ",
                  "The data are assumed to be grouped by dilution level, with each group providing ",
                  "the total number of replicates and the number of positive results."
                ),
                
                h4("Input data and columns"),
                tags$ul(
                  tags$li(
                    strong("Number of replicates (N): "),
                    "total number of tests performed at each dilution."
                  ),
                  tags$li(
                    strong("Number of positives (Y): "),
                    "number of positive results at each dilution."
                  ),
                  tags$li(
                    strong("Dilution (log10): "),
                    "log10 dilution value. ",
                    "You can either supply this column in the data (",
                    code("Use dilution column in the data"),
                    " checked), or set an initial dilution and let the app generate a sequence."
                  )
                ),
                p(
                  "The app derives a copy number for each dilution assuming: ",
                  "copy number at the initial dilution × ",
                  HTML("10<sup>(log dilution − initial dilution)</sup>"),
                  ", scaled by the user-specified initial copy number."
                ),
                
                h4("Statistical method"),
                p(
                  "For each dilution level, the proportion of positive results is modelled using a probit regression:"
                ),
                tags$ul(
                  tags$li(
                    strong("Model: "),
                    "a generalized linear model with binomial response ",
                    code("cbind(Y, N − Y)"),
                    " and probit link is fitted using ",
                    code("glm(..., family = binomial(link = 'probit'), method = 'brglmFit')"),
                    " from the ",
                    code("brglm2"),
                    " package to reduce bias, especially at extreme proportions (0% or 100%)."
                  ),
                  tags$li(
                    strong("Dose–response curve: "),
                    "the fitted probit model defines the relationship between log10 dilution and ",
                    "the probability of detection."
                  ),
                  tags$li(
                    strong("Percentiles of interest: "),
                    "the app uses ",
                    code("dose.p"),
                    " to calculate the dilution (and corresponding copy number) ",
                    "at user-selected detection probabilities (e.g. 1%, 50%, 95%, 99%)."
                  )
                ),
                
                h4("Description of outputs"),
                h5("Input data"),
                p(
                  "Shows the processed dataset used in the analysis. ",
                  "This includes the dilution values used in the model and derived quantities such as copy number."
                ),
                
                h5("LOD: dilution"),
                p("This table reports, for each probability level (percentile) evaluated:"),
                tags$ul(
                  tags$li(
                    strong("Estimate: "),
                    "the estimated log10 dilution at which that fraction of replicates would be positive."
                  ),
                  tags$li(
                    strong("LowerCI, UpperCI: "),
                    "the 95% confidence interval for the dilution associated with that percentile."
                  )
                ),
                p(
                  "The percentile selected in the sidebar (e.g. 95%) is highlighted in blue."
                ),
                
                h5("LOD: copy number"),
                p(
                  "This table converts the estimated dilutions and confidence intervals ",
                  "into estimated copy numbers using the user-specified initial copy number and dilution. ",
                  "The columns mirror the dilution table but are on the copy-number scale."
                ),
                
                h5("Percentile plot"),
                p(
                  "The plot shows point estimates and 95% confidence intervals for dilution across all percentiles. ",
                  "The selected percentile (e.g. 95%) is indicated by red dashed lines and labelled with its estimate and CI."
                ),
                
                h4("Interpretation"),
                tags$ul(
                  tags$li(
                    strong("LOD at a given probability (e.g. 95%): "),
                    "the estimated dilution or copy number where approximately that proportion of replicates ",
                    "would test positive under the probit model."
                  ),
                  tags$li(
                    strong("Uncertainty: "),
                    "the confidence interval around the LOD reflects uncertainty from the dose–response fit, ",
                    "not just simple binomial variability at a single dilution."
                  ),
                  tags$li(
                    strong("Use in validation reports: "),
                    "the estimated LOD (e.g. at 95% detection) and its confidence interval can be reported ",
                    "as the analytical sensitivity of the assay."
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

lodAppServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    
    getNames <- reactive({
      inFile1 <- input$file1
      if (is.null(inFile1)) return(NULL)
      theD <- read.csv(inFile1$datapath)
      theNames <- names(which(sapply(theD, is.numeric)))
      list(theNames = theNames, N = nrow(theD))
    })
    
    Data <- reactive({
      inFile1 <- input$file1
      if (is.null(inFile1)) return(NULL)
      theD <- read.csv(inFile1$datapath)
      
      table  <- matrix(NA)
      tableC <- matrix(NA)
      theDD  <- data.frame(Error = "Please make sure you select the right columns")
      index  <- NA_integer_
      
      nrowD <- nrow(theD)
      if (isTRUE(input$manD)) {
        theD$dilution <- theD[, input$D]
      } else {
        theD$dilution <- seq(input$dilution, input$dilution - nrowD + 1)
      }
      
      theD$N  <- theD[, input$N]
      theD$Y  <- theD[, input$Y]
      theD$copyNumber <- 10^theD$dilution * input$copy * 10^abs(max(theD$dilution))
      
      if (!sum(theD$N, na.rm = TRUE) < sum(theD$Y, na.rm = TRUE)) {
        theD$NY <- theD$N - theD$Y
        
        r <- glm(
          cbind(Y, NY) ~ dilution,
          data   = theD,
          family = binomial(link = "probit"),
          method = "brglmFit"
        )
        
        rr <- dose.p(
          r,
          p = c(
            seq(0.01, 0.10, by = 0.01),
            seq(0.15, 0.85, by = 0.05),
            seq(0.90, 0.99, by = 0.01)
          )
        )
        
        tableU <- as.matrix(rr) + 1.96 * attr(rr, "SE")
        tableL <- as.matrix(rr) - 1.96 * attr(rr, "SE")
        table  <- round(as.data.frame(cbind(as.matrix(rr), tableL, tableU)), 3)
        names(table) <- c("Estimate", "LowerCI", "UpperCI")
        table$label  <- as.numeric(substr(rownames(table), 5, 8))
        
        theDD <- theD[, -which(names(theD) %in% c("NY", "N", "Y"))]
        table <- table[order(table$label, decreasing = TRUE), ]
        index <- which(table$label %in% input$p)
        
        tableC <- round(10^table * input$copy * 10^abs(max(theD$dilution)), input$Ncopy)
      }
      
      list(dd = theDD, tt = table, cc = tableC, index = index)
    })
    
    output$namesN <- renderUI({
      theNames <- getNames()$theNames
      if (is.null(theNames)) return(NULL)
      theL <- as.list(theNames); names(theL) <- theNames
      selectInput(
        session$ns("N"),
        label   = h5("Select column for number of replicates"),
        choices = theL,
        selected = theL[1]
      )
    })
    
    output$namesY <- renderUI({
      theNames <- getNames()$theNames
      if (is.null(theNames)) return(NULL)
      theL <- as.list(theNames); names(theL) <- theNames
      selectInput(
        session$ns("Y"),
        label   = h5("Select column for number of positives"),
        choices = theL,
        selected = theL[2]
      )
    })
    
    output$namesD <- renderUI({
      theNames <- getNames()$theNames
      if (is.null(theNames)) return(NULL)
      theL <- as.list(theNames); names(theL) <- theNames
      selectInput(
        session$ns("D"),
        label   = h5("Select column for log dilution"),
        choices = theL,
        selected = theL[2]
      )
    })
    
    output$plot <- renderPlot({
      if (is.null(input$file1)) return(NULL)
      table <- Data()$tt
      if (is.null(table) || all(dim(table) == c(1, 1))) return(NULL)
      table <- round(table, 3)
      
      y <- table$Estimate[table$label == input$p]
      rownames(table) <-
        paste(table$Estimate, " (", table$LowerCI, ", ", table$UpperCI, ")", sep = "")
      
      ggplot(
        data = table,
        aes(
          x      = label,
          y      = Estimate,
          ymin   = LowerCI,
          ymax   = UpperCI,
          label  = rownames(table)
        )
      ) +
        geom_pointrange(
          fill    = '#60A5FA',
          color   = '#9CA3AF',
          shape   = 21,
          fatten  = 1,
          size    = 2,
          alpha   = 0.7,
          linewidth = 1
        ) +
        geom_point(
          color = '#60A5FA',
          shape = 21,
          size  = 2.5,
          alpha = 0.9
        ) +
        geom_hline(
          aes(yintercept = y),
          color   = "#F97316",
          linetype = "dashed",
          linewidth = 1
        ) +
        geom_vline(
          aes(xintercept = as.numeric(input$p)),
          color   = "#F97316",
          linetype = "dashed",
          linewidth = 1
        ) +
        coord_flip() +
        geom_text(
          data  = table[table$label == input$p, ],
          size  = 5,
          hjust = 1.5,
          vjust = 1.5,
          color = "#FBBF24",
          aes(label, Estimate, label = rownames(table[table$label == input$p, ]))
        ) +
        xlab("Percentile") +
        ylab("Estimate (95% CI)") +
        theme_minimal(base_size = 13) +
        theme(
          plot.background  = element_rect(fill = "#111827", colour = NA),
          panel.background = element_rect(fill = "#111827", colour = NA),
          text             = element_text(colour = "#E5E7EB"),
          axis.text        = element_text(colour = "#D1D5DB"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(color = "#374151")
        ) +
        scale_x_continuous(breaks = unique(table$label))
    })
    
    output$table <- DT::renderDataTable({
      if (is.null(input$file1)) {
        return(
          datatable(
            data.frame(Message = "Upload a LOD csv file to get started."),
            options = list(dom = "t"),
            rownames = FALSE
          )
        )
      }
      temp1 <- Data()$dd
      if (is.null(temp1)) temp1 <- data.frame()
      if (nrow(temp1) == 0) {
        datatable(data.frame(No_Data = NULL),
                  options = list(dom = "t"),
                  rownames = FALSE)
      } else {
        datatable(
          data.frame(temp1),
          rownames   = FALSE,
          extensions = "Buttons",
          options    = list(
            dom        = "Bfrtip",
            buttons    = c("copy", "csv", "excel"),
            pageLength = 15,
            lengthMenu = c(10, 15, 30, 50)
          )
        )
      }
    })
    
    output$tableRR <- DT::renderDataTable({
      if (is.null(input$file1)) {
        return(
          datatable(
            data.frame(Message = "Upload a LOD csv file to see dilution-based results."),
            options = list(dom = "t"),
            rownames = FALSE
          )
        )
      }
      temp  <- Data()
      temp1 <- temp$tt
      if (!is.null(temp1)) {
        temp1 <- data.frame(temp1)
        temp1$label <- NULL
      } else {
        temp1 <- data.frame()
      }
      if (nrow(temp1) == 0) {
        return(datatable(data.frame(No_Data = NULL),
                         options = list(dom = "t"),
                         rownames = FALSE))
      }
      index <- temp$index - 1
      js <- paste0(
        "function(row, data) {",
        "$(this.api().row(", index, ").node()).css({'background-color': '#1F2933'});",
        "}"
      )
      
      datatable(
        temp1,
        rownames   = TRUE,
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons    = c("copy", "csv", "excel"),
          pageLength = 35,
          lengthMenu = c(10, 35, 50, 100),
          drawCallback = JS(js)
        )
      )
    })
    
    output$tableRRC <- DT::renderDataTable({
      if (is.null(input$file1)) {
        return(
          datatable(
            data.frame(Message = "Upload a LOD csv file to see copy-number results."),
            options = list(dom = "t"),
            rownames = FALSE
          )
        )
      }
      temp  <- Data()
      temp1 <- temp$cc
      if (!is.null(temp1)) {
        temp1 <- data.frame(temp1)
        temp1$label <- NULL
      } else {
        temp1 <- data.frame()
      }
      if (nrow(temp1) == 0) {
        return(datatable(data.frame(No_Data = NULL),
                         options = list(dom = "t"),
                         rownames = FALSE))
      }
      index <- temp$index - 1
      js <- paste0(
        "function(row, data) {",
        "$(this.api().row(", index, ").node()).css({'background-color': '#1F2933'});",
        "}"
      )
      
      datatable(
        temp1,
        rownames   = TRUE,
        extensions = "Buttons",
        options    = list(
          dom        = "Bfrtip",
          buttons    = c("copy", "csv", "excel"),
          pageLength = 35,
          lengthMenu = c(10, 35, 50, 100),
          drawCallback = JS(js)
        )
      )
    })
  })
}

## ---------- Top-level UI & server (dark theme) --------

ui <- fluidPage(
  theme = bs_theme(
    version = 4,
    bootswatch = "darkly",
    base_font = font_google("Roboto"),
    heading_font = font_google("Roboto Slab"),
    bg  = "#020617",
    fg  = "#E5E7EB",
    primary   = "#3B82F6",
    secondary = "#F97316"
  ),
  
  tags$head(
    tags$style(HTML("
      body {
        background-color: #020617;
      }

      .diag-input-table th, .diag-input-table td {
        text-align: center;
        vertical-align: middle;
      }

      /* ===========================
         TOP-LEVEL APP TABS (main_tabs)
         Big pill buttons
         =========================== */
      /* ===========================
   TOP-LEVEL APP TABS (main_tabs)
   Big pill buttons, app-specific colors
   =========================== */
/* ===========================
   TOP-LEVEL APP TABS (main_tabs)
   Big pill buttons, app-specific colors
   =========================== */

/* base style: all three top-level tabs */
/* TOP-LEVEL APP TABS (main_tabs) */
/* Base style for ALL three tabs */
ul.nav[data-tabsetid='main_tabs'] > li > a {
  background-color: #111827;   /* dark pill on dark bg */
  color: #E5E7EB;
  border-radius: 999px;
  margin-right: 8px;
  border: 1px solid #4B5563;
  padding: 8px 22px;
  font-weight: 500;
  text-transform: uppercase;
  letter-spacing: 0.03em;
}

ul.nav[data-tabsetid='main_tabs'] > li > a:hover {
  background-color: #1F2937;
  color: #F9FAFB;
  border-color: #6B7280;
}

/* Tab 1 active: Ct Precision (blue pill) */
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(1).active > a,
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(1).active > a:focus,
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(1).active > a:hover {
  background-color: #2563EB;   /* blue */
  color: #F9FAFB;
  border-color: #60A5FA;
  box-shadow: 0 6px 14px rgba(37,99,235,0.6);
}

/* Tab 2 active: LOD (orange pill) */
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(2).active > a,
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(2).active > a:focus,
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(2).active > a:hover {
  background-color: #EA580C;   /* orange */
  color: #FFF7ED;
  border-color: #FDBA74;
  box-shadow: 0 6px 14px rgba(234,88,12,0.6);
}

/* Tab 3 active: Diagnostic 2×2 (teal pill) */
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(3).active > a,
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(3).active > a:focus,
ul.nav[data-tabsetid='main_tabs'] > li:nth-child(3).active > a:hover {
  background-color: #0D9488;   /* teal */
  color: #ECFEFF;
  border-color: #5EEAD4;
  box-shadow: 0 6px 14px rgba(13,148,136,0.6);
}


#main_tabs .nav-pills > li > a:hover {
  background-color: #1F2937;
  color: #F9FAFB;
  border-color: #6B7280;
}

/* Tab 1: Ct Precision (blue) */
#main_tabs .nav-pills > li:nth-child(1).active > a,
#main_tabs .nav-pills > li:nth-child(1).active > a:focus,
#main_tabs .nav-pills > li:nth-child(1).active > a:hover {
  background-color: #2563EB;   /* blue */
  color: #F9FAFB;
  border-color: #60A5FA;
  box-shadow: 0 6px 14px rgba(37,99,235,0.6);
}

/* Tab 2: LOD Estimation (orange) */
#main_tabs .nav-pills > li:nth-child(2).active > a,
#main_tabs .nav-pills > li:nth-child(2).active > a:focus,
#main_tabs .nav-pills > li:nth-child(2).active > a:hover {
  background-color: #EA580C;   /* orange */
  color: #FFF7ED;
  border-color: #FDBA74;
  box-shadow: 0 6px 14px rgba(234,88,12,0.6);
}

/* Tab 3: Diagnostic 2x2 (teal) */
#main_tabs .nav-pills > li:nth-child(3).active > a,
#main_tabs .nav-pills > li:nth-child(3).active > a:focus,
#main_tabs .nav-pills > li:nth-child(3).active > a:hover {
  background-color: #0D9488;   /* teal */
  color: #ECFEFF;
  border-color: #5EEAD4;
  box-shadow: 0 6px 14px rgba(13,148,136,0.6);
}


      /* ===========================
         INNER TABS (within each app)
         Slim, underlined style
         =========================== */

      /* Common base style for all inner tabsets */
      #ct-ct_tabs .nav-pills > li > a,
      #lod-lod_tabs .nav-pills > li > a,
      #diag2-diag_tabs .nav-pills > li > a {
        background-color: transparent;
        color: #9CA3AF;
        border-radius: 0;
        border: none;
        border-bottom: 2px solid transparent;
        margin-right: 18px;
        padding: 4px 2px 6px 2px;
        font-weight: 400;
        text-transform: none;
        letter-spacing: 0.01em;
      }

      #ct-ct_tabs .nav-pills > li > a:hover,
      #lod-lod_tabs .nav-pills > li > a:hover,
      #diag2-diag_tabs .nav-pills > li > a:hover {
        color: #E5E7EB;
        border-bottom-color: #4B5563;
        background-color: transparent;
      }

      /* Ct app inner tabs: blue accent */
      #ct-ct_tabs .nav-pills > li.active > a,
      #ct-ct_tabs .nav-pills > li.active > a:focus,
      #ct-ct_tabs .nav-pills > li.active > a:hover {
        background-color: transparent;
        color: #BFDBFE;
        border-bottom-color: #3B82F6;
        font-weight: 500;
      }

      /* LOD app inner tabs: orange accent */
      #lod-lod_tabs .nav-pills > li.active > a,
      #lod-lod_tabs .nav-pills > li.active > a:focus,
      #lod-lod_tabs .nav-pills > li.active > a:hover {
        background-color: transparent;
        color: #FED7AA;
        border-bottom-color: #F97316;
        font-weight: 500;
      }

      /* Diagnostic 2x2 inner tabs: teal accent */
      #diag2-diag_tabs .nav-pills > li.active > a,
      #diag2-diag_tabs .nav-pills > li.active > a:focus,
      #diag2-diag_tabs .nav-pills > li.active > a:hover {
        background-color: transparent;
        color: #A7F3D0;
        border-bottom-color: #14B8A6;
        font-weight: 500;
      }

      /* ===========================
         Side panels styling
         =========================== */

      .ct-panel .well {
        border-left: 5px solid #3B82F6;
        background-color: #020617;
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.6);
        border: 1px solid #1F2937;
      }
      .ct-panel h4 {
        color: #60A5FA;
        font-weight: 600;
        margin-top: 0;
      }

      .lod-panel .well {
        border-left: 5px solid #F97316;
        background-color: #020617;
        border-radius: 10px;
        box-shadow: 0 2px 6px rgba(0,0,0,0.6);
        border: 1px solid #1F2937;
      }
      .lod-panel h4 {
        color: #FB923C;
        font-weight: 600;
        margin-top: 0;
      }

      .ct-panel .form-group label,
      .lod-panel .form-group label {
        font-weight: 500;
        color: #E5E7EB;
      }

      /* ===========================
         Hero cards on top
         =========================== */
      .hero-card {
        border-radius: 12px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.7);
        border: 1px solid #1F2937;
        background: radial-gradient(circle at top left, #111827, #020617);
        color: #E5E7EB;
      }
      .hero-card h4 {
        margin-top: 0;
        margin-bottom: 6px;
        font-weight: 600;
      }
      .hero-icon {
        margin-right: 8px;
      }
      .hero-card p {
        color: #9CA3AF;
      }


      /* ===========================
         Hero cards as app navigation
         =========================== */
      .hero-card.clickable {
        cursor: pointer;
        transition: transform 120ms ease, box-shadow 120ms ease, border-color 120ms ease;
      }
      .hero-card.clickable:hover {
        transform: translateY(-2px);
        border-color: #6B7280;
      }
      .hero-card.active {
        border-color: #93C5FD;
      }
      #hero_ct.active {
        box-shadow: 0 0 0 1px rgba(59,130,246,0.35), 0 10px 22px rgba(37,99,235,0.25);
      }
      #hero_lod.active {
        box-shadow: 0 0 0 1px rgba(249,115,22,0.35), 0 10px 22px rgba(234,88,12,0.25);
      }
      #hero_diag2.active {
        box-shadow: 0 0 0 1px rgba(20,184,166,0.35), 0 10px 22px rgba(13,148,136,0.25);
      }

      /* ===========================
         Plots & text colors
         =========================== */
      .hero-card,
      .well {
        color: #E5E7EB;
      }

      /* ===========================
         DataTables / table styling
         =========================== */
      table.dataTable {
        color: #E5E7EB !important;
      }
      .dataTables_wrapper .dataTables_length,
      .dataTables_wrapper .dataTables_filter,
      .dataTables_wrapper .dataTables_info,
      .dataTables_wrapper .dataTables_processing,
      .dataTables_wrapper .dataTables_paginate {
        color: #9CA3AF !important;
      }
  ")),
    
    tags$script(HTML("  (function(){    function setActive(id){      $('#hero_ct, #hero_lod, #hero_diag2').removeClass('active');      $(id).addClass('active');    }    $(document).on('click', '#hero_ct', function(){      setActive('#hero_ct');      Shiny.setInputValue('app_choice', 'ct', {priority: 'event'});    });    $(document).on('click', '#hero_lod', function(){      setActive('#hero_lod');      Shiny.setInputValue('app_choice', 'lod', {priority: 'event'});    });    $(document).on('click', '#hero_diag2', function(){      setActive('#hero_diag2');      Shiny.setInputValue('app_choice', 'diag2', {priority: 'event'});    });    $(document).on('shiny:connected', function(){      setActive('#hero_ct');    });  })();"))
  ),
  
  
  titlePanel(
    div(
      style = "display:flex; align-items:center; gap:20px;",
      # img(
      #   src   = "PHO.png",
      #   height = 96,
      #   width  = 288,
      #   style  = "flex-shrink:0; margin-right:10px;"
      # ),
      h2("Ct Precision & LOD Tool",
         style = "margin-bottom:0; font-weight:600; color:#F9FAFB;")
    ),
    windowTitle = "Ct Precision & LOD Tool"
  ),
  
  tags$hr(style="margin-top:0; margin-bottom:15px; border-color:#1F2937;"),
  
  withMathJax(),
  
  br(),
  fluidRow(
    column(
      4,
      wellPanel(
        id = "hero_ct",
        class = "hero-card hero-ct clickable",
        h4(
          span(icon("vial", class = "hero-icon")),
          "Ct regression & variance"
        ),
        p(
          "Model Ct as a continuous outcome; decompose measurement variability into within- and between-experiment SD; ",
          "compare mixed-model and ANOVA-based variance components."
        )
      )
    ),
    column(
      4,
      wellPanel(
        id = "hero_lod",
        class = "hero-card hero-lod clickable",
        h4(
          span(icon("bullseye", class = "hero-icon")),
          "Probit LOD (binary detection)"
        ),
        p(
          "Fit probit dose–response models for detection/non-detection by dilution; ",
          "estimate LOD at chosen detection probabilities on both dilution and copy-number scales."
        )
      )
    ),
    column(
      4,
      wellPanel(
        id = "hero_diag2",
        class = "hero-card hero-diag2 clickable",
        h4(
          span(icon("table", class = "hero-icon")),
          "Diagnostic 2×2"
        ),
        p(
          "Summarize 2×2 diagnostic test data (TP, FP, TN, FN); ",
          "compute prevalence, sensitivity, specificity, predictive values, and likelihood ratios ",
          "with exact binomial confidence intervals."
        )
      )
    )
  ),
  br(),
  
  uiOutput("selected_app_ui"),
  
  tags$hr(style="border-color:#1F2937;"),
  div(
    style = "font-size:0.85em; color:#6B7280; text-align:right; margin-top:5px;",
    "Version 1.0 • For questions, contact: ",
    tags$a("Lennon Li", href = "mailto:lennon.li@oahpp.ca")
  )
)

server <- function(input, output, session) {
  selected_app <- reactiveVal("ct")
  
  observeEvent(input$app_choice, {
    req(input$app_choice)
    selected_app(as.character(input$app_choice))
  }, ignoreInit = TRUE)
  
  output$selected_app_ui <- renderUI({
    switch(
      selected_app(),
      ct    = ctAppUI("ct"),
      lod   = lodAppUI("lod"),
      diag2 = diag2x2UI("diag2"),
      ctAppUI("ct")
    )
  })
  
  ctAppServer("ct")
  lodAppServer("lod")
  diag2x2Server("diag2")
}

shinyApp(ui, server)
