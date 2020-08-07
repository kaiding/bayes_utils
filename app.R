library(shiny)
library(shinybusy)
library(markdown)
library(mosaic)
library(flextable)
library(officer)
library(Rcpp)
library(RcppDist)
library(RcppNumerical)
#library(RcppEigen)

Rcpp::sourceCpp("Bayes_Bin_Two.cpp")
Rcpp::sourceCpp("Bayes_Bin_One.cpp")
Rcpp::sourceCpp("Bayes_Normal_.cpp")
source("Bayes_Utils.R")


# Define UI for dataset viewer app ----
ui <- navbarPage(title = "Bayesian Predictive App",
                 
                 tabPanel(title = span("Introduction",
                                       style = "color: orange; font-size: 24px"),
                          tabsetPanel(
                              tabPanel("About",
                                       fluidRow(
                                           column(10,includeMarkdown("about.Rmd"))
                                       )
                              ),
                              
                              tabPanel("Bayesian Inference Graphical Illustration - Beta Binomial",
                                       sidebarLayout(
                                           sidebarPanel(
                                               p("Beta prior (in blue)", style="color: blue"),
                                               p("Binomial likelihood (in green)", style="color: green"),
                                               p("Beta posterior (in red)", style="color: red"),
                                               numericInput("alphaG", "alpha, shape parameter for Beta prior:", 
                                                            value = 1,min = 0.01, max = 100),
                                               numericInput("betaG", "beta, shape parameter for Beta prior:", 
                                                            value = 1,min = 0.01, max = 100),
                                               numericInput("rG", "Number of observed responses:", 15,
                                                            min = 1, max = 100),
                                               numericInput("nG", "Number of participants enrolled:", 20,
                                                            min = 1, max = 100),
                                               numericInput("pG", "Target response rate:", 0.5,
                                                            min = 0, max = 1),
                                               actionButton(inputId = "plotG",
                                                            label = "Draw")
                                           ),
                                           
                                           mainPanel(plotOutput("plotG"),
                                                     
                                                     textOutput("description_G"))
                                       )
                              ),
                              
                              tabPanel("Estimating beta parameters for prior",
                                       sidebarLayout(
                                           sidebarPanel(
                                               p("Beta prior B(a, b)", style="color: blue"),

                                               
                                               selectInput(inputId = "method_beta_par",
                                                           label = "Method for estimating Beta parameters, from historical data",
                                                           choices = c(
                                                               "Given mean and variance of response rates" = "1",
                                                               "Given location and interval of response rates" = "2"
                                                           )),
                                               
                                               conditionalPanel(
                                                   condition = "input.method_beta_par == '1'",
                                                   
                                                   p("What is the mean of response rates?"),
                                                   
                                                   textInput(inputId = "beta_par_1_mean",
                                                             label = "Mean:",
                                                             value = "0.4"),
                                                   p("What is the variance of response rates?"),
                                                   
                                                   textInput(inputId = "beta_par_1_var",
                                                             label = "Variance:",
                                                             value = "0.02"),
                                                   
                                                   p("Warning: (a,b) could be negative using this method!", style="color: red"),
                                                   
                                                   actionButton(inputId = "beta_par_1_GO",
                                                                label = "Calculate")
                                               ),
                                               
                                               conditionalPanel(
                                                   condition = "input.method_beta_par == '2'",
                                                   
                                                   p("What is the response rate that is mostly likely to occur? (eg., mode)"),
                                                   
                                                   textInput(inputId = "beta_par_2_mode",
                                                             label = "Mode:",
                                                             value = "0.95"),
                                                   
                                                   p("What is the chance the response rate falls between two values (L and U)? 
                                                     (eg., AUC of the Beta distribution between L and U)"),
                                                   
                                                   textInput(inputId = "beta_par_2_AUC",
                                                             label = "AUC:",
                                                             value = "0.0001"),
                                                   
                                                   textInput(inputId = "beta_par_2_L",
                                                             label = "L:",
                                                             value = "0"),
                                                   
                                                   textInput(inputId = "beta_par_2_U",
                                                             label = "U:",
                                                             value = "0.8"),
                                                   
                                                   actionButton(inputId = "beta_par_2_GO",
                                                                label = "Calculate")
                                               )
                                           ),
                                           
                                           mainPanel(
                                               conditionalPanel(
                                                   condition = "input.method_beta_par == '1'",
                                                   withMathJax(helpText("Solving for $$Mean=\\frac{a}{a+b}$$", 
                                                                        "and,",
                                                                        "$$Variance=\\frac{ab}{\\left(a+b\\right)^2
                                                                        \\left(a+b+1\\right)}$$")),
                                                   uiOutput("beta_par_1")
                                                   ),
                                               
                                               conditionalPanel(
                                                   condition = "input.method_beta_par == '2'",
                                                   withMathJax(helpText("Solving for $$Mode=\\frac{a-1}{a+b-2}$$", 
                                                                        "and,",
                                                                        "$$AUC=\\int_L^U\\frac{\\Gamma(a+b)}{\\Gamma(a)\\Gamma(b)}
                                                                        r^{a-1}(1-r)^{b-1}dr$$")),
                                                   uiOutput("beta_par_2")
                                               )
                                                   )
                                           )
                                       )
                              ,
                              
                              tabPanel("Test",
                                       tabsetPanel(
                                           tabPanel("Under Development"),
                                           tabPanel("Under Development")
                                       )
                                       )
                              )
      
                 ),
                 
                 tabPanel(title = span("Interim Monitoring",
                                       style = "color: orange; font-size: 24px"),
                          tabsetPanel(
                              tabPanel(title = "Binomial Endpoint",
                                       tabsetPanel(
                                           tabPanel(title = "One-Sample",
                                                    sidebarLayout(
                                                      sidebarPanel(
                                                        textInput(inputId = "sampleSize_IA_Bin_One",
                                                                  label = "Sample size per arm",
                                                                  value = "50"),
                                                        
                                                        textInput(inputId = "delta_IA_Bin_One",
                                                                  label = "Clinically important difference (eg., p1 - p0 > delta)",
                                                                  value = "0"),
                                                    
                                                        textInput(inputId = "neta_IA_Bin_One",
                                                                  label = "Threshold of posterior probability for declaring efficacy 
                                                   (eg., Pr(p1 - p0 > delta|X) > neta)",
                                                                  value = "0.95"),
                                                        
                                                        textInput(inputId = "alpha_IA_Bin_One",
                                                                  label = "One-sided alpha for predictive/conditional power approach",
                                                                  value = "0.05"),
                                                        
                                                        textInput(inputId = "p0_IA_Bin_One",
                                                                  label = "Hypothesized proportion under the null hypothesis",
                                                                  value = "0.55"),
                                                        
                                                        textInput(inputId = "p1_IA_Bin_One",
                                                                  label = "Hypothesized proportion under the alternative hypothesis",
                                                                  value = "0.65"),
                                                        
                                                        textInput(inputId = "a_IA_Bin_One",
                                                                  label = "Parameter a for Beta prior B(a, b)
                                                                  (semicolon-delimited for different priors.)",
                                                                  value = "1;5"),
                                                        
                                                        textInput(inputId = "b_IA_Bin_One",
                                                                  label = "Parameter b for Beta prior B(a, b)
                                                                  (semicolon-delimited for different priors.)",
                                                                  value = "1;5"),

                                                        textInput(inputId = "n_IA_Bin_One",
                                                                  label = "Number of subjects at interim",
                                                                  value = "25"),
                                                        

                                                        textInput(inputId = "r_IA_Bin_One",
                                                                  label = "Number of responders  at interim",
                                                                  value = "17"),
                                                    
                                                        actionButton(inputId = "go_IA_Bin_One",
                                                                     label = "Calculate")
                                                      ),
                                                      mainPanel(
                                                        textOutput("description_IA_Bin_One"),
                                                        
                                                        tabPanel("Table",
                                                                 uiOutput("res_IA_Bin_One")
                                                        )
                                                      )
                                                    )
                                                    ),
                                           
                                           tabPanel(title = "Two-Sample",
                                                    sidebarLayout(
                                                        sidebarPanel(
                                                            textInput(inputId = "sampleSize_IA_Bin",
                                                                      label = "Sample size per arm",
                                                                      value = "212"),
                                                            
                                                            textInput(inputId = "delta_IA_Bin",
                                                                      label = "Clinically important difference (eg., p1 - p2 > delta)",
                                                                      value = "0"),
                                                            
                                                            textInput(inputId = "neta_IA_Bin",
                                                                      label = "Threshold of posterior probability for declaring efficacy 
                                                   (eg., Pr(p1 - p2 > 0|X) > neta)",
                                                                      value = "0.9"),
                                                            
                                                            textInput(inputId = "alpha_IA_Bin",
                                                                      label = "One-sided alpha for predictive/conditional power approach",
                                                                      value = "0.1"),
                                                            
                                                            textInput(inputId = "es_IA_Bin",
                                                                      label = "Effective size under althernative hypothesis",
                                                                      value = "0.2065"),
                                                            
                                                            textInput(inputId = "p1_IA_Bin",
                                                                      label = "Prior rate in group 1 
                                                   (Specify if Beta priori is not entered. Semicolon-delimited for different priors.
                                                         Beta parameters will be used in case both are entered.)",
                                                                      value = "0.5;0.79;0.79;0.7"),
                                                            
                                                            textInput(inputId = "p2_IA_Bin",
                                                                      label = "Hypothesized rate in group 2 
                                                   (Specify if Beta priori is not entered. Semicolon-delimited for different priors.
                                                         Beta parameters will be used in case both are entered.)",
                                                                      value = "0.5;0.7;0.7;0.7"),
                                                            
                                                            textInput(inputId = "tau_IA_Bin",
                                                                      label = "Hypothesized common coefficient of variation in both arms 
                                        (Specify if Beta priori is not entered. Semicolon-delimited for different priors.
                                                         Beta parameters will be used in case both are entered.)",
                                                                      value = "0.5773503;0.1;0.05;0.05"),
                                                            
                                                            textInput(inputId = "a1_IA_Bin",
                                                                      label = "Parameter a for Beta priori B(a, b) in group 1
                                        (Specify if group 1 rate is not entered, semicolon-delimited for different priors.)",
                                                                      value = ""),
                                                            
                                                            textInput(inputId = "b1_IA_Bin",
                                                                      label = "Parameter b for Beta priori B(a, b) in group 1
                                        (Specify if group 1 rate is not entered, semicolon-delimited for different priors.)",
                                                                      value = ""),
                                                            
                                                            textInput(inputId = "a2_IA_Bin",
                                                                      label = "Parameter a for Beta priori B(a, b) in group 2 
                                        (Specify if group 2 rate is not entered, semicolon-delimited for different priors.)",
                                                                      value = ""),
                                                            
                                                            textInput(inputId = "b2_IA_Bin",
                                                                      label = "Parameter b for Beta priori B(a, b) in group 2 
                                        (Specify if group 2 rate is not entered, semicolon-delimited for different priors.)",
                                                                      value = ""),
                                                            
                                                            textInput(inputId = "n1_IA_Bin",
                                                                      label = "Number of subjects in group 1 at interim",
                                                                      value = "55"),
                                                            
                                                            textInput(inputId = "n2_IA_Bin",
                                                                      label = "Number of subjects in group 2 at interim",
                                                                      value = "45"),
                                                            
                                                            textInput(inputId = "r1_IA_Bin",
                                                                      label = "Number of responders in group 1 at interim",
                                                                      value = "33"),
                                                            
                                                            textInput(inputId = "r2_IA_Bin",
                                                                      label = "Number of responders in group 2 at interim",
                                                                      value = "29"),
                                                            
                                                            actionButton(inputId = "go_IA_Bin",
                                                                         label = "Calculate")
                                                        ),
                                                        mainPanel(
                                                            textOutput("description_IA_Bin"),
                                                            
                                                            tabPanel("Table",
                                                                     uiOutput("res_IA_Bin")
                                                            )
                                                        )
                                                    )
                                                    )
                                           )
                                       ),
                              
                              tabPanel(title = "Normal Endpoint",
                                       tabsetPanel(
                                           tabPanel(title = "Under Development"),
                                           
                                           tabPanel(title = "Two-Sample",
                                                    sidebarLayout(
                                                        sidebarPanel(
                                                            textInput(inputId = "sampleSize_IA_Norm",
                                                                      label = "Sample size per group",
                                                                      value = "41"),
                                                            
                                                            textInput(inputId = "delta_IA_Norm",
                                                                      label = "Clinically important difference (eg., mu1 - mu2 > delta)",
                                                                      value = "0"),
                                                            
                                                            textInput(inputId = "neta_IA_Norm",
                                                                      label = "Threshold of posterior probability for declaring efficacy 
                                                   (eg., Pr(mu1 - mu2 > delta|X) > neta)",
                                                                      value = "0.9"),
                                                            
                                                            textInput(inputId = "alpha_IA_Norm",
                                                                      label = "One-sided alpha for predictive/conditional power",
                                                                      value = "0.1"),
                                                            
                                                            textInput(inputId = "es_IA_Norm",
                                                                      label = "Effective size under althernative hypothesis",
                                                                      value = "0.4706"),
                                                            
                                                            textInput(inputId = "mu1_IA_Norm",
                                                                      label = "Prior means in group 1, semicolon-delimited for different priors.",
                                                                      value = "12;12;12;8"),
                                                            
                                                            textInput(inputId = "mu2_IA_Norm",
                                                                      label = "Prior means in group 2, semicolon-delimited for different priors.",
                                                                      value = "8;8;8;8"),
                                                            
                                                            textInput(inputId = "SD_IA_Norm",
                                                                      label = "Prior common standard deviations, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "1000;10;2;2"),
                                                            
                                                            textInput(inputId = "n1_IA_Norm",
                                                                      label = "Number of subjects in group 1 at interim",
                                                                      value = "10"),
                                                            
                                                            textInput(inputId = "n2_IA_Norm",
                                                                      label = "Number of subjects in group 2 at interim",
                                                                      value = "11"),
                                                            
                                                            textInput(inputId = "x1bar_IA_Norm",
                                                                      label = "Means in group 1 at interim",
                                                                      value = "9.2"),
                                                            
                                                            textInput(inputId = "sd1_IA_Norm",
                                                                      label = "SD's in group 1 at interim",
                                                                      value = "7.3"),
                                                            
                                                            textInput(inputId = "x2bar_IA_Norm",
                                                                      label = "Means in group 2 at interim",
                                                                      value = "8.4"),
                                                            
                                                            textInput(inputId = "sd2_IA_Norm",
                                                                      label = "SD's in group 2 at interim",
                                                                      value = "6.4"),
                                                            
                                                            selectInput(inputId = "methodpb_IA_Norm",
                                                                        label = "Method for Bayesian Predictive Probability",
                                                                        choices = c(
                                                                            "Analytical" = "1",
                                                                            "Simulation" = "2"
                                                                        )),
                                                            
                                                            textInput(inputId = "nsim_IA_Norm",
                                                                      label = "Number of replica for simulation-based Bayesian Predictive Probability",
                                                                      value = "10000"),
                                                            
                                                            actionButton(inputId = "go_IA_Norm",
                                                                         label = "Calculate")
                                                        ),
                                                        mainPanel(
                                                            textOutput("description_IA_Norm"),
                                                            
                                                            tabPanel("Table",
                                                                     uiOutput("res_IA_Norm")
                                                            )
                                                        )
                                                    ))
                                       )
                                       )
                          )
                          ),
                 
                 
                 tabPanel(title = span("Operating Characteristics",
                                       style = "color: orange;  font-size: 24px"),
                          tabsetPanel(
                              tabPanel(title = "Binomial Endpoint",
                                       tabsetPanel(
                                           tabPanel(title = "One-Sample",
                                                    sidebarLayout(
                                                      sidebarPanel(
                                                        textInput(inputId = "sampleSize_OC_Bin_One",
                                                                  label = "Planned sample size",
                                                                  value = "40"),
                                                        
                                                        textInput(inputId = "n_OC_Bin_One",
                                                                  label = "Planned number of subjects 
                                                                  across each interims (semicolon-delimited)",
                                                                  value = "10;20;30"),
                                                        
                                                        textInput(inputId = "as_OC_Bin_One",
                                                                  label = 'Parameter a for Beta priori B(a, b) in simulation',
                                                                  value = "1"),
                                                        
                                                        textInput(inputId = "bs_OC_Bin_One",
                                                                  label = 'Parameter b for Beta priori B(a, b) in simulation',
                                                                  value = "1"),
                                                    
                                                        textInput(inputId = "nsim_OC_Bin_One",
                                                                  label = "Number of replicas for evaluation",
                                                                  value = "1000"),
                                                        
                                                        textInput(inputId = "delta_OC_Bin_One",
                                                                  label = "Clinically important difference (eg., p1 - p0 > delta)",
                                                                  value = "0"),
                                                        
                                                        textInput(inputId = "neta_OC_Bin_One",
                                                                  label = "Threshold of posterior probability for declaring efficacy 
                                                   (eg., Pr(p1 - p0 > delta|X) > neta)",
                                                                  value = "0.95"),
                                                        
                                                        textInput(inputId = "tau_OC_Bin_One",
                                                                  label = "Threshold for declaring futility 
                                                         (eg., declaring futility if predictive probabilit below it at interim)",
                                                                  value = "0.2"),
                                                        
                                                        textInput(inputId = "alpha_OC_Bin_One",
                                                                  label = "One-sided alpha for predictive/conditional power",
                                                                  value = "0.05"),
                                                        
                                                        textInput(inputId = "p0_OC_Bin_One",
                                                                  label = "Hypothesized proportion under the null hypothesis",
                                                                  value = "0.3"),
                                                        
                                                        textInput(inputId = "p1_OC_Bin_One",
                                                                  label = "Hypothesized proportion under the alternative hypothesis",
                                                                  value = "0.5"),
                                                        
                                                        textInput(inputId = "a_OC_Bin_One",
                                                                  label = "Prior a for Beta priori B(a, b) for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                  value = "1;5"),
                                                        
                                                        textInput(inputId = "b_OC_Bin_One",
                                                                  label = "Prior b for Beta priori B(a, b) for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                  value = "1;5"),

                                                        actionButton(inputId = "go_OC_Bin_One",
                                                                     label = "Simulate")
                                                      ),
                                                      mainPanel(
                                                        textOutput("description_OC_Bin_One"),
                                                        
                                                        tabPanel("Table",
                                                                 uiOutput("res_OC_Bin_One")
                                                        )
                                                      )
                                                      )
                                                    ),
                                           
                                           tabPanel(title = "Two-Sample",
                                                    sidebarLayout(
                                                        sidebarPanel(
                                                            textInput(inputId = "sampleSize_OC_Bin",
                                                                      label = "Planned sample size per group",
                                                                      value = "40"),
                                                            
                                                            textInput(inputId = "n_OC_Bin",
                                                                      label = "Planned number of subjects per group 
                                                   across each interims (semicolon-delimited)",
                                                                      value = "10;20;30"),
                                                            
                                                            textInput(inputId = "a1s_OC_Bin",
                                                                      label = 'Parameter a for Beta priori B(a, b) in simulation for group 1',
                                                                      value = "1"),
                                                            
                                                            textInput(inputId = "b1s_OC_Bin",
                                                                      label = 'Parameter b for Beta priori B(a, b) in simulation for group 1',
                                                                      value = "1"),
                                                            
                                                            textInput(inputId = "a2s_OC_Bin",
                                                                      label = 'Parameter a for Beta priori B(a, b) in simulation for group 2',
                                                                      value = "1"),
                                                            
                                                            textInput(inputId = "b2s_OC_Bin",
                                                                      label = 'Parameter b for Beta priori B(a, b) in simulation for group 2',
                                                                      value = "1"),
                                                            
                                                            
                                                            textInput(inputId = "nsim_OC_Bin",
                                                                      label = "Number of replicas for evaluation",
                                                                      value = "1000"),
                                                            
                                                            textInput(inputId = "delta_OC_Bin",
                                                                      label = "Clinically important difference (eg., p1 - p2 > delta)",
                                                                      value = "0"),
                                                            
                                                            textInput(inputId = "neta_OC_Bin",
                                                                      label = "Threshold of posterior probability for declaring efficacy 
                                                   (eg., Pr(p1 - p2 > delta|X) > neta)",
                                                                      value = "0.9"),
                                                            
                                                            textInput(inputId = "tau_OC_Bin",
                                                                      label = "Threshold for declaring futility 
                                                         (eg., declaring futility if predictive probabilit below it at interim)",
                                                                      value = "0.2"),
                                                            
                                                            textInput(inputId = "alpha_OC_Bin",
                                                                      label = "One-sided alpha for predictive/conditional power",
                                                                      value = "0.1"),
                                                            
                                                            textInput(inputId = "es_OC_Bin",
                                                                      label = "Effective size under althernative hypothesis",
                                                                      value = "0.2076"),
                                                            
                                                            textInput(inputId = "a1_OC_Bin",
                                                                      label = "Prior a for Beta priori B(a, b) in group 1 for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "1;20"),
                                                            
                                                            textInput(inputId = "b1_OC_Bin",
                                                                      label = "Prior b for Beta priori B(a, b) in group 1 for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "1;5"),
                                                            
                                                            textInput(inputId = "a2_OC_Bin",
                                                                      label = "Prior a for Beta priori B(a, b) in group 1 for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "1;29"),
                                                            
                                                            textInput(inputId = "b2_OC_Bin",
                                                                      label = "Prior b for Beta priori B(a, b) in group 1 for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "1;12"),
                                                            
                                                            selectInput(inputId = "methodpb_OC_Bin",
                                                                        label = "Method for Bayesian Predictive Probability",
                                                                        choices = c(
                                                                            "Analytical" = "1",
                                                                            "Simulation" = "2"
                                                                        )),
                                                            
                                                            textInput(inputId = "nsim_p_OC_Bin",
                                                                      label = "Number of replica for simulation-based Bayesian Predictive Probability",
                                                                      value = "1000"),
                                                            
                                                            actionButton(inputId = "go_OC_Bin",
                                                                         label = "Simulate")
                                                        ),
                                                        
                                                        mainPanel(
                                                            textOutput("description_OC_Bin"),
                                                            
                                                            tabPanel("Table",
                                                                     uiOutput("res_OC_Bin")
                                                            )
                                                            
                                                        )
                                                    )
                                                    )
                                           )
                                       ),
                              
                              tabPanel(title = "Normal Endpoint",
                                       tabsetPanel(
                                           tabPanel(title = "Under Development"),
                                           
                                           tabPanel(title = "Two-Sample",
                                                    sidebarLayout(
                                                        sidebarPanel(
                                                            textInput(inputId = "sampleSize_OC_Norm",
                                                                      label = "Planned sample size per group",
                                                                      value = "41"),
                                                            
                                                            textInput(inputId = "n_OC_Norm",
                                                                      label = "Planned number of subjects per group 
                                                   across each interims (semicolon-delimited)",
                                                                      value = "10;20;30"),
                                                            
                                                            textInput(inputId = "mu1s_OC_Norm",
                                                                      label = 'Prior mean in group 1',
                                                                      value = "12"),
                                                            
                                                            textInput(inputId = "mu2s_OC_Norm",
                                                                      label = "Prior mean in group 2",
                                                                      value = "8"),
                                                            
                                                            textInput(inputId = "ss_OC_Norm",
                                                                      label = "Prior common standard deviation",
                                                                      value = "2"),
                                                            
                                                            textInput(inputId = "sdcom_OC_Norm",
                                                                      label = "Common standard deviation for simulation",
                                                                      value = "8.499788"),
                                                            
                                                            textInput(inputId = "nsim_OC_Norm",
                                                                      label = "Number of replicas for evaluation",
                                                                      value = "1000"),
                                                            
                                                            textInput(inputId = "delta_OC_Norm",
                                                                      label = "Clinically important difference (eg., mu1 - mu2 > delta)",
                                                                      value = "0"),
                                                            
                                                            textInput(inputId = "neta_OC_Norm",
                                                                      label = "Threshold of posterior probability for declaring efficacy 
                                                   (eg., Pr(mu1 - mu2 > delta|X) > neta)",
                                                                      value = "0.9"),
                                                            
                                                            textInput(inputId = "tau_OC_Norm",
                                                                      label = "Threshold for declaring futility 
                                                         (eg., declaring futility if predictive probabilit below it at interim)",
                                                                      value = "0.2"),
                                                            
                                                            textInput(inputId = "alpha_OC_Norm",
                                                                      label = "One-sided alpha for predictive/conditional power",
                                                                      value = "0.1"),
                                                            
                                                            textInput(inputId = "es_OC_Norm",
                                                                      label = "Effective size under althernative hypothesis",
                                                                      value = "0.4706"),
                                                            
                                                            textInput(inputId = "mu1_OC_Norm",
                                                                      label = "Prior means in group 1 for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "12;12;12;8"),
                                                            
                                                            textInput(inputId = "mu2_OC_Norm",
                                                                      label = "Prior means in group 2 for power/probability calculation,
                                                         semicolon-delimited for different priors.",
                                                                      value = "8;8;8;8"),
                                                            
                                                            textInput(inputId = "SD_OC_Norm",
                                                                      label = "Prior common standard deviations for power/probability calculation, 
                                                         semicolon-delimited for different priors.",
                                                                      value = "1000;10;2;2"),
                                                            
                                                            selectInput(inputId = "methodpb_OC_Norm",
                                                                        label = "Method for Bayesian Predictive Probability",
                                                                        choices = c(
                                                                            "Analytical" = "1",
                                                                            "Simulation" = "2"
                                                                        )),
                                                            
                                                            textInput(inputId = "nsim_p_OC_Norm",
                                                                      label = "Number of replica for simulation-based Bayesian Predictive Probability",
                                                                      value = "1000"),
                                                            
                                                            actionButton(inputId = "go_OC_Norm",
                                                                         label = "Simulate")
                                                        ),
                                                        
                                                        mainPanel(
                                                            textOutput("description_OC_Norm"),
                                                            
                                                            tabPanel("Table",
                                                                     uiOutput("res_OC_Norm")
                                                            )
                                                            
                                                        )
                                                    ))
                                       )
                                       
                              )
                          )
                 ),
                 
                 tabPanel(title = span("Under Development",
                                       style = "color: orange;  font-size: 24px"))
                 )


# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
    
    ###Bayesian Graphica Illustration
    
    
    drawG <- eventReactive(input$plotG, {
        pG = input$pG
        
        priorG = list(input$alphaG, input$betaG)
        posteriorG = list((input$alphaG + input$rG), (input$betaG + input$nG - input$rG))
        
        xG = seq(0, 1, 0.01)
        max1 = max(dbeta(xG, priorG[[1]], priorG[[2]]))
        max2 = max(dbeta(xG, posteriorG[[1]], posteriorG[[2]]))
        maxy = min(8, max(max1, max2) + 0.1)
        
        # pesky detail to rescale the likelihood
        maxbinom = dbinom(input$rG, input$nG, input$rG/input$nG)
        mle = makeFun(0.99*maxy/maxbinom*choose(nG, rG)*p^rG*(1-p)^(nG-rG) ~ p, nG=input$nG, rG=input$rG)
        
        # display the results
        
        l1 = stat_function(fun = function(x) dbeta(x, posteriorG[[1]], posteriorG[[2]]), 
                          col = "red", lwd = 1.5)
        
        l2 = stat_function(fun = function(x) dbeta(x, posteriorG[[1]], posteriorG[[2]]), 
                          geom = "area", xlim = c(pG, 1), col = "red", fill = "red", alpha = 0.5)
            
        l3 = stat_function(fun = function(x) dbeta(x,priorG[[1]], priorG[[2]]), 
                          col = "blue", lwd = 1, alpha = 0.6)
        
        l4 =  stat_function(fun = function(x) mle(x), col = "green", lwd = 1, alpha = 0.6) 
        
        l5 = geom_vline(xintercept = pG, col = "black", lty = 2, lwd = 1.2)
        
        l6 = geom_text(aes(x = pG, y = maxy/2, label = "target response rate"),
                       col = "red", angle = 90, vjust = -0.4)
        
        ggplot(data.frame(x = c(0, 1)), ylim = c(0, maxy), aes(x)) + 
            theme_bw() + l1 + l2 + l3 + l4 + l5 + l6
    })
    
    output$plotG = renderPlot({drawG()})
    
    description_G = eventReactive(input$plotG, {
        pG = input$pG
        
        priorG = list(input$alphaG, input$betaG)
        posteriorG = list((input$alphaG + input$rG), (input$betaG + input$nG - input$rG))
        postPR = 1 - pbeta(pG, posteriorG[[1]], posteriorG[[2]])
        
        d1 = paste0("With prior Beta(", 
                    priorG[[1]], ",", priorG, "), and  ", 
                    input$rG, " responder(s) out of ", input$nG, " subjects, the posterior is Beta(",
                    posteriorG[[1]], ",", posteriorG[[2]], "). ",
                    "The posterior probability that the true response rate is greater than target rate (",
                    pG, ") is ", postPR, ".")
        
        d1
    })
    
    output$description_G = renderText({description_G()})
    
    ###Estimating beta parameters
    beta_par_1 = eventReactive(input$beta_par_1_GO,{
        m = as.numeric(input$beta_par_1_mean)
        s2 = as.numeric(input$beta_par_1_var)
        res_1 = Beta_par_method1(m, s2)
        
        res_beta_par_1 = withMathJax(
            sprintf("The estimated parameters are $$a =%f$$ and $$b =%f$$", res_1[1], res_1[2])
            )
        
        res_beta_par_1
    })
    
    output$beta_par_1 = renderUI({beta_par_1()})
    
    beta_par_2 = eventReactive(input$beta_par_2_GO,{
        MODE = as.numeric(input$beta_par_2_mode)
        AUC = as.numeric(input$beta_par_2_AUC)
        L = as.numeric(input$beta_par_2_L)
        U = as.numeric(input$beta_par_2_U)
        res_2 = Beta_par_method2(MODE, AUC, L, U)
        
        res_beta_par_2 = withMathJax(
            sprintf("The estimated parameters are $$a =%f$$ and $$b =%f$$", res_2[1], res_2[2])
        )
        
        res_beta_par_2
    })
    
    output$beta_par_2 = renderUI({beta_par_2()})
    
    ###Interim Monitoring - Binary - One-Sample
    description_IA_Bin_One = eventReactive(input$go_IA_Bin_One, {
      d1 = paste0("One-sample trial with binary endpoint, and ", 
                  input$sampleSize_IA_Bin_One, " planned subjects to enroll. ")
      
      d2 = paste0("Hypothesized proportions under null and alternative hypotheses are ", 
                  input$p0_IA_Bin_One, " and ", input$p1_IA_Bin_One,  ", respectively. ")
      
      d3 = paste0("Clinical meaningful treatment difference is set to ", input$delta_IA_Bin_One, ". ")
      
      d4 = paste0("Threshold of posterior probability for declaring efficacy is set to ", input$neta_IA_Bin_One, ". ")
      
      d5 = paste0("One-sided alpha for conditional power and predictive 
                power approaches is set to ", input$alpha_IA_Bin_One, ". ")
      
      d6 = paste0("At interim, there are ", input$r_IA_Bin_One, " responders out of ", input$n_IA_Bin_One, 
                  " subjects.")
      paste0(d1, d2, d3, d4, d5, d6)
    })
    
    output$description_IA_Bin_One = renderText({description_IA_Bin_One()})
    
    calculation_IA_Bin_One = eventReactive(input$go_IA_Bin_One, {
      show_modal_spinner()
      N = as.numeric(input$sampleSize_IA_Bin_One)
      neta = as.numeric(input$neta_IA_Bin_One)
      alpha = as.numeric(input$alpha_IA_Bin_One)
      p0 = as.numeric(input$p0_IA_Bin_One)
      p1 = as.numeric(input$p1_IA_Bin_One)
      delta = as.numeric(input$delta_IA_Bin_One)

      a = as.numeric(stringr::str_split(input$a_IA_Bin_One, ";")[[1]])
      b = as.numeric(stringr::str_split(input$b_IA_Bin_One, ";")[[1]])

      beta = cbind(a, b)
      
      Prior = apply(round(beta, digits = 3), 1, Beta_string)

      n = as.numeric(input$n_IA_Bin_One)
      r = as.numeric(input$r_IA_Bin_One)
      
      res_IA_Bin_One = IA_Bin_One_cpp(n, r, N, neta, p1, p0, alpha, delta, a, b)
      
      remove_modal_spinner()
      
      nrow_IA_Binpb_One = max(3, length(a))
      nrow_IA_Binpp_One = max(2, length(a))
      
      res_IA_Bin_pb_One = cbind(c("Predictive Probability", paste0("delta=", input$delta_IA_Bin_One), 
                                  paste0("neta=", input$neta_IA_Bin_One), rep("", nrow_IA_Binpb_One-3)),
                            c(Prior, rep("", nrow_IA_Binpb_One-length(a))), 
                            c(round(res_IA_Bin_One$PredProb, digits = 3),
                                  rep("", nrow_IA_Binpb_One-length(a))))
      
      res_IA_Bin_pp_One = cbind(c("Predictive Power", paste0("alpha=", input$alpha_IA_Bin_One), 
                                  rep("", nrow_IA_Binpp_One-2)),
                                c(Prior, rep("", nrow_IA_Binpp_One-length(a))), 
                                c(round(res_IA_Bin_One$PredPower, digits = 3),
                                  rep("", nrow_IA_Binpp_One-length(a))))
      

      
      res_IA_Bin_cp_One = cbind(c("Conditional Power", paste0("alpha=", input$alpha_IA_Bin_One), ""), 
                            rep("", 3), c(round(res_IA_Bin_One$CondPower, digits = 3), "", ""))
      
      
      res_IA_Bin_One = rbind(res_IA_Bin_pb_One, res_IA_Bin_pp_One, res_IA_Bin_cp_One)
      
      colnames(res_IA_Bin_One) = c("Method", "Prior", 
                               "Probability of a positive trial outcome given interim")
      
      parts_b_IA_Bin_One = fp_border(width = 3)
      
      as.data.frame(res_IA_Bin_One) %>%
        flextable() %>%
        theme_booktabs() %>%
        hline(i = c(nrow_IA_Binpb_One, nrow_IA_Binpb_One + nrow_IA_Binpp_One), border = parts_b_IA_Bin_One) %>%
        autofit() %>%
        fix_border_issues() %>%
        htmltools_value()
      
    })
    
    output$res_IA_Bin_One = renderUI({
      calculation_IA_Bin_One()
    })
    
    ###Interim Monitoring - Binary - Two-Sample
    description_IA_Bin = eventReactive(input$go_IA_Bin, {
        d1 = paste0("Two-sample trial with binary endpoint, and ", 
                    input$sampleSize_IA_Bin, " subjects per group. ")
        d2 = paste0("Clinical meaningful treatment difference is set to ", input$delta_IA_Bin, ". ")
        d3 = paste0("Threshold of posterior probability for declaring efficacy is set to ", input$neta_IA_Bin, ". ")
        d4 = paste0("One-sided alpha for conditional power and predictive 
                power approaches is set to ", input$alpha_IA_Bin, ". ")
        
        d5 = paste0("At interim, there are ", input$r1_IA_Bin, " responders out of ", input$n1_IA_Bin, 
                    " subjects in group 1, and ", input$r2_IA_Bin, " responders out of ", input$n2_IA_Bin, " subjects in group 2.")
        paste0(d1, d2, d3, d4, d5)
    })
    
    output$description_IA_Bin = renderText({description_IA_Bin()})
    
    calculation_IA_Bin = eventReactive(input$go_IA_Bin, {
        show_modal_spinner()
        N = as.numeric(input$sampleSize_IA_Bin)
        delta = as.numeric(input$delta_IA_Bin)
        neta = as.numeric(input$neta_IA_Bin)
        alpha = as.numeric(input$alpha_IA_Bin)
        es = as.numeric(input$es_IA_Bin)
        if(input$a1_IA_Bin == "" | input$b1_IA_Bin =="" |input$a2_IA_Bin == "" | input$b2_IA_Bin ==""){
            p1 = as.numeric(stringr::str_split(input$p1_IA_Bin, ";")[[1]])
            p2 = as.numeric(stringr::str_split(input$p2_IA_Bin, ";")[[1]])
            tau = as.numeric(stringr::str_split(input$tau_IA_Bin, ";")[[1]])
            beta1 = Beta_par_convert(p1, tau)
            beta2 = Beta_par_convert(p2, tau)
            a1 = beta1[,1]
            b1 = beta1[,2]
            a2 = beta2[,1]
            b2 = beta2[,2]
        }else{
            a1 = as.numeric(stringr::str_split(input$a1_IA_Bin, ";")[[1]])
            b1 = as.numeric(stringr::str_split(input$b2_IA_Bin, ";")[[1]])
            a2 = as.numeric(stringr::str_split(input$a2_IA_Bin, ";")[[1]])
            b2 = as.numeric(stringr::str_split(input$b2_IA_Bin, ";")[[1]])
            beta1 = cbind(a1, b1)
            beta2 = cbind(a2, b2)
        }
        Prior1 = apply(round(beta1, digits = 3), 1, Beta_string)
        Prior2 = apply(round(beta2, digits = 3), 1, Beta_string)
        
        n1 = as.numeric(input$n1_IA_Bin)
        n2 = as.numeric(input$n2_IA_Bin)
        r1 = as.numeric(input$r1_IA_Bin)
        r2 = as.numeric(input$r2_IA_Bin)
        
        res_IA_Bin = IA_Bin_cpp(n1, n2, r1, r2, N, 1, 1000, 
                                delta, neta, es, alpha, 
                                a1, b1, a2, b2)
        
        
        remove_modal_spinner()
        
        nrow_IA_Binpb = max(3, length(a1))
        nrow_IA_Binpp = max(2, length(a1))
        
        res_IA_Bin_pb = cbind(c("Predictive Probability", paste0("delta=", input$delta_IA_Bin),
                                 paste0("neta=", input$neta_IA_Bin), rep("", nrow_IA_Binpb-3)),
                              c(Prior1, rep("", nrow_IA_Binpb-length(a1))),
                              c(Prior2, rep("", nrow_IA_Binpb-length(a1))),
                              c(round(res_IA_Bin$PredProb, digits = 3),
                                rep("", nrow_IA_Binpb-length(a1))))
        
        res_IA_Bin_pp = cbind(c("Predictive Power", paste0("alpha=", input$alpha_IA_Bin),
                                 rep("", nrow_IA_Binpp-2)),
                              c(Prior1, rep("", nrow_IA_Binpp-length(a1))), 
                              c(Prior2, rep("", nrow_IA_Binpp-length(a1))), 
                              c(round(res_IA_Bin$PredPower, digits = 3),
                                rep("", nrow_IA_Binpp-length(a1))))
        
        res_IA_Bin_cp = cbind(c("Conditional Power", paste0("alpha=", input$alpha_IA_Bin), 
                                 paste0("effect size under alternative: ", input$es_IA_Bin)), 
                               rep("", 3), rep("", 3), rbind(round(res_IA_Bin$CondPower, digits = 3),
                                                 matrix(rep("", 2*length(res_IA_Bin$CondPower)), nrow = 2)))
        
        res_IA_Bin = rbind(res_IA_Bin_pb, res_IA_Bin_pp, res_IA_Bin_cp)
        
        colnames(res_IA_Bin) = c("Method", "Group 1", "Group 2", 
                                 "Probability of a positive trial outcome given interim")
        
        parts_b_IA_Bin = fp_border(width = 3)
        
        as.data.frame(res_IA_Bin) %>%
            flextable() %>%
            add_header_row(values = c("",
                                      rep("Priors", 2),
                                      "")) %>%
            merge_at(i=1, j=2:3, part = "header") %>%
            theme_booktabs() %>%
            hline(i = c(nrow_IA_Binpb, nrow_IA_Binpb+nrow_IA_Binpp), border = parts_b_IA_Bin) %>%
            autofit() %>%
            fix_border_issues() %>%
            htmltools_value()
    })
    
    output$res_IA_Bin = renderUI({
        calculation_IA_Bin()
    })
    
    ###Interim Monitoring - Normal
    description_IA_Norm = eventReactive(input$go_IA_Norm, {
        d1 = paste0("Two-sample trial with normal endpoint, and ", 
                    input$sampleSize_IA_Norm, " subjects per group. ")
        d2 = paste0("Clinical meaningful treatment difference is set to ", input$delta_IA_Norm, ". ")
        d3 = paste0("Threshold of posterior probability for declaring efficacy is set to ", input$neta_IA_Norm, ". ")
        d4 = paste0("One-sided alpha for conditional power and predictive 
                power approaches is set to ", input$alpha_IA_Norm, ". ")
        d5 = paste0("At interim, responses (Mean (SD)) are , ", input$x1bar_IA_Norm, " (", input$sd1_IA_Norm,
                    ") in group 1, and ", input$x2bar_IA_Norm, " (", input$sd2_IA_Norm,") in group 2.")
       
        paste0(d1, d2, d3, d4, d5)
    })
    
    output$description_IA_Norm = renderText({description_IA_Norm()})
    
    calculation_IA_Norm = eventReactive(input$go_IA_Norm, {
        show_modal_spinner()
        N = as.numeric(input$sampleSize_IA_Norm)
        delta = as.numeric(input$delta_IA_Norm)
        neta = as.numeric(input$neta_IA_Norm)
        alpha = as.numeric(input$alpha_IA_Norm)
        es = as.numeric(input$es_IA_Norm)
        
        
        m1star = as.numeric(stringr::str_split(input$mu1_IA_Norm, ";")[[1]])
        m2star = as.numeric(stringr::str_split(input$mu2_IA_Norm, ";")[[1]])
        sstar = as.numeric(stringr::str_split(input$SD_IA_Norm, ";")[[1]])

        Prior_IA_Norm = cbind(m1star, m2star, sstar)
        Prior_IA_Norm = apply(round(Prior_IA_Norm, digits = 3), 1, Norm_string)

        n1= as.numeric(input$n1_IA_Norm)
        n2= as.numeric(input$n2_IA_Norm)
        x1bar= as.numeric(input$x1bar_IA_Norm)
        x2bar= as.numeric(input$x2bar_IA_Norm)
        s1= as.numeric(input$sd1_IA_Norm)
        s2= as.numeric(input$sd2_IA_Norm)
        
        methodpb = as.numeric(input$methodpb_IA_Norm)
        nsim_p = as.numeric(input$nsim_IA_Norm)

        
        
        res_IA_Norm = IA_Normal_cpp(n1, n2, N, nsim_p, methodpb, 1,
                                    delta, neta, es, alpha, 
                                    x1bar, x2bar, s1, s2,
                                    m1star, m2star, sstar)
        
        remove_modal_spinner()
        
        nrow_IA_Normpb = max(3, length(m1star))
        nrow_IA_Normpp = max(2, length(m1star))
        
        res_IA_Norm_pb = cbind(c("Predictive Probability", paste0("delta=", input$delta_IA_Norm),
                                 paste0("neta=", input$neta_IA_Norm), rep("", nrow_IA_Normpb-3)),
                               c(Prior_IA_Norm, rep("", nrow_IA_Normpb - length(m1star))), 
                               c(round(res_IA_Norm$PredProb, digits = 3),
                                 rep("", nrow_IA_Normpb - length(m1star))))
        
        res_IA_Norm_pp = cbind(c("Predictive Power", paste0("alpha=", input$alpha_IA_Norm),
                                 rep("", nrow_IA_Normpp-2)),
                               c(Prior_IA_Norm, rep("", nrow_IA_Normpp - length(m1star))), 
                               c(round(res_IA_Norm$PredPower, digits = 3),
                               rep("", nrow_IA_Normpp - length(m1star))))
        
        res_IA_Norm_cp = cbind(c("Conditional Power", paste0("alpha=", input$alpha_IA_Norm), 
                                 paste0("effect size under alternative: ", input$es_IA_Norm)), 
                               rep("", 3), rbind(round(res_IA_Norm$CondPower, digits = 3),
                                                 matrix(rep("", 2*length(res_IA_Norm$CondPower)), nrow = 2)))
        
        res_IA_Norm = rbind(res_IA_Norm_pb, res_IA_Norm_pp, res_IA_Norm_cp)
        
        colnames(res_IA_Norm) = c("Method", "Priors: Mean1, Mean2, SD", 
                                 "Probability of a positive trial outcome given interim")
        
        parts_b_IA_Norm = fp_border(width = 3)
        
        as.data.frame(res_IA_Norm) %>%
            flextable() %>%
            theme_booktabs() %>%
            hline(i = c(nrow_IA_Normpb, nrow_IA_Normpb+nrow_IA_Normpp), border = parts_b_IA_Norm) %>%
            autofit() %>%
            fix_border_issues() %>%
            htmltools_value()
        
    })
    
    output$res_IA_Norm = renderUI({
        calculation_IA_Norm()
    })
    
    ###OC - Binomial - One-Sample
    description_OC_Bin_One = eventReactive(input$go_OC_Bin_One, {
      
      n_IA_OC_Bin_One = stringr::str_c(stringr::str_split(input$n_OC_Bin_One, ";")[[1]], collapse = ", ")
      
      d1 = paste0("One-sample trial of binary endpoint, and ", 
                  input$sampleSize_OC_Bin, " subjects planned per group, and ",
                  n_IA_OC_Bin_One, " at each interim. ")
      
      d2 = paste0("Operating characteristics are evaluated with ", input$nsim_OC_Bin_One, 
                  " simulations, with priors from Beta distributions B(",
                  input$as_OC_Bin_One, ",", input$bs_OC_Bin_One, "). ")
      
      d3 = paste0("Hypothesized proportions under null and alternative hypotheses are ", 
                  input$p0_OC_Bin_One, " and ", input$p1_OC_Bin_One,  ", respectively. ")
      
      d4 = paste0("Clinical meaningful treatment difference is set to ", input$delta_OC_Bin_One, ". ")
      
      d5 = paste0("Threshold of posterior probability for declaring efficacy is set to ", input$neta_OC_Bin_One, ". ")
      
      d6 = paste0("One-sided alpha for conditional power and predictive 
                power approaches is set to ", input$alpha_OC_Bin_One, ". ")
      
      d7 = paste0("Threshold for declaring futility is set to ", input$tau_OC_Bin_One, ". ")
      
      paste0(d1, d2, d3, d4, d5, d6, d7)
    })
    
    output$description_OC_Bin_One = renderText({description_OC_Bin_One()})
    
    calculation_OC_Bin_One = eventReactive(input$go_OC_Bin_One, {
      show_modal_spinner()
      N = as.numeric(input$sampleSize_OC_Bin_One)
      n= as.numeric(stringr::str_split(input$n_OC_Bin_One, ";")[[1]])
      n = c(n, N)
      delta = as.numeric(input$delta_OC_Bin_One)
      neta = as.numeric(input$neta_OC_Bin_One)
      alpha = as.numeric(input$alpha_OC_Bin_One)
      p0 = as.numeric(input$p0_OC_Bin_One)
      p1 = as.numeric(input$p1_OC_Bin_One)
      tau = as.numeric(input$tau_OC_Bin_One)
      
      as = as.numeric(input$as_OC_Bin_One)
      bs = as.numeric(input$bs_OC_Bin_One)
      
      a = as.numeric(stringr::str_split(input$a_OC_Bin_One, ";")[[1]])
      b = as.numeric(stringr::str_split(input$b_OC_Bin_One, ";")[[1]])

      Prior = apply(round(cbind(a, b), digits = 3), 1, Beta_string)

      nsim = as.numeric(input$nsim_OC_Bin_One)
      res_OC_Bin_One = OC_Bin_cpp_One(nsim, as, bs, p1, p0, n, delta, neta, alpha, tau, a, b)
      
      remove_modal_spinner()
      
      nrow_OC_Binpb_One = max(3, length(a))
      nrow_OC_Binpp_One = max(2, length(a))
      
      res_OC_Bin_pb_One = cbind(c("Predictive Probability", paste0("delta=", input$delta_OC_Bin_One),
                              paste0("neta=", input$neta_OC_Bin_One), rep("", nrow_OC_Binpb_One-3)),
                            c(Prior, rep("", nrow_OC_Binpb_One-length(a))),
                            rbind(round(res_OC_Bin_One$PredProb, digits = 3),
                                  matrix(rep("", (length(n)+2)*(nrow_OC_Binpb_One-length(a))), ncol = length(n)+2)))
      
      res_OC_Bin_pp_One = cbind(c("Predictive Power", paste0("alpha=", input$alpha_OC_Bin_One),
                              rep("", nrow_OC_Binpp_One-2)),
                            c(Prior, rep("", nrow_OC_Binpp_One-length(a))), 
                            rbind(round(res_OC_Bin_One$PredPower, digits = 3),
                                  matrix(rep("", (length(n)+2)*(nrow_OC_Binpp_One-length(a))), ncol = length(n)+2)))
      
      res_OC_Bin_cp_One = cbind(c("Conditional Power", paste0("alpha=", input$alpha_OC_Bin_One)), 
                            rep("", 2), rbind(round(res_OC_Bin_One$CondPower, digits = 3),
                                                          matrix(rep("", 1*length(res_OC_Bin_One$CondPower)), nrow = 1)))
      
      res_OC_Bin_One = rbind(res_OC_Bin_pb_One, res_OC_Bin_pp_One, res_OC_Bin_cp_One)
      
      res_OC_interims_Bin_One = sapply(1:(length(n)-1), Norm_IA_string)
      
      colnames(res_OC_Bin_One) = c("Method", "Prior", res_OC_interims_Bin_One, 
                               "Not rejecting H0", "Rejecting H0", "Average sample size per group")
      
      parts_b_OC_Bin_One = fp_border(width = 3)
      
      as.data.frame(res_OC_Bin_One) %>%
        flextable() %>%
        add_header_row(values = c("",
                                  "",
                                  rep("Probability of early stopping at interim", length(n)-1),
                                  rep("Cumulative probability at final", 2), ""))%>%
        merge_at(i=1, j=3:(length(n)+1), part = "header") %>%
        merge_at(i=1, j= c(length(n)+2, length(n)+3), part = "header") %>%
        theme_booktabs() %>%
        hline(i = c(nrow_OC_Binpb_One, nrow_OC_Binpb_One+nrow_OC_Binpp_One), border = parts_b_OC_Bin_One) %>%
        autofit() %>%
        fix_border_issues() %>%
        htmltools_value()
    })
    
    output$res_OC_Bin_One = renderUI({
      calculation_OC_Bin_One()
    })
    
    ###OC - Binomial - Two-Sample
    description_OC_Bin = eventReactive(input$go_OC_Bin, {
        
        n_IA_OC_Bin = stringr::str_c(stringr::str_split(input$n_OC_Bin, ";")[[1]], collapse = ", ")
        
        d1 = paste0("Two-sample trial of binary endpoint, and ", 
                    input$sampleSize_OC_Bin, " subjects planned per group, and ",
                    n_IA_OC_Bin, " planned per arm at each interim. ")
        
        d2 = paste0("Operating characteristics are evaluated with ", input$nsim_OC_Bin, 
                    " simulations, with priors from Beta distributions B(",
                    input$a1s_OC_Bin, ",", input$b1s_OC_Bin, ")", "and B(", 
                    input$a2s_OC_Bin, ",", input$b2s_OC_Bin, ")", 
                    " for group 1 and 2, respectively")
        
        d3 = paste0("Clinical meaningful treatment difference is set to ", input$delta_OC_Bin, ". ")
        
        d4 = paste0("Threshold of posterior probability for declaring efficacy is set to ", input$neta_OC_Bin, ". ")
        
        d5 = paste0("One-sided alpha for conditional power and predictive 
                power approaches is set to ", input$alpha_OC_Bin, ". ")
        
        d6 = paste0("Threshold for declaring futility is set to ", input$tau_OC_Bin, ". ")
        
        paste0(d1, d2, d3, d4, d5, d6)
    })
    
    output$description_OC_Bin = renderText({description_OC_Bin()})
    
    calculation_OC_Bin = eventReactive(input$go_OC_Bin, {
        show_modal_spinner()
        N = as.numeric(input$sampleSize_OC_Bin)
        n= as.numeric(stringr::str_split(input$n_OC_Bin, ";")[[1]])
        n = c(n, N)
        delta = as.numeric(input$delta_OC_Bin)
        neta = as.numeric(input$neta_OC_Bin)
        alpha = as.numeric(input$alpha_OC_Bin)
        es = as.numeric(input$es_OC_Bin)
        tau = as.numeric(input$tau_OC_Bin)
        
        a1s = as.numeric(input$a1s_OC_Bin)
        b1s = as.numeric(input$b1s_OC_Bin)
        a2s = as.numeric(input$a2s_OC_Bin)
        b2s = as.numeric(input$b2s_OC_Bin)
        
        a1 = as.numeric(stringr::str_split(input$a1_OC_Bin, ";")[[1]])
        b1 = as.numeric(stringr::str_split(input$b1_OC_Bin, ";")[[1]])
        a2 = as.numeric(stringr::str_split(input$a2_OC_Bin, ";")[[1]])
        b2 = as.numeric(stringr::str_split(input$a2_OC_Bin, ";")[[1]])
        
        Prior1 = apply(round(cbind(a1, b1), digits = 3), 1, Beta_string)
        Prior2 = apply(round(cbind(a2, b2), digits = 3), 1, Beta_string)
        
        nsim = as.numeric(input$nsim_OC_Bin)
        nsim_p = as.numeric(input$nsim_p_OC_Bin)
        methodpb = as.numeric(input$methodpb_OC_Bin)
        
        res_OC_Bin = OC_Bin_cpp(nsim, a1s, b1s, a2s, b2s, 
                                 n, methodpb, nsim_p, 
                                 delta, neta, es, alpha, tau,
                                 a1, b1, a2, b2)
        
        remove_modal_spinner()
        
        nrow_OC_Bin = max(3, length(a1))

        res_OC_Bin_pb = cbind(c("Predictive Probability", paste0("delta=", input$delta_OC_Bin),
                                paste0("neta=", input$neta_OC_Bin), rep("", nrow_OC_Bin-3)),
                              c(Prior1, rep("", nrow_OC_Bin-length(a1))), 
                              c(Prior2, rep("", nrow_OC_Bin-length(a1))), 
                              rbind(round(res_OC_Bin$PredProb, digits = 3),
                              matrix(rep("", (length(n)+2)*(nrow_OC_Bin-length(a1))), ncol = length(n)+2)))
        
        res_OC_Bin_pp = cbind(c("Predictive Power", paste0("alpha=", input$alpha_OC_Bin),
                                rep("", nrow_OC_Bin-2)),
                              c(Prior1, rep("", nrow_OC_Bin-length(a1))), 
                              c(Prior2, rep("", nrow_OC_Bin-length(a1))), 
                              rbind(round(res_OC_Bin$PredPower, digits = 3),
                                    matrix(rep("", (length(n)+2)*(nrow_OC_Bin-length(a1))), ncol = length(n)+2)))
        
        res_OC_Bin_cp = cbind(c("Conditional Power", paste0("alpha=", input$alpha_OC_Bin), 
                                paste0("effect size under alternative: ", input$es_OC_Bin)), 
                              rep("", 3), rep("", 3), rbind(round(res_OC_Bin$CondPower, digits = 3),
                                                            matrix(rep("", 2*length(res_OC_Bin$CondPower)), nrow = 2)))
        
        res_OC_Bin = rbind(res_OC_Bin_pb, res_OC_Bin_pp, res_OC_Bin_cp)
        
        res_OC_interims_Bin = sapply(1:(length(n)-1), Norm_IA_string)
        
        colnames(res_OC_Bin) = c("Method", "Group 1", "Group 2", res_OC_interims_Bin, 
                                 "Not rejecting H0", "Rejecting H0", "Average sample size per group")
        
        parts_b_OC_Bin = fp_border(width = 3)
        
        as.data.frame(res_OC_Bin) %>%
            flextable() %>%
            add_header_row(values = c("",
                                      rep("Priors", 2),
                                      rep("Probability of early stopping at interim", length(n)-1),
                                      rep("Cumulative probability at final", 2), ""))%>%
            merge_at(i=1, j=2:3, part = "header") %>%
            merge_at(i=1, j=4:(length(n)+2), part = "header") %>%
            merge_at(i=1, j= c(length(n)+3, length(n)+4), part = "header") %>%
            theme_booktabs() %>%
            hline(i = (1:2)*nrow_OC_Bin, border = parts_b_OC_Bin) %>%
            autofit() %>%
            fix_border_issues() %>%
            htmltools_value()
    })
    
    output$res_OC_Bin = renderUI({
        calculation_OC_Bin()
    })
    
    ###OC - Normal
    description_OC_Norm = eventReactive(input$go_OC_Norm, {
        n_IA_OC_Norm = stringr::str_c(stringr::str_split(input$n_OC_Norm, ";")[[1]], collapse = ", ")
        
        d1 = paste0("Two-sample trial of normal endpoint, and", 
                    input$sampleSize_OC_Norm, " subjects planned per group, and ",
                    n_IA_OC_Norm, " planned per arm in each interim. ")
        
        d2 = paste0("Operating characteristics are evaluated with ", input$nsim_OC_Norm, 
                    " simulations, with common standard deviation of ",
                    input$sdcom_OC_Norm, "for both groups, and means from N(",
                    input$mu1s_OC_Norm, ",", input$ss_OC_Norm, "^2)", "and N(", 
                    input$mu2s_OC_Norm, ",", input$ss_OC_Norm, "^2)", 
                    " for group 1 and 2, respectively")
        
        d3 = paste0("Clinical meaningful treatment difference is set to ", 
                    input$delta_OC_Norm, ". ")
        
        d4 = paste0("Threshold of posterior probability for declaring efficacy is set to ", 
                    input$neta_OC_Norm, ". ")
        
        d5 = paste0("One-sided alpha for conditional power and predictive 
                power approaches is set to ", input$alpha_OC_Norm, ". ")
        
        d6 = paste0("Threshold for declaring futility is set to ", input$tau_OC_Norm, ". ")

        paste0(d1, d2, d3, d4, d5, d6)
    })
    
    output$description_OC_Norm = renderText({description_OC_Norm()})
    
    calculation_OC_Norm = eventReactive(input$go_OC_Norm, {
        show_modal_spinner()
        N = as.numeric(input$sampleSize_OC_Norm)
        n= as.numeric(stringr::str_split(input$n_OC_Norm, ";")[[1]])
        n = c(n, N)
        delta = as.numeric(input$delta_OC_Norm)
        neta = as.numeric(input$neta_OC_Norm)
        alpha = as.numeric(input$alpha_OC_Norm)
        es = as.numeric(input$es_OC_Norm)
        tau = as.numeric(input$tau_OC_Norm)
        
        mu1star = as.numeric(input$mu1s_OC_Norm)
        mu2star = as.numeric(input$mu2s_OC_Norm)
        sdstar = as.numeric(input$ss_OC_Norm)
        sdcom = as.numeric(input$sdcom_OC_Norm)
        
        m1star = as.numeric(stringr::str_split(input$mu1_OC_Norm, ";")[[1]])
        m2star = as.numeric(stringr::str_split(input$mu2_OC_Norm, ";")[[1]])
        sstar = as.numeric(stringr::str_split(input$SD_OC_Norm, ";")[[1]])
        
        Prior_OC_Norm = cbind(m1star, m2star, sstar)
        Prior_OC_Norm = apply(round(Prior_OC_Norm, digits = 3), 1, Norm_string)
        
        nsim = as.numeric(input$nsim_OC_Norm)
        nsim_p = as.numeric(input$nsim_p_OC_Norm)
        methodpb = as.numeric(input$methodpb_OC_Norm)
        
        res_OC_Norm = OC_Normal_cpp(nsim, mu1star, mu2star, sdstar,
                                    n, sdcom, nsim_p, methodpb, 1,
                                    delta, neta, es, alpha, tau,
                                    m1star, m2star, sstar)
        
        remove_modal_spinner() 
        
        nrow_OC_Norm = max(3, length(m1star))
        
        res_OC_Norm_pb = cbind(c("Predictive Probability", paste0("delta=", input$delta_OC_Norm),
                                        paste0("neta=", input$neta_OC_Norm), rep("", nrow_OC_Norm-3)),
                              c(Prior_OC_Norm, rep("", nrow_OC_Norm-length(m1star))), 
                              rbind(round(res_OC_Norm$PredProb, digits = 3),
                                    matrix(rep("", (length(n)+2)*(nrow_OC_Norm-length(m1star))), ncol = length(n)+2)))
        
        res_OC_Norm_pp = cbind(c("Predictive Power", paste0("alpha=", input$alpha_OC_Norm),
                                 rep("", nrow_OC_Norm-2)),
                               c(Prior_OC_Norm, rep("", nrow_OC_Norm-length(m1star))), 
                               rbind(round(res_OC_Norm$PredPower, digits = 3),
                                     matrix(rep("", (length(n)+2)*(nrow_OC_Norm-length(m1star))), ncol = length(n)+2)))
        
        res_OC_Norm_cp = cbind(c("Conditional Power", paste0("alpha=", input$alpha_OC_Norm), 
                                      paste0("effect size under alternative: ", input$es_OC_Norm)), 
                               rep("", 3), rbind(round(res_OC_Norm$CondPower, digits = 3),
                                         matrix(rep("", 2*length(res_OC_Norm$CondPower)), nrow = 2)))
        
    
        res_OC_interims = sapply(1:(length(n)-1), Norm_IA_string)
        
        res_OC_Norm = rbind(res_OC_Norm_pb, res_OC_Norm_pp,  res_OC_Norm_cp)
        
        colnames(res_OC_Norm) = c("Method", "Priors: Mean1, Mean2, SD", 
                                  res_OC_interims, "Not rejecting H0",
                                  "Rejecting H0", "Average sample size per group")

        parts_b_OR_Norm = fp_border(width = 3)
        
        as.data.frame(res_OC_Norm) %>%
            flextable() %>%
            add_header_row(values = c("Method/Case", "Method/Case",
                                  rep("Probability of early stopping at interim", length(n)-1),
                                  rep("Cumulative probability at final", 2), "")) %>%
            merge_at(i=1, j=1:2, part = "header") %>%
            merge_at(i=1, j=3:(length(n)+1), part = "header") %>%
            merge_at(i=1, j= c(length(n)+2, length(n)+3), part = "header") %>%
            theme_booktabs() %>%
            hline(i = (1:2)*nrow_OC_Norm, border = parts_b_OR_Norm) %>%
            autofit() %>%
            fix_border_issues() %>%
            htmltools_value()
    })
    
    output$res_OC_Norm = renderUI({
        calculation_OC_Norm()
    })
}

# Create Shiny app ----
shinyApp(ui, server)