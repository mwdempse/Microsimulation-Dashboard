# Libraries
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinydashboard)
library(waiter)
library(ggplot2)
library(reshape2)
library(ggsci)
# library(thematic)
# library(bslib)

# Define UI for application
ui <- fluidPage(
  # theme = bslib::bs_theme(bootswatch = 'darkly'),
  waiter::use_waiter(),
  titlePanel("Health Microsimulation Model"),
  sidebarLayout(
    sidebarPanel(
      position = 'left',
      fluid = T,
      h4("Introduction"),
      p("This tool performs a microsimulation based on the Sick-Sicker Model, where individuals move between the states healthy, sick, sicker, and dead. For the given inputs, the model generates an incremental cost-effectiveness analysis as well as helpful diagrams on the proportion of individuals in different health states and the distribution of costs and quality-adjusted life years (QALYs) across individuals."),
      tabsetPanel(
        tabPanel("General", 
                 numericInput("n_i", "Number of Individuals:", 10000, min = 1000, max = 100000, step = 1000),
                 numericInput("n_t", "Number of Cycles:", 30, min = 10, step = 1),
                 numericInput("d_c", "Discount Rate for Costs:", 0.03, min = 0, step = 0.01),
                 numericInput("d_e", "Discount Rate for QALYs:", 0.03, min = 0, step = 0.01)
        ),
        tabPanel("Transition Probabilities", 
                 numericInput("p_HD", "Probability to Die when Healthy:", 0.005, min = 0, max = 1, step = 0.001),
                 numericInput("p_HS1", "Probability to Become Sick when Healthy:", 0.15, min = 0, max = 1, step = 0.01),
                 numericInput("p_S1H", "Probability to Become Healthy when Sick:", 0.5, min = 0, max = 1, step = 0.01),
                 numericInput("p_S1S2", "Probability to Become Sicker when Sick:", 0.105, min = 0, max = 1, step = 0.01),
                 numericInput("rr_S1", "Rate Ratio of Death when Sick vs Healthy:", 3.0, min = 0, step = 0.1),
                 numericInput("rr_S2", "Rate Ratio of Death when Sicker vs Healthy:", 10.0, min = 0, step = 0.1),
                 numericInput("rp_S1S2", "Mortality Rate Increase with Every Additional Year Being Sick:", 0.2, min = 0, step = 0.01)
        ),
        tabPanel("Costs & Utilities", 
                 numericInput("c_H", "Cost in $ of Remaining Healthy per Cycle:", 2000, min = 0, step = 100),
                 numericInput("c_S1", "Cost in $ of Being Sick per Cycle:", 4000, min = 0, step = 100),
                 numericInput("c_S2", "Cost in $ of Being Sicker per Cycle:", 15000, min = 0, step = 500),
                 numericInput("c_Trt", "Cost in $ of Treatment per Cycle:", 12000, min = 0, step = 500),
                 numericInput("u_H", "Utility when Healthy:", 1, min = 0, max = 1, step = 0.01),
                 numericInput("u_S1", "Utility when Sick:", 0.75, min = 0, max = 1, step = 0.01),
                 numericInput("u_S2", "Utility when Sicker:", 0.5, min = 0, max = 1, step = 0.01),
                 numericInput("u_Trt", "Utility when Treated:", 0.95, min = 0, max = 1, step = 0.01)
        ),
        tabPanel("Advanced Settings", 
                 numericInput("ru_S1S2", "Decrease in Utility per Year Being Sick (ru_S1S2):", 0.03, min = 0, step = 0.001),
                 numericInput("x_lower", "Lower Bound for Individual Effect Modifier (x_lower):", 0.95, min = 0, max = 1, step = 0.01),
                 numericInput("x_upper", "Upper Bound for Individual Effect Modifier (x_upper):", 1.05, min = 1, max = 10, step = 0.01)
        )
      ),
      actionButton("run", "Run Simulation")
    ),
    mainPanel(
      h3("Incremental Cost-Effectiveness Analysis"),
      tableOutput("table_micro"),
      p("QALYs = Quality-Adjusted Life Years"),
      p("ICER = Incremental Costs Effectiveness Ratio"),
      p("MCSE = Monte-Carlo Standard Error"),
      h3("Diagrams"),
      plotOutput("proportions_plot"),
      fluidRow(
        column(
          width = 6,
          selectInput(
            inputId = "n_breaks_thc",
            label = "Bins for Treatment: Healthcare Costs",
            choices = c(10, 20, 30, 40, 50),
            selected = 40
          )
        ),
        column(
          width = 6,
          selectInput(
            inputId = "n_breaks_tqly",
            label = "Bins for Treatment: QALYs",
            choices = c(10, 20, 30, 40, 50),
            selected = 40
          )
        ),
        
        plotOutput("hist_trt", inline = T, width = "1000px", height = "600px"),
        br(),
        column(width = 6,
               selectInput(
                 inputId = "n_breaks_nthc",
                 label = "Bins for No Treatment: Healthcare Costs",
                 choices = c(10, 20, 30, 40, 50),
                 selected = 40
               )
        ),
        column(width = 6,
               selectInput(
                 inputId = "n_breaks_ntqly",
                 label = "Bins for No Treatment: QALYs",
                 choices = c(10, 20, 30, 40, 50),
                 selected = 40
               )
        ),
        plotOutput("hist_no_trt", inline = T, width = "1000px", height = "600px")
      )
    )
  )
)



# Define server logic
server <- function(input, output) {
  # thematic::thematic_shiny()
  observeEvent(input$run, {
    showPageSpinner(
      background = getOption("page.spinner.background", default = "darkgrey"),
      type = getOption("page.spinner.type", default = 1),
      color = getOption("page.spinner.color", default = "#0275D8"),
      caption = 'Calculation in progress. This may take a while...'
    )
    ##################################### Model input #######################################
    n_i   <- input$n_i                # number of simulated individuals
    n_t   <- input$n_t                # time horizon, 30 cycles
    v_n   <- c("H","S1","S2","D")     # the model states: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
    n_s   <- length(v_n)              # the number of health states
    v_M_1 <- rep("H", n_i)            # everyone begins in the healthy state 
    d_c   <- input$d_c                # discounting of costs
    d_e   <- input$d_e                # discounting of QALYs
    v.Trt <- c("No Treatment", "Treatment") # store the strategy names
    
    # Transition probabilities (per cycle)
    p_HD    <- input$p_HD             # probability to die when healthy
    p_HS1   <- input$p_HS1            # probability to become sick when healthy
    p_S1H   <- input$p_S1H            # probability to become healthy when sick
    p_S1S2  <- input$p_S1S2           # probability to become sicker when sick
    rr_S1   <- input$rr_S1            # rate ratio of death when sick vs healthy
    rr_S2   <- input$rr_S2            # rate ratio of death when sicker vs healthy 
    r_HD    <- -log(1 - p_HD)         # rate of death when healthy 
    r_S1D   <- rr_S1 * r_HD           # rate of death when sick
    r_S2D   <- rr_S2 * r_HD           # rate of death when sicker
    p_S1D   <- 1 - exp(- r_S1D)       # probability to die when sick
    p_S2D   <- 1 - exp(- r_S2D)       # probability to die when sicker
    rp_S1S2 <- input$rp_S1S2          # increase of the mortality rate with every additional year being sick
    
    # Cost and utility inputs 
    c_H     <- input$c_H              # cost of remaining one cycle healthy
    c_S1    <- input$c_S1             # cost of remaining one cycle sick
    c_S2    <- input$c_S2             # cost of remaining one cycle sicker
    c_Trt   <- input$c_Trt            # cost of treatment (per cycle)
    
    u_H     <- input$u_H              # utility when healthy 
    u_S1    <- input$u_S1             # utility when sick 
    u_S2    <- input$u_S2             # utility when sicker 
    u_Trt   <- input$u_Trt            # utility when sick(er) and being treated
    ru_S1S2 <- input$ru_S1S2          # decrease in utility of treated sick individuals with every additional year being sick/sicker
    x_lower <- input$x_lower          # lower bound for the individuals' effect modifier at baseline
    x_upper <- input$x_upper          # upper bound for the individuals' effect modifier at baseline
    v_x     <- runif(n_i, x_lower, x_upper) # vector capturing individuals' effect modifier at baseline
    
    ##################################### Functions ###########################################
    
    MicroSim <- function(v_M_1, n_i, n_t, v_n, X = NULL, d_c, d_e, c_H, c_S1, c_S2, c_Trt,
                         u_H, u_S1, u_S2, u_Trt, ru_S1S2, p_HD, p_HS1, p_S1H, p_S1S2, 
                         rr_S1, rr_S2, rp_S1S2, x_lower, x_upper, TS.out = TRUE, TR.out = TRUE, Trt = FALSE, seed = 1) {
      
      v.dwc <- 1 / ((1 + d_c) ^ (0:n_t))   # calculate the cost discount weight based on the discount rate d_c
      v.dwe <- 1 / ((1 + d_e) ^ (0:n_t))   # calculate the QALY discount weight based on the discount rate d_e
      
      # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
      m_M <- m_C <- m_E <- matrix(nrow = n_i, ncol = n_t + 1,
                                  dimnames = list(paste("ind",   1:n_i, sep =" "),
                                                  paste("cycle", 0:n_t, sep =" "))) 
      
      m_M[, 1] <- v_M_1             # indicate the initial health state 
      
      for (i in 1:n_i) {
        set.seed(seed + i)          # set the seed for every individual for the random number generator
        
        # create the dur variable that stores the number of consecutive cycles the individual occupies either when sick or sicker
        dur <- 0                            # the individual start without history        
        m_C[i, 1] <- Costs(m_M[i, 1], Trt)  # estimate costs per individual for the initial health state conditional on treatment
        m_E[i, 1] <- Effs(m_M[i, 1], dur, Trt, X = X[i])  # estimate QALYs per individual for the initial health state conditional on treatment, duration of being sick/sicker and individual characteristics
        
        for (t in 1:n_t) {
          v.p <- Probs(m_M[i, t], dur)         # calculate the transition probabilities at cycle t conditional on the duration of being sick/sicker
          
          m_M[i, t + 1] <- sample(v_n, prob = v.p, size = 1)  # sample the new health state and store that state in matrix m_M 
          m_C[i, t + 1] <- Costs(m_M[i, t + 1], Trt)     # estimate the cost per individual during cycle t + 1 conditional on treatment
          m_E[i, t + 1] <-  Effs(m_M[i, t + 1], dur, Trt, X = X[i])    # estimate the utility per individual during cycle t + 1 conditional on treatment, duration of being sick/sicker and individual characteristics
          
          if (m_M[i, t + 1] == "S1" | m_M[i, t + 1] == "S2") {  # expression to identify sick/sicker individuals
            dur <- dur + 1   # update the duration of being sick/sicker
          } else {
            dur <- 0}        # reset duration variable 
          
        } # close the loop for the time points 
      } # close the loop for the individuals
      
      tc <- m_C %*% v.dwc       # total (discounted) cost per individual
      te <- m_E %*% v.dwe       # total (discounted) QALYs per individual 
      
      tc_hat <- mean(tc)        # average (discounted) cost 
      te_hat <- mean(te)        # average (discounted) QALYs
      
      if (TS.out == TRUE) {  # create a  matrix of transitions across states
        TS <- paste(m_M, cbind(m_M[, -1], NA), sep = "->")  # transitions from one state to the other ###
        TS <- matrix(TS, nrow = n_i)
        rownames(TS) <- paste("Cycle", 0:n_t, sep = " ")    # name the rows of the matrix
        colnames(TS) <- paste("Ind",   1:n_s, sep = " ")    # name the columns of the matrix
      } else {
        TS <- NULL
      }
      
      if (TR.out == TRUE) {  # create a trace from the individual trajectories
        TR <- t(apply(m_M, 2, function(x) table(factor(x, levels = v_n, ordered = TRUE))))
        TR <- TR / n_i                                    # create a distribution trace
        rownames(TR) <- paste("Cycle", 0:n_t, sep = " ")  # name the rows of the matrix
        colnames(TR) <- v_n                               # name the columns of the matrix
      } else {
        TR <- NULL
      }
      
      results <- list(m_M = m_M, m_C = m_C, m_E = m_E, tc = tc, te = te, tc_hat = tc_hat, te_hat = te_hat, TS = TS, TR = TR)   # store the results from the simulation in a list  
      return(results)   # return the results
    } # end of the MicroSim function
    
    #### Probability function
    Probs <- function(M_it, dur) { 
      # M_it:   health state occupied by individual i at cycle t (character variable)
      # dur:    the duration of being sick (sick/sicker)
      
      v_p_it <- rep(NA, n_s)     # create vector of state transition probabilities
      names(v_p_it) <- v_n       # name the vector
      
      # update probabilities of death after first converting them to rates and applying the rate ratio
      r_S1D <-  - log(1 - p_S1D)
      r_S2D <-  - log(1 - p_S2D)
      p_S1D <- 1 - exp(- r_S1D * (1 + dur * rp_S1S2)) # calculate p_S1D conditional on duration of being sick/sicker
      p_S2D <- 1 - exp(- r_S2D * (1 + dur * rp_S1S2)) # calculate p_S2D conditional on duration of being sick/sicker
      
      # update v_p_it with the appropriate probabilities   
      v_p_it[M_it ==  "H"] <- c(1 - p_HS1 - p_HD, p_HS1, 0, p_HD)                    # transition probabilities when healthy
      v_p_it[M_it == "S1"] <- c(p_S1H, 1- p_S1H - p_S1S2 - p_S1D, p_S1S2, p_S1D)     # transition probabilities when sick
      v_p_it[M_it == "S2"] <- c(0, 0, 1 - p_S2D, p_S2D)                              # transition probabilities when sicker
      v_p_it[M_it ==  "D"] <- c(0, 0, 0, 1)                                          # transition probabilities when dead
      ifelse(sum(v_p_it) == 1, return(v_p_it), print("Probabilities do not sum to 1")) # return the transition probabilities or produce an error
    }    
    
    ### Costs function
    Costs <- function (M_it, Trt = FALSE) {  
      # M_it: health state occupied by individual i at cycle t (character variable)
      # Trt:  is the individual being treated? (default is FALSE)
      
      c_it <- 0                                   # by default the cost for everyone is zero
      c_it[M_it == "H"]  <- c_H                   # update the cost if healthy
      c_it[M_it == "S1"] <- c_S1 + c_Trt * Trt    # update the cost if sick conditional on treatment
      c_it[M_it == "S2"] <- c_S2 + c_Trt * Trt    # update the cost if sicker conditional on treatment
      return(c_it)                                # return the costs
    }
    
    ### Health outcome function 
    Effs <- function (M_it, dur, Trt = FALSE, cl = 1, X = NULL) { 
      # M_it: health state occupied by individual i at cycle t (character variable)
      # dur:  the duration of being sick/sicker
      # Trt:  is the individual being treated? (default is FALSE)
      # cl:   the cycle length (default = 1 )
      # X:    the vector or matrix of individual characteristics (optional)
      
      u_it               <- 0        # by default the utility for everyone is zero
      u_it[M_it == "H"]  <- u_H      # update the utility if healthy 
      u_it[M_it == "S1"] <- X * Trt * (u_Trt - dur * ru_S1S2) + (1 - Trt) * u_S1 # update the utility if sick conditional on treatment and duration of being sick/sicker
      u_it[M_it == "S2"] <- u_S2     # update the utility if sicker
      QALYs <- u_it * cl             # calculate the QALYs during cycle t
      return(QALYs)                  # return the results
    }
    
    ##################################### Run the simulation ##################################
    sim_no_trt <- MicroSim(v_M_1, n_i, n_t, v_n, X = v_x, d_c, d_e, TS.out = FALSE, TR.out = TRUE, Trt = FALSE) # run for no treatment
    sim_trt    <- MicroSim(v_M_1, n_i, n_t, v_n, X = v_x, d_c, d_e, TS.out = FALSE, TR.out = TRUE, Trt = TRUE) # run for treatment
    
    # Generate proportion plot
    df_tr <- as.data.frame(sim_no_trt$TR)
    df_tr$Cycle <- 0:input$n_t
    df_tr_melt <- melt(df_tr, id.vars = "Cycle")
    
    # Update health state names in the legend
    df_tr_melt$variable <- factor(df_tr_melt$variable, levels = c("H", "S1", "S2", "D"),
                                  labels = c("Healthy", "Sick", "Sicker", "Dead"))
    
    output$proportions_plot <- renderPlot({
      ggplot(df_tr_melt, aes(x = Cycle, y = value, fill = variable)) +
        geom_area(alpha = 0.8) +
        # scale_fill_brewer(palette = "Dark2") +
        # scale_fill_viridis_d()+
        scale_fill_jco() +
        labs(title = "Proportion of Individuals in Each Health State Over Time",
             x = "Cycle",
             y = "Proportion of Individuals",
             fill = "Health State") +
        scale_y_continuous(label = scales::percent_format(),
                           limits = c(0,1),
                           expand = c(0,0)
        )+
        scale_x_continuous(expand = c(0,0))+
        theme_minimal() +
        theme(
          panel.grid.major = element_blank(),
          axis.line = element_line(),
          axis.text = element_text(face = 'bold',
                                   size = 12)
          # axis.text.x = element_text()
        )
    }#, width = "1000px", height = "600px"
    )
    
    # Generate histograms
    output$hist_trt <- renderPlot({
      par(mfrow = c(1, 2))
      hist(
        sim_trt$tc,
        breaks = as.numeric(input$n_breaks_thc),
        col = "skyblue",
        main = "Treatment: Healthcare Costs",
        xlab = "Total Costs in $",
        ylab = "Frequency",
        border = "black",
        xaxt = "n",
        labels = FALSE,
        xlim = c(0, max(sim_trt$tc) + 100000)
      )
      axis(1,
           at = axTicks(1),
           labels = format(axTicks(1), scientific = FALSE))  # Remove scientific notation
      hist(
        sim_trt$te,
        breaks = as.numeric(input$n_breaks_tqly),
        col = "gold",
        main = "Treatment: QALYs",
        xlab = "Total QALYs",
        ylab = "Frequency",
        border = "black",
        labels = FALSE,
        xlim = c(0, max(sim_trt$te) + 5)
      )
    }, width = 1000, height = 600
    )
    
    output$hist_no_trt <- renderPlot({
      par(mfrow = c(1, 2))
      hist(
        sim_no_trt$tc,
        breaks = as.numeric(input$n_breaks_nthc),
        col = "gray",
        main = "No Treatment: Healthcare Costs",
        xlab = "Total Costs in $",
        ylab = "Frequency",
        border = "black",
        xaxt = "n",
        labels = FALSE,
        xlim = c(0, max(sim_no_trt$tc) + 100000)
      )
      axis(1,
           at = axTicks(1),
           labels = format(axTicks(1), scientific = FALSE))  # Remove scientific notation
      hist(
        sim_no_trt$te,
        breaks = as.numeric(input$n_breaks_ntqly),
        col = "tomato",
        main = "No Treatment: QALYs",
        xlab = "Total QALYs",
        ylab = "Frequency",
        border = "black",
        labels = FALSE,
        xlim = c(0, max(sim_no_trt$te) + 5)
      )
    }, width = 1000, height = 600
    )
    
    ################################# Cost-effectiveness analysis #############################
    
    # store the mean costs (and the MCSE)of each strategy in a new variable C (vector costs)
    v_C <- c(sim_no_trt$tc_hat, sim_trt$tc_hat) 
    se_C<- c(sd(sim_no_trt$tc), sd(sim_trt$tc)) / sqrt(n_i)
    # store the mean QALYs (and the MCSE) of each strategy in a new variable E (vector health outcomes)
    v_E <- c(sim_no_trt$te_hat, sim_trt$te_hat)
    se_E<- c(sd(sim_no_trt$te), sd(sim_trt$te)) / sqrt(n_i)
    
    delta_C <- v_C[2] - v_C[1]                   # calculate incremental costs
    delta_E <- v_E[2] - v_E[1]                   # calculate incremental QALYs
    se_delta_E <- sd(sim_trt$te - sim_no_trt$te) / sqrt(n_i) # Monte Carlo squared error (MCSE) of incremental QALYS
    se_delta_C <- sd(sim_trt$tc - sim_no_trt$tc) / sqrt(n_i) # Monte Carlo squared error (MCSE) of incremental costs
    ICER <- delta_C / delta_E                    # calculate the ICER
    results <- c(delta_C, delta_E, ICER)         # store the values in a new variable
    
    # Create full incremental cost-effectiveness analysis table
    table_micro <- data.frame(
      c(round(v_C, 0),  ""),           # costs per arm
      c(round(se_C, 0), ""),           # MCSE for costs
      c(round(v_E, 3),  ""),           # health outcomes per arm
      c(round(se_E, 3), ""),           # MCSE for health outcomes
      c("", round(delta_C, 0),   ""),  # incremental costs
      c("", round(se_delta_C, 0),""),  # MCSE for incremental costs
      c("", round(delta_E, 3),   ""),  # incremental QALYs 
      c("", round(se_delta_E, 3),""),  # MCSE for health outcomes (QALYs) gained
      c("", round(ICER, 0),      "")   # ICER
    )
    
    rownames(table_micro) = c("No Treatment", "Treatment", "")  # name the rows
    colnames(table_micro) = c("Costs", "Cost MCSE",  "QALYs", "MCSE: QUALYs", "Incremental Costs", "MCSE: Incremental Costs", "QALYs Gained", "MCSE: QALYs Gained", "ICER") # name the columns 
    
    output$table_micro <- renderTable({ table_micro }, rownames = TRUE, colnames = TRUE)
    
    hidePageSpinner()
  })
}

# Run the application
shinyApp(ui = ui, server = server)
