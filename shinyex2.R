library(shiny)
library(igraph)
library(tidyverse)
source("~/Documents/shiny/data.sample.R")


# coeffa, coeffb, coeffc, effect, ntot, nobs
# Define UI ----
ui <- fluidPage(
  titlePanel("Causal Inference for a Time-Varying Treatment"),
  sidebarLayout(
    sidebarPanel(
      sliderInput("coeffa", withMathJax(helpText("$$a$$")),
                  min = 0, max = 1, value = 0.3),
      sliderInput("coeffb", withMathJax(helpText("$$b$$")),
                  min = -2, max = 0, value = -0.75),
      sliderInput("coeffc", withMathJax(helpText("$$c$$")),
                  min = 0, max = 1, value = 0.3),
      sliderInput("coeffd", withMathJax(helpText("$$d$$")),
                  min = 0, max = 1, value = 0.75),
      sliderInput("effect", withMathJax(helpText("$$\\theta$$")),
                  min = -1, max = 0, value = -0.25),
      sliderInput("su", withMathJax(helpText("$$\\sigma_u$$")),
                  min = 0, max = 5, value = 0),
      withMathJax(helpText("$$G$$ 
                           $$\\sigma \\sim $$"))
    ),
    mainPanel(withMathJax(helpText("$$u \\sim N(0, \\sigma^2_u) $$
                                   $$X_1 \\sim N(0, \\sigma^2_x)$$ 
                                   $$g(P(Z_1 \\vert X_1)) = -0.1 + a x_1 $$
                                   $$P(X_2 \\vert X_1, Z_1) \\sim N(x_1 + b Z_1 + u, \\sigma^2_x/16) $$
                                   $$g(P(Z_2 \\vert Z_1, X_2)) = -0.1 + 2 z_1 + c x_2 $$
                                   $$g(P(Y \\vert Z_1, Z_2, Z_3, X_2)) = -0.1 + \\theta \\sum_{k=1}^3 z_k + d x_2 + u $$")),
    plotOutput("hist"),
    plotOutput("g"),
    h4("Special Cases"),
    p("- d=0; The effect, theta, is equal to -0.25. In this case, the effect estimates from traditional models (logistic regression in this case) and causal models can be compared.
      However, when d != 0, some of the effect of treatment is mediated through X2. In this case, the causal effect includes the effect transmitted through X2. Therefore, casual effect estimates can't be directly compared with 
      logistic regression parameters since the latter are conditioned on X2."),
    p("- U=0; no unmeasured confounders. When d=0, logistic regression and causal methods provide unbiased estimates of the true effect (-0.25). When d != 0, logistic regression model provides a conditional causal effect and causal models provide a marginal causal effect."),
    p("- U != 0; this is when causal models should be used! Conditioning on X2 induces a non-causal association between Z1 and Y (through U).")
    
      
    
    )
  )
)

# Define server logic ----
server <- function(input, output) {
  output$g = renderPlot({
    n = 6
    
    e = c(1,2, # X0 - A1
          1,3, # X0 
          2,3,
          2,5,
          3,4, # tvc
          3,5,
          4,5,
          6,3,
          6,5)
    
    gt = make_graph(n, edges = e)
    
    V(gt)$name = c("X1", "Z1", "X2", "Z2", "Y","U")
    V(gt)$size <- 30
    V(gt)$frame.color = c(rep("black", 5), "white")
    E(gt)$weight = 1
    
    
    vertex.color= c( rep(c("dodgerblue","red4"),2), "yellow", "grey")
    set.seed(12)
    plot(gt, vertex.label.color= c("black","white"), 
         vertex.label.dist=0,
         vertex.color = vertex.color,
         edge.label = c("a", "", "b", "", "c", "d", "", "" ),
         vertex.label.cex= 1.5, vertex.shape = c("circle", "square"))
    
  })
  
  output$hist = renderPlot( {

    a = data.fit(nsim = 10, coeffa=input$coeffa, coeffb = input$coeffb, coeffc = input$coeffc,
                 coeffd = input$coeffd,
                 su = input$su, effect = input$effect, ntot = 10000, nobs = 1000)
    
    ggplot(a, aes(x = sim, y = theta.hat, group = mod, color = mod)) + 
      geom_errorbar(aes(ymin = theta.hat - 2 * se, ymax = theta.hat + 2 * se), 
                    position = position_dodge()) + 
      geom_hline(yintercept=mean(a[,5])) + theme_classic() + 
      scale_color_manual(values=c('#999999','#E69F00'))
    # hist(rnorm(10))
    } )

}

# Run the app ----
shinyApp(ui = ui, server = server)