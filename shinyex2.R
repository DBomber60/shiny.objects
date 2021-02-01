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
    mainPanel("here is some text explaining",
              withMathJax(helpText("$$u \\sim N(0, \\sigma^2_u) $$
                                   $$X_1 \\sim N(0, \\sigma^2_x)$$ 
                                   $$g(P(Z_1 \\vert X_1)) = -0.1 + a x_1 $$
                                   $$P(X_2 \\vert X_1, Z_1) \\sim N(x_1 + b Z_1 + u, \\sigma^2_x/16) $$
                                   $$g(P(Z_2 \\vert Z_1, X_2)) = -0.1 + 2 z_1 + c x_2 $$
                                   $$g(P(Y \\vert Z_1, Z_2, Z_3, X_2)) = -0.1 + \\theta \\sum_{k=1}^3 z_k + d x_2 + u $$")),
              plotOutput("hist"),
    plotOutput("g"),
    textOutput("text"))
  )
)

# Define server logic ----
server <- function(input, output) {
  output$g = renderPlot({
    n=7
    gt = graph(n, edges = c(1,2))
    V(gt)$name = c("X0", "A1", "X1", "A2", "X3", "A3", "Y")
    plot(gt)
    V(gt)$size <- 30
    E(gt)$weight = 1
    
    
    vertex.color= c( rep(c("dodgerblue","red4"),3), "yellow")
    
    plot(gt, vertex.label.color= c("black","white"), 
         vertex.label.dist=0,
         vertex.color = vertex.color,
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
  
  output$text <- renderText({ 
    "General Setup"
  })

}

# Run the app ----
shinyApp(ui = ui, server = server)