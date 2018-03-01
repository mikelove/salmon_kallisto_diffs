library(shiny)
library(ggplot2)
library(MASS)

tests <- read.csv(file.path("geuvadis_results","tests.csv"))
tests$num <- rownames(tests)
boxdata <- read.csv(file.path("geuvadis_results","boxdata.csv"))
boxdata$center <- factor(boxdata$center, levels=boxdata$center[1:7])
levels(boxdata$center)[1] <- "CNAG"
boxdata$bias <- ifelse(as.integer(boxdata$center) < 4, "high", "low")
boxdata$bias <- factor(boxdata$bias, levels=c("high","low","none"))

# make a density on the scatter
# http://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
tests$density <- get_density(tests$SalmonF, tests$kallistoF, n=400)

tab <- with(tests, 
            table(factor(ifelse(SalmonF > 2*kallistoF, "Salmon", 
                         ifelse(kallistoF > 2*SalmonF, "kallisto", "both")),
                         levels=c("kallisto","both","Salmon"))))
tab.dat <- data.frame(x=c(500,1000,1000),y=c(1500,1500,100),num=as.numeric(tab))

### shiny app:

ui <- fluidPage(mainPanel(fluidRow(
  column(6, plotOutput("plot1", click="plot_click", height=900)),
  column(6, plotOutput("plot2", height=900)
  )),width=12))
th <- theme(axis.title=element_text(size=20), 
            strip.text=element_text(size=16),
            legend.position="none")
server <- function(input, output) {
  output$plot1 <- renderPlot({
    ggplot(tests, aes(x=SalmonF, y=kallistoF, col=density)) + 
      geom_point() +
      geom_abline(intercept=0,slope=c(.5,2)) + 
      geom_text(data=tab.dat, aes(x,y,label=num),size=12,col="red") + 
      xlab("Salmon F test on center") + ylab("kallisto F test on center") + th
  })
  output$plot2 <- renderPlot({
    near <- nearPoints(tests,input$plot_click,maxpoints=1,"SalmonF","kallistoF",threshold=10)
    if (nrow(near) == 0) return()
    boxsub <-subset(boxdata, gene==near$gene)
    boxsub$bias[boxsub$txp != near$txp]  <- "none"
    ggplot(boxsub, aes(x=center,color=bias,fill=bias)) +
      geom_boxplot(aes(ymin=X0,lower=X25,middle=X50,upper=X75,ymax=X100),size=1,stat="identity") +
      ylab("normalized counts") + 
      expand_limits(y=0) + 
      facet_grid(txp~method, scale="free_y") +
      scale_color_brewer(palette="Dark2") + 
      scale_fill_manual(values=c("darkseagreen1","wheat","white")) + th
  })
}
shinyApp(ui = ui, server = server)
