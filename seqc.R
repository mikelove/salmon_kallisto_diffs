library(shiny)
library(ggplot2)
library(MASS)

## replace with "B" for other SEQC samples

dataset <- "A"
#dataset <- "B"

tests <- read.csv(file.path("seqc_results",paste0("seqc_",dataset,"_tests.csv")))

# NA p-values go to 1 (log10 of 0)
# Inf p-values go to the max finite value plus jitter for plotting
tests$SalmonMinusLog10P[is.na(tests$SalmonMinusLog10P)] <- 0
tests$kallistoMinusLog10P[is.na(tests$kallistoMinusLog10P)] <- 0
fin <- is.finite(tests$SalmonMinusLog10P)
tests$SalmonMinusLog10P[!fin] <- max(tests$SalmonMinusLog10P[fin]) + 
  rnorm(sum(!fin),0,2.2)
fin <- is.finite(tests$kallistoMinusLog10P)
tests$kallistoMinusLog10P[!fin] <- max(tests$kallistoMinusLog10P[fin]) + 
  rnorm(sum(!fin),0,2.5)

# make a density on the scatter
# http://slowkow.com/notes/ggplot2-color-by-density/
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
tests$density <- get_density(tests$SalmonMinusLog10P, tests$kallistoMinusLog10P, n=400)

tab <- with(tests, 
            table(factor(ifelse(SalmonMinusLog10P > 2*kallistoMinusLog10P, "Salmon", 
                         ifelse(kallistoMinusLog10P > 2*SalmonMinusLog10P, "kallisto", "both")),
                         levels=c("kallisto","both","Salmon"))))
tab.dat <- data.frame(x=c(100,250,250),y=c(275,250,50),num=as.numeric(tab))

ncounts <- list()
for (m in c("Salmon","kallisto")) {
  ncounts[[m]] <- as.matrix(read.csv(
    file.path("seqc_results",paste0("seqc_",dataset,"_",m,"_ncts.csv")), 
    row.names=1))
}

### shiny app:

ui <- fluidPage(mainPanel(fluidRow(
  column(6, plotOutput("plot1", click="plot_click", height=600)),
  column(6, plotOutput("plot2", height=600)
  )),width=12))
th <- theme(axis.title=element_text(size=20), 
            strip.text=element_text(size=16),
            legend.position="none")
server <- function(input, output) {
  output$plot1 <- renderPlot({
    ggplot(tests, aes(x=SalmonMinusLog10P, y=kallistoMinusLog10P, col=density)) + 
      geom_point() +
      geom_abline(intercept=0,slope=c(.5,2)) + 
      geom_text(data=tab.dat, aes(x,y,label=num),size=8,col="red") + 
      xlab("Salmon -log10(p) on center") + ylab("kallisto -log10(p) on center") + th
    })
  output$plot2 <- renderPlot({
    near <- nearPoints(tests,input$plot_click,maxpoints=1,
                       "SalmonMinusLog10P","kallistoMinusLog10P",threshold=10)
    if (nrow(near) == 0) return()
    this.txp <- as.character(near$txp)
    txps <- as.character(tests$txp[tests$gene == near$gene])
    ntxps <- length(txps)
    # this is a vector 2 * ntxps * 12 
    counts <- c(as.vector(t(ncounts[["Salmon"]][txps,])),
                as.vector(t(ncounts[["kallisto"]][txps,])))
    center.lvls <- c("CNL","BGI","MAY")
    center <- factor(rep(rep(c("BGI","CNL","MAY"),each=4), 2 * ntxps), center.lvls)
    method <- factor(rep(c("Salmon","kallisto"),each=12 * ntxps))
    txp <- factor(rep(rep(txps, each=12), 2))
    bias <- factor(ifelse(center == "CNL", "high", "low"), levels=c("high","low","other"))
    bias[txp != this.txp] <- "other"
    dat <- data.frame(counts, center, method, txp, bias)
    ggplot(dat, aes(x=center, y=counts, col=bias, fill=bias)) +
      geom_jitter(size=3, width=0.1, height=0, shape=21, stroke=2) +
      expand_limits(y=0) + 
      ylab("normalized counts") + 
      scale_color_brewer(palette="Dark2") + 
      scale_fill_manual(values=c("darkseagreen1","wheat","white")) + 
      facet_grid(txp~method, scale="free_y") + th
  })
}
shinyApp(ui = ui, server = server)
