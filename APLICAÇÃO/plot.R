###########################################################################latex
if (!require("pacman")) install.packages("pacman")
pacman::p_load('tidyverse','vroom',"openxlsx", "knitr", "kableExtra", "xtable")

dados <- read.xlsx("C:\\Users\\matheus.braga\\Documents\\rascunho.xlsx")


xtable(dados)


dados %>%
  kbl(caption = "Example Table", format = "latex", booktabs = TRUE)


#############################################################modificações plot()
library(gamlss)
library(tidyverse)
library(Cairo)

data(abdom)
abd10<-gamlss(y~pb(x),sigma.fo=~pb(x),data=abdom,family=BCT,
              trace=FALSE)

#plot(abd10)

####alterações
newpar <-par(mfrow = c(1,1), #disposição dos gráficos
            mar = par("mar") + c(0,1,0,0),
            col.axis = "black",
            col = "darkgreen", 
            col.main = "black",
            col.lab = "black",
            pch = 1,
            #cex = 0.45, 
            cex.lab = 0.8, 
            cex.axis = 0.8, 
            cex.main = 0.8, 
            bg = "white")

plot(abd10,par = newpar) 


#
# getwd()
# 
# png(filename = str_c("plot2", ".png"))
# plot(abd10,par = newpar, main = "OI", which = 1)
# dev.off()