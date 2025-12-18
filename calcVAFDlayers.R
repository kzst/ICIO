library(readr)
sector_LUT <- read_delim("sector_LUT.csv",
                         delim = ";", escape_double = FALSE, trim_ws = TRUE)
View(sector_LUT)
country_LUT_SP <- read_delim("country_LUT_SP.csv",
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
load("data/ICIO2023econVA.RData")
load("data/ICIO2023econFD.RData")

library(igraph)
library(stringr)
library(leidenAlg)

COLS<-colnames(ICIO2023econVA)
VAcent<-matrix(0,nrow=nrow(ICIO2023econVA),ncol=44)

for (i in 1:nrow(ICIO2023econVA)){
  VA<-matrix(0,ncol=nrow(country_LUT_SP),nrow=nrow(sector_LUT))
  for (j in 1:ncol(ICIO2023econVA)){
    CTRY<-substr(COLS[j],1,3)
    SECT<-paste("D",substr(COLS[j],5,str_length(COLS[j])),sep="")
    ctry<-as.vector(as.matrix(country_LUT_SP[country_LUT_SP$Code %in% CTRY,1]))
    sect<-as.vector(as.matrix(sector_LUT[sector_LUT$`Old code` %in% SECT,1]))
    if (!is.numeric(ctry)) ctry<-77
    if ((sect>0)&&(ctry>0)){
      VA[sect,ctry]<-ICIO2023econVA[i,j]
    }
  }
  VAcent[i,]<-rowSums(VA)[1:44]
}
