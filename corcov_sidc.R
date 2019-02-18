
##############################################################
#     multivariate time seres analysis
#     corcov_sidc.R
##############################################################

# Created   dd mmm yyyy, J. Riggs
# Modified 
#			dd mmm yyyy, R. Howe, _ns rhowe to work with sidc_gn data. http://stat.ethz.ch/R-manual/R-patched/library/stats/html/smooth.spline.html
#			08 May 2018, J. Riggs, modeified from wolf_ns_spline_arima_1990-1707.R for multivariate time series analysis


# Working directories
# -------------------
# \Users\howe\desktop

# Imports
# -------
# sidc_corcov.csv

# Exports
# -------
# sidc.RData



##############################################################
# Initialization
##############################################################

library(MASS)
library(stats)
library(graphics)
library(xtable)
library(ggplot2)
library(PerformanceAnalytics)
#library("ggfortify")
#library("reshape2")
#library(rgl)
#require(fUnitRoots)
require(fBasics)
#library(car)
#library(fGarch)
#library(rmgarch)
#library(curl)
#library(devtools)
#library(quantmod)
#library(FinTS)
#library(forecast)
#library(ggpubr)

options(digits=4)

#source('backtest.R')
#source('ccm.R')
#source ('mq.R')


Path <- "/Users/howe/desktop/sidc_gn/"
setwd(Path)
(WD <- getwd())


##########################################################
#     Functions
##########################################################

fetch <- function(fn,ext) {# fn <- Ex
	infile <- paste0(WD, "/", fn, ".", ext)
     X <- data.frame(read.csv(infile, header=TRUE))
	}

WriteCSV <- function(RdataSet, CSVfileName) {
	outfile <- paste(WD, CSVfileName, sep="/")
	write.csv(RdataSet, file=outfile, row.names=F)
	}

##########################################################
#     Import data
##########################################################

#( Ex <- "sidc_corcov")   # Specify the name of the data file
#(ver <- "Obs")    # only 24th cycle JD 2455168
( Ex <- "sidc_gn_1818")   # Specify the name of the data file
(ver <- "2009")

X <- fetch(Ex,"csv")   # import csv data file
summary(X)
nrow(X)
H <- X    # hold current month's raw data
# X <- H  # restore
#(Ex <- "sidc")

# rename data column names
names(H) <- c("date","gn","ISN","gnRSN","logISN","logRg")
summary(H)
nrow(H)
X <- H
summary(X)
nrow(X)

H <- X     # hold modified monthly data losing original raw data
# X <- H    # restore if needed

(outfile <- paste0(WD, "/", Ex, ver, ".RData"))
save(X, file=outfile, ascii=FALSE)


##########################################################
#     Univariate ts EDA
##########################################################

# additional libraries

#(Ex <- "Carr")
#(ver <- "Obs")    # only 24th cycle JD 2455168
#(ver <- "SDO")

(infile <- paste0(WD, "/", Ex, ver, ".RData"))
load(infile)     # loads as X
summary(X)
nrow(X)

H <- X    # hold current month's raw data
# X <- H  # restore

# subset to cycle 24 only

#X <- H[H$jd >= 2455167,]
X$date <- as.character(X$date)
str(X)

######################################################
     part <- "VAR"  # multivariate time series
######################################################

#(Ex <- "Carr")
#(ver <- "Obs")    # only 24th cycle JD 2455168

(infile <- paste0(WD, "/", Ex, ver, ".RData"))
load(infile)     # loads as X
summary(X)
nrow(X)

#H <- X    # hold current month's raw data
# X <- H  # restore

# subset to cycle 24 only

#X <- H[H$jd >= 2455167,]
X$date <- as.character(X$date)
str(X)


#require(MTS)
X$ymd <- as.character(paste(as.character(X$date),as.character(X$jd),sep="-"))
z <- data.frame(diff(X$gnRSN),diff(X$logRg),row.names=X$ymd[-1])
colnames(z) <- c("gnRSN","logRg")
str(z)

#m <- VAR(z,p=20)
#m3.1 <- m
#m <- refVAR(m, thres=1.96)
#str(m)
#MTSdiag(m)


cor.fun <- function(x){
  cor(x)[1,2]
  }

cov.fun <- function(x){
  cov(x)[1,2]
  }

#Differenced data


# pick months for the plots
width <- 13
#width <- 2* 15
#width <- 2* 30
#width <- 2*180
#width <- 2*360

roll.cov <- rollapply(zoo(z), FUN=cov.fun, width=width, by.column=FALSE, align="right")

roll.cor <- rollapply(zoo(z), FUN=cor.fun, width=width, by.column=FALSE, align="right")
 mean(roll.cor)

plot.ts(X[,c("gnRSN","logRg")],col=c(1,2),plot.type="multiple",main="GN Ratios Log ISN")
x11()

par(mfrow=c(2,1))
plot(roll.cov, main="Rolling Covariance",
      ylab="Covariance", lwd=2, col="blue")
grid()
abline(h=cov(z)[1,2], lwd=2, col="red")
plot(roll.cor, main="Rolling Correlations",
      ylab="Correlation", lwd=2, col="blue")
grid()
abline(h=cor(z)[1,2], lwd=2, col="red")
par(mfrow=c(1,1))

#x11()

a <- data.frame(z,ret=z[,1],row.names=row.names(z))
head(a)

# Rolling correlation can be used to examine how relationships between two assets change over time.
   yl <- "RollCorr"
   xl <- ""
  ilk <- width
(main <- paste0("GN Ratio, Log ISN ", ilk, "-Month Correlation ", xl))
 (loc <- paste0("/Users/howe/desktop/sidc_gn/", Ex, ver, part, yl, "By", xl, ilk, ".png"))
png(loc)
	chart.RollingCorrelation(a[, 1, drop=FALSE], a[, 2, drop=FALSE], colorset=rich8equal, legend.loc="bottomright", width=2*as.numeric(ilk), main = main)
	dev.off()

cor(z[,1:2])


######################################################
     part <- "Volatility"  # volatility models
######################################################

#(Ex <- "Carr")
#(ver <- "Obs")    # only 24th cycle JD 2455168

(infile <- paste0(WD, "/", Ex, ver, ".RData"))
load(infile)     # loads as X
summary(X)
nrow(X)

H <- X    # hold current month's raw data
# X <- H  # restore

# subset to cycle 24 only

#X <- H[H$jd >= 2455167,]
X$date <- as.character(X$date)
str(X)

summary(z)
# test the series for ARCH effects
# (mean or volatility models)

# scatterplot to see if there's any obvious correlation

   yl <- "gnRSN"
   xl <- "logRg"
  (is <- coef(lm(gnRSN ~ logRg, data = z)))
(main <- paste(part, yl,"vs", xl))
 (loc <- paste0("/Users/howe/desktop/sidc_gn/", Ex, ver, part, "Scatter", yl, "By", xl, ".png"))
   gp <- ggplot(z, aes(x=gnRSN, y=logRg)) +
         geom_abline(intercept=is[[1]], slope=is[[2]]) +
         geom_point(col='grey45') + 
#         ylim(0,4) +
#         xlim(0.75,2.05) +
         ggtitle(main) + 
         xlab(xl) + 
         ylab(yl) +
         geom_point(col='grey45') + 
         theme(axis.text=element_text(size=16)) +
         theme(axis.title=element_text(size=20)) +
         theme(plot.title=element_text(size=20,hjust = 0.5)) +
         theme(panel.background = element_rect(fill = "grey92"))
         ggsave(loc)


X$ymd <- as.character(paste(as.character(X$jd),as.character(X$Month),sep="-"))
z <- data.frame(diff(X$gnRSN),diff(X$logRg),row.names=X$ymd[-1])
colnames(z) <- c("gnRSN","logRg")
str(z)
