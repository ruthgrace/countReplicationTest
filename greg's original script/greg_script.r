rdirichlet <- function (n, alpha)
{
  if(length(n) > 1) n <- length(n)
  #if(length(n) == 0 || as.integer(n) == 0) return(numeric(0))
  #n <- as.integer(n)
  if(n < 0) stop("value(n) can not be negative in rtriang")

  if(is.vector(alpha)) alpha <- t(alpha)
  l <- dim(alpha)[2]
  x <- matrix(rgamma(l * n, t(alpha)), ncol = l, byrow=TRUE)  # Gere le recycling
  return(x / rowSums(x))
}

prior <- 0.5 #dirichlet and plotting prior

#original data
data <- read.table("/Users/ggloor/Documents/0_techrep_blog/analysis_exact/td_OTU_tag_mapped_exact.txt", comment.char="",sep="\t", header=T, check.names=F, skip=1, row.names=1)
newdR <- data[order(colnames(data))]

d.21 <- data.frame(newdR[,"CACTAC-CACTAC-1"], newdR[,"CACTAC-CACTAC-2"])

rownames(d.21) <- rownames(newdR)

colnames(d.21) <- c("Rep_A", "Rep_B")

d.21[d.21 == 0] <- prior

round.digits <- 1 # This affects the low end variation estimate strongly
# dirichlet
xa <- round(rdirichlet(100000, d.21[,"Rep_A"]) *  sum(d.21[,"Rep_A"]), round.digits)
xa[ xa==0 ] <- prior
ya <- apply(xa, 1, function(x){ abs(log2(x) - log2(d.21[,"Rep_A"]) ) })
ya.q <- apply(ya, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})
yaE <- apply(xa, 1, function(x){ log2(abs(x- d.21[,"Rep_A"]) ) })
yaE[ is.nan(yaE) ] <- 0
yaE.q <- apply(yaE, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

# rarefaction dataset
#rare1 <- matrix(data=NA, nrow=1000,ncol=nrow(newdR))
#rare2 <- matrix(data=NA, nrow=1000,ncol=nrow(newdR))
#rare_R_1 <- matrix(data=NA, nrow=1000,ncol=nrow(newdR))
#rare_R_2 <- matrix(data=NA, nrow=1000,ncol=nrow(newdR))
#for(i in 0:999){
#       file = paste("/Users/ggloor/Documents/0_techrep_blog/analysis_exact/qiime/jknife_txt/jknife",i,".txt", sep="")
#       d <- read.table(file,comment.char="", header=T, check.names=F, skip=1, row.names=1, sep="\t")
#       j <- i+1
#       rare1[j,] <- d[,"CACTAC-CACTAC-1"]
#       rare2[j,] <- d[,"CACTAC-CACTAC-2"]
#       file = paste("/Users/ggloor/Documents/0_techrep_blog/analysis_exact/qiime/rare_txt/jknife",i,".txt", sep="")
#       d <- read.table(file,comment.char="", header=T, check.names=F, skip=1, row.names=1, sep="\t")
#       j <- i+1
#       rare_R_1[j,] <- d[,"CACTAC-CACTAC-1"]
#       rare_R_2[j,] <- d[,"CACTAC-CACTAC-2"]
#}
#write.table(rare_R_1, file="rare_R_1.txt")
#write.table(rare_R_2, file="rare_R_2.txt")

rare1 <- read.table("~/Documents/0_techrep_blog/analysis_exact/rare_R_1.txt")
rare2 <- read.table("~/Documents/0_techrep_blog/analysis_exact/rare_R_2.txt")

#correction for number of reads in rarefied samples
rare2[rare2 == 0] <- prior
rare1[rare1 == 0] <- prior
ra <- apply(rare1, 1, function(x){ abs(log2(x/(10000/sum(d.21[,"Rep_A"]))) - log2(d.21[,"Rep_A"]) ) })
ra.q <- apply(ra, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

raE <- apply(rare1, 1, function(x){ log2(abs(x/(10000/sum(d.21[,"Rep_A"])) - d.21[,"Rep_A"]) ) })
raE[ is.nan(raE) ] <- 0
raE.q <- apply(raE, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

# resampling
rare2[rare2 == 0] <- prior
raR <- apply(rare2, 1, function(x){ abs(log2(x/(10000/sum(d.21[,"Rep_A"]))) - log2(d.21[,"Rep_A"]) ) })
raR.q <- apply(raR, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

raRE <- apply(rare2, 1, function(x){ log2(abs(x/(10000/sum(d.21[,"Rep_A"])) - d.21[,"Rep_A"]) ) })
raRE[ is.nan(raRE) ] <- 0
raRE.q <- apply(raRE, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})


#A estimated variance plotted vs A-B difference
mincount <- apply(d.21[,1:2], 1, min)
rat <- sum(d.21[,1]) / sum(d.21[,2])

par(mfrow=c(1,2))

plot( log2(mincount), abs(log2(abs(d.21[,"Rep_A"] - round( rat * d.21[,"Rep_B"], 1)))), xaxt="n", yaxt="n", main="Difference", pch=19, cex=1.1, col=rgb(0,0,0,0.3),  ylab="Absolute Difference", xlab="Feature Count", ylim=c(0,8), type="h", lwd=4)
points(loess.smooth(x=log2(d.21[,"Rep_A"]), y= raE.q[3,]), type="l", col="red", lwd=2)
points(loess.smooth(x=log2(d.21[,"Rep_A"]), y= raRE.q[3,]), type="l", col="orange", lwd=2)
points(loess.smooth(x=log2(d.21[,"Rep_A"]), y= yaE.q[3,]), type="l", col="blue", lwd=2)
axis(1, at=log2(c(1,10,100,1000,10000)), labels=c(1,10,100,"1e3","1e4"))
axis(2, at=log2(c(1,10,100)), labels=c(1,10,100))

#rect(log2(1e3), -1, log2(1e5), 10.5,  col=rgb(0,0,0,0.15), lwd=NULL, border = NA)
#rect(2, -1, 5, 10.5,  col=rgb(.64,.2,.06,0.15), lwd=NULL, border = NA)


plot( log2(mincount), abs(log2(d.21[,"Rep_A"]) - log2(round(rat * d.21[,"Rep_B"], 1))), xaxt="n", yaxt="n", main="Ratio ", pch=19, cex=1.1, col=rgb(0,0,0,0.3),  ylab="Replicate Ratio", xlab="Feature Count", ylim=c(0,3), type="h", lwd=4)
points(loess.smooth(x=log2(d.21[,"Rep_A"]), y= ra.q[3,]), type="l", col="red", lwd=2)#points(loess.smooth(x=log2(d.21[,"Rep_A"]), y= raR.q[3,]), type="l", col="orange", lwd=2)
#points(log2(d.21[,"Rep_A"]),raR.q[3,], pch=19,cex=0.5, col="orange")
points(loess.smooth(x=log2(d.21[,"Rep_A"]), y= ya.q[3,]), type="l", col="blue", lwd=2)
#points(log2(d.21[,"Rep_A"]), ya.q[3,], pch=19, cex=0.5, col=rgb(0,0,1,0.4))
axis(1, at=log2(c(1,10,100,1000,10000)), labels=c(1,10,100,"1e3","1e4"))
axis(2, at=log2(c(1,2,4,8)), labels=c(1,2,4,8))

#rect(-2, -1, 2, 5.5,  col=rgb(0,0,0,0.15), lwd=NULL, border = NA)
#rect(2, -1, 5, 5.5,  col=rgb(.64,.2,.06,0.15), lwd=NULL, border = NA)