options(error=recover)
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

rarefy <- function(features, counts, samples, replacement) {
  # pool has the name of each OTU duplicated as many times as the count of that OTU
  counts.no.prior <- counts
  counts.no.prior[counts.no.prior==prior] <- 0
  pool <- c(1:sum(counts.no.prior))
  index <- 0
  for (i in c(1:length(features))) {
    if(counts.no.prior[i] >= 1) {
      pool[(index+1):(index+counts.no.prior[i])] <- rep(features[i],counts.no.prior[i])
      index <- index + counts.no.prior[i]
    }
  }
  sample.pool <- sample(pool, samples, replace=replacement)
  newcounts <- counts.no.prior
  newcounts[] <- 0
  aggregate.pool <- table(sample.pool)
  newcounts[match(dimnames(aggregate.pool)[[1]],features)] <- aggregate.pool
  return(newcounts)
}

prior <- 0.5 #dirichlet and plotting prior
round.digits <- 1 # dirichlet values rounded to 1 digit after decimal point


data <- read.table("data/td_OTU_tag_mapped_lineage.txt", comment.char="",sep="\t", header=T, check.names=F, skip=1, row.names=1)

# three replicates: samples in the a and b replicate end with 'a' and 'b'
# samples in the C replicate end with 'C' and often have more detailed names
a <- data[,grep("a$", colnames(data), value=FALSE)]
b <- data[,grep("b$", colnames(data), value=FALSE)]
c <- data[,grep("C$", colnames(data), value=FALSE)]

a <- a[,order(colnames(a))]
b <- b[,order(colnames(b))]
c <- c[,order(colnames(c))]

# normalizing sample names
colnames(a) <- gsub("a$","",colnames(a))
colnames(b) <- gsub("b$","",colnames(b))
colnames(c) <- gsub("-Healthy.*C$","",colnames(c))
colnames(c) <- gsub("-NASH.*C$","",colnames(c))
colnames(c) <- gsub("-SS.*C$","",colnames(c))

sample.names <- unique(c(colnames(a),colnames(b),colnames(c)))
sample.names <- sample.names[which(sample.names%in%colnames(a))]
sample.names <- sample.names[which(sample.names%in%colnames(b))]
# ONLY RUN THE REPLICATION EXPERIMENT WITH A AND B - c has weird bias
#sample.names <- sample.names[which(sample.names%in%colnames(c))]

print("Retrieved sample names.")

# calculate proportion difference in read count (difference / average total read count)
readCountDiff <- function(sampleName,replicate.a,replicate.b) {
  sum.a <- sum(replicate.a[,which(colnames(replicate.a)==sampleName)])
  sum.b <- sum(replicate.b[,which(colnames(replicate.b)==sampleName)])
  return(abs(sum.a-sum.b)/mean(sum.a,sum.b))
}
difference <- sapply(sample.names, function(x) { return(readCountDiff(x,a,b)) } )

# the sample for the count replication study has the most closely replicated read count
# HLD-38 has a read count of 16397 in replicate a and 16389 in replicate b
best.sample <- sample.names[which(difference == min(difference))]

print(paste("Found best sample for replication experiment:", best.sample))

replicate.data <- data.frame(a[,which(colnames(a)==best.sample)], b[,which(colnames(b)==best.sample)])
colnames(replicate.data) <- c("Rep_A", "Rep_B")
# since replicates a and b were taken from the same table, the order of the rows (OTUs) should be the same
rownames(replicate.data) <- rownames(a)

rowsums <- replicate.data[,"Rep_A"] + replicate.data[,"Rep_B"]
replicate.data <- replicate.data[which(rowsums!=0),]

replicate.data[replicate.data==0] <- prior

dirichlet.data <- round(rdirichlet(100000, replicate.data[,"Rep_A"]) *  sum(replicate.data[,"Rep_A"]), round.digits)
colnames(dirichlet.data) <- rownames(replicate.data)
dirichlet.data[ dirichlet.data==0 ] <- prior

logtransform.dirichlet <- apply(dirichlet.data, 1, function(x){ abs(log2(x) - log2(replicate.data[,"Rep_A"]) ) })
logtransform.dirichlet.q <- apply(logtransform.dirichlet, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))}) # limit of variance we get withdirich
logproportions.dirichlet <- apply(dirichlet.data, 1, function(x){ log2(abs(x- replicate.data[,"Rep_A"]) ) })
logproportions.dirichlet[ is.nan(logproportions.dirichlet) ] <- 0
logproportions.dirichlet[ is.infinite(logproportions.dirichlet) ] <- 0
logproportions.dirichlet.q <- apply(logproportions.dirichlet, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

print("Completed dirichlet.")

# rarefy 1000 times
rarefy.replace.a <- matrix(data=NA,nrow=1000,ncol=nrow(replicate.data))
rarefy.replace.b <- matrix(data=NA,nrow=1000,ncol=nrow(replicate.data))
rarefy.a <- matrix(data=NA,nrow=1000,ncol=nrow(replicate.data))
rarefy.b <- matrix(data=NA,nrow=1000,ncol=nrow(replicate.data))

for (i in 1:1000) {
  rarefy.replace.a[i,] <- rarefy(rownames(replicate.data), replicate.data[,"Rep_A"], sum(replicate.data$Rep_A[replicate.data$Rep_A!=prior])/2, FALSE)
  rarefy.replace.a[i,] <- rarefy.replace.a[i,] * 2
  rarefy.replace.b[i,] <- rarefy(rownames(replicate.data), replicate.data[,"Rep_B"], sum(replicate.data$Rep_B[replicate.data$Rep_B!=prior])/2, FALSE)
  rarefy.replace.b[i,] <- rarefy.replace.b[i,] * 2
  # # Not using jackknife rarefaction
  # rarefy.a[i,] <- sample(replicate.data$Rep_A, size=length(replicate.data$Rep_A)/2, replace=FALSE)
  # rarefy.a[i,] <- rarefy.a[i,] * 2
  # rarefy.b[i,] <- sample(replicate.data$Rep_B, size=length(replicate.data$Rep_B)/2, replace=FALSE)
  # rarefy.b[i,] <- rarefy.b[i,] * 2
}

print("Completed 1000 rarefactions")

rarefy.replace.a[rarefy.replace.a == 0] <- prior
rarefy.replace.b[rarefy.replace.b == 0] <- prior

logtransform.rarefy.a <- apply(rarefy.replace.a, 1, function(x){ abs(log2(x) - log2(replicate.data[,"Rep_A"]) ) })
logtransform.rarefy.a.q <- apply(logtransform.rarefy.a, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

logproportions.rarefy.a <- apply(rarefy.replace.a, 1, function(x){ log2(abs(x - replicate.data[,"Rep_A"]) ) })
logproportions.rarefy.a[ is.nan(logproportions.rarefy.a) ] <- 0
logproportions.rarefy.a[ is.infinite(logproportions.rarefy.a) ] <- 0
logproportions.rarefy.a.q <- apply(logproportions.rarefy.a, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

logtransform.rarefy.b <- apply(rarefy.replace.b, 1, function(x){ abs(log2(x) - log2(replicate.data[,"Rep_B"]) ) })
logtransform.rarefy.b.q <- apply(logtransform.rarefy.b, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

logproportions.rarefy.b <- apply(rarefy.replace.b, 1, function(x){ log2(abs(x - replicate.data[,"Rep_B"]) ) })
logproportions.rarefy.b[ is.nan(logproportions.rarefy.b) ] <- 0
logproportions.rarefy.b[ is.infinite(logproportions.rarefy.b) ] <- 0
logproportions.rarefy.b.q <- apply(logproportions.rarefy.b, 1, function(x) {quantile(x,probs=c(0.05,0.4,0.99))})

#A estimated variance plotted vs A-B difference
mincount <- apply(replicate.data[,1:2], 1, min)
rat <- sum(replicate.data[,1]) / sum(replicate.data[,2])

print("Processed data for plotting")

pdf("count_replication.pdf")

plot( log2(mincount), abs(log2(abs(replicate.data[,"Rep_A"] - round( rat * replicate.data[,"Rep_B"], 1)))), xaxt="n", yaxt="n", main="Difference", pch=19, cex=1.1, col=rgb(0,0,0,0.3),  ylab="Absolute Difference", xlab="Feature Count", ylim=c(0,10), type="h", lwd=4)
points(loess.smooth(x=log2(replicate.data[,"Rep_A"]), y= logproportions.rarefy.a.q[3,]), type="l", col="red", lwd=2)
points(loess.smooth(x=log2(replicate.data[,"Rep_A"]), y= logproportions.rarefy.b.q[3,]), type="l", col="orange", lwd=2)
points(loess.smooth(x=log2(replicate.data[,"Rep_A"]), y= logproportions.dirichlet.q[3,]), type="l", col="blue", lwd=2)
axis(1, at=log2(c(1,10,100,1000,10000,100000)), labels=c(1,10,100,"1e3","1e4","1e5"))
axis(2, at=log2(c(1,10,100,1000)), labels=c(1,10,100,1000))

#rect(log2(1e3), -1, log2(1e5), 10.5,  col=rgb(0,0,0,0.15), lwd=NULL, border = NA)
#rect(2, -1, 5, 10.5,  col=rgb(.64,.2,.06,0.15), lwd=NULL, border = NA)


plot( log2(mincount), abs(log2(replicate.data[,"Rep_A"]) - log2(round(rat * replicate.data[,"Rep_B"], 1))), xaxt="n", yaxt="n", main="Ratio ", pch=19, cex=1.1, col=rgb(0,0,0,0.3),  ylab="Replicate Ratio", xlab="Feature Count", ylim=c(0,7), type="h", lwd=4)
points(loess.smooth(x=log2(replicate.data[,"Rep_A"]), y= logtransform.rarefy.a.q[3,]), type="l", col="red", lwd=2)
points(loess.smooth(x=log2(replicate.data[,"Rep_A"]), y= logtransform.rarefy.b.q[3,]), type="l", col="orange", lwd=2)
points(loess.smooth(x=log2(replicate.data[,"Rep_A"]), y= logtransform.dirichlet.q[3,]), type="l", col="blue", lwd=2)
#points(log2(d.21[,"Rep_A"]), ya.q[3,], pch=19, cex=0.5, col=rgb(0,0,1,0.4))
axis(1, at=log2(c(1,10,100,1000)), labels=c(1,10,100,"1e3"))
axis(2, at=log2(c(1,2,4,8,16,32,64,128)), labels=c(1,2,4,8,16,32, 64, 128))

dev.off()

print("Complete")
