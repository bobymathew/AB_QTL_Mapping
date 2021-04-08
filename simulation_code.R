
################################33

#simulation

################################33
#select 5 QTLs for simulation
qtl=c(122,411,1021,1700,2730)

train=setdiff(1:2745,qtl)



library(rrBLUP)
#read the orginal genotypes
dat=read.table("B22_geno.txt",header=T,sep="\t",na.strings="NA")

M=as.matrix(dat[,4:253])

mark=dim(M)[1]
line=250
#set the heritability to 0.6
h2 <- 0.6

u <- rep(0,mark) #marker effects


set.seed(10)
#simulate the qtls effects from a uniform distribution
u[qtl] <- runif(length(qtl),9,10)


g <- as.vector(crossprod(M,u))
set.seed(1000)
y <- g + rnorm(line,mean=0,sd=sqrt((1-h2)/h2*var(g)))
pheno <- data.frame(line=colnames(M),y=y)

scores <- GWAS(pheno,dat,plot=TRUE)	


new_score=scores[order(-scores$y),]
head(new_score,15)

