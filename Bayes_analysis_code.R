
#Set the number iterations for the MCMC chain with the burnin and thinning period 
ite=50000;burnin=10000;th=50

#read the Genotype file with the marker information 
geno=read.table("Z86_genotype.txt",header=T,sep="\t")

#Read the phenorype file
pheno=read.table("Z86_pheno.txt",header=T,sep="\t")

#standardize the phenotype file
ytot = (pheno[,2]-mean(pheno[,2]))/sd(pheno[,2])



#Marker information as matrix for MCMC
M = as.matrix(geno[,4:dim(geno)[2]])


map=geno[,1:3]

Sys.time()


mark=dim(M)[1]
line=dim(M)[2]


y=as.vector(ytot)
A=vector()
va=vector()
#marker varaince
#fvb=matrix(0,nrow=mark,ncol=ite)
t=1


A[1]=mean(y)
B=matrix(0,nrow=mark,ncol=2) #marker x 2
vb=matrix(0,nrow=mark,ncol=2) #marker x 2
# set first values as 0
B[,1]=rep(0,mark)
#initial values 
vb[,1]=0.8*rep(1,mark)
va[1]=0.5
#final values
fA=vector()
fva=vector()
#marker effect
#fB=matrix(0,nrow=mark,ncol=ite)

res=matrix(0,nrow=mark,ncol=1)


for (s in 1:ite ) 
	{
	#update A(mean)
	mA=mean(y-t(B[,t])%*%M)
	vA=va[t]/line
	A[t+1]=mA+sqrt(vA)*rnorm(1)

	#update B (marker effect)
	B[,t+1]=B[,t]
	Bmxj=t(B[,t+1])%*%M

	for (j in 1:mark ) 
		{	
		Bmxj=Bmxj-B[j,t+1]*M[j,] # to keep k not j
		e=sum(M[j,]*(y-A[t+1]-Bmxj))

		d1=sum(M[j,]^2) #Sig^2_ij

		mB=vb[j,t]*e/(vb[j,t]*d1+va[t])

		vB=vb[j,t]*va[t]/(vb[j,t]*d1+va[t])


		B[j,t+1]=mB+sqrt(vB)*rnorm(1)
		Bmxj=Bmxj+B[j,t+1]*M[j,]
		}
# update va the residual error
	c2=rchisq(1,line)
	e2=sum((y-A[t+1]-t(B[,t+1])%*%M)^2)
	va[t+1]=e2/c2
#update vb
	vb[,t+1]=(B[,t+1]^2)/rchisq(mark,1)
#store data
	fA[s]=A[t+1]
	fva[s]=va[t+1]
	if((s>=burnin) & (s %% th  ==0)) 
		{
	res=res+B[,t+1]
		}
	
	#fvb[,s]=vb[,t+1]
#next round initial values	
	A[t]=A[t+1]
	B[,t]=B[,t+1]
	va[t]=va[t+1]
	vb[,t]=vb[,t+1]

	}
Sys.time()


m_eff=res/((ite-burnin)/th)


