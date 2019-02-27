options(echo=TRUE)

user<-Sys.getenv("USER")
wkpath<-paste("/data/",user, "/fst_hmm",sep="")
setwd(wkpath)

# library to fit Hidden Markov Models
library(HiddenMarkov)

# library to speed up loading of big tables
library(data.table)

# Modified version of Mstep.norm for fitting HMM below
# standard deviation is fixed to values specified in pm$sd
# ------------------------------------------------------------------------------
Mstep.normc <- function(x, cond, pm, pn){
  nms <- sort(names(pm))
  n <- length(x)
  m <- ncol(cond$u)
  if (all(nms==c("mean", "sd"))){
    mean <- as.numeric(matrix(x, nrow = 1) %*% cond$u)/apply(cond$u, 
       MARGIN = 2, FUN = sum)
    sd <- pm$sd
    return(list(mean=mean, sd=sd))
  }
}
rnormc <- rnorm
dnormc <- dnorm
pnormc <- pnorm
qnormc <- qnorm
# ------------------------------------------------------------------------------


# Prepare data
# ------------------------------------------------------------------------------
# Load Fst if need be
fst.HVAxHVC<-fread("timemaHVAxHVC.fst.dsv", header=T, sep="\t")

# Use only one linkage group for this example
fst<-fst.HVAxHVC[grep("^lg01_", fst.HVAxHVC$locus),]$fst

# Remove missing data
fst<-fst[!is.na(fst)]

# Replace too low Fst values by tractable value
fst[fst<1e-5]<-1e-5

# logit transform Fst, we will assume a normal logit distribution for each state
lfst<-log(fst/(1-fst))

# nloci, mean and sd (used below)
nloci<-length(lfst)
mean.lfst<-mean(lfst)
sd.lfst<-sd(lfst)
# ------------------------------------------------------------------------------


# Fit homogenous HMM of 3 discrete states
# ------------------------------------------------------------------------------
# initial m x m transition probability matrix
# disable transitions from 1 to 3 and 3 to 1 (set value to zero)
# (direct transitions from low to high and high to low differentiation)
Pi<-matrix(c(0.9,0.1,0,
             0.05,0.9,0.05,
             0,0.1,0.9),nrow=3,byrow=T)
# it is strongly recommended to run the EM algorithm multiple times (hundreds
# or thousands) with different starting matrices to avoid maximum-likelihood
# local optima; e.g. code to produce 100 random matrices:
# random.matrix<-function(x){
#                   x<-runif(3)
#                   matrix(c(x[1],1-x[1],0,
#                            x[2]/2,1-x[2],x[2]/2,
#                            0,x[3],1-x[3]),nrow=3,byrow=T)}
# matrices<-lapply(1:100, random.matrix)

# initial marginal probability distribution of the m hidden states
delta<-c(0,1,0)

# initial list of parameter values associated with the distribution 
# of observed process
# they must coincide with the parameters of the function use
# we are using a normal and must give 3 means and 3 sds
pm<-list(mean=c(mean.lfst*0.9,mean.lfst,mean.lfst*1.1),
         sd=rep(sd.lfst,3))

# set discrete time HMM object and initial values
# normc is modified distribution so that variance is fixed
# to genome-wide estimate (sd.lfst)
init<-dthmm(x=lfst,
            Pi=Pi,
            delta=delta,
            distn="normc",
            pm=pm,
            discrete=FALSE)

# set maximum number of iterations and convergence criterion
ctrl<-bwcontrol(maxiter=1000, tol=1e-6)

# estimate parameters of HMM using Baum-Welch EM algorithm
# this may take a while (a few minutes)
paramHMM<-BaumWelch(init, ctrl)

# predict Markov states using Viterbi algorithm
fitHMM<-Viterbi(paramHMM)
# ------------------------------------------------------------------------------

# Summarize results
# ------------------------------------------------------------------------------
# mean Fst of states
fst.states<-round(1/(1+exp(-1 * paramHMM$pm$mean)),5)

# number of loci in each state
nloci.states<-table(fitHMM)
nloci.states

# proportion of loci in each state
proploci.states<-nloci.states/nloci
proploci.states

# number of contiguous regions for each state
nregions.states<-unlist(lapply(1:3, 
                    function(x) sum((fitHMM==x)[1:(nloci-1)]-(fitHMM==x)[2:nloci]==1)))
nregions.states

summary.table<-cbind(state=c(1,2,3),
               mean.fst=fst.states,
               nloci=nloci.states,
               proploci=proploci.states,
               nregions=nregions.states)
# show table
summary.table

# save to file
write.table(summary.table, file="fst_hmm_lg01_regions_summary.dsv", 
              quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------

# plot regions using differentiation legend
# Plot to a png file because the number of dots is very high
# drawing this kind of plot over the network is very slow
# also opening vectorial files with many objects is slow
# ------------------------------------------------------------------------------
png("fst_hmm_lg01.png", width=4000, height=1000)
par(mar=c(10,12,10,4)+0.1)
colours<-c("red","dark grey","blue")
plot(fst, col=colours[fitHMM], ylim=c(0,0.8), 
     xlab="relative position", ylab=expression("F"[ST]), main="regions of differentiation", 
     cex=1, cex.axis=4, cex.lab=4, cex.main=4, mgp=c(7,3,0))

# add horizontal line to illustrate top 0.1% FST values
abline(h=quantile(fst, prob=0.999,na.rm=T), col="black", lty=3, cex=4, lwd=3)
text (0,quantile(fst, prob=0.999,na.rm=T)+0.025, label=expression("">=" 0.1%"), adj=1, pos=2, cex=3)

# legend
legend("topright",legend=c("high", "medium", "low"), 
       pch=1, col=colours, title="differentiation", bty="n", cex=4)
dev.off()
# ------------------------------------------------------------------------------

# Calculate size of HMM regions in base pairs
# ============

# function to get size of regions
# ------------------------------------------------------------------------------
get.regions<-function(lb=NA,ub=NA, loci.ord=NA, loci.pos=NA, lgordscalen=NA){
  nregions<-length(lb)
  bpsize<-numeric(nregions)
  for (i in 1:nregions){
    if (is.na(loci.ord[lb[i]]) || is.na(loci.ord[ub[i]])){
      cat ("I: ", i, "\n")
    }
    if (loci.ord[lb[i]] == loci.ord[ub[i]]){ # same scaffold
      bpsize[i]<-loci.pos[ub[i]]-loci.pos[lb[i]]
    }
    else{
      bpsize[i]<-(lgordscalen[lgordscalen$ord==loci.ord[lb[i]],]$length
                  -loci.pos[lb[i]]
                  +sum(lgordscalen[(lgordscalen$ord>loci.ord[lb[i]] & 
                                    lgordscalen$ord<loci.ord[ub[i]]),]$length)
                  +loci.pos[ub[i]])
    }
  }
  
  regions<-data.frame(lb.ord=loci.ord[lb], ub.ord=loci.ord[ub],
                      lb.sca=loci.sca[lb], ub.sca=loci.sca[ub],
                      lb.pos=loci.pos[lb], ub.pos=loci.pos[ub],
                      length=bpsize)
  return(regions)
}
# ------------------------------------------------------------------------------

# load order and size of scaffolds in draft genome
# ------------------------------------------------------------------------------
lgordscalen<-read.table("data/lg_ord_sca_length.dsv",header=T,sep="\t")
lgordscalen<-lgordscalen[lgordscalen$lg==1,]

loci.pos<-as.numeric(gsub(".*:","",fst.HVAxHVC$locus))
loci.ord<-as.numeric(gsub(".*ord|_scaf.*","",fst.HVAxHVC$locus))
loci.sca<-as.numeric(gsub(".*_scaf|:.*","",fst.HVAxHVC$locus))
# ------------------------------------------------------------------------------

# high differentiation regions
# ------------------------------------------------------------------------------
lb.locindex<-which(fitHMM[2:nloci]==1 & fitHMM[1:(nloci-1)]!=1)
ub.locindex<-which(fitHMM[1:(nloci-1)]==1 & fitHMM[2:nloci]!=1)
if (length(lb.locindex)<length(ub.locindex)) lb.locindex<-c(1,lb.locindex)
if (length(lb.locindex)>length(ub.locindex)) ub.locindex<-c(ub.locindex,nloci)

hidif.regions<-get.regions(lb.locindex,ub.locindex,loci.ord,loci.pos,lgordscalen)
write.table(hidif.regions, file="fst_hmm_lg01_hidif_regions.dsv", quote=F, sep="\t")
# ------------------------------------------------------------------------------

# medium differentiation regions
# ------------------------------------------------------------------------------
lb.locindex<-which(fitHMM[2:nloci]==2 & fitHMM[1:(nloci-1)]!=2)
ub.locindex<-which(fitHMM[1:(nloci-1)]==2 & fitHMM[2:nloci]!=2)
if (length(lb.locindex)<length(ub.locindex)) lb.locindex<-c(1,lb.locindex)
if (length(lb.locindex)>length(ub.locindex)) ub.locindex<-c(ub.locindex,nloci)

medif.regions<-get.regions(lb.locindex,ub.locindex,loci.ord,loci.pos,lgordscalen)
write.table(medif.regions, file="fst_hmm_lg01_medif_regions.dsv", quote=F, sep="\t")
# ------------------------------------------------------------------------------

# low differentiation regions
# ------------------------------------------------------------------------------
lb.locindex<-which(fitHMM[2:nloci]==3 & fitHMM[1:(nloci-1)]!=3)
ub.locindex<-which(fitHMM[1:(nloci-1)]==3 & fitHMM[2:nloci]!=3)
if (length(lb.locindex)<length(ub.locindex)) lb.locindex<-c(1,lb.locindex)
if (length(lb.locindex)>length(ub.locindex)) ub.locindex<-c(ub.locindex,nloci)

lodif.regions<-get.regions(lb.locindex,ub.locindex,loci.ord,loci.pos,lgordscalen)
write.table(lodif.regions, file="fst_hmm_lg01_lodif_regions.dsv", quote=F, sep="\t")
# ------------------------------------------------------------------------------

# Compare distribution of sizes
# ------------------------------------------------------------------------------
png("fst_hmm_lg01_dif_regions_sizes.png", width=3000, height=1000)
par(mfrow=c(1,3))
par(mar=c(15,15,10,4)+0.1,  mgp=c(7,3,0))
hist(hidif.regions$length/1000,main="high differentiation regions", 
        xlim=c(0,300),xlab="size (Kbp)",breaks=50, cex.axis=4, cex.lab=4, cex.main=4)
hist(medif.regions$length/1000,main="medium differentiation regions", 
        xlim=c(0,300), xlab="size (Kbp)",breaks=50, cex=4, cex.axis=4, cex.lab=4, cex.main=4)
hist(lodif.regions$length/1000,main="low differentiation regions", 
        xlim=c(0,300), xlab="size (Kbp)", breaks=50, cex=4, cex.axis=4, cex.lab=4, cex.main=4)
dev.off()
# ------------------------------------------------------------------------------
# ==============================================================================
