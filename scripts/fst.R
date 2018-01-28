options(echo=TRUE)

user<-Sys.getenv("USER")
wkpath<-paste("/data/",user, "/fst_hmm",sep="")
setwd(wkpath)

# library to speed up loading of big tables
library(data.table)

# function to calculate Hudson's Fst
# --------------------------------------------------
hudsonFst<-function(locus=NA, p1=NA, p2=NA){
    numerator<-p1 * (1 - p1) + p2 * (1 - p2)
    denominator<-p1 * (1 - p2) + p2 * (1 - p1)
    fst<-1 - numerator/denominator
    out<-data.frame(locus,numerator,denominator,fst)
    return(out)
}
# --------------------------------------------------

# load file with allele frequencies
pops.af<-fread("timemaHVAxHVC.alfreq.txt", header=T, sep="\t")

# Calculate Fst for all loci
fst.HVAxHVC<-hudsonFst(locus=pops.af$locus, p1=pops.af$afHVA,p2=pops.af$afHVC)

# Write table to file
write.table(fst.HVAxHVC, file="timemaHVAxHVC.fst.dsv", 
                quote=F, row.names=F, sep="\t")

# plot FST histogram to png file
png(filename="timemaHVAxHVC.fst.png", width=1600, height=1200)
hist(fst.HVAxHVC$fst, breaks=100)
dev.off()

# Genome-wide FST statistics
# ------------------------------------------------------------------------------
# Calculate genome-wide FST as "ratio of averages" following Bhatia et al. 2013
nloci<-length(na.exclude(fst.HVAxHVC$fst))
fst.mean<-1-((sum(fst.HVAxHVC$numerator)/nloci)/(sum(fst.HVAxHVC$denominator)/nloci))
fst.mean

# "average of ratios" yields a noticeable lower value
mean(fst.HVAxHVC$fst,na.rm=T)

# median and upper quantiles for all loci
fst.quantile<-quantile(fst.HVAxHVC$fst, prob=c(0.5,0.95,0.99,0.999,0.9999), na.rm=T)
fst.quantile

fst.max<-max(fst.HVAxHVC$fst,na.rm=T)
fst.max

fst.genome<-c("genome",nloci, fst.mean,fst.quantile)
# ------------------------------------------------------------------------------

# FST stats by linkage group
# ------------------------------------------------------------------------------
# get linkage groups from loci names
lgs<-unique(gsub("_ord.*","", fst.HVAxHVC$locus))

# function to get fst (including median and quantiles) for a linkage group
fst.lg<-function(lg){
              fst<-fst.HVAxHVC[grep(paste("^",lg,"_",sep=""), fst.HVAxHVC$locus),]
              nloci<-length(na.exclude(fst$fst))
              fst.mean<-1-((sum(fst$numerator)/nloci)/(sum(fst$denominator)/nloci))
              fst.quantile<-quantile(fst$fst, prob=c(0.5,0.95,0.99,0.999,0.9999), names=F, na.rm=T)
              fst.max<-max(fst$fst, na.rm=T)
              return(c(nloci=nloci,mean=fst.mean,fst.quantile,fst.max))
            }

# calculate fst for each linkage group
fst.stats.lgs<-lapply(lgs, fst.lg)

# put them together in a data frame
fst.stats.lgs<-data.frame(lg=lgs,do.call("rbind",fst.stats.lgs),stringsAsFactors=F)
# ------------------------------------------------------------------------------

# summarize results
# ------------------------------------------------------------------------------
# put together whole genome and lgs estimates
fst.stats<-rbind(fst.genome,fst.stats.lgs)

# put names of columns
colnames(fst.stats)<-c("lg","nloci","mean","median","q95%","q99%","q99.9%","q99.99%","max")

# Write summary table to file
write.table(fst.stats,file="timemaHVAxHVC.lgs.fst.dsv", 
            quote=F, row.names=F, sep="\t")
# ------------------------------------------------------------------------------
