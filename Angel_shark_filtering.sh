## New filtering with minor allele filter removed from SOL tree ##

#Split Multiallelic Freebayes calls
{```{bash}```
vcfallelicprimitives -k -g TotalRawSNPs.vcf | sed 's:\.|\.:\.\/\.:g' > TRS.prim
} #notepad cleanup

#Remove indels
{```{bash}```
vcftools --vcf TRS.prim --recode-INFO-all --recode --out SNP.TRS --remove-indels

#After filtering, kept 69 out of 69 Individuals
#After filtering, kept 601882 out of a possible 641156 Sites
#Run Time = 232.00 seconds
} #notepad cleanup

#Basic quality filters of data (Minimum of 2 alleles, Depth of at least 5, Minimum quality of at least 20)
{```{bash}```
vcftools --vcf SNP.TRS.recode.vcf --out SNP.TRS.QC --recode --recode-INFO-all --min-alleles 2 --minDP 5 --minQ 20 --min-meanDP 20

#After filtering, kept 69 out of 69 Individuals
#After filtering, kept 169844 out of a possible 601882 Sites
#Run Time = 78.00 seconds

charts.sh SNP.TRS.QC.recode.vcf &
} #notepad cleanup

#Filters allelic balance, quality vs depth, strand representation and paired read representation
{```{bash}```
dDocent_filters SNP.TRS.QC.recode.vcf SNP.TRS.dDocent
##Interactive responses
#no
#100000

#Remaining sites
# 59639
} #notepad cleanup

#Filtering singleton and doubleton loci for depth (Singletons >20 reads; Doubletons > 10 reads)
{```{bash}```
vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out out --singletons

awk ' $3=="S" {print $1, $2}' out.singletons > sing.loci
awk ' $3=="D" {print $1, $2}' out.singletons > doub.loci

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01a --recode --recode-INFO-all --exclude-positions sing.loci
vcftools --vcf SNP.TRS.F01a.recode.vcf --out SNP.TRS.F01b --recode --recode-INFO-all --exclude-positions doub.loci

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01.sing --recode --recode-INFO-all --positions sing.loci
vcftools --vcf SNP.TRS.F01.sing.recode.vcf --out SNP.TRS.F02.sing --recode --recode-INFO-all --minDP 20

vcftools --vcf SNP.TRS.dDocent.FIL.recode.vcf --out SNP.TRS.F01.doub --recode --recode-INFO-all --positions doub.loci
vcftools --vcf SNP.TRS.F01.doub.recode.vcf --out SNP.TRS.F02.doub --recode --recode-INFO-all --minDP 10

vcf-concat SNP.TRS.F02.sing.recode.vcf SNP.TRS.F02.doub.recode.vcf SNP.TRS.F01b.recode.vcf > SNP.TRS.F02.vcf

rm out.singletons
} #notepad cleanup

#Filter loci that have high variation in depth across a locus with an individual
{```{bash}```
vcftools --vcf SNP.TRS.F02.vcf --out out --geno-depth
#After filtering, kept 69 out of 69 Individuals
#After filtering, kept 59663 out of a possible 59663 Sites
#Run Time = 2.00 seconds
} #notepad cleanup
{```{R}```
gdepth<-read.table(file="out.gdepth", head=T)
gdepth[gdepth==-1]<-NA

for (i in 3:dim(gdepth)[2]) {
temp<-aggregate(gdepth[,i],by=list(gdepth[,1]), sd)
if(i==3){indv.site.sd<-data.frame(temp,row.names=1)} else
{indv.site.sd[,(i-2)]<-temp[,2]}
}
colnames(indv.site.sd)<-colnames(gdepth[3:dim(gdepth)[2]])
tmp<-apply(indv.site.sd, 1, mean, na.rm=T)
tmp2<-unique(c(names(which(tmp>25))))
length(tmp)
length(tmp2)
write.table(tmp2,file="bad.loci.sd", quote=F, col.names=F, row.names=F)
q("no")
}} #Notepad cleanup
{```{bash}```
grep "dDocent" SNP.TRS.F02.vcf | cut -f 1,2 | uniq | tail -n +2 > contigs.txt
grep -wf bad.loci.sd contigs.txt > bad.loci
vcftools --vcf SNP.TRS.F02.vcf  --out SNP.TRS.F03 --exclude-positions bad.loci --recode-INFO-all --recode
#After filtering, kept 69 out of 69 Individuals
#After filtering, kept 59365 out of a possible 59663 Sites
#Run Time = 23.00 seconds

charts.sh SNP.TRS.F03.recode.vcf &
} #notepad cleanup

#Filtering Tree without MAF
{```{bash}```
cd /home/afields/Workspace/Angels/filtering_Sep2023
cp ../analysis/Aug2020/SNP.TRS.F03.recode.vcf ../analysis/Aug2020/popmap .
mkdir SOL
cd SOL 
cp -s ../SNP.TRS.F03.recode.vcf .
SOL.filter1.no_mac.sh SNP.TRS.F03.recode.vcf
#Looked at in excel for best filter balance
cp vcf/B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf ..
cd ..

charts.sh B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf
} #notepad cleanup

#Removing loci with super high depth
{```{R}```
dat <- read.table("B.1.3.3.2.3.2.1.2.SNP.finalb.1.2023-09-20/loci.csv", head=T, sep=",")
hist(dat$MEAN_DEPTH, breaks=100, col="red4")
#The central portion of this graphic looks like it should be symetrical

#finding the highest point in the density plot
peak <- which.max(density(dat$MEAN_DEPTH)$y)
peak
#[1] 253
density(dat$MEAN_DEPTH)$x[peak]
#[1] 64.67747
abline(v = density(dat$MEAN_DEPTH)$x[peak])

## Method 1 ##
#Min after the peak 
tmp.min <- min(density(dat$MEAN_DEPTH)$y[peak:400])

#set cutoff
cutoff <- density(dat$MEAN_DEPTH)$x[density(dat$MEAN_DEPTH)$y == tmp.min]
cutoff
#[1] 379
abline(v = cutoff)

## Method 2 ##
#Symetrical around the peak

#Minimum depth
min(dat$MEAN_DEPTH)
#[1] 20.4615

#set cutoff
cutoff <- 2*density(dat$MEAN_DEPTH)$x[peak] - min(dat$MEAN_DEPTH)
cutoff
#[1] 108.8934
abline(v = cutoff)

write.table(dat[dat$MEAN_DEPTH > cutoff, c("CHROM","POS")], file="high_depth_loci.txt", quote=F, col.names=F, row.names=F, sep="\t")
} #notepad cleanup

{```{bash}```
vcftools --vcf B.1.3.3.2.3.2.1.2.SNP.finalb.1.recode.vcf --out SNP.TRS.F04 --recode --recode-INFO-all --exclude-positions high_depth_loci.txt
#After filtering, kept 65 out of 65 Individuals
#After filtering, kept 35311 out of a possible 35323 Sites
#Run Time = 14.00 seconds

charts.sh SNP.TRS.F04.recode.vcf
} #notepad cleanup

#MAF filtering
{```{R}```
#Load Libraries
{```{R}```
library('vcfR')
library('adegenet')
library('dartR')
library('related')
library('VennDiagram')
} #notepad cleanup

#Load data
{```{R}```
#Load sample information
strata <- read.csv(file = "SNP.TRS.F04.2023-09-22/indv.csv", header = TRUE)
#Adding Library, Index and Well to the sample information
strata$INDV <- as.character(as.matrix(strata$INDV))
indv.list <- strsplit(strata$INDV, "_")
tmp.df <- data.frame(matrix(ncol=5))
for(i in 1:nrow(strata)){tmp.df <- rbind(tmp.df, c(strata$INDV[i],indv.list[[i]]))}
tmp.df <- tmp.df[-1, ]
for(i in 13:16){strata[,i] <- tmp.df[match(strata$INDV, tmp.df[,1]),(i-11)]}
#Checking the sample information
head(strata)

#Loading the vcf file
vcf <- read.vcfR(file="SNP.TRS.F04.recode.vcf")
gen.vcf<-vcfR2genind(vcf)
#Assing the sample information to the genind object
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
#Checking the sample information of the genind
head(strata(gen.vcf))
} #notepad cleanup

#PCA to compare groups
{```{R}```
X <- scaleGen(gen.vcf, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
par(mfrow=c(1,1))
setPop(gen.vcf) <- ~POP
ade4::s.class(pca1$li, pop(gen.vcf), col=c("red4","mediumblue","darkgreen", "chocolate"), cstar=0, cellipse=0)

setPop(gen.vcf) <- ~Lib
ade4::s.class(pca1$li, pop(gen.vcf), col=c("red4","mediumblue","darkgreen", "chocolate"), cstar=0, cellipse=0)
} #notepad cleanup

#Removing Library effects
{```{R}```
setPop(gen.vcf) <- ~Lib
tmp<-gl.outflank(gen.vcf, qthreshold = 0.05)
length(which(tmp$outflank$results[15]=="TRUE"))	
#2428

outliers <- tmp$outflank$results[which(tmp$outflank$results[15]=="TRUE"),1]
out.list <- dput(matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1])

set.loc <- locNames(gen.vcf)[which(!(locNames(gen.vcf) %in% out.list))]
gen2.vcf <- gen.vcf[, loc=set.loc, drop=T]
set.loc <- locNames(gen.vcf)[which((locNames(gen.vcf) %in% out.list))]
gen.out.vcf <- gen.vcf[, loc=set.loc, drop=T]

X <- scaleGen(gen.vcf, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.net <- scaleGen(gen2.vcf, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.out <- scaleGen(gen.out.vcf, NA.method="mean", scale=F)
pca.out <- dudi.pca(X.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen.vcf) <- setPop(gen2.vcf) <- setPop(gen.out.vcf) <- ~Lib
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen.vcf), col=funky(4), cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen.vcf))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen2.vcf), col=funky(4), cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen2.vcf))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen.out.vcf), col=funky(4), cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.out.vcf))," loci)", sep=""), 3, 2, adj = 0.95)
} #notepad cleanup

#PCA to compare groups
{```{R}```
par(mfrow=c(1,1))
setPop(gen2.vcf) <- ~POP
ade4::s.class(pca.net$li, pop(gen2.vcf), col=funky(4), cstar=0, axesell=F)
} #notepad cleanup

#Removing duplicates
{```{R}```
gl <- gi2gl(gen.vcf)
test.cov <- gl2related(gl, save=F)
dim(test.cov)
#[1]   65 70623
wang.list <- coancestry(test.cov, wang =1)

par(mfrow=c(3,1))
hist(wang.list$relatedness$wang, breaks=200, main="All Relationships", col="red4")
hist(wang.list$relatedness$wang, breaks=200, main="All Relationships", ylim=c(0,1000), col="red4")

dups.rel2 <- wang.list$relatedness[wang.list$relatedness$wang > 0.4,2:3]
strata[grepl(dups.rel2[1,1],strata$INDV) | grepl(dups.rel2[1,2],strata$INDV), ]
{ #Results
                       INDV N_SITES.x MEAN_DEPTH N_DATA N_GENOTYPES_FILTERED
12  Atl.AS.124_Lib1_I12_D05     35315    63.6532  35323                    0
13 Atl.AS.124._Lib2_I12_A03     35312   179.5530  35323                    0
   N_MISS      F_MISS O.HOM.  E.HOM. N_SITES.y      Fis POP      Sample  Lib
12     12 0.000339722  26056 27415.4     33323 -0.23012 Atl  Atl.AS.124 Lib1
13     14 0.000396342  26459 27413.5     33321 -0.16157 Atl Atl.AS.124. Lib2
   Index Row
12   I12 D05
13   I12 A03
}
strata[grepl(dups.rel2[2,1],strata$INDV) | grepl(dups.rel2[2,2],strata$INDV), ]
{ #Results
                       INDV N_SITES.x MEAN_DEPTH N_DATA N_GENOTYPES_FILTERED
43 EAST.AS.9.2_Lib1_I12_D04     35313    65.6581  35323                    0
44 EAST.AS.9.2_Lib2_I12_A02     35314   154.2450  35323                    0
   N_MISS      F_MISS O.HOM.  E.HOM. N_SITES.y      Fis  POP      Sample  Lib
43     22 0.000622824  25941 27409.2     33315 -0.24859 EAST EAST.AS.9.2 Lib1
44     10 0.000283102  26481 27417.2     33325 -0.15847 EAST EAST.AS.9.2 Lib2
   Index Row
43   I12 D04
44   I12 A02
}

rm.indv <- c("Atl.AS.124_Lib1_I12_D05","EAST.AS.9.2_Lib1_I12_D04")
set.indv <- indNames(gen.vcf)[which(!(indNames(gen.vcf) %in% rm.indv))]
set.loc <- locNames(gen.vcf)[which(!(locNames(gen.vcf) %in% out.list))]
gen2.vcf <- gen.vcf[set.indv, loc=set.loc, drop=T]
} #notepad cleanup

#PCA to compare groups
{```{R}```
par(mfrow=c(1,1))
X2 <- scaleGen(gen2.vcf, NA.method="mean", scale=F)
pca2 <- dudi.pca(X2,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen2.vcf) <- ~POP
ade4::s.class(pca2$li, pop(gen2.vcf), col=funky(4), cstar=0, axesell=F)
} #notepad cleanup

#Adding kmeans groups to strata file
{```{R}```
setPop(gen2.vcf) <- ~POP
grp4 <- find.clusters(gen2.vcf, max.n.clust=40, n.pca=400, n.clust =4, method="kmeans")
table(pop(gen2.vcf), grp4$grp)
{ #Results
        1  2  3  4
  Atl  17  0  0  0
  EAST  0 24  1  0
  Pac   0  0  0  1
  WEST  0  0 20  0
}

tmp <- data.frame(grp4$grp)[match(gen2.vcf@strata$INDV,names(grp4$grp)),]
gen2.vcf@strata$G4 <- as.numeric(tmp)
} #notepad cleanup

#PCA to compare groups
{```{R}```
setPop(gen2.vcf) <- ~G4
par(mfrow=c(1,1))
ade4::s.class(pca2$li, pop(gen2.vcf), col=c("red4","mediumblue","darkgreen","chocolate"), cstar=0, cellipse=0)
} #notepad cleanup

#Looking for minor alleles in each group
{```{R}```
ESU <- seppop(gen2.vcf)

minor.all <- list()
for(i in 1:length(ESU)){
gen.tmp <- ESU[[i]]
tmp.v <- apply(gen.tmp@tab, 2, function(x) sum(x,na.rm=T)/(2*length(which(!is.na(x)))))
minor.all[[i]] <- names(tmp.v)[tmp.v < 0.01]
}
minor.all <- lapply(minor.all, function(x) x[!is.na(x)])

lapply(minor.all, length)
{ #Results
[[1]]	#Atlantic
[1] 9494

[[2]]	#East Gulf
[1] 5028

[[3]]	#Pacific
[1] 25086

[[4]]	#West Gulf
[1] 7162
}
} #notepad cleanup

#Plotting a Venn diagram of the overlap in minor alleles
{```{R}```
POP.tab <- lapply(ESU, function(x) table(x$strata$POP))
list.names <- unlist(lapply(POP.tab, function(x) names(x)[x == max(x)]))
venn.diagram(minor.all,
category.names=c(
paste(list.names[1],"\nloci\n(n=", length(minor.all[[1]]), ")", sep = ""),
paste(list.names[2],"\nloci\n(n=", length(minor.all[[2]]), ")", sep = ""),
paste(list.names[3],"\nloci\n(n=", length(minor.all[[3]]), ")", sep = ""),
paste(list.names[4],"\nloci\n(n=", length(minor.all[[4]]), ")", sep = "")
), "Venn_min_all.tiff", main="Minor alleles by species", fill=c("red4", "mediumblue", "darkgreen", "chocolate"), alpha=rep(0.5,length(list.names)), cat.cex=rep(0.8,length(list.names))
)
} #notepad cleanup

#Getting the overlap between lists of minor alleles
#Needs to be scaled to the number of groups
#currently set for 3 groups
{```{R}```
x <- minor.all
A <- x[[1]]
B <- x[[2]]
C <- x[[3]]
D <- x[[4]]
nab <- intersect(A, B)
nabc <- intersect(nab, C)
length(nabc)
#1005
#Group 4 is Pacific and therefore not the same question
nabcd <- intersect(nabc, D)
length(nabcd)
#0

min_all.m <- matrix(unlist(strsplit(nabc, "[.]")), ncol=2, byrow=T)
loc.parts <- matrix(unlist(strsplit(min_all.m[,1], "_")), ncol=5, byrow=T)
loc.df <- data.frame(CHROM=apply(loc.parts[,1:3],1,function(x) paste(x, collapse="_")), POS=loc.parts[,4], Allele=min_all.m[,2], UKN=loc.parts[,5])

vcf.rec <- vcf

#Replacing low maf values with missing data indicators
for(i in 1:length(nabc)){
#Record of progress
if(i/50 == round(i/50,0)){print(paste("Processing locus number", i, "out of", length(nabc)))}
#Pulling Locus
tmp.loc <- which(vcf@fix[,"CHROM"] == loc.df$CHROM[i] & vcf@fix[,"POS"] == loc.df$POS[i])
#Pulling genotype calls
tmp.gt <- vcf@gt[tmp.loc,]
#Seperating the vcf information for each individual
tmp.list <- strsplit(tmp.gt,":")
#Replacing minor alleles
for(j in 2:length(tmp.list)){if(length(grep(loc.df$Allele[i], tmp.list[[j]][1])) > 0){tmp.list[[j]][1] <- ".|."}}
#Reassembling vcf information
tmp.gt <- lapply(tmp.list, function(x) paste(x,collapse=":"))
#Putting it back in the vcf object
vcf@gt[tmp.loc,] <- unlist(tmp.gt)[match(colnames(vcf@gt), names(tmp.gt))]
}

write.vcf(vcf, file="SNP.TRS.F05.vcf.gz")
}}} #notepad cleanup
} #notepad cleanup

#Applying an MAC with vcftools
{```{bash}```
vcftools --vcf SNP.TRS.F04.recode.vcf --out SNP.TRS.F05.mac --recode --recode-INFO-all --mac 3
#After filtering, kept 65 out of 65 Individuals
#After filtering, kept 27245 out of a possible 35311 Sites
#Run Time = 11.00 seconds


gzip -d SNP.TRS.F05.vcf.gz
charts.sh SNP.TRS.F05.vcf
#
charts.sh SNP.TRS.F05.mac.recode.vcf
#

#Removing non-variable SNPs after removing minor alleles?
vcftools --vcf SNP.TRS.F05.vcf --out SNP.TRS.F05a --recode --recode-INFO-all --min-alleles 2
#After filtering, kept 65 out of 65 Individuals
#After filtering, kept 35311 out of a possible 35311 Sites
#Run Time = 13.00 seconds
} #notepad cleanup
{```{R}```
dat4 <- read.table("SNP.TRS.F04.loci")
dat5m <- read.table("SNP.TRS.F05.mac.loci")
dat5a <- read.table("SNP.TRS.F05a.loci")

nrow(dat4)
#[1] 9406
nrow(dat5m)
#[1] 8688
nrow(dat5a)
#[1] 9199
length(which(dat5a[,1] %in% dat4[,1]))
#[1] 9199
length(which(dat5m[,1] %in% dat4[,1]))
#[1] 8688
length(which(dat5m[,1] %in% dat5a[,1]))
#[1] 8670
dat5m[which(!dat5m[,1] %in% dat5a[,1]),]
{ #Results
 [1] dDocent_Contig_10891 dDocent_Contig_12993 dDocent_Contig_13845
 [4] dDocent_Contig_13864 dDocent_Contig_13961 dDocent_Contig_16524
 [7] dDocent_Contig_17283 dDocent_Contig_4197  dDocent_Contig_49269
[10] dDocent_Contig_5805  dDocent_Contig_6606  dDocent_Contig_7289
[13] dDocent_Contig_7773  dDocent_Contig_7929  dDocent_Contig_8838
[16] dDocent_Contig_8979  dDocent_Contig_9442  dDocent_Contig_9611
}
dat5a[which(!dat5a[,1] %in% dat5m[,1]),]
{ #Results
  [1] dDocent_Contig_10003 dDocent_Contig_10022 dDocent_Contig_10031
  [4] dDocent_Contig_10061 dDocent_Contig_10077 dDocent_Contig_10079
  [7] dDocent_Contig_10109 dDocent_Contig_10124 dDocent_Contig_10129
 [10] dDocent_Contig_10182 dDocent_Contig_10254 dDocent_Contig_10277
 [13] dDocent_Contig_10319 dDocent_Contig_10348 dDocent_Contig_10355
 [16] dDocent_Contig_10366 dDocent_Contig_10404 dDocent_Contig_10418
 [19] dDocent_Contig_10429 dDocent_Contig_10432 dDocent_Contig_10437
 [22] dDocent_Contig_10448 dDocent_Contig_10449 dDocent_Contig_10469
 [25] dDocent_Contig_10493 dDocent_Contig_10501 dDocent_Contig_10507
 [28] dDocent_Contig_10510 dDocent_Contig_10564 dDocent_Contig_10694
 [31] dDocent_Contig_10695 dDocent_Contig_10750 dDocent_Contig_10807
 [34] dDocent_Contig_10816 dDocent_Contig_10832 dDocent_Contig_10837
 [37] dDocent_Contig_10848 dDocent_Contig_10850 dDocent_Contig_10875
 [40] dDocent_Contig_10877 dDocent_Contig_10887 dDocent_Contig_10889
 [43] dDocent_Contig_10894 dDocent_Contig_10943 dDocent_Contig_10992
 [46] dDocent_Contig_10994 dDocent_Contig_11003 dDocent_Contig_11014
 [49] dDocent_Contig_11021 dDocent_Contig_11023 dDocent_Contig_11024
 [52] dDocent_Contig_11032 dDocent_Contig_11055 dDocent_Contig_11103
 [55] dDocent_Contig_11106 dDocent_Contig_11153 dDocent_Contig_11181
 [58] dDocent_Contig_11182 dDocent_Contig_11216 dDocent_Contig_11218
 [61] dDocent_Contig_11219 dDocent_Contig_11220 dDocent_Contig_11230
 [64] dDocent_Contig_11242 dDocent_Contig_11244 dDocent_Contig_11287
 [67] dDocent_Contig_11300 dDocent_Contig_11301 dDocent_Contig_11311
 [70] dDocent_Contig_11330 dDocent_Contig_11379 dDocent_Contig_11384
 [73] dDocent_Contig_11385 dDocent_Contig_11389 dDocent_Contig_114
 [76] dDocent_Contig_11412 dDocent_Contig_11416 dDocent_Contig_11472
 [79] dDocent_Contig_11491 dDocent_Contig_11493 dDocent_Contig_11531
 [82] dDocent_Contig_11561 dDocent_Contig_11578 dDocent_Contig_11579
 [85] dDocent_Contig_11582 dDocent_Contig_11592 dDocent_Contig_11630
 [88] dDocent_Contig_11657 dDocent_Contig_11682 dDocent_Contig_11767
 [91] dDocent_Contig_11813 dDocent_Contig_11867 dDocent_Contig_11885
 [94] dDocent_Contig_11891 dDocent_Contig_11916 dDocent_Contig_12001
 [97] dDocent_Contig_12024 dDocent_Contig_12025 dDocent_Contig_12032
[100] dDocent_Contig_12041 dDocent_Contig_12043 dDocent_Contig_12048
[103] dDocent_Contig_12051 dDocent_Contig_12089 dDocent_Contig_12106
[106] dDocent_Contig_12111 dDocent_Contig_12140 dDocent_Contig_12146
[109] dDocent_Contig_12149 dDocent_Contig_12187 dDocent_Contig_12295
[112] dDocent_Contig_12309 dDocent_Contig_12327 dDocent_Contig_12412
[115] dDocent_Contig_12416 dDocent_Contig_12422 dDocent_Contig_12448
[118] dDocent_Contig_12450 dDocent_Contig_12456 dDocent_Contig_12461
[121] dDocent_Contig_12467 dDocent_Contig_12475 dDocent_Contig_12480
[124] dDocent_Contig_12481 dDocent_Contig_12485 dDocent_Contig_12534
[127] dDocent_Contig_12543 dDocent_Contig_12585 dDocent_Contig_12679
[130] dDocent_Contig_12687 dDocent_Contig_12697 dDocent_Contig_12722
[133] dDocent_Contig_12746 dDocent_Contig_12752 dDocent_Contig_12772
[136] dDocent_Contig_12832 dDocent_Contig_12929 dDocent_Contig_12943
[139] dDocent_Contig_12956 dDocent_Contig_12981 dDocent_Contig_12998
[142] dDocent_Contig_13022 dDocent_Contig_13047 dDocent_Contig_13151
[145] dDocent_Contig_13166 dDocent_Contig_13170 dDocent_Contig_13184
[148] dDocent_Contig_13315 dDocent_Contig_13333 dDocent_Contig_13406
[151] dDocent_Contig_13586 dDocent_Contig_13687 dDocent_Contig_13834
[154] dDocent_Contig_13839 dDocent_Contig_13925 dDocent_Contig_13936
[157] dDocent_Contig_13963 dDocent_Contig_14088 dDocent_Contig_14263
[160] dDocent_Contig_14355 dDocent_Contig_14365 dDocent_Contig_14366
[163] dDocent_Contig_14395 dDocent_Contig_14443 dDocent_Contig_14457
[166] dDocent_Contig_14462 dDocent_Contig_14505 dDocent_Contig_14514
[169] dDocent_Contig_14540 dDocent_Contig_14553 dDocent_Contig_14721
[172] dDocent_Contig_14792 dDocent_Contig_14798 dDocent_Contig_14879
[175] dDocent_Contig_14889 dDocent_Contig_14890 dDocent_Contig_14897
[178] dDocent_Contig_15203 dDocent_Contig_15307 dDocent_Contig_15430
[181] dDocent_Contig_15491 dDocent_Contig_15567 dDocent_Contig_15633
[184] dDocent_Contig_15669 dDocent_Contig_15754 dDocent_Contig_15872
[187] dDocent_Contig_16084 dDocent_Contig_16116 dDocent_Contig_16198
[190] dDocent_Contig_16472 dDocent_Contig_16570 dDocent_Contig_16714
[193] dDocent_Contig_16744 dDocent_Contig_16753 dDocent_Contig_16794
[196] dDocent_Contig_16970 dDocent_Contig_17097 dDocent_Contig_17107
[199] dDocent_Contig_17297 dDocent_Contig_17514 dDocent_Contig_17536
[202] dDocent_Contig_17553 dDocent_Contig_17638 dDocent_Contig_17854
[205] dDocent_Contig_18538 dDocent_Contig_18542 dDocent_Contig_18581
[208] dDocent_Contig_18729 dDocent_Contig_18746 dDocent_Contig_18753
[211] dDocent_Contig_19160 dDocent_Contig_19243 dDocent_Contig_1955
[214] dDocent_Contig_1959  dDocent_Contig_19619 dDocent_Contig_19715
[217] dDocent_Contig_19736 dDocent_Contig_19854 dDocent_Contig_20
[220] dDocent_Contig_20031 dDocent_Contig_20123 dDocent_Contig_20283
[223] dDocent_Contig_20321 dDocent_Contig_2036  dDocent_Contig_20428
[226] dDocent_Contig_20627 dDocent_Contig_20810 dDocent_Contig_21098
[229] dDocent_Contig_21171 dDocent_Contig_21205 dDocent_Contig_21361
[232] dDocent_Contig_21803 dDocent_Contig_21901 dDocent_Contig_22002
[235] dDocent_Contig_22034 dDocent_Contig_22063 dDocent_Contig_22239
[238] dDocent_Contig_22259 dDocent_Contig_22293 dDocent_Contig_2274
[241] dDocent_Contig_22826 dDocent_Contig_2290  dDocent_Contig_22931
[244] dDocent_Contig_22964 dDocent_Contig_2304  dDocent_Contig_23445
[247] dDocent_Contig_23467 dDocent_Contig_23531 dDocent_Contig_23586
[250] dDocent_Contig_23690 dDocent_Contig_23711 dDocent_Contig_2403
[253] dDocent_Contig_2412  dDocent_Contig_2414  dDocent_Contig_24343
[256] dDocent_Contig_24374 dDocent_Contig_2467  dDocent_Contig_2474
[259] dDocent_Contig_24777 dDocent_Contig_24998 dDocent_Contig_25116
[262] dDocent_Contig_25154 dDocent_Contig_25807 dDocent_Contig_26450
[265] dDocent_Contig_26580 dDocent_Contig_27048 dDocent_Contig_27484
[268] dDocent_Contig_27540 dDocent_Contig_2766  dDocent_Contig_2788
[271] dDocent_Contig_2795  dDocent_Contig_28141 dDocent_Contig_29014
[274] dDocent_Contig_29516 dDocent_Contig_29918 dDocent_Contig_30098
[277] dDocent_Contig_30840 dDocent_Contig_31618 dDocent_Contig_32257
[280] dDocent_Contig_32712 dDocent_Contig_33208 dDocent_Contig_3323
[283] dDocent_Contig_33581 dDocent_Contig_3390  dDocent_Contig_3393
[286] dDocent_Contig_3395  dDocent_Contig_34346 dDocent_Contig_374
[289] dDocent_Contig_39697 dDocent_Contig_41182 dDocent_Contig_41724
[292] dDocent_Contig_4188  dDocent_Contig_4190  dDocent_Contig_4204
[295] dDocent_Contig_4212  dDocent_Contig_4213  dDocent_Contig_4215
[298] dDocent_Contig_4228  dDocent_Contig_4241  dDocent_Contig_4245
[301] dDocent_Contig_4260  dDocent_Contig_4261  dDocent_Contig_4296
[304] dDocent_Contig_4390  dDocent_Contig_43903 dDocent_Contig_4399
[307] dDocent_Contig_4444  dDocent_Contig_4466  dDocent_Contig_4473
[310] dDocent_Contig_4499  dDocent_Contig_4532  dDocent_Contig_4533
[313] dDocent_Contig_4563  dDocent_Contig_4632  dDocent_Contig_4634
[316] dDocent_Contig_4642  dDocent_Contig_4645  dDocent_Contig_4648
[319] dDocent_Contig_4661  dDocent_Contig_4676  dDocent_Contig_4695
[322] dDocent_Contig_48250 dDocent_Contig_48992 dDocent_Contig_4996
[325] dDocent_Contig_4998  dDocent_Contig_5014  dDocent_Contig_5056
[328] dDocent_Contig_5093  dDocent_Contig_5097  dDocent_Contig_5110
[331] dDocent_Contig_5176  dDocent_Contig_5205  dDocent_Contig_5217
[334] dDocent_Contig_52782 dDocent_Contig_5495  dDocent_Contig_5500
[337] dDocent_Contig_5568  dDocent_Contig_5573  dDocent_Contig_5626
[340] dDocent_Contig_5664  dDocent_Contig_57353 dDocent_Contig_5890
[343] dDocent_Contig_62464 dDocent_Contig_6257  dDocent_Contig_6263
[346] dDocent_Contig_6266  dDocent_Contig_6283  dDocent_Contig_6287
[349] dDocent_Contig_6330  dDocent_Contig_6395  dDocent_Contig_6489
[352] dDocent_Contig_6490  dDocent_Contig_6492  dDocent_Contig_6499
[355] dDocent_Contig_6534  dDocent_Contig_6535  dDocent_Contig_65424
[358] dDocent_Contig_6554  dDocent_Contig_6555  dDocent_Contig_6573
[361] dDocent_Contig_6574  dDocent_Contig_6594  dDocent_Contig_6615
[364] dDocent_Contig_6855  dDocent_Contig_6993  dDocent_Contig_7051
[367] dDocent_Contig_7085  dDocent_Contig_7088  dDocent_Contig_7120
[370] dDocent_Contig_7123  dDocent_Contig_7139  dDocent_Contig_715
[373] dDocent_Contig_7219  dDocent_Contig_7220  dDocent_Contig_7223
[376] dDocent_Contig_7224  dDocent_Contig_7226  dDocent_Contig_7232
[379] dDocent_Contig_7245  dDocent_Contig_7247  dDocent_Contig_7251
[382] dDocent_Contig_7258  dDocent_Contig_7276  dDocent_Contig_7331
[385] dDocent_Contig_7368  dDocent_Contig_7446  dDocent_Contig_7454
[388] dDocent_Contig_7466  dDocent_Contig_7494  dDocent_Contig_7607
[391] dDocent_Contig_7620  dDocent_Contig_7688  dDocent_Contig_7692
[394] dDocent_Contig_7705  dDocent_Contig_7769  dDocent_Contig_7775
[397] dDocent_Contig_7776  dDocent_Contig_7786  dDocent_Contig_7787
[400] dDocent_Contig_7792  dDocent_Contig_7801  dDocent_Contig_7809
[403] dDocent_Contig_7811  dDocent_Contig_7835  dDocent_Contig_7840
[406] dDocent_Contig_7857  dDocent_Contig_7861  dDocent_Contig_7874
[409] dDocent_Contig_7920  dDocent_Contig_7935  dDocent_Contig_7939
[412] dDocent_Contig_7994  dDocent_Contig_8000  dDocent_Contig_80103
[415] dDocent_Contig_8021  dDocent_Contig_8058  dDocent_Contig_8060
[418] dDocent_Contig_8140  dDocent_Contig_8176  dDocent_Contig_8218
[421] dDocent_Contig_8223  dDocent_Contig_8226  dDocent_Contig_8248
[424] dDocent_Contig_8253  dDocent_Contig_8255  dDocent_Contig_8272
[427] dDocent_Contig_8299  dDocent_Contig_8300  dDocent_Contig_8302
[430] dDocent_Contig_8310  dDocent_Contig_8314  dDocent_Contig_8318
[433] dDocent_Contig_8343  dDocent_Contig_8344  dDocent_Contig_8354
[436] dDocent_Contig_8362  dDocent_Contig_8384  dDocent_Contig_8408
[439] dDocent_Contig_8415  dDocent_Contig_8434  dDocent_Contig_8440
[442] dDocent_Contig_8446  dDocent_Contig_8447  dDocent_Contig_8448
[445] dDocent_Contig_8458  dDocent_Contig_8459  dDocent_Contig_8468
[448] dDocent_Contig_8470  dDocent_Contig_8484  dDocent_Contig_8494
[451] dDocent_Contig_8540  dDocent_Contig_8548  dDocent_Contig_8553
[454] dDocent_Contig_8559  dDocent_Contig_8570  dDocent_Contig_8606
[457] dDocent_Contig_8611  dDocent_Contig_8613  dDocent_Contig_8651
[460] dDocent_Contig_8653  dDocent_Contig_8668  dDocent_Contig_8684
[463] dDocent_Contig_8688  dDocent_Contig_8701  dDocent_Contig_8718
[466] dDocent_Contig_8745  dDocent_Contig_8752  dDocent_Contig_8767
[469] dDocent_Contig_8768  dDocent_Contig_8771  dDocent_Contig_8773
[472] dDocent_Contig_8794  dDocent_Contig_8810  dDocent_Contig_8818
[475] dDocent_Contig_8819  dDocent_Contig_8833  dDocent_Contig_8927
[478] dDocent_Contig_8930  dDocent_Contig_8949  dDocent_Contig_8951
[481] dDocent_Contig_8961  dDocent_Contig_8967  dDocent_Contig_9096
[484] dDocent_Contig_9107  dDocent_Contig_9110  dDocent_Contig_9172
[487] dDocent_Contig_9193  dDocent_Contig_9216  dDocent_Contig_9228
[490] dDocent_Contig_9237  dDocent_Contig_9267  dDocent_Contig_9286
[493] dDocent_Contig_9316  dDocent_Contig_9327  dDocent_Contig_9353
[496] dDocent_Contig_9363  dDocent_Contig_9374  dDocent_Contig_9377
[499] dDocent_Contig_9379  dDocent_Contig_9397  dDocent_Contig_9416
[502] dDocent_Contig_9431  dDocent_Contig_9447  dDocent_Contig_9497
[505] dDocent_Contig_9501  dDocent_Contig_9510  dDocent_Contig_9513
[508] dDocent_Contig_9515  dDocent_Contig_9524  dDocent_Contig_9553
[511] dDocent_Contig_9600  dDocent_Contig_9602  dDocent_Contig_9612
[514] dDocent_Contig_9626  dDocent_Contig_9658  dDocent_Contig_9660
[517] dDocent_Contig_9687  dDocent_Contig_9688  dDocent_Contig_9699
[520] dDocent_Contig_9723  dDocent_Contig_9727  dDocent_Contig_9753
[523] dDocent_Contig_9757  dDocent_Contig_9781  dDocent_Contig_9851
[526] dDocent_Contig_9863  dDocent_Contig_9887  dDocent_Contig_9966
[529] dDocent_Contig_9995
}
} #notepad cleanup

#Removing monomorphic SNPs
{```{R}```
#Load Libraries
library('vcfR')
library('adegenet')

#Loading the vcf file
vcf <- read.vcfR(file="SNP.TRS.F05.vcf")
gen.vcf<-vcfR2genind(vcf)
length(which(gen.vcf@loc.n.all == 1))
#[1] 2528
rm.loci <- names(gen.vcf@loc.n.all)[gen.vcf@loc.n.all == 1]
loc.m <- matrix(unlist(strsplit(rm.loci,"_")),ncol=5,byrow=T)
loc.df <- data.frame(CHROM=apply(loc.m[,1:3],1,function(x) paste(x, collapse="_")), POS=loc.m[,4])
write.table(loc.df, file="mono_loci.txt", col.names=F, row.names=F, quote=F, sep="\t")
} #notepad cleanup
{```{bash}```
vcftools --vcf SNP.TRS.F05.vcf --out SNP.TRS.F05a --recode --recode-INFO-all --exclude-positions mono_loci.txt
#After filtering, kept 65 out of 65 Individuals
#After filtering, kept 32783 out of a possible 35311 Sites
#Run Time = 13.00 seconds

charts.sh SNP.TRS.F05a
} #notepad cleanup

#Removing paralogs and indivduals with high missing data
{```{bash}```
mkdir haplotyping
cd haplotyping
cp -s ../../mapping/reference.fasta .
cp -s ../popmap .

ls *.bam | while read i; do echo "Processing $i"; samtools index -@ 20 $i; done
cp -s ../SNP.TRS.F05a.recode.vcf .

rad_haplotyper.pl -v SNP.TRS.F05a.recode.vcf -r reference.fasta -p popmap -x 20 -m 0.8 -o SNP.TRS.F06.vcf -g SNP.TRS.F06.gen -a SNP.TRS.F06.ima

cd ..
cp -s haplotyping/SNP.TRS.F06.vcf .
charts.sh SNP.TRS.F06.vcf
} #notepad cleanup

#Removing high Heterozygosity and high depth per locus
{```{bash}```
cd SNP.TRS.F06.2023-09-25
} #notepad cleanup
{```{R}```
indv <- read.table("indv.csv", head=T, sep=",")
loci <- read.table("loci.csv", head=T, sep=",")

loci_high <- loci[loci$PER_HET > 0.5, ]
loci_low <- loci[loci$PER_HET <= 0.5, ]

loci_high.contig <- unique(loci_high$CHROM)
loci_low.contig <- unique(loci_low$CHROM)

length(loci_high.contig)
#[1] 528
length(loci_low.contig)
#[1] 8494

length(which(loci_high.contig %in% loci_low.contig))
#[1] 446

nrow(loci[loci$CHROM %in% loci_high.contig,])
#[1] 2147

write.table(loci[loci$CHROM %in% loci_high.contig, c("CHROM","POS")], file="high_het.txt", sep="\t", col.names=F, row.names=F, quote=F)
write.table(loci[loci$MEAN_DEPTH > 83.5, c("CHROM","POS")], file="high_depth.txt", sep="\t", col.names=F, row.names=F, quote=F)
}
{```{bash}```
cd ..
vcftools --vcf SNP.TRS.F06.vcf --out SNP.TRS.F06a --recode --recode-INFO-all --exclude-positions SNP.TRS.F06.2023-09-25/high_het.txt
#After filtering, kept 65 out of 65 Individuals
#After filtering, kept 25127 out of a possible 27274 Sites
#Run Time = 10.00 seconds
vcftools --vcf SNP.TRS.F06a.recode.vcf --out SNP.TRS.F07 --recode --recode-INFO-all --exclude-positions SNP.TRS.F06.2023-09-25/high_depth.txt
#After filtering, kept 65 out of 65 Individuals
#After filtering, kept 25006 out of a possible 25127 Sites
#Run Time = 9.00 seconds

charts.sh SNP.TRS.F07.recode.vcf
} #notepad cleanup

#Exploring the increased missing data per locus
{```{R}```
library('vcfR')
library('adegenet')

#Loading the vcf file
vcf <- read.vcfR(file="SNP.TRS.F05.vcf")
gen.vcf <- vcfR2genind(vcf)
apply(gen.vcf@tab[1:4, grep("dDocent_Contig_10017",colnames(gen.vcf@tab))], 2, function(x) sum(is.na(x)))
dDocent_Contig_10017_311_5253.0 dDocent_Contig_10017_311_5253.1
                              0                               0
dDocent_Contig_10017_316_5254.0 dDocent_Contig_10017_316_5254.1
                              0                               0
dDocent_Contig_10017_339_5255.0 dDocent_Contig_10017_339_5255.1
                              0                               0
dDocent_Contig_10017_381_5256.0 dDocent_Contig_10017_381_5256.1
                              0                               0

locus.m <- matrix(unlist(strsplit(locNames(gen.vcf),"_")),ncol=5,byrow=T)
locus.df <- data.frame(CHROM=apply(locus.m[,1:3], 1, function(x) paste(x,collapse="_")), POS=locus.m[,4], ORDER=locus.m[,5])
locus.tab <- table(locus.df$CHROM)
test.list <- names(locus.tab)[locus.tab > 1]

weird.df <- data.frame(matrix(ncol=3))
colnames(weird.df) <- c("CHROM1", "CHROM2", "Mismatch")

for(i in test.list){
	if(which(test.list == i)/100 == round(which(test.list == i)/100,0)){print(paste("processing", which(test.list == i), "of", length(test.list)))}
	POS <- locNames(gen.vcf)[grep(i, locNames(gen.vcf))]
	tmp.list <- list()
	for(j in POS){tmp.list[[j]] <- which(is.na(gen.vcf@tab[grep(i,colnames(gen.vcf@tab))[1]]))}
	tmp.m <- t(combn(names(tmp.list), 2))
	tmp.df <- data.frame(CHROM1=tmp.m[,1], CHROM2=tmp.m[,2], Mismatch=rep("NA", nrow(tmp.m)))
	tmp.df$Mismatch <- as.numeric(tmp.df$Mismatch)
	for(j in 1:nrow(tmp.df)){tmp.df$Mismatch[j] <- length(which(!tmp.list[[as.character(as.matrix(tmp.df[j,1]))]] %in% tmp.list[[as.character(as.matrix(tmp.df[j,2]))]]))}
	if(sum(tmp.df$Mismatch) > 0){weird.df <- rbind(weird.df, tmp.df)}
}
nrow(weird.df)
#1	#Only the NA row from creating the matrix
}}} #notepad cleanup