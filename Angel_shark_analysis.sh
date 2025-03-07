#Loading Libraries
{```{R}```
library('devtools')
library('pegas')
library('adegenet')
library('vcfR') 
library('rospca')
library('dartR')
library('zvau')
library('geosphere')
library('stringr')
library('ggmap')
library('ggcompoplot')
library('vegan')
library('spdep')
library('adespatial')
library('igraph')
library('poppr') 
library('smatr')
library('radiator')
library('related')
library('ggcompoplot')
library('scales')
library('hierfstat')
library('VennDiagram')
library('sdmpredictors') 
library('ggord')
} #notepad cleanup

#Add functions
{```{R}```
#Plot multiple ggplot objects in the same window or save file
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
} #http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/

}}

#Function from Forester to get alleles for selecting outliers along RDA axes
#https://popgen.nescent.org/2018-03-27_RDA_GEA.html
outliers <- function(x,z){
   lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading
   x[x < lims[1] | x > lims[2]]               # locus names in these tails
 }

#Function to split out loci which are driving the difference between two PCA groups
Locus_hunt <- function(GEN, DAPC, THRES=0.8, GROUP="~POP", AXIS=1, ALPHA=0.05, MIN=0.5, ABS=0.001, MODEL=c("normal", "outlier")){
#Objects to put into the function

#GEN <- The genind object used to make a DAPC
#DAPC <- The DAPC oject made from the genind object
#THRES <- The precent of loci to remove during the 1st round
#GROUP <- The strata group that should be used
#AXIS <- The axis you would like to search along
#ALPHA <- The significance value you would like to work with 
#MIN <- The lower limit that the threshold can approach, but never reach
#ABS <- The difference between the thresholds which determines if the script continues to refine the loci

#Fills in missing values
if(missing(THRES)){THRES=0.8}
if(missing(GROUP)){GROUP="~POP"}
if(missing(AXIS)){AXIS=1}
if(missing(ALPHA)){ALPHA=0.05}
if(missing(MIN)){MIN=0.5}
if(missing(MODEL)){"normal"}

# Sets up the varibles tested
if(THRES<=MIN){THRES <- MIN + 0.05}
if(THRES>=0.9975){THRES <- 0.95}

# Loop to look for optimal value
if(MODEL=="normal"){
GOOD <- as.vector(MIN)
BAD <- as.vector(1)
while(abs(THRES-max(GOOD)) > ABS){
if(length(GOOD)>20){cat(paste("Threshold not found after 20 iterations\nPlease try a new starting point\nBest Threshold tried was",max(GOOD),"\n")); break}
print(paste("Processing threshold", THRES))
# Getting the loci to test
contrib <- rownames(DAPC$var.contr)[which(DAPC$var.contr[,AXIS] > quantile(DAPC$var.contr[,AXIS], prob=THRES))]
pc.loci <- unique(matrix(unlist(strsplit(contrib, "[.]")),ncol=2,byrow=T)[,1])

# Select out the data sets
set.loc<-locNames(GEN)[which(!(locNames(GEN) %in% pc.loci))]
tmp.keep<-GEN[, loc=set.loc]
set.loc<-locNames(GEN)[which((locNames(GEN) %in% pc.loci))]
tmp.rm<-GEN[, loc=set.loc]

# Attribute Variance
tmp.keep.result <- poppr.amova(tmp.keep, as.formula(GROUP), quiet=T)
tmp.rm.result <- poppr.amova(tmp.rm, as.formula(GROUP), quiet=T)

# Test for significance
tmp.keep.test <- randtest(tmp.keep.result)
tmp.rm.test <- randtest(tmp.rm.result)

#Tables for comparison
tmp.ktable <- data.frame(cbind(tmp.keep.test$names, tmp.keep.test$pvalue),row.names=1)
names(tmp.ktable) <- "Kept"
tmp.rtable <- data.frame(cbind(tmp.rm.test$names, tmp.rm.test$pvalue),row.names=1)
names(tmp.rtable) <- "Removed"
print(cbind(tmp.ktable, tmp.rtable))

#Decision of what to do next
if(as.numeric(as.matrix(tmp.ktable[nrow(tmp.ktable),1])) < ALPHA |	as.numeric(as.matrix(tmp.rtable[nrow(tmp.rtable),1])) > ALPHA){
	BAD <- c(BAD, THRES)
} else if(as.numeric(as.matrix(tmp.ktable[nrow(tmp.ktable),1])) > ALPHA | as.numeric(as.matrix(tmp.rtable[nrow(tmp.rtable),1])) < ALPHA){
	GOOD <- c(GOOD, THRES)
	STORE.GOOD.gen.keep <- tmp.keep
	STORE.GOOD.gen.rm <- tmp.rm
	STORE.GOOD.tab.keep <- tmp.ktable
	STORE.GOOD.tab.rm <- tmp.rtable
}
THRES <- mean(c(max(GOOD),min(BAD)))
}
THRES <- max(GOOD)
} else if(MODEL=="outlier"){
GOOD <- as.vector(1)
BAD <- vector()
while(abs(THRES-min(GOOD)) > ABS){
if(length(GOOD)>20){cat(paste("Threshold not found after 20 iterations\nPlease try a new starting point\nBest Threshold tried was",max(GOOD),"\n")); break}
print(paste("Processing threshold", THRES))
# Getting the loci to test
contrib <- rownames(DAPC$var.contr)[which(DAPC$var.contr[,AXIS] > quantile(DAPC$var.contr[,AXIS], prob=THRES))]
pc.loci <- unique(matrix(unlist(strsplit(contrib, "[.]")),ncol=2,byrow=T)[,1])

# Select out the data sets
set.loc<-locNames(GEN)[which(!(locNames(GEN) %in% pc.loci))]
tmp.keep<-GEN[, loc=set.loc]
set.loc<-locNames(GEN)[which((locNames(GEN) %in% pc.loci))]
tmp.rm<-GEN[, loc=set.loc]

# Attribute Variance
tmp.keep.result <- poppr.amova(tmp.keep, as.formula(GROUP), quiet=T)
tmp.rm.result <- poppr.amova(tmp.rm, as.formula(GROUP), quiet=T)

# Test for significance
tmp.keep.test <- randtest(tmp.keep.result)
tmp.rm.test <- randtest(tmp.rm.result)

#Tables for comparison
tmp.ktable <- data.frame(cbind(tmp.keep.test$names, tmp.keep.test$pvalue),row.names=1)
names(tmp.ktable) <- "Kept"
tmp.rtable <- data.frame(cbind(tmp.rm.test$names, tmp.rm.test$pvalue),row.names=1)
names(tmp.rtable) <- "Removed"
print(cbind(tmp.ktable, tmp.rtable))

#Decision of what to do next
if(as.numeric(as.matrix(tmp.rtable[nrow(tmp.rtable),1])) > ALPHA){
	BAD <- c(BAD, THRES)
	THRES <- mean(c(min(GOOD),max(BAD)))
} else if(as.numeric(as.matrix(tmp.rtable[nrow(tmp.rtable),1])) < ALPHA){
	GOOD <- c(GOOD, THRES)
	STORE.GOOD.gen.keep <- tmp.keep
	STORE.GOOD.gen.rm <- tmp.rm
	STORE.GOOD.tab.keep <- tmp.ktable
	STORE.GOOD.tab.rm <- tmp.rtable
	THRES <- mean(c(min(GOOD),MIN))
}
}
THRES <- min(GOOD)
}

#Outputs
OUTPUT.list <- list()
OUTPUT.list[[1]] <- STORE.GOOD.gen.keep
OUTPUT.list[[2]] <- STORE.GOOD.gen.rm
OUTPUT.list[[3]] <- data.frame(Summary = c(paste("The Best Threshold tested was ", THRES, sep=""), 
paste("This removed ", length(locNames(STORE.GOOD.gen.rm)), " loci (", round(length(locNames(STORE.GOOD.gen.rm))/length(locNames(GEN)),3)*100,"%)",sep=""),
paste("Kept p-value is ", STORE.GOOD.tab.keep[nrow(STORE.GOOD.tab.keep),1],sep=""),
paste("Removed p-value is ", STORE.GOOD.tab.rm[nrow(STORE.GOOD.tab.rm),1], sep=""),
paste("Data has been output to the 1st two objects of this list"),
paste("Item 1 contains the kept genind object"),
paste("Item 2 contains the removed genind object")))
OUTPUT.list[[3]] <- format(OUTPUT.list[[3]], justify = "left")
return(OUTPUT.list)
}

}}}}}}}}}}}}}}}}#Cleanup

#Pulls the unique loci names from alleles selected from genind object names (e.g. dDocent_Contig_23960.001)
Loci_names <- function(NAMES, SEP="[.]", REMOVE=1){
COL <- length(strsplit(head(NAMES,n=1), SEP)[[1]])
TMP_U <- unique(matrix(unlist(strsplit(NAMES,SEP)),ncol=COL,byrow=T)[,1:(COL-REMOVE)])
if(is.matrix(TMP_U)){TMP_DF <- data.frame(TMP_U)
} else {TMP_DF <- data.frame(matrix(TMP_U, ncol=(COL-REMOVE), byrow=T))}
return(tidyr::unite(TMP_DF, "loci", 1:ncol(TMP_DF), sep=SEP))
}

} #notepad cleanup

#Removing loci for genind object from a vcf; VCF locNames need to have a "dDocent_Contig_#_Pos" format
VCF_remove <- function(genind.vcf, loci.list){
tmp.loci <- matrix(unlist(strsplit(locNames(genind.vcf),"_")), ncol=4, byrow=T)
if(sum(is.na(tail(tmp.loci, n=1)))){print("locNames(genind) in incorrect format"); break}
tmp.loci <- data.frame(Locus=paste(tmp.loci[,1], tmp.loci[,2], tmp.loci[,3],sep="_"), Pos=tmp.loci[,4])

tmp <- vector()
for(i in loci.list){
ROW <- which(tmp.loci$Locus == i)
tmp <- append(tmp,paste(tmp.loci$Locus, tmp.loci$Pos, sep="_")[ROW])
}

set.loc <- subset(locNames(genind.vcf), !locNames(genind.vcf) %in% tmp)
return(genind.vcf[ ,loc=set.loc])
}

} #Notepad cleanup

#updating ordiR2Stap
ordiR2step <- function (object, scope, Pin = 0.05, R2scope = TRUE, permutations = how(nperm = 499),
    trace = TRUE, R2permutations = 1000, VIF = 3, ...)
{
    if (is.null(object$terms))
        stop("ordination model must be fitted using formula")
    if (missing(scope))
        stop("needs scope")
    if (inherits(scope, "cca"))
        scope <- delete.response(formula(scope))
    if (!inherits(scope, "formula"))
        scope <- reformulate(scope)
    if (is.null(object$CCA))
        R2.0 <- 0
    else R2.0 <- RsquareAdj(object, permutations = R2permutations,
        ...)$adj.r.squared
    if (is.list(scope) && length(scope) <= 2L)
        scope <- scope$upper
    if (is.null(scope) || !length(add.scope(object, scope)))
        stop("needs upper 'scope': no terms can be added")
    if (R2scope)
        R2.all <- RsquareAdj(update(object, delete.response(formula(scope))), permutations = R2permutations, ...)
    else R2.all <- list(adj.r.squared = NA)
    if (is.na(R2.all$adj.r.squared) && R2scope)
        stop("the upper scope cannot be fitted (too many terms?)")
    R2.all <- R2.all$adj.r.squared
    anotab <- list()
    R2.previous <- R2.0
	high.vif <- matrix(ncol=3)
    repeat {
        if (trace) {
            cat("Step: R2.adj=", R2.previous, "\n")
            cat(pasteCall(formula(object)), "\n")
        }
		#print("object:")
		#print(as.character(object$call)[2])
		#print("scope")
		#print(scope)
        adds <- add.scope(object, scope)
        if (length(adds) == 0)
            break
        R2.adds <- numeric(length(adds))
        adds <- paste("+", adds)
        names(R2.adds) <- adds
        for (trm in seq_along(R2.adds)) {
            fla <- paste(". ~ .", names(R2.adds[trm]))
            R2.tmp <- RsquareAdj(update(object, fla), permutations = R2permutations)$adj.r.squared
            if (!length(R2.tmp))
                R2.tmp <- 0
            R2.adds[trm] <- R2.tmp
        }
        best <- which.max(R2.adds)
        if (trace) {
            out <- sort(c(`<All variables>` = R2.all, `<none>` = R2.previous,
                R2.adds), decreasing = TRUE)
            out <- as.matrix(out)
            colnames(out) <- "R2.adjusted"
            print(out)
            cat("\n")
        }
        if (R2.adds[best] > R2.previous && (!R2scope || R2scope &&
            R2.adds[best] <= R2.all)) {
            adds <- names(sort(R2.adds, decreasing=T))
            vif.test <- VIF+1
            while(max(vif.test, na.rm=T) > VIF){
            object_tmp <- object
            fla <- paste("~  .", adds[1])
            object_tmp <- update(object_tmp, fla)
            vif.test <- vif.cca(object_tmp)
            if (max(vif.test, na.rm=T) > VIF) {
                if (trace) {
                    print(paste("Max VIF =", max(vif.test, na.rm=T)))
                    print(paste(str_sub(adds[1], start=3), "removed from consideration"))
                    }
                high.vif<- rbind(high.vif, c(str_sub(adds[1], start=3),max(vif.test, na.rm=T),as.character(object$call)[2]))
                scope <- terms.formula(scope)
                if(length(attr(scope,"term.labels")) == 2){break}
                scope <- reformulate(attr(scope, "term.labels")[!attr(scope, "term.labels") %in% str_sub(adds[1], start=3)])
                adds <- adds[-1]
				if(length(adds) == 0){break}
				#print(scope)
                }
            }
			if(length(adds) == 0){break}
            tst <- add1(object, scope = adds[1], test = "permu", permutations = permutations, alpha = Pin, trace = FALSE)
            if (trace) {
                print(tst[-1, ])
                cat("\n")
            }
			if (is.na(tst[, "Pr(>F)"][2])){
				print(paste(adds[1], "has no constrained component"))
                scope <- terms.formula(scope)
                if(length(attr(scope,"term.labels")) == 2){break}
                scope <- reformulate(attr(scope, "term.labels")[!attr(scope, "term.labels") %in% str_sub(adds[1], start=3)])
                adds <- adds[-1]
			} else if (tst[, "Pr(>F)"][2] <= Pin) {
                fla <- paste("~  .", adds[1])
                object <- update(object, fla)
            } else break
        }
        else {
            break
        }
		R2.previous <- RsquareAdj(object, permutations = R2permutations)$adj.r.squared
        anotab <- rbind(anotab, cbind(R2.adj = R2.previous, tst[2, ]))
        }
    if (NROW(anotab)) {
        if (R2scope)
            anotab <- rbind(anotab, `<All variables>` = c(R2.all, rep(NA, 4)))
        class(anotab) <- c("anova", class(anotab))
        object$anova <- anotab
    }
	if(nrow(high.vif)>0){colnames(high.vif) <- c("Rejected_variable","Max_VIF","Formula")}
    attr(object,'VIF_remove_list') <- na.omit(high.vif)
    object
}

}}}}}}}}}}}}}}}}}}}}}}}}}} #notepad cleanup
} #notepad cleanup

#Importing Data
{```{R}```
strata <- read.csv(file = "../../filtering_Sep2023/SNP.TRS.F07.2023-09-26/indv.csv", header = TRUE)
head(strata)

vcf<-read.vcfR(file="SNP.TRS.F07.recode.vcf")
gen.vcf<-vcfR2genind(vcf)
strata(gen.vcf) <- strata[match(indNames(gen.vcf),strata$INDV),]
head(strata(gen.vcf))
rm(vcf)

set.loc <- Loci_names(locNames(gen.vcf), SEP="_", REMOVE=1)
gen <- read.genepop(file = "SNP.TRS.F06.gen", ncode=3L, quiet = FALSE)[, loc=set.loc[,1], drop=T]
strata(gen) <- strata[match(indNames(gen),strata$INDV),]
head(strata(gen))

gen@other$lat<-gen@strata[,c("Lon","Lat")]
gen@other$lat[,1]<-as.numeric(as.character(gen@strata$Lon))
gen@other$lat[,2]<-as.numeric(as.character(gen@strata$Lat))	#As numeric
gen@other$xy<-gen@strata[,c("Lon","Lat")]	#As factor

#Adding related ID to strata
gen.indv <- substr(gen@strata$INDV, 1, 20)
gen.indv <- data.frame(cbind(as.character(gen@strata$INDV), gen.indv))
names(gen.indv) <- c("gen.ID", "relate.ID")

gen@strata$relate.ID <- gen.indv[match(indNames(gen),gen.indv$gen.ID),"relate.ID"]
head(strata(gen))

save(gen, file="gen.gz", compress=T)
#load("gen.gz")
} #notepad cleanup

#Setting color schemes
{```{R}```
#           "Atl"   "EAST"    "Pac"     "WEST"
col.All<-c("red4","dodgerblue","grey50","mediumblue")

c2 <- c("mediumblue", "red4")
c3 <- c("dodgerblue", "red4", "mediumblue")
c4 <- c("lightslateblue", "mediumblue", "red4", "dodgerblue")
c5 <- c("dodgerblue", "mediumblue", "royalblue1", "lightslateblue", "red4")
c6 <- c("mediumblue", "lightslateblue", "royalblue1", "red4","cornflowerblue", "dodgerblue")
} #notepad cleanup

#Removing duplates
{```{R}```
remove.ind <- c("Atl.AS.124_Lib1_I12_D05", "EAST.AS.9.2_Lib1_I12_D04")
set.ind <- subset(indNames(gen), !indNames(gen) %in% remove.ind)

gen2 <- gen[set.ind, , drop=T]
gen2.vcf <- gen.vcf[set.ind, , drop=T]

save(gen2, file="gen2_wPac.gz", compress=T)
save(gen2.vcf, file="gen2.vcf_wPac.gz", compress=T)
} #notepad cleanup

#Looking for rough relatedness between individuals
{```{R}```
gl <- gi2gl(gen2)
test.cov <- gl2related(gl, save=F)
dim(test.cov)
#[1]    63 16021
test.out <- coancestry(test.cov, wang =1)
save(test.out, file="F07.related_wang.gz", compress=T)
#load("F07.related_wang.gz")
} #notepad cleanup

#Removing duplicates based upon relatedness
#Removing 1 of those which have a relationship value larger than 1
{```{R}```
dups.rel <- test.out$relatedness[test.out$relatedness$dyadml > 1,2:3]
dups.rel
} #notepad cleanup

#PCA with Pacific sample
{```{R}```
#Standard
X <- scaleGen(gen2, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen2)<-~POP
par(mfrow=c(1,1))
ade4::s.class(pca1$li, pop(gen2), col=col.All, cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
legend(x=-25.5, y=-22,legend=c("Atlantic (n=17)", "East GOM (n=25)", "West GOM (n=20)", "Pacific (n=1)"),col=c("red4","dodgerblue", "mediumblue","grey50"), cex=1.5, pch=16, bty="n", ncol=1)

tiff("PCA_POP_all_loci_wPac.tif", res=300, height =2000, width=2000)
ade4::s.class(pca1$li, pop(gen2), col=col.All, cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
legend(x=-25.5, y=-22,legend=c("Atlantic (n=17)", "East GOM (n=25)", "West GOM (n=20)", "Pacific (n=1)"),col=c("red4","dodgerblue", "mediumblue","grey50"), cex=1.5, pch=16, bty="n", ncol=1)
dev.off()
} #notepad cleanup

#Removing Pacific sample
{```{R}```
remove.ind <- c("Atl.AS.124_Lib1_I12_D05", "EAST.AS.9.2_Lib1_I12_D04", "Pac.AS.T437_Lib2_I12_D01")
set.ind <- subset(indNames(gen), !indNames(gen) %in% remove.ind)

gen2 <- gen[set.ind, , drop=T]
gen2.vcf <- gen.vcf[set.ind, , drop=T]

save(gen2, file="gen2_noPac.gz", compress=T)
save(gen2.vcf, file="gen2.vcf_noPac.gz", compress=T)
} #notepad cleanup

#Library Bias (Outflank)
{```{R}```
setPop(gen2.vcf)<-~Lib
tmp<-gl.outflank(gen2.vcf, qthreshold = 0.1)	#Can use plot=F option if you want
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))
#283

outliers <- tmp$outflank$results[which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"),1]
out.list <- matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1]
tmp <- matrix(unlist(strsplit(as.vector(out.list),"[_]")),ncol=4,byrow=T)
#head(tmp)
out.list <- unique(apply(tmp[,1:3], 1, function(x) paste(x, collapse="_")))
length(out.list) 
#90
gen.out.list<-(out.list)

set.loc<-locNames(gen2)[which(!(locNames(gen2) %in% gen.out.list))]
gen.net<-gen2[, loc=set.loc]
set.loc<-locNames(gen2)[which((locNames(gen2) %in% gen.out.list))]
gen.out<-gen2[, loc=set.loc]

X <- scaleGen(gen2, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.out <- scaleGen(gen.out, NA.method="mean", scale=F)
pca.out <- dudi.pca(X.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen2) <- setPop(gen.net) <- setPop(gen.out) <- ~POP
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen2), col=col.All, cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen2))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen.net), col=col.All, cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen.out), col=col.All, cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.out))," loci)", sep=""), 3, 2, adj = 0.95)

setPop(gen2) <- setPop(gen.net) <- setPop(gen.out) <- ~Lib
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen2), col=col.All, cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen2))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen.net), col=col.All, cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen.out), col=col.All, cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.out))," loci)", sep=""), 3, 2, adj = 0.95)

write.table(gen.out.list, "Outflank_Library_Outliers.list", quote=F, col.names=F, row.names=F)
gen.out.list <- read.table("Outflank_Library_Outliers.list", head=F)[,1]
gen.out.list <- as.character(as.matrix(gen.out.list))
gen3 <- gen.net
gen3.vcf <- VCF_remove(gen2.vcf, gen.out.list)

save(gen3, file="gen3_noPac.gz", compress=T)
save(gen3.vcf, file="gen3.vcf_noPac.gz", compress=T)
} #notepad cleanup

#Library Bias (Bayescan)
{```{R}```
#Only 2 groups present, therefore cannot run Bayescan
} #notepad cleanup

#Finding monomorphic loci
{```{R}```
mono.list <- names(which(gen2@loc.n.all==1))
length(mono.list)
#[1] 21
} #notepad cleanup

#Looking for those with extreme numbers of alleles
#Experimental
{```{R}```
tmp_dat <- apply(gen2@tab, 2, function(x) sum(x,na.rm=T))
tmp_dat <- tmp_dat[tmp_dat > 0]
tmp_names <- matrix(unlist(strsplit(names(tmp_dat),"[.]")),ncol=2, byrow=T)[,1]
tmp_tab <- table(tmp_names)
par(mfrow=c(2,1))
hist(tmp_tab, breaks=100, col="red4")
hist(tmp_tab, breaks=100, col="red4", ylim=c(0,100))
Increased_SNPs <- names(boxplot(tmp_tab,plot=F, range=3)$out)
table(gen2@loc.n.all[Increased_SNPs])
{ #Results
15 16 17 18 19 20 22 23 24 35 44
11  2  4  2  2  2  1  1  1  1  1
}
length(Increased_SNPs)
#28

pois2.v <- dpois(0:45,2)
SCALE <- 2000/max(pois2.v)
hist(tmp_tab, breaks=100, col="red4")
lines(pois2.v*SCALE)

pois3.v <- dpois(0:45,3)
SCALE <- 2000/max(pois3.v)
lines(pois3.v*SCALE, col="mediumblue")

pois25.v <- dpois(0:45,2.5)
SCALE <- 2000/max(pois25.v)
lines(pois25.v*SCALE, col="darkgreen")

pois25.v*SCALE
{ #Results
 [1] 6.400000e+02 1.600000e+03 2.000000e+03 1.666667e+03 1.041667e+03
 [6] 5.208333e+02 2.170139e+02 7.750496e+01 2.422030e+01 6.727861e+00
[11] 1.681965e+00 3.822648e-01 7.963851e-02 1.531510e-02 2.734839e-03
[16] 4.558065e-04 7.121976e-05 1.047349e-05 1.454652e-06 1.914016e-07
[21] 2.392520e-08 2.848238e-09 3.236634e-10 3.518080e-11 3.664667e-12
[26] 3.664667e-13 3.523718e-14 3.262702e-15 2.913127e-16 2.511316e-17
[31] 2.092764e-18 1.687713e-19 1.318525e-20 9.988829e-22 7.344727e-23
[36] 5.246234e-24 3.643218e-25 2.461634e-26 1.619496e-27 1.038138e-28
[41] 6.488365e-30 3.956320e-31 2.354952e-32 1.369158e-33 7.779309e-35
[46] 4.321838e-36
}

length(which(gen2@loc.n.all > 11))
#72
gen2@loc.n.all[gen2@loc.n.all == 16]
#dDocent_Contig_65662 dDocent_Contig_23976 dDocent_Contig_32935 dDocent_Contig_30954 dDocent_Contig_36070
#                  16                   16                   16                   16                   16
} #notepad cleanup

#Removing Library outliers, monomorphic loci, high SNP count
{```{R}```
OF.out <- read.table("Outflank_Library_Outliers.list", head=F)
gen.out.list <- unique(c(as.character(as.matrix(OF.out)),mono.list, Increased_SNPs))
set.loc <- locNames(gen2)[which(!(locNames(gen2) %in% gen.out.list))]

gen3 <- gen2[, loc=set.loc]
gen3.vcf <- VCF_remove(gen2.vcf, gen.out.list)

save(gen3, file="gen3_noPac.gz", compress=T)
save(gen3.vcf, file="gen3.vcf_noPac.gz", compress=T)
#load("gen3_noPac.gz")
#load("gen3.vcf_noPac.gz")
} #notepad cleanup

#PCA of filtered data
{```{R}```
X <- scaleGen(gen3, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen3) <- ~G3

col_pts <- as.character(as.matrix(gen3@strata$G3))
col_pts[col_pts == "Atlantic"] <- "chocolate2"
col_pts[col_pts == "East"] <- "green"
col_pts[col_pts == "West"] <- "mediumblue"

#Axes 1 & 2
png("PCA_POP_all_loci.png", res=600, width=6000, height=6000)
plot(pca1$li[,1:2], pch=19, col=col_pts, xlab=paste("PC 1 (variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%)",sep=""), cex.lab=1)
abline(v=0,h=0,col="black")
legend(-12, -5, legend=c("Atlantic (n=17)","East (n=21)","West (n=24)"), pch=19, col=c("chocolate2","green","mediumblue"), bty="n", cex=1)
dataEllipse(as.matrix(pca1$li[gen3@strata$G3 == "Atlantic", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2")
dataEllipse(as.matrix(pca1$li[gen3@strata$G3 == "East", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green")
dataEllipse(as.matrix(pca1$li[gen3@strata$G3 == "West", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
#mtext(paste("PC1 variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""), 1, adj=0.05, line=-3)
#mtext(paste("PC2 variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""), 1, adj=0.05, line=-2)
dev.off()


png("PCA_POP_all_loci_apr.png", res=600, width=6000, height=6000)
plot(pca1$li[,1:2], pch=19, col=col_pts, xlab=paste("PC 1 (variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%)",sep=""), cex.lab=1)
abline(v=0,h=0,col="black")
legend(-12, -5, legend=c("Atlantic (n=17)","East (n=21)","West (n=24)"), pch=19, col=c("chocolate2","green","mediumblue"), bty="n", cex=1)
dataEllipse(as.matrix(pca1$li[gen3@strata$G3 == 3, 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2")
dataEllipse(as.matrix(pca1$li[gen3@strata$G3 == 2, 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green")
dataEllipse(as.matrix(pca1$li[gen3@strata$G3 == 1, 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
dev.off()
} #notepad cleanup

#Population Outliers (Outflank)
{```{R}```
setPop(gen3.vcf) <- ~G3
tmp<-gl.outflank(gen3.vcf, qthreshold = 0.1)	#Can use plot=F option if you want
length(which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"))
#70

outliers <- tmp$outflank$results[which(tmp$outflank$results[15]=="TRUE" & tmp$outflank$results[11]=="goodH"),1]
out.list <- matrix(unlist(strsplit(as.vector(outliers),"[.]")),ncol=2,byrow=T)[,1]
tmp <- matrix(unlist(strsplit(as.vector(out.list),"[_]")),ncol=4,byrow=T)
#head(tmp)
out.list <- unique(apply(tmp[,1:3], 1, function(x) paste(x, collapse="_")))
length(out.list)
#27
gen.out.list <- (out.list)

set.loc<-locNames(gen3)[which(!(locNames(gen3) %in% gen.out.list))]
gen.net<-gen3[, loc=set.loc]
set.loc<-locNames(gen3)[which((locNames(gen3) %in% gen.out.list))]
gen.out<-gen3[, loc=set.loc]

X <- scaleGen(gen3, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.out <- scaleGen(gen.out, NA.method="mean", scale=F)
pca.out <- dudi.pca(X.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen3) <- setPop(gen.net) <- setPop(gen.out) <- ~G3
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen3), col=col.All, cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen3))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen.net), col=col.All, cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen.out), col=col.All, cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.out))," loci)", sep=""), 3, 2, adj = 0.95)
legend(2.5, 1.5, legend=c("Atlantic", "East Gulf", "West Gulf"), pch=19, col=col.All, bty="n")

setPop(gen3) <- setPop(gen.net) <- setPop(gen.out) <- ~Lib
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen3), col=col.All, cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen3))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen.net), col=col.All, cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen.out), col=col.All, cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen.out))," loci)", sep=""), 3, 2, adj = 0.95)

write.table(gen.out.list, "Outflank_Fst_Outliers.list", quote=F, col.names=F, row.names=F)
} #notepad cleanup

#Population Outliers (Bayescan)
{```{R}```
#Exporting data for Bayescan
setPop(gen3)<-~G3
writeGenPop(gen3, "G3_groups.gen", "NWA Angel shark data without dups by G3")
} #notepad cleanup
{```{bash}```
#Converting to BS format
java -jar /usr/local/bin/PGDSpider2.1.1.3-cli.jar -inputfile G3_groups.gen -inputformat GENEPOP -outputfile G3_groups.BS -outputformat BAYESCAN -spid /home/afields/bin/genepop_to_BS.spid

#Running Bayescan
bayescan_2.1 G3_groups.BS -od ./Bayescan -all-trace -threads 24 -thin 100 -nbp 30 -pr_odds 100
} #notepad cleanup
{```{R}```
#Checking convergance
#http://evomics.org/wp-content/uploads/2016/01/BayeScan_BayeScEnv_exercises.pdf
library('coda')
library('adegenet')
#Load data
chain <- read.table("G3_group.sel",header=TRUE)
chain <- chain[-c(1)]
chain <- mcmc(chain,thin=10)
#Plot
plot(chain)
#Summary
summary(chain)
effectiveSize(chain)
#Checking Convergance
#Geweke’s convergence diagnostic (values should be between -1.96 and 1.96 for alpha = 0.05)
geweke.diag(chain, frac1=0.1, frac2=0.5)
#Heidelberg and Welch’s convergence diagnostic
heidel.diag(chain, eps=0.1, pvalue=0.05)

#Plotting the Fst relative to the ancestor
chain <- read.table("G3_group.sel",header=TRUE)
chain <- chain[-c(1)]

png("Population_outlier.png", res=200, width=2000, height=2000)
plot(chain$Fst1, pch=19, cex=0.5, col=funky(1), ylim=c(0,0.04))
points(chain$Fst2, pch=19, cex=0.5, col=tail(funky(2),n=1))
points(chain$Fst3, pch=19, cex=0.5, col=funky(3)[2])
dev.off()
} #notepad cleanup
{```{bash}```
#Analyzing results
head -n2 ../G3_groups.gen | tail -n 1 | sed 's/, /\n/g' > contigs
echo "Contigs prob" | cat - contigs > tmp; mv tmp contigs
paste contigs G3_group_fst.txt | awk '{$2=""; print}' > fst.txt

awk 'NR==1{next;} $4 < 0.05{print $0}' fst.txt | wc -l
#11
awk 'NR==1{next;} $4 < 0.05{print $1}' fst.txt | sort > Bayescan_Outliers.list
awk 'NR==1{print $0; next} $4 < 0.05{print $0}' fst.txt | less
} #notepad cleanup

#Looking at effects of Bayescan outliers
{```{R}```
tmp.out <- as.character(as.matrix(read.csv("Bayescan/Bayescan_Outliers.list", head=F)))

set.loc<-locNames(gen3)[which(!(locNames(gen3) %in% tmp.out))]
gen3.net<-gen3[, loc=set.loc]
set.loc<-locNames(gen3)[which((locNames(gen3) %in% tmp.out))]
gen3.out<-gen3[, loc=set.loc]

X <- scaleGen(gen3, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.net <- scaleGen(gen3.net, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.out <- scaleGen(gen3.out, NA.method="mean", scale=F)
pca.out <- dudi.pca(X.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen3) <- setPop(gen3.net) <- setPop(gen3.out) <- ~G3
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen3), col=col.All, cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen3))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen3.net), col=col.All, cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen3.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen3.out), col=col.All, cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen3.out))," loci)", sep=""), 3, 2, adj = 0.95)

setPop(gen3) <- setPop(gen3.net) <- setPop(gen3.out) <- ~Lib
par(mfrow=c(3,1))
ade4::s.class(pca1$li, pop(gen3), col=col.All, cstar=0, axesell=F)
mtext(paste("All data (n=",length(locNames(gen3))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.net$li, pop(gen3.net), col=col.All, cstar=0, axesell=F)
mtext(paste("Neutral data (n=",length(locNames(gen3.net))," loci)", sep=""), 3, 2, adj = 0.95)
ade4::s.class(pca.out$li, pop(gen3.out), col=col.All, cstar=0, axesell=F)
mtext(paste("Outlier data (n=",length(locNames(gen3.out))," loci)", sep=""), 3, 2, adj = 0.95)
} #notepad cleanup

#Getting the neutral and outlier dataset
{```{R}```
OF.out <- read.table("Outflank_Fst_Outliers.list", head=F)
BS.out <- read.table("Bayescan/Bayescan_Outliers.list", head=F)
gen.out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out))))
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% gen.out.list))]

gen4 <- gen3[, loc=set.loc]
gen4.vcf <- VCF_remove(gen3.vcf, gen.out.list)

save(gen4, file="gen4_noPac.gz", compress=T)
save(gen4.vcf, file="gen4.vcf_noPac.gz", compress=T)
#load("gen4_noPac.gz")
#load("gen4.vcf_noPac.gz")
} #notepad cleanup

#RDA
#MEMS do not really make sense due to each group being independent and there is a large distance gap between Atl and Gulf
{```{R}```
#### GeoDist, Environmental factors, strata and RDAs ####

#Getting the geographical distance between samples along the gradient
{```{R}```
## Getting the network
#Setting up data
{```{R}```
#gen3
gen3@other$lat <- gen3@strata[,c("Lon","Lat")]
gen3@other$lat[,1] <- as.numeric(as.matrix(gen3@strata$Lon))	#As numeric
gen3@other$lat[,2] <- as.numeric(as.matrix(gen3@strata$Lat))	#As numeric
rownames(gen3@other$lat) <- indNames(gen3)

tmp.lat <- rbind(gen3@other$lat, data.frame(Lon=c(-81.141374, -79.880862, -80.579324, -83.001458), Lat=c(31.020688, 26.698130, 24.831460, 27.640040)))
} #notepad cleanup
#Getting distance matrix
{```{R}```
dist.mat <- distm(tmp.lat[,1:2], fun=distHaversine)/1000
colnames(dist.mat) <- rownames(dist.mat) <- rownames(tmp.lat)
} #notepad cleanup
#Getting FL corrected distance matrix
{```{R}```
#Getting Unique locations
uniq.lats <- unique(tmp.lat)
#Getting Minimum spanning network
tmp.cn <- chooseCN(uniq.lats, type=4)
#Converting network to matrix
tmp.nb <-neig2mat(nb2neig(tmp.cn))
colnames(tmp.nb) <- rownames(tmp.nb) <- rownames(uniq.lats)
#Removing cross FL connections
tmp.nb[rownames(tmp.nb) %in% gen3@strata$INDV[gen3@strata$G2 == "Gulf"],colnames(tmp.nb) %in% gen3@strata$INDV[gen3@strata$G2 == "Atlantic"]] <- 0
tmp.nb[rownames(tmp.nb) %in% gen3@strata$INDV[gen3@strata$G2 == "Atlantic"],colnames(tmp.nb) %in% gen3@strata$INDV[gen3@strata$G2 == "Gulf"]] <- 0
tmp.nb[rownames(tmp.nb) == 1, rownames(tmp.nb) == 4] <- tmp.nb[rownames(tmp.nb) == 4, rownames(tmp.nb) == 1] <- 0
tmp.nb[rownames(tmp.nb) == 2, rownames(tmp.nb) == 4] <- tmp.nb[rownames(tmp.nb) == 4, rownames(tmp.nb) == 2] <- 0
#Connecting Gulf to Atlantic through Keys
tmp.nb[rownames(tmp.nb) == 1, rownames(tmp.nb) == 2] <- tmp.nb[rownames(tmp.nb) == 2, rownames(tmp.nb) == 1] <- 1
tmp.nb[rownames(tmp.nb) == 2, rownames(tmp.nb) == 3] <- tmp.nb[rownames(tmp.nb) == 3, rownames(tmp.nb) == 2] <- 1
tmp.nb[rownames(tmp.nb) == 3, rownames(tmp.nb) == 4] <- tmp.nb[rownames(tmp.nb) == 4, rownames(tmp.nb) == 3] <- 1

#Plot check
plot(mat2listw(tmp.nb), uniq.lats, bg="red4", pch=21)
mtext("Minimum spanning tree network", 3, cex=2, font=2, line=0.25)
grid(col="grey60")
box()
} #notepad cleanup
#Getting closest distance between the Gulf and Atlantic
{```{R}```
min.tmp <- min(dist.mat.cor[indNames(gen3)[gen3@strata$G2 == "Gulf"], indNames(gen3)[gen3@strata$G2 == "Atlantic"]])
tmp.tab <- apply(dist.mat.cor, 2, function(x) sum(x == min.tmp))
ATL.link <- names(tmp.tab)[which(tmp.tab>0 & names(tmp.tab) %in% indNames(gen3)[gen3@strata$G2 == "Atlantic"])][1]
Gulf.link <- names(tmp.tab)[which(tmp.tab>0 & names(tmp.tab) %in% indNames(gen3)[gen3@strata$G2 == "Gulf"])][1]
} #notepad cleanup

## Getting the distances
#Getting the most NE point
{```{R}```
uniq.lats <- unique(gen3@other$lat[indNames(gen3)[gen3@strata$G3 == "Atlantic"],])
tmp.Start <- rownames(uniq.lats)[which(uniq.lats$Lat == max(uniq.lats$Lat))][1]
} #notepad cleanup
#Making the igraph object
{```{R}```
tmp_mat <- as.matrix(tmp.nb)
tmp_mat[tmp_mat > 0] <- 1
tmp_graph <- graph.adjacency(tmp_mat, mode="undirected", weighted=T)
} #notepad cleanup
#Getting the distance from the NE
{```{R}```
uniq.lats <- unique(tmp.lat)

DIST_m <- NULL
for(i in indNames(gen3)){
if(i == tmp.Start){next}
COORD <- gen3@other$lat[i,]
INDV <- rownames(uniq.lats)[which(uniq.lats[,1] == COORD[,1] & uniq.lats[,2] == COORD[,2])]
ROUTE <- attr(shortest_paths(tmp_graph, tmp.Start, as.character(INDV))$vpath[[1]], "names")
DIST <- 0
if(length(ROUTE) == 1){DIST <- dist.mat[tmp.Start, INDV]
} else {for(j in 1:(length(ROUTE)-1)){
DIST <- DIST + dist.mat[ROUTE[j], ROUTE[j+1]]
}}
DIST_m <- rbind(DIST_m, c(i,DIST))
}
DIST_m <- rbind(DIST_m, c(tmp.Start, 0))
}}}} #notepad cleanup
save.image("RDA_prep.RData.gz", compress=T)
#load("RDA_prep.RData.gz")

#Environmental with geo_dist
{```{R}```
DIST_df <- as.data.frame(DIST_m)

RDA_data <- cbind(env_dat_filter, DIST_df[match(rownames(env_dat_filter),DIST_df[,1]),2])
colnames(RDA_data)[ncol(RDA_data)] <- "Geo_dist"
RDA_data$Geo_dist <- as.numeric(as.matrix(RDA_data$Geo_dist))
RDA_data_scale <- as.data.frame(scale(RDA_data[match(indNames(gen3), rownames(RDA_data)),]))

X <- scaleGen(gen3, NA.method="mean", scale=F)
m1<-rda(X ~ ., RDA_data_scale)
m0<-rda(X ~ 1, RDA_data_scale)
set.seed(1235)
m.ord_All <- ordiR2step(m0, scope = formula(m1), Pin=0.05, Pout=0.1, permutations = 999, parallel=10, R2scope=F)
{ #Results
	R2.adjusted	
+	MS_sst09_5m	1.77E-02
+	Geo_dist	1.75E-02
+	BO_parmean	1.75E-02
+	BO22_parmean	1.75E-02
+	Lat	1.74E-02
+	MS_sst10_5m	1.73E-02
+	BO2_tempmax_ss	1.73E-02
+	BO21_tempmax_ss	1.73E-02
+	BO22_tempmax_ss	1.73E-02
+	MS_sst11_5m	1.72E-02
+	BO2_ppmean_bdmean	1.72E-02
+	BO2_pprange_bdmean	1.72E-02
+	BO21_ppmean_bdmean	1.72E-02
+	BO21_pprange_bdmean	1.72E-02
+	BO22_ppmean_bdmean	1.72E-02
+	BO22_pprange_bdmean	1.72E-02
+	BO2_ppltmax_bdmean	1.72E-02
+	BO21_ppltmax_bdmean	1.72E-02
+	BO22_ppltmax_bdmean	1.72E-02
+	BO2_carbonphytomax_bdmean	1.71E-02
+	BO21_carbonphytomax_bdmean	1.71E-02
+	BO22_carbonphytomax_bdmean	1.71E-02
+	BO2_carbonphytoltmax_bdmean	1.71E-02
+	BO21_carbonphytoltmax_bdmean	1.71E-02
+	BO22_carbonphytoltmax_bdmean	1.71E-02
+	BO2_ppmax_bdmean	1.70E-02
+	BO21_ppmax_bdmean	1.70E-02
+	BO22_ppmax_bdmean	1.70E-02
+	BO22_chlomax_bdmean	1.70E-02
+	BO_silicate	1.70E-02
+	BO22_chloltmax_bdmean	1.69E-02
+	BO2_carbonphytomax_bdmin	1.68E-02
+	BO21_carbonphytomax_bdmin	1.68E-02
+	BO22_carbonphytomax_bdmin	1.68E-02
+	BO2_carbonphytoltmax_bdmin	1.68E-02
+	BO21_carbonphytoltmax_bdmin	1.68E-02
+	BO22_carbonphytoltmax_bdmin	1.68E-02
+	MS_sst06_5m	1.68E-02
+	BO2_templtmax_ss	1.67E-02
+	BO21_templtmax_ss	1.67E-02
+	BO22_templtmax_ss	1.67E-02
+	BO2_tempmin_ss	1.67E-02
+	BO21_tempmin_ss	1.67E-02
+	BO22_tempmin_ss	1.67E-02
+	BO_sstmean	1.67E-02
+	BO2_carbonphytomean_bdmean	1.66E-02
+	BO2_carbonphytorange_bdmean	1.66E-02
+	BO21_carbonphytomean_bdmean	1.66E-02
+	BO21_carbonphytorange_bdmean	1.66E-02
+	BO22_carbonphytomean_bdmean	1.66E-02
+	BO22_carbonphytorange_bdmean	1.66E-02
+	MS_biogeo15_sst_max_5m	1.66E-02
+	MS_sst08_5m	1.66E-02
+	BO22_chlomax_bdmin	1.66E-02
+	BO2_ppltmax_bdmin	1.66E-02
+	BO21_ppltmax_bdmin	1.66E-02
+	BO22_ppltmax_bdmin	1.66E-02
+	MS_biogeo13_sst_mean_5m	1.66E-02
+	MS_sst05_5m	1.65E-02
+	MS_biogeo14_sst_min_5m	1.65E-02
+	BO_sstmin	1.65E-02
+	BO22_chloltmax_bdmin	1.64E-02
+	BO2_ppltmin_bdmean	1.64E-02
+	BO21_ppltmin_bdmean	1.64E-02
+	BO22_ppltmin_bdmean	1.64E-02
+	BO2_ppmax_bdmin	1.64E-02
+	BO21_ppmax_bdmin	1.64E-02
+	BO22_ppmax_bdmin	1.64E-02
+	BO_sstmax	1.64E-02
+	MS_sst02_5m	1.64E-02
+	BO2_tempmean_ss	1.63E-02
+	BO21_tempmean_ss	1.63E-02
+	BO22_tempmean_ss	1.63E-02
+	BO2_carbonphytomin_ss	1.63E-02
+	BO21_carbonphytomin_ss	1.63E-02
+	BO22_carbonphytomin_ss	1.63E-02
+	BO2_chlomax_ss	1.63E-02
+	BO22_chlomax_ss	1.63E-02
+	BO2_chloltmax_ss	1.61E-02
+	BO22_chloltmax_ss	1.61E-02
+	MS_sst07_5m	1.61E-02
+	MS_sst03_5m	1.60E-02
+	BO2_dissoxmax_bdmax	1.60E-02
+	BO2_dissoxmax_bdmean	1.60E-02
+	BO2_dissoxmax_bdmin	1.60E-02
+	BO2_dissoxmean_bdmax	1.60E-02
+	BO2_dissoxmin_bdmax	1.60E-02
+	BO2_dissoxrange_bdmax	1.60E-02
+	BO2_dissoxmax_ss	1.60E-02
+	BO21_dissoxmax_bdmax	1.60E-02
+	BO21_dissoxmax_bdmean	1.60E-02
+	BO21_dissoxmax_bdmin	1.60E-02
+	BO21_dissoxmax_ss	1.60E-02
+	BO21_dissoxmean_bdmax	1.60E-02
+	BO21_dissoxmin_bdmax	1.60E-02
+	BO21_dissoxrange_bdmax	1.60E-02
+	BO22_dissoxmax_bdmax	1.60E-02
+	BO22_dissoxmax_bdmean	1.60E-02
+	BO22_dissoxmax_bdmin	1.60E-02
+	BO22_dissoxmax_ss	1.60E-02
+	BO22_dissoxmean_bdmax	1.60E-02
+	BO22_dissoxmin_bdmax	1.60E-02
+	BO22_dissoxrange_bdmax	1.60E-02
+	BO2_chlorange_ss	1.60E-02
+	BO2_temprange_ss	1.60E-02
+	BO21_temprange_ss	1.60E-02
+	BO22_temprange_ss	1.60E-02
+	Lon	1.60E-02
+	BO2_chlomean_ss	1.59E-02
+	BO22_chlomean_ss	1.59E-02
+	MS_biogeo16_sst_range_5m	1.59E-02
+	BO2_carbonphytomin_bdmean	1.59E-02
+	BO21_carbonphytomin_bdmean	1.59E-02
+	BO22_carbonphytomin_bdmean	1.59E-02
+	BO22_chlomax_bdmax	1.58E-02
+	BO22_chlomean_bdmax	1.58E-02
+	BO22_chlomin_bdmax	1.58E-02
+	BO22_chlorange_bdmax	1.58E-02
+	BO22_chloltmax_bdmax	1.58E-02
+	BO22_chloltmin_bdmax	1.58E-02
+	BO2_dissoxrange_ss	1.58E-02
+	BO21_dissoxrange_ss	1.58E-02
+	BO22_dissoxrange_ss	1.58E-02
+	BO_sstrange	1.57E-02
+	BO2_carbonphytoltmin_bdmean	1.57E-02
+	BO21_carbonphytoltmin_bdmean	1.57E-02
+	BO22_carbonphytoltmin_bdmean	1.57E-02
+	BO22_chlomean_bdmean	1.56E-02
+	BO22_chlorange_bdmean	1.56E-02
+	MS_sst04_5m	1.56E-02
+	BO2_templtmin_ss	1.56E-02
+	BO21_templtmin_ss	1.56E-02
+	BO22_templtmin_ss	1.56E-02
+	BO2_carbonphytoltmax_ss	1.55E-02
+	BO21_carbonphytoltmax_ss	1.55E-02
+	BO22_carbonphytoltmax_ss	1.55E-02
+	BO2_chlomin_ss	1.55E-02
+	BO22_chlomin_ss	1.55E-02
+	BO2_carbonphytomax_bdmax	1.55E-02
+	BO2_carbonphytomean_bdmax	1.55E-02
+	BO2_carbonphytomin_bdmax	1.55E-02
+	BO2_carbonphytorange_bdmax	1.55E-02
+	BO21_carbonphytomax_bdmax	1.55E-02
+	BO21_carbonphytomean_bdmax	1.55E-02
+	BO21_carbonphytomin_bdmax	1.55E-02
+	BO21_carbonphytorange_bdmax	1.55E-02
+	BO22_carbonphytomax_bdmax	1.55E-02
+	BO22_carbonphytomean_bdmax	1.55E-02
+	BO22_carbonphytomin_bdmax	1.55E-02
+	BO22_carbonphytorange_bdmax	1.55E-02
+	BO_parmax	1.55E-02
+	BO22_parmax	1.55E-02
+	BO_damean	1.55E-02
+	BO22_damean	1.55E-02
+	BO2_carbonphytoltmin_ss	1.55E-02
+	BO21_carbonphytoltmin_ss	1.55E-02
+	BO22_carbonphytoltmin_ss	1.55E-02
+	BO2_carbonphytoltmax_bdmax	1.54E-02
+	BO2_carbonphytoltmin_bdmax	1.54E-02
+	BO21_carbonphytoltmax_bdmax	1.54E-02
+	BO21_carbonphytoltmin_bdmax	1.54E-02
+	BO22_carbonphytoltmax_bdmax	1.54E-02
+	BO22_carbonphytoltmin_bdmax	1.54E-02
+	BO2_ppmax_bdmax	1.54E-02
+	BO2_ppmean_bdmax	1.54E-02
+	BO2_ppmin_bdmax	1.54E-02
+	BO2_pprange_bdmax	1.54E-02
+	BO21_ppmax_bdmax	1.54E-02
+	BO21_ppmean_bdmax	1.54E-02
+	BO21_ppmin_bdmax	1.54E-02
+	BO21_pprange_bdmax	1.54E-02
+	BO22_ppmax_bdmax	1.54E-02
+	BO22_ppmean_bdmax	1.54E-02
+	BO22_ppmin_bdmax	1.54E-02
+	BO22_pprange_bdmax	1.54E-02
+	MS_sst01_5m	1.53E-02
+	BO2_carbonphytomax_ss	1.52E-02
+	BO21_carbonphytomax_ss	1.52E-02
+	BO22_carbonphytomax_ss	1.52E-02
+	BO2_ppltmax_bdmax	1.52E-02
+	BO2_ppltmin_bdmax	1.52E-02
+	BO21_ppltmax_bdmax	1.52E-02
+	BO21_ppltmin_bdmax	1.52E-02
+	BO22_ppltmax_bdmax	1.52E-02
+	BO22_ppltmin_bdmax	1.52E-02
+	BO2_ppmin_bdmean	1.52E-02
+	BO21_ppmin_bdmean	1.52E-02
+	BO22_ppmin_bdmean	1.52E-02
+	BO2_carbonphytomean_ss	1.52E-02
+	BO21_carbonphytomean_ss	1.52E-02
+	BO22_carbonphytomean_ss	1.52E-02
+	BO2_salinityltmax_bdmin	1.51E-02
+	BO21_salinityltmax_bdmin	1.51E-02
+	BO22_salinityltmax_bdmin	1.51E-02
+	BO_ph	1.51E-02
+	BO22_ph	1.51E-02
+	BO_calcite	1.50E-02
+	BO22_calcite	1.50E-02
+	BO2_dissoxmean_bdmin	1.49E-02
+	BO2_dissoxmin_bdmean	1.49E-02
+	BO2_dissoxmin_bdmin	1.49E-02
+	BO2_dissoxrange_bdmin	1.49E-02
+	BO2_dissoxmin_ss	1.49E-02
+	BO21_dissoxmean_bdmin	1.49E-02
+	BO21_dissoxmin_bdmean	1.49E-02
+	BO21_dissoxmin_bdmin	1.49E-02
+	BO21_dissoxmin_ss	1.49E-02
+	BO21_dissoxrange_bdmin	1.49E-02
+	BO22_dissoxmean_bdmin	1.49E-02
+	BO22_dissoxmin_bdmean	1.49E-02
+	BO22_dissoxmin_bdmin	1.49E-02
+	BO22_dissoxmin_ss	1.49E-02
+	BO22_dissoxrange_bdmin	1.49E-02
+	BO2_ppltmin_bdmin	1.48E-02
+	BO21_ppltmin_bdmin	1.48E-02
+	BO22_ppltmin_bdmin	1.48E-02
+	MS_sss11_5m	1.48E-02
+	BO2_salinitymax_bdmin	1.48E-02
+	BO21_salinitymax_bdmin	1.48E-02
+	BO22_salinitymax_bdmin	1.48E-02
+	BO2_carbonphytorange_ss	1.48E-02
+	BO21_carbonphytorange_ss	1.48E-02
+	BO22_carbonphytorange_ss	1.48E-02
+	BO_damin	1.47E-02
+	BO22_damin	1.47E-02
+	MS_sss02_5m	1.47E-02
+	BO2_salinitymean_bdmin	1.45E-02
+	BO2_salinitymin_bdmin	1.45E-02
+	BO2_salinityrange_bdmin	1.45E-02
+	BO21_salinitymean_bdmin	1.45E-02
+	BO21_salinitymin_bdmin	1.45E-02
+	BO21_salinityrange_bdmin	1.45E-02
+	BO22_salinitymean_bdmin	1.45E-02
+	BO22_salinitymin_bdmin	1.45E-02
+	BO22_salinityrange_bdmin	1.45E-02
+	BO2_carbonphytomean_bdmin	1.44E-02
+	BO2_carbonphytomin_bdmin	1.44E-02
+	BO2_carbonphytorange_bdmin	1.44E-02
+	BO21_carbonphytomean_bdmin	1.44E-02
+	BO21_carbonphytomin_bdmin	1.44E-02
+	BO21_carbonphytorange_bdmin	1.44E-02
+	BO22_carbonphytomean_bdmin	1.44E-02
+	BO22_carbonphytomin_bdmin	1.44E-02
+	BO22_carbonphytorange_bdmin	1.44E-02
+	BO2_salinityltmin_bdmin	1.44E-02
+	BO21_salinityltmin_bdmin	1.44E-02
+	BO22_salinityltmin_bdmin	1.44E-02
+	BO_salinity	1.44E-02
+	BO2_chloltmin_ss	1.43E-02
+	BO22_chloltmin_ss	1.43E-02
+	MS_biogeo10_sss_max_5m	1.43E-02
+	BO2_dissoxltmin_bdmean	1.42E-02
+	BO2_dissoxltmin_bdmin	1.42E-02
+	BO2_dissoxltmin_ss	1.42E-02
+	BO21_dissoxltmin_bdmean	1.42E-02
+	BO21_dissoxltmin_bdmin	1.42E-02
+	BO21_dissoxltmin_ss	1.42E-02
+	BO22_dissoxltmin_bdmean	1.42E-02
+	BO22_dissoxltmin_bdmin	1.42E-02
+	BO22_dissoxltmin_ss	1.42E-02
+	BO2_carbonphytoltmin_bdmin	1.41E-02
+	BO21_carbonphytoltmin_bdmin	1.41E-02
+	BO22_carbonphytoltmin_bdmin	1.41E-02
+	BO_damax	1.38E-02
+	BO22_damax	1.38E-02
+	MS_biogeo17_sst_variance_5m	1.36E-02
+	BO2_ppmean_ss	1.35E-02
+	BO21_ppmean_ss	1.35E-02
+	BO22_ppmean_ss	1.35E-02
+	BO_dissox	1.34E-02
+	BO2_ppltmin_ss	1.32E-02
+	BO21_ppltmin_ss	1.32E-02
+	BO22_ppltmin_ss	1.32E-02
+	MS_sss09_5m	1.32E-02
+	MS_biogeo05_dist_shore_5m	1.30E-02
+	BO2_ppmax_ss	1.30E-02
+	BO21_ppmax_ss	1.30E-02
+	BO22_ppmax_ss	1.30E-02
+	MS_sss10_5m	1.29E-02
+	BO2_lightbotmean_bdmin	1.28E-02
+	BO2_lightbotmin_bdmin	1.28E-02
+	BO2_lightbotrange_bdmin	1.28E-02
+	BO21_lightbotmean_bdmin	1.28E-02
+	BO21_lightbotmin_bdmin	1.28E-02
+	BO21_lightbotrange_bdmin	1.28E-02
+	BO22_lightbotmean_bdmin	1.28E-02
+	BO22_lightbotmin_bdmin	1.28E-02
+	BO22_lightbotrange_bdmin	1.28E-02
+	BO2_ppmean_bdmin	1.28E-02
+	BO2_ppmin_bdmin	1.28E-02
+	BO2_pprange_bdmin	1.28E-02
+	BO21_ppmean_bdmin	1.28E-02
+	BO21_ppmin_bdmin	1.28E-02
+	BO21_pprange_bdmin	1.28E-02
+	BO22_ppmean_bdmin	1.28E-02
+	BO22_ppmin_bdmin	1.28E-02
+	BO22_pprange_bdmin	1.28E-02
+	BO2_phosphateltmin_ss	1.28E-02
+	BO21_phosphateltmin_ss	1.28E-02
+	BO22_phosphateltmin_ss	1.28E-02
+	BO2_dissoxmean_bdmean	1.27E-02
+	BO2_dissoxrange_bdmean	1.27E-02
+	BO2_dissoxmean_ss	1.27E-02
+	BO21_dissoxmean_bdmean	1.27E-02
+	BO21_dissoxmean_ss	1.27E-02
+	BO21_dissoxrange_bdmean	1.27E-02
+	BO22_dissoxmean_bdmean	1.27E-02
+	BO22_dissoxmean_ss	1.27E-02
+	BO22_dissoxrange_bdmean	1.27E-02
+	MS_biogeo08_sss_mean_5m	1.27E-02
+	BO2_dissoxltmax_bdmax	1.27E-02
+	BO2_dissoxltmax_bdmean	1.27E-02
+	BO2_dissoxltmax_bdmin	1.27E-02
+	BO2_dissoxltmin_bdmax	1.27E-02
+	BO2_dissoxltmax_ss	1.27E-02
+	BO21_dissoxltmax_bdmax	1.27E-02
+	BO21_dissoxltmax_bdmean	1.27E-02
+	BO21_dissoxltmax_bdmin	1.27E-02
+	BO21_dissoxltmax_ss	1.27E-02
+	BO21_dissoxltmin_bdmax	1.27E-02
+	BO22_dissoxltmax_bdmax	1.27E-02
+	BO22_dissoxltmax_bdmean	1.27E-02
+	BO22_dissoxltmax_bdmin	1.27E-02
+	BO22_dissoxltmax_ss	1.27E-02
+	BO22_dissoxltmin_bdmax	1.27E-02
+	BO2_ppmin_ss	1.23E-02
+	BO21_ppmin_ss	1.23E-02
+	BO22_ppmin_ss	1.23E-02
+	MS_biogeo12_sss_variance_5m	1.23E-02
+	BO2_nitratemax_bdmax	1.22E-02
+	BO2_nitratemean_bdmax	1.22E-02
+	BO2_nitratemin_bdmax	1.22E-02
+	BO2_nitraterange_bdmax	1.22E-02
+	BO21_nitratemax_bdmax	1.22E-02
+	BO21_nitratemean_bdmax	1.22E-02
+	BO21_nitratemin_bdmax	1.22E-02
+	BO21_nitraterange_bdmax	1.22E-02
+	BO22_nitratemax_bdmax	1.22E-02
+	BO22_nitratemean_bdmax	1.22E-02
+	BO22_nitratemin_bdmax	1.22E-02
+	BO22_nitraterange_bdmax	1.22E-02
+	BO2_lightbotltmin_bdmin	1.21E-02
+	BO21_lightbotltmin_bdmin	1.21E-02
+	BO22_lightbotltmin_bdmin	1.21E-02
+	BO2_phosphatemax_bdmax	1.21E-02
+	BO2_phosphatemean_bdmax	1.21E-02
+	BO2_phosphatemin_bdmax	1.21E-02
+	BO2_phosphaterange_bdmax	1.21E-02
+	BO21_phosphatemax_bdmax	1.21E-02
+	BO21_phosphatemean_bdmax	1.21E-02
+	BO21_phosphatemin_bdmax	1.21E-02
+	BO21_phosphaterange_bdmax	1.21E-02
+	BO22_phosphatemax_bdmax	1.21E-02
+	BO22_phosphatemean_bdmax	1.21E-02
+	BO22_phosphatemin_bdmax	1.21E-02
+	BO22_phosphaterange_bdmax	1.21E-02
+	MS_sss03_5m	1.21E-02
+	BO2_ironmax_bdmin	1.20E-02
+	BO21_ironmax_bdmin	1.20E-02
+	BO22_ironmax_bdmin	1.20E-02
+	BO2_ironmax_bdmean	1.19E-02
+	BO21_ironmax_bdmean	1.19E-02
+	BO22_ironmax_bdmean	1.19E-02
+	BO2_ppltmax_ss	1.19E-02
+	BO21_ppltmax_ss	1.19E-02
+	BO22_ppltmax_ss	1.19E-02
+	BO2_phosphatemin_ss	1.19E-02
+	BO21_phosphatemin_ss	1.19E-02
+	BO22_phosphatemin_ss	1.19E-02
+	MS_biogeo11_sss_range_5m	1.18E-02
+	BO2_phosphateltmax_bdmax	1.18E-02
+	BO2_phosphateltmin_bdmax	1.18E-02
+	BO21_phosphateltmax_bdmax	1.18E-02
+	BO21_phosphateltmin_bdmax	1.18E-02
+	BO22_phosphateltmax_bdmax	1.18E-02
+	BO22_phosphateltmin_bdmax	1.18E-02
+	BO2_pprange_ss	1.17E-02
+	BO21_pprange_ss	1.17E-02
+	BO22_pprange_ss	1.17E-02
+	BO2_nitrateltmax_bdmax	1.17E-02
+	BO2_nitrateltmin_bdmax	1.17E-02
+	BO21_nitrateltmax_bdmax	1.17E-02
+	BO21_nitrateltmin_bdmax	1.17E-02
+	BO22_nitrateltmax_bdmax	1.17E-02
+	BO22_nitrateltmin_bdmax	1.17E-02
+	MS_sss12_5m	1.16E-02
+	MS_sss06_5m	1.16E-02
+	BO22_chloltmin_bdmean	1.15E-02
+	BO2_salinitymin_bdmean	1.13E-02
+	BO21_salinitymin_bdmean	1.13E-02
+	BO22_salinitymin_bdmean	1.13E-02
+	BO2_tempmax_bdmax	1.13E-02
+	BO2_tempmean_bdmax	1.13E-02
+	BO2_tempmin_bdmax	1.13E-02
+	BO2_temprange_bdmax	1.13E-02
+	BO21_tempmax_bdmax	1.13E-02
+	BO21_tempmean_bdmax	1.13E-02
+	BO21_tempmin_bdmax	1.13E-02
+	BO21_temprange_bdmax	1.13E-02
+	BO22_tempmax_bdmax	1.13E-02
+	BO22_tempmean_bdmax	1.13E-02
+	BO22_tempmin_bdmax	1.13E-02
+	BO22_temprange_bdmax	1.13E-02
+	BO22_chlomin_bdmean	1.12E-02
+	BO2_templtmax_bdmax	1.11E-02
+	BO2_templtmin_bdmax	1.11E-02
+	BO21_templtmax_bdmax	1.11E-02
+	BO21_templtmin_bdmax	1.11E-02
+	BO22_templtmax_bdmax	1.11E-02
+	BO22_templtmin_bdmax	1.11E-02
+	BO2_ironltmax_bdmin	1.11E-02
+	BO21_ironltmax_bdmin	1.11E-02
+	BO22_ironltmax_bdmin	1.11E-02
+	BO2_ironltmax_bdmean	1.10E-02
+	BO21_ironltmax_bdmean	1.10E-02
+	BO22_ironltmax_bdmean	1.10E-02
+	MS_sss01_5m	1.10E-02
+	BO2_phosphatemax_bdmean	1.10E-02
+	BO21_phosphatemax_bdmean	1.10E-02
+	BO22_phosphatemax_bdmean	1.10E-02
+	BO2_tempmax_bdmean	1.08E-02
+	BO21_tempmax_bdmean	1.08E-02
+	BO22_tempmax_bdmean	1.08E-02
+	BO2_salinityltmin_bdmean	1.08E-02
+	BO21_salinityltmin_bdmean	1.08E-02
+	BO22_salinityltmin_bdmean	1.08E-02
+	BO2_templtmax_bdmean	1.05E-02
+	BO21_templtmax_bdmean	1.05E-02
+	BO22_templtmax_bdmean	1.05E-02
+	BO2_tempmean_bdmin	1.05E-02
+	BO2_tempmin_bdmin	1.05E-02
+	BO2_temprange_bdmin	1.05E-02
+	BO21_tempmean_bdmin	1.05E-02
+	BO21_tempmin_bdmin	1.05E-02
+	BO21_temprange_bdmin	1.05E-02
+	BO22_tempmean_bdmin	1.05E-02
+	BO22_tempmin_bdmin	1.05E-02
+	BO22_temprange_bdmin	1.05E-02
+	BO2_nitratemax_bdmean	1.05E-02
+	BO21_nitratemax_bdmean	1.05E-02
+	BO22_nitratemax_bdmean	1.05E-02
+	BO2_salinitymean_bdmean	1.04E-02
+	BO2_salinityrange_bdmean	1.04E-02
+	BO21_salinitymean_bdmean	1.04E-02
+	BO21_salinityrange_bdmean	1.04E-02
+	BO22_salinitymean_bdmean	1.04E-02
+	BO22_salinityrange_bdmean	1.04E-02
+	BO2_phosphateltmax_bdmean	1.04E-02
+	BO21_phosphateltmax_bdmean	1.04E-02
+	BO22_phosphateltmax_bdmean	1.04E-02
+	BO_cloudmax	1.02E-02
+	BO22_cloudmax	1.02E-02
+	BO2_ironmean_bdmean	1.02E-02
+	BO2_ironrange_bdmean	1.02E-02
+	BO21_ironmean_bdmean	1.02E-02
+	BO21_ironrange_bdmean	1.02E-02
+	BO22_ironmean_bdmean	1.02E-02
+	BO22_ironrange_bdmean	1.02E-02
+	BO2_salinityltmax_bdmean	1.00E-02
+	BO21_salinityltmax_bdmean	1.00E-02
+	BO22_salinityltmax_bdmean	1.00E-02
+	BO2_ironltmin_bdmin	9.95E-03
+	BO21_ironltmin_bdmin	9.95E-03
+	BO22_ironltmin_bdmin	9.95E-03
+	BO_bathymax	9.87E-03
+	BO2_nitrateltmax_bdmean	9.78E-03
+	BO21_nitrateltmax_bdmean	9.78E-03
+	BO22_nitrateltmax_bdmean	9.78E-03
+	BO_nitrate	9.77E-03
+	BO2_salinitymax_bdmean	9.73E-03
+	BO21_salinitymax_bdmean	9.73E-03
+	BO22_salinitymax_bdmean	9.73E-03
+	BO2_curvelmean_ss	9.68E-03
+	BO2_tempmax_bdmin	9.50E-03
+	BO21_tempmax_bdmin	9.50E-03
+	BO22_tempmax_bdmin	9.50E-03
+	BO2_phosphatemean_bdmean	9.40E-03
+	BO2_phosphaterange_bdmean	9.40E-03
+	BO21_phosphatemean_bdmean	9.40E-03
+	BO21_phosphaterange_bdmean	9.40E-03
+	BO22_phosphatemean_bdmean	9.40E-03
+	BO22_phosphaterange_bdmean	9.40E-03
+	BO2_ironmax_bdmax	9.28E-03
+	BO2_ironmean_bdmax	9.28E-03
+	BO2_ironmin_bdmax	9.28E-03
+	BO2_ironrange_bdmax	9.28E-03
+	BO21_ironmax_bdmax	9.28E-03
+	BO21_ironmean_bdmax	9.28E-03
+	BO21_ironmin_bdmax	9.28E-03
+	BO21_ironrange_bdmax	9.28E-03
+	BO22_ironmax_bdmax	9.28E-03
+	BO22_ironmean_bdmax	9.28E-03
+	BO22_ironmin_bdmax	9.28E-03
+	BO22_ironrange_bdmax	9.28E-03
+	BO_chlomean	9.27E-03
+	BO2_templtmax_bdmin	9.21E-03
+	BO21_templtmax_bdmin	9.21E-03
+	BO22_templtmax_bdmin	9.21E-03
+	BO_bathymean	9.07E-03
+	BO2_nitratemean_bdmean	9.06E-03
+	BO2_nitraterange_bdmean	9.06E-03
+	BO21_nitratemean_bdmean	9.06E-03
+	BO21_nitraterange_bdmean	9.06E-03
+	BO22_nitratemean_bdmean	9.06E-03
+	BO22_nitraterange_bdmean	9.06E-03
+	MS_bathy_5m	9.05E-03
+	BO2_phosphatemax_bdmin	9.03E-03
+	BO21_phosphatemax_bdmin	9.03E-03
+	BO22_phosphatemax_bdmin	9.03E-03
+	BO2_phosphatemean_ss	8.82E-03
+	BO21_phosphatemean_ss	8.82E-03
+	BO22_phosphatemean_ss	8.82E-03
+	BO2_ironmean_bdmin	8.69E-03
+	BO2_ironmin_bdmin	8.69E-03
+	BO2_ironrange_bdmin	8.69E-03
+	BO21_ironmean_bdmin	8.69E-03
+	BO21_ironmin_bdmin	8.69E-03
+	BO21_ironrange_bdmin	8.69E-03
+	BO22_ironmean_bdmin	8.69E-03
+	BO22_ironmin_bdmin	8.69E-03
+	BO22_ironrange_bdmin	8.69E-03
+	BO_chlomin	8.52E-03
+	BO2_ironltmax_bdmax	8.51E-03
+	BO2_ironltmin_bdmax	8.51E-03
+	BO21_ironltmax_bdmax	8.51E-03
+	BO21_ironltmin_bdmax	8.51E-03
+	BO22_ironltmax_bdmax	8.51E-03
+	BO22_ironltmin_bdmax	8.51E-03
+	BO2_salinityltmax_ss	8.42E-03
+	BO21_salinityltmax_ss	8.42E-03
+	BO22_salinityltmax_ss	8.42E-03
+	BO2_phosphateltmin_bdmean	8.26E-03
+	BO21_phosphateltmin_bdmean	8.26E-03
+	BO22_phosphateltmin_bdmean	8.26E-03
+	BO2_lightbotmin_bdmean	8.25E-03
+	BO21_lightbotmin_bdmean	8.25E-03
+	BO22_lightbotmin_bdmean	8.25E-03
+	BO2_ironltmin_bdmean	8.19E-03
+	BO21_ironltmin_bdmean	8.19E-03
+	BO22_ironltmin_bdmean	8.19E-03
+	BO2_nitrateltmin_bdmean	8.17E-03
+	BO21_nitrateltmin_bdmean	8.17E-03
+	BO22_nitrateltmin_bdmean	8.17E-03
+	BO2_phosphateltmax_bdmin	8.16E-03
+	BO21_phosphateltmax_bdmin	8.16E-03
+	BO22_phosphateltmax_bdmin	8.16E-03
+	BO_bathymin	8.14E-03
+	BO2_templtmin_bdmin	8.08E-03
+	BO21_templtmin_bdmin	8.08E-03
+	BO22_templtmin_bdmin	8.08E-03
+	BO_cloudmean	7.97E-03
+	BO22_cloudmean	7.97E-03
+	MS_biogeo02_aspect_NS_5m	7.43E-03
+	BO2_tempmin_bdmean	7.39E-03
+	BO21_tempmin_bdmean	7.39E-03
+	BO22_tempmin_bdmean	7.39E-03
+	BO_chlomax	7.36E-03
+	BO21_curvelltmin_ss	7.26E-03
+	BO22_curvelltmin_ss	7.26E-03
+	MS_sss08_5m	7.18E-03
+	BO2_nitratemax_bdmin	6.95E-03
+	BO21_nitratemax_bdmin	6.95E-03
+	BO22_nitratemax_bdmin	6.95E-03
+	BO2_salinitymax_ss	6.91E-03
+	BO21_salinitymax_ss	6.91E-03
+	BO22_salinitymax_ss	6.91E-03
+	BO2_nitratemin_bdmean	6.77E-03
+	BO21_nitratemin_bdmean	6.77E-03
+	BO22_nitratemin_bdmean	6.77E-03
+	BO2_phosphatemin_bdmean	6.77E-03
+	BO21_phosphatemin_bdmean	6.77E-03
+	BO22_phosphatemin_bdmean	6.77E-03
+	BO2_curvelmax_bdmax	6.49E-03
+	BO2_curvelmean_bdmax	6.49E-03
+	BO2_curvelmin_bdmax	6.49E-03
+	BO2_curvelrange_bdmax	6.49E-03
+	BO2_nitrateltmax_bdmin	6.36E-03
+	BO21_nitrateltmax_bdmin	6.36E-03
+	BO22_nitrateltmax_bdmin	6.36E-03
+	BO2_curvelltmax_bdmax	6.35E-03
+	BO2_curvelltmin_bdmax	6.35E-03
+	BO2_salinitymean_ss	6.25E-03
+	BO21_salinitymean_ss	6.25E-03
+	BO22_salinitymean_ss	6.25E-03
+	MS_sss04_5m	6.22E-03
+	BO22_chloltmin_bdmin	6.15E-03
+	MS_sss07_5m	6.13E-03
+	MS_sss05_5m	6.11E-03
+	BO2_lightbotltmax_bdmax	5.84E-03
+	BO2_lightbotltmin_bdmax	5.84E-03
+	BO21_lightbotltmax_bdmax	5.84E-03
+	BO21_lightbotltmin_bdmax	5.84E-03
+	BO22_lightbotltmax_bdmax	5.84E-03
+	BO22_lightbotltmin_bdmax	5.84E-03
+	BO2_phosphateltmin_bdmin	5.79E-03
+	BO21_phosphateltmin_bdmin	5.79E-03
+	BO22_phosphateltmin_bdmin	5.79E-03
+	BO2_phosphateltmax_ss	5.78E-03
+	BO21_phosphateltmax_ss	5.78E-03
+	BO22_phosphateltmax_ss	5.78E-03
+	BO2_ironmin_bdmean	5.77E-03
+	BO21_ironmin_bdmean	5.77E-03
+	BO22_ironmin_bdmean	5.77E-03
+	BO2_lightbotltmin_bdmean	5.72E-03
+	BO21_lightbotltmin_bdmean	5.72E-03
+	BO22_lightbotltmin_bdmean	5.72E-03
+	BO22_chlomean_bdmin	5.64E-03
+	BO22_chlomin_bdmin	5.64E-03
+	BO22_chlorange_bdmin	5.64E-03
+	BO_chlorange	5.58E-03
+	BO2_ironmin_ss	5.48E-03
+	BO21_ironmin_ss	5.48E-03
+	BO22_ironmin_ss	5.48E-03
+	BO2_nitrateltmin_bdmin	5.38E-03
+	BO21_nitrateltmin_bdmin	5.38E-03
+	BO22_nitrateltmin_bdmin	5.38E-03
+	MS_biogeo09_sss_min_5m	5.33E-03
+	MS_biogeo06_bathy_slope_5m	5.12E-03
+	BO2_lightbotltmax_bdmin	4.71E-03
+	BO21_lightbotltmax_bdmin	4.71E-03
+	BO22_lightbotltmax_bdmin	4.71E-03
+	BO2_salinityltmax_bdmax	4.56E-03
+	BO2_salinityltmin_bdmax	4.56E-03
+	BO21_salinityltmax_bdmax	4.56E-03
+	BO21_salinityltmin_bdmax	4.56E-03
+	BO22_salinityltmax_bdmax	4.56E-03
+	BO22_salinityltmin_bdmax	4.56E-03
+	BO2_salinitymax_bdmax	4.41E-03
+	BO2_salinitymean_bdmax	4.41E-03
+	BO2_salinitymin_bdmax	4.41E-03
+	BO2_salinityrange_bdmax	4.41E-03
+	BO21_salinitymax_bdmax	4.41E-03
+	BO21_salinitymean_bdmax	4.41E-03
+	BO21_salinitymin_bdmax	4.41E-03
+	BO21_salinityrange_bdmax	4.41E-03
+	BO22_salinitymax_bdmax	4.41E-03
+	BO22_salinitymean_bdmax	4.41E-03
+	BO22_salinitymin_bdmax	4.41E-03
+	BO22_salinityrange_bdmax	4.41E-03
+	BO2_ironltmin_ss	4.34E-03
+	BO21_ironltmin_ss	4.34E-03
+	BO22_ironltmin_ss	4.34E-03
+	BO2_nitraterange_ss	4.27E-03
+	BO21_nitraterange_ss	4.27E-03
+	BO22_nitraterange_ss	4.27E-03
+	BO2_phosphatemean_bdmin	4.23E-03
+	BO2_phosphatemin_bdmin	4.23E-03
+	BO2_phosphaterange_bdmin	4.23E-03
+	BO21_phosphatemean_bdmin	4.23E-03
+	BO21_phosphatemin_bdmin	4.23E-03
+	BO21_phosphaterange_bdmin	4.23E-03
+	BO22_phosphatemean_bdmin	4.23E-03
+	BO22_phosphatemin_bdmin	4.23E-03
+	BO22_phosphaterange_bdmin	4.23E-03
+	BO2_curvelltmax_bdmean	4.22E-03
+	BO2_templtmin_bdmean	4.19E-03
+	BO21_templtmin_bdmean	4.19E-03
+	BO22_templtmin_bdmean	4.19E-03
+	BO21_curvelmean_ss	4.14E-03
+	BO22_curvelmean_ss	4.14E-03
+	BO2_nitratemax_ss	4.12E-03
+	BO21_nitratemax_ss	4.12E-03
+	BO22_nitratemax_ss	4.12E-03
+	BO2_curvelmax_bdmean	4.08E-03
+	BO_phosphate	4.06E-03
+	BO2_curvelmean_bdmean	4.00E-03
+	BO2_curvelrange_bdmean	4.00E-03
+	BO2_nitratemean_bdmin	4.00E-03
+	BO2_nitratemin_bdmin	4.00E-03
+	BO2_nitraterange_bdmin	4.00E-03
+	BO21_nitratemean_bdmin	4.00E-03
+	BO21_nitratemin_bdmin	4.00E-03
+	BO21_nitraterange_bdmin	4.00E-03
+	BO22_nitratemean_bdmin	4.00E-03
+	BO22_nitratemin_bdmin	4.00E-03
+	BO22_nitraterange_bdmin	4.00E-03
+	BO2_phosphatemax_ss	3.97E-03
+	BO21_phosphatemax_ss	3.97E-03
+	BO22_phosphatemax_ss	3.97E-03
+	BO2_phosphaterange_ss	3.94E-03
+	BO21_phosphaterange_ss	3.94E-03
+	BO22_phosphaterange_ss	3.94E-03
+	BO2_lightbotmax_bdmax	3.90E-03
+	BO2_lightbotmean_bdmax	3.90E-03
+	BO2_lightbotmin_bdmax	3.90E-03
+	BO2_lightbotrange_bdmax	3.90E-03
+	BO21_lightbotmax_bdmax	3.90E-03
+	BO21_lightbotmean_bdmax	3.90E-03
+	BO21_lightbotmin_bdmax	3.90E-03
+	BO21_lightbotrange_bdmax	3.90E-03
+	BO22_lightbotmax_bdmax	3.90E-03
+	BO22_lightbotmean_bdmax	3.90E-03
+	BO22_lightbotmin_bdmax	3.90E-03
+	BO22_lightbotrange_bdmax	3.90E-03
+	BO2_curvelltmin_bdmean	3.78E-03
+	BO2_nitrateltmax_ss	3.70E-03
+	BO21_nitrateltmax_ss	3.70E-03
+	BO22_nitrateltmax_ss	3.70E-03
+	BO2_silicatemax_bdmin	3.43E-03
+	BO21_silicatemax_bdmin	3.43E-03
+	BO22_silicatemax_bdmin	3.43E-03
+	BO21_curvelltmax_ss	3.28E-03
+	BO22_curvelltmax_ss	3.28E-03
+	BO2_lightbotmean_bdmean	3.22E-03
+	BO2_lightbotrange_bdmean	3.22E-03
+	BO21_lightbotmean_bdmean	3.22E-03
+	BO21_lightbotrange_bdmean	3.22E-03
+	BO22_lightbotmean_bdmean	3.22E-03
+	BO22_lightbotrange_bdmean	3.22E-03
+	BO2_nitratemean_ss	3.13E-03
+	BO21_nitratemean_ss	3.13E-03
+	BO22_nitratemean_ss	3.13E-03
+	BO21_curvelmean_bdmin	3.10E-03
+	BO21_curvelmin_bdmin	3.10E-03
+	BO21_curvelrange_bdmin	3.10E-03
+	BO22_curvelmean_bdmin	3.10E-03
+	BO22_curvelmin_bdmin	3.10E-03
+	BO22_curvelrange_bdmin	3.10E-03
+	BO2_salinityltmin_ss	3.07E-03
+	BO21_salinityltmin_ss	3.07E-03
+	BO22_salinityltmin_ss	3.07E-03
+	BO2_silicateltmax_bdmin	2.82E-03
+	BO21_silicateltmax_bdmin	2.82E-03
+	BO22_silicateltmax_bdmin	2.82E-03
+	BO2_salinitymin_ss	2.80E-03
+	BO21_salinitymin_ss	2.80E-03
+	BO22_salinitymin_ss	2.80E-03
+	BO21_curvelrange_ss	2.80E-03
+	BO22_curvelrange_ss	2.80E-03
+	BO21_curvelmax_ss	2.75E-03
+	BO22_curvelmax_ss	2.75E-03
+	BO2_silicaterange_ss	2.71E-03
+	BO21_silicaterange_ss	2.71E-03
+	BO22_silicaterange_ss	2.71E-03
+	BO2_curvelmin_bdmean	2.67E-03
+	BO2_curvelmax_bdmin	2.63E-03
+	BO2_lightbotmax_bdmin	2.62E-03
+	BO21_lightbotmax_bdmin	2.62E-03
+	BO22_lightbotmax_bdmin	2.62E-03
+	BO2_curvelltmax_bdmin	2.43E-03
+	BO2_curvelltmin_bdmin	2.38E-03
+	BO2_salinityrange_ss	2.30E-03
+	BO21_salinityrange_ss	2.30E-03
+	BO22_salinityrange_ss	2.30E-03
+	BO2_nitrateltmin_ss	2.26E-03
+	BO21_nitrateltmin_ss	2.26E-03
+	BO22_nitrateltmin_ss	2.26E-03
+	BO2_silicatemax_ss	2.22E-03
+	BO21_silicatemax_ss	2.22E-03
+	BO22_silicatemax_ss	2.22E-03
+	BO2_curvelltmin_ss	2.09E-03
+	BO2_lightbotltmax_bdmean	2.08E-03
+	BO21_lightbotltmax_bdmean	2.08E-03
+	BO22_lightbotltmax_bdmean	2.08E-03
+	BO21_curvelmax_bdmax	1.96E-03
+	BO21_curvelmean_bdmax	1.96E-03
+	BO21_curvelmin_bdmax	1.96E-03
+	BO21_curvelrange_bdmax	1.96E-03
+	BO22_curvelmax_bdmax	1.96E-03
+	BO22_curvelmean_bdmax	1.96E-03
+	BO22_curvelmin_bdmax	1.96E-03
+	BO22_curvelrange_bdmax	1.96E-03
+	BO2_tempmean_bdmean	1.94E-03
+	BO2_temprange_bdmean	1.94E-03
+	BO21_tempmean_bdmean	1.94E-03
+	BO21_temprange_bdmean	1.94E-03
+	BO22_tempmean_bdmean	1.94E-03
+	BO22_temprange_bdmean	1.94E-03
+	BO2_silicateltmax_ss	1.90E-03
+	BO21_silicateltmax_ss	1.90E-03
+	BO22_silicateltmax_ss	1.90E-03
+	MS_biogeo01_aspect_EW_5m	1.88E-03
+	BO2_ironrange_ss	1.81E-03
+	BO21_ironrange_ss	1.81E-03
+	BO22_ironrange_ss	1.81E-03
+	BO21_curvelmin_ss	1.66E-03
+	BO22_curvelmin_ss	1.66E-03
+	BO_cloudmin	1.65E-03
+	BO22_cloudmin	1.65E-03
+	BO2_nitratemin_ss	1.65E-03
+	BO21_nitratemin_ss	1.65E-03
+	BO22_nitratemin_ss	1.65E-03
+	BO21_curvelltmax_bdmin	1.63E-03
+	BO22_curvelltmax_bdmin	1.63E-03
+	BO2_curvelrange_ss	1.58E-03
+	BO2_ironmean_ss	1.57E-03
+	BO21_ironmean_ss	1.57E-03
+	BO22_ironmean_ss	1.57E-03
+	BO21_curvelltmin_bdmin	1.55E-03
+	BO22_curvelltmin_bdmin	1.55E-03
+	BO2_lightbotmax_bdmean	1.50E-03
+	BO21_lightbotmax_bdmean	1.50E-03
+	BO22_lightbotmax_bdmean	1.50E-03
+	BO2_silicatemax_bdmax	1.49E-03
+	BO2_silicatemean_bdmax	1.49E-03
+	BO2_silicatemin_bdmax	1.49E-03
+	BO2_silicaterange_bdmax	1.49E-03
+	BO21_silicatemax_bdmax	1.49E-03
+	BO21_silicatemean_bdmax	1.49E-03
+	BO21_silicatemin_bdmax	1.49E-03
+	BO21_silicaterange_bdmax	1.49E-03
+	BO22_silicatemax_bdmax	1.49E-03
+	BO22_silicatemean_bdmax	1.49E-03
+	BO22_silicatemin_bdmax	1.49E-03
+	BO22_silicaterange_bdmax	1.49E-03
+	BO2_silicateltmax_bdmax	1.47E-03
+	BO2_silicateltmin_bdmax	1.47E-03
+	BO21_silicateltmax_bdmax	1.47E-03
+	BO21_silicateltmin_bdmax	1.47E-03
+	BO22_silicateltmax_bdmax	1.47E-03
+	BO22_silicateltmin_bdmax	1.47E-03
+	BO2_ironmax_ss	1.45E-03
+	BO21_ironmax_ss	1.45E-03
+	BO22_ironmax_ss	1.45E-03
+	BO2_ironltmax_ss	1.28E-03
+	BO21_ironltmax_ss	1.28E-03
+	BO22_ironltmax_ss	1.28E-03
+	BO2_silicatemin_bdmean	1.20E-03
+	BO21_silicatemin_bdmean	1.20E-03
+	BO22_silicatemin_bdmean	1.20E-03
+	BO2_silicatemean_ss	1.19E-03
+	BO21_silicatemean_ss	1.19E-03
+	BO22_silicatemean_ss	1.19E-03
+	BO2_silicateltmin_bdmean	1.18E-03
+	BO21_silicateltmin_bdmean	1.18E-03
+	BO22_silicateltmin_bdmean	1.18E-03
+	BO21_curvelmax_bdmin	7.97E-04
+	BO22_curvelmax_bdmin	7.97E-04
+	BO21_curvelltmin_bdmean	7.92E-04
+	BO22_curvelltmin_bdmean	7.92E-04
+	BO21_curvelltmax_bdmax	6.82E-04
+	BO21_curvelltmin_bdmax	6.82E-04
+	BO22_curvelltmax_bdmax	6.82E-04
+	BO22_curvelltmin_bdmax	6.82E-04
+	BO2_curvelmean_bdmin	5.91E-04
+	BO2_curvelmin_bdmin	5.91E-04
+	BO2_curvelrange_bdmin	5.91E-04
+	MS_biogeo07_concavity_5m	5.37E-04
+	BO2_silicatemin_ss	5.10E-04
+	BO21_silicatemin_ss	5.10E-04
+	BO22_silicatemin_ss	5.10E-04
+	BO2_curvelmax_ss	4.98E-04
+	BO21_curvelltmax_bdmean	4.85E-04
+	BO22_curvelltmax_bdmean	4.85E-04
+	BO21_curvelmean_bdmean	4.38E-04
+	BO21_curvelrange_bdmean	4.38E-04
+	BO22_curvelmean_bdmean	4.38E-04
+	BO22_curvelrange_bdmean	4.38E-04
+	BO2_curvelmin_ss	4.21E-04
+	BO2_silicateltmin_ss	3.71E-04
+	BO21_silicateltmin_ss	3.71E-04
+	BO22_silicateltmin_ss	3.71E-04
+	BO2_silicatemean_bdmean	3.45E-04
+	BO2_silicaterange_bdmean	3.45E-04
+	BO21_silicatemean_bdmean	3.45E-04
+	BO21_silicaterange_bdmean	3.45E-04
+	BO22_silicatemean_bdmean	3.45E-04
+	BO22_silicaterange_bdmean	3.45E-04
+	BO21_curvelmax_bdmean	2.17E-04
+	BO22_curvelmax_bdmean	2.17E-04
+	MS_biogeo04_profile_curvature_5m	1.00E-04
+	MS_biogeo03_plan_curvature_5m	9.16E-05
+	BO2_curvelltmax_ss	1.39E-05
<none>	0.00E+00	
+	BO2_silicateltmin_bdmin	-9.22E-05
+	BO21_silicateltmin_bdmin	-9.22E-05
+	BO22_silicateltmin_bdmin	-9.22E-05
+	BO2_silicatemax_bdmean	-1.08E-04
+	BO21_silicatemax_bdmean	-1.08E-04
+	BO22_silicatemax_bdmean	-1.08E-04
+	BO2_silicateltmax_bdmean	-1.18E-04
+	BO21_silicateltmax_bdmean	-1.18E-04
+	BO22_silicateltmax_bdmean	-1.18E-04
+	BO21_curvelmin_bdmean	-1.37E-04
+	BO22_curvelmin_bdmean	-1.37E-04
+	BO2_silicatemean_bdmin	-1.99E-04
+	BO2_silicatemin_bdmin	-1.99E-04
+	BO2_silicaterange_bdmin	-1.99E-04
+	BO21_silicatemean_bdmin	-1.99E-04
+	BO21_silicatemin_bdmin	-1.99E-04
+	BO21_silicaterange_bdmin	-1.99E-04
+	BO22_silicatemean_bdmin	-1.99E-04
+	BO22_silicatemin_bdmin	-1.99E-04
+	BO22_silicaterange_bdmin	-1.99E-04

}
m.ord_Env <- m.ord_All

save(m.ord_All, file="RDA_geo_env.gz", compress=T)
#load("RDA_geo_env.gz")

m.ord_All
{ #Results
Call: rda(formula = X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)

                Inertia Proportion Rank
Total         1.264e+03  1.000e+00
Constrained   8.868e+01  7.016e-02    3
Unconstrained 1.175e+03  9.298e-01   58
Inertia is variance

Eigenvalues for constrained axes:
 RDA1  RDA2  RDA3
43.02 24.27 21.39

Eigenvalues for unconstrained axes:
   PC1    PC2    PC3    PC4    PC5    PC6    PC7    PC8
23.713 23.292 23.067 22.905 22.788 22.545 22.364 22.323
(Showing 8 of 58 unconstrained eigenvalues)
}
attributes(m.ord_All)$VIF_remove_list
{ #Results
      Rejected_variable          Max_VIF
 [1,] "BO2_carbonphytomean_ss"   "7.3615684058146"
 [2,] "BO21_carbonphytomean_ss"  "7.3615684058146"
 [3,] "BO22_carbonphytomean_ss"  "7.3615684058146"
 [4,] "BO2_carbonphytoltmax_ss"  "8.7601191679018"
 [5,] "BO21_carbonphytoltmax_ss" "8.7601191679018"
 [6,] "BO22_carbonphytoltmax_ss" "8.7601191679018"
 [7,] "MS_sst03_5m"              "31.3768003998775"
 [8,] "MS_sst04_5m"              "24.9736737483256"
 [9,] "MS_biogeo13_sst_mean_5m"  "55.8898158309897"
      Formula
 [1,] "X ~ MS_sst09_5m + BO_nitrate"
 [2,] "X ~ MS_sst09_5m + BO_nitrate"
 [3,] "X ~ MS_sst09_5m + BO_nitrate"
 [4,] "X ~ MS_sst09_5m + BO_nitrate"
 [5,] "X ~ MS_sst09_5m + BO_nitrate"
 [6,] "X ~ MS_sst09_5m + BO_nitrate"
 [7,] "X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss"
 [8,] "X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss"
 [9,] "X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss"
}
m.ord_All$anova
{ #Results
                      R2.adj Df    AIC      F Pr(>F)
+ MS_sst09_5m       0.017669  1 443.67 2.0972  0.001 ***
+ BO_nitrate        0.020889  1 444.42 1.1973  0.001 ***
+ BO2_nitratemin_ss 0.022064  1 445.29 1.0709  0.005 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
RsquareAdj(m.ord_All)
{ #Results
$r.squared
[1] 0.07017061

$adj.r.squared
[1] 0.02207598
}
m.ord_All$CCA$eig[1:4]/m.ord_All$tot.chi
{ #Results
      RDA1       RDA2       RDA3       <NA>
0.03402604 0.01921955 0.01692501         NA
}
anova(m.ord_All, parallel=20, permutations = 9999)
{ #Results
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 9999

Model: rda(formula = X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
         Df Variance     F Pr(>F)
Model     3    88.74 1.459  1e-04 ***
Residual 58  1175.92
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
tmp_anova <- anova(m.ord_All, by="axis", parallel=20, permutations = 9999)
tmp_anova
{ #Results
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 9999

Model: rda(formula = X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
         Df Variance      F Pr(>F)
RDA1      1    43.03 2.1224 0.0001 ***
RDA2      1    24.31 1.1989 0.0001 ***
RDA3      1    21.40 1.0557 0.0172 *
Residual 58  1175.92
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
tmp_anova$Variance[1:4]/sum(tmp_anova$Variance)
{ #Results
[1] 0.03402604 0.01921955 0.01692501 0.92982939
}
anova_MARG <- anova(m.ord_All, by="margin", parallel=20)
anova_MARG
{ #Results
Permutation test for rda under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

Model: rda(formula = X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
                  Df Variance      F Pr(>F)
MS_sst09_5m        1    32.19 1.5877  0.001 ***
BO_nitrate         1    23.72 1.1701  0.005 **
BO2_nitratemin_ss  1    21.72 1.0713  0.042 *
Residual          58  1175.92
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
cbind(rownames(anova_MARG),round(anova_MARG$Variance[1:length(anova_MARG$Variance)]/sum(anova_MARG$Variance),6))
{ #Results
     [,1]                [,2]
[1,] "MS_sst09_5m"       "0.025679"
[2,] "BO_nitrate"        "0.018925"
[3,] "BO2_nitratemin_ss" "0.017327"
[4,] "Residual"          "0.938069"
}
#What are these factors?
env_data_list[env_data_list$layer_code == "MS_sst09_5m", 1:6]
{ #Results
    dataset_code  layer_code                                name
133      MARSPEC MS_sst09_5m Sea surface temperature (september)
                                                 description terrestrial marine
133 Average sea surface temperature for the month september.       FALSE   TRUE
}
env_data_list[env_data_list$layer_code == "BO_nitrate", 1:6]
{ #Results
   dataset_code layer_code    name
81   Bio-ORACLE BO_nitrate Nitrate
                                                                                                                                                                                                                                                                                                                                             description
81 This layer contains both [NO3] and [NO3+NO2] data. By this we mean chemically reactive dissolved inorganic nitrate and nitrate or nitrite. (It is important to note that data reported as [NO3] in the WOD09 should be used with caution because it is difficult to verify that the [NO3] (nitrate) data are [NO3+NO2] or [NO3]. (Boyer et al. 2009))
   terrestrial marine
81       FALSE   TRUE
}
env_data_list[env_data_list$layer_code == "BO2_nitratemin_ss", 1:6]
{ #Results
    dataset_code        layer_code                            name
421   Bio-ORACLE BO2_nitratemin_ss Nitrate concentration (minimum)
                                                 description terrestrial marine
421 Minimum mole concentration of nitrate at the sea surface       FALSE   TRUE
}
} #notepad cleanup

#Scale measurements
{```{R}```
#Pulling measurements
meas <- gen3@strata[,c(1, 12:16, 114, 29:112)]
for(i in 9:ncol(meas)){meas[ ,i] <- as.numeric(as.matrix(meas[ ,i])) }
meas[meas == 0] <- NA
rownames(meas) <- meas$INDV
#meas <- meas[, - c(86,87)]

#Making a matrix for the allometric scaling data
tmp <- t(utils::combn(colnames(meas)[9:ncol(meas)],2))
m.results <- matrix(nrow=nrow(tmp), ncol=25)
colnames(m.results) <- c("i","j","b_pvalue", "a_pvalue","Atl_b_min", "Atl_b", "Atl_b_max", "GOM_b_min", "GOM_b", "GOM_b_max", "Combine_b_min", "Combine_b", "Combine_b_max", "Atl_a_min", "Atl_a", "Atl_a_max", "GOM_a_min", "GOM_a", "GOM_a_max", "Combine_a_min", "Combine_a", "Combine_a_max", "Mean_combo_i", "Mean_Atl_i", "Mean_Gulf_i")
m.results[,1:2] <- tmp

#Loop for scaling the data
#The coefficent (a values) 

for(k in c(1:1629, 1632:3397,3399:nrow(m.results))){
	print(paste("starting", m.results[k,1], "and", m.results[k,2]))
	i <- m.results[k,1]
	j <- m.results[k,2]

		if(nrow(na.omit(meas[meas$G2 == "Atlantic",c(i,j)])) < 4){next}		#Skips if there are not enough atlantic values to use	#Probablly should just make a combine then#
		if(nrow(na.omit(meas[meas$G2 == "Gulf",c(i,j)])) < 4){next}		#Skips if there are not enough Gulf values to use
		if(nrow(meas) - sum(meas[,i] == meas[,j], na.rm=T) < 4){next}		#Skips if the two columns are too similar	#sma can't handle highly similar values
		tmp.meas <- meas[,c("G2",i,j)]
		tmp.meas <- na.omit(tmp.meas)
		tmp.err <- FALSE
		tryCatch(tmp <- sma(tmp.meas[,i] ~ tmp.meas[,j]*G2, tmp.meas, log="xy", na.action=na.omit),error=function(e){cat("Error :", conditionMessage(e)); tmp.err <- TRUE})
		if(tmp.err){next}
#		if(tmp$commoncoef$p > 0.05){tryCatch({ tmp <- sma(meas[,i] ~ meas[,j]+Region, meas, log="xy", na.action=na.omit)},error=function(e){cat("Error :", conditionMessage(e), "\n"); next})}
		m.results[k, "b_pvalue"] <- tmp$commoncoef$p
		m.results[k, c("Atl_b", "Atl_b_min", "Atl_b_max")] <- as.matrix(tmp$coef$Atlantic["slope",])
		m.results[k, c("GOM_b", "GOM_b_min", "GOM_b_max")] <- as.matrix(tmp$coef$Gulf["slope",])
		m.results[k, "Combine_b"] <- tmp$commoncoef$b
		m.results[k, c("Combine_b_min", "Combine_b_max")] <- as.matrix(t(tmp$commoncoef$ci))

		tmp <- sma(meas[ ,i] ~ meas[ ,j] + G2, meas, log="xy", type="elevation", na.action=na.omit)
		m.results[k, "Combine_a"] <- tmp$gtr$a
		m.results[k, c("Combine_a_min", "Combine_a_max")] <- as.matrix(t(tmp$gtr$ci))
		m.results[k, "Mean_combo_i"] <- mean(meas[,i], na.rm=T)

		if(tmp$commoncoef$p > 0.05) {next}
		m.results[k, "a_pvalue"] <- tmp$gtr$p
		m.results[k, c("Atl_a", "Atl_a_min", "Atl_a_max")] <- as.matrix(tmp$coef$Atlantic["elevation",])
		m.results[k, c("GOM_a", "GOM_a_min", "GOM_a_max")] <- as.matrix(tmp$coef$Gulf["elevation",])
		m.results[k, "Mean_Atl_i"] <- mean(meas[meas$G2 == "Atlantic", i], na.rm=T)
		m.results[k, "Mean_Gulf_i"] <- mean(meas[meas$G2 == "Gulf", i], na.rm=T)
}

}}}}}} #notepad cleanup

#Converting the data into a data frame with numeric values
{```{R}```
m.results.df <- data.frame(m.results)
for(i in 3:ncol(m.results.df)){m.results.df[ ,i] <- as.numeric(as.matrix(m.results.df[ ,i]))}

#Preparing matrix for converted data based upon the allometric scaling
est.meas <- matrix(nrow=nrow(m.results), ncol=(nrow(meas)+2))
est.meas[ ,1:2] <- as.matrix(m.results[ ,1:2])
colnames(est.meas) <- c("i", "j", as.character(meas[,1]))

#Loop for converting values based upon allometric scaling	#TRY3#
for(L in 1:nrow(est.meas)){
	i <- est.meas[L,1]
	j <- est.meas[L,2]
	print(paste("starting", i, "and", j))
	if(sum(is.na(m.results.df[L, "b_pvalue"])) > 0){next}
		for(k in 1:nrow(meas)){
			if(m.results.df[L, "b_pvalue"] < 0.05){
				reg <- as.character(meas[k,"G2"])
					if(reg=="Atlantic"){
					est.meas[L,as.character(meas[k,1])] <- ((mean(meas[,j], na.rm=T)/meas[k,j])^m.results.df[L,"Atl_b"])*meas[k,i]
					}
					if(reg=="Gulf"){
					est.meas[L,as.character(meas[k,1])] <- ((mean(meas[,j], na.rm=T)/meas[k,j])^m.results.df[L,"GOM_b"])*meas[k,i]
					}
				}
				else{est.meas[L,as.character(meas[k,1])] <- ((mean(meas[,j], na.rm=T)/meas[k,j])^m.results.df[L,"Combine_b"])*meas[k,i]
				}
	}
}

}}}} # notepad cleanup

# Convert to data frame with numerical columns
{```{R}```
est.meas.df <- as.data.frame(est.meas[!is.na(est.meas[,1]),])
for(i in 3:ncol(est.meas.df)){est.meas.df[,i] <- as.numeric(as.character(est.meas.df[,i]))}

#Determine the number of conversions with different slopes between basins
j=7
j <- colnames(meas)[j]
tmp.b <- m.results[m.results[,2]==j & !(is.na(m.results[,2])),]
b_scores <- tmp.b[tmp.b[,3] < 0.05, 1]
length(b_scores)

g_splits <- NULL
for(j in unique(c(m.results[,1],m.results[,2]))){
tmp.b <- m.results.df[m.results.df[,1]==j | m.results.df[,2]==j,]
b_scores <- unlist(tmp.b[tmp.b[,3] < 0.05, 1:2])
b_scores <- as.character(b_scores[which(b_scores != j)])
g_splits<- c(g_splits, length(b_scores))
}
names(g_splits) <- unique(c(m.results[,1],m.results[,2]))
#Making graphs
j=10
#tiff("ln_PCL_correlations.tif", height = 6000, width=8000, res = 200)
j <- colnames(meas)[j]
par(mfrow = c(3,4))
for(i in b_scores){print(m.results[m.results[,1]==i & m.results[,2]==j & !(is.na(m.results[,2])),1:3])
tmp <- na.omit(meas[,c("INDV", "G2", j, i)])
tmp <- tmp[ which(tmp[,3] != 0 & tmp[,4] != 0), ]
tmp$ln_X <- log(tmp[,j])
tmp$ln_Y <- log(tmp[,i])
tmp_A.lm <- lm(tmp[tmp[,"G2"] == "Atlantic", "ln_Y"] ~ tmp[tmp[,"G2"] == "Atlantic", "ln_X"])
tmp_B.lm <- lm(tmp[tmp[,"G2"] == "Gulf", "ln_Y"] ~ tmp[tmp[,"G2"] == "Gulf", "ln_X"])
plot(tmp[tmp[,"G2"] == "Gulf", c("ln_X", "ln_Y")], col = "blue", pch=19, ylim = c(min(tmp$ln_Y, na.rm=T)*0.9,max(tmp$ln_Y, na.rm=T)*1.1), xlim = c(min(tmp$ln_X)*0.9,max(tmp$ln_X)*1.1), xlab="", ylab = "")
abline(tmp_B.lm, col="blue")
par(new=T)
plot(tmp[tmp[,"G2"] == "Atlantic", c("ln_X", "ln_Y")], col = "red", pch=19, ylim = c(min(tmp$ln_Y, na.rm=T)*0.9,max(tmp$ln_Y, na.rm=T)*1.1), xlim = c(min(tmp$ln_X)*0.9,max(tmp$ln_X)*1.1), xlab = paste("ln(", j, ")", sep=""), ylab = paste("ln(", i, ")", sep=""), main = i)
abline(tmp_A.lm, col="red")
}
#dev.off()
rm(tmp, tmp_A.lm, tmp_B.lm, tmp.b)

#Putting the scaling estimates into the genind object
tmp.est <- data.frame(row.names = as.character(meas[,1]))
tmp.est <- t(est.meas.df[est.meas.df[,1] == "PCL" | est.meas.df[,2] == "PCL", 3:ncol(est.meas.df)])
colnames(tmp.est) <- as.character(unique(c(as.matrix(est.meas.df)[,1], as.matrix(est.meas.df[,2]))))[-2]
sum(indNames(gen3) == rownames(tmp.est))
#[1] 62
gen3@other$scaled.size <- as.data.frame(tmp.est)

#Looking for missing data
na.col <- apply(gen3@other$scaled.size, 2, function(x) sum(is.na(x)))
na.row <- apply(gen3@other$scaled.size, 1, function(x) sum(is.na(x)))

table(na.col)
#35 36 37 62
#67  7  4  4

table(na.row)
# 4  5  6  7  8 82
#19  4  2  1  1 35

#Removing measurements and individuals with no measurement values
rm.meas <- names(na.col)[na.col == 62]
rm.indv <- names(na.row)[na.row == 82]

tmp.meas <- gen3@other$scaled.size[!rownames(gen3@other$scaled.size) %in% rm.indv, !colnames(gen3@other$scaled.size) %in% rm.meas]

#Measurements with 0 values
na.col <- apply(tmp.meas, 2, function(x) sum(is.na(x)))
na.row <- apply(tmp.meas, 1, function(x) sum(is.na(x)))

table(na.col)
# 0  1  2
#67  7  4

table(na.row)
# 0  1  2  3  4
#19  4  2  1  1

#Filling in NA values from density plots
tmp.rank <- apply(tmp.meas, 2, rank)
tmp.rank[is.na(tmp.meas)] <- NA
med.rank <- rank(apply(tmp.rank,1,function(x) median(x, na.rm=T)))

for(i in which(na.col > 0)){
INDV <- which(is.na(tmp.meas[,i]))
meas.val <- quantile(tmp.meas[- INDV,i], probs = med.rank[INDV]/nrow(tmp.rank))
tmp.meas[INDV,i] <- meas.val
}

na.col <- apply(tmp.meas, 2, function(x) sum(is.na(x)))
na.row <- apply(tmp.meas, 1, function(x) sum(is.na(x)))

table(na.col)
# 0
#78

table(na.row)
# 0
#27
} #notepad cleanup

#Scaled measurements RDA
{```{R}```
set.indv <- rownames(tmp.meas)
gen.tmp <- gen3[set.indv,]
gen.tmp@other$scaled.size <- as.data.frame(tmp.meas)

DIST_df <- as.data.frame(DIST_m)

RDA_data <- cbind(env_dat_filter, DIST_df[match(rownames(env_dat_filter),DIST_df[,1]),2])
colnames(RDA_data)[ncol(RDA_data)] <- "Geo_dist"
RDA_data$Geo_dist <- as.numeric(as.matrix(RDA_data$Geo_dist))
RDA_data <- RDA_data[set.indv, ]
RDA_data <- cbind(RDA_data, tmp.meas[match(rownames(RDA_data),rownames(tmp.meas))])

RDA_data_scale <- as.data.frame(scale(RDA_data[match(indNames(gen.tmp), rownames(RDA_data)),]))

X <- scaleGen(gen.tmp, NA.method="mean", scale=F)
m1<-rda(X ~ ., RDA_data_scale)
m0<-rda(X ~ 1, RDA_data_scale)
set.seed(1235)
m.ord_All <- ordiR2step(m0, scope = formula(m1), Pin=0.05, Pout=0.1, permutations = 999, parallel=10, R2scope=F)
{ #Results
                                     R2.adjusted
+ Geo_dist                          1.330049e-02
+ MS_sst09_5m                       1.316670e-02
+ Lat                               1.293012e-02
+ MS_sst10_5m                       1.292323e-02
+ MS_sst06_5m                       1.291236e-02
+ BO_chlomin                        1.268854e-02
+ MS_sst11_5m                       1.265417e-02
+ BO_damin                          1.258135e-02
+ BO22_damin                        1.258135e-02
+ Lon                               1.237149e-02
+ BO2_tempmax_ss                    1.233704e-02
+ BO21_tempmax_ss                   1.233704e-02
+ BO22_tempmax_ss                   1.233704e-02
+ MS_sst08_5m                       1.227517e-02
+ MS_biogeo15_sst_max_5m            1.227121e-02
+ BO2_ppmean_bdmean                 1.224406e-02
+ BO2_pprange_bdmean                1.224406e-02
+ BO21_ppmean_bdmean                1.224406e-02
+ BO21_pprange_bdmean               1.224406e-02
+ BO22_ppmean_bdmean                1.224406e-02
+ BO22_pprange_bdmean               1.224406e-02
+ BO2_ppmax_bdmean                  1.219809e-02
+ BO21_ppmax_bdmean                 1.219809e-02
+ BO22_ppmax_bdmean                 1.219809e-02
+ MS_biogeo13_sst_mean_5m           1.218824e-02
+ BO_sstmean                        1.215599e-02
+ BO2_ppltmax_bdmean                1.213148e-02
+ BO21_ppltmax_bdmean               1.213148e-02
+ BO22_ppltmax_bdmean               1.213148e-02
+ BO_parmean                        1.212568e-02
+ BO22_parmean                      1.212568e-02
+ BO2_ppltmin_bdmean                1.211247e-02
+ BO21_ppltmin_bdmean               1.211247e-02
+ BO22_ppltmin_bdmean               1.211247e-02
+ BO2_ppltmax_bdmax                 1.209159e-02
+ BO2_ppltmin_bdmax                 1.209159e-02
+ BO21_ppltmax_bdmax                1.209159e-02
+ BO21_ppltmin_bdmax                1.209159e-02
+ BO22_ppltmax_bdmax                1.209159e-02
+ BO22_ppltmin_bdmax                1.209159e-02
+ BO2_ppmax_bdmax                   1.204173e-02
+ BO2_ppmean_bdmax                  1.204173e-02
+ BO2_ppmin_bdmax                   1.204173e-02
+ BO2_pprange_bdmax                 1.204173e-02
+ BO21_ppmax_bdmax                  1.204173e-02
+ BO21_ppmean_bdmax                 1.204173e-02
+ BO21_ppmin_bdmax                  1.204173e-02
+ BO21_pprange_bdmax                1.204173e-02
+ BO22_ppmax_bdmax                  1.204173e-02
+ BO22_ppmean_bdmax                 1.204173e-02
+ BO22_ppmin_bdmax                  1.204173e-02
+ BO22_pprange_bdmax                1.204173e-02
+ BO_silicate                       1.201705e-02
+ MS_sst05_5m                       1.187602e-02
+ BO2_carbonphytomax_bdmean         1.184101e-02
+ BO21_carbonphytomax_bdmean        1.184101e-02
+ BO22_carbonphytomax_bdmean        1.184101e-02
+ BO2_carbonphytoltmax_bdmean       1.179730e-02
+ BO21_carbonphytoltmax_bdmean      1.179730e-02
+ BO22_carbonphytoltmax_bdmean      1.179730e-02
+ BO_sstmax                         1.176743e-02
+ BO_sstmin                         1.174910e-02
+ BO22_chlomax_bdmax                1.171580e-02
+ BO22_chlomean_bdmax               1.171580e-02
+ BO22_chlomin_bdmax                1.171580e-02
+ BO22_chlorange_bdmax              1.171580e-02
+ BO22_chlomax_bdmean               1.171035e-02
+ BO2_carbonphytomax_bdmax          1.168528e-02
+ BO2_carbonphytomean_bdmax         1.168528e-02
+ BO2_carbonphytomin_bdmax          1.168528e-02
+ BO2_carbonphytorange_bdmax        1.168528e-02
+ BO21_carbonphytomax_bdmax         1.168528e-02
+ BO21_carbonphytomean_bdmax        1.168528e-02
+ BO21_carbonphytomin_bdmax         1.168528e-02
+ BO21_carbonphytorange_bdmax       1.168528e-02
+ BO22_carbonphytomax_bdmax         1.168528e-02
+ BO22_carbonphytomean_bdmax        1.168528e-02
+ BO22_carbonphytomin_bdmax         1.168528e-02
+ BO22_carbonphytorange_bdmax       1.168528e-02
+ MS_biogeo14_sst_min_5m            1.167010e-02
+ BO2_ppmin_bdmean                  1.166146e-02
+ BO21_ppmin_bdmean                 1.166146e-02
+ BO22_ppmin_bdmean                 1.166146e-02
+ BO2_chloltmax_ss                  1.165011e-02
+ BO22_chloltmax_ss                 1.165011e-02
+ BO22_chloltmax_bdmax              1.164288e-02
+ BO22_chloltmin_bdmax              1.164288e-02
+ BO2_chlomax_ss                    1.163168e-02
+ BO22_chlomax_ss                   1.163168e-02
+ MS_sst07_5m                       1.159147e-02
+ BO2_phosphatemean_ss              1.158859e-02
+ BO21_phosphatemean_ss             1.158859e-02
+ BO22_phosphatemean_ss             1.158859e-02
+ BO_damean                         1.158778e-02
+ BO22_damean                       1.158778e-02
+ BO2_carbonphytoltmax_bdmax        1.157653e-02
+ BO2_carbonphytoltmin_bdmax        1.157653e-02
+ BO21_carbonphytoltmax_bdmax       1.157653e-02
+ BO21_carbonphytoltmin_bdmax       1.157653e-02
+ BO22_carbonphytoltmax_bdmax       1.157653e-02
+ BO22_carbonphytoltmin_bdmax       1.157653e-02
+ MS_sst03_5m                       1.157597e-02
+ BO2_lightbotmin_bdmean            1.153853e-02
+ BO21_lightbotmin_bdmean           1.153853e-02
+ BO22_lightbotmin_bdmean           1.153853e-02
+ BO22_chloltmax_bdmean             1.149911e-02
+ BO2_carbonphytoltmax_bdmin        1.146213e-02
+ BO21_carbonphytoltmax_bdmin       1.146213e-02
+ BO22_carbonphytoltmax_bdmin       1.146213e-02
+ BO2_carbonphytomax_bdmin          1.144587e-02
+ BO21_carbonphytomax_bdmin         1.144587e-02
+ BO22_carbonphytomax_bdmin         1.144587e-02
+ BO2_templtmax_ss                  1.142109e-02
+ BO21_templtmax_ss                 1.142109e-02
+ BO22_templtmax_ss                 1.142109e-02
+ BO2_chlorange_ss                  1.140956e-02
+ MS_sst02_5m                       1.140620e-02
+ BO2_tempmin_ss                    1.134845e-02
+ BO21_tempmin_ss                   1.134845e-02
+ BO22_tempmin_ss                   1.134845e-02
+ MS_sst01_5m                       1.129277e-02
+ BO2_ppmax_bdmin                   1.128644e-02
+ BO21_ppmax_bdmin                  1.128644e-02
+ BO22_ppmax_bdmin                  1.128644e-02
+ BO2_ironmax_bdmean                1.125669e-02
+ BO21_ironmax_bdmean               1.125669e-02
+ BO22_ironmax_bdmean               1.125669e-02
+ BO2_lightbotltmin_bdmean          1.120706e-02
+ BO21_lightbotltmin_bdmean         1.120706e-02
+ BO22_lightbotltmin_bdmean         1.120706e-02
+ BO2_ironltmax_bdmean              1.119225e-02
+ BO21_ironltmax_bdmean             1.119225e-02
+ BO22_ironltmax_bdmean             1.119225e-02
+ BO2_ppltmin_bdmin                 1.115572e-02
+ BO21_ppltmin_bdmin                1.115572e-02
+ BO22_ppltmin_bdmin                1.115572e-02
+ BO_calcite                        1.115211e-02
+ BO22_calcite                      1.115211e-02
+ BO2_tempmean_ss                   1.113096e-02
+ BO21_tempmean_ss                  1.113096e-02
+ BO22_tempmean_ss                  1.113096e-02
+ BO22_chlomax_bdmin                1.110773e-02
+ BO2_carbonphytomean_bdmean        1.105597e-02
+ BO2_carbonphytorange_bdmean       1.105597e-02
+ BO21_carbonphytomean_bdmean       1.105597e-02
+ BO21_carbonphytorange_bdmean      1.105597e-02
+ BO22_carbonphytomean_bdmean       1.105597e-02
+ BO22_carbonphytorange_bdmean      1.105597e-02
+ BO2_ppltmax_bdmin                 1.102291e-02
+ BO21_ppltmax_bdmin                1.102291e-02
+ BO22_ppltmax_bdmin                1.102291e-02
+ BO22_chloltmax_bdmin              1.096067e-02
+ BO2_carbonphytomin_ss             1.074443e-02
+ BO21_carbonphytomin_ss            1.074443e-02
+ BO22_carbonphytomin_ss            1.074443e-02
+ BO2_carbonphytomin_bdmean         1.073776e-02
+ BO21_carbonphytomin_bdmean        1.073776e-02
+ BO22_carbonphytomin_bdmean        1.073776e-02
+ BO2_dissoxmax_bdmax               1.070300e-02
+ BO2_dissoxmax_bdmean              1.070300e-02
+ BO2_dissoxmax_bdmin               1.070300e-02
+ BO2_dissoxmean_bdmax              1.070300e-02
+ BO2_dissoxmin_bdmax               1.070300e-02
+ BO2_dissoxrange_bdmax             1.070300e-02
+ BO2_dissoxmax_ss                  1.070300e-02
+ BO21_dissoxmax_bdmax              1.070300e-02
+ BO21_dissoxmax_bdmean             1.070300e-02
+ BO21_dissoxmax_bdmin              1.070300e-02
+ BO21_dissoxmax_ss                 1.070300e-02
+ BO21_dissoxmean_bdmax             1.070300e-02
+ BO21_dissoxmin_bdmax              1.070300e-02
+ BO21_dissoxrange_bdmax            1.070300e-02
+ BO22_dissoxmax_bdmax              1.070300e-02
+ BO22_dissoxmax_bdmean             1.070300e-02
+ BO22_dissoxmax_bdmin              1.070300e-02
+ BO22_dissoxmax_ss                 1.070300e-02
+ BO22_dissoxmean_bdmax             1.070300e-02
+ BO22_dissoxmin_bdmax              1.070300e-02
+ BO22_dissoxrange_bdmax            1.070300e-02
+ BO2_lightbotmean_bdmean           1.068031e-02
+ BO2_lightbotrange_bdmean          1.068031e-02
+ BO21_lightbotmean_bdmean          1.068031e-02
+ BO21_lightbotrange_bdmean         1.068031e-02
+ BO22_lightbotmean_bdmean          1.068031e-02
+ BO22_lightbotrange_bdmean         1.068031e-02
+ BO2_chlomean_ss                   1.062860e-02
+ BO22_chlomean_ss                  1.062860e-02
+ BO_sstrange                       1.056141e-02
+ MS_sst04_5m                       1.055793e-02
+ BO2_chlomin_ss                    1.053257e-02
+ BO22_chlomin_ss                   1.053257e-02
+ BO2_dissoxmean_bdmin              1.042536e-02
+ BO2_dissoxmin_bdmean              1.042536e-02
+ BO2_dissoxmin_bdmin               1.042536e-02
+ BO2_dissoxrange_bdmin             1.042536e-02
+ BO2_dissoxmin_ss                  1.042536e-02
+ BO21_dissoxmean_bdmin             1.042536e-02
+ BO21_dissoxmin_bdmean             1.042536e-02
+ BO21_dissoxmin_bdmin              1.042536e-02
+ BO21_dissoxmin_ss                 1.042536e-02
+ BO21_dissoxrange_bdmin            1.042536e-02
+ BO22_dissoxmean_bdmin             1.042536e-02
+ BO22_dissoxmin_bdmean             1.042536e-02
+ BO22_dissoxmin_bdmin              1.042536e-02
+ BO22_dissoxmin_ss                 1.042536e-02
+ BO22_dissoxrange_bdmin            1.042536e-02
+ BO2_carbonphytoltmin_bdmean       1.033928e-02
+ BO21_carbonphytoltmin_bdmean      1.033928e-02
+ BO22_carbonphytoltmin_bdmean      1.033928e-02
+ MS_biogeo16_sst_range_5m          1.032229e-02
+ BO2_dissoxrange_ss                1.020540e-02
+ BO21_dissoxrange_ss               1.020540e-02
+ BO22_dissoxrange_ss               1.020540e-02
+ BO2_lightbotmean_bdmin            1.015061e-02
+ BO2_lightbotmin_bdmin             1.015061e-02
+ BO2_lightbotrange_bdmin           1.015061e-02
+ BO21_lightbotmean_bdmin           1.015061e-02
+ BO21_lightbotmin_bdmin            1.015061e-02
+ BO21_lightbotrange_bdmin          1.015061e-02
+ BO22_lightbotmean_bdmin           1.015061e-02
+ BO22_lightbotmin_bdmin            1.015061e-02
+ BO22_lightbotrange_bdmin          1.015061e-02
+ BO2_carbonphytoltmin_ss           1.008410e-02
+ BO21_carbonphytoltmin_ss          1.008410e-02
+ BO22_carbonphytoltmin_ss          1.008410e-02
+ BO2_chloltmin_ss                  1.006596e-02
+ BO22_chloltmin_ss                 1.006596e-02
+ BO2_temprange_ss                  9.988272e-03
+ BO21_temprange_ss                 9.988272e-03
+ BO22_temprange_ss                 9.988272e-03
+ BO2_lightbotltmax_bdmean          9.951589e-03
+ BO21_lightbotltmax_bdmean         9.951589e-03
+ BO22_lightbotltmax_bdmean         9.951589e-03
+ BO_chlomean                       9.931654e-03
+ BO2_ppmean_bdmin                  9.926599e-03
+ BO2_ppmin_bdmin                   9.926599e-03
+ BO2_pprange_bdmin                 9.926599e-03
+ BO21_ppmean_bdmin                 9.926599e-03
+ BO21_ppmin_bdmin                  9.926599e-03
+ BO21_pprange_bdmin                9.926599e-03
+ BO22_ppmean_bdmin                 9.926599e-03
+ BO22_ppmin_bdmin                  9.926599e-03
+ BO22_pprange_bdmin                9.926599e-03
+ BO2_pprange_ss                    9.911131e-03
+ BO21_pprange_ss                   9.911131e-03
+ BO22_pprange_ss                   9.911131e-03
+ BO22_chlomean_bdmean              9.905495e-03
+ BO22_chlorange_bdmean             9.905495e-03
+ BO2_salinityltmin_bdmin           9.895821e-03
+ BO21_salinityltmin_bdmin          9.895821e-03
+ BO22_salinityltmin_bdmin          9.895821e-03
+ BO2_ppmax_ss                      9.867094e-03
+ BO21_ppmax_ss                     9.867094e-03
+ BO22_ppmax_ss                     9.867094e-03
+ BO2_dissoxltmin_bdmean            9.818194e-03
+ BO2_dissoxltmin_bdmin             9.818194e-03
+ BO2_dissoxltmin_ss                9.818194e-03
+ BO21_dissoxltmin_bdmean           9.818194e-03
+ BO21_dissoxltmin_bdmin            9.818194e-03
+ BO21_dissoxltmin_ss               9.818194e-03
+ BO22_dissoxltmin_bdmean           9.818194e-03
+ BO22_dissoxltmin_bdmin            9.818194e-03
+ BO22_dissoxltmin_ss               9.818194e-03
+ BO2_carbonphytoltmax_ss           9.760284e-03
+ BO21_carbonphytoltmax_ss          9.760284e-03
+ BO22_carbonphytoltmax_ss          9.760284e-03
+ BO2_ironmean_bdmean               9.736745e-03
+ BO2_ironrange_bdmean              9.736745e-03
+ BO21_ironmean_bdmean              9.736745e-03
+ BO21_ironrange_bdmean             9.736745e-03
+ BO22_ironmean_bdmean              9.736745e-03
+ BO22_ironrange_bdmean             9.736745e-03
+ BO2_lightbotltmin_bdmin           9.703660e-03
+ BO21_lightbotltmin_bdmin          9.703660e-03
+ BO22_lightbotltmin_bdmin          9.703660e-03
+ BO_nitrate                        9.703193e-03
+ BO2_salinityltmax_bdmin           9.647313e-03
+ BO21_salinityltmax_bdmin          9.647313e-03
+ BO22_salinityltmax_bdmin          9.647313e-03
+ BO2_salinitymean_bdmin            9.636117e-03
+ BO2_salinitymin_bdmin             9.636117e-03
+ BO2_salinityrange_bdmin           9.636117e-03
+ BO21_salinitymean_bdmin           9.636117e-03
+ BO21_salinitymin_bdmin            9.636117e-03
+ BO21_salinityrange_bdmin          9.636117e-03
+ BO22_salinitymean_bdmin           9.636117e-03
+ BO22_salinitymin_bdmin            9.636117e-03
+ BO22_salinityrange_bdmin          9.636117e-03
+ BO2_templtmin_ss                  9.584630e-03
+ BO21_templtmin_ss                 9.584630e-03
+ BO22_templtmin_ss                 9.584630e-03
+ BO_ph                             9.548917e-03
+ BO22_ph                           9.548917e-03
+ BO2_curvelltmax_bdmax             9.323049e-03
+ BO2_curvelltmin_bdmax             9.323049e-03
+ BO2_carbonphytomax_ss             9.317949e-03
+ BO21_carbonphytomax_ss            9.317949e-03
+ BO22_carbonphytomax_ss            9.317949e-03
+ BO2_lightbotltmax_bdmax           9.288606e-03
+ BO2_lightbotltmin_bdmax           9.288606e-03
+ BO21_lightbotltmax_bdmax          9.288606e-03
+ BO21_lightbotltmin_bdmax          9.288606e-03
+ BO22_lightbotltmax_bdmax          9.288606e-03
+ BO22_lightbotltmin_bdmax          9.288606e-03
+ BO2_salinitymax_bdmin             9.183647e-03
+ BO21_salinitymax_bdmin            9.183647e-03
+ BO22_salinitymax_bdmin            9.183647e-03
+ BO2_carbonphytomean_ss            9.121860e-03
+ BO21_carbonphytomean_ss           9.121860e-03
+ BO22_carbonphytomean_ss           9.121860e-03
+ BO2_ironmax_bdmin                 9.118550e-03
+ BO21_ironmax_bdmin                9.118550e-03
+ BO22_ironmax_bdmin                9.118550e-03
+ BO2_phosphateltmax_ss             8.852648e-03
+ BO21_phosphateltmax_ss            8.852648e-03
+ BO22_phosphateltmax_ss            8.852648e-03
+ BO2_curvelmax_bdmax               8.814542e-03
+ BO2_curvelmean_bdmax              8.814542e-03
+ BO2_curvelmin_bdmax               8.814542e-03
+ BO2_curvelrange_bdmax             8.814542e-03
+ BO_dissox                         8.805378e-03
+ MS_biogeo11_sss_range_5m          8.764432e-03
+ BO2_ironltmax_bdmin               8.763136e-03
+ BO21_ironltmax_bdmin              8.763136e-03
+ BO22_ironltmax_bdmin              8.763136e-03
+ BO2_carbonphytorange_ss           8.688107e-03
+ BO21_carbonphytorange_ss          8.688107e-03
+ BO22_carbonphytorange_ss          8.688107e-03
+ BO2_curvelltmin_bdmean            8.679164e-03
+ BO2_ppltmax_ss                    8.595024e-03
+ BO21_ppltmax_ss                   8.595024e-03
+ BO22_ppltmax_ss                   8.595024e-03
+ BO2_carbonphytomean_bdmin         8.505678e-03
+ BO2_carbonphytomin_bdmin          8.505678e-03
+ BO2_carbonphytorange_bdmin        8.505678e-03
+ BO21_carbonphytomean_bdmin        8.505678e-03
+ BO21_carbonphytomin_bdmin         8.505678e-03
+ BO21_carbonphytorange_bdmin       8.505678e-03
+ BO22_carbonphytomean_bdmin        8.505678e-03
+ BO22_carbonphytomin_bdmin         8.505678e-03
+ BO22_carbonphytorange_bdmin       8.505678e-03
+ MS_sss11_5m                       8.469547e-03
+ BO2_carbonphytoltmin_bdmin        8.195287e-03
+ BO21_carbonphytoltmin_bdmin       8.195287e-03
+ BO22_carbonphytoltmin_bdmin       8.195287e-03
+ BO2_phosphateltmin_ss             8.190572e-03
+ BO21_phosphateltmin_ss            8.190572e-03
+ BO22_phosphateltmin_ss            8.190572e-03
+ BO_salinity                       8.091392e-03
+ BO_damax                          8.005907e-03
+ BO22_damax                        8.005907e-03
+ BO2_ironltmin_bdmin               7.924984e-03
+ BO21_ironltmin_bdmin              7.924984e-03
+ BO22_ironltmin_bdmin              7.924984e-03
+ MS_sss02_5m                       7.883844e-03
+ BO_phosphate                      7.875638e-03
+ BO2_curvelmean_bdmean             7.850093e-03
+ BO2_curvelrange_bdmean            7.850093e-03
+ BO_parmax                         7.790042e-03
+ BO22_parmax                       7.790042e-03
+ MS_biogeo12_sss_variance_5m       7.787819e-03
+ BO2_ppmean_ss                     7.740698e-03
+ BO21_ppmean_ss                    7.740698e-03
+ BO22_ppmean_ss                    7.740698e-03
+ BO2_curvelltmax_bdmean            7.656159e-03
+ MS_biogeo17_sst_variance_5m       7.628631e-03
+ BO2_dissoxmean_bdmean             7.567293e-03
+ BO2_dissoxrange_bdmean            7.567293e-03
+ BO2_dissoxmean_ss                 7.567293e-03
+ BO21_dissoxmean_bdmean            7.567293e-03
+ BO21_dissoxmean_ss                7.567293e-03
+ BO21_dissoxrange_bdmean           7.567293e-03
+ BO22_dissoxmean_bdmean            7.567293e-03
+ BO22_dissoxmean_ss                7.567293e-03
+ BO22_dissoxrange_bdmean           7.567293e-03
+ BO2_ironmax_bdmax                 7.566523e-03
+ BO2_ironmean_bdmax                7.566523e-03
+ BO2_ironmin_bdmax                 7.566523e-03
+ BO2_ironrange_bdmax               7.566523e-03
+ BO21_ironmax_bdmax                7.566523e-03
+ BO21_ironmean_bdmax               7.566523e-03
+ BO21_ironmin_bdmax                7.566523e-03
+ BO21_ironrange_bdmax              7.566523e-03
+ BO22_ironmax_bdmax                7.566523e-03
+ BO22_ironmean_bdmax               7.566523e-03
+ BO22_ironmin_bdmax                7.566523e-03
+ BO22_ironrange_bdmax              7.566523e-03
+ MS_biogeo10_sss_max_5m            7.424605e-03
+ BO2_curvelmean_ss                 7.383193e-03
+ BO2_ironltmax_bdmax               7.360692e-03
+ BO2_ironltmin_bdmax               7.360692e-03
+ BO21_ironltmax_bdmax              7.360692e-03
+ BO21_ironltmin_bdmax              7.360692e-03
+ BO22_ironltmax_bdmax              7.360692e-03
+ BO22_ironltmin_bdmax              7.360692e-03
+ BO2_ppltmin_ss                    7.104266e-03
+ BO21_ppltmin_ss                   7.104266e-03
+ BO22_ppltmin_ss                   7.104266e-03
+ BO2_lightbotltmax_bdmin           7.091680e-03
+ BO21_lightbotltmax_bdmin          7.091680e-03
+ BO22_lightbotltmax_bdmin          7.091680e-03
+ BO22_chloltmin_bdmean             6.859314e-03
+ BO2_dissoxltmax_bdmax             6.851720e-03
+ BO2_dissoxltmax_bdmean            6.851720e-03
+ BO2_dissoxltmax_bdmin             6.851720e-03
+ BO2_dissoxltmin_bdmax             6.851720e-03
+ BO2_dissoxltmax_ss                6.851720e-03
+ BO21_dissoxltmax_bdmax            6.851720e-03
+ BO21_dissoxltmax_bdmean           6.851720e-03
+ BO21_dissoxltmax_bdmin            6.851720e-03
+ BO21_dissoxltmax_ss               6.851720e-03
+ BO21_dissoxltmin_bdmax            6.851720e-03
+ BO22_dissoxltmax_bdmax            6.851720e-03
+ BO22_dissoxltmax_bdmean           6.851720e-03
+ BO22_dissoxltmax_bdmin            6.851720e-03
+ BO22_dissoxltmax_ss               6.851720e-03
+ BO22_dissoxltmin_bdmax            6.851720e-03
+ BO2_curvelmax_bdmean              6.819288e-03
+ BO22_chlomin_bdmean               6.756463e-03
+ MS_biogeo05_dist_shore_5m         6.740734e-03
+ BO2_nitratemax_bdmax              6.685816e-03
+ BO2_nitratemean_bdmax             6.685816e-03
+ BO2_nitratemin_bdmax              6.685816e-03
+ BO2_nitraterange_bdmax            6.685816e-03
+ BO21_nitratemax_bdmax             6.685816e-03
+ BO21_nitratemean_bdmax            6.685816e-03
+ BO21_nitratemin_bdmax             6.685816e-03
+ BO21_nitraterange_bdmax           6.685816e-03
+ BO22_nitratemax_bdmax             6.685816e-03
+ BO22_nitratemean_bdmax            6.685816e-03
+ BO22_nitratemin_bdmax             6.685816e-03
+ BO22_nitraterange_bdmax           6.685816e-03
+ MS_sss09_5m                       6.480763e-03
+ BO2_phosphatemax_bdmax            6.330168e-03
+ BO2_phosphatemean_bdmax           6.330168e-03
+ BO2_phosphatemin_bdmax            6.330168e-03
+ BO2_phosphaterange_bdmax          6.330168e-03
+ BO21_phosphatemax_bdmax           6.330168e-03
+ BO21_phosphatemean_bdmax          6.330168e-03
+ BO21_phosphatemin_bdmax           6.330168e-03
+ BO21_phosphaterange_bdmax         6.330168e-03
+ BO22_phosphatemax_bdmax           6.330168e-03
+ BO22_phosphatemean_bdmax          6.330168e-03
+ BO22_phosphatemin_bdmax           6.330168e-03
+ BO22_phosphaterange_bdmax         6.330168e-03
+ BO2_ironmean_bdmin                6.164595e-03
+ BO2_ironmin_bdmin                 6.164595e-03
+ BO2_ironrange_bdmin               6.164595e-03
+ BO21_ironmean_bdmin               6.164595e-03
+ BO21_ironmin_bdmin                6.164595e-03
+ BO21_ironrange_bdmin              6.164595e-03
+ BO22_ironmean_bdmin               6.164595e-03
+ BO22_ironmin_bdmin                6.164595e-03
+ BO22_ironrange_bdmin              6.164595e-03
+ BO2_phosphatemax_ss               6.093764e-03
+ BO21_phosphatemax_ss              6.093764e-03
+ BO22_phosphatemax_ss              6.093764e-03
+ BO2_phosphaterange_ss             6.072906e-03
+ BO21_phosphaterange_ss            6.072906e-03
+ BO22_phosphaterange_ss            6.072906e-03
+ BO2_nitrateltmax_bdmax            6.017129e-03
+ BO2_nitrateltmin_bdmax            6.017129e-03
+ BO21_nitrateltmax_bdmax           6.017129e-03
+ BO21_nitrateltmin_bdmax           6.017129e-03
+ BO22_nitrateltmax_bdmax           6.017129e-03
+ BO22_nitrateltmin_bdmax           6.017129e-03
+ BO2_phosphateltmax_bdmax          5.822542e-03
+ BO2_phosphateltmin_bdmax          5.822542e-03
+ BO21_phosphateltmax_bdmax         5.822542e-03
+ BO21_phosphateltmin_bdmax         5.822542e-03
+ BO22_phosphateltmax_bdmax         5.822542e-03
+ BO22_phosphateltmin_bdmax         5.822542e-03
+ BO2_curvelltmax_bdmin             5.524075e-03
+ BO2_ironltmin_bdmean              5.516451e-03
+ BO21_ironltmin_bdmean             5.516451e-03
+ BO22_ironltmin_bdmean             5.516451e-03
+ BO2_lightbotmax_bdmax             5.446115e-03
+ BO2_lightbotmean_bdmax            5.446115e-03
+ BO2_lightbotmin_bdmax             5.446115e-03
+ BO2_lightbotrange_bdmax           5.446115e-03
+ BO21_lightbotmax_bdmax            5.446115e-03
+ BO21_lightbotmean_bdmax           5.446115e-03
+ BO21_lightbotmin_bdmax            5.446115e-03
+ BO21_lightbotrange_bdmax          5.446115e-03
+ BO22_lightbotmax_bdmax            5.446115e-03
+ BO22_lightbotmean_bdmax           5.446115e-03
+ BO22_lightbotmin_bdmax            5.446115e-03
+ BO22_lightbotrange_bdmax          5.446115e-03
+ BO2_lightbotmax_bdmean            5.329859e-03
+ BO21_lightbotmax_bdmean           5.329859e-03
+ BO22_lightbotmax_bdmean           5.329859e-03
+ MS_sss10_5m                       5.312932e-03
+ BO2_templtmax_bdmax               5.205509e-03
+ BO2_templtmin_bdmax               5.205509e-03
+ BO21_templtmax_bdmax              5.205509e-03
+ BO21_templtmin_bdmax              5.205509e-03
+ BO22_templtmax_bdmax              5.205509e-03
+ BO22_templtmin_bdmax              5.205509e-03
+ BO2_curvelmax_bdmin               5.183772e-03
+ BO2_tempmax_bdmax                 5.133291e-03
+ BO2_tempmean_bdmax                5.133291e-03
+ BO2_tempmin_bdmax                 5.133291e-03
+ BO2_temprange_bdmax               5.133291e-03
+ BO21_tempmax_bdmax                5.133291e-03
+ BO21_tempmean_bdmax               5.133291e-03
+ BO21_tempmin_bdmax                5.133291e-03
+ BO21_temprange_bdmax              5.133291e-03
+ BO22_tempmax_bdmax                5.133291e-03
+ BO22_tempmean_bdmax               5.133291e-03
+ BO22_tempmin_bdmax                5.133291e-03
+ BO22_temprange_bdmax              5.133291e-03
+ BO2_nitratemax_bdmean             4.894053e-03
+ BO21_nitratemax_bdmean            4.894053e-03
+ BO22_nitratemax_bdmean            4.894053e-03
+ BO2_ppmin_ss                      4.823715e-03
+ BO21_ppmin_ss                     4.823715e-03
+ BO22_ppmin_ss                     4.823715e-03
+ MS_biogeo08_sss_mean_5m           4.716971e-03
+ BO2_phosphatemax_bdmean           4.696205e-03
+ BO21_phosphatemax_bdmean          4.696205e-03
+ BO22_phosphatemax_bdmean          4.696205e-03
+ BO_cloudmean                      4.653165e-03
+ BO22_cloudmean                    4.653165e-03
+ BO_cloudmax                       4.605666e-03
+ BO22_cloudmax                     4.605666e-03
+ BO2_tempmax_bdmean                4.589377e-03
+ BO21_tempmax_bdmean               4.589377e-03
+ BO22_tempmax_bdmean               4.589377e-03
+ MS_sss03_5m                       4.510331e-03
+ BO2_templtmax_bdmean              4.501385e-03
+ BO21_templtmax_bdmean             4.501385e-03
+ BO22_templtmax_bdmean             4.501385e-03
+ BO2_curvelltmin_bdmin             4.448117e-03
+ BO2_phosphatemin_ss               4.408080e-03
+ BO21_phosphatemin_ss              4.408080e-03
+ BO22_phosphatemin_ss              4.408080e-03
+ BO21_curvelmean_bdmin             4.232078e-03
+ BO21_curvelmin_bdmin              4.232078e-03
+ BO21_curvelrange_bdmin            4.232078e-03
+ BO22_curvelmean_bdmin             4.232078e-03
+ BO22_curvelmin_bdmin              4.232078e-03
+ BO22_curvelrange_bdmin            4.232078e-03
+ BO2_nitrateltmax_bdmean           4.158901e-03
+ BO21_nitrateltmax_bdmean          4.158901e-03
+ BO22_nitrateltmax_bdmean          4.158901e-03
+ BO21_curvelltmin_ss               4.097441e-03
+ BO22_curvelltmin_ss               4.097441e-03
+ BO2_phosphateltmax_bdmean         4.026644e-03
+ BO21_phosphateltmax_bdmean        4.026644e-03
+ BO22_phosphateltmax_bdmean        4.026644e-03
+ BO2_lightbotmax_bdmin             4.006678e-03
+ BO21_lightbotmax_bdmin            4.006678e-03
+ BO22_lightbotmax_bdmin            4.006678e-03
+ BO2_salinitymin_bdmean            3.895560e-03
+ BO21_salinitymin_bdmean           3.895560e-03
+ BO22_salinitymin_bdmean           3.895560e-03
+ BO_bathymax                       3.840640e-03
+ HDW                               3.834029e-03
+ MS_sss12_5m                       3.829136e-03
+ BO2_tempmean_bdmin                3.771373e-03
+ BO2_tempmin_bdmin                 3.771373e-03
+ BO2_temprange_bdmin               3.771373e-03
+ BO21_tempmean_bdmin               3.771373e-03
+ BO21_tempmin_bdmin                3.771373e-03
+ BO21_temprange_bdmin              3.771373e-03
+ BO22_tempmean_bdmin               3.771373e-03
+ BO22_tempmin_bdmin                3.771373e-03
+ BO22_temprange_bdmin              3.771373e-03
+ BO2_curvelrange_ss                3.753249e-03
+ BO2_tempmax_bdmin                 3.702627e-03
+ BO21_tempmax_bdmin                3.702627e-03
+ BO22_tempmax_bdmin                3.702627e-03
+ MS_sss07_5m                       3.620328e-03
+ MS_sss08_5m                       3.527455e-03
+ BO2_templtmax_bdmin               3.513149e-03
+ BO21_templtmax_bdmin              3.513149e-03
+ BO22_templtmax_bdmin              3.513149e-03
+ BO21_curvelltmax_bdmin            3.485123e-03
+ BO22_curvelltmax_bdmin            3.485123e-03
+ BO2_silicaterange_ss              3.358419e-03
+ BO21_silicaterange_ss             3.358419e-03
+ BO22_silicaterange_ss             3.358419e-03
+ BO2_silicatemax_ss                3.307172e-03
+ BO21_silicatemax_ss               3.307172e-03
+ BO22_silicatemax_ss               3.307172e-03
+ BO2_salinityltmin_bdmean          3.305628e-03
+ BO21_salinityltmin_bdmean         3.305628e-03
+ BO22_salinityltmin_bdmean         3.305628e-03
+ BO2_nitratemean_bdmean            3.305160e-03
+ BO2_nitraterange_bdmean           3.305160e-03
+ BO21_nitratemean_bdmean           3.305160e-03
+ BO21_nitraterange_bdmean          3.305160e-03
+ BO22_nitratemean_bdmean           3.305160e-03
+ BO22_nitraterange_bdmean          3.305160e-03
+ BO2_nitraterange_ss               3.292757e-03
+ BO21_nitraterange_ss              3.292757e-03
+ BO22_nitraterange_ss              3.292757e-03
+ BO2_silicateltmax_ss              3.230627e-03
+ BO21_silicateltmax_ss             3.230627e-03
+ BO22_silicateltmax_ss             3.230627e-03
+ BO2_phosphatemean_bdmean          3.211249e-03
+ BO2_phosphaterange_bdmean         3.211249e-03
+ BO21_phosphatemean_bdmean         3.211249e-03
+ BO21_phosphaterange_bdmean        3.211249e-03
+ BO22_phosphatemean_bdmean         3.211249e-03
+ BO22_phosphaterange_bdmean        3.211249e-03
+ BO2_nitratemax_ss                 3.183336e-03
+ BO21_nitratemax_ss                3.183336e-03
+ BO22_nitratemax_ss                3.183336e-03
+ BO2_nitrateltmax_ss               3.172271e-03
+ BO21_nitrateltmax_ss              3.172271e-03
+ BO22_nitrateltmax_ss              3.172271e-03
+ MS_sss01_5m                       3.143754e-03
+ BO2_curvelmean_bdmin              3.137701e-03
+ BO2_curvelmin_bdmin               3.137701e-03
+ BO2_curvelrange_bdmin             3.137701e-03
+ BO2_phosphatemax_bdmin            2.976562e-03
+ BO21_phosphatemax_bdmin           2.976562e-03
+ BO22_phosphatemax_bdmin           2.976562e-03
+ BO21_curvelmin_ss                 2.971881e-03
+ BO22_curvelmin_ss                 2.971881e-03
+ MS_sss06_5m                       2.934520e-03
+ BO2_nitratemax_bdmin              2.910155e-03
+ BO21_nitratemax_bdmin             2.910155e-03
+ BO22_nitratemax_bdmin             2.910155e-03
+ MS_sss05_5m                       2.832750e-03
+ BO2_silicatemean_ss               2.780217e-03
+ BO21_silicatemean_ss              2.780217e-03
+ BO22_silicatemean_ss              2.780217e-03
+ BO_bathymean                      2.771905e-03
+ MS_bathy_5m                       2.770459e-03
+ BO2_salinitymax_ss                2.762308e-03
+ BO21_salinitymax_ss               2.762308e-03
+ BO22_salinitymax_ss               2.762308e-03
+ MS_biogeo01_aspect_EW_5m          2.728418e-03
+ BO21_curvelltmin_bdmean           2.723720e-03
+ BO22_curvelltmin_bdmean           2.723720e-03
+ MOL                               2.686557e-03
+ BO2_salinityltmax_ss              2.661179e-03
+ BO21_salinityltmax_ss             2.661179e-03
+ BO22_salinityltmax_ss             2.661179e-03
+ BO2_salinitymean_bdmean           2.579902e-03
+ BO2_salinityrange_bdmean          2.579902e-03
+ BO21_salinitymean_bdmean          2.579902e-03
+ BO21_salinityrange_bdmean         2.579902e-03
+ BO22_salinitymean_bdmean          2.579902e-03
+ BO22_salinityrange_bdmean         2.579902e-03
+ BO2_ironmin_bdmean                2.563716e-03
+ BO21_ironmin_bdmean               2.563716e-03
+ BO22_ironmin_bdmean               2.563716e-03
+ BO2_nitratemean_ss                2.516646e-03
+ BO21_nitratemean_ss               2.516646e-03
+ BO22_nitratemean_ss               2.516646e-03
+ HDH                               2.482398e-03
+ BO21_curvelmax_bdmax              2.458235e-03
+ BO21_curvelmean_bdmax             2.458235e-03
+ BO21_curvelmin_bdmax              2.458235e-03
+ BO21_curvelrange_bdmax            2.458235e-03
+ BO22_curvelmax_bdmax              2.458235e-03
+ BO22_curvelmean_bdmax             2.458235e-03
+ BO22_curvelmin_bdmax              2.458235e-03
+ BO22_curvelrange_bdmax            2.458235e-03
+ BO2_nitrateltmin_bdmean           2.457617e-03
+ BO21_nitrateltmin_bdmean          2.457617e-03
+ BO22_nitrateltmin_bdmean          2.457617e-03
+ MS_sss04_5m                       2.450543e-03
+ PP1A                              2.409841e-03
+ BO2_silicateltmin_ss              2.393081e-03
+ BO21_silicateltmin_ss             2.393081e-03
+ BO22_silicateltmin_ss             2.393081e-03
+ BO2_phosphateltmin_bdmean         2.360217e-03
+ BO21_phosphateltmin_bdmean        2.360217e-03
+ BO22_phosphateltmin_bdmean        2.360217e-03
+ BO21_curvelltmax_bdmean           2.348665e-03
+ BO22_curvelltmax_bdmean           2.348665e-03
+ BO2_phosphateltmax_bdmin          2.316381e-03
+ BO21_phosphateltmax_bdmin         2.316381e-03
+ BO22_phosphateltmax_bdmin         2.316381e-03
+ BO_chlomax                        2.306780e-03
+ BO2_nitrateltmax_bdmin            2.299542e-03
+ BO21_nitrateltmax_bdmin           2.299542e-03
+ BO22_nitrateltmax_bdmin           2.299542e-03
+ BO2_salinityrange_ss              2.269943e-03
+ BO21_salinityrange_ss             2.269943e-03
+ BO22_salinityrange_ss             2.269943e-03
+ BO2_salinitymin_ss                2.233010e-03
+ BO21_salinitymin_ss               2.233010e-03
+ BO22_salinitymin_ss               2.233010e-03
+ BO_bathymin                       2.189411e-03
+ BO2_curvelmin_ss                  2.139999e-03
+ BO21_curvelmean_bdmean            2.032659e-03
+ BO21_curvelrange_bdmean           2.032659e-03
+ BO22_curvelmean_bdmean            2.032659e-03
+ BO22_curvelrange_bdmean           2.032659e-03
+ EYL                               1.963380e-03
+ BO2_salinityltmin_ss              1.929975e-03
+ BO21_salinityltmin_ss             1.929975e-03
+ BO22_salinityltmin_ss             1.929975e-03
+ BO2_salinityltmax_bdmean          1.900390e-03
+ BO21_salinityltmax_bdmean         1.900390e-03
+ BO22_salinityltmax_bdmean         1.900390e-03
+ OHW                               1.878881e-03
+ BO21_curvelmax_bdmean             1.856501e-03
+ BO22_curvelmax_bdmean             1.856501e-03
+ BO2_nitrateltmin_ss               1.837123e-03
+ BO21_nitrateltmin_ss              1.837123e-03
+ BO22_nitrateltmin_ss              1.837123e-03
+ TL                                1.802423e-03
+ BO22_chloltmin_bdmin              1.793695e-03
+ UAW                               1.749412e-03
+ BO2_templtmin_bdmin               1.708450e-03
+ BO21_templtmin_bdmin              1.708450e-03
+ BO22_templtmin_bdmin              1.708450e-03
+ LLA                               1.683104e-03
+ BO2_salinitymax_bdmean            1.653968e-03
+ BO21_salinitymax_bdmean           1.653968e-03
+ BO22_salinitymax_bdmean           1.653968e-03
+ BO21_curvelmax_bdmin              1.649642e-03
+ BO22_curvelmax_bdmin              1.649642e-03
+ BO2_silicatemin_ss                1.578779e-03
+ BO21_silicatemin_ss               1.578779e-03
+ BO22_silicatemin_ss               1.578779e-03
+ BO2_salinitymean_ss               1.557497e-03
+ BO21_salinitymean_ss              1.557497e-03
+ BO22_salinitymean_ss              1.557497e-03
+ BO21_curvelltmin_bdmin            1.513892e-03
+ BO22_curvelltmin_bdmin            1.513892e-03
+ BO2_curvelltmax_ss                1.508282e-03
+ BO2_nitratemin_bdmean             1.476268e-03
+ BO21_nitratemin_bdmean            1.476268e-03
+ BO22_nitratemin_bdmean            1.476268e-03
+ BO2_nitratemin_ss                 1.414419e-03
+ BO21_nitratemin_ss                1.414419e-03
+ BO22_nitratemin_ss                1.414419e-03
+ BO21_curvelmin_bdmean             1.380927e-03
+ BO22_curvelmin_bdmean             1.380927e-03
+ BO2_phosphatemin_bdmean           1.376039e-03
+ BO21_phosphatemin_bdmean          1.376039e-03
+ BO22_phosphatemin_bdmean          1.376039e-03
+ BO21_curvelltmax_bdmax            1.351009e-03
+ BO21_curvelltmin_bdmax            1.351009e-03
+ BO22_curvelltmax_bdmax            1.351009e-03
+ BO22_curvelltmin_bdmax            1.351009e-03
+ BO22_chlomean_bdmin               1.314719e-03
+ BO22_chlomin_bdmin                1.314719e-03
+ BO22_chlorange_bdmin              1.314719e-03
+ NOW                               1.308626e-03
+ MS_biogeo02_aspect_NS_5m          1.297198e-03
+ BO2_curvelmax_ss                  1.129555e-03
+ MS_biogeo06_bathy_slope_5m        1.123282e-03
+ BO21_curvelmean_ss                1.089191e-03
+ BO22_curvelmean_ss                1.089191e-03
+ PD1                               1.084543e-03
+ BO2_ironrange_ss                  1.053225e-03
+ BO21_ironrange_ss                 1.053225e-03
+ BO22_ironrange_ss                 1.053225e-03
+ MS_biogeo03_plan_curvature_5m     8.593626e-04
+ EYH                               8.458937e-04
+ BO21_curvelmax_ss                 7.911272e-04
+ BO22_curvelmax_ss                 7.911272e-04
+ INO                               7.886964e-04
+ MS_biogeo04_profile_curvature_5m  7.845104e-04
+ ESD                               7.837651e-04
+ BO_chlorange                      7.706994e-04
+ BO2_nitrateltmin_bdmin            6.915339e-04
+ BO21_nitrateltmin_bdmin           6.915339e-04
+ BO22_nitrateltmin_bdmin           6.915339e-04
+ SP1A                              6.859124e-04
+ BO2_tempmin_bdmean                6.515468e-04
+ BO21_tempmin_bdmean               6.515468e-04
+ BO22_tempmin_bdmean               6.515468e-04
+ BO2_phosphateltmin_bdmin          6.497912e-04
+ BO21_phosphateltmin_bdmin         6.497912e-04
+ BO22_phosphateltmin_bdmin         6.497912e-04
+ BO2_curvelmin_bdmean              6.117008e-04
+ POB                               5.836919e-04
+ BO2_curvelltmin_ss                5.811129e-04
+ BO_cloudmin                       5.699433e-04
+ BO22_cloudmin                     5.699433e-04
+ BO2_ironmax_ss                    4.956840e-04
+ BO21_ironmax_ss                   4.956840e-04
+ BO22_ironmax_ss                   4.956840e-04
+ BO2_ironltmax_ss                  4.262698e-04
+ BO21_ironltmax_ss                 4.262698e-04
+ BO22_ironltmax_ss                 4.262698e-04
+ BO2_silicatemax_bdmax             3.907383e-04
+ BO2_silicatemean_bdmax            3.907383e-04
+ BO2_silicatemin_bdmax             3.907383e-04
+ BO2_silicaterange_bdmax           3.907383e-04
+ BO21_silicatemax_bdmax            3.907383e-04
+ BO21_silicatemean_bdmax           3.907383e-04
+ BO21_silicatemin_bdmax            3.907383e-04
+ BO21_silicaterange_bdmax          3.907383e-04
+ BO22_silicatemax_bdmax            3.907383e-04
+ BO22_silicatemean_bdmax           3.907383e-04
+ BO22_silicatemin_bdmax            3.907383e-04
+ BO22_silicaterange_bdmax          3.907383e-04
+ BO2_tempmean_bdmean               3.069312e-04
+ BO2_temprange_bdmean              3.069312e-04
+ BO21_tempmean_bdmean              3.069312e-04
+ BO21_temprange_bdmean             3.069312e-04
+ BO22_tempmean_bdmean              3.069312e-04
+ BO22_temprange_bdmean             3.069312e-04
+ PD2                               2.435006e-04
+ SHW                               2.058506e-04
+ BO21_curvelltmax_ss               1.541196e-04
+ BO22_curvelltmax_ss               1.541196e-04
+ MS_biogeo07_concavity_5m          9.807926e-05
+ PSP                               5.924354e-05
+ BO2_silicateltmax_bdmax           3.931279e-05
+ BO2_silicateltmin_bdmax           3.931279e-05
+ BO21_silicateltmax_bdmax          3.931279e-05
+ BO21_silicateltmin_bdmax          3.931279e-05
+ BO22_silicateltmax_bdmax          3.931279e-05
+ BO22_silicateltmin_bdmax          3.931279e-05
+ BO21_curvelrange_ss               1.635399e-05
+ BO22_curvelrange_ss               1.635399e-05
<none>                              0.000000e+00
+ BO2_ironmin_ss                   -5.218986e-05
+ BO21_ironmin_ss                  -5.218986e-05
+ BO22_ironmin_ss                  -5.218986e-05
+ PG1                              -1.731057e-04
+ ANF                              -1.838903e-04
+ UAH                              -1.951250e-04
+ BO2_ironmean_ss                  -2.724383e-04
+ BO21_ironmean_ss                 -2.724383e-04
+ BO22_ironmean_ss                 -2.724383e-04
+ BO2_nitratemean_bdmin            -2.969745e-04
+ BO2_nitratemin_bdmin             -2.969745e-04
+ BO2_nitraterange_bdmin           -2.969745e-04
+ BO21_nitratemean_bdmin           -2.969745e-04
+ BO21_nitratemin_bdmin            -2.969745e-04
+ BO21_nitraterange_bdmin          -2.969745e-04
+ BO22_nitratemean_bdmin           -2.969745e-04
+ BO22_nitratemin_bdmin            -2.969745e-04
+ BO22_nitraterange_bdmin          -2.969745e-04
+ MS_biogeo09_sss_min_5m           -3.162320e-04
+ BO2_phosphatemean_bdmin          -3.290049e-04
+ BO2_phosphatemin_bdmin           -3.290049e-04
+ BO2_phosphaterange_bdmin         -3.290049e-04
+ BO21_phosphatemean_bdmin         -3.290049e-04
+ BO21_phosphatemin_bdmin          -3.290049e-04
+ BO21_phosphaterange_bdmin        -3.290049e-04
+ BO22_phosphatemean_bdmin         -3.290049e-04
+ BO22_phosphatemin_bdmin          -3.290049e-04
+ BO22_phosphaterange_bdmin        -3.290049e-04
+ PP1                              -3.838145e-04
+ BO2_ironltmin_ss                 -3.922944e-04
+ BO21_ironltmin_ss                -3.922944e-04
+ BO22_ironltmin_ss                -3.922944e-04
+ MOW                              -4.365201e-04
+ PP2                              -4.489243e-04
+ PEY                              -4.827085e-04
+ BO2_templtmin_bdmean             -7.903115e-04
+ BO21_templtmin_bdmean            -7.903115e-04
+ BO22_templtmin_bdmean            -7.903115e-04
+ HDL                              -8.610524e-04
+ BO2_salinitymax_bdmax            -8.805052e-04
+ BO2_salinitymean_bdmax           -8.805052e-04
+ BO2_salinitymin_bdmax            -8.805052e-04
+ BO2_salinityrange_bdmax          -8.805052e-04
+ BO21_salinitymax_bdmax           -8.805052e-04
+ BO21_salinitymean_bdmax          -8.805052e-04
+ BO21_salinitymin_bdmax           -8.805052e-04
+ BO21_salinityrange_bdmax         -8.805052e-04
+ BO22_salinitymax_bdmax           -8.805052e-04
+ BO22_salinitymean_bdmax          -8.805052e-04
+ BO22_salinitymin_bdmax           -8.805052e-04
+ BO22_salinityrange_bdmax         -8.805052e-04
+ BO2_salinityltmax_bdmax          -8.878034e-04
+ BO2_salinityltmin_bdmax          -8.878034e-04
+ BO21_salinityltmax_bdmax         -8.878034e-04
+ BO21_salinityltmin_bdmax         -8.878034e-04
+ BO22_salinityltmax_bdmax         -8.878034e-04
+ BO22_salinityltmin_bdmax         -8.878034e-04
+ BO2_silicateltmin_bdmean         -1.016931e-03
+ BO21_silicateltmin_bdmean        -1.016931e-03
+ BO22_silicateltmin_bdmean        -1.016931e-03
+ BO2_silicatemax_bdmean           -1.085877e-03
+ BO21_silicatemax_bdmean          -1.085877e-03
+ BO22_silicatemax_bdmean          -1.085877e-03
+ BO2_silicatemean_bdmean          -1.112205e-03
+ BO2_silicaterange_bdmean         -1.112205e-03
+ BO21_silicatemean_bdmean         -1.112205e-03
+ BO21_silicaterange_bdmean        -1.112205e-03
+ BO22_silicatemean_bdmean         -1.112205e-03
+ BO22_silicaterange_bdmean        -1.112205e-03
+ BO2_silicatemin_bdmean           -1.157677e-03
+ BO21_silicatemin_bdmean          -1.157677e-03
+ BO22_silicatemin_bdmean          -1.157677e-03
+ BO2_silicateltmax_bdmean         -1.243803e-03
+ BO21_silicateltmax_bdmean        -1.243803e-03
+ BO22_silicateltmax_bdmean        -1.243803e-03
+ BO2_silicatemax_bdmin            -1.309299e-03
+ BO21_silicatemax_bdmin           -1.309299e-03
+ BO22_silicatemax_bdmin           -1.309299e-03
+ BO2_silicatemean_bdmin           -1.460843e-03
+ BO2_silicatemin_bdmin            -1.460843e-03
+ BO2_silicaterange_bdmin          -1.460843e-03
+ BO21_silicatemean_bdmin          -1.460843e-03
+ BO21_silicatemin_bdmin           -1.460843e-03
+ BO21_silicaterange_bdmin         -1.460843e-03
+ BO22_silicatemean_bdmin          -1.460843e-03
+ BO22_silicatemin_bdmin           -1.460843e-03
+ BO22_silicaterange_bdmin         -1.460843e-03
+ BO2_silicateltmin_bdmin          -1.489194e-03
+ BO21_silicateltmin_bdmin         -1.489194e-03
+ BO22_silicateltmin_bdmin         -1.489194e-03
+ BO2_silicateltmax_bdmin          -1.691780e-03
+ BO21_silicateltmax_bdmin         -1.691780e-03
+ BO22_silicateltmax_bdmin         -1.691780e-03
}

m.ord_Meas <- m.ord_All

save(m.ord_All, file="RDA_meas_geo_env.gz", compress=T)
#load("RDA_meas_geo_env.gz")

m.ord_All
{ #Results
Call: rda(formula = X ~ Geo_dist + BO_ph, data = RDA_data_scale)

                Inertia Proportion Rank
Total         1.260e+03  1.000e+00
Constrained   1.174e+02  9.318e-02    2
Unconstrained 1.143e+03  9.068e-01   24
Inertia is variance

Eigenvalues for constrained axes:
 RDA1  RDA2
64.92 52.49

Eigenvalues for unconstrained axes:
  PC1   PC2   PC3   PC4   PC5   PC6   PC7   PC8
52.96 51.27 50.97 50.65 50.37 50.14 49.61 49.14
(Showing 8 of 24 unconstrained eigenvalues)
}
attributes(m.ord_All)$VIF_remove_list
{ #Results
     Rejected_variable Max_VIF            Formula
[1,] "BO_nitrate"      "3.3780165780482"  "X ~ Geo_dist"
[2,] "BO_parmean"      "4.30734884628966" "X ~ Geo_dist"
[3,] "BO22_parmean"    "4.30734884628966" "X ~ Geo_dist"
}
RsquareAdj(m.ord_All)
{ #Results
$r.squared
[1] 0.09318018

$adj.r.squared
[1] 0.01761186
}
m.ord_All$CCA$eig[1:4]/m.ord_All$tot.chi
{ #Results
      RDA1       RDA2       <NA>       <NA>
0.05152537 0.04165481         NA         NA
}
anova(m.ord_All, parallel=20, permutations = 9999)
{ #Results
Permutation test for rda under reduced model
Permutation: free
Number of permutations: 9999

Model: rda(formula = X ~ Geo_dist + BO_ph, data = RDA_data_scale)
         Df Variance      F Pr(>F)
Model     2   117.41 1.2331  1e-04 ***
Residual 24  1142.61
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
tmp_anova <- anova(m.ord_All, by="axis", parallel=20, permutations = 9999)
tmp_anova
{ #Results
Permutation test for rda under reduced model
Forward tests for axes
Permutation: free
Number of permutations: 9999

Model: rda(formula = X ~ Geo_dist + BO_ph, data = RDA_data_scale)
         Df Variance      F Pr(>F)
RDA1      1    64.92 1.3637  1e-04 ***
RDA2      1    52.49 1.1024  2e-04 ***
Residual 24  1142.61
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
tmp_anova$Variance[1:4]/sum(tmp_anova$Variance)
{ #Results
[1] 0.05152537 0.04165481 0.90681982         NA
}
anova_MARG <- anova(m.ord_All, by="margin", parallel=20)
anova_MARG
{ #Results
Permutation test for rda under reduced model
Marginal effects of terms
Permutation: free
Number of permutations: 999

Model: rda(formula = X ~ Geo_dist + BO_ph, data = RDA_data_scale)
         Df Variance      F Pr(>F)
Geo_dist  1    57.37 1.2051  0.001 ***
BO_ph     1    52.82 1.1095  0.001 ***
Residual 24  1142.61
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
}
cbind(rownames(anova_MARG),round(anova_MARG$Variance[1:length(anova_MARG$Variance)]/sum(anova_MARG$Variance),6))
{ #Results
     [,1]       [,2]      )
[1,] "Geo_dist" "0.045796"
[2,] "BO_ph"    "0.042163"
[3,] "Residual" "0.912041"
}
#What are these factors?
env_data_list[env_data_list$layer_code == "BO_ph", 1:6]
{ #Results
   dataset_code layer_code name                      description terrestrial
84   Bio-ORACLE      BO_ph   pH Measure of acidity in the ocean.       FALSE
   marine
84   TRUE
}
} #notepad cleanup

#Identiying the loci associated with SST
{```{R}```
#Getting the RDA and the loadings
tmp.RDA <- rda(X ~ MS_sst09_5m + Condition(BO_nitrate + BO2_nitratemin_ss), data = RDA_data_scale)
RDA.load <- abs(tmp.RDA$CCA$v[,1])
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]

#Getting the initial p-value
tmp.anova <- anova(tmp.RDA, parallel=20)
iteration <- 0
if(exists("loc.rm")){rm(loc.rm)}
if(exists("Loci")){rm(Loci)}

while(tmp.anova$"Pr(>F)"[1] < 0.05){
#Get Loci to remove
if(!exists("Loci")){Loci <- names(RDA.load[1])
} else {Loci <- c(Loci, names(RDA.load[1]))}
loc.rm <- unlist(Loci_names(Loci,"[.]"))
#Remove loci from genind object
set.loc <- locNames(gen3)[!locNames(gen3) %in% loc.rm]
gen.keep <- gen3[,loc = set.loc, drop=T]
#Scaling genind object
X <- scaleGen(gen.keep, NA.method="mean", scale=F)
tmp.rda <- rda(X ~ MS_sst09_5m + Condition(BO_nitrate + BO2_nitratemin_ss), data = RDA_data_scale)
tmp.anova <- anova(tmp.rda, parallel=20)
print(paste("Loci removed:", length(loc.rm), '	', "AMOVA p-value", tmp.anova$"Pr(>F)"[1]))
RDA.load <- abs(tmp.rda$CCA$v[,1])
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]
}
length(loc.rm)
#[1] 772

SST.loc <- loc.rm
write.table(loc.rm, file="SST_loci.txt", quote=F, col.names=F, row.names=F)
#SST.loc <- read.table("SST_loci.txt", head=F)[,1]
}}}} #notepad cleanup

#Identiying the loci associated with BO_nitrate
{```{R}```
#Getting the RDA and the loadings
tmp.RDA <- rda(X ~ BO_nitrate + Condition(MS_sst09_5m + BO2_nitratemin_ss), data = RDA_data_scale)
RDA.load <- abs(tmp.RDA$CCA$v[,1])
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]

#Getting the initial p-value
tmp.anova <- anova(tmp.RDA, parallel=20)
tmp.anova
iteration <- 0
if(exists("loc.rm")){rm(loc.rm)}
if(exists("Loci")){rm(Loci)}

while(tmp.anova$"Pr(>F)"[1] < 0.05){
#Get Loci to remove
if(!exists("Loci")){Loci <- names(RDA.load[1])
} else {Loci <- c(Loci, names(RDA.load[1]))}
loc.rm <- unlist(Loci_names(Loci,"[.]"))
#Remove loci from genind object
set.loc <- locNames(gen3)[!locNames(gen3) %in% loc.rm]
gen.keep <- gen3[,loc = set.loc, drop=T]
#Scaling genind object
X <- scaleGen(gen.keep, NA.method="mean", scale=F)
tmp.rda <- rda(X ~ BO_nitrate + Condition(MS_sst09_5m + BO2_nitratemin_ss), data = RDA_data_scale)
tmp.anova <- anova(tmp.rda, parallel=20)
print(paste("Loci removed:", length(loc.rm), '	', "AMOVA p-value", tmp.anova$"Pr(>F)"[1]))
RDA.load <- abs(tmp.rda$CCA$v[,1])
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]
}
length(loc.rm)
#[1] 127

Nitrate.loc <- loc.rm
write.table(loc.rm, file="Nitrate_loci_Jan2024.txt", quote=F, col.names=F, row.names=F)
#Nitrate.loc <- read.table("Nitrate_loci.txt", head=F)[,1]
}}}} #notepad cleanup

#Identiying the loci associated with nitrate
{```{R}```
#Getting the RDA and the loadings
tmp.RDA <- rda(X ~ BO2_nitratemin_ss + Condition(MS_sst09_5m + BO_nitrate), data = RDA_data_scale)
RDA.load <- abs(tmp.RDA$CCA$v[,1])
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]

#Getting the initial p-value
tmp.anova <- anova(tmp.RDA, parallel=20)
tmp.anova
iteration <- 0
if(exists("loc.rm")){rm(loc.rm)}
if(exists("Loci")){rm(Loci)}

while(tmp.anova$"Pr(>F)"[1] < 0.05){
#Get Loci to remove
if(!exists("Loci")){Loci <- names(RDA.load[1])
} else {Loci <- c(Loci, names(RDA.load[1]))}
loc.rm <- unlist(Loci_names(Loci,"[.]"))
#Remove loci from genind object
set.loc <- locNames(gen3)[!locNames(gen3) %in% loc.rm]
gen.keep <- gen3[,loc = set.loc, drop=T]
#Scaling genind object
X <- scaleGen(gen.keep, NA.method="mean", scale=F)
tmp.rda <- rda(X ~ BO2_nitratemin_ss + Condition(MS_sst09_5m + BO_nitrate), data = RDA_data_scale)
tmp.anova <- anova(tmp.rda, parallel=20)
print(paste("Loci removed:", length(loc.rm), '	', "AMOVA p-value", tmp.anova$"Pr(>F)"[1]))
RDA.load <- abs(tmp.rda$CCA$v[,1])
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]
}
length(loc.rm)
#[1] 19
#[1] 18	v2

Nitrate_ss.loc <- loc.rm
write.table(loc.rm, file="Nitrate_ss_loci_v2.txt", quote=F, col.names=F, row.names=F)
#Nitrate_ss.loc <- read.table("Nitrate_ss_loci.txt", head=F)[,1]
}}}} #notepad cleanup

#Identiying the loci associated with the RDA (Multiple axes v2)
{```{R}```
#Getting the RDA and the loadings
set.seed(1235)
X <- scaleGen(gen3, NA.method="mean", scale=F)
tmp.rda <- rda(X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
RDA.load <- apply(tmp.rda$CCA$v, 1, function(x) sum(abs(x)))
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]

#Getting the initial p-value
tmp.anova <- anova(tmp.rda, parallel=20)
tmp.anova
iteration <- 0
if(exists("loc.rm")){rm(loc.rm)}
if(exists("Loci")){rm(Loci)}
if(exists("taken_loci")){rm(taken_loci)}

while(tmp.anova$"Pr(>F)"[1] < 0.05){
#Testing loci
tmp_loci <- names(RDA.load)[1]
#print(test_loci)
for(i in 1:ncol(tmp.rda$CCA$v)){
tmp.load <- abs(tmp.rda$CCA$v[,i])
tmp.load <- tmp.load[order(tmp.load, decreasing = T)]
tmp_loci <- c(tmp_loci, names(tmp.load)[1])
#print(names(tmp.load)[1])
}
test_loci <- unique(unlist(Loci_names(tmp_loci,"[.]")))

tmp.var <- vector()
for(i in test_loci){
gen.test <- gen3[,loc = i, drop=T]
X <- scaleGen(gen.test, NA.method="mean", scale=F)
tmp.rda <- rda(X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
tmp.var <- c(tmp.var, sum(tmp.rda$CCA$eig/tmp.rda$tot.chi))
}

#Get Loci to remove
if(!exists("Loci")){Loci <- test_loci[which(tmp.var == max(tmp.var))]
} else {Loci <- c(Loci, test_loci[which(tmp.var == max(tmp.var))])}
#Record which column it came from
if(!exists("taken_loci")){taken_loci <- list(grep(tail(Loci,n=1),tmp_loci))
}else{taken_loci <- append(taken_loci, (grep(tail(Loci,n=1),tmp_loci)))}

#Remove loci from genind object
loc.rm <- Loci
set.loc <- locNames(gen3)[!locNames(gen3) %in% loc.rm]
gen.keep <- gen3[,loc = set.loc, drop=T]
#Scaling genind object
X <- scaleGen(gen.keep, NA.method="mean", scale=F)
tmp.rda <- rda(X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
tmp.anova <- anova(tmp.rda, parallel=20)
print(paste("Loci removed:", length(loc.rm), '	', "AMOVA p-value", tmp.anova$"Pr(>F)"[1]))
RDA.load <- apply(tmp.rda$CCA$v, 1, function(x) sum(abs(x)))
RDA.load <- RDA.load[order(RDA.load, decreasing = T)]
}
length(loc.rm)
#[1] 2665

All.loc <- loc.rm
write.table(loc.rm, file="RDA_loci_v3.txt", quote=F, col.names=F, row.names=F)
#All.loc <- read.table("RDA_loci_v3.txt", head=F)[,1]

#How many loci are associated with the locus
table(unlist(lapply(taken_loci,length)))
#   1
#1273

#Which categories are most and least seen
table(unlist(taken_loci))
#  1   2   3   4
#243 650 216 164

#1 is sum of all the values
#2 is RDA axis 1
#3 is RDA axis 2
#4 is RDA axis 3
}}}}}} #notepad cleanup

#Comparing loci between lists
{```{R}```
#tmp.m <- t(combn(c("SST","Nitrate","Nitrate_ss","All"),2))
#tmp.m <- t(combn(c("SST","Nitrate","Nitrate_ss","OF","BS"),2))

tmp.m <- t(combn(c("SST","Nitrate","Nitrate_ss"),2))
loc_comp <- data.frame(matrix(ncol=3, nrow=nrow(tmp.m)))
colnames(loc_comp) <- c("Fac_1", "Fac_2", "Shared")

for(i in 1:nrow(loc_comp)){
tmp.n <- length(which(get(paste(tmp.m[i,1],"loc",sep=".")) %in% get(paste(tmp.m[i,2],"loc",sep="."))))
loc_comp[i,] <- c(tmp.m[i,1], tmp.m[i,2], tmp.n)
}

loc_comp
{ #Results
       Fac_1      Fac_2 Shared
1        SST    Nitrate     42
2        SST Nitrate_ss      2
3        SST        All    595
4    Nitrate Nitrate_ss      2
5    Nitrate        All    127
6 Nitrate_ss        All     19
}

 Nitrate_ss.loc[which(Nitrate_ss.loc %in% Nitrate.loc)]
[1] dDocent_Contig_22146 dDocent_Contig_21276
19 Levels: dDocent_Contig_1039 dDocent_Contig_11855 ... dDocent_Contig_5435
> Nitrate_ss.loc[which(Nitrate_ss.loc %in% SST.loc)]
[1] dDocent_Contig_21276 dDocent_Contig_37176
19 Levels: dDocent_Contig_1039 dDocent_Contig_11855 ... dDocent_Contig_5435
} #notepad cleanup

#Making Venn diagram of the RDA loci
{```{R}```
venn.diagram(list(SST = SST.loc, Full=All.loc, Nitrate= Nitrate.loc, Iron = Nitrate_ss.loc),
category.names=c(paste(" SST\n(n=", length(SST.loc), ")", sep = ""), 
paste("Full\n(n=", length(All.loc), ")", sep = ""),
paste("Nitrate\n(n=", length(Nitrate.loc), ")", sep = ""), 
paste("Nitrate_ss\n(n=", length(Nitrate_ss.loc), ")", sep = "")),
"Venn_sig_sites.tiff", fill=c("red4", "darkgreen", "mediumblue", "yellow3"))
} #notepad cleanup

#Making PCA of the RDA loci
{```{R}```
gen.SST <- gen3[,loc=SST.loc]
gen.Nit <- gen3[,loc=Nitrate.loc]
gen.Nit_ss <- gen3[,loc=Nitrate_ss.loc]

X.SST <- scaleGen(gen.SST, NA.method="mean", scale=F)
pca.SST <- dudi.pca(X.SST,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.Nit <- scaleGen(gen.Nit, NA.method="mean", scale=F)
pca.Nit <- dudi.pca(X.Nit,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.Nit_ss <- scaleGen(gen.Nit_ss, NA.method="mean", scale=F)
pca.Nit_ss <- dudi.pca(X.Nit_ss,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen.SST) <- setPop(gen.Nit) <- setPop(gen.Nit_ss) <- ~G3
par(mfrow=c(3,1))
ade4::s.class(pca.SST$li, pop(gen.SST), col=c("chocolate2", "green", "mediumblue"), cstar=0, axesell=F, clabel=0)
#mtext(paste("x-axis variation: ",format(round(100*pca.SST$eig[1]/sum(pca.SST$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
#mtext(paste("y-axis variation: ",format(round(100*pca.SST$eig[2]/sum(pca.SST$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("SST loci (n=772)", adj=0.05)
legend(8,5,legend=levels(pop(gen.Nit_ss)), pch=19, col=c("chocolate2", "green", "mediumblue"), bty="n")

ade4::s.class(pca.Nit$li, pop(gen.Nit), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("Nitrate loci", adj=0.05)
legend(-1,-2,legend=levels(pop(gen.Nit_ss)), pch=19, col=c3[c(2,1,3)], bty="n")

ade4::s.class(pca.Nit_ss$li, pop(gen.Nit_ss), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("ss Nitrate loci", adj=0.05)
legend(2,-0.5,legend=levels(pop(gen.Nit_ss)), pch=19, col=c3[c(2,1,3)], bty="n")

png("PCA_SST_loci.png", res=200, width=2000, height=2000)
ade4::s.class(pca.SST$li, pop(gen.SST), col=c("chocolate2", "green", "mediumblue"), cstar=0, axesell=F, clabel=0)
#mtext(paste("x-axis variation: ",format(round(100*pca.SST$eig[1]/sum(pca.SST$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
#mtext(paste("y-axis variation: ",format(round(100*pca.SST$eig[2]/sum(pca.SST$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext(paste("SST loci (n=",nLoc(gen.SST),")",sep=""), adj=0.05, line=2)
legend(-8.5,7,legend=levels(pop(gen.SST)), pch=19, col=c("chocolate2", "green", "mediumblue"), bty="n")
dev.off()

png("PCA_Nit_loci.png", res=200, width=2000, height=2000)
ade4::s.class(pca.Nit$li, pop(gen.Nit), col=c("chocolate2", "green", "mediumblue"), cstar=0, axesell=F, clabel=0)
#mtext(paste("x-axis variation: ",format(round(100*pca.Nit$eig[1]/sum(pca.Nit$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
#mtext(paste("y-axis variation: ",format(round(100*pca.Nit$eig[2]/sum(pca.Nit$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext(paste("Nitrate loci (n=",nLoc(gen.Nit),")",sep=""), adj=0.05, line=2)
legend(-1,-2,legend=levels(pop(gen.Nit)), pch=19, col=c("chocolate2", "green", "mediumblue"), bty="n")
dev.off()

png("PCA_Nit_ss_loci.png", res=200, width=2000, height=2000)
ade4::s.class(pca.Nit_ss$li, pop(gen.Nit_ss), col=c("chocolate2", "green", "mediumblue"), cstar=0, axesell=F, clabel=0)
#mtext(paste("x-axis variation: ",format(round(100*pca.Nit_ss$eig[1]/sum(pca.Nit_ss$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
#mtext(paste("y-axis variation: ",format(round(100*pca.Nit_ss$eig[2]/sum(pca.Nit_ss$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext(paste("ss Nitrate loci (n=",nLoc(gen.Nit_ss),")",sep=""), adj=0.05, line=2)
legend(-2.5,-1,legend=levels(pop(gen.Nit_ss)), pch=19, col=c("chocolate2", "green", "mediumblue"), bty="n")
dev.off()
}

#Biplots
{```{R}```
library('car')
X <- scaleGen(gen3, NA.method="mean", scale=F)
Final_rda <- rda(formula = X ~ MS_sst09_5m + BO_nitrate + BO2_nitratemin_ss, data = RDA_data_scale)
#Biplot of all the full model
Atl.indv <- as.character(as.matrix(gen3@strata$INDV[gen3@strata$G3=="Atlantic"]))
East.indv <- as.character(as.matrix(gen3@strata$INDV[gen3@strata$G3=="East"]))
West.indv <- as.character(as.matrix(gen3@strata$INDV[gen3@strata$G3=="West"]))

col_pts <- rownames(Final_rda$CCA$wa)
col_pts[col_pts %in% Atl.indv] <- "chocolate2"
col_pts[col_pts %in% East.indv] <- "green"
col_pts[col_pts %in% West.indv] <- "mediumblue"

#Axes 1 & 2
MIN.x <- min(c(Final_rda$CCA$wa[,1]*14.5, Final_rda$CCA$biplot[,1]*2))*1.01
MAX.x <- max(c(Final_rda$CCA$wa[,1]*14.5, Final_rda$CCA$biplot[,1]*2))*1.01
MIN.y <- min(c(Final_rda$CCA$wa[,2]*14.5, Final_rda$CCA$biplot[,2]*2))*1.01
MAX.y <- max(c(Final_rda$CCA$wa[,2]*14.5, Final_rda$CCA$biplot[,2]*2))*1.01

png("Biplot_RDA_full_A12.png", res=, width=2000, height=2000)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(Final_rda$CCA$wa*14.5, pch=19, col=col_pts, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), main="Full RDA model", xlab="RDA Axis 1", ylab="RDA Axis 2", cex.lab=2, cex.main=2)
abline(v=0,h=0,col="black")
for(i in 1:nrow(Final_rda$CCA$biplot)){arrows(0,0,Final_rda$CCA$biplot[i,1]*1.5, Final_rda$CCA$biplot[i,2]*1.5, col="black")}
legend(-3, -2, legend=c("Atlantic (n=17)","East (n=21)","West (n=24)"), pch=19, col=c("chocolate2","green","mediumblue"), bty="n", cex=2)
text(Final_rda$CCA$biplot[1,1]*1.8, Final_rda$CCA$biplot[1,2]*2.5, labels=rownames(Final_rda$CCA$biplot)[1], cex=2)
text(Final_rda$CCA$biplot[2,1]*0.7, Final_rda$CCA$biplot[2,2]*1.1, labels=rownames(Final_rda$CCA$biplot)[2], cex=2)
text(Final_rda$CCA$biplot[3,1]*3, Final_rda$CCA$biplot[3,2]*0.9, labels=rownames(Final_rda$CCA$biplot)[3], cex=2)
dataEllipse(Final_rda$CCA$wa[Atl.indv, 1:2]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2")
dataEllipse(Final_rda$CCA$wa[East.indv, 1:2]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green")
dataEllipse(Final_rda$CCA$wa[West.indv, 1:2]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
dev.off()

#Axes 1 & 3
MIN.x <- min(c(Final_rda$CCA$wa[,1]*14.5, Final_rda$CCA$biplot[,1]*2))*1.01
MAX.x <- max(c(Final_rda$CCA$wa[,1]*14.5, Final_rda$CCA$biplot[,1]*2))*1.01
MIN.y <- min(c(Final_rda$CCA$wa[,3]*14.5, Final_rda$CCA$biplot[,3]*2))*1.01
MAX.y <- max(c(Final_rda$CCA$wa[,3]*14.5, Final_rda$CCA$biplot[,3]*2))*1.01

png("Biplot_RDA_full_A13.png", res=, width=2000, height=2000)
par(mar = c(5.1, 5.1, 4.1, 2.1))
plot(Final_rda$CCA$wa[,c(1,3)]*14.5, pch=19, col=col_pts, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), main="Full RDA model", xlab="RDA Axis 1", ylab="RDA Axis 3", cex.lab=2, cex.main=2)
abline(v=0,h=0,col="black")
for(i in 1:nrow(Final_rda$CCA$biplot)){arrows(0,0,Final_rda$CCA$biplot[i,1]*1.5, Final_rda$CCA$biplot[i,3]*1.5, col="black")}
legend(-3, 5, legend=c("Atlantic (n=17)","East (n=21)","West (n=24)"), pch=19, col=c("chocolate2","green","mediumblue"), bty="n", cex=2)
text(Final_rda$CCA$biplot[1,1]*1.9, Final_rda$CCA$biplot[1,3]*2, labels=rownames(Final_rda$CCA$biplot)[1], cex=2)
text(Final_rda$CCA$biplot[2,1]*0.6, Final_rda$CCA$biplot[2,3]*1.1, labels=rownames(Final_rda$CCA$biplot)[2], cex=2)
text(Final_rda$CCA$biplot[3,1]*2.2, Final_rda$CCA$biplot[3,3]*1.6, labels=rownames(Final_rda$CCA$biplot)[3], cex=2)
dataEllipse(Final_rda$CCA$wa[Atl.indv, c(1,3)]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2")
dataEllipse(Final_rda$CCA$wa[East.indv, c(1,3)]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green")
dataEllipse(Final_rda$CCA$wa[West.indv, c(1,3)]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
dev.off()

#Biplot of SST
SST_rda <- rda(formula = X ~ MS_sst09_5m + Condition(BO_nitrate + BO2_nitratemin_ss), data = RDA_data_scale)

col_pts <- rownames(SST_rda$CCA$wa)
col_pts[col_pts %in% Atl.indv] <- "red4"
col_pts[col_pts %in% East.indv] <- "dodgerblue"
col_pts[col_pts %in% West.indv] <- "mediumblue"

MIN.x <- min(c(SST_rda$CCA$wa[,1]*14.5, SST_rda$CCA$biplot[,1]*2))*1.01
MAX.x <- max(c(SST_rda$CCA$wa[,1]*14.5, SST_rda$CCA$biplot[,1]*2))*1.01
MIN.y <- min(SST_rda$CA$u[,1]*14.5)*1.01
MAX.y <- max(SST_rda$CA$u[,1]*14.5)*1.01

png("Biplot_RDA_SST.png", res=, width=2000, height=2000)
plot(SST_rda$CCA$wa[,1]*14.5, SST_rda$CA$u[,1]*14.5, pch=19, col=col_pts, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), main="SST RDA model", xlab="RDA Axis 1", ylab="PC Axis 1")
abline(v=0,h=0,col="black")
arrows(0,0,SST_rda$CCA$biplot[1,1]*1.5, SST_rda$CCA$biplot[1,2]*1.5, col="black")
legend(-4, 7, legend=c("Atlantic","East","West"), pch=19, col=c("red4","dodgerblue","mediumblue"), bty="n")
text(SST_rda$CCA$biplot[1,1]*1.5, SST_rda$CCA$biplot[1,2]*3, labels=rownames(SST_rda$CCA$biplot)[1])
dataEllipse(SST_rda$CCA$wa[Atl.indv, 1]*14.5, SST_rda$CA$u[Atl.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="red4")
dataEllipse(SST_rda$CCA$wa[East.indv, 1]*14.5, SST_rda$CA$u[East.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="dodgerblue")
dataEllipse(SST_rda$CCA$wa[West.indv, 1]*14.5, SST_rda$CA$u[West.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
dev.off()

#Biplot of Nitrate
Nit_rda <- rda(formula = X ~ BO_nitrate + Condition(MS_sst09_5m + BO2_nitratemin_ss), data = RDA_data_scale)

col_pts <- rownames(Nit_rda$CCA$wa)
col_pts[col_pts %in% gen3@strata$INDV[gen3@strata$G3=="Atlantic"]] <- "red4"
col_pts[col_pts %in% gen3@strata$INDV[gen3@strata$G3=="East"]] <- "dodgerblue"
col_pts[col_pts %in% gen3@strata$INDV[gen3@strata$G3=="West"]] <- "mediumblue"

MIN.x <- min(c(Nit_rda$CCA$wa[,1]*14.5, Nit_rda$CCA$biplot[,1]*2))*1.01
MAX.x <- max(c(Nit_rda$CCA$wa[,1]*14.5, Nit_rda$CCA$biplot[,1]*2))*1.01
MIN.y <- min(Nit_rda$CA$u[,1]*14.5)*1.01
MAX.y <- max(Nit_rda$CA$u[,1]*14.5)*1.01

png("Biplot_RDA_Nit.png", res=, width=2000, height=2000)
plot(Nit_rda$CCA$wa[,1]*14.5, Nit_rda$CA$u[,1]*14.5, pch=19, col=col_pts, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), main="Nitrate RDA model", xlab="RDA Axis 1", ylab="PC Axis 1")
abline(v=0,h=0,col="black")
arrows(0,0,Nit_rda$CCA$biplot[1,1]*1.5, 0, col="black")
legend(3, -3, legend=c("Atlantic","East","West"), pch=19, col=c("red4","dodgerblue","mediumblue"), bty="n")
text(Nit_rda$CCA$biplot[1,1]*0.7, -0.2, labels=rownames(Nit_rda$CCA$biplot)[1])
dataEllipse(Nit_rda$CCA$wa[Atl.indv, 1]*14.5, Nit_rda$CA$u[Atl.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="red4")
dataEllipse(Nit_rda$CCA$wa[East.indv, 1]*14.5, Nit_rda$CA$u[East.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="dodgerblue")
dataEllipse(Nit_rda$CCA$wa[West.indv, 1]*14.5, Nit_rda$CA$u[West.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
dev.off()

#Biplot of Nitreate ss
Nit_ss_rda <- rda(formula = X ~ BO2_nitratemin_ss + Condition(MS_sst09_5m + BO_nitrate), data = RDA_data_scale)

col_pts <- rownames(Nit_ss_rda$CCA$wa)
col_pts[col_pts %in% gen3@strata$INDV[gen3@strata$G3=="Atlantic"]] <- "red4"
col_pts[col_pts %in% gen3@strata$INDV[gen3@strata$G3=="East"]] <- "dodgerblue"
col_pts[col_pts %in% gen3@strata$INDV[gen3@strata$G3=="West"]] <- "mediumblue"

MIN.x <- min(c(Nit_ss_rda$CCA$wa[,1]*14.5, Nit_ss_rda$CCA$biplot[,1]*2))*1.01
MAX.x <- max(c(Nit_ss_rda$CCA$wa[,1]*14.5, Nit_ss_rda$CCA$biplot[,1]*2))*1.01
MIN.y <- min(Nit_ss_rda$CA$u[,1]*14.5)*1.01
MAX.y <- max(Nit_ss_rda$CA$u[,1]*14.5)*1.01

png("Biplot_RDA_Nit_ss.png", res=, width=2000, height=2000)
plot(Nit_ss_rda$CCA$wa[,1]*14.5, Nit_ss_rda$CA$u[,1]*14.5, pch=19, col=col_pts, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), main="Min Sea Surface Nitrate RDA model", xlab="RDA Axis 1", ylab="PC Axis 1")
abline(v=0,h=0,col="black")
arrows(0,0,Nit_ss_rda$CCA$biplot[1,1]*1.5, 0, col="black")
legend(-5.5, -2.5, legend=c("Atlantic","East","West"), pch=19, col=c("red4","dodgerblue","mediumblue"), bty="n")
text(Nit_ss_rda$CCA$biplot[1,1]*0.7, -0.2, labels=rownames(Nit_ss_rda$CCA$biplot)[1])
dataEllipse(Nit_ss_rda$CCA$wa[Atl.indv, 1]*14.5, Nit_ss_rda$CA$u[Atl.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="red4")
dataEllipse(Nit_ss_rda$CCA$wa[East.indv, 1]*14.5, Nit_ss_rda$CA$u[East.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="dodgerblue")
dataEllipse(Nit_ss_rda$CCA$wa[West.indv, 1]*14.5, Nit_ss_rda$CA$u[West.indv,1]*14.5, levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue")
dev.off()
}

#Comparing loci between outlier and env lists
{```{R}```
SST.loc <- read.table("SST_loci.txt", head=F)[,1]
Nitrate.loc <- read.table("Nitrate_loci.txt", head=F)[,1]
Nitrate_ss.loc <- read.table("Nitrate_ss_loci.txt", head=F)[,1]

env.loc <- unique(c(as.character(as.matrix(SST.loc)), as.character(as.matrix(Nitrate.loc)), as.character(as.matrix(Nitrate_ss.loc))))
length(env.loc)
#[1] 873

OF.loc <- read.table("Outflank_Fst_Outliers.list", head=F)[,1]
BS.loc <- read.table("Bayescan/Bayescan_Outliers.list", head=F)[,1]

tmp.m <- t(combn(c("env","OF","BS"),2))
loc_comp <- data.frame(matrix(ncol=3, nrow=nrow(tmp.m)))
colnames(loc_comp) <- c("Fac_1", "Fac_2", "Shared")

for(i in 1:nrow(loc_comp)){
tmp.n <- length(which(get(paste(tmp.m[i,1],"loc",sep=".")) %in% get(paste(tmp.m[i,2],"loc",sep="."))))
loc_comp[i,] <- c(tmp.m[i,1], tmp.m[i,2], tmp.n)
}

loc_comp
{ #Results
  Fac_1 Fac_2 Shared
1   env    OF      2
2   env    BS      3
3    OF    BS     10
}
} #notepad cleanup

#Dividing out the outliers
{```{R}```
OF.out <- read.table("Outflank_Fst_Outliers.list", head=F)
BS.out <- read.table("Bayescan/Bayescan_Outliers.list", head=F)

gen.out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out)), env.loc))
out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out))))
length(gen.out.list)
#898
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% gen.out.list))]

gen.net <- gen3[, loc=set.loc]
gen.net.vcf <- VCF_remove(gen3.vcf, gen.out.list)
nLoc(gen.net)
#6973
dim(Loci_names(locNames(gen.net.vcf), SEP="_"))
#6973

gen.out <- gen3[, loc=out.list]
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% out.list))]
gen.out.vcf <- VCF_remove(gen3.vcf, set.loc)
nLoc(gen.out)
#28
dim(Loci_names(locNames(gen.out.vcf), SEP="_"))
#28

gen.env <- gen3[, loc=env.loc]
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% env.loc))]
gen.env.vcf <- VCF_remove(gen3.vcf, set.loc)
nLoc(gen.env)
#873
dim(Loci_names(locNames(gen.env.vcf), SEP="_"))
#873

save(gen.net, file="gen.net_noPac.gz", compress=T)
save(gen.net.vcf, file="gen.net.vcf_noPac.gz", compress=T)
save(gen.out, file="gen.out_noPac.gz", compress=T)
save(gen.out.vcf, file="gen.out.vcf_noPac.gz", compress=T)
save(gen.env, file="gen.env_noPac.gz", compress=T)
save(gen.env.vcf, file="gen.env.vcf_noPac.gz", compress=T)

write.table(locNames(gen.net.vcf), file="gen.net_SNPs.txt", col.names=F, row.names=F, quote=F)
write.table(locNames(gen.out.vcf), file="gen.out_SNPs.txt", col.names=F, row.names=F, quote=F)
write.table(locNames(gen.env.vcf), file="gen.env_SNPs.txt", col.names=F, row.names=F, quote=F)
write.table(gen.net@strata$INDV[gen.net@strata$G3 == "Atlantic"], file="Atlantic_indv.txt", col.names=F, row.names=F, quote=F)
write.table(gen.net@strata$INDV[gen.net@strata$G3 == "East"], file="East_indv.txt", col.names=F, row.names=F, quote=F)
write.table(gen.net@strata$INDV[gen.net@strata$G3 == "West"], file="West_indv.txt", col.names=F, row.names=F, quote=F)
} #notepad cleanup

#Getting the neutral loci with the Pacific sample for Arlequin analysis
{```{R}```
load("gen2_wPac.gz")
load("gen3_noPac.gz")
load("gen.net_noPac.gz")
gen3P <- gen2[, loc=locNames(gen3)]
G3 <- as.character(as.matrix(gen.net@strata$G3))
names(G3) <- gen.net@strata$INDV

SST.loc <- read.table("SST_loci.txt", head=F)[,1]
Nitrate.loc <- read.table("Nitrate_loci.txt", head=F)[,1]
Nitrate_ss.loc <- read.table("Nitrate_ss_loci.txt", head=F)[,1]
env.loc <- unique(c(as.character(as.matrix(SST.loc)), as.character(as.matrix(Nitrate> oc)), as.character(as.matrix(Nitrate_ss.loc))))

OF.out <- read.table("Outflank_Fst_Outliers.list", head=F)
BS.out <- read.table("Bayescan/Bayescan_Outliers.list", head=F)
gen.out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out)), env.loc))
out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out))))
length(gen.out.list)
#898

set.loc <- locNames(gen3P)[which(!(locNames(gen3P) %in% gen.out.list))]
gen.net <- gen3P[, loc=set.loc]
gen.net@strata$INDV[grep("Pac",gen.net@strata$INDV)]
#[1] Pac.AS.T437_Lib2_I12_D01
G3["Pac.AS.T437_Lib2_I12_D01"] <- "Pacific"
gen.net@strata$G3 <- G3[match(gen.net@strata$INDV,names(G3))]

c4 <- c("red4", "mediumblue", "grey40", "chocolate2")

X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.net) <- ~G3
par(mfrow=c(1,1))
ade4::s.class(pca.net$li, pop(gen.net), col=c4, cstar=0, axesell=F, clabel=0)

setPop(gen.net)<-~G3
writeGenPop(gen.net, "Neutral_G3_groups_Pac.gen", "Neutral loci in Angel shark data without dups by G3 including Pacific")

all_pop <- seppop(gen.net)

# Function to modify the numbers in each vector
modify_numbers <- function(vec, VALUE=1) {
  modified_vec <- sub("^\\d", VALUE, vec)
  return(modified_vec)
}

for(i in 1:length(all_pop)){
modified_list <- lapply(all_pop[[i]]@all.names, function(x) modify_numbers(x,i))
all_pop[[i]]@all.names <- modified_list
}

#repool the geninds back together
gen2.net <- repool(all_pop)

setPop(gen2.net)<-~G3
writeGenPop(gen2.net, "Neutral_G3_groups_Pac_Fmax.gen", "Neutral loci in Angel shark data by G3 with Pacific and unique alleles")
} #notepad cleanup
{```{bash}```
java -jar /usr/local/bin/PGDSpider2.1.1.3-cli.jar -inputfile Neutral_G3_groups_Pac.gen -inputformat GENEPOP -outputfile Neutral_G3_groups_Pac.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java -jar /usr/local/bin/PGDSpider2.1.1.3-cli.jar -inputfile Neutral_G3_groups_Pac_Fmax.gen -inputformat GENEPOP -outputfile Neutral_G3_groups_Pac_Fmax.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
} #notepad cleanup
#Exporting data for Arlequin
{```{R}```
setPop(gen.net)<-~G3
writeGenPop(gen.net, "Neutral_G3_groups.gen", "Neutral loci in Angel shark data without dups by G3")
setPop(gen.out)<-~G3
writeGenPop(gen.out, "Outlier_G3_groups.gen", "Outlier loci in Angel shark data without dups by G3")
setPop(gen.env)<-~G3
writeGenPop(gen.env, "Environ_G3_groups.gen", "Environmental loci in Angel shark data without dups by G3")
} #notepad cleanup
{```{bash}```
java -jar /usr/local/bin/PGDSpider2.1.1.3-cli.jar -inputfile Neutral_G3_groups.gen -inputformat GENEPOP -outputfile Neutral_G3_groups.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java -jar /usr/local/bin/PGDSpider2.1.1.3-cli.jar -inputfile Outlier_G3_groups.gen -inputformat GENEPOP -outputfile Outlier_G3_groups.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
java -jar /usr/local/bin/PGDSpider2.1.1.3-cli.jar -inputfile Environ_G3_groups.gen -inputformat GENEPOP -outputfile Environ_G3_groups.arp -outputformat Arlequin -spid ~/bin/genepop_to_arlequin_STD.spid
} #notepad cleanup

#PCA of datasets
{```{R}```
#Standard
X.net <- scaleGen(gen.net, NA.method="mean", scale=F)
pca.net <- dudi.pca(X.net,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.out <- scaleGen(gen.out, NA.method="mean", scale=F)
pca.out <- dudi.pca(X.out,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.env <- scaleGen(gen.env, NA.method="mean", scale=F)
pca.env <- dudi.pca(X.env,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen.net) <- setPop(gen.out) <- setPop(gen.env) <- ~G3
par(mfrow=c(3,1))
ade4::s.class(pca.net$li, pop(gen.net), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("Neutral loci", adj=0.05)

pca.net$li[grep(64, rownames(pca.net$li)),1:2]

ade4::s.class(pca.out$li, pop(gen.out), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("Outlier loci", adj=0.05)

ade4::s.class(pca.env$li, pop(gen.env), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
legend(8, -2, legend=c(levels(pop(gen.env))), pch=19, col=c3[c(2,1,3)], bty="n")
mtext("Environmental loci", adj=0.05)

tiff("PCA_G3_split_loci_v2.tif", res=300, height =6000, width=2000)
par(mfrow=c(3,1))
ade4::s.class(pca.net$li, pop(gen.net), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca.net$eig[1]/sum(pca.net$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca.net$eig[2]/sum(pca.net$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("Neutral loci", adj=0.05)

ade4::s.class(pca.out$li, pop(gen.out), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca.out$eig[1]/sum(pca.out$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca.out$eig[2]/sum(pca.out$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("Outlier loci", adj=0.05)

ade4::s.class(pca.env$li, pop(gen.env), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca.env$eig[1]/sum(pca.env$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca.env$eig[2]/sum(pca.env$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
legend(8, -2, legend=c(levels(pop(gen.env))), pch=19, col=c3[c(2,1,3)], bty="n")
mtext("Environmental loci", adj=0.05)
legend(2, -5, legend=c(levels(pop(gen.env))), pch=19, col=c3[c(2,1,3)], bty="n")
dev.off()


setPop(gen.net) <- ~net_G3
tiff("PCA_G3_neutral_loci.tif", res=300, height =6000, width=2000)
ade4::s.class(pca.net$li, pop(gen.net), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca.net$eig[1]/sum(pca.net$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca.net$eig[2]/sum(pca.net$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext(paste("Neutral loci (n=",nLoc(gen.net),")",sep=""), adj=0.05, line=2)
legend(-15,-7.5, legend=c("Atlantic (n=17)", "East GOM (n=25)", "West GOM (n=20)"), pch=19, col=c3[c(2,1,3)], bty="n")
dev.off()

tiff("PCA_G3_outlier_loci.tif", res=300, height =6000, width=2000)
setPop(gen.out) <- ~out_G3
ade4::s.class(pca.out$li, pop(gen.out), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
mtext("Outlier loci", adj=0.05)
dev.off()

tiff("PCA_G3_env_loci.tif", res=300, height =6000, width=2000)
setPop(gen.env) <- ~env_G3
ade4::s.class(pca.env$li, pop(gen.env), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("x-axis variation: ",format(round(100*pca1$eig[1]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=0.8, adj=0.15, cex=1.5)
mtext(paste("y-axis variation: ",format(round(100*pca1$eig[2]/sum(pca1$eig),3),nsmall=3),"%",sep=""),1, line=2, adj=0.15, cex=1.5)
legend(8, -2, legend=c(levels(pop(gen.env))), pch=19, col=c3[c(2,1,3)], bty="n")
mtext("Environmental loci", adj=0.05)
dev.off()

pca.out$li[pca.out$li[,1] < 1 & rownames(pca.out$li) %in% gen.out@strata$INDV[gen.out@strata$G3 == "West"],1:2]
} #notepad cleanup

#Grouping Neutral data
{```{R}```
grp <- find.clusters(gen.net,  max.n.clust=20, n.pca=400, method = "kmeans")
#2
grp$Kstat
{ #Results
     K=1      K=2      K=3      K=4      K=5      K=6      K=7      K=8
432.7147 435.2906 438.2236 441.1931 444.1411 447.0897 450.0160 452.9500
     K=9     K=10     K=11     K=12     K=13     K=14     K=15     K=16
455.8222 458.7262 461.5647 464.3969 467.1998 469.9835 472.7311 475.4876
    K=17     K=18     K=19     K=20
478.1784 480.8600 483.4990 486.1163
}

#Getting clustering
for(i in 2:5){assign(paste("grp",i,sep=""), find.clusters(gen.net, max.n.clust=40, n.pca=400, n.clust =i, method="kmeans"))}

#Putting assignments into strata
for(i in 2:5){
tmp.grp <- get(paste("grp",i,sep=""))
tmp.v <- data.frame(tmp.grp$grp)[match(gen.net@strata$INDV,names(tmp.grp$grp)),]
gen.net@strata[,paste("net_G",i,sep="")] <- as.numeric(tmp.v)
}

#scaling genind objects
X <- scaleGen(gen.net, NA.method="mean", scale=F)

#DAPC cross validation
par(mfrow=c(4,1))
xval.2 <- xvalDapc(X, gen.net@strata$net_G2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.3 <- xvalDapc(X, gen.net@strata$net_G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.4 <- xvalDapc(X, gen.net@strata$net_G4, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.5 <- xvalDapc(X, gen.net@strata$net_G5, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)

xval.G3 <- xvalDapc(X, gen.net@strata$G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.G3[2:6]	#Success: 1.0000000
{ #Results
$`Median and Confidence Interval for Random Chance`
     2.5%       50%     97.5%
0.2167454 0.3223623 0.4450047

$`Mean Successful Assignment by Number of PCs of PCA`
        5        10        15        20        25        30        35        40
0.9900000 1.0000000 0.9966667 1.0000000 0.9766667 0.9766667 0.9666667 0.9633333
       45        50
0.9233333 0.8000000

$`Number of PCs Achieving Highest Mean Success`
[1] "10"

$`Root Mean Squared Error by Number of PCs of PCA`
         5         10         15         20         25         30         35
0.04082483 0.00000000 0.02357023 0.00000000 0.06236096 0.06236096 0.08164966
        40         45         50
0.08498366 0.11303883 0.26246693

$`Number of PCs Achieving Lowest MSE`
[1] "20"
}

#Outputting the crossvalidation data
results.df <- data.frame(matrix(ncol=5))
colnames(results.df) <- c("K", "Success", "Success_PCs", "RMSE", "RMSE_PCs")
for(i in 2:5){
xval.tmp <- get(paste("xval",i, sep="."))
tmp.row <- c(i, max((xval.tmp[[3]])), xval.tmp[[4]], min((xval.tmp[[5]])), xval.tmp[[6]])
results.df <- rbind(results.df, tmp.row)
}
results.df <- results.df[-1,]
results.df
{ #Results
  K           Success Success_PCs               RMSE RMSE_PCs
2 2             0.985          15 0.0612372435695794       15
3 3 0.986666666666667          15 0.0471404520791032       15
4 4 0.964444444444444           5 0.0812707713240433       15
5 5 0.910833333333333          10  0.133982793256779       10
}
table(gen.net@strata$G3, gen.net@strata$net_G2)
{ #Results
            1  2
  Atlantic 17  0
  East      1 20
  West      0 24
}
table(gen.net@strata$G3, gen.net@strata$net_G3)
{ #Results
            1  2  3
  Atlantic  0  0 17
  East      5 16  0
  West     20  4  0
}
save(gen.net, file="gen.net_noPac.gz", compress=T)

setPop(gen.net) <- ~G3
p2.comp <- ggcompoplot(xval.2$DAPC, gen.net, pal = funky(2), cols=2) + theme(axis.text.x = element_blank()) + ggtitle("K-means 2 Groups")
p3.comp <- ggcompoplot(xval.3$DAPC, gen.net, pal = funky(3), cols=2) + theme(axis.text.x = element_blank()) + ggtitle("K-means 3 Groups")
p4.comp <- ggcompoplot(xval.4$DAPC, gen.net, pal = funky(4), cols=2) + theme(axis.text.x = element_blank()) + ggtitle("K-means 4 Groups")
} #notepad cleanup

#Identifying loci Driving Neutral structure
{```{R}```
#scaling genind objects
X <- scaleGen(gen.net, NA.method="mean", scale=F)
xval.G3 <- xvalDapc(X, gen.net@strata$G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
xval.G3[2:6]	#Success: 1.0000000
{ #Results
$`Median and Confidence Interval for Random Chance`
     2.5%       50%     97.5%
0.2167454 0.3223623 0.4450047

$`Mean Successful Assignment by Number of PCs of PCA`
        5        10        15        20        25        30        35        40
0.9900000 1.0000000 0.9966667 1.0000000 0.9766667 0.9766667 0.9666667 0.9633333
       45        50
0.9233333 0.8000000

$`Number of PCs Achieving Highest Mean Success`
[1] "10"

$`Root Mean Squared Error by Number of PCs of PCA`
         5         10         15         20         25         30         35
0.04082483 0.00000000 0.02357023 0.00000000 0.06236096 0.06236096 0.08164966
        40         45         50
0.08498366 0.11303883 0.26246693

$`Number of PCs Achieving Lowest MSE`
[1] "20"
}

par(mfrow=c(1,1))
scatter(xval.G3$DAPC, cstar=0, col=c3[c(2,1,3)], posi.da="bottomright")

#Old result
#Neutral_select <- Locus_hunt(gen.net, xval.G3$DAPC, THRES = 0.95, GROUP="~G3", AXIS=1, ALPHA=0.05, MIN=0.8)
#print(Neutral_select[[3]], row.names = FALSE)
{ #Results 
                                                 Summary
 The Best Threshold tested was 0.88671875
 This removed 2230 loci (28.4%)
 Kept p-value is 0.06
 Removed p-value is 0.01
 Data has been output to the 1st two objects of this list
 Item 1 contains the kept genind object
 Item 2 contains the removed genind object
}

Neutral_select <- Locus_hunt(gen.net, xval.G3$DAPC, THRES = 0.95, GROUP="~G3", AXIS=1, ALPHA=0.05, MIN=0.8, ABS=0.00001, MODEL="normal")
print(Neutral_select[[3]], row.names = FALSE)
{ #Results
                                                  Summary
 The Best Threshold tested was 0.93153076171875
 This removed 1234 loci (17.7%)
 Kept p-value is 0.09
 Removed p-value is 0.01
 Data has been output to the 1st two objects of this list
 Item 1 contains the kept genind object
 Item 2 contains the removed genind object
}

xval.G2 <- xvalDapc(X, gen.net@strata$G2, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
par(mfrow=c(1,1))
scatter(xval.G2$DAPC, cstar=0, col=c2, posi.da="bottomright")

Neutral2_select <- Locus_hunt(gen.net, xval.G2$DAPC, THRES = 0.95, GROUP="~G2", AXIS=1, ALPHA=0.05, MIN=0.8)
print(Neutral2_select[[3]], row.names = FALSE)
{ #Results
                                                  Summary
 The Best Threshold tested was 0.94296875
 This removed 1051 loci (15%)
 Kept p-value is 0.08
 Removed p-value is 0.01
 Data has been output to the 1st two objects of this list
 Item 1 contains the kept genind object
 Item 2 contains the removed genind object
}

Neutral2_select <- Locus_hunt(gen.net, xval.G2$DAPC, THRES = 0.95, GROUP="~G2", AXIS=1, ALPHA=0.05, MIN=0.94, ABS=0.00001, MODEL="normal")
print(Neutral2_select[[3]], row.names = FALSE)
{ #Results
 The Best Threshold tested was 0.9425
 This removed 1061 loci (15.2%)
 Kept p-value is 0.09
 Removed p-value is 0.01
 Data has been output to the 1st two objects of this list
 Item 1 contains the kept genind object
 Item 2 contains the removed genind object
}

Neutral2_select <- Locus_hunt(gen.net, xval.G2$DAPC, THRES = 0.95, GROUP="~G2", AXIS=1, ALPHA=0.05, MIN=0.94, ABS=0.00001, MODEL="normal")
print(Neutral2_select[[3]], row.names = FALSE)
{ #Results
                                                  Summary
 The Best Threshold tested was 0.9442578125
 This removed 1026 loci (14.7%)
 Kept p-value is 0.05
 Removed p-value is 0.01
 Data has been output to the 1st two objects of this list
 Item 1 contains the kept genind object
 Item 2 contains the removed genind object
}

save(Neutral2_select, file="Gul_Atl_loci.RData.gz", compress=T)

X.Net1 <- scaleGen(Neutral2_select[[1]], NA.method="mean", scale=F)
pca.Net1 <- dudi.pca(X.Net1,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.Net2 <- scaleGen(Neutral2_select[[2]], NA.method="mean", scale=F)
pca.Net2 <- dudi.pca(X.Net2,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

png("Gulf_Atlantic_loci.png", res=200, height=4000, width=2000)
setPop(Neutral2_select[[1]]) <- setPop(Neutral2_select[[2]]) <- ~G3
par(mfrow=c(2,1))
ade4::s.class(pca.Net1$li, pop(Neutral2_select[[1]]), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("Kept loci ( n =", nLoc(Neutral2_select[[1]]),")"), adj=0.05)

ade4::s.class(pca.Net2$li, pop(Neutral2_select[[2]]), col=c3[c(2,1,3)], cstar=0, axesell=F, clabel=0)
mtext(paste("Removed loci ( n =", nLoc(Neutral2_select[[2]]),")"), adj=0.05)
dev.off()

setPop(gen.net) <- ~G3
temp <- seppop(gen.net)
east <- temp$`East`
west <- temp$`West`
pooled_gens <- repool(east,west)

X <- scaleGen(pooled_gens, NA.method="mean", scale=F)
xval.Gulf <- xvalDapc(X, pooled_gens@strata$G3, n.pca.max = 300, training.set = 0.9, result = "groupMean", center = TRUE, scale = FALSE, n.pca = NULL, n.rep = 50, xval.plot = TRUE)
scatter(xval.Gulf$DAPC, cstar=0, col=c2, posi.da="bottomright")

Gulf_select <- Locus_hunt(pooled_gens, xval.Gulf$DAPC, THRES = 0.95, GROUP="~G3", AXIS=1, ALPHA=0.05, MIN=0.8)
print(Gulf_select[[3]], row.names = FALSE)
{ #Results
                                                  Summary
 The Best Threshold tested was 0.9953125
 This removed 93 loci (1.3%)
 Kept p-value is 0.07
 Removed p-value is 0.01
 Data has been output to the 1st two objects of this list
 Item 1 contains the kept genind object
 Item 2 contains the removed genind object
}

save(Gulf_select, file="Gul_split_loci.RData.gz", compress=T)

X.Gulf1 <- scaleGen(Gulf_select[[1]], NA.method="mean", scale=F)
pca.Gulf1 <- dudi.pca(X.Gulf1,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

X.Gulf2 <- scaleGen(Gulf_select[[2]], NA.method="mean", scale=F)
pca.Gulf2 <- dudi.pca(X.Gulf2,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(Gulf_select[[1]]) <- setPop(Gulf_select[[2]]) <- ~G3
par(mfrow=c(2,1))
ade4::s.class(pca.Gulf1$li, pop(Gulf_select[[1]]), col=c2, cstar=0, axesell=F, clabel=0)
mtext("Kept loci", adj=0.05)

ade4::s.class(pca.Gulf2$li, pop(Gulf_select[[2]]), col=c2, cstar=0, axesell=F, clabel=0)
mtext("Removed loci", adj=0.05)


length(which(locNames(Gulf_select[[2]]) %in% locNames(Neutral2_select[[2]])))
#32
loc_overlap <- locNames(Gulf_select[[2]])[which(locNames(Gulf_select[[2]]) %in% locNames(Neutral2_select[[2]]))]
set.loci <- locNames(gen.net)[!locNames(gen.net) %in% loc_overlap]
gen.tmp <- gen.net[,loc=set.loci]

tmp <- poppr.amova(gen.net, ~G3, quiet=T)
randtest(tmp, nrepet = 10000)
{ #Results
class: krandtest lightkrandtest
Monte-Carlo tests
Call: randtest.amova(xtest = tmp, nrepet = 10000)

Number of tests:   3

Adjustment method for multiple comparisons:   none
Permutation number:   10000
                        Test        Obs   Std.Obs   Alter     Pvalue
1  Variations within samples 1729.44981 -2.277862    less 0.01109889
2 Variations between samples   43.51087  1.779529 greater 0.03809619
3      Variations between G3   14.30693 19.489061 greater 0.00009999
}
tmp2 <- poppr.amova(gen.tmp, ~G3, quiet=T)
randtest(tmp2, nrepet = 10000)
{ #Results
class: krandtest lightkrandtest
Monte-Carlo tests
Call: randtest.amova(xtest = tmp2, nrepet = 10000)

Number of tests:   3

Adjustment method for multiple comparisons:   none
Permutation number:   10000
                        Test        Obs   Std.Obs   Alter     Pvalue
1  Variations within samples 1716.22007 -2.271795    less 0.01129887
2 Variations between samples   43.34852  1.800639 greater 0.03289671
3      Variations between G3   13.36880 19.039696 greater 0.00009999
}
} #notepad cleanup

#New Diversity Analysis
{```{R}```

#Preparing the final data frame for Neutral data
#Allelic Richness
{```{R}```
tmp.v <- gen.net@strata$G3
tmp.v[tmp.v == "1"] <- "Atlantic"
tmp.v[tmp.v == "2"] <- "East_Gulf"
tmp.v[tmp.v == "3"] <- "West_Gulf"
gen.net@strata$G3 <- as.factor(tmp.v)

setPop(gen.net) <- ~G3
Ar_data.df <- allelic.richness(gen.net)$Ar


#Expected Heterozygosity
He_data.df <- read.table("Neutral_He.csv", head=T, sep="\t")
He_data.df <- He_data.df[! He_data.df[,7] == 0.00000,]
He_data.df <- He_data.df[,2:4]

colnames(He_data.df)[1] <- "Atlantic"
colnames(He_data.df)[2] <- "East_Gulf"
colnames(He_data.df)[3] <- "West_Gulf"
} #notepad cleanup

#Running the Friedmans on Neutral data
{```{R}```
for(i in 1:2){
VAR <- c("Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

#Removing rows with NA
na.test <- apply(data.df, 1, function(x) sum(is.na(x)))
data.df <- data.df[which(na.test == 0),]

#"Gathering" data
data.m <- matrix(ncol=3)
for(j in 1:length(locNames(gen.net)[which(na.test == 0)])){
tmp.m <- matrix(c(rep(locNames(gen.net)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP2","VAR")

#Friedman's test
tmp.dat <- friedman.test(VAR ~ POP2 | Locus, data=test.df)
print(paste("Friedman's test on", VAR, "( n=", nrow(data.df),")"))
print(tmp.dat)
print(apply(get(paste(VAR, "_data.df", sep="")), 2, function(x) mean(x, na.rm=T)))
}
{ #Results
[1] "Friedman's test on Ar ( n= 6973 )"

        Friedman rank sum test

data:  VAR and POP2 and Locus
Friedman chi-squared = 144.56, df = 2, p-value < 2.2e-16

 Atlantic West_Gulf East_Gulf
 2.567658  2.705672  2.622757
[1] "Friedman's test on He ( n= 6970 )"

        Friedman rank sum test

data:  VAR and POP2 and Locus
Friedman chi-squared = 69.591, df = 2, p-value = 7.736e-16

 Atlantic East_Gulf West_Gulf
0.2870911 0.2950072 0.2879220
}
} #notepad cleanup
#While Ar is out of order, the file and result have been confirmed

#Saving image
{```{R}```
save.image("Fried_G3_net.RData.gz", compress=T)
#load("Fried_G3_net.RData.gz")
}

#Running the Wilcoxon Test
{```{R}```
error.list <- NULL

for(i in 1:2){
VAR <- c("Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

data.m <- matrix(ncol=3)
for(j in 1:length(locNames(gen.net))){
tmp.m <- matrix(c(rep(locNames(gen.net)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP2","VAR")

test.df$POP2 <- factor(test.df$POP2, levels=c(c("Atlantic","West_Gulf","East_Gulf")))
test.df$VAR <- as.numeric(as.matrix(test.df$VAR))
result_m <- t(combn(levels(test.df$POP2), 2))
result_m <- cbind(result_m, matrix(ncol=3, nrow=nrow(result_m)))

for(k in 1:nrow(result_m)){
tmp.df <- test.df[test.df$POP2 %in% result_m[k,1:2],]
tmp.df$POP2 <- gdata::drop.levels(tmp.df$POP2)
if(length(table(tmp.df$POP2)) != 2){next}
rm.loci <- as.character(as.matrix(unique(tmp.df$Locus[is.na(tmp.df$VAR)])))
tmp.df <- tmp.df[!tmp.df$Locus %in% rm.loci,]
tmp.df$Locus <- gdata::drop.levels(tmp.df$Locus)
possibleError <- tryCatch(tmp.out <- wilcoxsign_test(VAR ~ POP2 | Locus, tmp.df), error=function(e) e)
if(inherits(possibleError, "error")){error.list <- c(error.list,c(VAR,result_m[k,1:2])); next}

tmp.df$POP2 <- factor(tmp.df$POP2, levels=c(result_m[k,1], result_m[k,2]))
tmp.out <- wilcoxsign_test(VAR ~ POP2 | Locus,tmp.df)
tmp.stat <- statistic(tmp.out)
tmp.p <- pvalue(tmp.out)
result_m[k,3:4] <- c(tmp.stat, tmp.p)

tmp.td <- tidyr::spread(tmp.df, POP2, VAR)
tmp.td$diff <- tmp.td[,2] - tmp.td[,3]
tmp.td$rank <- rank(abs(tmp.td$diff))
tmp.td$sgn <- tmp.td$diff
tmp.td$sgn[tmp.td$sgn > 0] <- 1
tmp.td$sgn[tmp.td$sgn < 0] <- -1

result_m[k,5] <- sum(tmp.td$sgn*tmp.td$rank)
}

result_df <- as.data.frame(result_m)
colnames(result_df) <- c("Loc1", "Loc2", "stat", "pvalue", "T_stat")
for(j in 3:5){result_df[,j] <- as.numeric(as.matrix(result_df[,j]))}
result_df$pvalue_adj <- p.adjust(result_df$pvalue, method="fdr")
print(paste("Wilcox test on", VAR))
print(paste("Number of significant results:", length(which(result_df$pvalue < 0.05))))
print(paste("Number of significant results adj:", length(which(result_df$pvalue_adj < 0.05))))
#print(result_df[result_df$pvalue < 0.05])
assign(paste(VAR,"_wilcox.df",sep=""), result_df)
}
{ #Results
[1] "Wilcox test on Ar"
[1] "Number of significant results: 3"
[1] "Number of significant results adj: 3"
[1] "Wilcox test on He"
[1] "Number of significant results: 2"
[1] "Number of significant results adj: 2"

for(i in 1:6){
VAR <- c("Nuc_div","Hap_div","Seg_sites","theta.s", "Ar", "He")[i]
print(VAR)
print(get(paste(VAR,"_wilcox.df",sep="")))
}
{ #Results
[1] "Nuc_div"
       Loc1      Loc2       stat       pvalue   T_stat   pvalue_adj
1  Atlantic West_Gulf   3.226107 1.254865e-03  1084656 1.254865e-03
2  Atlantic East_Gulf  -5.016826 5.253213e-07 -1686725 7.879820e-07
3 West_Gulf East_Gulf -10.492484 9.353651e-26 -3527705 2.806095e-25
[1] "Hap_div"
       Loc1      Loc2      stat       pvalue   T_stat   pvalue_adj
1  Atlantic West_Gulf -6.161071 7.225470e-10 -2071425 1.083820e-09
2  Atlantic East_Gulf -1.322448 1.860191e-01  -444622 1.860191e-01
3 West_Gulf East_Gulf  6.261501 3.812894e-10  2105191 1.083820e-09
[1] "Seg_sites"
       Loc1      Loc2       stat        pvalue   T_stat    pvalue_adj
1  Atlantic West_Gulf -25.271806 6.523189e-141 -8111386 1.956957e-140
2  Atlantic East_Gulf -20.880199  8.105309e-97 -6664952  1.215796e-96
3 West_Gulf East_Gulf   6.746751  1.511924e-11  2069457  1.511924e-11
[1] "theta.s"
       Loc1      Loc2       stat       pvalue   T_stat   pvalue_adj
1  Atlantic West_Gulf -5.5579009 2.730383e-08 -1868038 4.095574e-08
2  Atlantic East_Gulf  0.1036074 9.174809e-01    34828 9.174809e-01
3 West_Gulf East_Gulf 31.8251086 0.000000e+00 10697493 0.000000e+00
[1] "Ar"
       Loc1      Loc2       stat       pvalue   T_stat   pvalue_adj
1  Atlantic West_Gulf -13.403247 5.787467e-41 -4506327 8.681200e-41
2  Atlantic East_Gulf  -5.482747 4.187723e-08 -1843361 4.187723e-08
3 West_Gulf East_Gulf  11.314872 0.000000e+00  3804191 0.000000e+00
[1] "He"
       Loc1      Loc2      stat       pvalue   T_stat   pvalue_adj
1  Atlantic West_Gulf -1.136332 2.558177e-01  -381801 2.558177e-01
2  Atlantic East_Gulf -5.969294 2.382827e-09 -2005651 3.574241e-09
3 West_Gulf East_Gulf -6.409130 1.463527e-10 -2153433 4.390582e-10
}

}}}}#clean-up notepad

#Saving image
{```{R}```
save.image("Fried_G3_net.RData.gz", compress=T)
#load("Fried_G3_net.RData.gz")
}

#Preparing the data tables (Outlier data)
{```{R}```
#Allelic Richness
tmp.v <- GENIND@strata$G3
tmp.v[tmp.v == "1"] <- "Atlantic"
tmp.v[tmp.v == "2"] <- "East_Gulf"
tmp.v[tmp.v == "3"] <- "West_Gulf"
GENIND@strata$G3 <- as.factor(tmp.v)

setPop(GENIND) <- ~G3
Ar_data.df <- allelic.richness(GENIND)$Ar

#Expected Heterozygosity
He_data.df <- read.table("Outlier_He.txt", head=T, sep="\t")
He_data.df <- He_data.df[! He_data.df[,7] == 0.00000,]
He_data.df <- He_data.df[,2:4]

colnames(He_data.df)[1] <- "Atlantic"
colnames(He_data.df)[2] <- "East_Gulf"
colnames(He_data.df)[3] <- "West_Gulf"
} #notepad cleanup

#Running the Friedmans on Outlier data
{```{R}```
for(i in 1:2){
VAR <- c("Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

#Removing rows with NA
na.test <- apply(data.df, 1, function(x) sum(is.na(x)))
data.df <- data.df[which(na.test == 0),]

#"Gathering" data
data.m <- matrix(ncol=3)
for(j in 1:length(locNames(GENIND)[which(na.test == 0)])){
tmp.m <- matrix(c(rep(locNames(GENIND)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP2","VAR")

#Friedman's test
tmp.dat <- friedman.test(VAR ~ POP2 | Locus, data=test.df)
print(paste("Friedman's test on", VAR, "( n=", nrow(data.df),")"))
print(tmp.dat)
print(apply(get(paste(VAR, "_data.df", sep="")), 2, function(x) mean(x, na.rm=T)))
}
{ #Results
[1] "Friedman's test on Ar ( n= 28 )"

        Friedman rank sum test

data:  VAR and POP2 and Locus
Friedman chi-squared = 3.7818, df = 2, p-value = 0.1509

 Atlantic West_Gulf East_Gulf
 3.140421  3.210856  2.775940
[1] "Friedman's test on He ( n= 28 )"

        Friedman rank sum test

data:  VAR and POP2 and Locus
Friedman chi-squared = 10.571, df = 2, p-value = 0.005063

 Atlantic East_Gulf West_Gulf
0.4690057 0.4345143 0.3035414
}
} #notepad cleanup
#While Ar is out of order, the file and result have been confirmed

#Saving image
{```{R}```
save.image("Fried_G3_out.RData.gz", compress=T)
#load("Fried_G3_out.RData.gz")
}

#Running the Wilcox Test
{```{R}```
error.list <- NULL

for(i in 1:2){
VAR <- c("Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

data.m <- matrix(ncol=3)
for(j in 1:length(locNames(GENIND))){
tmp.m <- matrix(c(rep(locNames(GENIND)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP2","VAR")

test.df$POP2 <- factor(test.df$POP2, levels=c(c("Atlantic","West_Gulf","East_Gulf")))
test.df$VAR <- as.numeric(as.matrix(test.df$VAR))
result_m <- t(combn(levels(test.df$POP2), 2))
result_m <- cbind(result_m, matrix(ncol=3, nrow=nrow(result_m)))

for(k in 1:nrow(result_m)){
tmp.df <- test.df[test.df$POP2 %in% result_m[k,1:2],]
tmp.df$POP2 <- gdata::drop.levels(tmp.df$POP2)
if(length(table(tmp.df$POP2)) != 2){next}
rm.loci <- as.character(as.matrix(unique(tmp.df$Locus[is.na(tmp.df$VAR)])))
tmp.df <- tmp.df[!tmp.df$Locus %in% rm.loci,]
tmp.df$Locus <- gdata::drop.levels(tmp.df$Locus)
possibleError <- tryCatch(tmp.out <- wilcoxsign_test(VAR ~ POP2 | Locus, tmp.df), error=function(e) e)
if(inherits(possibleError, "error")){error.list <- c(error.list,c(VAR,result_m[k,1:2])); next}

tmp.df$POP2 <- factor(tmp.df$POP2, levels=c(result_m[k,1], result_m[k,2]))
tmp.out <- wilcoxsign_test(VAR ~ POP2 | Locus,tmp.df)
tmp.stat <- statistic(tmp.out)
tmp.p <- pvalue(tmp.out)
result_m[k,3:4] <- c(tmp.stat, tmp.p)

tmp.td <- tidyr::spread(tmp.df, POP2, VAR)
tmp.td$diff <- tmp.td[,2] - tmp.td[,3]
tmp.td$rank <- rank(abs(tmp.td$diff))
tmp.td$sgn <- tmp.td$diff
tmp.td$sgn[tmp.td$sgn > 0] <- 1
tmp.td$sgn[tmp.td$sgn < 0] <- -1

result_m[k,5] <- sum(tmp.td$sgn*tmp.td$rank)
}

result_df <- as.data.frame(result_m)
colnames(result_df) <- c("Loc1", "Loc2", "stat", "pvalue", "T_stat")
for(j in 3:5){result_df[,j] <- as.numeric(as.matrix(result_df[,j]))}
result_df$pvalue_adj <- p.adjust(result_df$pvalue, method="fdr")
print(paste("Wilcox test on", VAR))
print(paste("Number of significant results:", length(which(result_df$pvalue < 0.05))))
print(paste("Number of significant results adj:", length(which(result_df$pvalue_adj < 0.05))))
#print(result_df[result_df$pvalue < 0.05])
assign(paste(VAR,"_wilcox.df",sep=""), result_df)
}
{ #Results
[1] "Wilcox test on Ar"
[1] "Number of significant results: 1"
[1] "Number of significant results adj: 0"
[1] "Wilcox test on He"
[1] "Number of significant results: 2"
[1] "Number of significant results adj: 2"
}

for(i in 1:6){
VAR <- c("Ar", "He")[i]
print(VAR)
print(get(paste(VAR,"_wilcox.df",sep="")))
}
{ #Results
[1] "Ar"
       Loc1      Loc2      stat     pvalue T_stat pvalue_adj
1  Atlantic West_Gulf 0.6831427 0.49451667     60 0.49451667
2  Atlantic East_Gulf 1.6514635 0.09864397    145 0.14796595
3 West_Gulf East_Gulf 2.2315995 0.02564145    196 0.07692434
[1] "He"
       Loc1      Loc2       stat       pvalue T_stat   pvalue_adj
1  Atlantic West_Gulf  2.8008851 0.0050962664    246 0.0076443996
2  Atlantic East_Gulf  0.7059141 0.4802415427     62 0.4802415427
3 West_Gulf East_Gulf -3.5978849 0.0003208155   -316 0.0009624466
}

}}}} #clean-up notepad

#Saving image
{```{R}```
save.image("Fried_G3_out.RData.gz", compress=T)
#load("Fried_G3_out.RData.gz")
}

#Preparing the Environmental data tables
{```{R}```
#Allelic Richness
tmp.v <- GENIND@strata$G3
tmp.v[tmp.v == "1"] <- "Atlantic"
tmp.v[tmp.v == "2"] <- "East_Gulf"
tmp.v[tmp.v == "3"] <- "West_Gulf"
GENIND@strata$G3 <- as.factor(tmp.v)

setPop(GENIND) <- ~G3
Ar_data.df <- allelic.richness(GENIND)$Ar

#Expected Heterozygosity
He_data.df <- read.table("Environ_He.txt", head=T, sep="\t")
He_data.df <- He_data.df[! He_data.df[,7] == 0.00000,]
He_data.df <- He_data.df[,2:4]

colnames(He_data.df)[1] <- "Atlantic"
colnames(He_data.df)[2] <- "East_Gulf"
colnames(He_data.df)[3] <- "West_Gulf"
} #Notepad cleanup

#Running the Friedmans on Environmental data
{```{R}```
for(i in 1:2){
VAR <- c("Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

#Removing rows with NA
na.test <- apply(data.df, 1, function(x) sum(is.na(x)))
data.df <- data.df[which(na.test == 0),]

#"Gathering" data
data.m <- matrix(ncol=3)
for(j in 1:length(locNames(GENIND)[which(na.test == 0)])){
tmp.m <- matrix(c(rep(locNames(GENIND)[j],ncol(danvta.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP2","VAR")

#Friedman's test
tmp.dat <- friedman.test(VAR ~ POP2 | Locus, data=test.df)
print(paste("Friedman's test on", VAR, "( n=", nrow(data.df),")"))
print(tmp.dat)
print(apply(get(paste(VAR, "_data.df", sep="")), 2, function(x) mean(x, na.rm=T)))
}
{ #Results
[1] "Friedman's test on Ar ( n= 873 )"

        Friedman rank sum test

data:  VAR and POP2 and Locus
Friedman chi-squared = 9.5518, df = 2, p-value = 0.008431

 Atlantic West_Gulf East_Gulf
 3.493793  3.622685  3.510904
[1] "Friedman's test on He ( n= 873 )"

        Friedman rank sum test

data:  VAR and POP2 and Locus
Friedman chi-squared = 5.4176, df = 2, p-value = 0.06662

 Atlantic East_Gulf West_Gulf
0.5098087 0.5041053 0.4940786
}
} #notepad cleanup
#While Ar is out of order, the file and result have been confirmed

#Saving image
{```{R}```
save.image("Fried_G3_env.RData.gz", compress=T)
#load("Fried_G3_env.RData.gz")
}

#Running the Wilcox Test
{```{R}```
error.list <- NULL

for(i in 1:2){
VAR <- c("Ar", "He")[i]
data.df <- get(paste(VAR, "_data.df", sep=""))

data.m <- matrix(ncol=3)
for(j in 1:length(locNames(GENIND))){
tmp.m <- matrix(c(rep(locNames(GENIND)[j],ncol(data.df)), colnames(data.df), as.numeric(data.df[j,])),ncol=3,byrow=F)
data.m <- rbind(data.m, tmp.m)
}
test.df <- as.data.frame(data.m[-1, ])
colnames(test.df) <- c("Locus","POP2","VAR")

test.df$POP2 <- factor(test.df$POP2, levels=c(c("Atlantic","West_Gulf","East_Gulf")))
test.df$VAR <- as.numeric(as.matrix(test.df$VAR))
result_m <- t(combn(levels(test.df$POP2), 2))
result_m <- cbind(result_m, matrix(ncol=3, nrow=nrow(result_m)))

for(k in 1:nrow(result_m)){
tmp.df <- test.df[test.df$POP2 %in% result_m[k,1:2],]
tmp.df$POP2 <- gdata::drop.levels(tmp.df$POP2)
if(length(table(tmp.df$POP2)) != 2){next}
rm.loci <- as.character(as.matrix(unique(tmp.df$Locus[is.na(tmp.df$VAR)])))
tmp.df <- tmp.df[!tmp.df$Locus %in% rm.loci,]
tmp.df$Locus <- gdata::drop.levels(tmp.df$Locus)
possibleError <- tryCatch(tmp.out <- wilcoxsign_test(VAR ~ POP2 | Locus, tmp.df), error=function(e) e)
if(inherits(possibleError, "error")){error.list <- c(error.list,c(VAR,result_m[k,1:2])); next}

tmp.df$POP2 <- factor(tmp.df$POP2, levels=c(result_m[k,1], result_m[k,2]))
tmp.out <- wilcoxsign_test(VAR ~ POP2 | Locus,tmp.df)
tmp.stat <- statistic(tmp.out)
tmp.p <- pvalue(tmp.out)
result_m[k,3:4] <- c(tmp.stat, tmp.p)

tmp.td <- tidyr::spread(tmp.df, POP2, VAR)
tmp.td$diff <- tmp.td[,2] - tmp.td[,3]
tmp.td$rank <- rank(abs(tmp.td$diff))
tmp.td$sgn <- tmp.td$diff
tmp.td$sgn[tmp.td$sgn > 0] <- 1
tmp.td$sgn[tmp.td$sgn < 0] <- -1

result_m[k,5] <- sum(tmp.td$sgn*tmp.td$rank)
}

result_df <- as.data.frame(result_m)
colnames(result_df) <- c("Loc1", "Loc2", "stat", "pvalue", "T_stat")
for(j in 3:5){result_df[,j] <- as.numeric(as.matrix(result_df[,j]))}
result_df$pvalue_adj <- p.adjust(result_df$pvalue, method="fdr")
print(paste("Wilcox test on", VAR))
print(paste("Number of significant results:", length(which(result_df$pvalue < 0.05))))
print(paste("Number of significant results adj:", length(which(result_df$pvalue_adj < 0.05))))
#print(result_df[result_df$pvalue < 0.05])
assign(paste(VAR,"_wilcox.df",sep=""), result_df)
}
{ #Results
[1] "Wilcox test on Ar"
[1] "Number of significant results: 2"
[1] "Number of significant results adj: 2"
[1] "Wilcox test on He"
[1] "Number of significant results: 2"
[1] "Number of significant results adj: 2"
}

for(i in 1:2){
VAR <- c("Ar", "He")[i]
print(VAR)
print(get(paste(VAR,"_wilcox.df",sep="")))
}
{ #Results
[1] "Ar"
       Loc1      Loc2        stat       pvalue T_stat   pvalue_adj
1  Atlantic West_Gulf -3.00481920 0.0026573876 -44787 0.0039860815
2  Atlantic East_Gulf -0.07125098 0.9431980060  -1062 0.9431980060
3 West_Gulf East_Gulf  3.85781813 0.0001144037  57501 0.0003432112
[1] "He"
       Loc1      Loc2       stat     pvalue T_stat pvalue_adj
1  Atlantic West_Gulf  2.3132416 0.02070936  34479 0.03106404
2  Atlantic East_Gulf  0.9436394 0.34535394  14065 0.34535394
3 West_Gulf East_Gulf -2.5436333 0.01097062 -37913 0.03106404
}

}}}} #clean-up notepad

#Saving image
{```{R}```
save.image("Fried_G3_env.RData.gz", compress=T)
#load("Fried_G3_env.RData.gz")
}
}

## Moments analysis ##
#Preparing vcf file 
{```{R}```
loci_names <- Loci_names(locNames(gen.net.vcf), SEP="_")

rand_loci <- data.frame(matrix(ncol=1,nrow=length(loci_names)))
for(i in 1:nrow(loci_names)){
set.loc <- locNames(gen.net.vcf)[grep(paste(as.character(loci_names[i,]),"_",sep=""), locNames(gen.net.vcf))]
tmp.gen <- gen.net.vcf[,loc=set.loc]
rand_loci[i,] <- sample(locNames(tmp.gen),1)
if(i/500==round(i/500,0)){print(paste("iteration", i, "finished"))}
}

write.table(rand_loci, "Random_loci_for_dadi.txt", col.names=F, row.names=F, quote=F)
#rand_loci <- read.table("Random_loci_for_dadi.txt", head=F)
write.table(indNames(gen.net.vcf), "sample_for_dadi.txt", col.names=F, row.names=F, quote=F)
write.table(gen.net@strata[,c("INDV","G3")], "G3_strata.txt", col.names=F, row.names=F, quote=F, sep="\t")
}} #notepad cleanup
{```{bash}```
sed -i 's/_/\t/3' Random_loci_for_dadi.txt
vcftools --vcf SNP.TRS.F07.recode.vcf --out SNP.TRS.F07_DADI --recode --recode-INFO-all --positions Random_loci_for_dadi.txt --keep sample_for_dadi.txt

python3 ~/bin/easySFS.py -i SNP.TRS.F07_DADI.recode.vcf -p G3_strata.txt --preview
{ #Results
Atlantic
(2, 898.0)      (3, 1345.0)     (4, 1679.0)     (5, 1956.0)     (6, 2195.0)     (7, 2406.0)     (8, 2595.0)     (9, 2766.0)     (10, 2923.0)    (11, 3066.0)    (12, 3199.0)    (13, 3322.0)    (14, 3436.0)   (15, 3543.0)     (16, 3644.0)    (17, 3738.0)    (18, 3827.0)    (19, 3911.0)    (20, 3990.0)    (21, 4065.0)    (22, 4136.0)    (23, 4203.0)    (24, 4267.0)    (25, 4322.0)    (26, 4380.0)    (27, 4428.0)   (28, 4481.0)     (29, 4496.0)    (30, 4544.0)    (31, 4419.0)    (32, 4462.0)    (33, 3605.0)    (34, 3637.0)

East
(2, 931.0)      (3, 1393.0)     (4, 1743.0)     (5, 2036.0)     (6, 2291.0)     (7, 2519.0)     (8, 2725.0)     (9, 2912.0)     (10, 3085.0)    (11, 3245.0)    (12, 3394.0)    (13, 3533.0)    (14, 3663.0)   (15, 3786.0)     (16, 3902.0)    (17, 4012.0)    (18, 4116.0)    (19, 4215.0)    (20, 4309.0)    (21, 4399.0)    (22, 4485.0)    (23, 4567.0)    (24, 4645.0)    (25, 4721.0)    (26, 4793.0)    (27, 4862.0)   (28, 4929.0)     (29, 4992.0)    (30, 5054.0)    (31, 5111.0)    (32, 5168.0)    (33, 5214.0)    (34, 5268.0)    (35, 5302.0)    (36, 5351.0)    (37, 5340.0)    (38, 5386.0)    (39, 5155.0)    (40, 5196.0)   (41, 3955.0)     (42, 3984.0)

West
(2, 898.0)      (3, 1346.0)     (4, 1681.0)     (5, 1958.0)     (6, 2199.0)     (7, 2413.0)     (8, 2606.0)     (9, 2782.0)     (10, 2943.0)    (11, 3092.0)    (12, 3230.0)    (13, 3359.0)    (14, 3480.0)   (15, 3593.0)     (16, 3700.0)    (17, 3801.0)    (18, 3896.0)    (19, 3987.0)    (20, 4073.0)    (21, 4155.0)    (22, 4233.0)    (23, 4307.0)    (24, 4379.0)    (25, 4447.0)    (26, 4512.0)    (27, 4574.0)   (28, 4634.0)     (29, 4691.0)    (30, 4747.0)    (31, 4800.0)    (32, 4851.0)    (33, 4900.0)    (34, 4947.0)    (35, 4991.0)    (36, 5034.0)    (37, 5068.0)    (38, 5109.0)    (39, 5132.0)    (40, 5170.0)   (41, 5164.0)     (42, 5199.0)    (43, 5083.0)    (44, 5115.0)    (45, 4524.0)    (46, 4550.0)    (47, 2909.0)    (48, 2925.0)
}

python3 ~/bin/easySFS.py -f -i SNP.TRS.F07_DADI.recode.vcf -o output_AEW -p AEW_popmap --proj 32,38,42
python3 ~/bin/easySFS.py -f -i SNP.TRS.F07_DADI.recode.vcf -o output_WEA -p WEA_popmap --proj 32,38,42
}

############## Making Pretty figures ##############
# Figure 1
{```{R}```
load("DAPC.RData.gz")

library('sp')
library('rnaturalearth')
library('dichromat')
library('ggtern')
library('prettymapr')

#rivers10 <- ne_load(scale=10,type = "rivers_lake_centerlines", category = "physical", destdir="~/Workspace/Spinner/analysis")

for(i in c("Lat", "Lon")){gen3@strata[,i] <- as.numeric(as.matrix(gen3@strata[,i]))}
pt.col <- as.character(as.matrix(gen3@strata$G3))
pt.col[pt.col == 1] <- "mediumblue"
pt.col[pt.col == 2] <- "green"
pt.col[pt.col == 3] <- "chocolate2"

select_states <- c("Texas", "Louisiana", "Alabama", "Mississippi", "Florida", "Georgia", "South Carolina", "North Carolina", "Virginia", "Maryland")
tmp.state <- ne_states(geounit="United States of America")
state.df <- data.frame(matrix(ncol=4, nrow=length(select_states)))
names(state.df) <- c("State","abbr","Lat","Lon")
state.df$State <- select_states
for(i in select_states){
state.df$abbr[state.df$State == i] <- tmp.state$postal[tmp.state$name == i]
state.df$Lat[state.df$State == i] <- tmp.state$latitude[tmp.state$name == i]
state.df$Lon[state.df$State == i] <- st_centroid(tmp.state[tmp.state$name == i])
}

plot(st_geometry(st_centroid(nc)), pch = 3, col = 'red', add = TRUE)

<- ne_states(geounit="United States of America")


png("Fig1_May2024.png", res=600, width=6000, height=6000)
par(mfrow=c(1,1))
plot(st_geometry(ne_countries(type = "countries", scale = "large")[,1]), xlim=c(-100, -74), ylim=c(24,38.3), bg=scales::alpha("grey60",0.5), col="white", axes=T, main="")
sp::plot(ne_states(geounit="United States of America"), add=T, col="white")
points(gen3@strata$Lon, gen3@strata$Lat, pch=21, cex=1, bg=pt.col)
legend(-79.8, 30, legend=c("Atlantic Ocean (n=17)", "eastern Gulf (n=21)", "western Gulf (n=24)"), title="Samples", pch=21, pt.bg=c("chocolate2","green","mediumblue"), bg="white")
text(-90, 24.5, "Gulf of\nMexico")
text(-77, 31.5, "Atlantic\nOcean")
addnortharrow(pos = "topleft", scale=0.75)
for(i in select_states){if(i == "Louisiana"){text(st_coordinates(st_centroid(tmp.state[tmp.state$name ==i,])), label=state.df$abbr[state.df$State == i], adj=1)
}else if(i == "Florida"){text(st_coordinates(st_centroid(tmp.state[tmp.state$name ==i,])), label=state.df$abbr[state.df$State == i], adj=-0.2)}
else{text(st_coordinates(st_centroid(tmp.state[tmp.state$name ==i,])), label=state.df$abbr[state.df$State == i], adj=0.5)}}
box()
dev.off()

library('terra')
pt.col <- as.character(as.matrix(gen3@strata$G3))
pt.col[pt.col == 3] <- "mediumblue"
pt.col[pt.col == 1] <- "green"
pt.col[pt.col == 2] <- "chocolate2"

png("Fig1v2_May2024.png", res=600, width=6000, height=4000)
par(mfrow=c(1,1))
plot(st_geometry(ne_countries(type = "countries", scale = "large")[,1]), xlim=c(-100, -74), ylim=c(24,38.3), bg=scales::alpha("grey60",0.5), col="white", axes=T, main="")
sp::plot(ne_states(geounit="United States of America"), add=T, col="white")
points(gen3@strata$Lon, gen3@strata$Lat, pch=21, cex=1, bg=pt.col)
legend(-78.5, 27, legend=c("Atlantic Ocean (n=17)", "eastern Gulf (n=21)", "western Gulf (n=24)"), title="Samples", pch=21, pt.bg=c("chocolate2","green","mediumblue"), bg="white")
text(-90, 24.5, "Gulf of\nMexico")
text(-77, 31.5, "Atlantic\nOcean")
north(xy=c(-100, 38.5), label="", type=2, adj=0.75)
text(-100, 37.5, "N")
text(st_coordinates(st_centroid(tmp.state[tmp.state$name == "Texas",])), label="Texas", adj=0.5)
tmp.cood <- st_coordinates(st_centroid(tmp.state[tmp.state$name == "Louisiana",]))
text(tmp.cood[1], 30.5, label="Louisiana", adj=0.5)
text(-82, 29, label="Florida", srt=290, adj=0)
tmp.cood <- st_coordinates(st_centroid(tmp.state[tmp.state$name == "Mississippi",]))
text(st_coordinates(st_centroid(tmp.state[tmp.state$name == "Mississippi",])), label="Mississippi", srt=270, adj=0.5)
text(st_coordinates(st_centroid(tmp.state[tmp.state$name == "Alabama",])), label="Alabama", srt=270, adj=0.5)
text(st_coordinates(st_centroid(tmp.state[tmp.state$name == "Georgia",])), label="Georgia", adj=0.5)
tmp.cood <- st_coordinates(st_centroid(tmp.state[tmp.state$name == "South Carolina",]))
text(tmp.cood[1], tmp.cood[2]*1.01, label="South", adj=0.5)
text(tmp.cood[1], tmp.cood[2]*0.995, label="Carolina", adj=0.5)
tmp.cood <- st_coordinates(st_centroid(tmp.state[tmp.state$name == "North Carolina",]))
text(tmp.cood[1], tmp.cood[2]*1.005, label="North Carolina", adj=0.5)
tmp.cood <- st_coordinates(st_centroid(tmp.state[tmp.state$name == "Virginia",]))
text(tmp.cood[1], tmp.cood[2]*0.99, label="Virginia", adj=0.5)
text(-100.5, 25, "Mexico")
box()
dev.off()
}}} #notepad cleanup

# Figure 2
{```{R}```
load("DAPC.RData.gz")

X <- scaleGen(gen3, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen3) <- ~G3

col_pts <- as.character(as.matrix(gen3@strata$G3))
col_pts[col_pts == 3] <- "chocolate2"
col_pts[col_pts == 2] <- "green"
col_pts[col_pts == 1] <- "mediumblue"

png("PCA_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca1$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dev.off()
}

# Figure 3
{```{R}```
#Getting datasets
{```{R}```
load("DAPC.RData.gz")
SST.loc <- read.table("SST_loci.txt", head=F)[,1]
Nitrate.loc <- read.table("Nitrate_loci.txt", head=F)[,1]
Nitrate_ss.loc <- read.table("Nitrate_ss_loci.txt", head=F)[,1]

env.loc <- unique(c(as.character(as.matrix(SST.loc)), as.character(as.matrix(Nitrate.loc)), as.character(as.matrix(Nitrate_ss.loc))))
length(env.loc)
#[1] 873

OF.out <- read.table("Outflank_Fst_Outliers.list", head=F)
BS.out <- read.table("Bayescan/Bayescan_Outliers.list", head=F)

gen.out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out)), env.loc))
out.list <- unique(c(as.character(as.matrix(OF.out)),as.character(as.matrix(BS.out))))
length(gen.out.list)
#898
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% gen.out.list))]

gen.net <- gen3[, loc=set.loc]
nLoc(gen.net)
#6973

gen.out <- gen3[, loc=out.list]
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% out.list))]
nLoc(gen.out)
#28

gen.env <- gen3[, loc=env.loc]
set.loc <- locNames(gen3)[which(!(locNames(gen3) %in% env.loc))]
nLoc(gen.env)
#873
}
#Making plots
{```{R}```
X <- scaleGen(gen.net, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.net) <- ~G3

col_pts <- as.character(as.matrix(gen.net@strata$G3))
col_pts[col_pts == 3] <- "chocolate2"
col_pts[col_pts == 2] <- "green"
col_pts[col_pts == 1] <- "mediumblue"

png("PCA_neutral_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca1$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dev.off()

X <- scaleGen(gen.out, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.out) <- ~G3

png("PCA_outlier_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca1$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dev.off()

X <- scaleGen(gen.env, NA.method="mean", scale=F)
pca1 <- dudi.pca(X,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
setPop(gen.env) <- ~G3

png("PCA_envir_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca1$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca1$eig[1]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca1$eig[2]/sum(pca1$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dev.off()
}
}

# Figure S3
{```{R}```
library('plotrix')
load("RDA_geo_env.gz")

col_pts <- as.character(as.matrix(gen3@strata$G3))
col_pts[col_pts == 3] <- "chocolate2"
col_pts[col_pts == 2] <- "green"
col_pts[col_pts == 1] <- "mediumblue"

#Axes 1 & 2
MIN.x <- min(c(m.ord_All$CCA$wa[,1]*14.5, m.ord_All$CCA$biplot[,1]*2))*1.01
MAX.x <- max(c(m.ord_All$CCA$wa[,1]*14.5, m.ord_All$CCA$biplot[,1]*2))*1.01
MIN.y <- min(c(m.ord_All$CCA$wa[,2]*14.5, m.ord_All$CCA$biplot[,2]*2))*1.01
MAX.y <- max(c(m.ord_All$CCA$wa[,2]*14.5, m.ord_All$CCA$biplot[,2]*2))*1.01

png("RDA_A12_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(m.ord_All$CCA$wa*14.5, pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), xlab="RDA Axis 1", ylab="RDA Axis 2", cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
plotrix::boxed.labels(m.ord_All$CCA$biplot[1,1]*1.8, m.ord_All$CCA$biplot[1,2]*3.5, rownames(m.ord_All$CCA$biplot)[1], cex = 1, border = NA, bg ="white", xpad = 1.1, ypad = 1.4)
plotrix::boxed.labels(m.ord_All$CCA$biplot[2,1]*0.7, m.ord_All$CCA$biplot[2,2]*1.3, rownames(m.ord_All$CCA$biplot)[2], cex = 1, border = NA, bg ="white", xpad = 1.1, ypad = 1.4)
plotrix::boxed.labels(m.ord_All$CCA$biplot[3,1]*2.0, m.ord_All$CCA$biplot[3,2]*2.0, rownames(m.ord_All$CCA$biplot)[3], cex = 1, border = NA, bg ="white", xpad = 1.1, ypad = 1.4)
for(i in 1:nrow(m.ord_All$CCA$biplot)){arrows(0,0,m.ord_All$CCA$biplot[i,1]*1.5, m.ord_All$CCA$biplot[i,2]*1.5, col="black")}
dev.off()

png("RDA_A13_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(m.ord_All$CCA$wa[,c(1,3)]*14.5, pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlim=c(MIN.x, MAX.x), ylim=c(MIN.y, MAX.y), xlab="RDA Axis 1", ylab="RDA Axis 3", cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
plotrix::boxed.labels(0.8, -0.25, rownames(m.ord_All$CCA$biplot)[1], cex = 1, border = NA, bg ="white", xpad = 1.1, ypad = 1.4)
plotrix::boxed.labels(0.45, 0.2, rownames(m.ord_All$CCA$biplot)[2], cex = 1, border = NA, bg ="white", xpad = 1.1, ypad = 1.4)
plotrix::boxed.labels(m.ord_All$CCA$biplot[3,1]*2.0, m.ord_All$CCA$biplot[3,3]*1.6, rownames(m.ord_All$CCA$biplot)[3], cex = 1, border = NA, bg ="white", xpad = 1.1, ypad = 1.4)
for(i in 1:nrow(m.ord_All$CCA$biplot)){arrows(0,0,m.ord_All$CCA$biplot[i,1]*1.5, m.ord_All$CCA$biplot[i,3]*1.5, col="black")}
dev.off()
}

# Figure S4
{```{R}```
gen.SST <- gen3[,loc=SST.loc]
gen.Nit <- gen3[,loc=Nitrate.loc]
gen.Nit_ss <- gen3[,loc=Nitrate_ss.loc]

X.SST <- scaleGen(gen.SST, NA.method="mean", scale=F)
X.Nit <- scaleGen(gen.Nit, NA.method="mean", scale=F)
X.Nit_ss <- scaleGen(gen.Nit_ss, NA.method="mean", scale=F)

pca.SST <- dudi.pca(X.SST,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
pca.Nit <- dudi.pca(X.Nit,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)
pca.Nit_ss <- dudi.pca(X.Nit_ss,cent=FALSE,scale=FALSE,scannf=FALSE,nf=5000)

setPop(gen.SST) <- ~G3
setPop(gen.Nit) <- ~G3
setPop(gen.Nit_ss) <- ~G3

col_pts <- as.character(as.matrix(gen3@strata$G3))
col_pts[col_pts == "Atlantic"] <- "chocolate2"
col_pts[col_pts == "East"] <- "green"
col_pts[col_pts == "West"] <- "mediumblue"

png("PCA_SST_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca.SST$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca.SST$eig[1]/sum(pca.SST$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca.SST$eig[2]/sum(pca.SST$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dataEllipse(as.matrix(pca.SST$li[gen.SST@strata$G3 == "Atlantic", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2", lwd=0.5)
dataEllipse(as.matrix(pca.SST$li[gen.SST@strata$G3 == "East", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green", lwd=0.5)
dataEllipse(as.matrix(pca.SST$li[gen.SST@strata$G3 == "West", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue", lwd=0.5)
dev.off()

png("PCA_Nit_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca.Nit$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca.Nit$eig[1]/sum(pca.Nit$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca.Nit$eig[2]/sum(pca.Nit$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dataEllipse(as.matrix(pca.Nit$li[gen.Nit@strata$G3 == "Atlantic", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2", lwd=0.5, xpd=F)
dataEllipse(as.matrix(pca.Nit$li[gen.Nit@strata$G3 == "East", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green", lwd=0.5, xpd=F)
dataEllipse(as.matrix(pca.Nit$li[gen.Nit@strata$G3 == "West", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue", lwd=0.5, xpd=F)
dev.off()

png("PCA_Nit_ss_May2024.png", res=600, width=6000, height=6000)
par(mar=c(8.1, 4.6, 4.1, 2.1), xpd=T)
plot(pca.Nit_ss$li[,1:2], pch=21, col="black", bg=col_pts, lwd=0.5, cex=1.5, xlab=paste("PC 1 (variation = ",format(round(100*pca.Nit_ss$eig[1]/sum(pca.Nit_ss$eig),2),nsmall=3),"%)",sep=""), 
ylab=paste("PC 2 (variation = ",format(round(100*pca.Nit_ss$eig[2]/sum(pca.Nit_ss$eig),2),nsmall=3),"%)",sep=""), cex.lab=1.5)
abline(v=0,h=0,col="black", lty=2, xpd=F)
legend("bottom", inset=c(0,-0.17), x.intersp=1, xjust=0, legend=c("Atlantic (n=17)","east Gulf (n=21)","west Gulf (n=24)"), pt.cex=1.25, pch=21, col="black", pt.bg=c("chocolate2","green","mediumblue"), bty="n", cex=1, ncol=3)
dataEllipse(as.matrix(pca.Nit_ss$li[gen.Nit_ss@strata$G3 == "Atlantic", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="chocolate2", lwd=0.5, xpd=F)
dataEllipse(as.matrix(pca.Nit_ss$li[gen.Nit_ss@strata$G3 == "East", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="green", lwd=0.5, xpd=F)
dataEllipse(as.matrix(pca.Nit_ss$li[gen.Nit_ss@strata$G3 == "West", 1:2]), levels=c(0.75), center.pch=F, add=T, plot.points=F, col="mediumblue", lwd=0.5, xpd=F)
dev.off()
}



