##########################################################################################
# R SCRIPT FOR MORPHOMETRIC ANALYSES (having snakes in mind...)
# version 2021.06.11 (reduced from Barbo et al., 2021 - Systematics and Biodiversity)
# Designed to analyze data for Barbo et al., 2021 (Canadian Journal of Zoology)
##########################################################################################
# Felipe G. Grazziotin - fgrazziotin@gmail.com
##########################################################################################

# "All organisms contain an array of biological information, some of which is useful in a comparative or taxonomic sense, some of which is useless, and some of which may even be misleading."
# Judith E. Winston, 1999 (Describing Species, p. 53)

# "... Molecular Biology is beginning to be seen as a restricted science, narrowing our vision... the need for taxonomists to draw attention to the enormous diversity and variation of this earth's biota becomes more and more pressing."
# Tod. F. Stuessy, 1990 (Plant Taxonomy, p. xvii)

# "Is it not extraordinary that young taxonomists are trained like performing monkeys, almost wholly by imitation, and that in only the rarest cases are they given any instruction in taxonomic theory?
# A. J. Cain, 1959 (quoted in Simpson, 1961 - Principles of Animal Taxonomy, p. vii)

# Make sure you have identified all of the possible confounding variables in your study... always work hard digging deep trying to identify them...

# "Check the distribution, check the outliers, check the residuals."

# "In God we trust, all others bring data."
# William Edwards Deming (1900-1993) as quoted in Hastie, Tibshirani, Friedman (2017) - The Elements of Statistical Learning, pg. vii.

##########################################################################################
# General script to analyze snake morphological data 
# To configure your dataset, please follow these procedures:
# 1 - Have a column named "Sex" with variables "F", "M" or "" (for sex not defined);
# 2 - Have one or more columns named as "OTUs", "OTUs01", "OTUs02", etc. In these columns you should include your preliminary taxonomic hypotheses like: "lineage01", "spnov01", etc;
# 3 - Geographic coordinates should be in decimal degrees and split in two columns named "Lat" for latitude and "Long" for longitude;
# 4 - Name all your morphological variables with the prefixes: "cont_" for continuous variables, "cat_" for categorical variables, rat_ for ratios between variables, and "disc_" for discrete variables;
# 5 - Use the following strict convention to name the following variables: "SVL" for snout vent length, "TaL" for tail length, "HL" for head length, "HW" for head length, ...
# 6 - Have a column named "ID" with a specific label for each row.

# X - Your table needs to be in csv format.

# OBS01 - I like loops... I know it is not "R style" and can perform badly when compared to apply functions... feel free to apply "apply" fuctions when they can be applied...
# OBS02 - Please, let me know when you find errors or make substantial improvements.
# OBS03 - This is a modular general script combining original functions and tricks for data handling with several other commands derived from blogs, forums, books and regular papers. I deeply acknowledge all clever brains that originally created several pieces of the following code.
# OBS04 - Please, if you are using this script or parts of it, cite in your publication the packages and their associated references

##########################################################################################

#######################################

##########################################################################################
# INSTALL/LOAD PACKAGES
##########################################################################################

#######################################

# install.packages("nortest")
# needed for lillie.test
library(nortest)

# install.packages("MASS")
# needed for LDA and MDS
# library(MASS)

# install.packages("corrplot")
# needed for corrplot
# library(corrplot)

# install.packages("coin")
# needed for independence_test
# library(coin)

# install.packages("wPerm")
# needed for perm.ind.test
# library(wPerm)

# install.packages("caret")
# needed for PLSDA and RF
# library(caret)

# install.packages("devtools")
# needed for to install ggbiplot
# library(devtools)

# install.packages("ggplot2")
# needed for all plots
library(ggplot2)

# install.packages("gridExtra")
# needed to plot in grid
library(gridExtra)

# install.packages("rfPermute")
# library(rfPermute)

# install.packages("psych")
# needed for some summary statistics
# library(psych)

# install.packages("vegan")
# needed for Permanova
# library(vegan)

# install.packages("RVAideMemoire")
# needed for pairwise Post-hoc Permanova and PLSDA.VIP
# library(RVAideMemoire)

# install.packages("doSNOW")
# needed to run CV in parallel
# library(doSNOW)


#######################################

##########################################################################################
# READ AND SORT DATA 
##########################################################################################

#######################################

# read .csv file
x<-read.csv("dataset.csv",header=T)

#######################

# read OTUs and cat variables as factor
for(i in seq(grep("OTUs",colnames(x)))){
	x[,grep("OTUs",colnames(x))[i]]<-as.factor(x[,grep("OTUs",colnames(x))[i]])
}
for(i in seq(grep("cat",colnames(x)))){
	x[,grep("cat",colnames(x))[i]]<-as.factor(x[,grep("cat",colnames(x))[i]])
}

#######################

# create a backup for data
x_backup<-x

# x<-x_backup

#######################

# sort variables
x_disc<-x[,grep("disc",colnames(x))]
x_cat<-x[,grep("cat",colnames(x))]
# transform x_cat in numeric
for(i in seq(x_cat[1,])) {
	x_cat[,i]<-as.numeric(x_cat[,i])
}
x_cont<-x[,grep("cont",colnames(x))]
x_rate<-x[,grep("rat",colnames(x))]

#######################

# recombine variables in different dataframes
x_all<-cbind(x_disc,x_rate,x_cat,x_cont)
x_multi<-cbind(x_disc,x_rate,x_cat)

#######################

# set hypotheses of OTU for clustering tests
Hindex<-grep("OTUs",colnames(x))
Hlabel<-colnames(x)[grep("OTUs",colnames(x))]
OTUs<-list()
for (k in seq(Hlabel)) {
	OTUs[Hlabel[k]]<-list(levels(x[,Hindex[k]]))
}

#######################

##########################################################################################
# CHECK, TRANSFORM AND MANIPULATE DATA 
##########################################################################################

#######################

# check subsets
str(x_disc)
str(x_cat)
str(x_cont)
str(x_rate)
str(x_all)
str(x_multi)

#######################################

# create directories for results
dir.create("./01_tables",showWarnings = FALSE)
dir.create("./02_histograms",showWarnings = FALSE)
dir.create("./03_QQplots",showWarnings = FALSE)
dir.create("./04_boxplots",showWarnings = FALSE)
# dir.create("./05_correlograms",showWarnings = FALSE)
# dir.create("./06_PCAs",showWarnings = FALSE)
# dir.create("./07_PLSDAs",showWarnings = FALSE)
# dir.create("./08_RF",showWarnings = FALSE)
dir.create("./09_Combined",showWarnings = FALSE)
#######################################

##########################################################################################
# HISTOGRAMS TO CHECK FOR ERRORS 
##########################################################################################

#######################################

# Plot histograms for ALL variables NO groups
# OBS: Only to check for erros, weird distributions and outliers
dir.create("./02_histograms/01_all_variables-no_groups",showWarnings = FALSE)

# histogram for ALL variables
for(i in seq(x_all[1,])){
	pdf(paste("./02_histograms/01_all_variables-no_groups/histo_",colnames(x_all)[i],".pdf",sep=""), family="ArialMT", width=9, height=12)
		a<-hist(x_all[,i],freq=F,xlab=NULL,,main=colnames(x_all)[i])
	dev.off()
}

#######################################

# Plot histograms for ALL variables grouped by SEX
# OBS: Only to check for erros and weird distributions
dir.create("./02_histograms/02_all_variables-sex",showWarnings = FALSE)

#histogram for ALL variables
for(i in seq(x_all[1,])){
	pdf(paste("./02_histograms/02_all_variables-sex/histo_",colnames(x_all)[i],".pdf",sep=""), family="ArialMT", width=9, height=12)
	layout(matrix(c(1,2),2,1))
		a<-hist(x_all[,i][x$Sex=="F"],freq=F,xlab="F",,main=colnames(x_all)[i],xlim=c(min(na.omit(x_all[,i])),max(na.omit(x_all[,i]))))
		a<-hist(x_all[,i][x$Sex=="M"],freq=F,xlab="M",,main=colnames(x_all)[i],xlim=c(min(na.omit(x_all[,i])),max(na.omit(x_all[,i]))))
	dev.off()
}

#######################################

##########################################################################################
# BOXPLOT + HISTOGRAM TO CHECK FOR ERRORS - first attempt
##########################################################################################

#######################################

# Plot Boxplot + histogram for ALL variables - complete dataset
# Only to check for outliers

dir.create("./04_boxplots/01_all_variables-no-groups-check")

for(i in 1:ncol(x_all)){
	#get histogram data
	b1<-hist(x_all[,i],breaks=length(x_all[,i]),plot=FALSE)
	b2<-hist(x_all[,i],plot=FALSE)
	#open pdf file
	pdf(paste("./04_boxplots/01_all_variables-no-groups-check/boxplot_",colnames(x_all)[i],".pdf",sep=""),family="ArialMT", width=9, height=12)
		#remove zeros from counts
		q1<-b1$counts
		q1[q1==0]<-NA
		q2<-b2$counts
		q2[q2==0]<-NA
		#transform counts to boxplot coordinates
		m1<-(0.8-1.2)/(max(na.omit(q1))-min(na.omit(q1)))
		y1<-(m1*q1)+(1.2-(m1*min(na.omit(q1))))
		m2<-(0.8-1.2)/(max(na.omit(q2))-min(na.omit(q2)))
		y2<-(m2*q2)+(1.2-(m2*min(na.omit(q2))))
		#plot boxplot
		boxplot(x_all[,i],main=colnames(x_all)[i],outline=F,ylim=c(min(na.omit(x_all[,i])),max(na.omit(x_all[,i]))),lwd=1.5)
		#plot points from histogram
		points(y1,b1$mids,col="blue",cex=0.8)
		#plot bars from histogram
		rect(y2,b2$breaks[-length(b2$breaks)],rep(1.2,length(y2)),b2$breaks[-1],border="blue",col="blue", lty=3,lwd=0.3,density = 10)
		#plot mean as a line
		segments(0.8,mean(na.omit(x_all[,i])),1.2,mean(na.omit(x_all[,i])),col="red")
		#plot N for variable
		text(1,max(na.omit(x_all[,i])), labels=paste("N = ",length(na.omit(x_all[,i])),sep=""),pos=3)
		#get range amplitude
		k<-max(na.omit(x_all[,i]))-min(na.omit(x_all[,i]))
		#plot axis for histogram
		segments(0.8,min(na.omit(x_all[,i]))-(0.01*k),1.2,min(na.omit(x_all[,i]))-(0.01*k),col="grey",adj=0)
		#plot max line for histogram
		segments(0.8,min(na.omit(x_all[,i])),0.8,max(na.omit(x_all[,i])),lty=3,col="grey")
		segments(0.8,min(na.omit(x_all[,i]))-(0.02*k),0.8,min(na.omit(x_all[,i])),col="grey")
		#plot min line for histogram
		segments(1.2,min(na.omit(x_all[,i])),1.2,max(na.omit(x_all[,i])),lty=3,col="grey")
		segments(1.2,min(na.omit(x_all[,i]))-(0.02*k),1.2,min(na.omit(x_all[,i])),col="grey")
		#plot labels for the histogram axis
		text(1,min(na.omit(x_all[,i]))-(0.025*k),labels="Count",cex=0.8)
		text(0.8,min(na.omit(x_all[,i]))-(0.03*k),labels=max(na.omit(q1)),cex=0.8)
		text(1.2,min(na.omit(x_all[,i]))-(0.03*k),labels=min(na.omit(q1)),cex=0.8)
		
	dev.off()
}

#######################################

##########################################################################################
# CLEANING BAD VARIABLES
##########################################################################################

#######################################

# adjust dataset by removing non-informative variables
# OBS: nonrandom missing data (e.g. missing variable for all individuals from an OTU) can have pervasive effects in multivariate analyses

# remove all variables that are completely missing for an OTU
to_remove<-list()
list_chars<-vector()
# loop for each taxon set
for(k in seq(Hlabel)){
	# loop for each variable
	for (i in 1:ncol(x_all)){
		for (j in seq(OTUs[[k]])) {
			if(length(na.omit(x_all[,i][x[,Hindex[k]]==OTUs[[k]][j]]))<1){
				list_chars<-c(list_chars,colnames(x_all[i]))
				to_remove[[k]]<-unique(list_chars)
			}
		}
	}
}

# write list of missing variables
write.table(to_remove,file=paste("./01_tables/Removed_Characters.txt",sep=""),sep="\t",col.names=F)

# Check the differences among taxon sets...
# If the missing variables are the same, follow the script. Otherwise, choose the variables to remove.
# Ex:
# to_remove2<-c("cat_LLd","cat_LLe",cat_con_fx_posoc")

to_remove2<-c()

to_remove2<-c(to_remove2,unique(unlist(to_remove)))

# cleaning datasets
x_cln0<-x[,!colnames(x)%in%to_remove2]
x_all_cln0<-x_all[,!colnames(x_all)%in%to_remove2]
x_multi_cln0<-x_multi[,!colnames(x_multi)%in%to_remove2]

#######################################

##########################################################################################
# BOXPLOT + HISTOGRAM TO CHECK FOR ERRORS - second attempt
##########################################################################################

#######################################

# Plot Boxplot + histograms for ALL variables grouped by OTUs
# Only to check for outliers

# loop for each taxon set
for(k in seq(Hlabel)){
	dir.create(paste("./04_boxplots/02_all_variables-",Hlabel[k],"-check",sep=""),showWarnings = FALSE)
	#loop for each variable
	for(i in 1:ncol(x_all_cln0)){
		#open pdf file
		pdf(paste("./04_boxplots/02_all_variables-",Hlabel[k],"-check/boxplot_",colnames(x_all_cln0)[i],".pdf",sep=""),family="ArialMT", width=length(OTUs[[k]]), height=12)
			#plot boxplot
			par(mar = c(7, 4, 4, 2) + 0.1)
			boxplot(x_all_cln0[,i] ~ x[,Hindex[k]],main=colnames(x_all_cln0)[i],outline=F,ylim=c(min(na.omit(x_all_cln0[,i])),max(na.omit(x_all_cln0[,i]))),boxwex = 0.8,xaxt = "n", xlab= Hlabel[k], ylab ="")
	# 		axis(1, labels = FALSE)
			xlabels <- levels(x[,Hindex[k]])
			text(1:length(levels(x[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE)
			#loop for each OTUs
			for(j in 1:length(levels(x[,Hindex[k]]))){
				if(length(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))==0){
					#plot N for variable
					text(j,mean(na.omit(x_all_cln0[,i])), labels=paste("N = ",length(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
				}
				else{
					#get histogram data
					b<-hist(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]],breaks=length(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]),plot=FALSE)
					#remove zeros from counts
					q<-b$counts
					q[q==0]<-NA
					#transform counts to boxplot coordinates
					m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
					if(m!=-Inf){
						y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
					}
					else{
						if(length(q)>1){
							y=q
							y[!is.na(y)]<-j+0.4
						}
						else{
							y=j+0.4
						}
					}
					#plot points from histogram
					if(length(b$mids)<2){
						points(y,unique(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="blue")
					}
					else{
						points(y,b$mids,col="blue")
					}
					#plot mean as a line
					segments(j-0.4,mean(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="red")
					#plot N for variable
					text(j,max(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					#get range amplitude
					rango<-max(na.omit(x_all_cln0[,i]))-min(na.omit(x_all_cln0[,i]))
					#plot axis for histogram
					segments(j-0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
					#plot max line for histogram
					segments(j-0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j-0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="grey")
					#plot mean line for histogram
					segments(j,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="grey")
					#plot min line for histogram
					segments(j+0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j+0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="grey")
					#plot labels for the histogram axis
					text(j,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
					text(j-0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
					text(j+0.4,min(na.omit(x_all_cln0[,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
				}	
			}
		dev.off()
	}
}

#######################################

##########################################################################################
# CLEANING BAD TERMINALS
##########################################################################################

#######################################

# adjust dataset by removing outliers
# OBS: consider removing outliers completelly if there are doubts about collected data or too many missing variables for them

# include in the vector "to_remove01" all IDs for outliers that you want to completelly remove from your analyses.
# EX:
# to_remove01<-c("1381","2","1337","820","1272","410")

to_remove01<-c()

x_cln<-x_cln0[!x_cln0$ID%in%to_remove01,]
x_all_cln<-x_all_cln0[!x_cln0$ID%in%to_remove01,]
x_multi_cln<-x_multi_cln0[!x_cln0$ID%in%to_remove01,]


# write list of previouly removed terminals
removed_terminals<-cbind(x_cln0$ID[which(x_cln0$ID%in%to_remove01)],x_cln0$Voucher[which(x_cln0$ID%in%to_remove01)])
write.table(removed_terminals,file=paste("./01_tables/Removed_Terminals.txt",sep=""),sep="\t",col.names=F)


#######################################

# removing outliers based on data values - removing a specific value only when it is out of 2*(interquartile range) for that particular variable

out_check<-list()

# loop for each taxon set
for(k in seq(Hlabel)){
	#loop for each variable
	for(i in 1:ncol(x_all_cln)){
		outlier_values <- boxplot(x_all_cln[,i] ~ x_cln[,Hindex[k]],range=2,plot=F)$stats
		outliers_ind<-vector()
		for (j in seq(OTUs[[k]])) {
			outliers_ind<-c(outliers_ind,x_cln$ID[x_cln[,Hindex[k]]==OTUs[[k]][j]][x_all_cln[x_cln[,Hindex[k]]==OTUs[[k]][j],i]<outlier_values[1,j]|x_all_cln[x_cln[,Hindex[k]]==OTUs[[k]][j],i]>outlier_values[5,j]])
		}
		out_check[[Hlabel[k]]][[colnames(x_all_cln)[i]]]<-unique(na.omit(outliers_ind))
		# names(out_check[k])[i]<-colnames(x_all_cln)[i]
	}
}

# write list of removed outliers
sink("./01_tables/list_of_removed_outliers.txt")
	print(out_check)
sink()

# cleaning dataset from outliers

x_cln_noOut<-list()
for(k in seq(Hlabel)){
	x_cln_noOut[[Hlabel[k]]]<-x_cln
	for(i in seq(out_check[[k]])){
		x_cln_noOut[[k]][which(x_cln_noOut[[k]]$ID%in%out_check[[k]][[i]]),grep(names(out_check[[k]][i]),colnames(x_cln_noOut[[k]]))]<-NA
	}
}

x_all_cln_noOut<-list()
for(k in seq(Hlabel)){
	x_all_cln_noOut[[Hlabel[k]]]<-x_all_cln
	for(i in seq(out_check[[k]])){
		x_all_cln_noOut[[k]][which(x_cln_noOut[[k]]$ID%in%out_check[[k]][[i]]),grep(names(out_check[[k]][i]),colnames(x_all_cln_noOut[[k]]))]<-NA
	}
}

x_multi_cln_noOut<-list()
for(k in seq(Hlabel)){
	x_multi_cln_noOut[[Hlabel[k]]]<-x_multi_cln
	for(i in seq(out_check[[k]])){
		x_multi_cln_noOut[[k]][which(x_cln_noOut[[k]]$ID%in%out_check[[k]][[i]]),grep(names(out_check[[k]][i]),colnames(x_multi_cln_noOut[[k]]))]<-NA
	}
}


#######################################

# removing zero-variance variables after correcting for outliers

zero_var<-list()
x_all_cln_noOut_noZ<-list()
for(k in seq(Hlabel)){
	zero_var[[Hlabel[k]]]<-which(apply(na.omit(x_all_cln_noOut[[k]]), 2, var) == 0)
	x_all_cln_noOut_noZ[[Hlabel[k]]]<-x_all_cln_noOut[[k]][-as.numeric(which(apply(na.omit(x_all_cln_noOut[[k]]), 2, var) == 0))]
}

# write list of zero-variance variables
sink("./01_tables/zero-variance_variables.txt")
	print(zero_var)
sink()

zero_var<-list()
x_multi_cln_noOut_noZ<-list()
for(k in seq(Hlabel)){
	zero_var[[Hlabel[k]]]<-which(apply(na.omit(x_multi_cln_noOut[[k]]), 2, var) == 0)
	x_multi_cln_noOut_noZ[[Hlabel[k]]]<-x_multi_cln_noOut[[k]][-as.numeric(which(apply(na.omit(x_multi_cln_noOut[[k]]), 2, var) == 0))]
}

#######################################

##########################################################################################
# PLOT QQPLOTS
##########################################################################################

#######################################

# Plot Quantiles-Quantiles normal plot for ALL variables
# Only to check distribution. 
# OBS: Discrete and categorical variables are not Normal distributed by definition (they can follow Poisson, binomial, etc). Only use this to check for remaining outliers and possible problems

dir.create("./03_QQplots/01_all_variables-no-groups",showWarnings = FALSE)

#qqnorm for ALL variables
for(k in seq(Hlabel)){
	dir.create(paste("./03_QQplots/01_all_variables-no-groups/",Hlabel[k],sep=""),showWarnings = FALSE)
	for (i in 1:length(x_all_cln_noOut_noZ[[k]][1,])){
	pdf(paste("./03_QQplots/01_all_variables-no-groups/",Hlabel[k],"/qqnorm_",colnames(x_all_cln_noOut_noZ[[k]])[i],".pdf",sep=""), family="ArialMT", width=9, height=12)
		qqnorm(x_all_cln_noOut_noZ[[k]][,i], main = colnames(x_all_cln_noOut_noZ[[k]][i]))
		qqline(x_all_cln_noOut_noZ[[k]][,i], col=2)
		if (length(x_all_cln_noOut_noZ[[k]][,i])<=50){
			text(-1,max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))-(diff(range(na.omit(x_all_cln_noOut_noZ[[k]][,i])))*0.1),paste("Shapiro-Wilk - p.value = ",round(shapiro.test(x_all_cln_noOut_noZ[[k]][,i])$p.value,4)))
		} else {
			text(-0.9,max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))-(diff(range(na.omit(x_all_cln_noOut_noZ[[k]][,i])))*0.1),paste("Lilliefors - p.value = ",round(lillie.test(na.omit(x_all_cln_noOut_noZ[[k]][,i]))$p.value,4)))
		}	
	dev.off()
	}
}

#######################################

##########################################################################################
# TEST FOR SEXUAL DIMORPHISM 
##########################################################################################

#######################################

# Statistical test for sexual dimorphism for ALL variables
# The t-test will be used if the variable is continuous or some ratios, normally distributed and the variances are equal, otherwise the Wilcoxon test will be applied.
# OBS: This test will not evaluate SD within OTUs to avoid small sample sizes, but if your sampling is comprehensive enough or your OTUs are phylogenetically distant consider changing the script.

#new matrix
dimorphism<-list()

for(k in seq(Hlabel)){
	dimorphism[[Hlabel[k]]]<-matrix(nrow=ncol(x_all_cln_noOut_noZ[[k]]),ncol=2)
	#row names
	rownames(dimorphism[[k]])<-colnames(x_all_cln_noOut_noZ[[k]])
	colnames(dimorphism[[k]])<-c("Test","p.value")
	#loop to fill table
	for (i in 1:ncol(x_all_cln_noOut_noZ[[k]])){
		if(length(grep("cont",colnames(x_all_cln_noOut_noZ[[k]])[i]))==1|length(grep("rat_",colnames(x_all_cln_noOut_noZ[[k]])[i]))==1){
			if (fligner.test(x_all_cln_noOut_noZ[[k]][,i]~x_cln$Sex,data=x_all_cln_noOut_noZ[[k]])$p.value>=0.05) {
				dimorphism[[k]][i,1]<-"T.test"
				dimorphism[[k]][i,2]<-round(t.test(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i],x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i])$p.value,4)
			} else {
				dimorphism[[k]][i,1]<-"Wilcoxon.test"
				dimorphism[[k]][i,2]<-round(wilcox.test(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i],x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i])$p.value,4)
			}
		} else {
			dimorphism[[k]][i,1]<-"Wilcoxon.test"
			dimorphism[[k]][i,2]<-round(wilcox.test(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i],x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i])$p.value,4)
		}
	}
	rownames(dimorphism[[k]])<-gsub("disc_","",rownames(dimorphism[[k]]))
	rownames(dimorphism[[k]])<-gsub("cont_","",rownames(dimorphism[[k]]))
	rownames(dimorphism[[k]])<-gsub("rat_","",rownames(dimorphism[[k]]))
	rownames(dimorphism[[k]])<-gsub("cat_","",rownames(dimorphism[[k]]))
	write.table(dimorphism[[k]],file=paste("./01_tables/dimorphism-no_groups-",Hlabel[k],".txt",sep=""),sep="\t")
}

# dimorphism<-read.table("./01_tables/dimorphism-no_groups.txt")

#######################################

##########################################################################################
# BOXPLOT + HISTOGRAM - FINAL
##########################################################################################

#######################################

# Boxplot + histogram for ALL variables - per OTUs
# Final

for(k in 1:length(Hlabel)){
	dir.create(paste("./04_boxplots/03_all_variables-",Hlabel[k],sep=""),showWarnings = FALSE)
	#loop for each variable
	for(i in 1:ncol(x_all_cln_noOut_noZ[[k]])){
		#open pdf file
		pdf(paste("./04_boxplots/03_all_variables-",Hlabel[k],"/boxplot_",colnames(x_all_cln_noOut_noZ[[k]])[i],".pdf",sep=""),family="ArialMT", width=9, height=12)
			#plot boxplot
			par(mar = c(7, 4, 4, 2) + 0.1)
			boxplot(x_all_cln_noOut_noZ[[k]][,i] ~ x[,Hindex[k]],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= "", ylab ="")
	# 		axis(1, labels = FALSE)
			xlabels <- levels(x[,Hindex[k]])
			text(1:length(levels(x[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1.1, labels = xlabels, xpd = TRUE)
			#loop for each OTUs
			for(j in 1:length(levels(x[,Hindex[k]]))){
				if(length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))==0){
					#plot N for variable
					text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
				}
				else{
					#transform counts to boxplot coordinates
					q<-as.numeric(table(x_all_cln_noOut_noZ[[k]][x[,Hindex[k]]==xlabels[j],i]))
					m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
					if(m!=-Inf){
						y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
					}else{
						if(length(q)>1){
							y=q
							y[!is.na(y)]<-j+0.4
						}else{
							y=j+0.4
						}
					}
					#plot points from histogram
					if(length(q)<2){
						points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x[,Hindex[k]]==xlabels[j],i])),col="blue",cex=0.8)
					}else{
						# points(y,b$mids,col="blue")
						points(y,as.numeric(attr(table(x_all_cln_noOut_noZ[[k]][x[,Hindex[k]]==xlabels[j],i]),"dimnames")[[1]]),col="blue",cex=0.8)
					}
					#plot mean as a line
					segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="red")
					#plot N for variable
					text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					#get range amplitude
					rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][,i]))
					#plot axis for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
					#plot max line for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="grey")
					#plot mean line for histogram
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="grey")
					#plot min line for histogram
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]])),col="grey")
					#plot labels for the histogram axis
					text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
					text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
					text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x[,Hindex[k]]==levels(x[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
				}	
			}
		dev.off()
	}
}

#######################################

# Boxplot + histogram for ALL variables - per Sex (if dimorphic)- per OTUs
# Need to run/load test for sexual dimorphism first

for(k in 1:length(Hlabel)){
	dir.create(paste("./04_boxplots/04_all_variables-",Hlabel[k],"-sex",sep=""),showWarnings = FALSE)
	#loop for each variable
	for(i in 1:ncol(x_all_cln_noOut_noZ[[k]])){
		if(as.numeric(dimorphism[[k]][i,2])>=0.05){
			#open pdf file
			pdf(paste("./04_boxplots/04_all_variables-",Hlabel[k],"-sex/boxplot_",colnames(x_all_cln_noOut_noZ[[k]])[i],".pdf",sep=""),family="ArialMT", width=3+length(OTUs[[k]]), height=14)
				#plot boxplot
				boxplot(x_all_cln_noOut_noZ[[k]][,i] ~ x_cln_noOut[[k]][,Hindex[k]],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= "", ylab ="")
				xlabels <- levels(x_cln_noOut[[k]][,Hindex[k]])
				text(1:length(levels(x_cln_noOut[[k]][,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1.1, labels = xlabels, xpd = TRUE)
				#loop for each OTUs
				for(j in 1:length(levels(x_cln_noOut[[k]][,Hindex[k]]))){
					if(length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))==0){
						#plot N for variable
						text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					}
					else{
						q<-as.numeric(table(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]][,Hindex[k]]==xlabels[j],i]))
						#transform counts to boxplot coordinates
						m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
						if(m!=-Inf){
							y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
						}else{
							if(length(q)>1){
								y=q
								y[!is.na(y)]<-j+0.4
							}else{
								y=j+0.4
							}
						}
						#plot points from histogram
						if(length(q)<2){
							points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]][,Hindex[k]]==xlabels[j],i])),col="blue",cex=0.8)
						}else{
							# points(y,b$mids,col="blue")
							points(y,as.numeric(attr(table(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]][,Hindex[k]]==xlabels[j],i]),"dimnames")[[1]]),col="blue",cex=0.8)
						}
						#plot mean as a line
						segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="red")
						#plot N for variable
						text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
						#get range amplitude
						rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][,i]))
						#plot axis for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
						#plot max line for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot mean line for histogram
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot min line for histogram
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot labels for the histogram axis
						text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
						text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
						text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln_noOut[[k]][,Hindex[k]]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
					}
				}
			dev.off()
		} else {
			#open pdf file
			pdf(paste("./04_boxplots/04_all_variables-",Hlabel[k],"-sex/boxplot_",colnames(x_all_cln_noOut_noZ[[k]])[i],".pdf",sep=""),family="ArialMT", width=10+length(OTUs[[k]]), height=14)
				layout(matrix(c(1,2),1,2))
				#plot boxplot for females
				boxplot(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i] ~ x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= "F", ylab ="")
				xlabels <- levels(x_cln_noOut[[k]][,Hindex[k]])
				text(1:length(levels(x_cln_noOut[[k]][,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1.1, labels = xlabels, xpd = TRUE)
				#loop for each OTUs
				for(j in 1:length(levels(x_cln_noOut[[k]][,Hindex[k]]))){
					if(length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))==0){
						#plot N for variable
						text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					}else{
						q<-as.numeric(table(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))
						# q[q==0]<-NA
						#transform counts to boxplot coordinates
						m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
						if(m!=-Inf){
							y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
						}else{
							if(length(q)>1){
								y=q
								y[!is.na(y)]<-j+0.4
							}else{
								y=j+0.4
							}
						}
						#plot points from histogram
						if(length(q)<2){
							points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="blue",cex=0.8)
						}else{
							# points(y,b$mids,col="blue")
							points(y,as.numeric(attr(table(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]),"dimnames")[[1]]),col="blue",cex=0.8)
						}
						#plot mean as a line
						segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="red")
						#plot N for variable
						text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
						#get range amplitude
						rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i]))
						#plot axis for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
						#plot max line for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot mean line for histogram
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot min line for histogram
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot labels for the histogram axis
						text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
						text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
						text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="F",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
					}
				}
				#plot boxplot for males
				boxplot(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i] ~ x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= "M", ylab ="")
				xlabels <- levels(x_cln_noOut[[k]][,Hindex[k]])
				text(1:length(levels(x_cln_noOut[[k]][,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1.1, labels = xlabels, xpd = TRUE)
				#loop for each OTUs
				for(j in 1:length(levels(x_cln_noOut[[k]][,Hindex[k]]))){
					if(length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))==0){
						#plot N for variable
						text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					}else{
						q<-as.numeric(table(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))
						#transform counts to boxplot coordinates
						m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
						if(m!=-Inf){
							y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
						}else{
							if(length(q)>1){
								y=q
								y[!is.na(y)]<-j+0.4
							}else{
								y=j+0.4
							}
						}
						#plot points from histogram
						if(length(q)<2){
							points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="blue",cex=0.8)
						}else{
							# points(y,b$mids,col="blue")
							points(y,as.numeric(attr(table(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]),"dimnames")[[1]]),col="blue",cex=0.8)
						}
						#plot mean as a line
						segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="red")
						#plot N for variable
						text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
						#get range amplitude
						rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i]))
						#plot axis for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
						#plot max line for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot mean line for histogram
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot min line for histogram
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]])),col="grey")
						#plot labels for the histogram axis
						text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
						text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
						text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln_noOut[[k]]$Sex=="M",i][x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"]==levels(x_cln_noOut[[k]][,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
					}
				}
			dev.off()
		}
	}
}

#######################################


##########################################################################################
# composed pdf for figure
##########################################################################################

#######################################


# Boxplot + histogram for VIP - per Sex (if dimorphic)- per OTUs
# Need to run/load test for sexual dimorphism first

vip<-c("disc_VEN","disc_DOA","disc_DOM")

for(k in 1:length(Hlabel)){
	dir.create(paste("./09_Combined/01_combined-",Hlabel[k],"-sex",sep=""),showWarnings = FALSE)
	vipInd<-which(colnames(x_all_cln_noOut_noZ[[k]])%in%vip)
	#open pdf file
	pdf(paste("./09_Combined/01_combined-",Hlabel[k],"-sex/boxplot.pdf",sep=""),family="ArialMT", width=8, height=11)
		layout(matrix(c(seq(length(vip)*2)),length(vip),2,byrow=T))
		par(mar=c(4,2,2,2))
		#loop for each variable
		for(i in seq(vipInd)){
			#plot boxplot for females
			boxplot(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]] ~ x_cln[,Hindex[k]][x_cln$Sex=="F"],main=paste(gsub("disc_","",colnames(x_all_cln_noOut_noZ[[k]])[vipInd[i]]),"- Female",sep=" "), outline=F, ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]]))-(diff(range(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]])))*0.02),max(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]]))+(diff(range(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]])))*0.02)),boxwex = 0.8,xaxt = "n", xlab= "", ylab ="")
			xlabels <- levels(x_cln[,Hindex[k]])
			text(1:length(levels(x_cln[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE, cex=0.8)
			#loop for each OTUs
			for(j in 1:length(levels(x_cln[,Hindex[k]]))){
				if(length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))==0){
					#plot N for variable
					text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
				}else{
					#get histogram data
					# b<-hist(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]],breaks=length(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]),plot=FALSE)
					
					#remove zeros from counts
					# q<-b$counts
					q<-as.numeric(table(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))
					# q[q==0]<-NA
					#transform counts to boxplot coordinates
					m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
					if(m!=-Inf){
						y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
					}else{
						if(length(q)>1){
							y=q
							y[!is.na(y)]<-j+0.4
						}else{
							y=j+0.4
						}
					}
					#plot points from histogram
					if(length(q)<2){
						points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="blue",cex=0.8)
					}else{
						# points(y,b$mids,col="blue")
						points(y,as.numeric(attr(table(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]),"dimnames")[[1]]),col="blue",cex=0.8)
					}
					#plot mean as a line
					segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="red")
					#plot N for variable
					text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					#get range amplitude
					rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]]))-min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]]))
					#plot axis for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
					#plot max line for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot mean line for histogram
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot min line for histogram
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot labels for the histogram axis
					text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
					text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
					text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
				}
			}
			#plot boxplot for males
			boxplot(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]] ~ x_cln[,Hindex[k]][x_cln$Sex=="M"],main=paste(gsub("disc_","",colnames(x_all_cln_noOut_noZ[[k]])[vipInd[i]]),"- Male",sep=" "),outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]]))-(diff(range(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]])))*0.02),max(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]]))+(diff(range(na.omit(x_all_cln_noOut_noZ[[k]][,vipInd[i]])))*0.02)),boxwex = 0.8,xaxt = "n", xlab= "", ylab ="")
			xlabels <- levels(x_cln[,Hindex[k]])
			text(1:length(levels(x_cln[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE, cex=0.8)
			#loop for each OTUs
			for(j in 1:length(levels(x_cln[,Hindex[k]]))){
				if(length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))==0){
					#plot N for variable
					text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
				}else{
					#get histogram data
					# b<-hist(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]],breaks=length(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]),plot=FALSE)
					#remove zeros from counts
					# q<-b$counts
					# q[q==0]<-NA
					q<-as.numeric(table(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))
					#transform counts to boxplot coordinates
					m<-((j-0.4)-(j+0.4))/(max(na.omit(q))-min(na.omit(q)))
					if(m!=-Inf){
						y<-(m*q)+((j+0.4)-(m*min(na.omit(q))))
					}else{
						if(length(q)>1){
							y=q
							y[!is.na(y)]<-j+0.4
						}else{
							y=j+0.4
						}
					}
					#plot points from histogram
					if(length(q)<2){
						points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="blue",cex=0.8)
					}else{
						# points(y,b$mids,col="blue")
						points(y,as.numeric(attr(table(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]),"dimnames")[[1]]),col="blue",cex=0.8)
					}
					#plot mean as a line
					segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="red")
					#plot N for variable
					text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					#get range amplitude
					rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]]))-min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]]))
					#plot axis for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
					#plot max line for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot mean line for histogram
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot min line for histogram
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot labels for the histogram axis
					text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
					text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
					text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",vipInd[i]][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
				}
			}
		}
	dev.off()
}

#######################################

##########################################################################################
# SUMMARIES
##########################################################################################

#######################################

# make summary table
for(k in 1:length(Hlabel)){
	dir.create(paste("./01_tables/01_summaries-",Hlabel[k],sep=""),showWarnings = FALSE)
	OTUs<-unique(x[,Hindex[k]])	
	data_summary_up<-character()
	data_summary_down<-character()
	for (i in 1:length(OTUs)){
		all_up<-character()
		all_down<-character()
		for(j in 1:length(x_multi_cln[1,])){
			aA<-round(mean(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i],j])),2)
			bA<-round(sd(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i],j])),2)
			cA<-round(min(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i],j])),2)
			dA<-round(max(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i],j])),2)
			eA<-length(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i],j]))
			mA<-round(median(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i],j])),2)
			fA<-(paste("range=",cA,"-",dA,", n=",eA,sep=""))
			gA<-c("x=",aA,", s=",bA,", m=",mA)
			hA<-paste(gA,collapse ="")
			all_up<-c(all_up,fA)
			all_down<-c(all_down,hA)
		}
		data_summary_up<-cbind(data_summary_up,all_up)
		data_summary_down<-cbind(data_summary_down,all_down)
	}
	rownames(data_summary_up)<-colnames(x_multi_cln)
	rownames(data_summary_up)<-gsub("disc_","",rownames(data_summary_up))
	rownames(data_summary_up)<-gsub("cont_","",rownames(data_summary_up))
	rownames(data_summary_up)<-gsub("rat_","",rownames(data_summary_up))
	OTUs_labels<-as.character(OTUs)
	l <- list(a=data_summary_up,b=data_summary_down)
	data_summary <- do.call(rbind, l)[order(sequence(sapply(l, nrow))), ]
	data_summary<-rbind(OTUs=OTUs_labels,data_summary)
	write.table(data_summary,file=paste("./01_tables/01_summaries-",Hlabel[k],"/summary-ALL-variables.txt",sep=""),sep="\t",col.names=F)
}

#########################################

# make summary table for females and males
for(k in 1:length(Hlabel)){
	dir.create(paste("./01_tables/01_summaries-",Hlabel[k],sep=""),showWarnings = FALSE)
	OTUs<-unique(x[,Hindex[k]])	
	data_summary_up<-character()
	data_summary_down<-character()
	for (i in 1:length(OTUs)){
		females_up<-character()
		females_down<-character()
		males_up<-character()
		males_down<-character()
		for(j in 1:length(x_multi_cln[1,])){
			aF<-round(mean(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="F",j])),2)
			bF<-round(sd(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="F",j])),2)
			cF<-round(min(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="F",j])),2)
			dF<-round(max(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="F",j])),2)
			eF<-length(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="F",j]))
			mF<-round(median(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="F",j])),2)
			fF<-(paste("range=",cF,"-",dF,", n=",eF,sep=""))
			gF<-c("x=",aF,", s=",bF,", m=",mF)
			hF<-paste(gF,collapse ="")
			females_up<-c(females_up,fF)
			females_down<-c(females_down,hF)
			aM<-round(mean(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="M",j])),2)
			bM<-round(sd(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="M",j])),2)
			cM<-round(min(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="M",j])),2)
			dM<-round(max(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="M",j])),2)
			eM<-length(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="M",j]))
			mM<-round(median(na.omit(x_multi_cln[x_cln[,Hindex[k]]==OTUs[i]&x_cln$Sex=="M",j])),2)
			fM<-(paste("range=",cM,"-",dM,", n=",eM,sep=""))
			gM<-c("x=",aM,", s=",bM,", m=",mM)
			hM<-paste(gM,collapse ="")
			males_up<-c(males_up,fM)
			males_down<-c(males_down,hM)			
		}
		data_summary_up<-cbind(data_summary_up,females_up,males_up)
		data_summary_down<-cbind(data_summary_down,females_down,males_down)
	}
	rownames(data_summary_up)<-colnames(x_multi_cln)
	rownames(data_summary_up)<-gsub("disc_","",rownames(data_summary_up))
	rownames(data_summary_up)<-gsub("cont_","",rownames(data_summary_up))
	rownames(data_summary_up)<-gsub("rat_","",rownames(data_summary_up))
	rownames(data_summary_down)<-rep("",length(x_multi_cln[1,]))
	sex<-rep(c("Female","Male"),length(OTUs))
	OTUs_labels<-c(rbind(as.character(OTUs),rep("",length(OTUs))))
	l <- list(a=data_summary_up,b=data_summary_down)
	data_summary <- do.call(rbind, l)[order(sequence(sapply(l, nrow))), ]
	data_summary<-rbind(OTUs=OTUs_labels,Sex=sex,data_summary)
	write.table(data_summary,file=paste("./01_tables/01_summaries-",Hlabel[k],"/summary-ALL-variables-F&M.txt",sep=""),sep="\t",col.names=F)
}

#########################################

##########################################################################################
##########################################################################################
# END
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################
