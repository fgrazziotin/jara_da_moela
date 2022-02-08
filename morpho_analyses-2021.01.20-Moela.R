##########################################################################################
# R SCRIPT FOR MORPHOMETRIC ANALYSES (having snakes in mind...)
# version 2021.02.02
# Designed to analyze data for Barbo et al., 2022 - Systematics and Biodiversity
##########################################################################################
# Felipe G. Grazziotin - fgrazziotin@gmail.com
##########################################################################################

# "All organisms contain an array of biological information, some of which is useful in a comparative or taxonomic sense, some of which is useless, and some of which may even be misleading."
# Judith E. Winston, 1999 (Describing Species, p. 53)

# "... Molecular Biology is beginning to be seen as a restricted science, narrowing our vision... the need for taxonomists to draw attention to the enormous diversity and variation of this earth's biota becomes more and more pressing."
# Tod. F. Stuessy, 1990 (Plant Taxonomy, p. xvii)

# "Is it not extraordinary that young taxonomists are trained like performing monkeys, almost wholly by imitation, and that in only the rarest cases are they given any instruction in taxonomic theory?
# A. J. Cain, 1959 (quoted in Simpson, 1961 - Principles of Animal Taxonomy, p. vii)

# Make sure you have identified all of the possible confounding variables in your study... always get dirty hands digging deep trying to identify them...

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
library(MASS)

# install.packages("corrplot")
# needed for corrplot
library(corrplot)

# install.packages("coin")
# needed for independence_test
library(coin)

# install.packages("wPerm")
# needed for perm.ind.test
library(wPerm)

# install.packages("caret")
# needed for PLSDA and RF
library(caret)

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
library(rfPermute)

# install.packages("psych")
# needed for some summary statistics
# library(psych)

# install.packages("vegan")
# needed for Permanova
library(vegan)

# install.packages("RVAideMemoire")
# needed for pairwise Post-hoc Permanova and PLSDA.VIP
library(RVAideMemoire)

# install.packages("doSNOW")
# needed to run CV in parallel
library(doSNOW)

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

x<-x_backup

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

# set hypotheses of OTU clustering
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
dir.create("./05_correlograms",showWarnings = FALSE)
dir.create("./06_PCAs",showWarnings = FALSE)
dir.create("./07_PLSDAs",showWarnings = FALSE)
dir.create("./08_RF",showWarnings = FALSE)
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
# The t-test will be used if the variable is continuous or some rate, normally distributed and the variances are equal, otherwise the Wilcoxon test will be applied.
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
		pdf(paste("./04_boxplots/03_all_variables-",Hlabel[k],"/boxplot_",colnames(x_all_cln_noOut_noZ[[k]])[i],".pdf",sep=""),family="ArialMT", width=3+length(OTUs[[k]]), height=12)
			#plot boxplot
			par(mar = c(7, 4, 4, 2) + 0.1)
			boxplot(x_all_cln_noOut_noZ[[k]][,i] ~ x_cln[,Hindex[k]],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= Hlabel[k], ylab ="")
	# 		axis(1, labels = FALSE)
			xlabels <- levels(x_cln[,Hindex[k]])
			text(1:length(levels(x_cln[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE)
			#loop for each OTUs
			for(j in 1:length(levels(x_cln[,Hindex[k]]))){
				if(length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))==0){
					#plot N for variable
					text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
				}
				else{
					#get histogram data
					b<-hist(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]],breaks=length(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]),plot=FALSE)
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
						points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="blue")
					}
					else{
						points(y,b$mids,col="blue")
					}
					#plot mean as a line
					segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="red")
					#plot N for variable
					text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					#get range amplitude
					rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][,i]))
					#plot axis for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
					#plot max line for histogram
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot mean line for histogram
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot min line for histogram
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
					segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
					#plot labels for the histogram axis
					text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
					text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
					text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
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
				boxplot(x_all_cln_noOut_noZ[[k]][,i] ~ x_cln[,Hindex[k]],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= Hlabel[k], ylab ="")
				xlabels <- levels(x_cln[,Hindex[k]])
				text(1:length(levels(x_cln[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE)
				#loop for each OTUs
				for(j in 1:length(levels(x_cln[,Hindex[k]]))){
					if(length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))==0){
						#plot N for variable
						text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					}
					else{
						#get histogram data
						b<-hist(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]],breaks=length(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]),plot=FALSE)
						#remove zeros from counts
						q<-b$counts
						q[q==0]<-NA
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
						if(length(b$mids)<2){
							points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="blue")
						}else{
							points(y,b$mids,col="blue")
						}
						#plot mean as a line
						segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="red")
						#plot N for variable
						text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
						#get range amplitude
						rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][,i]))
						#plot axis for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
						#plot max line for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot mean line for histogram
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot min line for histogram
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot labels for the histogram axis
						text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
						text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
						text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][,i][x_cln[,Hindex[k]]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
					}
				}
			dev.off()
		} else {
			#open pdf file
			pdf(paste("./04_boxplots/04_all_variables-",Hlabel[k],"-sex/boxplot_",colnames(x_all_cln_noOut_noZ[[k]])[i],".pdf",sep=""),family="ArialMT", width=10+length(OTUs[[k]]), height=14)
				layout(matrix(c(1,2),1,2))
				#plot boxplot for females
				boxplot(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i] ~ x_cln[,Hindex[k]][x_cln$Sex=="F"],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= paste(Hlabel[k]," Females",sep=""), ylab ="")
				xlabels <- levels(x_cln[,Hindex[k]])
				text(1:length(levels(x_cln[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE)
				#loop for each OTUs
				for(j in 1:length(levels(x_cln[,Hindex[k]]))){
					if(length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))==0){
						#plot N for variable
						text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					}else{
						#get histogram data
						b<-hist(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]],breaks=length(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]),plot=FALSE)
						#remove zeros from counts
						q<-b$counts
						q[q==0]<-NA
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
						if(length(b$mids)<2){
							points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="blue")
						}else{
							points(y,b$mids,col="blue")
						}
						#plot mean as a line
						segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="red")
						#plot N for variable
						text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
						#get range amplitude
						rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i]))
						#plot axis for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
						#plot max line for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot mean line for histogram
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot min line for histogram
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot labels for the histogram axis
						text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
						text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
						text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="F",i][x_cln[,Hindex[k]][x_cln$Sex=="F"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
					}
				}
				#plot boxplot for males
				boxplot(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i] ~ x_cln[,Hindex[k]][x_cln$Sex=="M"],main=colnames(x_all_cln_noOut_noZ[[k]])[i],outline=F,ylim=c(min(na.omit(x_all_cln_noOut_noZ[[k]][,i])),max(na.omit(x_all_cln_noOut_noZ[[k]][,i]))),boxwex = 0.8,xaxt = "n", xlab= paste(Hlabel[k]," Males",sep=""), ylab ="")
				xlabels <- levels(x_cln[,Hindex[k]])
				text(1:length(levels(x_cln[,Hindex[k]])), par("usr")[3], offset=1, srt = 45, adj = 1, labels = xlabels, xpd = TRUE)
				#loop for each OTUs
				for(j in 1:length(levels(x_cln[,Hindex[k]]))){
					if(length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))==0){
						#plot N for variable
						text(j,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
					}else{
						#get histogram data
						b<-hist(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]],breaks=length(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]),plot=FALSE)
						#remove zeros from counts
						q<-b$counts
						q[q==0]<-NA
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
						if(length(b$mids)<2){
							points(y,unique(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="blue")
						}else{
							points(y,b$mids,col="blue")
						}
						#plot mean as a line
						segments(j-0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,mean(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="red")
						#plot N for variable
						text(j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])), labels=paste("N = ",length(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),sep=""),cex=0.8,pos=3)
						#get range amplitude
						rango<-max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i]))-min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i]))
						#plot axis for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.01*rango),col="grey",adj=0)
						#plot max line for histogram
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j-0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot mean line for histogram
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot min line for histogram
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),j+0.4,max(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),lty=3,col="grey")
						segments(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.02*rango),j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]])),col="grey")
						#plot labels for the histogram axis
						text(j,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.025*rango),labels="Count",cex=0.8)
						text(j-0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=max(na.omit(q)),cex=0.8)
						text(j+0.4,min(na.omit(x_all_cln_noOut_noZ[[k]][x_cln$Sex=="M",i][x_cln[,Hindex[k]][x_cln$Sex=="M"]==levels(x_cln[,Hindex[k]])[j]]))-(0.03*rango),labels=min(na.omit(q)),cex=0.8)
					}
				}
			dev.off()
		}
	}
}

#######################################

##########################################################################################
# CLEANING DATA FOR MULTIVARIATE
##########################################################################################
# dataset cannot have missing entries for most of the multivariate analyses.

#######################################

# prepare dataset for multivariate for ALL variables - But NO OUTLIERS 
# fill missing entries with overal median
x_all_cln_noOut_noZ_cln<-x_all_cln_noOut_noZ

for (k in seq(Hlabel)) {
	for (i in 1:length(x_all_cln_noOut_noZ_cln[[k]][1,])) {
		m <- median(na.omit(x_all_cln_noOut_noZ_cln[[k]][,i]))
		x_all_cln_noOut_noZ_cln[[k]][,i][is.na(x_all_cln_noOut_noZ_cln[[k]][,i])] <- m
	}
}

#######################################

# prepare dataset for multivariate for RATIOS, transformed CAT and DISC variables
# fill missing entries with overal median
x_multi_cln_noOut_noZ_cln<-x_multi_cln_noOut_noZ

for (k in seq(Hlabel)) {
	for (i in 1:length(x_multi_cln_noOut_noZ_cln[[k]][1,])) {
		m <- median(na.omit(x_multi_cln_noOut_noZ_cln[[k]][,i]))
		x_multi_cln_noOut_noZ_cln[[k]][,i][is.na(x_multi_cln_noOut_noZ_cln[[k]][,i])] <- m
	}
}

#######################################

# prepare dataset for multivariate for RATIOS, transformed CAT and DISC variables
# fill missing entries by using proximity results from Random Forest

# x_multi_cln_noOut_noZ_Imp <- list()

# for (k in seq(Hlabel)) {
# 	x_multi_cln_noOut_noZ_Imp[[Hlabel[[k]]]] <- rfImpute(x=x_multi_cln_noOut_noZ[[Hlabel[[k]]]], y = x_cln_noOut[[k]][,Hindex[k]], iter=6)[,-1]
# }

#######################################

#########################################################################################
# Collinearity tests
#########################################################################################

#######################################

# Remove variables with collinearity
# OBS: Size measurements are frequently correlated and affect linear models

#create directory for correlograms
dir.create("./05_correlograms",showWarnings = FALSE)

# correlogram for ALL variables
corr_matrix<-list()
for (k in seq(Hlabel)) {
	corr_matrix[[Hlabel[k]]] <- cor(na.omit(x_all_cln_noOut_noZ_cln[[k]]))

	pdf(paste("./05_correlograms/corrplot-",Hlabel[k],"-ALL.pdf",sep=""), family="ArialMT", width=20, height=20)
		corrplot.mixed(corr_matrix[[Hlabel[k]]],tl.pos="lt")
	dev.off()
}

#######################################

# correlogram for ALL MULTI (discrete, categorical and RATIOS)
corr_matrix_multi<-list()
for (k in seq(Hlabel)) {
	corr_matrix_multi[[Hlabel[k]]] <- cor(na.omit(x_multi_cln_noOut_noZ_cln[[k]]))
	pdf(paste("./05_correlograms/corrplot-",Hlabel[k],"-MULTI.pdf",sep=""), family="ArialMT", width=20, height=20)
		corrplot.mixed(corr_matrix_multi[[Hlabel[k]]],tl.pos="lt")
	dev.off()
}

#######################################

##########################################################################################
# CLEANING Collinear VARIABLES
##########################################################################################

#######################################

# adjust dataset by removing variables presenting collinearity (if needed)

# manually
# to_remove01<-c("rat_RW2RW1",)

# automatically remove collinear variables

list_corVar<-list()
for (k in seq(Hlabel)) {
	temp<-which(corr_matrix_multi[[Hlabel[k]]]>=0.7, arr.ind=TRUE)
	corVar<-rownames(corr_matrix_multi[[Hlabel[k]]])[temp[!temp[,1]==temp[,2],][,1]]
	list_corVar[[Hlabel[k]]]<-cbind(corVar[seq(1,length(corVar),by=2)],corVar[seq(2,length(corVar),by=2)])
}

# it will randomly select which collinear variable to remove
x_all_cln_noOut_noZ_cln_noCor<-list()
x_multi_cln_noOut_noZ_cln_noCor<-list()
for(k in seq(Hlabel)){
	set.seed(1001)
	setcol<-sample(c(1,2),1)
	to_remove02<-list_corVar[[Hlabel[k]]][,setcol]
	x_all_cln_noOut_noZ_cln_noCor[[Hlabel[k]]]<-x_all_cln_noOut_noZ_cln[[Hlabel[k]]][,-which(colnames(x_all_cln_noOut_noZ_cln[[Hlabel[k]]])%in%c(to_remove02))]
	x_multi_cln_noOut_noZ_cln_noCor[[Hlabel[k]]]<-x_multi_cln_noOut_noZ_cln[[Hlabel[k]]][,-which(colnames(x_multi_cln_noOut_noZ_cln[[Hlabel[k]]])%in%c(to_remove02))]
}

#######################################

##########################################################################################
# SET COLOR AND SHAPES FOR PLOTS - MULTIVARIATE ANALYSES
##########################################################################################

#######################################
# create fixed color pallete for multivariated plots

# manually define colors

# get OTUs from all hypotheses
all_OTUs<-vector()
for(k in 1:length(Hindex)){
	all_OTUs<-c(all_OTUs,levels(x[,Hindex[k]]))
}

all_OTUs<-unique(all_OTUs)

# all_OTUs
# check them and order them as you wish

all_OTUs<-c("Mainland_NC","Mainland_NSC01","Queimada_Grande","Mainland_NSC02","Franceses","Mainland_NSC03","Sao_Sebastiao","Moela","Vitoria","Alcatrazes","Mainland_SC")

# define the colors
# myColors <- c("magenta4","darkgreen","PaleGreen4","NavyBlue","RoyalBlue","darkred","orange4","red","salmon3","OrangeRed1","goldenrod1")
myColors <- c("magenta4","darkgreen","mediumseagreen","NavyBlue","RoyalBlue","darkred","orange4","red","salmon3","darkorange1","goldenrod1")


names(myColors) <- all_OTUs

# set color and fill for all graphs
colScale <- scale_colour_manual(name = "OTUS",values = myColors)
fillScale <- scale_fill_manual(name = "OTUS",values = myColors)

#######################################

# manually define shapes 
# OBS: Here I'm using stars to represent the Islands and circles for Mainland. The UNICODE symbols can be messy with ggplot2, so I'm using inverted triangles to do the trick. If it's not you case you can skip this block or change it as you prefer.
# Island = stars (using two inverted filled triangles)
# Mainland = filled circle

Shape_numbers1<-list()
Shape_numbers2<-list()
for(k in 1:length(Hindex)){
	Shape_numbers1[[Hlabel[k]]]<-x_cln_noOut[[k]]$Geography
	Shape_numbers2[[Hlabel[k]]]<-x_cln_noOut[[k]]$Geography
	Shape_numbers2[[Hlabel[k]]][Shape_numbers2[[Hlabel[k]]]=="Island"]<-"Island2"
	Shape_numbers2[[Hlabel[k]]][Shape_numbers2[[Hlabel[k]]]=="Mainland"]<-NA
}

myShapes<-c('Island'=24, 'Mainland'=16, 'Island2'=25)

#######################################

##########################################################################################
# PCA
##########################################################################################

#######################################

# plot PCA for RATIOS, CATEGORICAL and DISC variables - polygon instead of ellipse 
for(k in seq(Hlabel)){
	dir.create(paste("./06_PCAs/01_PCA_multi_variables-",Hlabel[k],sep=""),showWarnings = FALSE)

	# plot pca for females - ALL variables
	pca_multi_F <- prcomp(x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln$Sex=="F",], scale.=T)
	dfscores_F <- as.data.frame(pca_multi_F$x)
	prop_pca_F = pca_multi_F$sdev^2/sum(pca_multi_F$sdev^2)
	forChull<-cbind(dfscores_F[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="F"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=`PC1`,y=`PC2`))
	g <- g + geom_polygon(data=afterChull, aes(x=`PC1`,y=`PC2`, fill = `OTUS` ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = paste("PC1 (", round(prop_pca_F[1]*100,1), "%)", sep=""), y = paste("PC2 (", round(prop_pca_F[2]*100,1), "%)", sep=""))
	pdf(paste("./06_PCAs/01_PCA_multi_variables-",Hlabel[k],"/PCA-plot-F.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()
	
	# plot pca for males - ALL variables
	pca_multi_M <- prcomp(x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln$Sex=="M",], scale.=T)
	dfscores_M <- as.data.frame(pca_multi_M$x)
	prop_pca_M = pca_multi_M$sdev^2/sum(pca_multi_M$sdev^2)
	forChull<-cbind(dfscores_M[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="M"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=`PC1`,y=`PC2`))
	g <- g + geom_polygon(data=afterChull, aes(x=`PC1`,y=`PC2`, fill = `OTUS` ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = paste("PC1 (", round(prop_pca_M[1]*100,1), "%)", sep=""), y = paste("PC2 (", round(prop_pca_M[2]*100,1), "%)", sep=""))
	pdf(paste("./06_PCAs/01_PCA_multi_variables-",Hlabel[k],"/PCA-plot-M.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()
	
	# plot pca - ALL variables
	pca_multi <- prcomp(x_multi_cln_noOut_noZ_cln_noCor[[k]], scale.=T)
	dfscores <- as.data.frame(pca_multi$x)
	prop_pca = pca_multi$sdev^2/sum(pca_multi$sdev^2)
	forChull<-cbind(dfscores[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]],Labels=x_cln_noOut[[k]]$ID,Shapes1=Shape_numbers1[[Hlabel[k]]],Shapes2=Shape_numbers2[[Hlabel[k]]])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=`PC1`,y=`PC2`))
	g <- g + geom_polygon(data=afterChull, aes(x=`PC1`,y=`PC2`, fill = `OTUS` ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = paste("PC1 (", round(prop_pca[1]*100,1), "%)", sep=""), y = paste("PC2 (", round(prop_pca[2]*100,1), "%)", sep=""))
	pdf(paste("./06_PCAs/01_PCA_multi_variables-",Hlabel[k],"/PCA-plot.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()	
}

#######################################

##########################################################################################
# Create training and testing datasets, as well multiple constant FOLDS for Machine Learning
##########################################################################################

# because our main goal is only to visually explore the morphological differences among OTUs, instead of strictly classify our samples, we will train our models twice: 
# 1) using the complete dataset - to get the best estimation of latent components (PSL) and proximity scores (RF); 
# 2) using 60% of our samples as training and 40% as testing - to briefly evaluate how our models are behaving (avoiding using the same samples for training and testing).
# OBS01: although both models will be different, it seems fair to suppose that the model build with the complete dataset will provide better classifications. Based on the same idea, by predicting the testing with a model built on only 60% of our dataset, we are proceeding in a more conservative way when reporting the parameters of model confidence.
# OBS02: this approach will generate two CMs, one based on the predictions provided by the CV and OTB methods (training with complete dataset) and another based on the predictions for the testing dataset. 

#######################################

# create folds and repeats for each OTU configuration for cross-validation - COMPLETE DATASET

myfolds_F <- list()
myfolds_M <- list()
myfolds <- list()
set.seed(1001)
tmp <- .Random.seed
save(tmp,file = "./01_tables/Fold_seeds.RData")
for(k in seq(Hlabel)){
	myfolds_F[[Hlabel[k]]] <- createMultiFolds(x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"], k = 2, times = 2)
	myfolds_M[[Hlabel[k]]] <- createMultiFolds(x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"], k = 2, times = 2)
	myfolds[[Hlabel[k]]] <- createMultiFolds(x_cln_noOut[[k]][,Hindex[k]], k = 2, times = 2)
}

#######################################

# create training and testing datasets keeping imbalance among OTUs for RF

training <- list()
testing <- list()
training_OTUs <- list()
testing_OTUs <- list()
training_F <- list()
testing_F <- list()
training_OTUs_F <- list()
testing_OTUs_F <- list()
training_M <- list()
testing_M <- list()
training_OTUs_M <- list()
testing_OTUs_M <- list()

for(k in seq(Hlabel)){
	indexes <- as.vector(createDataPartition(x_cln_noOut[[k]][,Hindex[k]], p = 0.6))
	training[[Hlabel[k]]] <- x_multi_cln_noOut_noZ_cln_noCor[[k]][indexes$Resample1,]
	testing[[Hlabel[k]]] <- x_multi_cln_noOut_noZ_cln_noCor[[k]][-indexes$Resample1,]
	training_OTUs[[Hlabel[k]]] <- x_cln_noOut[[k]][,Hindex[k]][indexes$Resample1]
	testing_OTUs[[Hlabel[k]]] <- x_cln_noOut[[k]][,Hindex[k]][-indexes$Resample1]
	indexes_F <- as.vector(createDataPartition(x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"], p = 0.6))
	training_F[[Hlabel[k]]] <- x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",][indexes_F$Resample1,]
	testing_F[[Hlabel[k]]] <- x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",][-indexes_F$Resample1,]
	training_OTUs_F[[Hlabel[k]]] <- x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"][indexes_F$Resample1]
	testing_OTUs_F[[Hlabel[k]]] <- x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"][-indexes_F$Resample1]
	indexes_M <- as.vector(createDataPartition(x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"], p = 0.6))
	training_M[[Hlabel[k]]] <- x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",][indexes_M$Resample1,]
	testing_M[[Hlabel[k]]] <- x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",][-indexes_M$Resample1,]
	training_OTUs_M[[Hlabel[k]]] <- x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"][indexes_M$Resample1]
	testing_OTUs_M[[Hlabel[k]]] <- x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"][-indexes_M$Resample1]
}

	# uncomment the following lines to check the proportion of samples among OTUs
	# round(100*table(x_cln_noOut[[k]][,Hindex[k]])/length(x_cln_noOut[[k]][,Hindex[k]]))
	# round(100*table(x_cln_noOut[[k]][,Hindex[k]][indexes$Resample1])/length(x_cln_noOut[[k]][,Hindex[k]][indexes$Resample1]))
	# table(x_cln_noOut[[k]][,Hindex[k]])
	# table(x_cln_noOut[[k]][,Hindex[k]][indexes$Resample1])

	# uncomment the following lines to check the proportion of samples among OTUs
	# round(100*table(testing_OTUs)/length(testing_OTUs))
	# round(100*table(training_OTUs)/length(training_OTUs))


#######################################

# create folds and repeats for each OTU configuration for cross-validation - TRAINING AND TESTING

myfolds_train_F <- list()
myfolds_train_M <- list()
myfolds_train <- list()
set.seed(1001)
tmp <- .Random.seed
save(tmp,file = "./01_tables/Fold_seeds.RData")
for(k in seq(Hlabel)){
	myfolds_train_F[[Hlabel[k]]] <- createMultiFolds(training_OTUs_F[[Hlabel[k]]], k = 10, times = 3)
	myfolds_train_M[[Hlabel[k]]] <- createMultiFolds(training_OTUs_M[[Hlabel[k]]], k = 10, times = 3)
	myfolds_train[[Hlabel[k]]] <- createMultiFolds(training_OTUs[[Hlabel[k]]], k = 10, times = 3)
}

##########################################################################################
# PLS-DA 
##########################################################################################

# if want to use imputed by RF proximities
# x_multi_cln_noOut_noZ_cln_noCor<-x_multi_cln_noOut_noZ_Imp

#######################################

# create lists for final tables
cmPls <- list()
cmPls_M <- list()
cmPls_F <- list()
cmPls_train <- list()
cmPls_train_M <- list()
cmPls_train_F <- list()
varImpPLS <- list()
varImpPLS_M <- list()
varImpPLS_F <- list()
plsda_multi <- list()
plsda_multi_M <- list()
plsda_multi_F <- list()
plsda_multi_train <- list()
plsda_multi_train_M <- list()
plsda_multi_train_F <- list()
plsPred <- list()
plsPred_M <- list()
plsPred_F <- list()
plsPred_train <- list()
plsPred_train_M <- list()
plsPred_train_F <- list()

# start parallel processes
cl <- makeCluster(4,type = "SOCK")
registerDoSNOW(cl)

# plot PLS-DA for MULTI variables - NO OUTLIERS - polygon instead of ellipse 
for(k in seq(Hlabel)){
	dir.create(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],sep=""),showWarnings = FALSE)

	# whole dataset
	control <- trainControl("repeatedcv", index = myfolds[[k]], selectionFunction = "oneSE")
	control_train<- trainControl("repeatedcv", index = myfolds_train[[k]], selectionFunction = "oneSE")
	plsda_multi[[Hlabel[k]]]<-train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]],y=x_cln_noOut[[k]][,Hindex[k]],method="pls",tuneLength=30,trControl=control,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)
	plsda_multi_train[[Hlabel[k]]] <- train(x=training[[k]],y=training_OTUs[[k]],method="pls",tuneLength=30,trControl=control_train,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)
	plsPred_train[[Hlabel[k]]] <- predict(plsda_multi_train[[Hlabel[k]]], newdata = testing[[k]])
	cmPls_train[[Hlabel[k]]] <- caret::confusionMatrix(data=plsPred_train[[Hlabel[k]]], testing_OTUs[[k]])
	plsPred[[Hlabel[k]]] <- predict(plsda_multi[[Hlabel[k]]], newdata = x_multi_cln_noOut_noZ_cln_noCor[[k]])
	cmPls[[Hlabel[k]]] <- caret::confusionMatrix(data=plsPred[[Hlabel[k]]], x_cln_noOut[[k]][,Hindex[k]])
	p<-ggplot(plsda_multi[[Hlabel[k]]])
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-CV-accuracy_by_Comp.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(p)
	dev.off()
	varImpPLS[[Hlabel[k]]]<-varImp(plsda_multi[[Hlabel[k]]])
	p_imp<-ggplot(varImpPLS[[Hlabel[k]]], top=10)
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-VarImp_by_OTU.pdf",sep=""),family="ArialMT", width=16, height=14)
		print(p_imp)
	dev.off()
	dfscores <- as.data.frame(plsda_multi[[Hlabel[k]]]$finalModel$scores[,])
	prop_var <- round(plsda_multi[[Hlabel[k]]]$finalModel$Xvar/sum(plsda_multi[[Hlabel[k]]]$finalModel$Xvar)*100, 1)
	forChull<-cbind(dfscores[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]],Labels=x_cln_noOut[[k]]$ID,Shapes1=Shape_numbers1[[Hlabel[k]]],Shapes2=Shape_numbers2[[Hlabel[k]]])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=`Comp 1`,y=`Comp 2`))
	g <- g + geom_polygon(data=afterChull, aes(x=`Comp 1`,y=`Comp 2`, fill = OTUS ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")	
	g <- g + labs(x = paste("Comp 1 (", prop_var[1], "%)", sep=""), y = paste("Comp 2 (", prop_var[2], "%)", sep=""))
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-plot.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()

	# only Males
	control_M <- trainControl("repeatedcv", index = myfolds_M[[k]], selectionFunction = "oneSE")
	control_train_M<- trainControl("repeatedcv", index = myfolds_train_M[[k]], selectionFunction = "oneSE")
	plsda_multi_M[[Hlabel[k]]]<-train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",],y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],method="pls",tuneLength=30,trControl=control_M,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)
	plsda_multi_train_M[[Hlabel[k]]] <- train(x=training_M[[k]],y=training_OTUs_M[[k]],method="pls",tuneLength=30,trControl=control_train_M,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)
	plsPred_train_M[[Hlabel[k]]] <- predict(plsda_multi_train_M[[Hlabel[k]]], newdata = testing_M[[k]])
	cmPls_train_M[[Hlabel[k]]] <- caret::confusionMatrix(data=plsPred_train_M[[Hlabel[k]]], testing_OTUs_M[[k]])
	plsPred_M[[Hlabel[k]]] <- predict(plsda_multi_M[[Hlabel[k]]], newdata = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",])
	cmPls_M[[Hlabel[k]]] <- caret::confusionMatrix(data=plsPred_M[[Hlabel[k]]], x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"])
	
	p<-ggplot(plsda_multi_M[[k]])
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-CV-accuracy_by_Comp-M.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(p)
	dev.off()
	varImpPLS_M[[Hlabel[k]]]<-varImp(plsda_multi_M[[k]])
	p_imp<-ggplot(varImpPLS_M[[k]], top=10)
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-VarImp_by_OTU-M.pdf",sep=""),family="ArialMT", width=16, height=14)
		print(p_imp)
	dev.off()
	dfscores_M <- as.data.frame(plsda_multi_M[[k]]$finalModel$scores[,])
	prop_var_M <- round(plsda_multi_M[[k]]$finalModel$Xvar/sum(plsda_multi_M[[k]]$finalModel$Xvar)*100, 1)
	forChull<-cbind(dfscores_M[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="M"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=`Comp 1`,y=`Comp 2`))
	g <- g + geom_polygon(data=afterChull, aes(x=`Comp 1`,y=`Comp 2`, fill = OTUS ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star.
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = paste("Comp 1 (", prop_var_M[1], "%)", sep=""), y = paste("Comp 2 (", prop_var_M[2], "%)", sep=""))
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-plot-M.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()

	# only females
	control_F <- trainControl("repeatedcv", index = myfolds_F[[k]], selectionFunction = "oneSE")
	control_train_F<- trainControl("repeatedcv", index = myfolds_train_F[[k]], selectionFunction = "oneSE")
	plsda_multi_F[[Hlabel[k]]]<-train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",],y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],method="pls",tuneLength=30,trControl=control_F,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)	
	plsda_multi_train_F[[Hlabel[k]]] <- train(x=training_F[[k]],y=training_OTUs_F[[k]],method="pls",tuneLength=30,trControl=control_train_F,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)
	plsPred_train_F[[Hlabel[k]]] <- predict(plsda_multi_train_F[[Hlabel[k]]], newdata = testing_F[[k]])
	cmPls_train_F[[Hlabel[k]]] <- caret::confusionMatrix(data=plsPred_train_F[[Hlabel[k]]], testing_OTUs_F[[k]])
	plsPred_F[[Hlabel[k]]] <- predict(plsda_multi_F[[Hlabel[k]]], newdata = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",])
	plsPred_M[[Hlabel[k]]] <- predict(plsda_multi_M[[Hlabel[k]]], newdata = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",])
	cmPls_F[[Hlabel[k]]] <- caret::confusionMatrix(data=plsPred_F[[Hlabel[k]]], x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"])
	p<-ggplot(plsda_multi_F[[Hlabel[k]]])
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-CV-accuracy_by_Comp-F.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(p)
	dev.off()
	varImpPLS_F[[k]]<-varImp(plsda_multi_F[[Hlabel[k]]])
	p_imp<-ggplot(varImpPLS_F[[Hlabel[k]]], top=10)
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-VarImp_by_OTU-F.pdf",sep=""),family="ArialMT", width=16, height=14)
		print(p_imp)
	dev.off()
	dfscores_F <- as.data.frame(plsda_multi_F[[Hlabel[k]]]$finalModel$scores[,])
	prop_var_F <- round(plsda_multi_F[[Hlabel[k]]]$finalModel$Xvar/sum(plsda_multi_F[[Hlabel[k]]]$finalModel$Xvar)*100, 1)
	forChull<-cbind(dfscores_F[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="F"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=`Comp 1`,y=`Comp 2`))
	g <- g + geom_polygon(data=afterChull, aes(x=`Comp 1`,y=`Comp 2`, fill = OTUS ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star.
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = paste("Comp 1 (", prop_var_F[1], "%)", sep=""), y = paste("Comp 2 (", prop_var_F[2], "%)", sep=""))
	pdf(paste("./07_PLSDAs/01_PLSDA_multi_variables-",Hlabel[k],"/PLSDA-plot-F.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()
}

# end parallel processes
stopCluster(cl)

# save Confusion Matrix and Table of VIP
for(k in seq(Hlabel)){
	write.table(cmPls_train[[k]]$table,paste("01_tables/cmPls-",Hlabel[k],"-testing.txt",sep=""),quote = F, sep = "\t", col.names=NA)
	write.table(cmPls[[k]]$table,paste("01_tables/cmPls-",Hlabel[k],".txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(cmPls_train_M[[k]]$table,paste("01_tables/cmPls_M-",Hlabel[k],"-testing.txt",sep=""),quote = F, sep = "\t", col.names=NA)
	write.table(cmPls_M[[k]]$table,paste("01_tables/cmPls_M-",Hlabel[k],".txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(cmPls_train_F[[k]]$table,paste("01_tables/cmPls_F-",Hlabel[k],"-testing.txt",sep=""),quote = F, sep = "\t", col.names=NA)
	write.table(cmPls_F[[k]]$table,paste("01_tables/cmPls_F-",Hlabel[k],".txt",sep=""),quote = F,sep="\t", col.names=NA)

	varImpPLS_print<-cbind(varImpPLS[[k]]$importance,overall=apply(varImpPLS[[k]]$importance, 1, max))
	varImpPLS_M_print<-cbind(varImpPLS_M[[k]]$importance,overall=apply(varImpPLS_M[[k]]$importance, 1, max))
	varImpPLS_F_print<-cbind(varImpPLS_F[[k]]$importance,overall=apply(varImpPLS_F[[k]]$importance, 1, max))

	write.table(varImpPLS_print[order(-varImpPLS_print$overall),],paste("01_tables/varImpPLS-",Hlabel[k],".txt",sep=""),quote = F, sep = "\t", col.names=NA)
	write.table(varImpPLS_M_print[order(-varImpPLS_M_print$overall),],paste("01_tables/varImpPLS_M-",Hlabel[k],".txt",sep=""),quote = F, sep = "\t", col.names=NA)
	write.table(varImpPLS_F_print[order(-varImpPLS_F_print$overall),],paste("01_tables/varImpPLS_F-",Hlabel[k],".txt",sep=""),quote = F, sep = "\t", col.names=NA)
}

#######################################

##########################################################################################
# RF
##########################################################################################

#######################################

# try to keep e^Nf < No (Nf = number of features, No = number of observations)
# if you have to many variables or few observations... consider to reduce the variables
# compare values and think about running RF or not
for(k in seq(Hlabel)){
	print(exp(ncol(x_multi_cln_noOut_noZ_cln_noCor[[k]])))
	print(ncol(x_multi_cln_noOut_noZ_cln_noCor[[k]])*nrow(x_multi_cln_noOut_noZ_cln_noCor[[k]]))
}


# estimate optimal maximum number of variables based on number of samples
opt<-ncol(x_multi_cln_noOut_noZ_cln_noCor[[k]])
while(opt*nrow(x_multi_cln_noOut_noZ_cln_noCor[[k]]) <= exp(opt) ){
	opt = opt-1
} 
print(paste("The optimal number of variables for ", nrow(x_multi_cln_noOut_noZ_cln_noCor[[k]])," samples is ",opt, sep=""))

# you can select the "most important" variables by using RF, PLSDA, PCA or other methods.

#######################################

# set Y or N if you want to retrain the models
RC<-"Y"
if(RC=="Y"){
	rf_multi_F<-list()
	prox_scores_F<-list()
	prox_dist_F<-list()
	prox_mds_F<-list()	
	rf_multi_M<-list()
	prox_scores_M<-list()
	prox_dist_M<-list()	
	prox_mds_M<-list()
	rf_multi<-list()
	prox_scores<-list()
	prox_dist<-list()
	prox_mds<-list()
	rfPred_F<-list()
	cmRF_F<-list()
	rfPred_M<-list()
	cmRF_M<-list()
	rfPred<-list()
	cmRF<-list()
	varImpRF_F<-list()
	varImpRF_M<-list()
	varImpRF<-list()
	rf_multi_train<-list()
	rf_multi_train_F<-list()
	rf_multi_train_M<-list()
	control_train<-list()
	control_train_F<-list()
	control_train_M<-list()
}

#######################################
# start parallel processes
cl <- makeCluster(8,type = "SOCK")
registerDoSNOW(cl)

# plot RF for MULTI variables - NO OUTLIERS - polygon instead of ellipse 
for(k in seq(Hlabel)){
	dir.create(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],sep=""),showWarnings = FALSE)

	# whole dataset
	control <- trainControl("repeatedcv", index = myfolds[[k]], search = "grid", verboseIter=T)
	control_train<- trainControl("repeatedcv", index = myfolds_train[[k]], search = "grid", verboseIter=T)
	if(RC=="Y"){
		rf_multi[[Hlabel[k]]] <- caret::train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]], y=x_cln_noOut[[k]][,Hindex[k]], method="rf", tuneLength=20, trControl=control, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
		rf_multi_train[[Hlabel[k]]] <- caret::train(x=training[[k]], y=training_OTUs[[k]], method="rf", tuneLength=20, trControl=control_train, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
		prox_scores[[Hlabel[k]]]<-rf_multi[[Hlabel[k]]]$finalModel$proximity
		prox_dist[[Hlabel[k]]]<-dist(1-prox_scores[[Hlabel[k]]])
		prox_mds[[Hlabel[k]]]<-isoMDS(prox_dist[[Hlabel[k]]], k=2)
		# prox_mds[[Hlabel[k]]] <- cmdscale(prox_dist[[Hlabel[k]]], eig=TRUE, x.ret=TRUE) # use this line if you prefer classical (Torgerson) MDS instead.
		rfPred[[Hlabel[k]]] <- predict(rf_multi_train[[Hlabel[k]]], newdata = testing[[k]])
		cmRF[[Hlabel[k]]] <- caret::confusionMatrix(data=rfPred[[Hlabel[k]]], testing_OTUs[[k]])
	}
	p<-ggplot(rf_multi[[Hlabel[k]]])
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-CV-accuracy_by_Comp.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(p)
	dev.off()
	varImpRF[[Hlabel[k]]]<-varImp(rf_multi[[Hlabel[k]]])
	p_imp<-ggplot(varImpRF[[k]], top=10)
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-VarImp_by_OTU.pdf",sep=""),family="ArialMT", width=16, height=14)
		print(p_imp)
	dev.off()

	forChull<-data.frame(Dim1=prox_mds[[Hlabel[k]]]$points[,1],Dim2=prox_mds[[Hlabel[k]]]$points[,2],OTUS=x_cln_noOut[[k]][,Hindex[k]],Labels=x_cln_noOut[[k]]$ID,Shapes1=Shape_numbers1[[Hlabel[k]]],Shapes2=Shape_numbers2[[Hlabel[k]]])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=Dim1,y=Dim2))
	g <- g + geom_polygon(data=afterChull, aes(x=Dim1,y=Dim2, fill = `OTUS` ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star.
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = "isoMDS1", y = "isoMDS2" )
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-proximityPlot.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()

	# only females
	control_F <- trainControl("repeatedcv", index = myfolds_F[[k]], search = "grid", verboseIter=T)
	control_train_F<- trainControl("repeatedcv", index = myfolds_train_F[[k]], search = "grid", verboseIter=T)
	if(RC=="Y"){
		rf_multi_F[[Hlabel[k]]] <- caret::train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",], y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"], method="rf", tuneLength=20, trControl=control_F, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
		rf_multi_train_F[[Hlabel[k]]] <- caret::train(x=training_F[[k]], y=training_OTUs_F[[k]], method="rf", tuneLength=20, trControl=control_train_F, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
		prox_scores_F[[Hlabel[k]]]<-rf_multi_F[[Hlabel[k]]]$finalModel$proximity
		prox_dist_F[[Hlabel[k]]]<-dist(1-prox_scores_F[[Hlabel[k]]])
		prox_mds_F[[Hlabel[k]]]<-isoMDS(prox_dist_F[[Hlabel[k]]], k=2)
		# prox_mds_F[[Hlabel[k]]] <- cmdscale(prox_dist_F[[Hlabel[k]]], eig=TRUE, x.ret=TRUE) # use this line if you prefer classical (Torgerson) MDS instead.
		rfPred_F[[Hlabel[k]]] <- predict(rf_multi_train_F[[Hlabel[k]]], newdata = testing_F[[k]])
		cmRF_F[[Hlabel[k]]] <- caret::confusionMatrix(data=rfPred_F[[Hlabel[k]]], testing_OTUs_F[[k]])
	}
	p<-ggplot(rf_multi_F[[Hlabel[k]]])
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-CV-accuracy_by_Comp-F.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(p)
	dev.off()
	varImpRF_F[[Hlabel[k]]]<-varImp(rf_multi_F[[Hlabel[k]]])
	p_imp<-ggplot(varImpRF_F[[k]], top=10)
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-VarImp_by_OTU-F.pdf",sep=""),family="ArialMT", width=16, height=14)
		print(p_imp)
	dev.off()
	forChull<-data.frame(Dim1=prox_mds_F[[Hlabel[k]]]$points[,1],Dim2=prox_mds_F[[Hlabel[k]]]$points[,2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="F"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=Dim1,y=Dim2))
	g <- g + geom_polygon(data=afterChull, aes(x=Dim1,y=Dim2, fill = `OTUS` ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star.
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = "isoMDS1", y = "isoMDS2" )
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-proximityPlot-F.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()

	# only males
	control_M <- trainControl("repeatedcv", index = myfolds_M[[k]], search = "grid", verboseIter=T)
	control_train_M<- trainControl("repeatedcv", index = myfolds_train_M[[k]], search = "grid", verboseIter=T)
	if(RC=="Y"){
		rf_multi_M[[Hlabel[k]]] <- caret::train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",], y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"], method="rf", tuneLength=20, trControl=control_M, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
		rf_multi_train_M[[Hlabel[k]]] <- caret::train(x=training_M[[k]], y=training_OTUs_M[[k]], method="rf", tuneLength=20, trControl=control_train_M, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
		prox_scores_M[[Hlabel[k]]]<-rf_multi_M[[Hlabel[k]]]$finalModel$proximity
		prox_dist_M[[Hlabel[k]]]<-dist(1-prox_scores_M[[Hlabel[k]]])
		prox_mds_M[[Hlabel[k]]]<-isoMDS(prox_dist_M[[Hlabel[k]]], k=2)
		# prox_mds_M[[Hlabel[k]]] <- cmdscale(prox_dist_M[[Hlabel[k]]], eig=TRUE, x.ret=TRUE) # use this line if you prefer classical (Torgerson) MDS instead.
		rfPred_M[[Hlabel[k]]] <- predict(rf_multi_M[[Hlabel[k]]], newdata = testing_M[[k]])
		cmRF_M[[Hlabel[k]]] <- caret::confusionMatrix(data=rfPred_M[[Hlabel[k]]], testing_OTUs_M[[k]])
	}
	p<-ggplot(rf_multi_M[[Hlabel[k]]])
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-CV-accuracy_by_Comp-M.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(p)
	dev.off()
	varImpRF_M[[Hlabel[k]]]<-varImp(rf_multi_M[[Hlabel[k]]])
	p_imp<-ggplot(varImpRF_M[[k]], top=10)
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-VarImp_by_OTU-M.pdf",sep=""),family="ArialMT", width=16, height=14)
		print(p_imp)
	dev.off()

	forChull<-data.frame(Dim1=prox_mds_M[[Hlabel[k]]]$points[,1],Dim2=prox_mds_M[[Hlabel[k]]]$points[,2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="M"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"])
	forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
	afterChull<-vector()
	for(i in 1:length(OTUs[[k]])){
		afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
	}
	g <- ggplot(data=forChull, aes(x=Dim1,y=Dim2))
	g <- g + geom_polygon(data=afterChull, aes(x=Dim1,y=Dim2, fill = `OTUS` ), alpha = 0.3) 
	g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes1, fill=OTUS), size=3)
	g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes2, fill=OTUS), size=3) # I'm plotting the points twice only to overlap the triagles to produce a star.
	# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
	g <- g + scale_shape_manual(values = myShapes)
	g <- g + scale_colour_manual(values = myColors)
	g <- g + scale_color_discrete(name = '')
	g <- g + colScale + fillScale
	g <- g + theme(legend.direction = 'horizontal', legend.position = 'top')
	g <- g + guides(shape="none")
	g <- g + labs(x = "isoMDS1", y = "isoMDS2" )
	pdf(paste("./08_RF/01_RF_multi_variables-",Hlabel[k],"/RF-proximityPlot-M.pdf",sep=""), family="ArialMT", width=16, height=14)
		print(g)
	dev.off()

}

# end parallel processes
stopCluster(cl)

# save Confusion Matrix and Table of VIP
for(k in seq(Hlabel)){
	write.table(cmRF[[k]]$table,paste("01_tables/cmRF-",Hlabel[k],"-testing.txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(rf_multi[[k]]$finalModel$confusion,paste("01_tables/cmRF-",Hlabel[k],".txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(cmRF_M[[k]]$table,paste("01_tables/cmRF_M-",Hlabel[k],"-testing.txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(rf_multi_M[[k]]$finalModel$confusion,paste("01_tables/cmRF_M-",Hlabel[k],".txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(cmRF_F[[k]]$table,paste("01_tables/cmRF_F-",Hlabel[k],"-testing.txt",sep=""),quote = F,sep="\t", col.names=NA)
	write.table(rf_multi_F[[k]]$finalModel$confusion,paste("01_tables/cmRF_F-",Hlabel[k],".txt",sep=""),quote = F,sep="\t", col.names=NA)
}


#######################################

##########################################################################################
# COMPOSED PDF WITH SELECTED MULTIVARIATED PLOTS
##########################################################################################
# Here I'm running the analyses again for consistency, but "g" objects can be previously saved and used here if you wish...
#######################################

# composed pdf
# start parallel processes
cl <- makeCluster(8,type = "SOCK")
registerDoSNOW(cl)


# plot PCA for RATIOS, CATEGORICAL and DISC variables - polygon instead of ellipse 
for(k in seq(Hlabel)){
	dir.create(paste("./09_Combined/01_combined-",Hlabel[k],sep=""),showWarnings = FALSE)

		# plot pca for females - ALL variables
		pca_multi_F <- prcomp(x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln$Sex=="F",], scale.=T)
		dfscores_F <- as.data.frame(pca_multi_F$x)
		prop_pca_F = pca_multi_F$sdev^2/sum(pca_multi_F$sdev^2)
		forChull<-cbind(dfscores_F[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="F"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"])
		forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
		afterChull<-vector()
		for(i in 1:length(OTUs[[k]])){
			afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
		}
		g <- ggplot(data=forChull, aes(x=`PC1`,y=`PC2`))
		g <- g + geom_polygon(data=afterChull, aes(x=`PC1`,y=`PC2`, fill = `OTUS` ), alpha = 0.3) 
		g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes1, fill=OTUS), size=1.5)
		g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes2, fill=OTUS), size=1.5) # I'm plotting the points twice only to overlap the triagles to produce a star
		# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
		g <- g + scale_shape_manual(values = myShapes)
		g <- g + scale_colour_manual(values = myColors)
		g <- g + scale_color_discrete(name = '')
		g <- g + colScale + fillScale
		g <- g + theme(legend.position = 'none')
		g <- g + guides(shape="none")
		g <- g + labs(x = paste("PC1 (", round(prop_pca_F[1]*100,1), "%)", sep=""), y = paste("PC2 (", round(prop_pca_F[2]*100,1), "%)", sep=""))
		g1<-g

		# PSLDA only females
		control_F <- trainControl("repeatedcv", index = myfolds_F[[k]], selectionFunction = "oneSE")
		control_train_F<- trainControl("repeatedcv", index = myfolds_train_F[[k]], selectionFunction = "oneSE")
		plsda_multi_F[[Hlabel[k]]]<-train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",],y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],method="pls",tuneLength=30,trControl=control_F,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)	
		dfscores_F <- as.data.frame(plsda_multi_F[[Hlabel[k]]]$finalModel$scores[,])
		prop_var_F <- round(plsda_multi_F[[Hlabel[k]]]$finalModel$Xvar/sum(plsda_multi_F[[Hlabel[k]]]$finalModel$Xvar)*100, 1)
		forChull<-cbind(dfscores_F[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="F"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"])
		forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
		afterChull<-vector()
		for(i in 1:length(OTUs[[k]])){
			afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
		}
		g <- ggplot(data=forChull, aes(x=`Comp 1`,y=`Comp 2`))
		g <- g + geom_polygon(data=afterChull, aes(x=`Comp 1`,y=`Comp 2`, fill = OTUS ), alpha = 0.3) 
		g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes1, fill=OTUS), size=1.5)
		g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes2, fill=OTUS), size=1.5) # I'm plotting the points twice only to overlap the triagles to produce a star.
		# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
		g <- g + scale_shape_manual(values = myShapes)
		g <- g + scale_colour_manual(values = myColors)
		g <- g + scale_color_discrete(name = '')
		g <- g + colScale + fillScale
		g <- g + theme(legend.position = 'none')
		g <- g + guides(shape="none")
		g <- g + labs(x = paste("Comp 1 (", prop_var_F[1], "%)", sep=""), y = paste("Comp 2 (", prop_var_F[2], "%)", sep=""))
		g2<-g

		# RF only females
		control_F <- trainControl("repeatedcv", index = myfolds_F[[k]], search = "grid", verboseIter=T)
		control_train_F<- trainControl("repeatedcv", index = myfolds_train_F[[k]], search = "grid", verboseIter=T)
		if(RC=="Y"){
			rf_multi_F[[Hlabel[k]]] <- caret::train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",], y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"], method="rf", tuneLength=20, trControl=control_F, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
			rf_multi_train_F[[Hlabel[k]]] <- caret::train(x=training_F[[k]], y=training_OTUs_F[[k]], method="rf", tuneLength=20, trControl=control_train_F, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
			prox_scores_F[[Hlabel[k]]]<-rf_multi_F[[Hlabel[k]]]$finalModel$proximity
			prox_dist_F[[Hlabel[k]]]<-dist(1-prox_scores_F[[Hlabel[k]]])
			prox_mds_F[[Hlabel[k]]]<-isoMDS(prox_dist_F[[Hlabel[k]]], k=2)
			# prox_mds_F[[Hlabel[k]]] <- cmdscale(prox_dist_F[[Hlabel[k]]], eig=TRUE, x.ret=TRUE) # use this line if you prefer classical (Torgerson) MDS instead.
			rfPred_F[[Hlabel[k]]] <- predict(rf_multi_train_F[[Hlabel[k]]], newdata = testing_F[[k]])
			cmRF_F[[Hlabel[k]]] <- caret::confusionMatrix(data=rfPred_F[[Hlabel[k]]], testing_OTUs_F[[k]])
		}
		forChull<-data.frame(Dim1=prox_mds_F[[Hlabel[k]]]$points[,1],Dim2=prox_mds_F[[Hlabel[k]]]$points[,2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="F"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="F"])
		forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
		afterChull<-vector()
		for(i in 1:length(OTUs[[k]])){
			afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
		}
		g <- ggplot(data=forChull, aes(x=Dim1,y=Dim2))
		g <- g + geom_polygon(data=afterChull, aes(x=Dim1,y=Dim2, fill = `OTUS` ), alpha = 0.3) 
		g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes1, fill=OTUS), size=1.5)
		g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes2, fill=OTUS), size=1.5) # I'm plotting the points twice only to overlap the triagles to produce a star.
		# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
		g <- g + scale_shape_manual(values = myShapes)
		g <- g + scale_colour_manual(values = myColors)
		g <- g + scale_color_discrete(name = '')
		g <- g + colScale + fillScale
		g <- g + theme(legend.position = 'none')
		g <- g + guides(shape="none")
		g <- g + labs(x = "isoMDS1", y = "isoMDS2" )
		g3<-g

		# plot pca for males - ALL variables
		pca_multi_M <- prcomp(x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln$Sex=="M",], scale.=T)
		dfscores_M <- as.data.frame(pca_multi_M$x)
		prop_pca_M = pca_multi_M$sdev^2/sum(pca_multi_M$sdev^2)
		forChull<-cbind(dfscores_M[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="M"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"])
		forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
		afterChull<-vector()
		for(i in 1:length(OTUs[[k]])){
			afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
		}
		g <- ggplot(data=forChull, aes(x=`PC1`,y=`PC2`))
		g <- g + geom_polygon(data=afterChull, aes(x=`PC1`,y=`PC2`, fill = `OTUS` ), alpha = 0.3) 
		g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes1, fill=OTUS), size=1.5)
		g <- g + geom_point(data=forChull, aes(x=`PC1`,y=`PC2`, color = OTUS, shape=Shapes2, fill=OTUS), size=1.5) # I'm plotting the points twice only to overlap the triagles to produce a star
		# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
		g <- g + scale_shape_manual(values = myShapes)
		g <- g + scale_colour_manual(values = myColors)
		g <- g + scale_color_discrete(name = '')
		g <- g + colScale + fillScale
		g <- g + theme(legend.position = 'none')
		g <- g + guides(shape="none")
		g <- g + labs(x = paste("PC1 (", round(prop_pca_M[1]*100,1), "%)", sep=""), y = paste("PC2 (", round(prop_pca_M[2]*100,1), "%)", sep=""))
		g4<-g
		
		# PSLDA only Males
		control_M <- trainControl("repeatedcv", index = myfolds_M[[k]], selectionFunction = "oneSE")
		control_train_M<- trainControl("repeatedcv", index = myfolds_train_M[[k]], selectionFunction = "oneSE")
		plsda_multi_M[[Hlabel[k]]]<-train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",],y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],method="pls",tuneLength=30,trControl=control_M,preProc=c("zv","center","scale"),metric="Accuracy", maxit = 10000, verbose=T)
		dfscores_M <- as.data.frame(plsda_multi_M[[k]]$finalModel$scores[,])
		prop_var_M <- round(plsda_multi_M[[k]]$finalModel$Xvar/sum(plsda_multi_M[[k]]$finalModel$Xvar)*100, 1)
		forChull<-cbind(dfscores_M[,1:2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="M"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"])
		forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
		afterChull<-vector()
		for(i in 1:length(OTUs[[k]])){
			afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
		}
		g <- ggplot(data=forChull, aes(x=`Comp 1`,y=`Comp 2`))
		g <- g + geom_polygon(data=afterChull, aes(x=`Comp 1`,y=`Comp 2`, fill = OTUS ), alpha = 0.3) 
		g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes1, fill=OTUS), size=1.5)
		g <- g + geom_point(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, shape=Shapes2, fill=OTUS), size=1.5) # I'm plotting the points twice only to overlap the triagles to produce a star.
		# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
		g <- g + scale_shape_manual(values = myShapes)
		g <- g + scale_colour_manual(values = myColors)
		g <- g + scale_color_discrete(name = '')
		g <- g + colScale + fillScale
		g <- g + theme(legend.position = 'none')
		g <- g + guides(shape="none")
		g <- g + labs(x = paste("Comp 1 (", prop_var_M[1], "%)", sep=""), y = paste("Comp 2 (", prop_var_M[2], "%)", sep=""))
		g5<-g

		# RF only males
		control_M <- trainControl("repeatedcv", index = myfolds_M[[k]], search = "grid", verboseIter=T)
		control_train_M<- trainControl("repeatedcv", index = myfolds_train_M[[k]], search = "grid", verboseIter=T)
		if(RC=="Y"){
			rf_multi_M[[Hlabel[k]]] <- caret::train(x=x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",], y=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"], method="rf", tuneLength=20, trControl=control_M, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
			rf_multi_train_M[[Hlabel[k]]] <- caret::train(x=training_M[[k]], y=training_OTUs_M[[k]], method="rf", tuneLength=20, trControl=control_train_M, preProc=c("zv","center","scale"), metric="Accuracy", maxit = 5000, proximity=T, verbose=T)
			prox_scores_M[[Hlabel[k]]]<-rf_multi_M[[Hlabel[k]]]$finalModel$proximity
			prox_dist_M[[Hlabel[k]]]<-dist(1-prox_scores_M[[Hlabel[k]]])
			prox_mds_M[[Hlabel[k]]]<-isoMDS(prox_dist_M[[Hlabel[k]]], k=2)
			# prox_mds_M[[Hlabel[k]]] <- cmdscale(prox_dist_M[[Hlabel[k]]], eig=TRUE, x.ret=TRUE) # use this line if you prefer classical (Torgerson) MDS instead.
			rfPred_M[[Hlabel[k]]] <- predict(rf_multi_M[[Hlabel[k]]], newdata = testing_M[[k]])
			cmRF_M[[Hlabel[k]]] <- caret::confusionMatrix(data=rfPred_M[[Hlabel[k]]], testing_OTUs_M[[k]])
		}
		forChull<-data.frame(Dim1=prox_mds_M[[Hlabel[k]]]$points[,1],Dim2=prox_mds_M[[Hlabel[k]]]$points[,2],OTUS=x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],Labels=x_cln_noOut[[k]]$ID[x_cln_noOut[[k]]$Sex=="M"],Shapes1=Shape_numbers1[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"],Shapes2=Shape_numbers2[[Hlabel[k]]][x_cln_noOut[[k]]$Sex=="M"])
		forChull$OTUS<-factor(forChull$OTUS,all_OTUs[all_OTUs%in%forChull$OTUS])
		afterChull<-vector()
		for(i in 1:length(OTUs[[k]])){
			afterChull<-rbind(afterChull,forChull[forChull[,3]==OTUs[[k]][i],][chull(forChull[forChull[,3]==OTUs[[k]][i],1:2]),])
		}
		g <- ggplot(data=forChull, aes(x=Dim1,y=Dim2))
		g <- g + geom_polygon(data=afterChull, aes(x=Dim1,y=Dim2, fill = `OTUS` ), alpha = 0.3) 
		g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes1, fill=OTUS), size=1.5)
		g <- g + geom_point(data=forChull, aes(x=Dim1, y=Dim2, color = OTUS, shape=Shapes2, fill=OTUS), size=1.5) # I'm plotting the points twice only to overlap the triagles to produce a star.
		# g <- g + geom_text(data=forChull, aes(x=`Comp 1`,y=`Comp 2`, color = OTUS, label=Labels), hjust = -0.5, nudge_x = 0.0) # uncomment this line if you want to check the lables associated to each point (IDs)
		g <- g + scale_shape_manual(values = myShapes)
		g <- g + scale_colour_manual(values = myColors)
		g <- g + scale_color_discrete(name = '')
		g <- g + colScale + fillScale
		g <- g + theme(legend.position = 'none')
		g <- g + guides(shape="none")
		g <- g + labs(x = "isoMDS1", y = "isoMDS2" )
		g6<-g

		pdf(paste("./09_Combined/01_combined-",Hlabel[k],"/Combined_plot.pdf",sep=""), family="ArialMT", width=8.3, height=10.7)
			grid.arrange(g1, g4, g2, g5, g3, g6, nrow = 3, ncol = 2)
		dev.off()
}

# end parallel processes
stopCluster(cl)

#######################################

##########################################################################################
# COMPOSED PDF WITH SELECTED BOXPLOTS
##########################################################################################

#######################################

# composed pdf
# Boxplot + histogram for VIP - per Sex (if dimorphic)- per OTUs
# Need to run/load test for sexual dimorphism first

vip<-c("disc_IOLd","disc_DOM","disc_VEN")

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
# EQUALIZE SAMPLE SIZE
##########################################################################################
# boostraping OTUs per Sex for ANOVA-like approaches
# OBS1: Each OTU for each Sex will be randomly sampled 30 times
# OBS2: the complete dataset will be the simple combination of the bootstraped dataset for each sex
#######################################

bootData<-list()
bootOTUS<-list()
runs<-c("run01","run02","run03")
for(q in seq(runs)){
	sexo<-c("F","M")
	set.seed(1001+q)
	for(k in seq(Hlabel)){
		id_comb<-list()
		for(j in seq(sexo)){
			id_zero<-vector()
			id_boot<-vector()
			for(i in seq(OTUs[[k]])){
				id_zero<-which(x_cln[,Hindex[k]]==OTUs[[k]][i]&x_cln$Sex==sexo[j])
				id_boot<-sample(as.character(id_zero),20,replace=TRUE)
				id_comb[[sexo[j]]]<-c(id_comb[[sexo[j]]],id_boot)
			}
		bootData[[runs[q]]][[Hlabel[k]]][[sexo[j]]]<-x_multi_cln_noOut_noZ_cln_noCor[[k]][as.integer(id_comb[[sexo[j]]]),]
		bootOTUS[[runs[q]]][[Hlabel[k]]][[sexo[j]]]<-x_cln[,Hindex[k]][as.integer(id_comb[[sexo[j]]])]
		}
		bootOTUS[[runs[q]]][[Hlabel[k]]][["ALL"]]<-x_cln[,Hindex[k]][c(as.integer(id_comb[["F"]]),as.integer(id_comb[["M"]]))]
		bootData[[runs[q]]][[Hlabel[k]]][["ALL"]]<-x_multi_cln_noOut_noZ_cln_noCor[[k]][c(as.integer(id_comb[["F"]]),as.integer(id_comb[["M"]])),]
	}
}

#######################################


##########################################################################################
# MANOVA NP
##########################################################################################
# OBS01: by combining discrete and continuous variables, your dataset probably cannot be analysed by MANOVA. Use PerMANOVA instead.
# OBS02: small/unequal sample size among OTUs can be devastating for any ANOVA-like approach. Consider using an equalized dataset.

#######################################

# backup_noCor<-x_multi_cln_noOut_noZ_cln_noCor
# x_multi_cln_noOut_noZ_cln_noCor<-x_multi_cln_noOut_noZ_cln
# x_multi_cln_noOut_noZ_cln_noCor<-backup_noCor

#######################################

# Manova NP for original dataset

out_permanova_multi <- list()
out_permanova_multi_PHoc2 <- list()
out_permanova_multi_F <- list()
out_permanova_multi_PHoc2_F <- list()
out_permanova_multi_M <- list()
out_permanova_multi_PHoc2_M <- list()

for(k in seq(Hlabel)){
	# PerMANOVA on RATIOS and DISC variables (Multi)
	out_permanova_multi[[Hlabel[k]]]<-adonis(formula = x_multi_cln_noOut_noZ_cln_noCor[[k]] ~ x_cln_noOut[[k]][,Hindex[k]], data = x_multi_cln_noOut_noZ_cln_noCor[[k]], permutations = 1000, method = "euclidean") 
	# pairwise Post-hoc PerMANOVA
	out_permanova_multi_PHoc2[[Hlabel[k]]]<-pairwise.perm.manova(dist(x_multi_cln_noOut_noZ_cln_noCor[[k]],"euclidean"),x[,Hindex[k]],nperm=1000, p.method = "bonferroni")
	write.table(out_peranova_multi[[Hlabel[k]]]$aov.tab,paste("01_tables/PerMANOVA_Multi-",Hlabel[k],".txt",sep=""),quote=F, sep="\t", col.names=NA)
	write.table(out_permanova_multi_PHoc2[[Hlabel[k]]]$p.value,paste("01_tables/PerMANOVA_Multi-PHoc2-",Hlabel[k],".txt",sep=""),quote=F, sep="\t", col.names=NA)
	#######################################
	# PerMANOVA on RATIOS and DISC variables (Multi) - Sex
	out_permanova_multi_F[[Hlabel[k]]]<-adonis(formula = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",] ~ x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"], data = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",], permutations = 1000, method = "euclidean") 
	out_permanova_multi_M[[Hlabel[k]]]<-adonis(formula = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",] ~ x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"], data = x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",], permutations = 1000, method = "euclidean") 
	# pairwise Post-hoc PerMANOVA
	out_permanova_multi_PHoc2_F[[Hlabel[k]]]<-pairwise.perm.manova(dist(x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="F",],"euclidean"),x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="F"],nperm=1000, p.method = "bonferroni")
	out_permanova_multi_PHoc2_M[[Hlabel[k]]]<-pairwise.perm.manova(dist(x_multi_cln_noOut_noZ_cln_noCor[[k]][x_cln_noOut[[k]]$Sex=="M",],"euclidean"),x_cln_noOut[[k]][,Hindex[k]][x_cln_noOut[[k]]$Sex=="M"],nperm=1000, p.method = "bonferroni")
	write.table(out_permanova_multi_F[[Hlabel[k]]]$aov.tab,paste("01_tables/PerMANOVA_Multi-F-",Hlabel[k],".txt",sep=""),quote=F, sep="\t", col.names=NA)
	write.table(out_permanova_multi_PHoc2_F[[Hlabel[k]]]$p.value,paste("01_tables/PerMANOVA_Multi-PHoc2-F-",Hlabel[k],".txt",sep=""),quote=F, sep="\t", col.names=NA)
	write.table(out_permanova_multi_M[[Hlabel[k]]]$aov.tab,paste("01_tables/PerMANOVA_Multi-M-",Hlabel[k],".txt",sep=""),quote=F, sep="\t", col.names=NA)
	write.table(out_permanova_multi_PHoc2_M[[Hlabel[k]]]$p.value,paste("01_tables/PerMANOVA_Multi-PHoc2-M-",Hlabel[k],".txt",sep=""),quote=F, sep="\t", col.names=NA)
}

#######################################

# Manova NP for equalized dataset

out_permanova_multi_eq <- list()
out_permanova_multi_PHoc2_eq <- list()

# running in triplicate to avoid biases when sampling OTUs with large N

runs<-c("run01","run02","run03")
for(q in seq(runs)){
	for(q in seq(runs)){
		sexo<-c("ALL","F","M")	
		for(k in seq(Hlabel)){
			for(j in seq(sexo)){
				# PerMANOVA on RATIOS and DISC variables (Multi)
				out_permanova_multi_eq[[runs[q]]][[Hlabel[k]]][[sexo[j]]]<-adonis(formula = bootData[[runs[q]]][[Hlabel[k]]][[sexo[j]]] ~ bootOTUS[[runs[q]]][[Hlabel[k]]][[sexo[j]]], data = bootData[[runs[q]]][[Hlabel[k]]][[sexo[j]]], permutations = 1000, method = "euclidean") 
				# pairwise Post-hoc PerMANOVA
				out_permanova_multi_PHoc2_eq[[runs[q]]][[Hlabel[k]]][[sexo[j]]]<-pairwise.perm.manova(dist(bootData[[runs[q]]][[Hlabel[k]]][[sexo[j]]],"euclidean"),bootOTUS[[runs[q]]][[Hlabel[k]]][[sexo[j]]],nperm=1000, p.method = "bonferroni")
			}
		}
	}
}


out_permanova_multi_PHoc2_eq_sex<-list()
out_permanova_multi_PHoc2_eq_matrix<-list()


for(k in seq(Hlabel)){
	for(j in seq(sexo)){
		for(q in 1:length(runs)){
			out_permanova_multi_PHoc2_eq_matrix[[runs[q]]]<-as.matrix(out_permanova_multi_PHoc2_eq[[runs[q]]][[Hlabel[k]]][[sexo[j]]][["p.value"]])
		}
	out_permanova_multi_PHoc2_eq_sex[[sexo[j]]]<-Reduce("+", out_permanova_multi_PHoc2_eq_matrix)/3
	}
	table_to_write<-do.call(rbind, out_permanova_multi_PHoc2_eq_sex)[order(sequence(sapply(out_permanova_multi_PHoc2_eq_sex, nrow))), ]
	write.table(table_to_write,paste("01_tables/PerMANOVA_Multi-PHoc2-eq-",Hlabel[k],"-2.txt",sep=""),quote=F, sep="\t", col.names=NA)
}

for(k in seq(Hlabel)){
	write.table(out_permanova_multi_PHoc2_eq[[runs[q]]][[Hlabel[k]]][[sexo[j]]]$aov.tab,paste("01_tables/PerMANOVA_Multi-eq.txt",sep=""),quote=F, sep="\t", col.names=NA)
	write.table(out_permanova_multi_PHoc2_eq[[runs[q]]][[Hlabel[k]]][[sexo[j]]]$p.value,paste("01_tables/PerMANOVA_Multi-PHoc2-eq.txt",sep=""),quote=F, sep="\t", col.names=NA)
}

#######################################

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
