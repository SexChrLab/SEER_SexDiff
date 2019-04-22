######################################################################
#This script subsets the SEER data sets by males and females for each cancer type and performs a linear regression for the cancer types that occur in both sexes. 
#
#created: October 29, 2018
#updated: March 11, 2019
#
#created by: Karen Funderburk 
# kmfunder@asu.edu
######################################################################

types = read.csv('cancer_types.csv', sep = ",", header = T)$Type		#reading in list of cancer types
types = as.character(types)

fit_results = matrix(nrow = length(types), ncol = 17)
rownames(fit_results) = types
colnames(fit_results) = c('Regression Slope', 'Regression P value', 'Residual Standard Error', 'R squared', 'Shapiro Wilk', 'Male beginning', 'Male ending', 'Female beginning', 'Female ending', 'Male Percent Difference', 'Female Percent Difference', 'Beginning difference', 'Ending difference', 'Males Only Slope', 'Males Only P Value', 'Females Only Slope', 'Females Only P Value')


pdf(file = "CancerTypePlotsSubsetBySexls.pdf")

par(mfrow = c(3, 2), mai = c(.65,0.5,0.5,.5))		#initializing plot area
for (type in types){
	data = read.csv(paste0('./Data_files/', type, '.csv'), skip = 2, header = T, as.is = T)	#reading in data files
	data = data[,1:7]			#removing extra columns

	data[data == '-'] = NA		#replacing dashes with NAs

	data[,4] = as.numeric(data[,4])
	data[,5] = as.numeric(data[,5]) 			#converting to numerics
	data[,6] = as.numeric(data[,6])
	data[,7] = as.numeric(data[,7])

	#subsetting data
	idata = subset(data, (Sex == 'Female' | Sex =='Male') & Race.Ethnicity == 'All Races (includes Hispanic)' & Age == 'All Ages', select = c('Year', 'Sex', 'Age.adjusted..AA..Rate.per.100.000', 'AA.Rate.Lower.Confidence.Interval', 'AA.Rate.Upper.Confidence.Interval'))
	idata = aggregate(idata[,3:5], by = list(Year = idata$Year, Sex = idata$Sex), mean)		#combining multiple entries per year
		
	y_upper = max(idata$Age.adjusted..AA..Rate.per.100.000, na.rm = T) 
		
	if (length(idata$Sex[idata$Sex=='Male'])==0 | type=='Breast'){		#female only cancers
		
		plot(idata$Year[idata$Sex=='Female'], idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female'],type = 'l', lwd=2, cex = 1.2, col = 'red', ylim = c(0, y_upper), ylab = 'Age Adjusted Rater per 100,000', xlab = 'Year', main = type)
		plot.new()
		
		fit_results[which(types==type), 8] = head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1)
		fit_results[which(types==type), 9] = tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1)
		fit_results[which(types==type), 11] = (tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1) - head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1))/tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1)
		
	} else if (length(idata$Sex[idata$Sex=='Female'])==0){		#male only cancers
		
		plot(idata$Year[idata$Sex=='Male'], idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male'], type = 'l', lwd = 2, col = 'blue', ylim = c(0, y_upper), ylab = 'Age Adjusted Rater per 100,000', xlab = 'Year', main = type)
		plot.new()
		
		fit_results[which(types==type), 6] = head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1)
		fit_results[which(types==type), 7] = tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1)
		fit_results[which(types==type), 10] = (tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1) - head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1))/tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1)
		
	} else {
	
	idata$adj.Year = idata$Year - 1975		#shifting years to start at 0 so intercept from fitting is more accurate
		
#	m_ratio = data.frame(Year = idata$Year[idata$Sex=='Male'], Ratio = (idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']/(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female'] + idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male'])))		#calculating male fraction
	
	m_ratio = data.frame(Year = idata$Year[idata$Sex=='Male'], Ratio = (idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']/(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female'] )))		#calculating male to female ratio

#	myfit = lm(log(m_ratio$Ratio/(1-m_ratio$Ratio)) ~ m_ratio$Year)  #logit regression
	myfit = lm(m_ratio$Ratio ~ m_ratio$Year)		#least squares regression

	#individual fits for males and females

	m_fit = lm(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male'] ~ idata$Year[idata$Sex=='Male'])
	f_fit = lm(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female'] ~ idata$Year[idata$Sex=='Female'])
	
	
	fit_results[which(types==type), 1] = summary(myfit)$coefficients[2,1]	#slope of regression line
	fit_results[which(types==type), 2] = summary(myfit)$coefficients[2,4]	#pvalue of regression line
	fit_results[which(types==type), 3] = summary(myfit)$sigma						#Residual standard error for fit
	fit_results[which(types==type), 4] = summary(myfit)$r.squared					#R squared for fit
	fit_results[which(types==type), 5] = shapiro.test(residuals(myfit))$p.value		#shapiro wilk for myfit
	fit_results[which(types==type), 6] = head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1)
	fit_results[which(types==type), 7] = tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1)
	fit_results[which(types==type), 8] = head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1)
	fit_results[which(types==type), 9] = tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1)
	fit_results[which(types==type), 10] = (tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1) - head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1))/tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1)
	fit_results[which(types==type), 11] = (tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1) - head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1))/tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1)
	fit_results[which(types==type), 12] = abs(head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1) - head(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1))
	fit_results[which(types==type), 13] = abs(tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male']), n=1) - tail(na.omit(idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female']), n=1))

	fit_results[which(types==type), 14] = summary(m_fit)$coefficients[2,1]	#slope of male regression line
	fit_results[which(types==type), 15] = summary(m_fit)$coefficients[2,4]	#pval of male regression line
	fit_results[which(types==type), 16] = summary(f_fit)$coefficients[2,1]	#slope of female regression line
	fit_results[which(types==type), 17] = summary(f_fit)$coefficients[2,4]	#pval of female regression line

	#plotting data with fit results
	plot(idata$Year[idata$Sex=='Female'], idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Female'], type = 'l', lwd = 2, cex = 1.2, col = 'red', ylim = c(0, y_upper), ylab = 'Age Adjusted Rater per 100,000', xlab = 'Year', main = type)
	lines(idata$Year[idata$Sex=='Male'], idata$Age.adjusted..AA..Rate.per.100.000[idata$Sex=='Male'], lwd=2, col = 'blue')
#	lines(idata$Year[idata$Sex=='Male'], m_fit$fitted.values, lwd = 1.5)
#	lines(idata$Year[idata$Sex=='Female'], f_fit$fitted.values, lwd = 1.5)	
	legend('topleft', legend = c('female', 'male'), col = c('red', 'blue'), pch = c(18, 16), bg = 'white', cex =.7)

ymin = min(m_ratio$Ratio, na.rm = TRUE)
ymax = max(m_ratio$Ratio, na.rm = TRUE)

	#plotting male to female ratio
	plot(m_ratio[is.na(m_ratio$Ratio) == F,], pch = 18, col = 'orange', ylim = c(ymin, ymax), xlab = "Year", ylab = "Male to Female Ratio", main = 'Male to Female Ratio', sub = paste0('Slope: ', formatC(summary(myfit)$coefficients[2,1], format = 'e'), ', pval: ', formatC(summary(myfit)$coefficients[2,4], format = 'e')))
	lines(m_ratio$Year[is.na(m_ratio$Ratio) == F], myfit$fitted.values)		#ls lines
#	lines(m_ratio$Year[is.na(m_ratio$Ratio) == F], exp(myfit$fitted.values)/(1+exp(myfit$fitted.values)))
	
# Residual Analysis
	
	# res_m = residuals(indiv_fit_m)
	# res_norm_m = shapiro.test(res_m)
	# fit_results[which(types==type), 15] = res_norm_m$p.value 		#p value of shapiro wilk test to test for normality of residuals

	# res_f = residuals(indiv_fit_f)
	# res_norm_f = shapiro.test(res_f)
	# fit_results[which(types==type), 16] = res_norm_f$p.value 		#p value of shapiro wilk test to test for normality of residuals

	# qqnorm(res_m, main = 'Male Normal Q-Q Plot', sub = paste('pval:', res_norm_m$p.value))			#qqplot for male residuals
	# qqline(res_m)
	
	# qqnorm(res_f, main = 'Female Normal Q-Q Plot', sub = paste('pval:', res_norm_f$p.value))			#qqplot for female residuals
	# qqline(res_f)
	

	#Levene's test to test for homoskedasticity
	# sex = c(rep('female', length(res_f)), rep('male', length(res_m)))
	# res_both = c(res_f, res_m)
	# res = residuals( reg_fit)
	# res$sex = as.factor(res$sex)
	#lt = leveneTest(res_both~sex, data = res)
	#fit_results[which(types==type), 13] = lt$P[1]

	#plotting residuals and fitted values
	# plot(indiv_fit_m$fitted.values, res_m, main = "Male Residuals vs Fitted Values", ylab = 'Residuals', xlab = "Fitted Values")
	# abline(h=0)
	# plot(indiv_fit_f$fitted.values, res_f, main = "Female Residuals vs Fitted Values", ylab = 'Residuals', xlab = "Fitted Values")
	# abline(h=0)


	}
	
		
}
dev.off()


pdf(file = 'ls_slopes.pdf', height = 10, width = 12)
	#sorting differences between slopes and making plot
	par(mai = c(2.5, 1, .5, .5), cex.axis = 0.85)
	fit_results = fit_results[order(fit_results[,'Regression Slope']),]
	signif = rep(0, nrow(fit_results))
	signif[fit_results[,'Regression P value']<0.05] = 1			#adding column to identify significant results
	fit_results = cbind(fit_results, signif)
	barplot(fit_results[which(fit_results[,'Regression Slope'] != 'NA' ),'Regression Slope'], las = 2, col = signif + 3, ylab = "Slope of Regression Line", main = "Slope of Regression Line for Each Cancer Type")
	legend('topleft', legend = c('significant difference', 'insignificant difference'), fill = c(4,3), bty = 'n')
	
dev.off()
	
#	formatC(fit_results, digits = 3, format = 'G')

 fit_results = fit_results[order(rownames(fit_results)),]


#printing data for latex table
	
for (i in 1:nrow(fit_results)){
	if(is.na(fit_results[i,1]) != 1 ){
		cat(rownames(fit_results)[i], ' & ', round(fit_results[i,6],2), ' & ', round(fit_results[i,8],2), ' & ',  round(fit_results[i,7],2), ' & ', round(fit_results[i,9],2), ' & ', round(fit_results[i,10]*100,1), ' & ', round(fit_results[i,11]*100,1), '\\', '\n')
	}
}




#printing regression results for latex table	
	
for (i in 1:nrow(fit_results)){
	if(is.na(fit_results[i,1]) != 1 ){
		cat(rownames(fit_results)[i], ' & ', round(fit_results[i,12], 3), ' & ', round(fit_results[i,13], 3), ' & ')
		
		if(fit_results[i,15] < 0.05){ cat('\\textbf{', round(fit_results[i,14], 4), '} & ')}
		else if (fit_results[i,15] > 0.05){cat(round(fit_results[i,14],4), ' & ')}
		
		if(fit_results[i,17] < 0.05){ cat('\\textbf{', round(fit_results[i,16], 4), '} & ')}
		else if (fit_results[i,17] > 0.05){cat(round(fit_results[i,16],4), ' & ')}
		
		if(fit_results[i,2] < 0.05){cat('\\textbf{', formatC(fit_results[i,1], digits = 3, format = 'G'), '} \\ ', '\n')} 
		else if (fit_results[i,2] > 0.05){cat(formatC(fit_results[i,1], digits = 3, format = 'G'), ' \\', '\n')}
	}
}

#fit_results = formatC(fit_results)

write.csv(fit_results, file = "LS_RegressionResults.csv")




