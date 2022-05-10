
library(truncnorm)
library(plotrix)

########################################################################
#                                                                      #
# Monte-Carlo simulations for:                                         #
# Prediction of prophylactic replacement of voice                      #
# prosthesis in laryngectomized patients –                             #
# A retrospective cohort study                                         #
#                                                                      #
# This is the R-code for the actual Monte-Carlo                        #
# simulations. It works stand alone.                                   #
#                                                                      #
# Copyright 2022                                                       #
# R.J.J.H. van Son & Netherlands Cancer Institute,                     #
# Plesmanlaan 121                                                      #
# 1066CX Amsterdam, The Netherlands                                    #
#                                                                      #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or    #
# (at your option) any later version.                                  #
#                                                                      #
# This program is distributed in the hope that it will be useful,      #
# but WITHOUT ANY WARRANTY; without even the implied warranty of       #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the        #
# GNU General Public License for more details.                         #
#                                                                      #
# You should have received a copy of the GNU General Public License    #
# along with this program.  If not, see http://www.gnu.org/licenses/   #
#                                                                      #
########################################################################

# Monte Carlo Simulation Parameters

# Desired Leak fraction
#ideal.leak <- 0.12;
ideal.leak <- 0.30;
#ideal.leak <- 0.05;

# The fixed cutoff for the fixed change policy
fixed.cutoff <- 16;

# Device life parameters
mean.num.days <- 60;
sd.num.days <- 60;
minimum.num.days <- 7;

# Simulation settings
# Initial look-back window
window <- 3;

# Maximum look back window
window.dyn <- 10;

# The MC parameters
# Number of devices to follow
Nprostheses <- 32;

# Number of patients to test
Npatient <- 5000;

######################################################
# 
# Select distribution of device life-times
# All distributions are definite positive
#
######################################################

# Select the distribution to test
# Poisson distribution, this should be the default!
#Distribution <- "Poisson";
#Distribution <- "Uniform";
Distribution <- "trunc. Normal";
#Distribution <- "Exponential";

# Poisson distribution, this should be the default!
# Separate adaptations for positive and negative
mean.sd.text <- "Lambda=mean ($\\pm$ Sqrt(Lambda)";
overshoot.exp <- 1.5;

# Second adaptive run
overshoot2.exp <- 1;
mean.sd.text <- "Lambda";

# Initialize distribution parameters
# Uniform distribution
if(Distribution == "Uniform"){
	mean.sd.text <- "mean [lambda/2, lambda*3/2]";
}

# Normal distribution
if(Distribution == "trunc. Normal"){
	standard.dev <- 1/2;
	#standard.dev <- 2/3;
	mean.sd.text <- "mean ($\\pm$mean/2)";
}

# Exponential distribution with SD = Mean
if(Distribution == "Exponential"){
	standard.dev <- 1;
	mean.sd.text <- "mean ($\\pm$ mean)";
}

standard.dev.2 <- 0.8;

# Reporting table
results.header <- c("none", "none.leak", "ideal", "ideal.leak", "ideal.cutoff", "fixed", "fixed.leak", "fixed.cutoff", "initial.tot", "initial.leak.tot", "initial.red", "initial.leak.red", "initial.cutoff", "adapt.tot", "adapt.leak.tot", "adapt.red", "adapt.leak.red", "adapt.cutoff", "adapt2.tot", "adapt2.leak.tot", "adapt2.red", "adapt2.leak.red", "adapt2.cutoff", "sample", "window");
results.matrix<-matrix(ncol=length(results.header), byrow=FALSE);
colnames(results.matrix) <- results.header;
results.matrix <- data.frame(results.matrix);

# Accumulators for plots
sum.sequence <- rep(0, Nprostheses);
sum.sequence2 <- rep(0, Nprostheses);
lambda.list <- c(0);

# The Monte-Carlo simulation
for(lambda in rexp(Npatient, rate=1/mean.num.days)){
	# Lower cutoff of the average life-time of a device
    while(lambda <= minimum.num.days)lambda <- rexp(Npatient, rate=1/mean.num.days)[[1]];
    
    ######################################################
    # 
	# Generate survival times
    #
    ######################################################
    if(Distribution == "Uniform"){
		survival <- runif(Nprostheses, min = lambda/2, max = lambda*3/2);
		standard.dev <- sd(survival)/mean(survival);
    } else if(Distribution == "trunc. Normal"){
		survival <- rtruncnorm(Nprostheses, a=0, b=Inf, mean = lambda, sd = lambda*standard.dev);
    } else if(Distribution == "Exponential"){
		survival <- rexp(Nprostheses, rate = 1/lambda);
		standard.dev <- 1;
	} else {			# Poisson
		survival <- rpois(Nprostheses, lambda=lambda);
		standard.dev <- sqrt(lambda)/lambda;
	};
    results <- c();
    lambda.list <- append(lambda.list, lambda);
    
    ######################################################
    # 
    # Wait-to-Leakage policy
    #
    ######################################################
	mean.none <- sum(survival)/length(survival);
	leak.none <- 1;
	results <- c(results, mean.none, 100*leak.none);
	
    ######################################################
    # 
	# Ideal change policy with known parameters
    #
    ######################################################
	cutoff.ideal <- quantile(survival, probs=ideal.leak);
	survival.ideal <- survival;
	survival.ideal[survival>=cutoff.ideal] <- cutoff.ideal;
	mean.ideal <- sum(survival.ideal)/length(survival.ideal);
	leak.ideal <- sum(survival.ideal<cutoff.ideal)/length(survival.ideal);
	results <- c(results, mean.ideal, 100*leak.ideal, cutoff.ideal);
	
    ######################################################
    # 
	# Fixed-time changes policy
    #
    ######################################################
	cutoff.fixed <- fixed.cutoff;
	survival.fixed <- survival;
	survival.fixed[survival>cutoff.fixed] <- cutoff.fixed;
	mean.fixed <- sum(survival.fixed)/length(survival.fixed);
	leak.fixed <- sum(survival.fixed<cutoff.fixed)/length(survival.fixed);
	results <- c(results, mean.fixed, 100*leak.fixed, cutoff.fixed);

    ######################################################
    # 
	# Initial window
    #
    ######################################################

	# Extract initial window and determine cutoff
	window.initial <- survival[1:window];
	#cutoff.initial <- quantile(window.initial, probs=ideal.leak);
	cutoff.initial <- qtruncnorm(ideal.leak, a=0, b=Inf, mean=mean(window.initial), sd=sd(survival[window.initial]))	;
	
	# Extract remainder and apply cutoff to reduced list
	survival.initial <- survival[(window+1):length(survival)];
	survival.initial[survival.initial>cutoff.initial] <- cutoff.initial;
	mean.initial.red <- sum(survival.initial)/length(survival.initial);
	leak.initial.red <- sum(survival.initial<cutoff.initial)/length(survival.initial);
	leak.initial.tot <- (window + sum(survival.initial<cutoff.initial))/(window+length(survival.initial));
	
	# Complete the list with the initial survivial times
	survival.initial <- c(window.initial, survival.initial);
	# Calculate mean survival
	mean.initial.tot <- sum(survival.initial)/length(survival.initial);

	results <- c(results, mean.initial.tot, 100*leak.initial.tot, mean.initial.red, 100*leak.initial.red, cutoff.initial);
	
    ######################################################
    # 
	# Adaptive 1
    #
    ######################################################
	current.window <- window;
	# Reduce the ideal leakage target to "undershoot" the first part
	current.leakage <- ideal.leak^1.5;
	
	# Extract initial window and determine cutoff
	cutoff.adapt <- survival;
	survival.adapt <- survival;
	cutoff.adapt[1:window] <- survival.adapt[1:window]+0.1;
	
	# We do not know the distribution, so we use a truncated normal for the initial window
	ideal.cutoff <- qtruncnorm(current.leakage, a=0, b=Inf, mean=mean(survival[1:window]), sd=sd(survival[1:window])) 
	sd.initial <- sd(survival[1:window]);
	if(sd.initial < 1)sd.initial <- 1;
	expected.negative <- mean(rtruncnorm(100, a=0, b=ideal.cutoff, mean=mean(survival[1:window]), sd=sd.initial)) 

	negative.trend <- (ideal.cutoff - expected.negative);
	positive.trend <- negative.trend;
	
	# Initialize
	current.cutoff <- ideal.cutoff;
	minus.values <- c(negative.trend);
	for(i in (current.window+1):length(survival.adapt)){
		cutoff.adapt[i] <- current.cutoff;
		
		# The end point still depends on the adapt coeff.
		if(survival.adapt[i] >= current.cutoff) {
			# Add the current value to the list
			survival.adapt[i] <- current.cutoff;
			
			# Adapt next positive cutoff to previous leakages
			if(length(minus.values) > 0)positive.trend <- minus.values[length(minus.values)];
			
			current.adapt <- rnorm(1, mean=positive.trend, sd=positive.trend/2.5);
		} else {
			current.adapt <- 0;
			minus.values <- c(minus.values, current.cutoff - survival.adapt[i]);
		};
		# New, adapted, cutoff
		current.cutoff <- quantile(survival.adapt[(i-current.window+1):(i)], probs=current.leakage);
		current.cutoff <- current.cutoff + current.adapt;
		
		# Lengthen the window
		if(current.window < window.dyn){
			current.window <- current.window+1;
		} else {
			current.leakage <- ideal.leak;
		};
		
	};

	# Collect statistics
	mean.adapt.tot <- sum(survival.adapt)/length(survival.adapt);
	leak.adapt.tot <- sum(survival.adapt < cutoff.adapt)/length(survival.adapt);
	mean.adapt.red <- sum(survival.adapt[round((Nprostheses+1)/2):Nprostheses])/length(survival.adapt[round((Nprostheses+1)/2):Nprostheses]);
	leak.adapt.red <- sum(survival.adapt[round((Nprostheses+1)/2):Nprostheses] < cutoff.adapt[round((Nprostheses+1)/2):Nprostheses])/length(survival.adapt[round((Nprostheses+1)/2):Nprostheses]);
	sum.sequence <- sum.sequence + as.numeric(survival.adapt < cutoff.adapt);;

	results <- c(results, mean.adapt.tot, 100*leak.adapt.tot, mean.adapt.red, 100*leak.adapt.red, mean(survival.adapt[round((Nprostheses+1)/2):Nprostheses]));
	
    ######################################################
    # 
	# Adaptive 2
    #
    ######################################################

	current.window <- window;
	# Reduce the ideal leakage target to "undershoot" the first part
	current.leakage <- ideal.leak;
	
	# Extract initial window and determine cutoff
	cutoff.adapt2 <- survival;
	survival.adapt2 <- survival;
	cutoff.adapt2[1:window] <- survival.adapt2[1:window]+0.1;
	
	# We do not know the distribution, so we use a truncated normal for the initial window
	ideal.cutoff <- qtruncnorm(current.leakage, a=0, b=Inf, mean=mean(survival[1:window]), sd=sd(survival[1:window])) 
	sd.initial <- sd(survival[1:window]);
	if(sd.initial < 1)sd.initial <- 1;
	expected.negative <- mean(rtruncnorm(100, a=0, b=ideal.cutoff, mean=mean(survival[1:window]), sd=sd.initial)) 

	negative.trend <- (ideal.cutoff - expected.negative);
	positive.trend <- negative.trend;
	
	# Initialize
	current.cutoff <- ideal.cutoff;
	minus.values <- c(negative.trend);
	for(i in (current.window+1):length(survival.adapt2)){
		cutoff.adapt2[i] <- current.cutoff;
		
		# The end point still depends on the adapt coeff.
		if(survival.adapt2[i] >= current.cutoff) {
			# Add the current value to the list
			survival.adapt2[i] <- current.cutoff;
			
			# Adapt next positive cutoff to previous leakages
			if(length(minus.values) > 0)positive.trend <- minus.values[length(minus.values)];
			
			current.adapt2 <- rnorm(1, mean=positive.trend, sd=positive.trend/2.5);
		} else {
			current.adapt2 <- 0;
			minus.values <- c(minus.values, current.cutoff - survival.adapt2[i]);
		};
		# New, adapted, cutoff
		current.cutoff <- quantile(survival.adapt2[(i-current.window+1):(i)], probs=current.leakage);
		current.cutoff <- current.cutoff + current.adapt2;
		
		# Lengthen the window
		if(current.window < window.dyn){
			current.window <- current.window+1;
		} else {
			current.leakage <- ideal.leak;
		};
		
	};

	# Collect statistics
	mean.adapt2.tot <- sum(survival.adapt2)/length(survival.adapt2);
	leak.adapt2.tot <- sum(survival.adapt2 < cutoff.adapt2)/length(survival.adapt2);
	mean.adapt2.red <- sum(survival.adapt2[round((Nprostheses+1)/2):Nprostheses])/length(survival.adapt2[round((Nprostheses+1)/2):Nprostheses]);
	leak.adapt2.red <- sum(survival.adapt2[round((Nprostheses+1)/2):Nprostheses] < cutoff.adapt2[round((Nprostheses+1)/2):Nprostheses])/length(survival.adapt2[round((Nprostheses+1)/2):Nprostheses]);
	sum.sequence2 <- sum.sequence2 + as.numeric(survival.adapt2 < cutoff.adapt2);;

	results <- c(results, mean.adapt2.tot, 100*leak.adapt2.tot, mean.adapt2.red, 100*leak.adapt2.red, mean(survival.adapt2[round((Nprostheses+1)/2):Nprostheses]));

    ######################################################
    # 
	# Final parameters
    #
    ######################################################
	results <- c(results, Nprostheses, window);
	
	results.matrix <- rbind(results, results.matrix, make.row.names = FALSE);
}

######################################################
# 
# Generate statistics
#
######################################################
results.matrix <- results.matrix[-c(length(results.matrix[,1])),];
averaged.results <- colMeans(results.matrix, na.rm=TRUE);
sd.results <- sapply(results.matrix, sd);

averaged.PercentLeak <- c(
100/(averaged.results["none"]/(averaged.results["none.leak"]/100)/averaged.results["none"]),
100/(averaged.results["ideal"]/(averaged.results["ideal.leak"]/100)/averaged.results["none"]),
100/(averaged.results["fixed"]/(averaged.results["fixed.leak"]/100)/averaged.results["none"]),
100/(averaged.results["initial.tot"]/(averaged.results["initial.leak.tot"]/100)/averaged.results["none"]),
100/(averaged.results["initial.red"]/(averaged.results["initial.leak.red"]/100)/averaged.results["none"]),
100/(averaged.results["adapt.tot"]/(averaged.results["adapt.leak.tot"]/100)/averaged.results["none"]),
100/(averaged.results["adapt.red"]/(averaged.results["adapt.leak.red"]/100)/averaged.results["none"]),
100/(averaged.results["adapt2.tot"]/(averaged.results["adapt2.leak.tot"]/100)/averaged.results["none"]),
100/(averaged.results["adapt2.red"]/(averaged.results["adapt2.leak.red"]/100)/averaged.results["none"])
);

averaged.PercentChanges <- c(
100*averaged.results["none"]/averaged.results["none"],
100*averaged.results["none"]/averaged.results["ideal"],
100*averaged.results["none"]/averaged.results["fixed"],
100*averaged.results["none"]/averaged.results["initial.tot"],
100*averaged.results["none"]/averaged.results["initial.red"],
100*averaged.results["none"]/averaged.results["adapt.tot"],
100*averaged.results["none"]/averaged.results["adapt.red"],
100*averaged.results["none"]/averaged.results["adapt2.tot"],
100*averaged.results["none"]/averaged.results["adapt2.red"]
);
names(averaged.PercentChanges) <- names(averaged.PercentLeak);

mean.lambda <- mean(lambda.list);
sd.lambda <- sd(lambda.list);

# Print results (when not using Rmd)
cat("Means\n");
print(averaged.results);

cat("\nStandard Deviation\n");
print(sd.results);

# Create graph data
mean.leakage <- sum.sequence / Npatient;
mean.leakage2 <- sum.sequence2 / Npatient;

######################################################
# 
# Sequence of leakage targets
# Using full knowledge of distribution (Ideal)
#
# This part is pretty pointless:
# Mean device life (% of mean) ~ Target leakage%
#
######################################################

cutoff.ideal.list <- c();
mean.ideal.list <- c();
leakage.ideal.list <- c();
for(leakage in 1:75/100){
	lambda <- mean.num.days;
    Nideal <- 10000;
    ######################################################
    # 
	# Generate survival times
    #
    ######################################################
    if(Distribution == "Uniform"){
		survival <- runif(Nideal, min = lambda/2, max = lambda*3/2);
    } else if(Distribution == "trunc. Normal"){
		survival <- rtruncnorm(Nideal, a=0, b=Inf, mean = lambda, sd = lambda/3);
    } else if(Distribution == "Exponential"){
		survival <- rexp(Nideal, rate = 1/lambda);
	} else {			# Poisson
		survival <- rpois(Nideal, lambda=lambda);
	};
	
    ######################################################
    # 
	# Ideal
    #
    ######################################################
	cutoff.ideal <- quantile(survival, probs=leakage);
	survival.ideal <- survival;
	survival.ideal[survival>=cutoff.ideal] <- cutoff.ideal;
	mean.ideal <- sum(survival.ideal)/length(survival.ideal);
	leak.ideal <- sum(survival.ideal<cutoff.ideal)/length(survival.ideal);
	
	cutoff.ideal.list <- c(cutoff.ideal.list, cutoff.ideal);
	mean.ideal.list <- c(mean.ideal.list, mean.ideal);
	leakage.ideal.list <- c(leakage.ideal.list, leak.ideal);
}
names(mean.ideal.list) <- names(cutoff.ideal.list);
names(leakage.ideal.list) <- names(cutoff.ideal.list);


######################################################
# 
# Sequence of standard deviations
# Using full knowledge of distribution (Ideal)
# and ideal.leak leakage as the target
#
######################################################

# Settings
norm.leakage <- ideal.leak;
lambda <- mean.num.days;

######################################################
# 
# Normal distribution
#
######################################################
Nideal.normalsd <- 1000000;

# Initialization
ideal.normalsd.list <- c(c(1,2,5)/100, c(1:10, 12, 15, 20)/10);
cutoff.normalsd.list <- c();
realt2f.normalsd.list <- c();
realsd.normalsd.list <- c();
mean.normalsd.list <- c();
leakage.normalsd.list <- c();
for(ideal.standard.dev in ideal.normalsd.list){
    ######################################################
    # 
	# Generate survival times
	# abs(Normal distribution)
    #
    ######################################################
	survival <- rtruncnorm(Nideal.normalsd, a=0, b=Inf, mean = lambda, sd = lambda*ideal.standard.dev);
	
    ######################################################
    # 
	# Ideal
    #
    ######################################################
	cutoff.ideal <- quantile(survival, probs=norm.leakage);
	names(cutoff.ideal) <- c(sprintf("%.0f%%", 100*ideal.standard.dev));
	
	survival.ideal <- survival;
	survival.ideal[survival>=cutoff.ideal] <- cutoff.ideal;
	mean.ideal <- sum(survival.ideal)/length(survival.ideal);
	leak.ideal <- sum(survival.ideal<cutoff.ideal)/length(survival.ideal);
	
	cutoff.normalsd.list <- c(cutoff.normalsd.list, cutoff.ideal);
	mean.normalsd.list <- c(mean.normalsd.list, mean.ideal);
	leakage.normalsd.list <- c(leakage.normalsd.list, leak.ideal);
	realt2f.normalsd.list <- c(realt2f.normalsd.list, mean(survival));
	realsd.normalsd.list <- c(realsd.normalsd.list, sd(survival));
}
names(mean.normalsd.list) <- names(cutoff.normalsd.list);
names(leakage.normalsd.list) <- names(cutoff.normalsd.list);
names(realt2f.normalsd.list) <- names(realt2f.normalsd.list);
names(realsd.normalsd.list) <- names(realsd.normalsd.list);


######################################################
# 
# Uniform distribution
#
######################################################
Nideal.uniform <- 1000000;

# Initialization
ideal.uniform.list <- c(c(1,2,5)/100, 1:10/10);
cutoff.uniform.list <- c();
realt2f.uniform.list <- c();
realsd.uniform.list <- c();
mean.uniform.list <- c();
leakage.uniform.list <- c();
for(ideal.standard.dev in ideal.uniform.list){
    ######################################################
    # 
	# Generate survival times
	# abs(Normal distribution)
    #
    ######################################################
	survival <- runif(Nideal.uniform, min = (1-ideal.standard.dev)*lambda, max = (1+ideal.standard.dev)*lambda);
	while(sum(survival <= 0) > 0){
		survival[survival <= 0] <- runif(Nideal.uniform, min = (1-ideal.standard.dev)*lambda, max = (1+ideal.standard.dev)*lambda);
	}
	
    ######################################################
    # 
	# Ideal
    #
    ######################################################
	cutoff.ideal <- quantile(survival, probs=norm.leakage);
	names(cutoff.ideal) <- c(sprintf("%.0f%%", 100*ideal.standard.dev));
	
	survival.ideal <- survival;
	survival.ideal[survival>=cutoff.ideal] <- cutoff.ideal;
	mean.ideal <- sum(survival.ideal)/length(survival.ideal);
	leak.ideal <- sum(survival.ideal<cutoff.ideal)/length(survival.ideal);
	
	cutoff.uniform.list <- c(cutoff.uniform.list, cutoff.ideal);
	mean.uniform.list <- c(mean.uniform.list, mean.ideal);
	leakage.uniform.list <- c(leakage.uniform.list, leak.ideal);
	realt2f.uniform.list <- c(realt2f.uniform.list, mean(survival));
	realsd.uniform.list <- c(realsd.uniform.list, sd(survival));
}
names(mean.uniform.list) <- names(cutoff.uniform.list);
names(leakage.uniform.list) <- names(cutoff.uniform.list);
names(realt2f.uniform.list) <- names(realt2f.uniform.list);
names(realsd.uniform.list) <- names(realsd.uniform.list);

######################################################
# 
# Poisson distribution
#
######################################################
Nideal.poisson <- 1000000;

# Initialization
ideal.poisson.list <- c(3:12)^2;
cutoff.poisson.list <- c();
realt2f.poisson.list <- c();
realsd.poisson.list <- c();
mean.poisson.list <- c();
leakage.poisson.list <- c();
for(lambda.poisson in ideal.poisson.list){
    ######################################################
    # 
	# Generate survival times
	# abs(Normal distribution)
    #
    ######################################################
	survival <- rpois(Nideal.poisson, lambda=lambda.poisson);
	
    ######################################################
    # 
	# Ideal
    #
    ######################################################
	cutoff.ideal <- quantile(survival, probs=norm.leakage);
	#names(cutoff.ideal) <- c(sprintf("%.0f%%", 100*standard.dev));
	
	survival.ideal <- survival;
	survival.ideal[survival>=cutoff.ideal] <- cutoff.ideal;
	mean.ideal <- sum(survival.ideal)/length(survival.ideal);
	leak.ideal <- sum(survival.ideal<cutoff.ideal)/length(survival.ideal);
	
	cutoff.poisson.list <- c(cutoff.poisson.list, cutoff.ideal);
	mean.poisson.list <- c(mean.poisson.list, mean.ideal);
	leakage.poisson.list <- c(leakage.poisson.list, leak.ideal);
	realt2f.poisson.list <- c(realt2f.poisson.list, mean(survival));
	realsd.poisson.list <- c(realsd.poisson.list, sd(survival));
}
names(mean.poisson.list) <- names(cutoff.poisson.list);
names(leakage.poisson.list) <- names(cutoff.poisson.list);
names(realt2f.poisson.list) <- names(realt2f.poisson.list);
names(realsd.poisson.list) <- names(realsd.poisson.list);

# Do the same for Exponential
survival <- rexp(Nideal.normalsd, rate = 1/lambda);
cutoff.ideal <- quantile(survival, probs=norm.leakage);
survival.ideal <- survival;
survival.ideal[survival>=cutoff.ideal] <- cutoff.ideal;

# Statistics
mean.normalsd.Exponential <- sum(survival.ideal)/length(survival.ideal);
leakage.normalsd.Exponential <- sum(survival.ideal<cutoff.ideal)/length(survival.ideal);
realt2f.normalsd.Exponential <- mean(survival);
realsd.normalsd.Exponential <- sd(survival);


######################################################
# 
# Graphics: Figure 1
#
######################################################

# Plot Figure 1
imageList <- c();
fileName.Fig.1 <- "Time-to-changeOverSD.png";
OutputWidth <-  8.267;
OutputHeight <- 8.267;

png(fileName.Fig.1, width=OutputWidth, height=OutputHeight, units="in", res=300);

par(cex=1.6);
cex <- par("cex");
cex.factor <- 1;
symbolsize <- 1 * cex.factor;
labelsize <- 1.0 * cex.factor;
axissize <- 1 * cex.factor;
namessize <- axissize
legendsize <- (1+axissize)/2;
cexlab <- cex;

par(mar=c(5, 5, 4, 5) + 0.1);
par(family="Helvetica")

labelsize <- 1 * cex.factor;
plot(realsd.normalsd.list/realt2f.normalsd.list, 100*mean.normalsd.list/(realt2f.normalsd.list*norm.leakage), pch=21, xlim=c(0,1), ylim=c(100,100/norm.leakage), ylab="Time-to-leakage (index) -> ", xlab="SD/Mean", cex.lab=labelsize, cex.axis=axissize, bg="gray70");
mtext("<- # Changes (index)", side=4, cex=1.5, line=2);
axis(1, at=c(1:10)/10, labels=FALSE, tck=-0.02);
axis(4, at=c(10000/(50*2:(floor(2/norm.leakage)))/(norm.leakage)), labels=c(paste(50*2:(floor(2/norm.leakage))))) ;
abline(h=100, lty=2, lwd=1);
abline(h=100/norm.leakage, lty=2, lwd=1);
abline(line(realsd.normalsd.list/realt2f.normalsd.list, 100*mean.normalsd.list/(realt2f.normalsd.list*norm.leakage), iter = 1), lty=3);
abline(line(realsd.uniform.list/realt2f.uniform.list, 100*mean.uniform.list/(realt2f.uniform.list*norm.leakage), iter = 1), lty=3);
abline(line(c(0, 1), c(100/norm.leakage, 100*mean.normalsd.Exponential/(realt2f.normalsd.Exponential*norm.leakage)), iter = 1), lty=3);
points(realsd.poisson.list/realt2f.poisson.list, 100*mean.poisson.list/(realt2f.poisson.list*norm.leakage), pch=24, col="red3", bg="red3", cex=0.8);
points(realsd.uniform.list/realt2f.uniform.list, 100*mean.uniform.list/(realt2f.uniform.list*norm.leakage), pch=23, col="green3", bg="green3", cex=0.8);
points(realsd.normalsd.list/realt2f.normalsd.list, 100*mean.normalsd.list/(realt2f.normalsd.list*norm.leakage), pch=21, col="black", bg="grey70", cex=0.8);
points(realsd.normalsd.Exponential/realt2f.normalsd.Exponential, 100*mean.normalsd.Exponential/(realt2f.normalsd.Exponential*norm.leakage), pch=22, col="blue3", bg="blue3");
text(0.1,100, labels=c(sprintf("N: 10^%.0f", log10(Nideal.normalsd))), pos=3, cex=0.9)

normalsd.line <- line(realsd.normalsd.list/realt2f.normalsd.list, 100*mean.normalsd.list/(realt2f.normalsd.list*norm.leakage), iter = 1);
y.leak.time <-  coefficients(normalsd.line)[2] * standard.dev + coefficients(normalsd.line)[1];
x.norm.sd <- (y.leak.time - coefficients(normalsd.line)[1])/coefficients(normalsd.line)[2];
y.leak.time.2 <-  coefficients(normalsd.line)[2] * standard.dev.2 + coefficients(normalsd.line)[1];
x.norm.sd.2 <- (y.leak.time.2 - coefficients(normalsd.line)[1])/coefficients(normalsd.line)[2];

if(y.leak.time > 100){
	text(0.1, y.leak.time, labels=c(sprintf("#leaks = %.0f%%", 100/(y.leak.time/100))), pos=1, cex=0.7, col="blue")
	segments(-0.1, y.leak.time, 1.1, y.leak.time, lty=3, lwd=2, col="blue");
	mtext(sprintf("%.0f",  y.leak.time), 2, at=c( y.leak.time), col="blue", cex=0.9);
	mtext(sprintf("%.0f",  100/(y.leak.time/100*norm.leakage)), 4, at=c( y.leak.time), col="blue", cex=0.9);
	
	# Crossing of the 0.5 SD/Mean line with NormalSD line
	segments(standard.dev, 0, standard.dev, y.leak.time, lty=3, lwd=2, col="blue");
	mtext(sprintf("%.2f", standard.dev), 1, at=c(standard.dev), col="blue", cex=0.8);
	
	# Crossing of the 0.8 SD/Mean line NormalSD line
	segments(standard.dev.2, 0, standard.dev.2, y.leak.time.2, lty=3, lwd=2, col="red");
	mtext(sprintf("%.2f", standard.dev.2), 1, at=c(standard.dev.2), col="blue", cex=0.8);
	draw.circle(standard.dev.2, y.leak.time.2, border="red", radius=0.075, lwd=2)
};

legend("topright", title = sprintf("Distribution; leak: %.0f%% (DL%.0f)", 100*norm.leakage, 100*(1-norm.leakage)), legend=c("Uniform", "Normal (truncated)", "Poisson (mean=sd^2)", "Exponential (mean=sd)"), col=c("green3", "black", "red3", "blue3"), lty=0, lwd=0, pch=c(23, 21, 24, 22), pt.bg=c("green3", "grey70", "red3", "blue3"), bty="o", cex=0.6, box.lty=0, bg="white");
segments(1, 100/norm.leakage, 1.1, 100/norm.leakage, lty=2, lwd=1);

dev.off();

######################################################
# 
# Graphics: Figure 2
#
######################################################

# Plot Figure 2
imageList <- c();
fileName.Fig.2 <- "LeakageOverChanges.png";
OutputWidth <-  11.692;
OutputHeight <- 8.267;

cex.factor <- 1;
symbolsize <- 1 * cex.factor;
labelsize <- 1 * cex.factor;
axissize <- 1 * cex.factor;
namessize <- axissize
legendsize <- (1+axissize)/2;
png(fileName.Fig.2, width=OutputWidth, height=OutputHeight, units="in", res=300);
par(cex=2);
cex=1.5;
cexlab <- 2*cex;

par(mar=c(5, 5, 4, 2) + 0.1);
par(family="Helvetica")
plot(100*mean.leakage2, type='p', pch=22, bg="grey", col="grey", lty=1, xlab="Change #", ylab="Leakge -> %", ylim=c(0,100), lwd=2, cex.lab=labelsize, cex.axis=axissize);
lines(100*mean.leakage2, type='l', pch=22, bg="grey", col="grey", lty=1, lwd=2, cex.lab=labelsize, cex.axis=axissize);
points(100*mean.leakage, type='p', pch=21, bg="black", col="black", lty=1, lwd=2, cex.lab=labelsize, cex.axis=axissize);
lines(100*mean.leakage, type='l', pch=21, bg="black", col="black", lty=1, lwd=2, cex.lab=labelsize, cex.axis=axissize);
axis(2, at=c(1:10)*10, labels=FALSE, tck=-0.02);
abline(h=100*ideal.leak, lty=2, lwd=2);
text(2,100*ideal.leak, labels=c(sprintf("%.0f%%", 100*ideal.leak)), pos=1);
text(27,0, labels=c(sprintf("N: %.0f", Npatient)), pos=2, cex=1)
text(0, 0, labels=c(Distribution), pos=4, cex=0.8)
legend("topright", title = sprintf("Overshoot exp. (%.0f#)", window.dyn), legend=c(sprintf("target^%.1f", 1), sprintf("target^%.1f", 1.5)), col=c("grey", "black"), lty=1, lwd=2, pch=c( 22,21), pt.bg=c("grey", "black"), bty="n");

dev.off();
