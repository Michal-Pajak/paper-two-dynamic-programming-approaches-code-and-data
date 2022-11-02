model_comparison_baseline_3_cases <- function()
{
	#Common and fixed
	t0 <- 100
	k1 <- 50
	k2 <- 50
	e <- 1

	#STATIONARY
	s_a<-1.4
	s_b<-0.1093
	s_c<-0.02
	s_d<-0.175
	s_x0<-11
	s_y0<-10

	#Discrete system - Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	s_d_trajx_0<-array(0,dim=t0+1)
	s_d_trajy_0<-array(0,dim=t0+1)

	s_d_trajx_0[1]<-s_x0
	s_d_trajy_0[1]<-s_y0

	for (t in 2:(t0+1)){
			x<-s_d_trajx_0[t-1]
			y<-s_d_trajy_0[t-1]
			xnew<-round(x+x*s_a*(1-x/k1)-x*s_b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-round(ynew+ynew*s_c*x*(1-e*ynew/k2)-s_d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			s_d_trajx_0[t]<-xnew
			s_d_trajy_0[t]<-ynew
	}

	#Continuous system - Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	s_c_trajx_0<-array(0,dim=t0+1)
	s_c_trajy_0<-array(0,dim=t0+1)

	s_c_trajx_0[1]<-s_x0
	s_c_trajy_0[1]<-s_y0

	for (t in 2:(t0+1)){
			x<-s_c_trajx_0[t-1]
			y<-s_c_trajy_0[t-1]
			xnew<-(x+x*s_a*(1-x/k1)-x*s_b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-(ynew+ynew*s_c*x*(1-e*ynew/k2)-s_d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			s_c_trajx_0[t]<-xnew
			s_c_trajy_0[t]<-ynew
	}

	#OSCILLATION
	o_a<-1.1
	o_b<-0.07
	o_c<-0.05
	o_d<-0.7
	o_x0<-35
	o_y0<-9


	#Discrete system - Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	o_d_trajx_0<-array(0,dim=t0+1)
	o_d_trajy_0<-array(0,dim=t0+1)

	o_d_trajx_0[1]<-o_x0
	o_d_trajy_0[1]<-o_y0

	for (t in 2:(t0+1)){
			x<-o_d_trajx_0[t-1]
			y<-o_d_trajy_0[t-1]
			xnew<-round(x+x*o_a*(1-x/k1)-x*o_b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-round(ynew+ynew*o_c*x*(1-e*ynew/k2)-o_d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			o_d_trajx_0[t]<-xnew
			o_d_trajy_0[t]<-ynew
	}

	#Continuous system - Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	o_c_trajx_0<-array(0,dim=t0+1)
	o_c_trajy_0<-array(0,dim=t0+1)

	o_c_trajx_0[1]<-o_x0
	o_c_trajy_0[1]<-o_y0

	for (t in 2:(t0+1)){
			x<-o_c_trajx_0[t-1]
			y<-o_c_trajy_0[t-1]
			xnew<-(x+x*o_a*(1-x/k1)-x*o_b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-(ynew+ynew*o_c*x*(1-e*ynew/k2)-o_d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			o_c_trajx_0[t]<-xnew
			o_c_trajy_0[t]<-ynew
	}

	#EXTINCTION
	e_a<-1.1
	e_b<-0.07
	e_c<-0.08
	e_d<-0.7
	e_x0<-32
	e_y0<-15

	#Discrete system - Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	e_d_trajx_0<-array(0,dim=t0+1)
	e_d_trajy_0<-array(0,dim=t0+1)

	e_d_trajx_0[1]<-e_x0
	e_d_trajy_0[1]<-e_y0

	for (t in 2:(t0+1)){
			x<-e_d_trajx_0[t-1]
			y<-e_d_trajy_0[t-1]
			xnew<-round(x+x*e_a*(1-x/k1)-x*e_b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-round(ynew+ynew*e_c*x*(1-e*ynew/k2)-e_d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			e_d_trajx_0[t]<-xnew
			e_d_trajy_0[t]<-ynew
	}

	#Continuous system - Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	e_c_trajx_0<-array(0,dim=t0+1)
	e_c_trajy_0<-array(0,dim=t0+1)

	e_c_trajx_0[1]<-e_x0
	e_c_trajy_0[1]<-e_y0

	for (t in 2:(t0+1)){
			x<-e_c_trajx_0[t-1]
			y<-e_c_trajy_0[t-1]
			xnew<-(x+x*e_a*(1-x/k1)-x*e_b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-(ynew+ynew*e_c*x*(1-e*ynew/k2)-e_d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			e_c_trajx_0[t]<-xnew
			e_c_trajy_0[t]<-ynew
	}

	png(filename="f_model_comparison.png",width = 15, height = 5, units="cm", res=300)

	par(mfrow=c(1,3))
	par(mar = c(3, 3, 1, 0.5))
    par(mgp=c(2,0.5,0))

	plot(s_c_trajx_0,s_c_trajy_0,type="l", ,xlab="Prey",ylab="Predator",ylim=c(0,50), xlim=c(0,50),col="CornflowerBlue", main="Case 1: Stationary", cex.main=0.75)
	points(s_c_trajx_0[1],s_c_trajy_0[1], pch=16, col="CornflowerBlue")
	points(s_c_trajx_0[t0+1],s_c_trajy_0[t0+1], pch=4, col="CornflowerBlue")
	lines(s_d_trajx_0,s_d_trajy_0,type="l", lty=1, col="black")
	points(s_d_trajx_0[1],s_d_trajy_0[1], pch=16, col="black")
	points(s_d_trajx_0[t0+1],s_d_trajy_0[t0+1], pch=4, col="black")
	legend(x="topleft",c("Continuous","Discrete"),lty=c(1,1), col=c("CornflowerBlue","black"), cex=0.75)

	plot(o_c_trajx_0,o_c_trajy_0,type="l", ,xlab="Prey",ylab="Predator",ylim=c(0,50), xlim=c(0,50),col="CornflowerBlue", main="Case 2: Oscillation", cex.main=0.75)
	points(o_c_trajx_0[1],o_c_trajy_0[1], pch=16, col="CornflowerBlue")
	points(o_c_trajx_0[t0+1],o_c_trajy_0[t0+1], pch=4, col="CornflowerBlue")
	lines(o_d_trajx_0,o_d_trajy_0,type="l", lty=1, col="black")
	points(o_d_trajx_0[1],o_d_trajy_0[1], pch=16, col="black")
	points(o_d_trajx_0[t0+1],o_d_trajy_0[t0+1], pch=4, col="black")
	legend(x="topleft",c("Continuous","Discrete"),lty=c(1,1), col=c("CornflowerBlue","black"), cex=0.75)

	plot(e_c_trajx_0,e_c_trajy_0,type="l", ,xlab="Prey",ylab="Predator",ylim=c(0,50), xlim=c(0,50),col="CornflowerBlue", main="Case 3: Extinction", cex.main=0.75)
	points(e_c_trajx_0[1],e_c_trajy_0[1], pch=16, col="CornflowerBlue")
	points(e_c_trajx_0[t0+1],e_c_trajy_0[t0+1], pch=4, col="CornflowerBlue")
	lines(e_d_trajx_0,e_d_trajy_0,type="l", lty=1, col="black")
	points(e_d_trajx_0[1],e_d_trajy_0[1], pch=16, col="black")
	points(e_d_trajx_0[t0+1],e_d_trajy_0[t0+1], pch=4, col="black")
	legend(x="topleft",c("Continuous","Discrete"),lty=c(1,1), col=c("CornflowerBlue","black"), cex=0.75)

	dev.off()
}

results_comparison_oscillation <- function()
{
	#Common and fixed parameters
	t0 <- 100
	k1 <- 50
	k2 <- 50
	e <- 1
	time<-1:(t0+1)

	#Data
	ds_data <- read.csv("oscillatory_trajectories.csv")	#Discrete system (all data)
	cs_so_data <- read.csv("scenario2_1_01.csv", sep=";", dec=",")	#Continuous system (social optimum)
	cs_ne_data <- read.csv("game2_1_01.csv", sep=";", dec=",")	#Continuous system (Nash equilibirum)

	png(filename="f_oscillation_comparison.png",width = 15, height = 10, units="cm", res=300)
	par(mfrow=c(2,3))
	par(mar = c(3, 3, 1, 0.5))
    par(mgp=c(2,0.5,0))

	#Plot: population dynamics under social optimum harvest - discrete system
	plot(time,ds_data$x_so,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Extensive search", col="darkgreen", cex.main=0.75)
	lines(time,ds_data$y_so, col="firebrick", lty=1)
	lines(time,ds_data$y_h_so, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - continuous system
	plot(time,cs_so_data$xt,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Approximate method", col="darkgreen", cex.main=0.75)
	lines(time,cs_so_data$yt, col="firebrick", lty=1)
	lines(time,cs_so_data$at, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Removing first and last N iterations from the dataset and calculate mean harvest for given period
	N<-10
	discrete_h_so <- ds_data$y_h_so
	discrete_h_so <- tail(discrete_h_so,length(discrete_h_so)-N)
	discrete_h_so <- head(discrete_h_so,length(discrete_h_so)-N)

	continuous_h_so <- cs_so_data$at
	continuous_h_so <- tail(continuous_h_so,length(continuous_h_so)-N)
	continuous_h_so <- head(continuous_h_so,length(continuous_h_so)-N)

	mean_discrete_h <- rep(mean(discrete_h_so), (t0+1))
	mean_continuous_h <- rep(mean(continuous_h_so), (t0+1))

	for (i in 1:N)
	{
		mean_discrete_h[i] <- NA
		mean_discrete_h[length(mean_discrete_h)+1-i] <- NA

		mean_continuous_h[i] <- NA
		mean_continuous_h[length(mean_continuous_h)+1-i] <- NA

	}

	#Plot: harvest and average harvest comparison
	plot(time,mean_discrete_h,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Mean harvest", col="black", cex.main=0.75)
	#lines(time,mean_discrete_h, col="black", lty=3)
	#lines(time,cs_so_data$at, col="CornflowerBlue", lty=1)
	lines(time,mean_continuous_h, col="CornflowerBlue", lty=1)
	legend(x="topleft",c("Extensive search","Approximate method"),lty=c(1,1), col=c("black","CornflowerBlue"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - discrete system
	plot(time,ds_data$x_ne,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium - Extensive search", col="darkgreen", cex.main=0.75)
	lines(time,ds_data$y_ne, col="firebrick", lty=1)
	lines(time,ds_data$y_h_ne, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - continuous system
	plot(time,cs_ne_data$xt,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium - Approximate method", col="darkgreen", cex.main=0.75)
	lines(time,cs_ne_data$yt, col="firebrick", lty=1)
	lines(time,cs_ne_data$at, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Removing first and last N iterations from the dataset and calculate mean harvest for given period
	N<-10
	discrete_h_ne <- ds_data$y_h_ne
	discrete_h_ne <- tail(discrete_h_ne,length(discrete_h_ne)-N)
	discrete_h_ne <- head(discrete_h_ne,length(discrete_h_ne)-N)

	continuous_h_ne <- cs_ne_data$at
	continuous_h_ne <- tail(continuous_h_ne,length(continuous_h_ne)-N)
	continuous_h_ne <- head(continuous_h_ne,length(continuous_h_ne)-N)

	mean_discrete_h <- rep(mean(discrete_h_ne), (t0+1))
	mean_continuous_h <- rep(mean(continuous_h_ne), (t0+1))

	for (i in 1:N)
	{
		mean_discrete_h[i] <- NA
		mean_discrete_h[length(mean_discrete_h)+1-i] <- NA

		mean_continuous_h[i] <- NA
		mean_continuous_h[length(mean_continuous_h)+1-i] <- NA

	}

	#Plot: harvest and average harvest comparison
	plot(time,mean_discrete_h,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash eqilibrium - Mean harvest", col="black", cex.main=0.75)
	#lines(time,mean_discrete_h, col="black", lty=3)
	#lines(time,cs_so_data$at, col="CornflowerBlue", lty=1)
	lines(time,mean_continuous_h, col="CornflowerBlue", lty=1)
	legend(x="topleft",c("Extensive search","Approximate method"),lty=c(1,1), col=c("black","CornflowerBlue"),cex=0.5)

	dev.off()

	#Creating csv value with mean data for both algorithms
	data_export <- cbind(mean(ds_data$x_so),mean(ds_data$y_so),mean(ds_data$y_h_so),mean(discrete_h_so),mean(cs_so_data$xt),mean(cs_so_data$yt),mean(cs_so_data$at),mean(continuous_h_so),mean(ds_data$x_ne),mean(ds_data$y_ne),mean(ds_data$y_h_ne),mean(discrete_h_ne),mean(cs_ne_data$xt),mean(cs_ne_data$yt),mean(cs_ne_data$at),mean(continuous_h_ne))
	colnames(data_export) <- c("d_mean_x_so","d_mean_y_so","d_mean_h_so","d_mean_h_so_mid","c_mean_x_so","c_mean_y_so","c_mean_h_so","c_mean_h_so_mid","d_mean_x_ne","d_mean_y_ne","d_mean_h_ne","d_mean_h_ne_mid","c_mean_x_ne","c_mean_y_ne","c_mean_h_ne","c_mean_h_ne_mid")
	data_frame <- as.data.frame(data_export)
	write.csv(data_frame,"oscillation_comparison.csv")
}

results_comparison_stationary <- function()
{
	#Common and fixed parameters
	t0 <- 100
	k1 <- 50
	k2 <- 50
	e <- 1
	time<-1:(t0+1)

	#Data
	ds_data <- read.csv("stationary_trajectories.csv")	#Discrete system (all data)
	cs_so_data <- read.csv("scenario1_1_01.csv", sep=";", dec=",")	#Continuous system (social optimum)
	cs_ne_data <- read.csv("game1_1_01.csv", sep=";", dec=",")	#Continuous system (Nash equilibirum)

	png(filename="f_stationary_comparison.png",width = 15, height = 10, units="cm", res=300)
	par(mfrow=c(2,3))
	par(mar = c(3, 3, 1, 0.5))
    par(mgp=c(2,0.5,0))

	#Plot: population dynamics under social optimum harvest - discrete system
	plot(time,ds_data$x_so,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Extensive search", col="darkgreen", cex.main=0.75)
	lines(time,ds_data$y_so, col="firebrick", lty=1)
	lines(time,ds_data$y_h_so, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - continuous system
	plot(time,cs_so_data$xt,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Approximate method", col="darkgreen", cex.main=0.75)
	lines(time,cs_so_data$yt, col="firebrick", lty=1)
	lines(time,cs_so_data$at, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Removing first and last N iterations from the dataset and calculate mean harvest for given period
	N<-10
	discrete_h_so <- ds_data$y_h_so
	discrete_h_so <- tail(discrete_h_so,length(discrete_h_so)-N)
	discrete_h_so <- head(discrete_h_so,length(discrete_h_so)-N)

	continuous_h_so <- cs_so_data$at
	continuous_h_so <- tail(continuous_h_so,length(continuous_h_so)-N)
	continuous_h_so <- head(continuous_h_so,length(continuous_h_so)-N)

	mean_discrete_h <- rep(mean(discrete_h_so), (t0+1))
	mean_continuous_h <- rep(mean(continuous_h_so), (t0+1))

	for (i in 1:N)
	{
		mean_discrete_h[i] <- NA
		mean_discrete_h[length(mean_discrete_h)+1-i] <- NA

		mean_continuous_h[i] <- NA
		mean_continuous_h[length(mean_continuous_h)+1-i] <- NA

	}

	#Plot: harvest and average harvest comparison
	plot(time,mean_discrete_h,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Mean harvest", col="black", cex.main=0.75)
	lines(time,mean_continuous_h, col="CornflowerBlue", lty=1)
	legend(x="topleft",c("Extensive search","Approximate method"),lty=c(1,1), col=c("black","CornflowerBlue"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - discrete system
	plot(time,ds_data$x_ne,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium - Extensive search", col="darkgreen", cex.main=0.75)
	lines(time,ds_data$y_ne, col="firebrick", lty=1)
	lines(time,ds_data$y_h_ne, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - continuous system
	plot(time,cs_ne_data$xt,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium - Approximate method", col="darkgreen", cex.main=0.75)
	lines(time,cs_ne_data$yt, col="firebrick", lty=1)
	lines(time,cs_ne_data$at, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Removing first and last N iterations from the dataset and calculate mean harvest for given period
	N<-10
	discrete_h_ne <- ds_data$y_h_ne
	discrete_h_ne <- tail(discrete_h_ne,length(discrete_h_ne)-N)
	discrete_h_ne <- head(discrete_h_ne,length(discrete_h_ne)-N)

	continuous_h_ne <- cs_ne_data$at
	continuous_h_ne <- tail(continuous_h_ne,length(continuous_h_ne)-N)
	continuous_h_ne <- head(continuous_h_ne,length(continuous_h_ne)-N)

	mean_discrete_h <- rep(mean(discrete_h_ne), (t0+1))
	mean_continuous_h <- rep(mean(continuous_h_ne), (t0+1))

	for (i in 1:N)
	{
		mean_discrete_h[i] <- NA
		mean_discrete_h[length(mean_discrete_h)+1-i] <- NA

		mean_continuous_h[i] <- NA
		mean_continuous_h[length(mean_continuous_h)+1-i] <- NA

	}

	#Plot: harvest and average harvest comparison
	plot(time,mean_discrete_h,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash eqilibrium - Mean harvest", col="black", cex.main=0.75)
	lines(time,mean_continuous_h, col="CornflowerBlue", lty=1)
	legend(x="topleft",c("Extensive search","Approximate method"),lty=c(1,1), col=c("black","CornflowerBlue"),cex=0.5)

	dev.off()

	#Creating csv value with mean data for both algorithms
	data_export <- cbind(mean(ds_data$x_so),mean(ds_data$y_so),mean(ds_data$y_h_so),mean(discrete_h_so),mean(cs_so_data$xt),mean(cs_so_data$yt),mean(cs_so_data$at),mean(continuous_h_so),mean(ds_data$x_ne),mean(ds_data$y_ne),mean(ds_data$y_h_ne),mean(discrete_h_ne),mean(cs_ne_data$xt),mean(cs_ne_data$yt),mean(cs_ne_data$at),mean(continuous_h_ne))
	colnames(data_export) <- c("d_mean_x_so","d_mean_y_so","d_mean_h_so","d_mean_h_so_mid","c_mean_x_so","c_mean_y_so","c_mean_h_so","c_mean_h_so_mid","d_mean_x_ne","d_mean_y_ne","d_mean_h_ne","d_mean_h_ne_mid","c_mean_x_ne","c_mean_y_ne","c_mean_h_ne","c_mean_h_ne_mid")
	data_frame <- as.data.frame(data_export)
	write.csv(data_frame,"stationary_comparison.csv")
}


results_comparison_extinction <- function()
{
	#Common and fixed parameters
	t0 <- 100
	k1 <- 50
	k2 <- 50
	e <- 1
	time<-1:(t0+1)

	#Data
	ds_data <- read.csv("extinction_trajectories.csv")	#Discrete system (all data)
	cs_so_data <- read.csv("scenario3_1_01.csv", sep=";", dec=",")	#Continuous system (social optimum)
	cs_ne_data <- read.csv("game3_1_01.csv", sep=";", dec=",")	#Continuous system (Nash equilibirum)

	png(filename="f_extinction_comparison.png",width = 15, height = 10, units="cm", res=300)
	par(mfrow=c(2,3))
	par(mar = c(3, 3, 1, 0.5))
    par(mgp=c(2,0.5,0))

	#Plot: population dynamics under social optimum harvest - discrete system
	plot(time,ds_data$x_so,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Extensive search", col="darkgreen", cex.main=0.75)
	lines(time,ds_data$y_so, col="firebrick", lty=1)
	lines(time,ds_data$y_h_so, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - continuous system
	plot(time,cs_so_data$xt,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Approximate method", col="darkgreen", cex.main=0.75)
	lines(time,cs_so_data$yt, col="firebrick", lty=1)
	lines(time,cs_so_data$at, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Removing first and last N iterations from the dataset and calculate mean harvest for given period
	N<-10
	discrete_h_so <- ds_data$y_h_so
	discrete_h_so <- tail(discrete_h_so,length(discrete_h_so)-N)
	discrete_h_so <- head(discrete_h_so,length(discrete_h_so)-N)

	continuous_h_so <- cs_so_data$at
	continuous_h_so <- tail(continuous_h_so,length(continuous_h_so)-N)
	continuous_h_so <- head(continuous_h_so,length(continuous_h_so)-N)

	mean_discrete_h <- rep(mean(discrete_h_so), (t0+1))
	mean_continuous_h <- rep(mean(continuous_h_so), (t0+1))

	for (i in 1:N)
	{
		mean_discrete_h[i] <- NA
		mean_discrete_h[length(mean_discrete_h)+1-i] <- NA

		mean_continuous_h[i] <- NA
		mean_continuous_h[length(mean_continuous_h)+1-i] <- NA

	}

	#Plot: harvest and average harvest comparison
	plot(time,mean_discrete_h,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum - Mean harvest", col="black", cex.main=0.75)
	#lines(time,mean_discrete_h, col="black", lty=3)
	#lines(time,cs_so_data$at, col="CornflowerBlue", lty=1)
	lines(time,mean_continuous_h, col="CornflowerBlue", lty=1)
	legend(x="topleft",c("Extensive search","Approximate method"),lty=c(1,1), col=c("black","CornflowerBlue"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - discrete system
	plot(time,ds_data$x_ne,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium - Extensive search", col="darkgreen", cex.main=0.75)
	lines(time,ds_data$y_ne, col="firebrick", lty=1)
	lines(time,ds_data$y_h_ne, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Plot: population dynamics under social optimum harvest - continuous system
	plot(time,cs_ne_data$xt,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium - Approximate method", col="darkgreen", cex.main=0.75)
	lines(time,cs_ne_data$yt, col="firebrick", lty=1)
	lines(time,cs_ne_data$at, col="black", lty=1, lwd = 0.5)
	legend(x="topleft",c("Prey","Predator","Harvest"),lty=c(1,1,1), lwd=c(1,1,0.5), col=c("darkgreen","firebrick","black"),cex=0.5)

	#Removing first and last N iterations from the dataset and calculate mean harvest for given period
	N<-10
	discrete_h_ne <- ds_data$y_h_ne
	discrete_h_ne <- tail(discrete_h_ne,length(discrete_h_ne)-N)
	discrete_h_ne <- head(discrete_h_ne,length(discrete_h_ne)-N)

	continuous_h_ne <- cs_ne_data$at
	continuous_h_ne <- tail(continuous_h_ne,length(continuous_h_ne)-N)
	continuous_h_ne <- head(continuous_h_ne,length(continuous_h_ne)-N)

	mean_discrete_h <- rep(mean(discrete_h_ne), (t0+1))
	mean_continuous_h <- rep(mean(continuous_h_ne), (t0+1))

	for (i in 1:N)
	{
		mean_discrete_h[i] <- NA
		mean_discrete_h[length(mean_discrete_h)+1-i] <- NA

		mean_continuous_h[i] <- NA
		mean_continuous_h[length(mean_continuous_h)+1-i] <- NA

	}

	#Plot: harvest and average harvest comparison
	plot(time,mean_discrete_h,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash eqilibrium - Mean harvest", col="black", cex.main=0.75)
	lines(time,mean_continuous_h, col="CornflowerBlue", lty=1)
	legend(x="topleft",c("Extensive search","Approximate method"),lty=c(1,1), col=c("black","CornflowerBlue"),cex=0.5)

	dev.off()

	#Creating csv value with mean data for both algorithms
	data_export <- cbind(mean(ds_data$x_so),mean(ds_data$y_so),mean(ds_data$y_h_so),mean(discrete_h_so),mean(cs_so_data$xt),mean(cs_so_data$yt),mean(cs_so_data$at),mean(continuous_h_so),mean(ds_data$x_ne),mean(ds_data$y_ne),mean(ds_data$y_h_ne),mean(discrete_h_ne),mean(cs_ne_data$xt),mean(cs_ne_data$yt),mean(cs_ne_data$at),mean(continuous_h_ne))
	colnames(data_export) <- c("d_mean_x_so","d_mean_y_so","d_mean_h_so","d_mean_h_so_mid","c_mean_x_so","c_mean_y_so","c_mean_h_so","c_mean_h_so_mid","d_mean_x_ne","d_mean_y_ne","d_mean_h_ne","d_mean_h_ne_mid","c_mean_x_ne","c_mean_y_ne","c_mean_h_ne","c_mean_h_ne_mid")
	data_frame <- as.data.frame(data_export)
	write.csv(data_frame,"ISDG2022_extinction_comparison.csv")
}

