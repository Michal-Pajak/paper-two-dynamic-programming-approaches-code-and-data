## Analysis scritps

predator_prey_full_analysis <- function(fname,parameters,line,no_processes)
{
	#Default script uses functions designed for parallel computations
	
	#Arguments
	#	parameters - dataframe with values of parameters
	#	line - position of parameter variables in the dataframe vectors
	#	no_processes - number of parallel processes run by the algorithms

	#Output
	#	set of files with outcomes of the analysis - for more information check the description of save_output_social_optimum_and_nash_equilibrium function

	output_so <- predator_prey_analysis_social_optimum_parallel(parameters,line,no_processes)
	output_ne <- predator_prey_analysis_nash_equilibrium_parallel(parameters,line,no_processes)
	save_output_social_optimum_and_nash_equilibrium(fname,output_so,output_ne,parameters,line)
}

predator_prey_analysis_social_optimum_parallel <-  function(parameters,line,no_processes)
{
	#Function finds the social optimum trajectory in the system using parallel computations

	#Arguments
	#	parameters - dataframe with values of parameters
	#	line - position of parameter variables in the dataframe vectors
	#	no_processes - number of parallel processes run by the algorithms

	#Output
	#	result - dataframe with trajectories (no harvest trajectory of the system and social optimum trajectory)

	#Setting up class used to save data from parallel computations
	resultClass <- function()
	{
		#Class stores data from parallel execution of the foreach function
		me<-list(
				rc_v_opt = array(0,dim=c(k2+1)),
				rc_p_opt = array(0,dim=c(k2+1)),
				rc_h_opt = array(list(c(0,0)),dim=c(k2+1)),
				rc_s_opt = array(list(c(0,0)),dim=c(k2+1)),
				rc_next_xy_opt = array(list(c(0,0)),dim=c(k2+1))
		)

		class(me) <- append(class(me),"resultClass")
		return(me)
	}

	#Setting up parameters based on the arguments
	k1<-parameters$k1[line]
	k2<-parameters$k2[line]
	max_x_s1 <-parameters$max_x_s1[line]
	max_x_s2 <-parameters$max_x_s2[line]
	max_y_s1 <-parameters$max_y_s1[line]
	max_y_s2 <-parameters$max_y_s2[line]
	HM_x_1 <- parameters$HM_x_1[line]
	HM_x_2 <- parameters$HM_x_2[line]
	HM_y_1 <- parameters$HM_y_1[line]
	HM_y_2 <- parameters$HM_y_2[line]
	PM_x_1 <- parameters$PM_x_1[line]
	PM_x_2 <- parameters$PM_x_2[line]
	PM_y_1 <- parameters$PM_y_1[line]
	PM_y_2 <- parameters$PM_y_2[line]
	HM_x <- round((HM_x_1+HM_x_2)/2)					#Setting up common harvest multiplier for prey population
	HM_y <- round((HM_y_1+HM_y_2)/2)					#Setting up common harvest multiplier for predator population
	PM_x <- round((PM_x_1+PM_x_2)/2)					#Setting up common payoff multiplier for prey population
	PM_y <- round((PM_y_1+PM_y_2)/2)					#Setting up common payoff multiplier for predator population
	max_x_ts <- min(max_x_s1+max_x_s2,round(k1/HM_x))	#Setting up common maximum decision for prey population - scaling decisions to the maximum possible harvest
	max_y_ts <- min(max_y_s1+max_y_s2,round(k2/HM_y))	#Setting up common maximum decision for predator population - scaling decisions to the maximum possible harvest
	##TODO TEST Scaling test when the values are not symmetric
	##TODO TEST Scaling test when the strategies are higher than the size of the resource
	a<-parameters$a[line]
	b<-parameters$b[line]
	c<-parameters$c[line]
	d<-parameters$d[line]
	e<-parameters$e[line]
	gamma<-parameters$gamma[line]
	t0<-parameters$t0[line]
	x0<-parameters$x0[line]
	y0<-parameters$y0[line]
	symmetric<-parameters$sym[line]

	#Setting up parallel procsses count
	cl <- makeCluster(no_processes)
	registerDoParallel(cl)

	#Setting up values for analysis and data structures for storing the results of calculations 
	time<-1:(t0+1)
	v_opt<-array(0,dim=c(k1+1,k2+1,t0+1))						#Value for optimal the strategies in a given state in given iteration
	p_opt<-array(0,dim=c(k1+1,k2+1,t0+1))						#Myopic payoff from the optimal strategies in a given state in given iteration
	h_opt<-array(list(c(-1,-1)),dim=c(k1+1,k2+1,t0+1))			#Set of optimal harvest in a given state in a given iteration
	s_opt<-array(list(c(-1,-1)),dim=c(k1+1,k2+1,t0+1))			#Set of optimal harvest strategies in a given state in given iteration in form (prey,predator)
	next_xy_opt <- array(list(c(-1,-1)),dim=c(k1+1,k2+1,t0+1))	#Next state for a given state in given iteration under social optimum strategy

	#Setting up data structures for trajectory - single path generated at the end of the analysis
	trajsx<-array(0,dim=c(t0+1))		#Optimal strategy for prey population
	trajsy<-array(0,dim=c(t0+1))		#Optimal strategy for predator population
	trajhx<-array(0,dim=c(t0+1))		#Optimal harvest for prey population
	trajhy<-array(0,dim=c(t0+1))		#Optimal harvest for predator population
	trajp<-array(0,dim=c(t0+1))			#Myopic payoff from the optimal strategies
	trajap<-array(0,dim=c(t0+1))		#Accumulated payoff from the optimal strategies
	trajx<-array(0,dim=c(t0+1))			#Size of the prey population under optimal control
	trajy<-array(0,dim=c(t0+1))			#Size of the predator population under optimal control
	trajv<-array(0,dim=c(t0+1))			#Value for the trajectory under optimal control - saved for testing purposes

	cat("\n")
	cat("Social optimum progress:\n")
	cat(paste0(round(1/(t0+2)*100),"%\n"))

	#For the final timestep of the analysis iterate over each possible state of the system for both population in parallel
	result_par <- foreach (i=1:(k1+1)) %dopar%
	{
		#Setting up the result class
		result_class <- resultClass()

		for (j in 1:(k2+1))
		{
			#Setting up placeholders for the payoff, harvest and strategy values
			temp_payoff <- -1
			test_payoff <- -1
			temp_h <- -1
			temp_s <- c(-1,-1)
			temp_next_xy <- c(-1,-1)

			#Adjusting maximum strategy/harvesting effort to match the size of the prey population
			temp_x_max_s <- min(max_x_ts,round((i-1)/HM_x))

			#Iterating over all possible harvest strategies in prey population
			for (sx in 1:(temp_x_max_s+1))
			{
				#Adjusting maximum strategy/harvesting effort to match the size of the predator population
				temp_y_max_s <- min(max_y_ts,round((j-1)/HM_y))

				#Iterating over all possible harvest strategies in predator population
				for (sy in 1:(temp_y_max_s+1))
				{
					#If the generated test_payoff is higher than the saved temp_payoff the new solution is saved as being the optimal one
					test_payoff <- (sx-1)*HM_x*PM_x + (sy-1)*HM_y*PM_y
					if (test_payoff > temp_payoff)
					{
						temp_payoff <- test_payoff
						temp_h <- c((sx-1)*HM_x,(sy-1)*HM_y)
						temp_s <- c(sx-1,sy-1)
						temp_next_xy <- c(0,0)
					}
				}
			}

			#Saving data for the optimal solution - solution with the highest payoff found for the given state of the system
			result_class$rc_v_opt[j]<-temp_payoff
			result_class$rc_p_opt[j]<-temp_payoff
			result_class$rc_h_opt[j]<-list(temp_h)
			result_class$rc_s_opt[j]<-list(temp_s)
			result_class$rc_next_xy_opt[j]<-list(temp_next_xy)
		}

		#Returning class with data
		return(result_class)
	}
	
	#Saving results of parallel computations in the matrices
	for (it_i in 1:length(result_par))
	{
		for (it_j in 1:(k2+1))
		{
			v_opt[it_i,it_j,1]<-result_par[[it_i]]$rc_v_opt[it_j]
			p_opt[it_i,it_j,1]<-result_par[[it_i]]$rc_p_opt[it_j]
			h_opt[it_i,it_j,1]<-result_par[[it_i]]$rc_h_opt[it_j]
			s_opt[it_i,it_j,1]<-result_par[[it_i]]$rc_s_opt[it_j]
			next_xy_opt[it_i,it_j,1]<-result_par[[it_i]]$rc_next_xy_opt[it_j]
		}
	}
	
	if (t0>0)
	{

		# Calculating the optimal trajectory through the backward induction
		for (t in 2:(t0+1))
		{
			cat(paste0(round(t/(t0+2)*100),"%\n"))

			#Iteration over all posible states of the system in parallel
			result_par <- foreach (i=1:(k1+1)) %dopar%
			{
				#Setting up the result class
				result_class <- resultClass()

				for (j in 1:(k2+1))
				{
					#Setting up placeholders for the payoff, harvest and strategy values
					temp_payoff <- -1
					test_payoff <- -1
					temp_h <- -1
					temp_s <- c(-1,-1)
					temp_next_xy <- c(-1,-1)

					#Adjusting maximum strategy/harvesting effort to match the size of the population
					temp_x_max_s <- min(max_x_ts,round((i-1)/HM_x))

					#Iterating over all possible harvesting strategies for prey population
					for (sx in 1:(temp_x_max_s+1))
					{
						#Adjusting maximum strategy/harvesting effort to match the size of the population
						temp_y_max_s <- min(max_y_ts,round((j-1)/HM_y))

						#Iterating over all possible harvesting strategies for predator population
						for (sy in 1:(temp_y_max_s+1))
						{	
							#Calculating system state for the next iteration based on analyzed harvesting strategies
							x_h <- (i-1)-(sx-1)*HM_x
							y_h <- (j-1)-(sy-1)*HM_y
							xnew<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
							xnew<-max(0,xnew)
							xnew<-min(k1,xnew)
							ynew<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
							ynew<-max(0,ynew)
							ynew<-min(k2,ynew)

							#If the generated test_payoff is higher than the saved temp_payoff the new solution is saved as being the optimal one
							test_payoff <- ((sx-1)*HM_x*PM_x+(sy-1)*HM_y*PM_y)+gamma*v_opt[xnew+1,ynew+1,t-1] 
							if (test_payoff > temp_payoff)
							{
								temp_payoff <- test_payoff
								temp_h <- c((sx-1)*HM_x,(sy-1)*HM_y)
								temp_s <- c(sx-1,sy-1)
								temp_next_xy <- c(xnew,ynew)
							}
						}
					}

					#Saving data for the optimal solution - solution with the highest payoff found for the given state of the system
					result_class$rc_v_opt[j]<-temp_payoff
					#result_class$rc_p_opt[j]<-(sx-1)*HM_x*PM_x+(sy-1)*HM_y*PM_y
					result_class$rc_p_opt[j]<-temp_h[1]*PM_x+temp_h[2]*PM_y
					result_class$rc_h_opt[j]<-list(temp_h)
					result_class$rc_s_opt[j]<-list(temp_s)
					result_class$rc_next_xy_opt[j]<-list(temp_next_xy)	
				}

				#Returning class with data
				return(result_class)
			}

			#Saving results of parallel computations in the matrices
			for (it_i in 1:length(result_par))
			{
				for (it_j in 1:(k2+1))
				{
					v_opt[it_i,it_j,t]<-result_par[[it_i]]$rc_v_opt[it_j]
					p_opt[it_i,it_j,t]<-result_par[[it_i]]$rc_p_opt[it_j]
					h_opt[it_i,it_j,t]<-result_par[[it_i]]$rc_h_opt[it_j]
					s_opt[it_i,it_j,t]<-result_par[[it_i]]$rc_s_opt[it_j]
					next_xy_opt[it_i,it_j,t]<-result_par[[it_i]]$rc_next_xy_opt[it_j]
				}
			}
		}
	}
	
	cat(paste0("100%\n"))
	cat("\n")

	#Calculating the optimal trajectory starting from the initial state of the system
	trajx[1]<-x0
	trajy[1]<-y0
	trajsx[1]<-s_opt[[x0+1,y0+1,t0+1]][1]
	trajsy[1]<-s_opt[[x0+1,y0+1,t0+1]][2]
	trajhx[1]<-h_opt[[x0+1,y0+1,t0+1]][1]
	trajhy[1]<-h_opt[[x0+1,y0+1,t0+1]][2]
	trajp[1]<-p_opt[x0+1,y0+1,t0+1]
	trajv[1]<-v_opt[x0+1,y0+1,t0+1]
	trajap[1]<-(trajp[1])

	if(t0>0)
	{
		for (t in 2:(t0+1))
		{
			x<-trajx[t-1]
			y<-trajy[t-1]
			hx<-trajhx[t-1]
			hy<-trajhy[t-1]
			x_h<-x-hx
			y_h<-y-hy
			xnew<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			trajx[t]<-xnew
			trajy[t]<-ynew
			trajsx[t]<-s_opt[[xnew+1,ynew+1,t0+2-t]][1]
			trajsy[t]<-s_opt[[xnew+1,ynew+1,t0+2-t]][2]
			trajhx[t]<-h_opt[[xnew+1,ynew+1,t0+2-t]][1]
			trajhy[t]<-h_opt[[xnew+1,ynew+1,t0+2-t]][2]
			trajp[t]<-p_opt[xnew+1,ynew+1,t0+2-t]*(gamma^(t-1))
			trajv[t]<-v_opt[xnew+1,ynew+1,t0+2-t]
			trajap[t]<-trajp[t]+trajap[t-1]
		}
	}

	#Calculating the trajectory for the system without harvesting
	h<-0							#Harvest set to 0
	trajx_0<-array(0,dim=t0+1)
	trajy_0<-array(0,dim=t0+1)

	trajx_0[1]<-x0
	trajy_0[1]<-y0

	for (t in 2:(t0+1)){
			x<-trajx_0[t-1]
			y<-trajy_0[t-1]
			xnew<-round(x+x*a*(1-x/k1)-x*b*y)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-y-h
			ynew<-round(ynew+ynew*c*x*(1-e*ynew/k2)-d*ynew)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			trajx_0[t]<-xnew
			trajy_0[t]<-ynew
	}

	#Creating and returning dataframe from variables
	data_export <- cbind(time,trajx_0,trajy_0,trajx,trajy,trajhx,trajhy,trajsx,trajsy,trajp,trajap,trajv)
	colnames(data_export) <- c("t","x_h0","y_h0","x_so","y_so","x_h_so","y_h_so","x_s_so","y_s_so","p_so","ap_so","v_so")
	data_frame <- as.data.frame(data_export)

	output <- list(data_frame,next_xy_opt,h_opt)

	return(output)
}

predator_prey_analysis_nash_equilibrium_parallel <- function(parameters, line, parallel_processes)
{
	#Function finds Nash equilibrium trajectory in the system using falsification method for symmetric solutions and matrix method for non-symmetric ones. Calculations are done as parallel processes.
	#	For each set of symmetric strategies analysis checks if it is possible for one player to deviate and obtain better payoff.
	#		If it is not possible solution is saved as a Nash equlibrium
	#	For non-symmetric solutions as given solution is verified the algorithm fills out the matrix with information whether the given solution is maximum/best response for the player given fixed strategy of the other player
	#		Subsequent loops first checks whether the matrix is already filled in for given solution for given player and already calculated data can be used.
	#		If not, algorithm conducts analysis filling out new rows and columns of the multi-dimensional matrix.
	#	If sym parameter is set to TRUE function checks the symmetric sets of strategies first. The solution with the set of strategies with the highest payoff is saved.
	#		If there are no symmetric solutions found the algorithm checks all possible strategy sets.
	#		The solution with the set of strategies with highest payoff and most symmetric payoff is saved.		
	#	If sym parameter is set to FALSE function just analyzes all possible strategy sets. The solution with the set of strategies with highest payoff and most symmetric payoff is saved.

	#NOTE:
	#	If both players decide to harvest more resource than it is available each player gets the payoff that is proportional to his strategy/harvesting effort
	#	Analysis of symmetric solutions is valid only if both harvest and payoff multiplayers for both players are equal

	#Arguments
	#	parameters - dataframe with values of parameters
	#	line - position of parameter variables in the dataframe vectors

	#Output
	#	result - dataframe with trajectory data

	#Setting up class used to save data from parallel computations
	resultClass <- function()
	{
		#Class stores data from parallel execution of the foreach function
		me<-list(
			rc_x_s_ne = array(list(c(0,0)),dim=c(k2+1)),
			rc_s_ne = array(list(c(0,0)),dim=c(k2+1)),
			rc_ne = array(list(c(0,0)),dim=c(k2+1)),
			rc_ne = array(list(c(0,0)),dim=c(k2+1)),
			rc_ne = array(list(c(0,0)),dim=c(k2+1)),
			rc_next_xy_ne = array(list(c(0,0)),dim=c(k2+1))
		)

		class(me) <- append(class(me),"resultClass")
		return(me)
	}

	#Setting up parameters based on the arguments
	k1<-parameters$k1[line]
	k2<-parameters$k2[line]
	max_x_s1 <-parameters$max_x_s1[line]
	max_x_s2 <-parameters$max_x_s2[line]
	max_y_s1 <-parameters$max_y_s1[line]
	max_y_s2 <-parameters$max_y_s2[line]
	max_x_ts <- min(max_x_s1+max_x_s2,k2)
	max_y_ts <- min(max_y_s1+max_y_s2,k2)
	HM_x_1 <- parameters$HM_x_1[line]
	HM_x_2 <- parameters$HM_x_2[line]
	HM_y_1 <- parameters$HM_y_1[line]
	HM_y_2 <- parameters$HM_y_2[line]
	PM_x_1 <- parameters$PM_x_1[line]
	PM_x_2 <- parameters$PM_x_2[line]
	PM_y_1 <- parameters$PM_y_1[line]
	PM_y_2 <- parameters$PM_y_2[line]
	a<-parameters$a[line]
	b<-parameters$b[line]
	c<-parameters$c[line]
	d<-parameters$d[line]
	e<-parameters$e[line]
	gamma<-parameters$gamma[line]
	t0<-parameters$t0[line]
	x0<-parameters$x0[line]
	y0<-parameters$y0[line]
	symmetric<-parameters$sym[line]

	#Setting up parallel procsses count
	cl <- makeCluster(parallel_processes)
	registerDoParallel(cl)

	time<-1:(t0+1)
	trajx_ne<-array(0,dim=c(t0+1))						#Trajectory of the prey population under Nash equilibrium
	trajy_ne<-array(0,dim=c(t0+1))						#Trajectory of the predator population under Nash equilibrium
	traj_x_h_ne<-array(0,dim=c(t0+1))						#Trajectory of harvest under Nash equilibrium
	traj_y_h_ne<-array(0,dim=c(t0+1))	
	traj_x_s1_ne<-array(0,dim=c(t0+1))						#Trajectory of strategy for Player 1 under Nash equilibrium
	traj_x_s2_ne<-array(0,dim=c(t0+1))
	traj_y_s1_ne<-array(0,dim=c(t0+1))						#Trajectory of strategy for Player 1 under Nash equilibrium
	traj_y_s2_ne<-array(0,dim=c(t0+1))							#Trajectory of strategy for Player 2 under Nash equilibrium
	trajp1_ne<-array(0,dim=c(t0+1))						#Trajectory of current payoff for Player 1 under Nash equilibrium
	trajp2_ne<-array(0,dim=c(t0+1))						#Trajectory of current payoff for Player 2 under Nash equilibrium
	trajap1_ne<-array(0,dim=c(t0+1))					#Trajectory of accumulated payoff for Player 1 under Nash equilibrium
	trajap2_ne<-array(0,dim=c(t0+1))					#Trajectory of accumulated payoff for Player 2 under Nash equilibrium
	trajv1_ne<-array(0,dim=c(t0+1))
	trajv2_ne<-array(0,dim=c(t0+1))
	x_s_ne<-array(list(c(0,0)),dim=c(k1+1,k2+1,t0+1))		#Strategies for Nash equilibrium for specific state of the system
	y_s_ne<-array(list(c(0,0)),dim=c(k1+1,k2+1,t0+1))		#Strategies for Nash equilibrium for specific state of the system
	p_ne<-array(list(c(0,0)),dim=c(k1+1,k2+1,t0+1))		#Payoff for Nash equilibrium for specific state of the system
	v_ne<-array(list(c(0,0)),dim=c(k1+1,k2+1,t0+1))	
	h_ne<-array(list(c(0,0)),dim=c(k1+1,k2+1,t0+1))
	next_xy_ne <- array(list(c(-1,-1)),dim=c(k1+1,k2+1,t0+1))	#Next state for a given state in given iteration under Nash equilibrium strategy

	cat("\n")
	cat("Nash equilibrium progress:\n")
	cat(paste0(round(1/(t0+2)*100),"%\n"))

	#Iterate over all possible states of the system in the final iteration
	#	Each process calculate all possible states of predator population for given prey population
	result_par <- foreach (i=1:(k1+1)) %dopar%
	{
		#Setting up the result class
		result_class <- resultClass()

		for (j in 1:(k2+1))
		{
			#Setting up data for the Nash equilibrium
			ne_found <- FALSE
			temp_ne_payoff <- c(-1,-1)
	
			#Setting up maximum harvesting range for i/x and j/y
			temp_max_x_s1 <- min(max_x_s1,round((i-1)/HM_x_1))
			temp_max_x_s2 <- min(max_x_s2,round((i-1)/HM_x_2))
			temp_max_y_s1 <- min(max_y_s1,round((j-1)/HM_y_1))
			temp_max_y_s2 <- min(max_y_s2,round((j-1)/HM_y_2))
			temp_next_xy <- c(-1,-1)

			#Analysis of symmetric strategies - falsification method
			if(symmetric==TRUE)
			{
				#Setting up maximum symmetric strategies and symmetric value for harvest multipliers
				temp_max_x_sym <- min(temp_max_x_s1,temp_max_x_s2)
				temp_max_y_sym <- min(temp_max_y_s1,temp_max_y_s2)
				HM_x <- round((HM_x_1+HM_x_2)/2)
				HM_y <- round((HM_y_1+HM_y_2)/2)

				#Iteration over all possible sests of symmetric strategies
				for (x_sym in 1:(temp_max_x_sym+1))
				{
					for (y_sym in 1:(temp_max_y_sym+1))
					{
						#Setting up temporary payoff value for P1 and P2 and their harvest
						temp_payoff <- c(-1,-1)
						
						#Calculating temporary payoff
						if ((2*(x_sym-1)*HM_x)<=(i-1))
						{
							temp_payoff <- c((x_sym-1)*HM_x*PM_x_1,(x_sym-1)*HM_x*PM_x_2)
						}
						else
						{
							temp_payoff <- c((x_sym-1)*(i-1)*HM_x/(2*(x_sym-1)*HM_x)*PM_x_1,(x_sym-1)*(i-1)*HM_x/(2*(x_sym-1)*HM_x)*PM_x_2)
						}

						if ((2*(y_sym-1)*HM_y)<=(j-1))
						{
							temp_payoff <- temp_payoff + c((y_sym-1)*HM_y*PM_y_1,(y_sym-1)*HM_y*PM_y_2)
						}
						else
						{
							temp_payoff <- temp_payoff + c((y_sym-1)*(j-1)*HM_y/(2*(y_sym-1)*HM_y)*PM_y_1,(y_sym-1)*(j-1)*HM_y/(2*(y_sym-1)*HM_y)*PM_y_2)
						}
								
						#Setting initial value for the temp_ne_payoff
						if (x_sym==1 & y_sym==1)
						{
							temp_ne_payoff<-temp_payoff
							temp_x_ne_s <- c(x_sym-1,x_sym-1)
							temp_y_ne_s <- c(y_sym-1,y_sym-1)
						}

						#Verification of the solution for P1
						real_ne <- TRUE		#initial assumption that the solution is Nash equilibrium

						#Iteration over all possible strategies for P1
						t_xs1_max <- min(max_x_s1,round((i-1)/HM_x_1))
								
						for (t_xs1 in 1:(t_xs1_max+1))
						{
							if(real_ne==FALSE){break}		#Exit the loop if solution is not a Nash equlibrium

							#t_ys1_max <- min(max_y_s1,j-1)
							t_ys1_max <- min(max_y_s1,round((j-1)/HM_y_1))

							for (t_ys1 in 1:(t_ys1_max+1))
							{
								#Setting up verification payoff
								if (((t_xs1-1)*HM_x_1+(x_sym-1)*HM_x)<=(i-1))
								{
									t_p <- (t_xs1-1)*HM_x_1*PM_x_1
								}
								else
								{
									t_p <- (t_xs1-1)*(i-1)*HM_x_1/((t_xs1-1)*HM_x_1+(x_sym-1)*HM_x)*PM_x_1
								}

								if (((t_ys1-1)*HM_y_1+(y_sym-1)*HM_y)<=(j-1))
								{
									t_p <- t_p + (t_ys1-1)*HM_y_1*PM_y_1
								}
								else
								{
									t_p <- t_p + (t_ys1-1)*(j-1)*HM_y_1/((t_ys1-1)*HM_y_1+(y_sym-1)*HM_y)*PM_y_1
								}

								#If verification payoff is higher than the initial payoff the solution is not a Nash equilibrium
								if (t_p>temp_payoff[1])
								{
									real_ne <- FALSE
									break
								}
							}
						}

						#Verification of the solution for the P2 is required only if the solution was not verified negatively
						if (real_ne == TRUE)
						{
							#Iteration over all possible strategies for P2
							t_xs2_max <- min(max_x_s2,round((i-1)/HM_x_2))
									
							for (t_xs2 in 1:(t_xs2_max+1))
							{
								if(real_ne==FALSE){break}		#Exit the loop if solution is not a Nash equlibrium

								t_ys2_max <- min(max_y_s2,round((j-1)/HM_y_2))

								for (t_ys2 in 1:(t_ys2_max+1))
								{
									#Caluclating verification payoff
									if (((x_sym-1)*HM_x+(t_xs2-1)*HM_x_2)<=(i-1))
									{
										t_p <- (t_xs2-1)*HM_x_2*PM_x_2
									}
									else
									{
										t_p <- (t_xs2-1)*(i-1)*HM_x_2/((x_sym-1)*HM_x+(t_xs2-1)*HM_x_2)*PM_x_2
									}

									if (((y_sym-1)*HM_y_2+(t_ys2-1)*HM_y)<=(j-1))
									{
										t_p <- t_p + (t_ys2-1)*HM_y_2*PM_y_2
									}
									else
									{
										t_p <- t_p + (t_ys2-1)*(j-1)*HM_y_2/((y_sym-1)*HM_y+(t_ys2-1)*HM_y_2)*PM_y_2
									}

									#If verification payoff is higher than the initial payoff the solution is not a Nash equilibrium
									if (t_p>temp_payoff[2])
									{
										real_ne <- FALSE
										break
									}
								}
							}
						}

						#If set of strategies generates Nash equilibrium it is compared against already found Nash equilibrium
						if (real_ne==TRUE)
						{
							ne_found <- TRUE

							#If the group payoff of analyzed Nash equilibrium yelds higher payoff than the one already found then new Nash equilibrium is saved
							if ((temp_payoff[1]+temp_payoff[2])>(temp_ne_payoff[1]+temp_ne_payoff[2]))
							{
								temp_ne_payoff<-temp_payoff
								temp_x_ne_s <- c(x_sym-1,x_sym-1)
								temp_y_ne_s <- c(y_sym-1,y_sym-1)
							}	
						}
					}
				}
			}
			
			#If there is no symmetric Nash equilibirum algorithm goes through all possible strategy sets - matrix method
			if(ne_found == FALSE)
			{
				#Best response/Nash equilibria matrix containing information about whether the combination of the strategies is the Nash equilibrium, size is based on the maximum possible harvesting effort
				#Information:
				#	-1 not checked
				#	 0 checked/not best response for specific player
				#	 1 checked/best response for specific player
				ne_matrix<-array(list(c(-1,-1)),dim=c(temp_max_y_s1+1,temp_max_x_s1+1,temp_max_y_s2+1,temp_max_x_s2+1))

				#Iterate over all possible strategies for i/x
				for (x_s1 in 1:(temp_max_x_s1+1))
				{
					for (x_s2 in 1:(temp_max_x_s2+1))
					{
						#Iterate over all possible strategies for j/y
						for (y_s1 in 1:(temp_max_y_s1+1))
						{
							for (y_s2 in 1:(temp_max_y_s2+1))
							{					
								#Setting up temporary payoff value for P1 and P2
								temp_payoff <- c(-1,-1)

								#Calculating temporary payoff
								if (((x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2)<=(i-1))
								{
									temp_payoff <- c((x_s1-1)*HM_x_1*PM_x_1,(x_s2-1)*HM_x_2*PM_x_2)
								}
								else
								{
									temp_payoff <- c((x_s1-1)*(i-1)*HM_x_1/((x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2)*PM_x_1,(x_s2-1)*(i-1)*HM_x_2/((x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2)*PM_x_2)
								}

								if (((y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2)<=(j-1))
								{
									temp_payoff <- temp_payoff + c((y_s1-1)*HM_y_1*PM_y_1,(y_s2-1)*HM_y_2*PM_y_2)
								}
								else
								{
									temp_payoff <- temp_payoff + c((y_s1-1)*(j-1)*HM_y_1/((y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2)*PM_y_1,(y_s2-1)*(j-1)*HM_y_2/((y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2)*PM_y_2)
								}

								#Setting initial value for the temp_ne_payoff
								if (x_s1==1 & x_s2==1 & y_s1==1 & y_s2==1)
								{
									temp_ne_payoff<-temp_payoff
									temp_x_ne_s <- c(x_s1-1,x_s2-1)
									temp_y_ne_s <- c(y_s1-1,y_s2-1)
								}

								#Verification for the player is required only 
								if (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][1]==-1)
								{							
									#Setting up the maximum payoff updated over the verification
									current_max_payoff <- temp_payoff[1]

									#Best response vectors, storing indices of the the best responses
									x_br<-vector()
									y_br<-vector()

									#Iteration over all possible strategies for P1
									t_xs1_max <- min(max_x_s1,i-1)
									
									for (t_xs1 in 1:(t_xs1_max+1))
									{
										t_ys1_max <- min(max_y_s1,j-1)

										for (t_ys1 in 1:(t_ys1_max+1))
										{
											#Calculating verification payoff
											if (((t_xs1-1)*HM_x_1+(x_s2-1)*HM_x_2)<=(i-1))
											{
												t_p <- (t_xs1-1)*HM_x_1*PM_x_1
											}
											else
											{
												t_p <- (t_xs1-1)*(i-1)*HM_x_1/((t_xs1-1)*HM_x_1+(x_s2-1)*HM_x_2)*PM_x_1
											}

											if (((t_ys1-1)*HM_y_1+(y_s2-1)*HM_y_2)<=(j-1))
											{
												t_p <- t_p + (t_ys1-1)*HM_y_1*PM_y_1
											}
											else
											{
												t_p <- t_p + (t_ys1-1)*(j-1)*HM_y_1/((t_ys1-1)*HM_y_1+(y_s2-1)*HM_y_2)*PM_y_1
											}

											#If verification payoff is higher than the current max payoff, the max payoff is updated and index of the solution is stored
											if (t_p>current_max_payoff)
											{
												current_max_payoff<-t_p
												x_br<-t_xs1
												y_br<-t_ys1
											}
											else
											{
												#If verification payoff is equal to the current max payoff, the index of the solution is added to the list
												if (t_p==current_max_payoff)
												{
													x_br<-c(x_br,t_xs1)
													y_br<-c(y_br,t_ys1)
												}
											}

											#Information that the specific solution was verified is added to the best response matrix
											ne_matrix[[t_ys1,t_xs1,y_s2,x_s2]][1]<-0
										}
									}

									#Setting up information about the best responses into best response matrix
									for (it in 1:length(x_br))
									{
										ne_matrix[[y_br[it],x_br[it],y_s2,x_s2]][1]<-1
									}
								}

								#Verification of the solution for the P2 is required only if the solution was not verified negatively
								if (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][1]==1)
								{
									#Verification is required only if the solution was not already analyzed
									if (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][2]==-1)
									{
										#Setting up the maximum payoff
										current_max_payoff<-temp_payoff[2]

										#Best response vectors, storing indexes of the the best responses
										x_br<-vector()
										y_br<-vector()

										#Iteration over all possible strategies for P2
										t_xs2_max <- min(max_x_s2,i-1)
										
										for (t_xs2 in 1:(t_xs2_max+1))
										{
											t_ys2_max <- min(max_y_s2,j-1)

											for (t_ys2 in 1:(t_ys2_max+1))
											{
												#Calculating verification payoff
												if (((x_s1-1)*HM_x_1+(t_xs2-1)*HM_x_2)<=(i-1))
												{
													t_p <- (t_xs2-1)*HM_x_2*PM_x_2
												}
												else
												{
													t_p <- (t_xs2-1)*(i-1)*HM_x_2/((x_s1-1)*HM_x_1+(t_xs2-1)*HM_x_2)*PM_x_2
												}

												if (((y_s1-1)*HM_y_1+(t_ys2-1)*HM_y_2)<=(j-1))
												{
													t_p <- t_p + (t_ys2-1)*HM_y_2*PM_y_2
												}
												else
												{
													t_p <- t_p + (t_ys2-1)*(j-1)*HM_y_2/((y_s1-1)*HM_y_1+(t_ys2-1)*HM_y_2)*PM_y_2
												}

												#If verification payoff is higher than the current max payoff, the max payoff is updated and index of the solution is stored
												if (t_p>current_max_payoff)
												{
													current_max_payoff<-t_p
													x_br<-t_xs2
													y_br<-t_ys2
												}
												else
												{
													#If verification payoff is equal to the current max payoff, the index of the solution is added to the list
													if (t_p==current_max_payoff)
													{
														x_br<-c(x_br,t_xs2)
														y_br<-c(y_br,t_ys2)
													}
												}

												#Information that the specific solution was verified is added to the best response matrix
												ne_matrix[[y_s1,x_s1,t_ys2,t_xs2]][2]<-0
											}
										}

										#Setting up information about the best responses into best response matrix
										for (it in 1:length(x_br))
										{
											ne_matrix[[y_s1,x_s1,y_br[it],x_br[it]]][2]<-1
										}
									}
								}

								#If set of strategies generates Nash equilibrium it is compared against already found Nash equilibrium
								if ((ne_matrix[[y_s1,x_s1,y_s2,x_s2]][1]==1) & (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][2]==1))
								{
									ne_found <- TRUE

									#If the group payoff of analyzed Nash equilibrium yelds higher payoff than the one already found then new Nash equilibrium is saved
									if ((temp_payoff[1]+temp_payoff[2])>(temp_ne_payoff[1]+temp_ne_payoff[2]))
									{
										temp_ne_payoff <- temp_payoff
										temp_x_ne_s <- c(x_s1-1,x_s2-1)
										temp_y_ne_s <- c(y_s1-1,y_s2-1)
									}

									#If the group payoff analyzed Nash equilibrium yelds the same payoff, the equilibrium with lower difference between payoffs is saved
									if ((temp_payoff[1]+temp_payoff[2])==(temp_ne_payoff[1]+temp_ne_payoff[2]))
									{
										if(abs(temp_payoff[1]-temp_payoff[2])<abs(temp_ne_payoff[1]-temp_ne_payoff[2]))
										{
											temp_ne_payoff <- temp_payoff
											temp_x_ne_s <- c(x_s1-1,x_s2-1)
											temp_y_ne_s <- c(y_s1-1,y_s2-1)
										}
									}	
								}
							}
						}
					}
				}
			}

			if (ne_found==TRUE)
			{
				#Setting up values for Nash equilibrium for the given state
				result_class$rc_x_s_ne[j] <- list(temp_x_ne_s)
				result_class$rc_y_s_ne[j] <- list(temp_y_ne_s)
				result_class$rc_p_ne[j] <- list(temp_ne_payoff)
				result_class$rc_v_ne[j] <- list(temp_ne_payoff)
				result_class$rc_h_ne[j] <- list(c(min(i-1,temp_x_ne_s[1]+temp_x_ne_s[2]),min(j-1,temp_y_ne_s[1]+temp_y_ne_s[2])))
				result_class$rc_next_xy_ne[j]<-list(temp_next_xy)	

				
				#result_class$rc_h_ne[j] <- list(c(min(i-1,temp_x_ne_s[1]*HM_x_1+temp_x_ne_s[2]*HM_x_2),min(j-1,temp_y_ne_s[1]*HM_y_1+temp_y_ne_s[2]*HM_y_2)))		#Harvest cannot be higher than the population of the species
			}
			else
			{
				#Even if no solution is found the algorithm needs to have some numeric data to run
				result_class$rc_x_s_ne[j] <- list(c(-99,-99))
				result_class$rc_y_s_ne[j] <- list(c(-99,-99))
				result_class$rc_p_ne[j] <- list(c(-99,-99))
				result_class$rc_v_ne[j] <- list(c(-99,-99))
				result_class$rc_h_ne[j] <- list(c(-99,-99))
				result_class$rc_next_xy_ne[j]<-list(c(-99,-99))	
			}

			gc()
		}
		
		#Returning class with data
		return(result_class)
	}

	##Saving the results of parallel computations in the matrices
	for (it_i in 1:length(result_par))
	{
		for (it_j in 1:(k2+1))
		{
			x_s_ne[it_i,it_j,1] <- result_par[[it_i]]$rc_x_s_ne[it_j]
			y_s_ne[it_i,it_j,1] <- result_par[[it_i]]$rc_y_s_ne[it_j]
			p_ne[it_i,it_j,1] <- result_par[[it_i]]$rc_p_ne[it_j]
			v_ne[it_i,it_j,1] <- result_par[[it_i]]$rc_v_ne[it_j]
			h_ne[it_i,it_j,1] <- result_par[[it_i]]$rc_h_ne[it_j]
			next_xy_ne[it_i,it_j,1]<-result_par[[it_i]]$rc_next_xy_ne[it_j]
		}
	}

	#Backward induction analysis is conducted only if the game has more than one/initial iteration
	if (t0>0)
	{
		for (t in 2:(t0+1))
		{
			cat(paste0(round(t/(t0+2)*100),"%\n"))
			
			#Iterate over all possible states of the system
			#	Each process calculate all possible states of predator population for given prey population
			result_par <- foreach (i=1:(k1+1)) %dopar%
			{
				#Setting up the result class
				result_class <- resultClass()

				for (j in 1:(k1+1))
				{
					#Setting up data for the Nash equilibrium
					ne_found <- FALSE			#Information whether the Nash equilibrium was found
					temp_ne_payoff <- c(-1,-1)

					#Adjusting maximum strategy/harvesting effort to match the size of the populations
					temp_max_x_s1 <- min(max_x_s1,round((i-1)/HM_x_1))
					temp_max_x_s2 <- min(max_x_s2,round((i-1)/HM_x_2))
					temp_max_y_s1 <- min(max_y_s1,round((j-1)/HM_y_1))
					temp_max_y_s2 <- min(max_y_s2,round((j-1)/HM_y_2))
					temp_next_xy <- c(-1,-1)

					
					#Analysis of symmetric strategies - falsification method
					if(symmetric==TRUE)
					{
						#Setting up maximum symmetric strategies and symmetric value for harvest multipliers
						temp_max_x_sym <- min(temp_max_x_s1,temp_max_x_s2)
						temp_max_y_sym <- min(temp_max_y_s1,temp_max_y_s2)
						HM_x <- round((HM_x_1+HM_x_2)/2)
						HM_y <- round((HM_y_1+HM_y_2)/2)

						#Iteration over all possible sests of symmetric strategies
						for (x_sym in 1:(temp_max_x_sym+1))
						{
							for (y_sym in 1:(temp_max_y_sym+1))
							{
								#Setting up temporary payoff value for P1 and P2
								temp_payoff <- c(-1,-1)

								#Calculating harvest and myopic temporary payoff
								if ((2*(x_sym-1)*HM_x)<=(i-1))
								{
									temp_payoff <- c((x_sym-1)*HM_x*PM_x_1,(x_sym-1)*HM_x*PM_x_2)
									temp_x_h <- 2*(x_sym-1)*HM_x
								}
								else
								{
									temp_payoff <- c((x_sym-1)*(i-1)*HM_x/(2*(x_sym-1)*HM_x)*PM_x_1,(x_sym-1)*(i-1)*HM_x/(2*(x_sym-1)*HM_x)*PM_x_2)
									temp_x_h <- i-1
								}

								if ((2*(y_sym-1)*HM_y)<=(j-1))
								{
									temp_payoff <- temp_payoff + c((y_sym-1)*HM_y*PM_y_1,(y_sym-1)*HM_y*PM_y_2)
									temp_y_h <- 2*(y_sym-1)*HM_y
								}
								else
								{
									temp_payoff <- temp_payoff + c((y_sym-1)*(j-1)*HM_y/(2*(y_sym-1)*HM_y)*PM_y_1,(y_sym-1)*(j-1)*HM_y/(2*(y_sym-1)*HM_y)*PM_y_2)
									temp_y_h <- j-1
								}

								#Calculating next state for temporary solution
								x_h <- (i-1)-temp_x_h
								y_h <- (j-1)-temp_y_h
								xnew<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
								xnew<-max(0,xnew)
								xnew<-min(k1,xnew)
								ynew<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
								ynew<-max(0,ynew)
								ynew<-min(k2,ynew)

								#Updating payoff to reflect the next state Nash equilibrium
								temp_payoff <- temp_payoff + gamma*v_ne[[xnew+1,ynew+1,t-1]]

								#Setting initial value for the temp_ne_payoff
								if (x_sym==1 & y_sym==1)
								{
									temp_ne_payoff<-temp_payoff
									temp_x_ne_s <- c(x_sym-1,x_sym-1)
									temp_y_ne_s <- c(y_sym-1,y_sym-1)
									temp_next_xy <- c(xnew,ynew)
								}

								#Verification of the solution for P1
								real_ne <- TRUE		#initial assumption that the solution is Nash equilibrium
								
								#Setting up verification harvest values
								t_x_h <- -1
								t_y_h <- -1

								#Iteration over all possible strategies for P1
								t_xs1_max <- min(max_x_s1,round((i-1)/HM_x_1))
								
								for (t_xs1 in 1:(t_xs1_max+1))
								{
									if(real_ne==FALSE){break}		#Exit the loop if solution is not a Nash equlibrium

									t_ys1_max <- min(max_y_s1,round((j-1)/HM_y_1))

									for (t_ys1 in 1:(t_ys1_max+1))
									{
										#Calculating myopic verification payoff
										if (((t_xs1-1)*HM_x_1+(x_sym-1)*HM_x)<=(i-1))
										{
											t_p <- (t_xs1-1)*HM_x_1*PM_x_1
											t_x_h <- (t_xs1-1)*HM_x_1+(x_sym-1)*HM_x
										}
										else
										{
											t_p <- (t_xs1-1)*(i-1)*HM_x_1/((t_xs1-1)*HM_x_1+(x_sym-1)*HM_x)*PM_x_1
											t_x_h <- i-1
										}

										if (((t_ys1-1)*HM_y_1+(y_sym-1)*HM_y)<=(j-1))
										{
											t_p <- t_p + (t_ys1-1)*HM_y_1*PM_y_1
											t_y_h <- (t_ys1-1)*HM_y_1+(y_sym-1)*HM_y
										}
										else
										{
											t_p <- t_p + (t_ys1-1)*(j-1)*HM_y_1/((t_ys1-1)*HM_y_1+(y_sym-1)*HM_y)*PM_y_1
											t_y_h <- j-1
										}

										#Calculating next state for the verified solution
										x_h <- (i-1)-t_x_h
										y_h <- (j-1)-t_y_h
										x_n<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
										x_n<-max(0,x_n)
										x_n<-min(k1,x_n)
										y_n<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
										y_n<-max(0,y_n)
										y_n<-min(k2,y_n)

										#Updating payoff to reflect the next state Nash equilibrium
										t_p <- t_p + gamma*v_ne[[x_n+1,y_n+1,t-1]][1]

										#If verification payoff is higher than the initial payoff the solution is not a Nash equilibrium
										if (t_p>temp_payoff[1])
										{
											real_ne <- FALSE
											break
										}
									}
								}

								#Verification of the solution for the P2 is required only if the solution was not verified negatively
								if (real_ne == TRUE)
								{
									#Iteration over all possible strategies for P2
									t_xs2_max <- min(max_x_s2,round((i-1)/HM_x_2))
									
									for (t_xs2 in 1:(t_xs2_max+1))
									{
										if(real_ne==FALSE){break}		#Exit the loop if solution is not a Nash equlibrium

										t_ys2_max <- min(max_y_s2,round((j-1)/HM_y_2))

										for (t_ys2 in 1:(t_ys2_max+1))
										{
											#Calculating verification payoff
											if (((x_sym-1)*HM_x+(t_xs2-1)*HM_x_2)<=(i-1))
											{
												t_p <- (t_xs2-1)*HM_x_2*PM_x_2
												t_x_h <- (x_sym-1)*HM_x+(t_xs2-1)*HM_x_2
											}
											else
											{
												t_p <- (t_xs2-1)*(i-1)*HM_x_2/((x_sym-1)*HM_x+(t_xs2-1)*HM_x_2)*PM_x_2
												t_x_h <- i-1
											}

											if (((y_sym-1)*HM_y_2+(t_ys2-1)*HM_y)<=(j-1))
											{
												t_p <- t_p + (t_ys2-1)*HM_y_2*PM_y_2
												t_y_h <- (y_sym-1)*HM_y_2+(t_ys2-1)*HM_y
											}
											else
											{
												t_p <- t_p + (t_ys2-1)*(j-1)*HM_y_2/((y_sym-1)*HM_y+(t_ys2-1)*HM_y_2)*PM_y_2
												t_y_h <- j-1
											}

											#Calculating next state for the verified solution
											x_h <- (i-1)-t_x_h
											y_h <- (j-1)-t_y_h
											x_n<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
											x_n<-max(0,x_n)
											x_n<-min(k1,x_n)
											y_n<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
											y_n<-max(0,y_n)
											y_n<-min(k2,y_n)

											#Updating payoff to reflect the next state Nash equilibrium
											t_p <- t_p + gamma*v_ne[[x_n+1,y_n+1,t-1]][2]

											#If verification payoff is higher than the initial payoff the solution is not a Nash equilibrium
											if (t_p>temp_payoff[2])
											{
												real_ne <- FALSE
												break
											}
										}
									}
								}

								#If set of strategies generates Nash equilibrium it is compared against already found Nash equilibrium
								if (real_ne==TRUE)
								{
									ne_found <- TRUE
									#If the group payoff of analyzed Nash equilibrium yelds higher payoff than the one already found then new Nash equilibrium is saved
									if ((temp_payoff[1]+temp_payoff[2])>(temp_ne_payoff[1]+temp_ne_payoff[2]))
									{
										temp_ne_payoff <- temp_payoff
										temp_x_ne_s <- c(x_sym-1,x_sym-1)
										temp_y_ne_s <- c(y_sym-1,y_sym-1)
										temp_next_xy <- c(xnew,ynew)
									}
								}
							}
						}
					}

					#If there is no symmetric Nash equilibirum algorithm goes through all possible strategy sets - matrix method
					if(ne_found == FALSE)
					{
						#Best response matrix containing information about whether the combination of the strategies is the Nash equilibrium, size is based on the maximum possible harvesting effort
						#Information:
						#	-1 not checked
						#	 0 checked/not best response for specific player
						#	 1 checked/best response for specific player
						ne_matrix<-array(list(c(-1,-1)),dim=c(max_y_s1+1,max_x_s1+1,max_y_s2+1,max_x_s2+1))
						
						#Iterate over all possible strategies for i/x
						for (x_s1 in 1:(temp_max_x_s1+1))
						{
							for (x_s2 in 1:(temp_max_x_s2+1))
							{
								#Iterate over all possible strategies for j/y
								for (y_s1 in 1:(temp_max_y_s1+1))
								{
									for (y_s2 in 1:(temp_max_y_s2+1))
									{
										#Setting up temporary payoff value for P1 and P2
										temp_payoff <- c(-1,-1)

										#Setting up temporary harvest values
										temp_x_h <- -1
										temp_y_h <- -1

										#Calculating temporary payoff and harvest
										if (((x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2)<=(i-1))
										{
											temp_payoff <- c((x_s1-1)*HM_x_1*PM_x_1,(x_s2-1)*HM_x_2*PM_x_2)
											temp_x_h <- (x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2
										}
										else
										{
											temp_payoff <- c((x_s1-1)*(i-1)*HM_x_1/((x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2)*PM_x_1,(x_s2-1)*(i-1)*HM_x_2/((x_s1-1)*HM_x_1+(x_s2-1)*HM_x_2)*PM_x_2)
											temp_x_h <- i-1
										}

										if (((y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2)<=(j-1))
										{
											temp_payoff <- temp_payoff + c((y_s1-1)*HM_y_1*PM_y_1,(y_s2-1)*HM_y_2*PM_y_2)
											temp_y_h <- (y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2
										}
										else
										{
											temp_payoff <- temp_payoff + c((y_s1-1)*(j-1)*HM_y_1/((y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2)*PM_y_1,(y_s2-1)*(j-1)*HM_y_2/((y_s1-1)*HM_y_1+(y_s2-1)*HM_y_2)*PM_y_2)
											temp_y_h <- j-1
										}

										#Calculating next state for temporary solution
										x_h <- (i-1)-temp_x_h
										y_h <- (j-1)-temp_y_h
										xnew<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
										xnew<-max(0,xnew)
										xnew<-min(k1,xnew)
										ynew<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
										ynew<-max(0,ynew)
										ynew<-min(k2,ynew)

										#Updating payoff to reflect the next state Nash equilibrium
										temp_payoff <- temp_payoff + gamma*v_ne[[xnew+1,ynew+1,t-1]]

										#Setting initial value for the temp_ne_payoff
										if (x_s1==1 & x_s2==1 & y_s1==1 & y_s2==1)
										{
											temp_ne_payoff<-temp_payoff
											temp_x_ne_s <- c(x_s1-1,x_s2-1)
											temp_y_ne_s <- c(y_s1-1,y_s2-1)
											temp_next_xy <- c(xnew,ynew)
										}

										#Setting up verification harvest values
										t_x_h <- -1
										t_y_h <- -1

										#Verification for P1 is required only if the best reponses for this set of strategies were not calculated
										if (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][1]==-1)
										{
											#Setting up the maximum payoff updated over the verification
											current_max_payoff<-temp_payoff[1]

											#Best response vectors, storing indices of the the best responses
											x_br<-vector()
											y_br<-vector()

											#Iteration over all possible strategies for P1
											t_xs1_max <- min(max_x_s1,i-1)
													
											for (t_xs1 in 1:(t_xs1_max+1))
											{
												t_ys1_max <- min(max_y_s1,j-1)

												for (t_ys1 in 1:(t_ys1_max+1))
												{
													#Calculating myopic verification payoff
													if (((t_xs1-1)*HM_x_1+(x_s2-1)*HM_x_2)<=(i-1))
													{
														t_p <- (t_xs1-1)*HM_x_1*PM_x_1
														t_x_h <- (t_xs1-1)*HM_x_1+(x_s2-1)*HM_x_2
													}
													else
													{
														t_p <- (t_xs1-1)*(i-1)*HM_x_1/((t_xs1-1)*HM_x_1+(x_s2-1)*HM_x_2)*PM_x_1
														t_x_h <- i-1
													}

													if (((t_ys1-1)*HM_y_1+(y_s2-1)*HM_y_2)<=(j-1))
													{
														t_p <- t_p + (t_ys1-1)*HM_y_1*PM_y_1
														t_y_h <- (t_ys1-1)*HM_y_1+(y_s2-1)*HM_y_2
													}
													else
													{
														t_p <- t_p + (t_ys1-1)*(j-1)*HM_y_1/((t_ys1-1)*HM_y_1+(y_s2-1)*HM_y_2)*PM_y_1
														t_y_h <- j-1
													}

													#Calculating next state for the verified solution
													x_h <- (i-1)-t_x_h
													y_h <- (j-1)-t_y_h
													x_n<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
													x_n<-max(0,x_n)
													x_n<-min(k1,x_n)
													y_n<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
													y_n<-max(0,y_n)
													y_n<-min(k2,y_n)

													#Updating payoff to reflect the next state Nash equilibrium
													t_p <- t_p + gamma*v_ne[[x_n+1,y_n+1,t-1]][1]

													#If verification payoff is higher than the initial payoff the solution is not a Nash equilibrium
													if (t_p > current_max_payoff)
													{												
														current_max_payoff<-t_p
														x_br<-t_xs1
														y_br<-t_ys1
													}
													else
													{
														#If verification payoff is equal to the current max payoff, the index of the solution is added to the list
														if (t_p==current_max_payoff)
														{
															x_br<-c(x_br,t_xs1)
															y_br<-c(y_br,t_ys1)
														}
													}

													#Information that the specific solution was verified is added to the best response matrix
													ne_matrix[[t_ys1,t_xs1,y_s2,x_s2]][1]<-0
												}
											}

											#Setting up information about the best responses into best response matrix
											for (it in 1:length(x_br))
											{
												ne_matrix[[y_br[it],x_br[it],y_s2,x_s2]][1]<-1
											}
										}

										#Verification for P2 is conducted only if given solution is already best response for P1. If not, the solution is not a Nash equlilibrium
										if (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][1]==1)
										{
											#Verification for P2 is required only if the best reponses for this set of strategies were not calculated
											if (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][2]==-1)
											{
												#Setting up the maximum payoff updated over the verification
												current_max_payoff<-temp_payoff[2]

												#Best response vectors, storing indices of the the best responses
												x_br<-vector()
												y_br<-vector()

												#Iteration over all possible strategies for P2
												t_xs2_max <- min(max_x_s2,i-1)
														
												for (t_xs2 in 1:(t_xs2_max+1))
												{
													t_ys2_max <- min(max_y_s2,j-1)

													for (t_ys2 in 1:(t_ys2_max+1))
													{				
														#Calculating myopic verification payoff
														if (((x_s1-1)*HM_x_1+(t_xs2-1)*HM_x_2)<=(i-1))
														{
															t_p <- (t_xs2-1)*HM_x_2*PM_x_2
															t_x_h <- (x_s1-1)*HM_x_1+(t_xs2-1)*HM_x_2
														}
														else
														{
															t_p <- (t_xs2-1)*(i-1)*HM_x_2/((x_s1-1)*HM_x_1+(t_xs2-1)*HM_x_2)*PM_x_2
															t_x_h <- i-1
														}

														if (((y_s1-1)*HM_y_1+(t_ys2-1)*HM_y_2)<=(j-1))
														{
															t_p <- t_p + (t_ys2-1)*HM_y_2*PM_y_2
															t_y_h <- (y_s1-1)*HM_y_1+(t_ys2-1)*HM_y_2
														}
														else
														{
															t_p <- t_p + (t_ys2-1)*(j-1)*HM_y_2/((y_s1-1)*HM_y_1+(t_ys2-1)*HM_y_2)*PM_y_2
															t_y_h <- j-1
														}

														#Calculating next state for the verified solution
														x_h <- (i-1)-t_x_h
														y_h <- (j-1)-t_y_h
														x_n<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
														x_n<-max(0,x_n)
														x_n<-min(k1,x_n)
														y_n<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
														y_n<-max(0,y_n)
														y_n<-min(k2,y_n)

														#Updating payoff to reflect the next state Nash equilibrium
														t_p <- t_p + gamma*v_ne[[x_n+1,y_n+1,t-1]][2]

														#If verification payoff is higher than the initial payoff the solution is not a Nash equilibrium
														if (t_p>current_max_payoff)
														{
															current_max_payoff <- t_p
															x_br<-t_xs2
															y_br<-t_ys2
														}
														else
														{
															#If verification payoff is equal to the current max payoff, the index of the solution is added to the list
															if (t_p==current_max_payoff)
															{
																x_br<-c(x_br,t_xs2)
																y_br<-c(y_br,t_ys2)
															}
														}

														ne_matrix[[y_s1,x_s1,t_ys2,t_xs2]][2]<-0
													}
												}

												#Information that the specific solution was verified is added to the best response matrix
												for (it in 1:length(x_br))
												{
													ne_matrix[[y_s1,x_s1,y_br[it],x_br[it]]][2]<-1
												}
											}
										}
									
										#If set of strategies generates Nash equilibrium it is compared against already found Nash equilibrium
										if ((ne_matrix[[y_s1,x_s1,y_s2,x_s2]][1]==1) & (ne_matrix[[y_s1,x_s1,y_s2,x_s2]][2]==1))
										{
											ne_found <- TRUE

											#If the group payoff of analyzed Nash equilibrium yelds higher payoff than the one already found then new Nash equilibrium is saved
											if (temp_payoff[1]+temp_payoff[2]>temp_ne_payoff[1]+temp_ne_payoff[2])
											{
												temp_ne_payoff <- temp_payoff
												temp_x_ne_s <- c(x_s1-1,x_s2-1)
												temp_y_ne_s <- c(y_s1-1,y_s2-1)
												temp_next_xy <- c(xnew,ynew)											
											}
											else
											{
												#If the group payoff of the analyzed Nash equilibrium yeld same payoff as the one already found then the Nash equilibrium with lower payoff difference is saved
												if ((temp_payoff[1]+temp_payoff[2])==(temp_ne_payoff[1]+temp_ne_payoff[2]))
												{
													if(abs(temp_payoff[1]-temp_payoff[2])<abs(temp_ne_payoff[1]-temp_ne_payoff[2]))
													{
														temp_ne_payoff <- temp_payoff
														temp_x_ne_s <- c(x_s1-1,x_s2-1)
														temp_y_ne_s <- c(y_s1-1,y_s2-1)
														temp_next_xy <- c(xnew,ynew)
													}
												}
											}
										}
									}
								}
							}
						}
					}	

					if (ne_found==TRUE)
					{
						#Setting up values for Nash equilibrium for the given state
						result_class$rc_x_s_ne[j] <- list(temp_x_ne_s)
						result_class$rc_y_s_ne[j] <- list(temp_y_ne_s)

						#Setting up harvest momentary payoff value based on the selected
						if ((temp_y_ne_s[1]*HM_y_1+temp_y_ne_s[2]*HM_y_2)<=j-1)
						{
							momentary_payoff <- c(temp_y_ne_s[1]*HM_y_1*PM_y_1,temp_y_ne_s[2]*HM_y_2*PM_y_2)
							y_harvest <- temp_y_ne_s[1]*HM_y_1 + temp_y_ne_s[2]*HM_y_2
						}
						else
						{
							momentary_payoff <- c(temp_y_ne_s[1]*(j-1)*HM_y_1/(temp_y_ne_s[1]*HM_y_1+temp_y_ne_s[2]*HM_y_2)*PM_y_1,temp_y_ne_s[2]*(j-1)*HM_y_2/(temp_y_ne_s[1]*HM_y_1+temp_y_ne_s[2]*HM_y_2)*PM_y_2)
							y_harvest <- j-1
						}

						if ((temp_x_ne_s[1]*HM_x_1+temp_x_ne_s[2]*HM_x_2)<=i-1)
						{
							momentary_payoff <- momentary_payoff+c(temp_x_ne_s[1]*HM_x_1*PM_x_1,temp_x_ne_s[2]*HM_x_2*PM_x_2)
							x_harvest<-temp_x_ne_s[1]*HM_x_1 + temp_x_ne_s[2]*HM_x_2
						}
						else
						{
							momentary_payoff <- momentary_payoff+c((temp_x_ne_s[1])*(i-1)*HM_x_1/(temp_x_ne_s[1]*HM_x_1+temp_x_ne_s[2]*HM_x_2)*PM_x_1,(temp_x_ne_s[2])*(i-1)*HM_x_2/(temp_x_ne_s[1]*HM_x_1+temp_x_ne_s[2]*HM_x_2)*PM_x_2)
							x_harvest <- i-1
						}

						result_class$rc_p_ne[j] <- list(momentary_payoff)
						result_class$rc_v_ne[j] <- list(temp_ne_payoff)
						result_class$rc_h_ne[j] <- list(c(x_harvest,y_harvest))
						result_class$rc_next_xy_ne[j]<-list(temp_next_xy)
					}
					else
					{
						#Even if no solution is found the algorithm needs to have some numeric data to run
						result_class$rc_x_s_ne[j] <- list(c(-99,-99))
						result_class$rc_y_s_ne[j] <- list(c(-99,-99))
						result_class$rc_p_ne[j] <- list(c(-99,-99))
						result_class$rc_v_ne[j] <- list(c(-99,-99))
						result_class$rc_h_ne[j] <- list(c(-99,-99))
						result_class$rc_next_xy_ne[j] <- list(c(-99,-99))
					}

					gc()
				}

				#Returning class with data
				return(result_class)	
			}

			##Saving the results of parallel computations in the matrices
			for (it_i in 1:length(result_par))
			{
				for (it_j in 1:(k2+1))
				{
					x_s_ne[it_i,it_j,t] <- result_par[[it_i]]$rc_x_s_ne[it_j]
					y_s_ne[it_i,it_j,t] <- result_par[[it_i]]$rc_y_s_ne[it_j]
					p_ne[it_i,it_j,t] <- result_par[[it_i]]$rc_p_ne[it_j]
					v_ne[it_i,it_j,t] <- result_par[[it_i]]$rc_v_ne[it_j]
					h_ne[it_i,it_j,t] <- result_par[[it_i]]$rc_h_ne[it_j]
					next_xy_ne[it_i,it_j,t]<-result_par[[it_i]]$rc_next_xy_ne[it_j]
				}
			}
		}
	}

	#Stopping setup for parallel computations
	stopCluster(cl)

	cat(paste0("100%\n"))
	cat("\n")

	#Save the Nash equlibrium for the inital state of the sytem
	x<-x0
	y<-y0
	trajx_ne[1]<-x
	trajy_ne[1]<-y
	traj_x_h_ne[1]<-h_ne[[x+1,y+1,t0+1]][1]
	traj_y_h_ne[1]<-h_ne[[x+1,y+1,t0+1]][2]
	traj_x_s1_ne[1]<-x_s_ne[[x+1,y+1,t0+1]][1]
	traj_y_s1_ne[1]<-y_s_ne[[x+1,y+1,t0+1]][1]
	traj_x_s2_ne[1]<-x_s_ne[[x+1,y+1,t0+1]][2]
	traj_y_s2_ne[1]<-y_s_ne[[x+1,y+1,t0+1]][2]
	trajp1_ne[1]<-p_ne[[x+1,y+1,t0+1]][1]
	trajp2_ne[1]<-p_ne[[x+1,y+1,t0+1]][2]
	trajap1_ne<-trajp1_ne[1]
	trajap2_ne<-trajp2_ne[1]
	trajv1_ne[1]<-v_ne[[x+1,y+1,t0+1]][1]
	trajv2_ne[1]<-v_ne[[x+1,y+1,t0+1]][2]

	#Construction of the Nash equlibrium trajectory is conducted only if the game has more than one/initial iteration
	if (t0>0)
	{
		#Construct Nash equilibrium trajectory going through all iterations of the game
		for (t in 2:(t0+1))
		{
			x<-trajx_ne[t-1]
			y<-trajy_ne[t-1]
			x_h<-x-traj_x_h_ne[t-1]
			y_h<-y-traj_y_h_ne[t-1]
			xnew<-round(x_h+x_h*a*(1-x_h/k1)-x_h*b*y_h)
			xnew<-max(0,xnew)
			xnew<-min(k1,xnew)
			ynew<-round(y_h+y_h*c*x_h*(1-e*y_h/k2)-d*y_h)
			ynew<-max(0,ynew)
			ynew<-min(k2,ynew)
			trajx_ne[t]<-xnew
			trajy_ne[t]<-ynew
			traj_x_h_ne[t]<-h_ne[[xnew+1,ynew+1,t0+2-t]][1]
			traj_y_h_ne[t]<-h_ne[[xnew+1,ynew+1,t0+2-t]][2]
			traj_x_s1_ne[t]<-x_s_ne[[xnew+1,ynew+1,t0+2-t]][1]
			traj_y_s1_ne[t]<-y_s_ne[[xnew+1,ynew+1,t0+2-t]][1]
			traj_x_s2_ne[t]<-x_s_ne[[xnew+1,ynew+1,t0+2-t]][2]
			traj_y_s2_ne[t]<-y_s_ne[[xnew+1,ynew+1,t0+2-t]][2]

			#Save the value connected with current strategy, taking into account the discounting
			trajp1_ne[t]<-p_ne[[xnew+1,ynew+1,t0+2-t]][1]*(gamma^(t-1))
			trajp2_ne[t]<-p_ne[[xnew+1,ynew+1,t0+2-t]][2]*(gamma^(t-1))
			trajap1_ne[t]<-trajp1_ne[t]+trajap1_ne[t-1]
			trajap2_ne[t]<-trajp2_ne[t]+trajap2_ne[t-1]
			trajv1_ne[t]<-v_ne[[xnew+1,ynew+1,t0+2-t]][1]
			trajv2_ne[t]<-v_ne[[xnew+1,ynew+1,t0+2-t]][2]
		}
	}

	#Creating and returning dataframe from variables
	data_export <- cbind(time,trajx_ne,trajy_ne,traj_x_h_ne, traj_y_h_ne,traj_x_s1_ne,traj_y_s1_ne,traj_x_s2_ne,traj_y_s2_ne,trajp1_ne,trajp2_ne,trajap1_ne,trajap2_ne,trajv1_ne,trajv2_ne)
	colnames(data_export) <- c("t","x_ne","y_ne","x_h_ne","y_h_ne","x_s1_ne","y_s1_ne","x_s2_ne","y_s2_ne","p1_ne","p2_ne","ap1_ne","ap2_ne","v1_ne","v2_ne")
	data_frame <- as.data.frame(data_export)

	output <- list(data_frame,next_xy_ne)

	return(output)
}


save_output_social_optimum_and_nash_equilibrium <- function(fname,output_so,output_ne,parameters,line)
{
	#Arguments:
	#	fname - base filename
	#	results - trajectories data
	#	parameters - base parameters used to calculate zero-growth isoclines
	#	line - position of parameter variables in the dataframe vectors

	#Output:
	#	Note: All files are saved in the Data folder

	#Setting up parameters values from parameters dataframe
	k1<-parameters$k1[line]
	k2<-parameters$k2[line]
	max_x_s1 <-parameters$max_x_s1[line]
	max_x_s2 <-parameters$max_x_s2[line]
	max_y_s1 <-parameters$max_y_s1[line]
	max_y_s2 <-parameters$max_y_s2[line]
	max_x_ts <- min(max_x_s1+max_x_s2,k2)
	max_y_ts <- min(max_y_s1+max_y_s2,k2)
	a<-parameters$a[line]
	b<-parameters$b[line]
	c<-parameters$c[line]
	d<-parameters$d[line]
	e<-parameters$e[line]
	gamma<-parameters$gamma[line]
	t0<-parameters$t0[line]
	x0<-parameters$x0[line]
	y0<-parameters$y0[line]

	#Separating data from dataframes
	results_so <- output_so[[1]]
	results_ne <- output_ne[[1]]

	#Calculating zero-growth isoclines
	#	Note: Potential division by 0 generates values that are not depicted on the phase diagram
	isox<-array(0,dim=k1+1)
	isoy<-array(0,dim=k2+1)
	
	for (i in 0:k1){
		isox[i+1]<-(i*a*(1-i/k1))/(b*i)
	}
	for (j in 0:k2){
		isoy[j+1]<-d/((1-(e*j)/k2)*c)
	}

	#Saving metadata for the analysis
	metadata_file<-file(paste0("Data/",fname,"_metadata.txt"))
	timestep <- as.character(as.POSIXct(Sys.time()))
	writeLines(c(
		paste0("Filename: ",fname),
		paste0("Created: ",timestep),
		"",
		"Datafiles:",
		paste0(fname,"_trajectories.csv - trajectory data (no harvest, social optimum maximum harvest, Nash equilibrium maximum harvest)"),
		paste0(fname,"_isoclines.csv - zero-growth isoclines data (no harvest scenario)"),
		paste0(fname,"_parameters.csv - description of the model parameters"),
		paste0(fname,"_plots.png - group plots (no harvest, social optimum maximum harvest, Nash equilibrium maximum harvest) and phase diagrams for all scenarios)")), metadata_file)
	close(metadata_file)

	#Saving trajectories data into file
	data_export <- cbind(results_so$t,results_so$x_h0,results_so$y_h0,results_so$x_so,results_so$y_so,results_so$x_h_so,results_so$y_h_so,results_so$p_so,results_so$ap_so,results_so$v_so,results_ne$x_ne,results_ne$y_ne,results_ne$x_h_ne,results_ne$y_h_ne,results_ne$x_s1_ne,results_ne$x_s2_ne,results_ne$y_s1_ne,results_ne$y_s2_ne,results_ne$p1_ne,results_ne$p2_ne,results_ne$ap1_ne,results_ne$ap2_ne,results_ne$v1_ne,results_ne$v2_ne)
	colnames(data_export) <- c("t","x_h0","y_h0","x_so","y_so","x_h_so","y_h_so","p_so","ap_so","v_so","x_ne","y_ne","x_h_ne","y_h_ne","x_s1_ne","x_s2_ne","y_s1_ne","y_s2_ne","p1_ne","p2_ne","ap1_ne","ap2_ne","v1_ne","v2_ne")
	write.csv(data_export,paste0("Data/",fname,"_trajectories.csv"), row.names = FALSE)

	#Saving zero-growth isoclines data into file
	isoclines <- cbind(isox,isoy)
	colnames(isoclines) <- c("x_iso","y_iso")
	write.csv(isoclines,paste0("Data/",fname,"_isoclines.csv"), row.names = FALSE)

	#Saving parameters data into file
	write.csv(parameters[line,],paste0("Data/",fname,"_parameters.csv"), row.names = FALSE)

	#Saving plots into data files
	png(filename=paste0("Data/",fname,"_plots.png"),width = 15, height = 15, units="cm", res=300)

	par(mfrow=c(3,3))
	par(mar = c(3, 3, 1, 0.5))
    par(mgp=c(2,0.5,0))

	time<-1:(t0+1)

	#Plot: population dynamics without harvest
	plot(time,results_so$x_h0,type="l", col = "darkgreen",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="No harvest")
	lines(time,results_so$y_h0, col = "firebrick", lty=1)
	legend(x="topleft",c("Prey","Predator"),lty=c(1,1), col=c("darkgreen","firebrick"))

	#Plot: population dynamics under social optimum harvest
	plot(time,results_so$x_so,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Social optimum", col="darkgreen")
	lines(time,results_so$y_so, col="firebrick", lty=1)
	legend(x="topleft",c("Prey","Predator"),lty=c(1,1), col=c("darkgreen","firebrick"))

	#Plot: population dynamics under Nash equilibrium harvest
	plot(time,results_ne$x_ne,type="l",xlab="Season",ylab="Population size",ylim=c(0,max(k1,k2)), main="Nash equilibrium", col="darkgreen")
	lines(time,results_ne$y_ne, col="firebrick", lty=1)
	legend(x="topleft",c("Prey","Predator"),lty=c(1,1), col=c("darkgreen","firebrick"))

	#Plot: empty panel
	plot.new()

	#Plot: harvest under social optimum
	plot(time,results_so$y_h_so,type="l",lty=1,xlab="Season",ylab="Harvest",ylim=c(0,max(k1,k2)), main="Social optimum", col="cornflowerblue")

	#Plot: harvest under Nash equilibrium
	plot(time,results_ne$y_h_ne,type="l",lty=1,xlab="Season",ylab="Harvest",ylim=c(0,max(k1,k2)), main="Nash equilibrium", col="blue4")

	#Plot: phase diagram of population dynamics under no harvest, social optimum and Nash equilibrium
	plot(results_so$x_h0,results_so$y_h0,type="l", ,xlab="Prey",ylab="Predator",ylim=c(0,50), xlim=c(0,50),col="black", main="Phase diagram")
	points(results_so$x_h0[1],results_so$y_h0[1], pch=16, col="black")
	points(results_so$x_h0[t0+1],results_so$y_h0[t0+1], pch=4, col="black")
	lines(results_so$x_so,results_so$y_so,col="cornflowerblue")
	points(results_so$x_so[1],results_so$y_so[1], pch=16, col="cornflowerblue")
	points(results_so$x_so[t0+1],results_so$y_so[t0+1], pch=4, col="cornflowerblue")
	lines(results_ne$x_ne,results_ne$y_ne,col="blue4")
	points(results_ne$x_ne[1],results_ne$y_ne[1], pch=16,col="blue4")
	points(results_ne$x_ne[t0+1],results_ne$y_ne[t0+1], pch=4,col="blue4")
	lines(c(0:k2),isox, lty=3)
	lines(isoy,c(0:k1), lty=3)
	legend(x="topleft",c("No harvest","Social optimum","Nash equilibrium"),lty=c(1,1,1), col=c("black","cornflowerblue","blue4"))

	#Plot: empty panel
	plot.new()

	#Plot: players' strategies under Nash equilibrium
	#Version with predator harvesting
	plot(time,results_ne$y_s1_ne,type="l",lty=1,xlab="Season",ylab="Strategy",ylim=c(0,max(k1,k2)), main="Nash equilibrium", col="blue4")
	lines(time,results_ne$y_s2_ne, col="cadetblue1", lty=3)
	legend(x="topleft",c("Player 1","Player 2"),lty=c(1,3), col=c("blue4","cadetblue1"))

	dev.off()
}

## Package management functions

load_packages <- function()
{
	load_parallel_packages()
	load_graphic_packages()
}

load_parallel_packages <- function()
{
	#Code install package and loads library neccessary for parallel computations

	install.packages("doParallel")
	library(doParallel)
}

load_graphic_packages <- function()
{
	#Code install package and loads library neccesary for plots

	install.packages("RColorBrewer")
	library(RColorBrewer)
}


## Parameters managment functions

save_parameters <- function(fname,parameters)
{
	#Function saves the parameters dataframe into csv file

	#Arguments
	#	fname - name of the csv file to write the parameters, extension and description of the file is added automatically
	#	parameters - dataframe with parameters

	#Output
	#	fname_parameters.csv - csv file with parameters, each line of the file represent the line of parameters, specific parameters are placed in the columns

	write.csv(parameters,paste0(fname,"_parameters.csv"), row.names = FALSE)
}

load_parameters <- function(fname)
{
	#Function reads the parameters dataframe from the csv file

	#Arguments
	#	fname - name of the csv file with the parameters

	#Output
	#	parameters - dataframe with default values of parameters

	parameters <- read.csv(fname)
	return(parameters)
}