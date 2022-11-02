##Note: Parameters have to be changed manualy

# optimal harvesting of predators
# k1 - carrying capacity (prey)
# k2 - carrying capacity (predators)
# a,b,c,d - parameters for birth death process
# t0- time horizon
# v[x,y,t] - optimal value starting in state (x*eps,y*eps) with t-1 moments to go
# eps - precision of measuring x,y
# hopt[x,y,t] - optimal action in state (x,y) with t-1 moments to go
# r - discount rate
## vv value function as a vector
k1<-50
k2<-50
e1<-1
e2<-0.01
# (xm+1), (ym+1) number of x,y values in the web
xm<-round(k1/e1)
ym<-round(k2/e1)
am<-round(1+k2/e2)
a<-1.1
b<-0.07
c<-0.7
d<-0.05
r<-1
# t0 no. of seasons
t0<-100
th<-t0+1
time<-1:th
# the trajectory under optimal control
xt<-array(0,dim=th)
yt<-array(0,dim=th)
at<-array(0,dim=th)
xa<-array(0,dim=th)
ya<-array(0,dim=th)
aa<-array(0,dim=th)
## v - array of values
xma<-xm+1
yma<-ym+1
v<-array(0,dim=c(xma,yma,th))
# hopt - array of optimal actions
# traj - trajectories of populations and harvest
hopt<-array(-1,dim=c(xma,yma,th))
trajx<-array(0,dim=c(xma,yma,th))
trajy<-array(0,dim=c(xma,yma,th))
trajh<-array(0,dim=c(xma,yma,th))
# populations and harvest in the middle of the period 
x50<-array(0,dim=c(xma,yma))
y50<-array(0,dim=c(xma,yma))
h50<-array(0,dim=c(xma,yma))
# web of x and y values
x<-(0:xm)*e1
y<-(0:ym)*e1
xt<-(0:xm)*e1
yt<-(0:ym)*e1
# in the final season it is optimal to harvest all the predators (unique)
for (i in 1:xma){
	for (j in 1:yma){
		v[i,j,1]<-y[j]
		hopt[i,j,1]<-y[j]
	}
}
for (i in 1:xma){
	for (j in 1:yma){
		c0<-r*y[j]*(1-c+d*x[i]-d*x[i]*y[j]/k2)
	c1<-1-r*(1-c)-r*d*x[i]+2*r*d*x[i]*y[j]/k2
	c2<-r*d*x[i]/k2
	# a0 - value of a for which the value function is maximized (by differentiation)
	a0<-0.5*c1/c2
	# if a0 is between 0 and the predator population, it is the optimal action
	if (a0<y[j] && a0>0){
		hopt[i,j,2]<-a0
		v[i,j,2]<-c0+c1*a0-c2*a0^2
	}
	# Future Expected Reward for maximum harvest
	else { # the optimal action must be extreme (i.e. harvest maximally or harvest nothing)
	  f1<-c0+c1*y[j]-c2*y[j]^2
	  	if (c0<y[j]){ # in this case, it is better to harvest all than nothing
			hopt[i,j,2]<-y[j]
			v[i,j,2]<-y[j]
			}
			else { # when it is better to harvest nothing than harvest everything
				hopt[i,j,2]<-0
				v[i,j,2]<-c0
			}
	}			
	}
}
# we now use quadratic interpolation in order to estimate the 
# future value function over the web of points
for (t in 3:t0){
  for (i in 1:xma){
    v[i,1,t]<-0
    hopt[i,1,t]<-0
  }
  for (j in 2:yma){
    v[1,j,t]<-y[j]
    hopt[1,j,t]<-y[j]
  }
for (i in 2:xma){
 for (j in 2:yma){
   # possible actions 
   na<-round(1+y[j]/e2)
   #v[i,j,t]<-0
   #hopt[i,j,t]<-0
   for (k in 1:na){
   act<-(k-1)*e2
   # calculate the next state given the actions
   z<-y[j]-act
   xp<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*z)
   xp<-max(0,xp)
   yp<-min(k2,z*(1-c)+d*x[i]*z*(1-z/k2))
   yp<-max(0,yp)
   ipl<-floor(1+xp/e1)
   ipl<-min(xma-1,ipl)
   iph<-ipl+1
   xe<-(xp-x[ipl])/e1
   jpl<-floor(1+yp/e1)
   jpl<-min(yma-1,jpl)
   jph<-jpl+1
   ye<-(yp-y[jpl])/e1
   v1<-v[ipl,jpl,t-1]
   v2<-v[iph,jpl,t-1]
   v3<-v[ipl,jph,t-1]
   v4<-v[iph,jph,t-1]
   vp<-act+r*(v1+xe*(v2-v1)+ye*(v3-v1)+xe*ye*(v4-v2-v3+v1))
   if (vp>v[i,j,t]){
     v[i,j,t]<-vp
     hopt[i,j,t]<-act
   }
   }
 }
}
}
# plot of optimal harvest and value when prey population is 25, predator population y varies
#plot(y, hopt[251,1:(ym+1),t], type="l", ylab="opt. harvest t steps left")
#lines(y, v[251,1:(ym+1),t],col="red")
# trace the approximate trajectory of the population for starting population (11,10)
xt[1]<-13
yt[1]<-10
for (iii in 1:t0){
  ie<-min(xma,1+xt[iii]/e1)
  ie<-max(1,ie)
  i<-trunc(ie)
  je<-min(yma,1+yt[iii]/e1)
  je<-max(1,je)
  j<-trunc(je)
  if (i<xma & j<yma){
    u0<-hopt[i,j,t0+1-iii]
    u1<-hopt[i+1,j,t0+1-iii]
    u2<-hopt[i,j+1,t0+1-iii]
    u3<-hopt[i+1,j+1,t0+1-iii]
    at[iii]<-u0+(ie-i)*(u1-u0)+(je-j)*(u2-u0)+(ie-i)*(je-j)*(u3-u2-u1+u0)
  }
  if (i<xma & j==xma){
    u0<-hopt[i,j,t0+1-iii]
    u1<-hopt[i+1,j,t0+1-iii]
    at[iii]<-u0+(ie-i)*(u1-u0)
  }
  if (i==xma & j<xma){
    u0<-hopt[i,j,t0+1-iii]
    u1<-hopt[i,j+1,t0+1-iii]
    at[iii]<-u0+(je-j)*(u2-u0)
  }
  if (i==xma & j==xma){
    at[iii]<-hopt[i,j,t0+1-iii]
  }
  xt[iii+1]<-xt[iii]+a*xt[iii]*(1-xt[iii]/k1)-b*xt[iii]*(yt[iii]-at[iii])
  yt[iii+1]<-(1-c)*(yt[iii]-at[iii])+d*xt[iii]*(yt[iii]-at[iii])*(1-(yt[iii]-at[iii])/k2)
}
# plot of trajectories
plot(1:t0,xt[1:t0], type="l", ylab="population",xlab="Season",ylim=c(0,70))
lines(1:t0,yt[1:t0],col="red")
lines(1:t0,at[1:t0],col="blue")
legend('topleft',legend=c("prey","predator","harvest"),lty=c(1,1,1),col=c("black","red","blue"),cex=0.5)
#for (ii in 1:xm){
# for (jj in 1:ym){
#  xa[1]<-x[ii]
# ya[1]<-y[jj]
#for (iii in 1:(t0/2)){
# i<-max(1,round(1+(xa[iii]-1)/eps))
#  i<-min(i,xm+1)
# j<-max(1,round(1+(ya[iii]-1)/eps))
#j<-min(j,ym+1)
#  aa[iii]<-hopt[i,j,t0+1-iii]
#  xa[iii+1]<-xa[iii]+a*xa[iii]*(1-xa[iii]/k1)-b*xa[iii]*(ya[iii]-aa[iii])
#  ya[iii+1]<-(1-c)*(ya[iii]-aa[iii])+d*xa[iii]*(ya[iii]-aa[iii])*(1-(ya[iii]-aa[iii])/k2)
#}
#   x50[ii,jj]<-xa[iii]
#  y50[ii,jj]<-ya[iii]
#  h50[ii,jj]<-aa[iii]
#  }
#}
# function defining the value of the regression function on the 
# basis of the regression coefficients
#regfun<-function(x,y,arg1,arg2,arg3,arg4,arg5,arg6) {
# regfun <- arg1+arg2*x+arg3*y+arg4*x*y+arg5*x^2+arg6*y^2
#return(regfun)
#}
df<-data.frame(xt,yt,at)
write.csv2(df,"C:\\Users\\david\\OneDrive\\Dokumenty\\pajak\\scenario0type2.csv")
