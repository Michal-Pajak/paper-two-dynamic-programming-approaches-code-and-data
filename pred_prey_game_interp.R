#NOTE: Parameters have to be changed manualy

# optimal harvesting of predators
# k1 - carrying capacity (prey)
# k2 - carrying capacity (predators)
# a,b,c,d - parameters for birth death process
# t0- time horizon
# v[x,y,t] - optimal value starting in state (x*eps,y*eps) with t-1 moments to go
# eps - precision of measuring x,y
# hopt[x,y,t] - optimal action in state (x,y) with t-1 moments to go
# num[t] - number of symmetric equilibria at time t
# r - discount rate
## vv value function as a vector
k1<-50
k2<-50
e1<-1
e2<-0.01
# (xm+1), (ym+1) number of x,y values in the ?web?
xm<-round(k1/e1)
ym<-round(k2/e1)
am<-round(1+k2/(2*e2))
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
v<-array(0,dim=c(xma,yma,t0))
num<-array(0,dim=c(xma,yma,t0))
num[1:xma,1:yma,1:2]<-1
# hopt - array of optimal actions
# traj - trajectories of populations and harvest
hequ<-array(-1,dim=c(xma,yma,t0))
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
		v[i,j,1]<-y[j]/2
		hequ[i,j,1]<-y[j]/2
	}
}
# in the penultimate season calculation of the equilibrium action a is 
# relatively easy - max{c0+c1*a-c2*a^2} - see notes
for (i in 1:xma){
	for (j in 1:yma){
		c0<-r*y[j]*(1-c+d*x[i]-d*x[i]*y[j]/k2)
	c1<-1-r*(1-c)-r*d*x[i]+2*r*d*x[i]*y[j]/k2
	s<-1-c1
	c2<-r*d*x[i]/k2
	# a0 - equilibrium value of a  (by differentiation)
	a0<-(2-s)/(4*c2)
	# if 2*a0 is between 0 and the predator population, it is the equilibrium action
	if ((a0<y[j]/2) && (a0>0)){
		hequ[i,j,2]<-a0
		v[i,j,2]<-(c0+2*c1*a0-4*c2*a0^2)/2
	}
	# Future Expected Reward for maximum harvest
	else { # the optimal action must be extreme (i.e. harvest maximally or harvest nothing)
	  	if (c0<y[j]){ # in this case, it is better to harvest all than nothing
			hequ[i,j,2]<-y[j]/2
			v[i,j,2]<-y[j]/2
			}
			else { # when it is better to harvest nothing than harvest everything
				hequ[i,j,2]<-0
				v[i,j,2]<-c0/2
			}
	}			
	}
}
# we now use interpolation in order to estimate the 
# future value function 
# We look for a local Nash equilibrium
for (t in 3:t0){
  for (i in 1:xma){
    v[i,1,t]<-0
    hequ[i,1,t]<-0
    num[i,1,t]<-1
  }
  for (j in 2:yma){
    v[1,j,t]<-y[j]/2
    hequ[1,j,t]<-y[j]/2
    num[1,j,t]<-1
  }
  for (i in 2:xma){
    for (j in 2:yma){
      # possible actions 
      na<-round(1.0000001+y[j]/(2*e2))
      vmax<-0
      #v[i,j,t]<-0
      #hequ[i,j,t]<-0
      # !!!!! check whether extreme harvests form an equilibrium !!!!!
      # calculate the next state given the pair of actions (0,0)
      z<-y[j]
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
      vp<-r*(v1+xe*(v2-v1)+ye*(v3-v1)+xe*ye*(v4-v2-v3+v1))
      # calculate the next state given the pair of actions (e2,0)
      zu<-y[j]-e2
      au<-e2
      xpu<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*zu)
      xpu<-max(0,xpu)
      ypu<-min(k2,zu*(1-c)+d*x[i]*zu*(1-zu/k2))
      ypu<-max(0,ypu)
      iplu<-floor(1+xpu/e1)
      iplu<-min(xma-1,iplu)
      iphu<-iplu+1
      xeu<-(xpu-x[iplu])/e1
      jplu<-floor(1+ypu/e1)
      jplu<-min(yma-1,jplu)
      jphu<-jplu+1
      yeu<-(ypu-y[jplu])/e1
      v1u<-v[iplu,jplu,t-1]
      v2u<-v[iphu,jplu,t-1]
      v3u<-v[iplu,jphu,t-1]
      v4u<-v[iphu,jphu,t-1]
      vpu<-au+r*(v1u+xeu*(v2u-v1u)+yeu*(v3u-v1u)+xeu*yeu*(v4u-v2u-v3u+v1u))
      # calculate the next state given the pair of actions (e2,e2)
      zu2<-y[j]-2*e2
      xpu2<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*zu2)
      xpu2<-max(0,xpu2)
      ypu2<-min(k2,zu2*(1-c)+d*x[i]*zu2*(1-zu2/k2))
      ypu2<-max(0,ypu2)
      iplu2<-floor(1+xpu2/e1)
      iplu2<-min(xma-1,iplu2)
      iphu2<-iplu2+1
      xeu2<-(xpu2-x[iplu2])/e1
      jplu2<-floor(1+ypu2/e1)
      jplu2<-min(yma-1,jplu2)
      jphu2<-jplu2+1
      yeu2<-(ypu2-y[jplu2])/e1
      v1u2<-v[iplu2,jplu2,t-1]
      v2u2<-v[iphu2,jplu2,t-1]
      v3u2<-v[iplu2,jphu2,t-1]
      v4u2<-v[iphu2,jphu2,t-1]
      vpu2<-au+r*(v1u2+xeu2*(v2u2-v1u2)+yeu2*(v3u2-v1u2)+xeu2*yeu2*(v4u2-v2u2-v3u2+v1u2))
      # check whether (0,0) is a pure equilibrium
      if (vp>=vpu){
        num[i,j,t]<-num[i,j,t]+1
        vmax<-vp
        v[i,j,t]<-vp
        hequ[i,j,t]<-0
      }
      # check whether (0,0) is part of a mixed equilibrium
      if ((vpu>vp)&&(vpu-e2>vpu2)){
        p<-(vpu-e2-vpu2)/(vpu-e2-vpu2+vpu-vp)
        num[i,j,t]<-num[i,j,t]+1
        vmax<-p*vp+(1-p)*(vpu-e2)
        v[i,j,t]<-vmax
        hequ[i,j,t]<-e2*(1-p)
      }
      for (k in 2:(na-1)){
        act<-(k-1)*e2
        ad<-act-e2
        au<-act+e2
        # calculate the next state given the pair of actions (act,act)
        z<-y[j]-2*act
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
        # calculate the next state given the pair of actions (act-e2,act)
        zd<-y[j]-2*act+e2
      xpd<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*zd)
       xpd<-max(0,xpd)
       ypd<-min(k2,zd*(1-c)+d*x[i]*zd*(1-zd/k2))
       ypd<-max(0,ypd)
        ipld<-floor(1+xpd/e1)
       ipld<-min(xma-1,ipld)
        iphd<-ipld+1
        xed<-(xpd-x[ipld])/e1
        jpld<-floor(1+ypd/e1)
        jpld<-min(yma-1,jpld)
        jphd<-jpld+1
        yed<-(ypd-y[jpld])/e1
        v1d<-v[ipld,jpld,t-1]
        v2d<-v[iphd,jpld,t-1]
        v3d<-v[ipld,jphd,t-1]
        v4d<-v[iphd,jphd,t-1]
        vpd<-ad+r*(v1d+xed*(v2d-v1d)+yed*(v3d-v1d)+xed*yed*(v4d-v2d-v3d+v1d))
        # calculate the next state given the pair of actions (act+e2,act)
        zu<-y[j]-2*act-e2
        xpu<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*zu)
        xpu<-max(0,xpu)
        ypu<-min(k2,zu*(1-c)+d*x[i]*zu*(1-zu/k2))
        ypu<-max(0,ypu)
        iplu<-floor(1+xpu/e1)
        iplu<-min(xma-1,iplu)
        iphu<-iplu+1
        xeu<-(xpu-x[iplu])/e1
        jplu<-floor(1+ypu/e1)
        jplu<-min(yma-1,jplu)
        jphu<-jplu+1
        yeu<-(ypu-y[jplu])/e1
        v1u<-v[iplu,jplu,t-1]
        v2u<-v[iphu,jplu,t-1]
        v3u<-v[iplu,jphu,t-1]
        v4u<-v[iphu,jphu,t-1]
        vpu<-au+r*(v1u+xeu*(v2u-v1u)+yeu*(v3u-v1u)+xeu*yeu*(v4u-v2u-v3u+v1u))
        # calculate the next state given the pair of actions (act+e2,act+e2)
        zu2<-y[j]-2*act-2*e2
        xpu2<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*zu2)
        xpu2<-max(0,xpu2)
        ypu2<-min(k2,zu2*(1-c)+d*x[i]*zu2*(1-zu2/k2))
        ypu2<-max(0,ypu2)
        iplu2<-floor(1+xpu2/e1)
        iplu2<-min(xma-1,iplu2)
        iphu2<-iplu2+1
        xeu2<-(xpu2-x[iplu2])/e1
        jplu2<-floor(1+ypu2/e1)
        jplu2<-min(yma-1,jplu2)
        jphu2<-jplu2+1
        yeu2<-(ypu2-y[jplu2])/e1
        v1u2<-v[iplu2,jplu2,t-1]
        v2u2<-v[iphu2,jplu2,t-1]
        v3u2<-v[iplu2,jphu2,t-1]
        v4u2<-v[iphu2,jphu2,t-1]
        vpu2<-au+r*(v1u2+xeu2*(v2u2-v1u2)+yeu2*(v3u2-v1u2)+xeu2*yeu2*(v4u2-v2u2-v3u2+v1u2))
        # check whether (act,act) is a pure equilibrium
                if ((vp>=vpu) && (vp>=vpd)){
                  num[i,j,t]<-num[i,j,t]+1
                  if (vp>vmax){
                    vmax<-vp
                  v[i,j,t]<-vp
                  hequ[i,j,t]<-act
                  }
                }
        # check whether (act,act) is part of a mixed equilibrium
        if ((vpu>vp)&&(vpu-e2>vpu2)){
          p<-(vpu-e2-vpu2)/(vpu-e2-vpu2+vpu-vp)
          num[i,j,t]<-num[i,j,t]+1
          veq<-p*vp+(1-p)*(vpu-e2)
          if (veq>vmax){
          vmax<-veq
          v[i,j,t]<-veq
          hequ[i,j,t]<-act+e2*(1-p)
        } 
        }
      }
      # check whether maximum harvest is an equilibrium
      act<-y[j]/2
      ad<-(na-2)*e2
      vp<-act
      # calculate the next state given the pair of actions (ad,act)
      zd<-act-(na-2)*e2
      xpd<-min(k1,x[i]*(1+a)-a*x[i]^2/k1-b*x[i]*zd)
      xpd<-max(0,xpd)
      ypd<-min(k2,zd*(1-c)+d*x[i]*zd*(1-zd/k2))
      ypd<-max(0,ypd)
      ipld<-floor(1+xpd/e1)
      ipld<-min(xma-1,ipld)
      iphd<-ipld+1
      xed<-(xpd-x[ipld])/e1
      jpld<-floor(1+ypd/e1)
      jpld<-min(yma-1,jpld)
      jphd<-jpld+1
      yed<-(ypd-y[jpld])/e1
      v1d<-v[ipld,jpld,t-1]
      v2d<-v[iphd,jpld,t-1]
      v3d<-v[ipld,jphd,t-1]
      v4d<-v[iphd,jphd,t-1]
      vpd<-ad+r*(v1d+xed*(v2d-v1d)+yed*(v3d-v1d)+xed*yed*(v4d-v2d-v3d+v1d))
      if (vp>vpd){
        num[i,j,t]<-num[i,j,t]+1
        if (vp>vmax){
          vmax<-vp
        v[i,j,t]<-vp
        hequ[i,j,t]<-act
        }
      }
      if (num[i,j,t]==0){
        if (v[i-1,j,t]>=v[i,j-1,t]){
        v[i,j,t]<-v[i-1,j,t]
        hequ[i,j,t]<-hequ[i-1,j,t]
        }
        if (v[i-1,j,t]<v[i,j-1,t]){
          v[i,j,t]<-v[i,j-1,t]
          hequ[i,j,t]<-hequ[i,j-1,t]
        }
      }
    }
  }
}
 # -----------------------
xt[1]<-11
yt[1]<-10
for (iii in 1:t0){
  ie<-min(xma,1+xt[iii]/e1)
  ie<-max(1,ie)
  i<-trunc(ie)
  je<-min(yma,1+yt[iii]/e1)
  je<-max(1,je)
  j<-trunc(je)
  if (i<xma & j<yma){
    u0<-hequ[i,j,t0+1-iii]
    u1<-hequ[i+1,j,t0+1-iii]
    u2<-hequ[i,j+1,t0+1-iii]
    u3<-hequ[i+1,j+1,t0+1-iii]
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
  xt[iii+1]<-min(k1,xt[iii]+a*xt[iii]*(1-xt[iii]/k1)-b*xt[iii]*(yt[iii]-2*at[iii]))
  xt[iii+1]<-max(0,xt[iii+1])
  yt[iii+1]<-min(k2,(1-c)*(yt[iii]-2*at[iii])+d*xt[iii]*(yt[iii]-2*at[iii])*(1-(yt[iii]-2*at[iii])/k2))
  yt[iii+1]<-max(0,yt[iii+1])
}
# plot of trajectories
plot(1:th,xt[1:th], type="l", ylab="population",xlab="Season",ylim=c(0,70))
lines(1:th,yt[1:th],col="red")
lines(1:th,2*at[1:th],col="blue")
legend('topleft',legend=c("prey","predator","total harvest"),lty=c(1,1,1),col=c("black","red","blue"),cex=0.5)

#NOTE: Filename has to be changed manually
df<-data.frame(xt,yt,at)
write.csv2(df,"Data/scenario0type2.csv")

