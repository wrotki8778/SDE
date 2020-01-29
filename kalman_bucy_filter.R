require(yuima)
T<-3
n<-2^10
a<-function(t){return(1/2)}
A<-function(t){return(1/2)}
b<-function(t){return(1)}
B<-function(t){return(1)}
samp <- setSampling(Terminal=T,n=n)
# Wielowymiarowe SDE
# A multi-dimensional (correlated) diffusion process.
# To describe the following model:
# X=(X1,X2,X3); dXt = U(t,Xt)dt + V(t)dWt
# For drift coeffcient
U <- c("a(t)*x1","A(t)*x1")
# coefficient matrix for diffusion term
V <- matrix( c( "b(t)",
                "0",
                "0",
                "B(t)"
), 2, 2)
# Model sde using "setModel" function
cor.mod <- setModel(drift = U, diffusion = V, solve.variable=c("x1","x2") )
str(cor.mod)
# Set the `yuima' object.
cor.samp <- setSampling(Terminal=T, n=n)
cor <- setYuima(model=cor.mod, sampling=cor.samp)
# Solve SDEs using Euler-Maruyama method.
set.seed(173)
cor <- simulate(cor,xinit=c(1,1))
plot(cor)
dane<-cor@data@original.data
czasy<-samp@grid[[1]]
sygnal<-dane[,1]
obserwacja<-dane[,2]

#delta=sqrt(4*a^2 + 4*c^2)
#beta<-(-2*a+delta)/(-2*c^2)
#alpha<-(-2*a-delta)/(-2*c^2)

#gamma<-(-beta)/(alpha)
#lambda<-c^2*(alpha-beta)

v<-c(0)
tmp=0
for(i in 2:length(czasy)){
  t=czasy[i-1]
  dt=czasy[i]-czasy[i-1]
  tmp=tmp+(2*a(t)*tmp+b(t)^2-((A(t)*tmp)/B(t))^2)*dt
  v<-append(v,tmp)
}

odszumiony<-c(1)
tmp=1
for(i in 2:length(czasy)){
  t=czasy[i-1]
  dt=czasy[i]-czasy[i-1]
  dYt=obserwacja[i]-obserwacja[i-1]
  tmp=tmp+(a(t)-v[i-1]*(A(t)/B(t))^2)*tmp*dt+v[i-1]*(A(t)/B(t)^2)*dYt
  odszumiony<-append(odszumiony,tmp)
}

nowe_dane<-ts(data=cbind(sygnal,obserwacja,odszumiony))
plot(nowe_dane)
plot(czasy,obserwacja,col='red',type='l',ylim=c(0,max(sygnal,odszumiony,obserwacja)))
lines(czasy,sygnal,col='blue',type='l')
lines(czasy,odszumiony,col='green',type='l')

plot(czasy,v,type='l')