                     #//For Poisson Model//

# Vector of observed number of cysts
observed <-c(65 ,14, 10, 6 ,4 ,2 ,2 ,2 ,1, 1, 1, 2,0,0,0,0,0,0,0,1)

# Relative Frequency
observed_prob<-observed/111

#Summarising the Dataset
table(observed)
summary(observed)

#Maximum Likelihood parameter estimation

likelihood<- rep(0,20)

g<-function(x)
{
  
  return(observed_prob[x+1])
  
}





Maximumlikelihood_parameter<- function(o)
{
  
  
  for(i in 0:19)
  {
    
    if(g(i)==0)
    {
      likelihood[i+1]<-0
      
    }
    else
    {
      
      likelihood[i+1]<-  g(i)*log10( g(i)*factorial(i) /(exp(-o)*(o^i) )  )
    }
    
    
    
    
  }
  
  t<-sum(likelihood)
  return(t)
  
}

optimize(Maximumlikelihood_parameter,c(0,10))
#1. Power Divergence Optimisation


#Considering values after 5 i.e Median(+/-) 1.5*IQR as outliers and optimising Power Divergence equation

#fixing the value of tuning parameter lambda

r<-1.5

#performing the optimisation for the expression of power divergence
 
term2 <-rep(0,20)
 g<-function(x)
  {
    
    return(observed_prob[x+1])
    
  }





Best_parameter<- function(o)
{
  
  
  for(i in 0:19)
  {
    
    
    
    term2[i+1]<-  (((g(i))^( 1/r) )* ( factorial(i)^( (1/r) -1  ) ) ) /( ((exp(-o))*(o^i) )^( (1/r) -1  )) 
    
    
    
    
    
  }
  
  t<-sum(term2)
  return(t)
  
}

optimize(Best_parameter,c(0,100),maximum = TRUE) #We get 0.73 as the optimised value of parameter

#*************************************************************************************************************#

#2. Density Power Divergence Optimisation


#fixing the value of tuning parameter lambda
r<- 1.2

#performing the optimisation for the expression of density power divergence

term3 <-rep(0,20)
g<-function(x)
{
  return(ifelse(x > 19,0,observed_prob[x+1]) )
  
}




Best_Density<- function(o)
{
  
  
  for(i in 0:100)
  {
    
    
    
    term3[i+1]<-  ( ( (  exp(-o)*(o^i) )/factorial(i)  )^( r )) - ( 1+(1/(r-1) ))*(   g(i)*( ((  exp(-o)*(o^i) )/factorial(i)  )^(r-1) )  )   + (1/(r-1) )*( (   g(i)   )^(r) )
    
    
    
    
    
  }
  
  t<-sum(term3)
  return(t)
  
}

optimize(Best_Density,c(0,10))


#3.Log Density Power Divergence Optimisation  

#fixing the value of tuning parameter lambda
r<- 1.5

#performing the optimisation for the expression of density power divergence

obj1 <-rep(0,20)
obj2 <- rep(0,20)
obj3 <- rep(0,20)

g<-function(x)
{
  
  return(ifelse(x > 19,0,observed_prob[x+1]) )
  
  
}




Log_Density<- function(o)
{
  
  for(i in 0:100)
  {
    
    obj1[i+1] <- ( (exp(-o)*(o^i) )/factorial(i) )^r
    
    obj2[i+1] <- (g(i) * ( ( (exp(-o)*(o^i) )/factorial(i) )^(r-1) )  )
    
    #obj3[i+1] <- ( g(i)^r ) #No use in optimisation of parameter 
    
  }
  
  
  t<- log10( sum(obj1) ) +( r/(1-r) )*( log10(sum(obj2)) )# + (1/(r-1) )*( log10(sum(obj3)) )
  return(t)
}

optimize(Log_Density,c(0,10))


#Estimation of fitted values with best parameter


estimated_prob <- rep(0,19)


fitted_observation<-function(y)
{
  
  
  for(i in 0:19)
  {
    
    estimated_prob[i+1]<- ( exp(-y)*(y^i) ) /factorial(i)
    
  }
  return(estimated_prob)
  
}

Esti_values<-fitted_observation(0.95)*111 #since function returns probability

#Histogram comparing Robustness for all plots 

PD_alpha_1.2<-dpois(0:19, lambda = 0.73)*111
DPD_alpha_1.2<-dpois(0:19, lambda = 0.44)*111
LDPD_alpha_1.2<-dpois(0:19, lambda = 0.36)*111

final_Boxplot<-data.frame(PD_alpha_1.2,DPD_alpha_1.2,LDPD_alpha_1.2)

par(mfrow=c(1,3))
barplot(PD_alpha_1.2,ylim = c(0,80) ,xlim = c(0,20),main = "PD_alpha_1.2",ylab = 'Fitted_Values',xlab = 'No of cysts in kidney')
barplot(DPD_alpha_1.2,ylim = c(0,80) ,xlim = c(0,20) ,main = "DPD_alpha_1.2",ylab = 'Fitted_Values',xlab = 'No of cysts in kidney')
barplot(LDPD_alpha_1.2,ylim = c(0,80) ,xlim = c(0,20) ,main = "LDPD_alpha_1.2",ylab = 'Fitted_Values',xlab = 'No of cysts in kidney')




              # \\ For Normal model \\


#Newcombs Dataset
a<-c(28,26,33,24,34,-44,27,16,40,-2,29,22,24,21,25,30,23,29,31,19,24,20,36,32,36,28,25,21,28,29,37,25,28,26,30,32,36,26,30,22,36,23,27,27,28,27,31,27,26,
     
     33,26,32,32,24,39,28,24,25,32,25,29,27,28,29,16,23)

#Summarising the Dataset
table(a)
summary(a)

#1.Density Power Divergence Optimisation


#fixing the value of tuning parameter lambda
r<- 2

#performing the optimisation for the expression of density power divergence
value_data<-rep(0,66)
s<-function(u)
{
  t<- function(x)
  {
    
    
    
    
    ( (1/(sqrt(2*3.14*(u[2])^2 ) )   )*(exp(-(x-u[1])^2 /(2*(u[2])^2) ) ) )^(r)
    
    
    
    
  }
  
  total_integrate<- integrate(t,-Inf,Inf)$value
  
  for(i in 1:66)
  {
    
    value_data[i]<-(  (1/(sqrt(2*3.14*(u[2])^2 ) )   ) *(exp(-(a[i]-u[1])^2/(2*(u[2])^2)) ) )^(r-1)
    
    
  }
  
  total_term2<- (1+1/(r-1) )*(1/66)* sum(value_data)
  
  
  final<-total_integrate-  total_term2
  return(final)
  
}
optim(par=c(22,3), fn=s) 

#2.Log Density Power Divergence Optimisation


#fixing the value of tuning parameter lambda
r<- 2
s<-function(u)
{
  
  
  for(i in 1:66)
  {
    
    value_data[i]<-(  (1/(sqrt(2*3.14*(u[2])^2 ) )   ) *(exp(-(a[i]-u[1])^2/(2*(u[2])^2)) ) )^(r-1)
    
    
  }
  
  total_term2<- (r/(r-1) )*log10( ifelse(sum(value_data)/66 <= 0,1,sum(value_data)/66) )
  
  
  final<-log10( 1/( (( sqrt(2*3.14)*u[2] )^(r-1) )*sqrt(r) ) ) - total_term2
  return(final)
  
}
optim(par=c(22,3), fn=s) 


