library(readxl)
library(ggplot2)

library(minpack.lm)
library(rootSolve)
library(coda)


# individual samples

df<-read_excel("appendix.xlsx", sheet="Sheet2",
               range="B3:O11",col_names = FALSE)
df<-matrix(unlist(df),ncol=14,byrow = FALSE)

maxv<-read_excel("appendix.xlsx", sheet= "Sheet2",
                  range= "B12:O12", col_names = FALSE)
maxv<-c(t(maxv))


dosepts<-c(t(read_excel("appendix.xlsx", sheet= "Sheet2",
                        range = "A3:A11",col_names = FALSE)))
x<-c(dosepts)
xpts=10^seq(-4,1.2,by=1e-2)
tend=6
mu=log(2)/2.5
x0=mu

# Fit corresponding parameters
control = list(maxiter =1000,maxfev=1000)

# all parameters together
start_values0a=c(n=1/2, m=2, g_max = 0.4, b_M=1.4, g_M=6)

for (ii in seq(6,13,by=1)) {
  head1 <- paste0("dataset_",ii,".png")
  headt <- paste0("dataset ",ii)
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  y<-df[,ii]; y<-abs(y)
  
  fit<-nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m/(g_M + x^m)+mu)) 
             + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m/(g_M + x^m)),
             start=start_values0a,
             data=data.frame(x,y), control=control)
  summary(fit)
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new1agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new1agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new1agg.dat", sep = " ", 
              dec = ".", 
              col.names = !file.exists("2026_myDF_new1agg.dat"), 
              row.names = FALSE, append = TRUE)
 
    # compute AIC
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new1agg.dat", sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new1agg.dat"), 
              row.names = FALSE, append = TRUE)
  
  png(file=head1,width=500,height=500,units="px")
  plot(x[2:9],y[2:9],log='x',xlim=c(1e-4,12.5),main=headt,
       xlab="concentration",ylab="response")
    lines(xpts,  (-b_max*xpts^out[1]/(out[4]+xpts^out[1]))*exp(-tend*(out[3]*xpts^out[2]/(out[5] + xpts^out[2])+mu)) 
        + (1+b_max*xpts^out[1]/(out[4]+xpts^out[1]))*exp(-tend*out[3]*xpts^out[2]/(out[5] + xpts^out[2]))
        ,lty=2)
  
  dev.off()
  
}


# with fixed m

start_values1=c(n=1/2, g_max = 0.4,b_M=1.4,g_M=6)

m=2
for (ii in seq(1,13,by=1)) {
  head1 <- paste0("dataset_",ii,".png")
  headt <- paste0("dataset ",ii)
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii]; y<-abs(y)
  
  fit<-nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m/(g_M + x^m)+mu)) 
             + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m/(g_M + x^m)),
             start=start_values1,
             data=data.frame(x,y), control=control)
  summary(fit)
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new2agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new2agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new2agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new2agg.dat"), 
              row.names = FALSE, append = TRUE)
  # compute AIC
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new2agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new2agg.dat"), 
              row.names = FALSE, append = TRUE)
}

m=3
for (ii in seq(1,13,by=1)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii]; y<-abs(y)
  
  fit<-nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m/(g_M + x^m)+mu)) 
             + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m/(g_M + x^m)),
             start=start_values1,
             data=data.frame(x,y), control=control)
  summary(fit)
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new3agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new3agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new3agg.dat", sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new3agg.dat"), 
              row.names = FALSE, append = TRUE)
  # compute AIC
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new3agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new3agg.dat"), 
              row.names = FALSE, append = TRUE)
  
}

m=4
for (ii in seq(11,13,by=1)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii]; y<-abs(y)
  
  fit<-nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m/(g_M + x^m)+mu)) 
             + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m/(g_M + x^m)),
             start=start_values1,
             data=data.frame(x,y), control=control)
  summary(fit)
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new3agg2.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new3agg2.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new3agg2.dat", sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new3agg2.dat"), 
              row.names = FALSE, append = TRUE)
  # compute AIC
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new3agg2.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new3agg2.dat"), 
              row.names = FALSE, append = TRUE)
  
}


# without g_M
# all parameters
start_values2=c(n=1/2, m=2.5,g_max = 0.4,b_M=1.4)

for (ii in seq(1,13,by=1)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii];y<-abs(y)
  
  fit <- nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m+mu)) 
               + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m), 
               start=start_values2, data=data.frame(x,y),control=control)  
  
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new4agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new4agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new4agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new4agg.dat"), 
              row.names = FALSE, append = TRUE)
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new4agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new4agg.dat"), 
              row.names = FALSE, append = TRUE)
}


# with fixed exponent
m=2
start_values2a=c(n=1/2, g_max = 0.04,b_M=1.4)

for (ii in seq(1,13,by=1)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii];y<-abs(y)
  
  if (ii==7){ x<-c(dosepts[1:5],dosepts[7:9]);
  y<-c(df[1:5,ii],df[7:9,ii]); y<-abs(y)}
  
  if (ii==11){ x<-c(dosepts[1:6],dosepts[8:9]);
  y<-c(df[1:6,ii],df[8:9,ii]); y<-abs(y)}
  
  fit <- nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m+mu)) 
               + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m), 
               start=start_values2a, data=data.frame(x,y),control=control)  
  
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new5agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new5agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new5agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new5agg.dat"), 
              row.names = FALSE, append = TRUE)
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new5agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new5agg.dat"), 
              row.names = FALSE, append = TRUE)
}

m=3

for (ii in seq(1,13,by=1)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii];y<-abs(y)
  
  if (ii==7){ x<-c(dosepts[1:5],dosepts[7:9]);
  y<-c(df[1:5,ii],df[7:9,ii]); y<-abs(y)}
  
  if (ii==11){ x<-c(dosepts[1:6],dosepts[8:9]);
  y<-c(df[1:6,ii],df[8:9,ii]); y<-abs(y)}
  
  fit <- nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m+mu)) 
               + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m), 
               start=start_values2a, data=data.frame(x,y),control=control)  
  
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new6agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new6agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new6agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new6agg.dat"), 
              row.names = FALSE, append = TRUE)
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new6agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new6agg.dat"), 
              row.names = FALSE, append = TRUE)
}

m=4

for (ii in seq(1,13,by=1)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  x<-dosepts
  y<-df[,ii]; y<-abs(y)
  
  if (ii==7){ x<-c(dosepts[1:5],dosepts[7:9]);
  y<-c(df[1:5,ii],df[7:9,ii]); y<-abs(y)}
  
  if (ii==11){ x<-c(dosepts[1:6],dosepts[8:9]);
  y<-c(df[1:6,ii],df[8:9,ii]); y<-abs(y)}
  
  fit <- nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m+mu)) 
               + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m), 
               start=start_values2a, data=data.frame(x,y),control=control)  
  
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new7agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new7agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new7agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new7agg.dat"), 
              row.names = FALSE, append = TRUE)
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new7agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new7agg.dat"), 
              row.names = FALSE, append = TRUE)
}


m=4

for (ii in c(7,11)) {
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  if (ii==7){ x<-c(dosepts[1:5],dosepts[7:9]);
    y<-c(df[1:5,ii],df[7:9,ii]); y<-abs(y)}
  
  if (ii==11){ x<-c(dosepts[1:6],dosepts[8:9]);
    y<-c(df[1:6,ii],df[8:9,ii]); y<-abs(y)}
  
  fit <- nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*(g_max*x^m+mu)) 
               + (1+b_max*x^n/(b_M+x^n))*exp(-tend*g_max*x^m), 
               start=start_values2a, data=data.frame(x,y),control=control)  
  
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new8agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new8agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new8agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new8agg.dat"), 
              row.names = FALSE, append = TRUE)
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new8agg.dat", 
              sep = " ", dec = ".",
              col.names = !file.exists("2026_myDF_new8agg.dat"), 
              row.names = FALSE, append = TRUE)
}


# fit monotone data sets
# all parameters together
start_values0=c(n=1/2,b_M=1.4)

for (ii in c(3,14) ) {
  head1 <- paste0("dataset_",ii,".png")
  headt <- paste0("dataset ",ii)
  b_max<-maxv[ii]
  b_max<-(b_max-1)/(1-exp(-mu*tend))
  
  y<-df[,ii]; y<-abs(y)
  
  fit<-nlsLM(y~(-b_max*x^n/(b_M+x^n))*exp(-tend*mu) 
             + (1+b_max*x^n/(b_M+x^n)),
             start=start_values0,
             data=data.frame(x,y), control=control)
  summary(fit)
  out<-summary(fit)$coefficients
  write.table(t(c(ii,mu,b_max)),"2026_myDF_new0agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new0agg.dat"), 
              row.names = FALSE, append = TRUE)
  write.table(out, "2026_myDF_new0agg.dat", sep = " ", 
              dec = ".", 
              col.names = !file.exists("2026_myDF_new0agg.dat"), 
              row.names = FALSE, append = TRUE)
  # compute AIC
  out1<-AIC(fit)
  write.table(out1, "2026_myDF_new0agg.dat", 
              sep = " ", dec = ".", 
              col.names = !file.exists("2026_myDF_new0agg.dat"), 
              row.names = FALSE, append = TRUE)
  png(file=head1,width=500,height=500,units="px")
  plot(x[2:9],y[2:9],log='x',xlim=c(1e-4,12.5),main=headt,
       xlab="concentration",ylab="response")
  lines(xpts,  (-b_max*xpts^out[1]/(out[2]+xpts^out[1]))*exp(-tend*mu) 
        + (1+b_max*xpts^out[1]/(out[2]+xpts^out[1]))
        ,lty=2)
  
  dev.off()
  
}
