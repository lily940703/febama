# A simple example of feature-based bayesian model averaging
* model 1  AR(1)  
* model 2  random walk 

```
x <- generate_ts(n.ts = 1, freq = 1, nComp = 2, n = 20)
y<-x$N1$x
#加初值
y<-c(0,y)
y<-c(0.00000,21.72881,26.80634,28.72464,33.37096,34.12180,40.02512,
43.12906,42.79923,43.28925,47.78483,50.43229,52.01138,56.37153,59.94831,
64.44778,64.94955,69.85086,72.23654,78.41327,80.10752)
```
