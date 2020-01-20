#Figure Log

set.seed(123)
simdata = sim.settlement(K1=1000,K2=1000)
nsim = 100
b = c(0,0.3,0.7)
r = c(0.1,0.3,0.7)
params=expand.grid(r=r,b=b,nsim=1:nsim)
params$pr=NA

for (i in 1:nrow(params))
{
  params$pr[i]=biasedsampling(simdata,r=params$r[i],b=params$b[i])
}

# Plot results
params$r=as.factor(params$r)
p <- ggboxplot(params, x = "b", y = "pr",
               color = "r", palette = "jco",
               add = "jitter",ggtheme=theme_grey(),alpha=0.9) + ylab("Percentage Change") + geom_hline(yintercept = 0, linetype="dashed") +theme(legend.position = 'top')

