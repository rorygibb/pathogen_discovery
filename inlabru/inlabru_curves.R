

setwd("C:/Users/roryj/Documents/PhD/202008_pathogendiscovery/code/pathogen_discovery/inlabru/")
library(inlabru)
library(INLA)
library(mgcv)

### fitting GLM
load("data/awards.RData")
head(awards)

# specify inlabru
cmp1 <- num_awards ~ math + Intercept
fit.glm.bru <- bru(cmp1, family = "poisson", data = awards)

# with individual level random effect
cmp3 <- num_awards ~ math + Intercept + rand.eff(map = 1:200, model = "iid", n = 200)
fit3.glmm.bru <- bru(cmp3, family = "poisson", data = awards )

# with PC prior specified
cmp5 <- num_awards ~ math + Intercept + rand.eff(map = 1:200, model = "iid", n = 200,
                                                 hyper=list(prec=list(param=c(10,0.1),
                                                                      prior="pc.prec")))
fit5.glmm.bru <- bru(cmp5, family = "poisson", data = awards )
summary(fit5.glmm.bru)
INLA:::summary.inla(fit5.glmm.bru)



### fitting to georeferenced data

init.tutorial()
colsc <- function(...) {
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(11,"RdYlBu")),
                       limits = range(..., na.rm=TRUE))
}

# load data
data(Poisson2_1D)
cd <- countdata2

# plot
ggplot(cd) + geom_point(aes(x, y = count)) + ylim(0, max(cd$count))

# fit 1dimensional SPDE to counts
x <- seq(0, 55, length = 50) # this sets mesh points - try others if you like
mesh1D <- inla.mesh.1d(x, boundary = "free")
ggplot() + gg(mesh1D) + xlim(0,55)

# specify spde
spde <- inla.spde2.pcmatern(mesh1D,
                            prior.range=c(1, 0.01),
                            prior.sigma=c(10, 0.01))

# fit to counts (effect is called "field" and we can call it what we want)
comp <- count ~ field(map = x, model = spde) + Intercept

# fit with E as exposure variable (i.e. offset)
fit2.bru <- bru(comp, cd, family = "poisson", options = list(E = cd$exposure))
summary(fit2.bru)

# predict at x locs#
# n.b expected true counts are stored in E_nc2
x4pred <- data.frame(x = x)
pred2.bru <- predict(fit2.bru, x4pred, x ~ exp(field+Intercept), n.samples=1000)
true.lambda <- data.frame(x = cd$x, y = E_nc2/cd$exposure)
ggplot() + 
  gg(pred2.bru) + 
  geom_point(data = cd, aes(x = x, y = count/cd$exposure), cex = 2) + 
  geom_point(data = true.lambda, aes(x, y), pch="_", cex = 9, col = "blue") +
  xlab("x") + ylab("Intensity")


### 1 dimensional point process data (i.e. log gaussian cox)

data(Poisson2_1D)

# locations along the x axis where events occur (i.e. not discretised)
head(pts2)
ggplot(pts2) + 
  geom_histogram(aes(x = x), binwidth = 55/20, boundary = 0, fill = NA, color = "black") +
  geom_point(aes(x), y = 0, pch = "|", cex = 4) + 
  coord_fixed(ratio = 1)

# create mesh - sets domain of model
# fit 1dimensional SPDE to counts
x <- seq(0, 55, length = 50) # this sets mesh points - try others if you like
mesh1D <- inla.mesh.1d(x, boundary = "free")
spde <- inla.spde2.pcmatern(mesh1D,
                            prior.range=c(1, 0.01),
                            prior.sigma=c(10, 0.01))

# model (here to the xs)
mdl <- x ~ spde1D(map = x, model = spde) + Intercept

# fit model
fit.spde <- lgcp(mdl, pts2, ips = ipoints(c(0,55), 50, name = "x"))

# examine posterior range
post.range = spde.posterior(fit.spde,name="spde1D",what="range"); plot(post.range)

# predict on exp intensity scale (note specifying the exponential)
predf <- data.frame(x = seq(0, 55, by = 1)) # Set up a data frame of explanatory values at which to predict
pred_spde <- predict(fit.spde, predf,  ~ exp(spde1D + Intercept))
plot(pred_spde, color = "red") + 
  geom_point(data = pts2, aes(x = x), y = 0, pch = "|", cex = 2) + 
  xlab("x") + ylab("Intensity")
