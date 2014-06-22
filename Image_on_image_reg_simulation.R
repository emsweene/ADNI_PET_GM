temp <- commandArgs(TRUE)

n <- as.numeric(temp[1])
dim_roi <- c(as.numeric(temp[2]), as.numeric(temp[3]), as.numeric(temp[4]))

print(temp)

library("nplargedb")
library("nplargela")
library("lineprof")
library("dplyr")

setwd('/dexter/disk1/smart/AD/ADNI/Greven_ADNI')
source("nplarge-regression-utils.R")



make_image <- function(dims){
  periods <- runif(length(dims), .5, 9)
  dim_vecs <- lapply(dims, function(x) seq_len(x)/x)
  waves <- mapply(function(dim_vec, period) {
    runif(1, 1, 2) * sin(period*(dim_vec + rnorm(1, 0, .5)))
  }, dim_vecs, periods, SIMPLIFY=FALSE)
  Reduce(`%o%`, waves)
}


################################################################################
# simulate data from 
#   y_i(v) = beta0 + z_i*beta_z(v) +  x_i(v)*beta_x(v) + e_iv
# where v is a multi-index over some ROI
################################################################################

set.seed(19450508L)
n_dims <- length(dim_roi)

print(dim_roi)

dims <- lapply(dim_roi, function(x) seq_len(x))
names(dims) <- paste0("dim", 1:n_dims)

snr <- 10

data <- data.frame(subject=factor(1:n),
                   z = scale(rnorm(n)))

# centered images:
x_raw <- replicate(n, make_image(dim_roi))
x_mean <- apply(x_raw, 1:n_dims, sum)/n 
x <- array(apply(x_raw, n_dims+1, function(x_slice) x_slice - x_mean), 
           dim=c(dim_roi, n)) ## stupid apply collapses first dims wtf...

x_cube <- tbl_cube(dimensions=append(dims, list(subject=1:n)), 
                   measures = list(x=x))
beta_x <- make_image(dim_roi) + make_image(dim_roi)

# eta_x = x * beta
eta_x <-  x_cube$mets$x * array(beta_x, dim=c(dim_roi, n))

beta_z <- make_image(dim_roi) + make_image(dim_roi)
z_array <- aperm(array(data$z, dim=c(n, dim_roi)), c(2:(n_dims+1), 1))
eta_z <-  z_array * array(beta_z, dim=c(dim_roi, n))

var_eta <- var(as.vector(eta_x + eta_z))
sd_e <- sqrt(var_eta/snr)

beta0 <- .1
eta <- beta0 + eta_x + eta_z
y <- eta + array(rnorm(prod(dim_roi) * n, sd=sd_e), dim=c(dim_roi, n))

y_cube <- tbl_cube(dimensions=append(dims, list(subject=1:n)), 
                   measures = list(y=y, eta=eta))


################################################################################
# import data, set up modelterms:
################################################################################                   

# import data                   
denv <- dataenv(data, features="z", dims="subject")
add_features(data=x_cube, dataenv=denv)
add_dims(data=x_cube, dataenv=denv)

# set up terms & pred
const_subj <- const(n, denv)
#const_dims <- lapply(dim_roi, function(dim) const(dim, denv))
#const_v <- do.call(kron, const_dims)
# slightly more efficent than the above:
const_v <- const(prod(dim_roi), denv)

bs_dims <- lapply(paste0("dim", 1:n_dims), 
                  function(dim) bspline(dim, denv, df=6, constraint="none"))
bs_v <- do.call(kron, bs_dims)
intercept <- kron(const_v, const_subj)
bs_z <- kron(bs_v, const(n, denv, by=data$z))
bs_x <- kron(bs_v, const_subj, by=as.vector(x_cube$mets$x))

# constant intercept, z_i*beta_z(v), x_i(v)*beta_x(v)
(pred <- concat(intercept, bs_z, bs_x))
resp <- as.vector(y_cube$mets$y)


################################################################################
# unpenalized fit:
################################################################################

mem.it.up <- lineprof(time.it.up <- system.time(m <- nplarge_lm_fit(pred, resp)))


print('boo')

################################################################################
# penalized fit:
################################################################################

pred$terms[[2]]$pen <- do.call(make_kron_sum,
lapply(pred$terms[[2]]$terms[1:n_dims],
make_diffpen))
pred$terms[[3]]$pen <- do.call(make_kron_sum,
lapply(pred$terms[[3]]$terms[1:n_dims],
make_diffpen))

print('the')
mem.it.BFGS <- lineprof(time.it.BFGS <- system.time(test_opt <- nplarge_pen_fit(x=pred, y=resp,
                                        optimizer="optim-BFGS",
                                        optimizer.control=list(trace=3))))

print('you')
mem.it.bobyqa <- lineprof(time.it.bobyqa <- system.time(test_minqa <- nplarge_pen_fit(x=pred, y=resp, 
                                          optimizer="bobyqa",
                                          optimizer.control=list(iprint=2, maxfun=1e3))))

print('hello')
mem.it.nlminb <- lineprof(mem.it.nlminb <- system.time(test_nlmin <- nplarge_pen_fit(x=pred, y=resp, 
                                          optimizer="nlminb",
                                          optimizer.control=list(trace=1, eval.max=1e3)))) 

print('hi')
setwd('/dexter/disk1/smart/AD/ADNI/Greven_ADNI/Time_Mem')
save(mem.it.up, time.it.up,
     mem.it.BFGS, time.it.BFGS,
     mem.it.bobyqa, time.it.bobyqa,
     mem.it.nlminb, mem.it.nlminb,
     file = paste('output', dim_roi[1], dim_roi[2], dim_roi[3], n, sep = '_') )
