##################################################################################
## Read in required libraries 
## For cluster: qrsh -l  mem_free=19G,h_vmem=20G
##################################################################################

library(oro.nifti)
library(RColorBrewer)
library(Matrix)
library(nplargedb)
library(nplargela)
library(dplyr)
library(lineprof)

##################################################################################
## Set working directory, load nplarge regression utilties 
##################################################################################

setwd('/dexter/disk1/smart/AD/ADNI/Greven_ADNI')
source("nplarge-regression-utils.R")

##################################################################################
## Get the subject ids 
##################################################################################

fileDir <- '/dexter/disk1/smart/AD/ADNI/Greven_ADNI'
files <- dir(fileDir, pattern = "*.nii.gz", full.names = TRUE)
subject.ids <- unlist(lapply(files, function(x) tail(unlist(strsplit(x, '/')), n=1)))
subject.ids <- unique(unlist(lapply(subject.ids, function(x) substr(x, 1, 10))))

##################################################################################
## Load the meta data for the subjects, create binary indicator for MCI and gender
## Assign 1 for MCI and females
##################################################################################

meta.data <- read.table('Patient_metadata_3yrs_with_fw-ups.csv', header = TRUE, sep = ',')
baseline.meta.data <- meta.data[meta.data$VISCODE == 'bl',]
baseline.subjects <- merge(data.frame(Subject = c(subject.ids)), baseline.meta.data)
baseline.subjects$MCI <- rep(0, length(baseline.subjects$DXCHANGE))
baseline.subjects$GENDER <- rep(0, length(baseline.subjects$PTGENDER))
baseline.subjects$MCI[baseline.subjects$DXCHANGE == 'Stable: MCI'] <- 1
baseline.subjects$GENDER[baseline.subjects$PTGENDER == '2'] <- 1

##################################################################################
## Calculate PET differenced from 1st and 2nd and 3rd and 4th visit (both have a 
## year in between)
##################################################################################


for(i in 1:length(subject.ids)){
  
  pet1 <- readNIfTI(paste0(subject.ids[i], '_PET_coreg_to_MRI_MNI_1.nii.gz'), reorient = FALSE)[,,]
  pet2 <- readNIfTI(paste0(subject.ids[i], '_PET_coreg_to_MRI_MNI_2.nii.gz'), reorient = FALSE)[,,]
  pet.diff.1 <- pet2 - pet1 
  
  pet3 <- readNIfTI(paste0(subject.ids[i], '_PET_coreg_to_MRI_MNI_3.nii.gz'), reorient = FALSE)[,,]
  pet4 <- readNIfTI(paste0(subject.ids[i], '_PET_coreg_to_MRI_MNI_4.nii.gz'), reorient = FALSE)[,,]
  pet.diff.2 <- pet4 - pet3
  
  writeNIfTI(pet.diff.1, paste0(subject.ids[i], '_pet_diff_2_1'))
  writeNIfTI(pet.diff.2, paste0(subject.ids[i], '_pet_diff_4_3'))
  
  print(i)
}

##################################################################################
## Load in the ROIs
##################################################################################

ROI.large <- readNIfTI('Precun_cube.nii', reorient = FALSE)[,,]
ROI.small <- readNIfTI('Precun_cube_small.nii', reorient = FALSE)[,,]

##################################################################################
## Do some EDA in the large ROI
##################################################################################

mean.grey.matter.1 <- c() 
mean.grey.matter.2 <- c() 
mean.grey.matter.3 <- c() 
mean.grey.matter.4 <- c() 
mean.pet.slopes.2.1 <- c()
mean.pet.slopes.4.3 <- c()  
grey.matter.1 <- c() 
grey.matter.2 <- c() 
grey.matter.3 <- c() 
grey.matter.4 <- c() 
pet.slopes.2.1 <- c()
pet.slopes.4.3 <- c()  
mci.status <- c()
age <- c()
edu.status <- c()
gender <- c()

ROI <- ROI.large 

for(i in 1:length(subject.ids)){
  PET21 <- readNIfTI( paste0(subject.ids[i], '_pet_diff_2_1'), reorient = FALSE)[,,]
  PET43 <- readNIfTI( paste0(subject.ids[i], '_pet_diff_4_3'), reorient = FALSE)[,,]  
  GM1 <- readNIfTI(paste0(subject.ids[i], '_MRI_segm_MNI_1.nii.gz'), reorient = FALSE)[,,]
  GM2 <- readNIfTI(paste0(subject.ids[i], '_MRI_segm_MNI_2.nii.gz'), reorient = FALSE)[,,]
  GM3 <- readNIfTI(paste0(subject.ids[i], '_MRI_segm_MNI_3.nii.gz'), reorient = FALSE)[,,]
  GM4 <- readNIfTI(paste0(subject.ids[i], '_MRI_segm_MNI_4.nii.gz'), reorient = FALSE)[,,]
  grey.matter.1 <- c(grey.matter.1, GM1[ROI == 1])
  mean.grey.matter.1 <- c(mean.grey.matter.1, mean(GM1[ROI == 1]))
  grey.matter.2 <- c(grey.matter.2, GM2[ROI == 1]) 
  mean.grey.matter.2 <- c(mean.grey.matter.2, mean(GM2[ROI == 1]))
  grey.matter.3 <- c(grey.matter.3, GM3[ROI == 1]) 
  mean.grey.matter.3 <- c(mean.grey.matter.3, mean(GM3[ROI == 1]))
  grey.matter.4 <- c(grey.matter.4, GM4[ROI == 1]) 
  mean.grey.matter.4 <- c(mean.grey.matter.4, mean(GM4[ROI == 1]))
  pet.slopes.2.1 <- c(pet.slopes.2.1, PET21[ROI == 1])
  mean.pet.slopes.2.1<- c(mean.pet.slopes.2.1, mean(PET21[ROI == 1]))  
  pet.slopes.4.3 <- c(pet.slopes.4.3, PET43[ROI == 1])  
  mean.pet.slopes.4.3<- c(mean.pet.slopes.4.3, mean(PET43[ROI == 1]))    
  mci.status <- c(mci.status, rep(baseline.subjects$MCI[i], sum(ROI)))
  age <- c(age, rep(baseline.subjects$AGE[i], sum(ROI)))
  edu.status <- c(edu.status, rep(baseline.subjects$PTEDUCAT[i], sum(ROI)))
  gender <- c(gender, rep(baseline.subjects$GENDER[i], sum(ROI)))
  print(i) 
}

large.roi <- data.frame(grey.matter.1, 
                        grey.matter.2,
                        grey.matter.3,
                        grey.matter.4, 
                        pet.slopes.2.1,
                        pet.slopes.4.3,  
                        mean.grey.matter.1,
                        mean.grey.matter.2,
                        mean.grey.matter.3, 
                        mean.grey.matter.4,
                        mean.pet.slopes.2.1,
                        mean.pet.slopes.4.3,  
                        mci.status,
                        age,
                        edu.status,
                        gender)

save(large.roi, file = 'largeROIexplore.rda')



##################################################################################
## Pairs plot with all data 
##################################################################################

downsample <- sample(1: dim(large.roi)[1], 5000, replace = F)
pdf('pairs_plot.pdf')
pairs( ~ grey.matter.1 +
      grey.matter.2 +
      grey.matter.3 +
      grey.matter.4 +
      pet.slopes.2.1 +
      pet.slopes.4.3 + 
      mci.status +
      age +
      edu.status +
      gender, data = large.roi[downsample,])
dev.off()


##################################################################################
## Do some EDA in the large ROI
##################################################################################

pdf('LargeROI_BaseGrey_PET43.pdf')
smooth.mci <- loess.smooth(grey.matter.1[mci.status == 1], pet.slopes.4.3[mci.status == 1],
                           data = large.roi) 
smooth.no.mci <- loess.smooth(grey.matter.1[mci.status == 0], pet.slopes.4.3[mci.status == 0], 
                              data = large.roi) 
plot(pet.slopes.4.3 ~ grey.matter.1, data = large.roi[downsample,],  pch = 20, 
     main = 'Large ROI', ylab = 'PET Difference visit 4 and 3', 
     xlab = c('Baseline Grey Matter'))
lines(smooth.mci, lwd = 3, col = 'red')
lines(smooth.no.mci, lwd = 3, col = 'blue')
legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), lwd = 3)
dev.off()


pdf('LargeROI_Grey2_PET43.pdf')
smooth.mci <- loess.smooth(grey.matter.2[mci.status == 1], pet.slopes.4.3[mci.status == 1],
                           data = large.roi) 
smooth.no.mci <- loess.smooth(grey.matter.2[mci.status == 0], pet.slopes.4.3[mci.status == 0], 
                              data = large.roi) 
plot(pet.slopes.4.3 ~ grey.matter.2, data = large.roi[downsample,],  pch = 20, 
     main = 'Large ROI', ylab = 'PET Difference visit 4 and 3', 
     xlab = c('Grey Matter Visit 2'))
lines(smooth.mci, lwd = 3, col = 'red')
lines(smooth.no.mci, lwd = 3, col = 'blue')
legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), lwd = 3)
dev.off()



pdf('LargeROI_Grey4_PET21.pdf')
smooth.mci <- loess.smooth(grey.matter.4[mci.status == 1], pet.slopes.2.1[mci.status == 1],
                           data = large.roi) 
smooth.no.mci <- loess.smooth(grey.matter.4[mci.status == 0], pet.slopes.2.1[mci.status == 0], 
                              data = large.roi) 
plot(pet.slopes.2.1 ~ grey.matter.4, data = large.roi[downsample,],  pch = 20, 
     main = 'Large ROI', ylab = 'PET Difference visit 2 and 1', 
     xlab = c('Grey Matter Visit 4'))
lines(smooth.mci, lwd = 3, col = 'red')
lines(smooth.no.mci, lwd = 3, col = 'blue')
legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), lwd = 3)
dev.off()



##################################################################################
## Plot mean pet slopes and mean gm -- see if we see any pattern in the region 
##################################################################################

pdf('mean_PET21_GM4.pdf')
use.colors <- c("blue", "red")
plot(mean.grey.matter.4, mean.pet.slopes.2.1, col = use.colors[c(baseline.subjects$MCI +1 )],
     pch = 20, ylab = 'PET Difference visit 2 and 1', 
     xlab = c('Grey Matter Visit 4'), main = 'Large ROI')
    legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), pch = 20, lwd = 3)
dev.off()

pdf('mean_PET21_GM3.pdf')
use.colors <- c("blue", "red")
plot(mean.grey.matter.3, mean.pet.slopes.2.1, col = use.colors[c(baseline.subjects$MCI +1 )],
     pch = 20, ylab = 'PET Difference visit 2 and 1', 
     xlab = c('Grey Matter Visit 3'), main = 'Large ROI')
legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), pch = 20, lwd = 3)
dev.off()

pdf('mean_PET34_GM1.pdf')
use.colors <- c("blue", "red")
plot(mean.grey.matter.1, mean.pet.slopes.4.3, col = use.colors[c(baseline.subjects$MCI +1 )],
     pch = 20, ylab = 'PET Difference visit 4 and 3', 
     xlab = c('Grey Matter Visit 1'), main = 'Large ROI')
legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), pch = 20, lwd = 3)
dev.off()

pdf('mean_PET34_GM2.pdf')
use.colors <- c("blue", "red")
plot(mean.grey.matter.2, mean.pet.slopes.4.3, col = use.colors[c(baseline.subjects$MCI +1 )],
     pch = 20, ylab = 'PET Difference visit 4 and 3', 
     xlab = c('Grey Matter Visit 2'), main = 'Large ROI')
legend("topleft",  c("MCI", "Controls"), col = c('red', 'blue'), pch = 20, lwd = 3)
dev.off()

fit.lm.grey4.petslope21 <- lm(mean.grey.matter.4 ~ mean.pet.slopes.2.1 + baseline.subjects$MCI + 
               baseline.subjects$GENDER + baseline.subjects$PTEDUCAT +
               baseline.subjects$AGE
               )

fit.lm.grey3.petslope21 <- lm(mean.grey.matter.3 ~ mean.pet.slopes.2.1 + baseline.subjects$MCI + 
                                baseline.subjects$GENDER + baseline.subjects$PTEDUCAT +
                                baseline.subjects$AGE
                                )


fit.lm.petslope43.grey1 <- lm( mean.pet.slopes.4.3 ~ mean.grey.matter.1 + baseline.subjects$MCI + 
                                baseline.subjects$GENDER + baseline.subjects$PTEDUCAT +
                                baseline.subjects$AGE
                                )

fit.lm.petslope43.grey2 <- lm( mean.pet.slopes.4.3 ~ mean.grey.matter.2 + baseline.subjects$MCI + 
                                 baseline.subjects$GENDER + baseline.subjects$PTEDUCAT +
                                 baseline.subjects$AGE
                                )

##################################################################################
## Preparing this data for use with the nplarge packages
##################################################################################

load('largeROIexplore.rda')

n <-  94
dim_roi <- c(7, 7, 7)
n_dims <- length(dim_roi) 
dims <- lapply(dim_roi, function(x) seq_len(x))
names(dims) <- paste0("dim", 1:n_dims)

transform.grey.matter <- large.roi$grey.matter.4
transform.grey.matter[transform.grey.matter  <= 0.001] <- .001
transform.grey.matter[transform.grey.matter  >= 0.999] <- 0.999
transform.grey.matter <- log(transform.grey.matter /( 1 - transform.grey.matter))

y <- array(transform.grey.matter, c(dim_roi, n))


y_cube <- tbl_cube(dimensions=append(dims, list(Subject=subject.ids)), 
                   measures = list(y=y))

x <- array(large.roi$pet.slopes.2.1, c(dim_roi, n))

x_cube <- tbl_cube(dimensions=append(dims, list(Subject=subject.ids)), 
                   measures = list(x=x))

keep <- c("AGE", "PTGENDER", "MCI", "PTEDUCAT", "Subject")
data <- baseline.subjects[keep]

denv <- dataenv(data,  dims="Subject")
add_features(data=x_cube, dataenv=denv)
add_dims(data=x_cube, dataenv=denv)


##Creat a compressed data representaion of data cube tbl
## denv <- dataenv(cov_cube)

## Create intercept term for subjects
const_subj <- const(n, denv)

##Create intercept term for the dimmensions of the ROI
const_v <- const(prod(dim_roi), denv)

##Create bspline basis
bs_dims <- lapply(paste0("dim", 1:n_dims), 
                  function(dim) bspline(dim, denv, df=6, constraint="none"))
bs_v <- do.call(kron, bs_dims)

##Create intercept array 
intercept <- kron(const_v, const_subj)

##create bspline basis for age, mci, education, and gender
bs_age <- kron(bs_v, const(n, denv, by=baseline.subjects$AGE))
bs_mci <- kron(bs_v, const(n, denv, by=baseline.subjects$MCI))
bs_edu <- kron(bs_v, const(n, denv, by=baseline.subjects$PTEDUCAT))
bs_gen <- kron(bs_v, const(n, denv, by=baseline.subjects$PTGENDER))

##create bspline basis for grey matter
bs_x <- kron(bs_v, const_subj, by=as.vector(x_cube$mets$x))

(pred <- concat(intercept, bs_x, bs_age, bs_mci, bs_edu, bs_gen))
resp <- as.vector(y_cube$mets$y)


##################################################################################
## Fit Unpenalized Model 
##################################################################################

system.time(m <- nplarge_lm_fit(pred, resp))

##################################################################################
## Visualize the Results (Analysis of Resdiuals)
##################################################################################

y_hat <- array(m$fitted.values, c(dim_roi, n))
my.colors <- colors()[sample(2:length(colors()), 94, replace = FALSE)]
my.colors <- c(apply(as.matrix(my.colors), 1, function(x) rep(x, prod(dim_roi))))
age <- c(apply(as.matrix(data$AGE), 1, function(x) rep(x, prod(dim_roi))))
edu <- c(apply(as.matrix(data$PTEDUCAT), 1, function(x) rep(x, prod(dim_roi))))

downsample <- sample(1: length( m$residuals), 5000, replace = FALSE)
pdf('grey_fit_residuals.pdf')
plot(m$fitted.values[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residuals', xlab = 'Fitted Transformed Gray Matter', 
     main = 'Residuals Colored by Subject')
dev.off()

pdf('age_vs_residuals.pdf')
plot(age[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'Age', 
     main = 'Residuals Colored by Subject')
dev.off()

pdf('edu_vs_residuals.pdf')
plot(edu[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'Education', 
     main = 'Residuals Colored by Subject')
dev.off()

pdf('pet_vs_residuals.pdf')
plot(large.roi$pet.slopes.2.1[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'PET Slopes', 
     main = 'Residuals Colored by Subject')
dev.off()

MSE <- mean(m$residuals^2)
print(MSE)
##  6.600595

##################################################################################
## Residuals versus other covariates not included in the model (mainly genetic 
## covariates)
##################################################################################

APGEN1 <- c(apply(as.matrix(baseline.subjects$APGEN1), 1, function(x) rep(x, prod(dim_roi))))
APGEN2 <- c(apply(as.matrix(baseline.subjects$APGEN2), 1, function(x) rep(x, prod(dim_roi))))
PTAU181P <- c(apply(as.matrix(baseline.subjects$PTAU181P), 1, function(x) rep(x, prod(dim_roi))))
Ab.status <- c(apply(as.matrix(baseline.subjects$Ab.status), 1, function(x) rep(x, prod(dim_roi))))


pdf('APGEN1_vs_residuals.pdf')
plot(APGEN1[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'APGEN1', 
     main = 'Residuals Colored by Subject')
dev.off()

pdf('APGEN2_vs_residuals.pdf')
plot(APGEN2[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'APGEN1', 
     main = 'Residuals Colored by Subject')
dev.off()

pdf('PTAU181P_vs_residuals.pdf')
plot(PTAU181P[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'PTAU181P', 
     main = 'Residuals Colored by Subject')
dev.off()

pdf('Abstatus_vs_residuals.pdf')
plot(Ab.status[downsample], m$residuals[downsample], 
     col = my.colors[downsample], ylab = 'Residual', xlab = 'PTAU181P', 
     main = 'Residuals Colored by Subject')
dev.off()


set.seed(33)

subjects <- sample(1:n, 10)
slices <- sample(1:dim_roi[3], 10, replace = TRUE)
clrs <- heat.colors(24)
brks <- seq(min(resp), max(resp), l=25) 

pdf('GMOUT_LARGEROI_scalar_reg_est.pdf')
for(i in 1:10){
  par(mfrow = c(1,2))
  image(y_cube$mets$y[, , slices[i], subjects[i]], col=clrs, breaks=brks,
        main="Grey Matter Probabilities")
  image(y_hat[, , slices[i], subjects[i]], col=clrs, breaks=brks, 
        main="Estimate")    
}
dev.off()

##################################################################################
## Visualize the Results (Coefficients)
##################################################################################


m$coef[1]
coef_pet <- m$coef[2 : cumsum(pred$c)[2]]
coef_age <- m$coef[c(cumsum(pred$c)[2] + 1) : cumsum(pred$c)[3]]
coef_mci <- m$coef[c(cumsum(pred$c)[3] + 1) : cumsum(pred$c)[4]]
coef_edu <- m$coef[c(cumsum(pred$c)[4] + 1) : cumsum(pred$c)[5]]
coef_gen <- m$coef[c(cumsum(pred$c)[5] + 1) : cumsum(pred$c)[6]]

beta_pet_hat <- array(pred$terms[[2]]$getM(use_by=FALSE) %*% coef_pet,
dim=c(dim_roi, n))
beta_age_hat <- array(pred$terms[[3]]$getM(use_by=FALSE) %*% coef_age,
                      dim=c(dim_roi, n))
beta_mci_hat <- array(pred$terms[[4]]$getM(use_by=FALSE) %*% coef_mci,
                      dim=c(dim_roi, n))
beta_edu_hat <- array(pred$terms[[5]]$getM(use_by=FALSE) %*% coef_edu,
                      dim=c(dim_roi, n))
beta_gen_hat <- array(pred$terms[[6]]$getM(use_by=FALSE) %*% coef_gen,
                      dim=c(dim_roi, n))

pdf('GMOUT_LARGEROI_beta_pet.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(3,3))

brks <- seq(min(beta_pet_hat), max(beta_pet_hat), l=25) 
for(slice in slices){
  image(beta_pet_hat[,,slice,1], main=bquote("beta_pet, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()


pdf('LARGEROI_beta_age.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(2,3))

brks <- seq(min(beta_age_hat), max(beta_age_hat), l=25) 
for(slice in slices){
  image(beta_age_hat[,,slice,1], main=bquote("beta_age, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()

pdf('LARGEROI_beta_mci.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(2,3))

brks <- seq(min(beta_mci_hat), max(beta_mci_hat), l=25) 
for(slice in slices){
  image(beta_mci_hat[,,slice,1], main=bquote("beta_mci, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()


pdf('LARGEROI_beta_edu.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(2,3))

brks <- seq(min(beta_edu_hat), max(beta_edu_hat), l=25) 
for(slice in slices){
  image(beta_edu_hat[,,slice,1], main=bquote("beta_edu, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()


pdf('LARGEROI_beta_gen.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(2,3))

brks <- seq(min(beta_gen_hat), max(beta_gen_hat), l=25) 
for(slice in slices){
  image(beta_gen_hat[,,slice,1], main=bquote("beta_gen, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()

gm <- array(large.roi$grey.matter.4, c(dim_roi, n))
pet <- array(large.roi$pet.slopes.2.1, c(dim_roi, n))

mean.pet <- apply(pet, c(1,2,3), mean)
sd.pet <- apply(pet, c(1,2,3), sd)
mean.gm <- apply(gm, c(1,2,3), mean)
sd.gm <- apply(gm, c(1,2,3), sd)

mean.pet.array <- array(mean.pet, c(dim_roi, n))
mean.gm.array <- array(mean.gm, c(dim_roi, n))

gm[1,1,1,]-mean.gm[1,1,1]
x <- gm - mean.gm.array
x[1,1,1,]

int.matrix <- (pet - mean.pet.array)*(gm - mean.gm.array)

cor.matrix <- apply(int.matrix, c(1,2,3), mean)
cor.matrix <- cor.matrix/(sd.pet*sd.gm)*93/94

pdf('corelations.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(2,3))

brks <- seq(min(cor.matrix), max(cor.matrix), l=25) 
for(slice in slices){
  image(cor.matrix[,,slice], main=bquote("cor_matrix, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()

pdf('histcorr.pdf')
hist(cor.matrix, breaks = 15)
dev.off()


##################################################################################
## Penalized Regression
##################################################################################


##add a difference penalty to all covariates (note that we do not add the penalty to the intercept)

pred$terms[[2]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[2]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[3]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[4]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[5]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[6]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))

system.time(test_opt <- nplarge_pen_fit(x=pred, y=resp, optimizer="optim-BFGS", 
                                        optimizer.control=list(trace=3)))

y_hat <- array(test_opt$fitted.values, c(dim_roi, n))


















##################################################################################
## Do regression in the opposite direction...
##################################################################################

load('largeROIexplore.rda')

n <-  94
dim_roi <- c(7, 7, 7)
n_dims <- length(dim_roi) 
dims <- lapply(dim_roi, function(x) seq_len(x))
names(dims) <- paste0("dim", 1:n_dims)


y <- array(large.roi$pet.slopes.2.1, c(dim_roi, n))

y_cube <- tbl_cube(dimensions=append(dims, list(Subject=subject.ids)), 
                   measures = list(y=y))

x <- array(large.roi$grey.matter.4, c(dim_roi, n))


x_cube <- tbl_cube(dimensions=append(dims, list(Subject=subject.ids)), 
                   measures = list(x=x))



keep <- c("AGE", "PTGENDER", "MCI", "PTEDUCAT", "Subject")
data <- baseline.subjects[keep]

denv <- dataenv(data,  dims="Subject")
add_features(data=x_cube, dataenv=denv)
add_dims(data=x_cube, dataenv=denv)


##Creat a compressed data representaion of data cube tbl
## denv <- dataenv(cov_cube)

## Create intercept term for subjects
const_subj <- const(n, denv)

##Create intercept term for the dimmensions of the ROI
const_v <- const(prod(dim_roi), denv)

##Create bspline basis
bs_dims <- lapply(paste0("dim", 1:n_dims), 
                  function(dim) bspline(dim, denv, df=6, constraint="none"))
bs_v <- do.call(kron, bs_dims)

##Create intercept array 
intercept <- kron(const_v, const_subj)

##create bspline basis for age, mci, education, and gender
bs_age <- kron(bs_v, const(n, denv, by=baseline.subjects$AGE))
bs_mci <- kron(bs_v, const(n, denv, by=baseline.subjects$MCI))
bs_edu <- kron(bs_v, const(n, denv, by=baseline.subjects$PTEDUCAT))
bs_gen <- kron(bs_v, const(n, denv, by=baseline.subjects$PTGENDER))

##create bspline basis for grey matter
bs_x <- kron(bs_v, const_subj, by=as.vector(x_cube$mets$x))

(pred <- concat(intercept, bs_x, bs_age, bs_mci, bs_edu, bs_gen))
resp <- as.vector(y_cube$mets$y)


##################################################################################
## Fit Unpenalized Model 
##################################################################################

system.time(m <- nplarge_lm_fit(pred, resp))

##144.527

y_hat <- array(m$fitted.values, c(dim_roi, n))

set.seed(33)

subjects <- sample(1:n, 10)
slices <- sample(1:dim_roi[3], 10, replace = TRUE)
clrs <- heat.colors(24)
brks <- seq(min(resp), max(resp), l=25) 

pdf('PETOUT__LARGEROI_scalar_reg_est.pdf')
for(i in 1:10){
  par(mfrow = c(1,2))
  image(y_cube$mets$y[, , slices[i], subjects[i]], col=clrs, breaks=brks,
        main="Change in PET")
  image(y_hat[, , slices[i], subjects[i]], col=clrs, breaks=brks, 
        main="Estimate")    
}
dev.off()



m$coef[1]
coef_pet <- m$coef[2 : cumsum(pred$c)[2]]
coef_age <- m$coef[c(cumsum(pred$c)[2] + 1) : cumsum(pred$c)[3]]
coef_mci <- m$coef[c(cumsum(pred$c)[3] + 1) : cumsum(pred$c)[4]]
coef_edu <- m$coef[c(cumsum(pred$c)[4] + 1) : cumsum(pred$c)[5]]
coef_gen <- m$coef[c(cumsum(pred$c)[5] + 1) : cumsum(pred$c)[6]]

beta_pet_hat <- array(pred$terms[[2]]$getM(use_by=FALSE) %*% coef_pet,
                      dim=c(dim_roi, n))
beta_age_hat <- array(pred$terms[[3]]$getM(use_by=FALSE) %*% coef_age,
                      dim=c(dim_roi, n))
beta_mci_hat <- array(pred$terms[[4]]$getM(use_by=FALSE) %*% coef_mci,
                      dim=c(dim_roi, n))
beta_edu_hat <- array(pred$terms[[5]]$getM(use_by=FALSE) %*% coef_edu,
                      dim=c(dim_roi, n))
beta_gen_hat <- array(pred$terms[[6]]$getM(use_by=FALSE) %*% coef_gen,
                      dim=c(dim_roi, n))

pdf('PETOUT_LARGEROI_beta_gm.pdf')
clrs <- heat.colors(24)
slices <- round(seq(1, dim_roi[3]))
par(mfrow = c(3,3))

brks <- seq(min(beta_pet_hat), max(beta_pet_hat), l=25) 
for(slice in slices){
  image(beta_pet_hat[,,slice,1], main=bquote("beta_gm, slice "~.(slice)),
        col=clrs, breaks=brks)
}
dev.off()


##################################################################################
## Penalized Regression
##################################################################################


##add a difference penalty to all covariates (note that we do not add the penalty to the intercept)

pred$terms[[2]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[2]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[3]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[4]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[5]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))
pred$terms[[6]]$pen <- do.call(make_kron_sum, 
                               lapply(pred$terms[[3]]$terms[1:n_dims], 
                                      make_diffpen))

system.time(test_opt <- nplarge_pen_fit(x=pred, y=resp, optimizer="optim-BFGS", 
                                        optimizer.control=list(trace=3)))

y_hat <- array(test_opt$fitted.values, c(dim_roi, n))

