
#library(survival)
#if(!require(flexsurv)){
#  install.packages("flexsurv", repos = "http://cran.us.r-project.org")
#}

library(flexsurv)
library(ggplot2)

#library(doFuture)
library(magrittr)

#library(foreach)
#library(doParallel)


library(data.table)


# This function generate random data from rweibull dist based on the given params
# Currently it saves the data to a specific location
generate_rweibull_data <- function(Beta = 1.9, theta0 = 0.1, theta1 = 0.2, n=20, N = 10000, data_path = ""){
  # generate a random initial design
  des = sample.int(10, n, replace = T)
  candidate = 1:10
  n <- length(des)  # Number of values in each simulated dataset
  
  # Beta <- 1.9 #shape
  # theta0 <- 0.1
  # theta1 <- 0.2
  
  # Parameter Gamma depends on the regression model. Whereas Beta is fixed.
  Gamma <- exp(theta0+theta1*des) 
  runs = N
  for(run in 1:runs){
    rand.data = data.frame(matrix(ncol = n, nrow = length(candidate)))
    colnames(rand.data) = des
    rownames(rand.data) = candidate
    for(row.num in 1:length(candidate)){
      tim<- rweibull(n, shape = Beta, scale = Gamma)
      rand.data[row.num,] <- tim
      
    }
    #print(min(rand.data))
    #print(max(rand.data))
    #print('---')
    
    write.csv(rand.data, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
  }
}

# This function loads existing or generates random data from rweibull dist based on the given params
# Currently it saves the data to a specific location, loads it into a 3D matrix
# It returns the 3D data matrix
load_or_generate_rweibull_data_matrix <- function(Beta = 1.9, theta0 = 0.1, theta1 = 0.2, n=20, N = 10000, data_path){
  # if data path is not there, create it!
  if(!dir.exists(data_path)){
    dir.create(data_path)
  }
  
  # check if folder already has csv files, load them, otherwise generate new data
  temp = list.files(data_path, full.names = T, pattern="*.csv")
  # rm = file.remove(temp)
  if(length(temp) == N){
    print(length(temp))
    print('Data already exists, load it!')
    myfiles = lapply(temp, read.csv)
    # this is the 3D matrix
    data.matrix = array(unlist(myfiles), dim=c(dim(myfiles[[1]]), length(myfiles)))
  } else {
    # generate a random initial design
    des = sample.int(10, n, replace = T)
    candidate = 1:10
    n <- length(des)  # Number of values in each simulated dataset
    
    # Parameter Gamma depends on the regression model. Whereas Beta is fixed.
    Gamma <- exp(theta0+theta1*des) 
    runs = N
    for(run in 1:runs){
      rand.data = data.frame(matrix(ncol = n, nrow = length(candidate)))
      colnames(rand.data) = des
      rownames(rand.data) = candidate
      for(row.num in 1:length(candidate)){
        tim<- rweibull(n, shape = Beta, scale = Gamma)
        rand.data[row.num,] <- tim
        
      }
      
      write.csv(rand.data, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
    }
    
    # this code loads the data from csv files into a 3D matrix
    temp = list.files(data_path, full.names = T, pattern="*.csv")
    
    myfiles = lapply(temp, read.csv)
    # this is the 3D matrix
    data.matrix = array(unlist(myfiles), dim=c(dim(myfiles[[1]]), length(myfiles)))
  }
  return(data.matrix)
}


generate_rweibull_data2 <- function(n = 20, candidate = 1:10 ,Beta = 1.9, theta0 = 0.1, theta1 = 0.2, N = 10000, data_path = ""){
  # generate a random initial design
  des = sample.int(10, n, replace = T)
  #candidate = 1:10
  #n <- length(des)  # Number of values in each simulated dataset
  
  # Beta <- 1.9 #shape
  # theta0 <- 0.1
  # theta1 <- 0.2
  
  # Parameter Gamma depends on the regression model. Whereas Beta is fixed.
  Gamma <- exp(theta0+theta1*des) 
  runs = N
  for(run in 1:runs){
    rand.data = data.frame(matrix(ncol = n, nrow = length(candidate)))
    colnames(rand.data) = des
    rownames(rand.data) = candidate
    for(row.num in 1:length(candidate)){
      tim<- rweibull(n, shape = Beta, scale = Gamma)
      rand.data[row.num,] <- tim
    }
    #print(min(rand.data))
    #print(max(rand.data))
    #print('---')
    
    write.csv(rand.data, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
  }
}

generate_rweibull_data3 <- function(n, candidate ,Beta = 1.9, theta0 = 0.1, theta1 = 0.2, N = 10000, data_path = ""){
  # generate a random initial design
  candidate_length = length(candidate)
  prob = c(0.4,0.4, rep(0.2/(candidate_length-2), (candidate_length-2)) )
  result = rep(candidate, round(n * prob))
  des = sample(result)
  
  #des = sample.int(10, n, replace = T)
  #candidate = 1:10
  #n <- length(des)  # Number of values in each simulated dataset
  
  # Beta <- 1.9 #shape
  # theta0 <- 0.1
  # theta1 <- 0.2
  
  # Parameter Gamma depends on the regression model. Whereas Beta is fixed.
  Gamma <- exp(theta0+theta1*des) 
  runs = N
  for(run in 1:runs){
    rand.data = data.frame(matrix(ncol = n, nrow = length(candidate)))
    #colnames(rand.data) = des
    #rownames(rand.data) = candidate
    for(row.num in 1:length(candidate)){
      tim<- rweibull(n, shape = Beta, scale = Gamma)
      rand.data[row.num,] <- tim
    }
    #print(min(rand.data))
    #print(max(rand.data))
    #print('---')
    
    write.csv(rand.data, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
  }
}

## We generate random data using the qweibull function
#The choice of (1:n)/N as the input to the qweibull() function in the 
#systematic sampling approach is based on the desire to obtain a set of evenly 
#spaced quantile values between 0 and 1.
#Systematic sampling intervals: The resulting quantile probabilities (1:n)/N
# create equally spaced intervals along the cumulative distribution function (CDF) 
# of the Weibull distribution. These intervals facilitate systematic sampling, 
# ensuring that samples are taken at regular intervals across the distribution.
#By using (1:n)/N as the input to the qweibull() function, 
# we obtain a systematic sample of quantile values from the Weibull distribution
# that covers the entire range of the distribution in a systematic and evenly 
# spaced manner.
#It's important to note that this choice assumes that the data is 
# well-represented by the Weibull distribution and that equally spaced intervals 
# are appropriate for your specific analysis. 
# If your data deviates significantly from the Weibull distribution or if there 
# are specific requirements or assumptions in your analysis, 
# you may need to adjust the systematic sampling approach accordingly.

generate_qweibull_data <- function(initial_design, candidate ,Beta = 1.9, theta0 = 0.1, theta1 = 0.2, N = 10000, data_path = ""){
  # generate a random initial design
  candidate_length = length(candidate)
  n = length(initial_design)
  #prob = c(0.4,0.4, rep(0.2/(candidate_length-2), (candidate_length-2)) )
  #result = rep(candidate, round(n * prob))
  #des = sample(result)
  
  # Parameter Gamma depends on the regression model. Whereas Beta is fixed.
  Gamma <- exp(theta0+theta1*initial_design) 
  runs = N
  for(run in 1:runs){
    rand.data = data.frame(matrix(ncol = n, nrow = length(candidate)))
    #colnames(rand.data) = des
    #rownames(rand.data) = candidate
    for(row.num in 1:length(candidate)){
      # here we change the intervals in qweibull by selecting a random number
      # between 1 and N
      # repeat until all data generated doesn't contain any NaN
      while(TRUE){
        tim<- qweibull((1:n)/(sample(N, 1)), shape = Beta, scale = Gamma)
        if(!any(is.nan(tim))){
          rand.data[row.num,] <- tim
          break
        }
      }
    }
    #print(min(rand.data))
    #print(max(rand.data))
    #print('---')
    
    write.csv(rand.data, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
  }
}

generate_bootstrapped_qweibull_data <- function(initial_design, candidate ,Beta = 1.9, theta0 = 0.1, theta1 = 0.2, N = 10000, data_path = ""){
  
  # generate a random initial design
  candidate_length = length(candidate)
  num_row_samples = N %/% candidate_length
  n = length(initial_design)
  #prob = c(0.4,0.4, rep(0.2/(candidate_length-2), (candidate_length-2)) )
  #result = rep(candidate, round(n * prob))
  #des = sample(result)
  
  # Parameter Gamma depends on the regression model. Whereas Beta is fixed.
  Gamma <- exp(theta0+theta1*initial_design) 
  #runs = N / num_row_samples
  # if it's the very first time, we generate random data based on the initial design
  rand.data = data.frame(matrix(ncol = n, nrow = length(candidate)))
  for(row.num in 1:length(candidate)){
    # here we change the intervals in qweibull by selecting a random number
    # between 1 and N
    # repeat until all data generated doesn't contain any NaN
    while(TRUE){
      tim<- qweibull((1:n)/(sample(N, 1)), shape = Beta, scale = Gamma)
      if(!any(is.nan(tim))){
        rand.data[row.num,] <- tim
        break
      }
    }
  }
  #write.csv(rand.data, file = paste(data_path,"original_.csv.t", sep=""), row.names = F)
  
  # if it's the second time or after, we take a bootstrap sample from the 
  # data matrix generated in the first step
  for(sn in 1:num_row_samples){
    for(rn in 1:length(candidate)){
      rand.subset = data.frame(matrix(ncol = n, nrow = length(candidate)))
      for(sample.row in 1:length(candidate)){
        rnd_row = sample(rand.data[rn,], replace = T)
        rand.subset[sample.row,] <- rnd_row
      }
      write.csv(rand.subset, file = paste(data_path,"sample_",as.character(rn),"_",as.character(sn),".csv", sep=""), row.names = F)
    }
  }
  
  #  for(run in 1:runs){
  
  #    if(run == 1){
  #colnames(rand.data) = des
  #rownames(rand.data) = candidate
  
  #print(min(rand.data))
  #print(max(rand.data))
  #print('---')
  
  #      write.csv(rand.data, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
  #    } else {
  
  #rand.data1 = sample_n(rand.data, nrow(rand.data), replace = T)
  #write.csv(rand.data1, file = paste(data_path,"sample_",as.character(run),".csv", sep=""), row.names = F)
  #    }
  #  }
}






generate_rweibull_data_matrix <- function(initial_design, candidate = 1:10, data_path, N, theta01, theta1, input_beta, censor, method, rep_num){
  # if data path is not there, create it!
  if(!dir.exists(data_path)){
    dir.create(data_path)
  }
  
  # check if folder already has csv files, load them, otherwise generate new data
  temp = list.files(data_path, full.names = T, pattern="*.csv")
  # rm = file.remove(temp)
  if(length(temp) == N){
    print(length(temp))
    print('Data already exists, load it!')
    myfiles = lapply(temp, read.csv)
    # this is the 3D matrix
    data.matrix = array(unlist(myfiles), dim=c(dim(myfiles[[1]]), length(myfiles)))
  } else {
    # Generate random data from weibull dist based on the given params
    #generate_qweibull_data(n,candidate, Beta=input_beta, theta1 = theta1, N = N, data_path = data_path)
    generate_bootstrapped_qweibull_data(initial_design,candidate, Beta=input_beta, theta0 = theta01, theta1 = theta1, N = N, data_path = data_path)
    
    #old
    #generate_rweibull_data3(n,candidate, Beta=input_beta, theta1 = theta1, N = N, data_path = data_path)
    
    # this code loads the data from csv files into a 3D matrix
    temp = list.files(data_path, full.names = T, pattern="*.csv")
    
    myfiles = lapply(temp, read.csv)
    # this is the 3D matrix
    data.matrix = array(unlist(myfiles), dim=c(dim(myfiles[[1]]), length(myfiles)))
  }
  return(data.matrix)
}

# plot density of data of optima design
# saves the plot to the given file_name
generate_density_plot <- function(optimal.design, data_matrix, file.name, N){
  # load data in 3D matrix
  #temp = list.files(data_path, full.names = T, pattern="*.csv")
  #print(data_path)
  #print(file.name)
  #myfiles = lapply(temp, read.csv)
  #print(length(myfiles))
  # this is the 3D matrix
  #data.matrix = array(unlist(myfiles), dim=c(dim(myfiles[[1]]), length(myfiles)))
  # load data for this best design
  #print(data_matrix)
  optimal.design.data = load_design_data(optimal.design, data_matrix, N) 
  #print('Generate density plot: data loaded')
  # get index values for 1 and 10 in  optimal.design
  loc_1 = which(optimal.design %in% 1)
  loc_10 = which(optimal.design %in% 10)
  
  # get corresponding rows for 1 and 10 from optimal.design.data
  data_1 = as.vector(optimal.design.data[loc_1,])
  data_10 = as.vector(optimal.design.data[loc_10,])
  
  df1 <- data.frame(x = data_1, label=rep('1 Data', length(data_1)))
  df2 <- data.frame(x = data_10, label=rep('10 Data', length(data_10)))
  df=rbind(df1, df2)
  #print('Generate density plot: data filtered, now creating plot!')
  myplot <- ggplot(df, aes(x=x, color=label, fill = label, ..scaled..)) + geom_density(alpha = 0.3)
  ggsave(file.name, device = "png")
  cat("Saved plot to", file.name, "\n")
  #unlink()
  #png(file.name)
  #print(myplot)
  #dev.off()
}

# plot density of data of optima design for competing risk
# we send it the min matrix and grab data from it!
# saves the plot to the given file_name
generate_comp_density_plot <- function(candidate, optimal.design, min.matrix, file.name, N){
  cat('Plotting .. Min Matrix Dims = ', dim(min.matrix), '\n')
  # load data for this best design
  #print(optimal.design)
  optimal.design.data = load_design_data2(candidate, optimal.design, min.matrix, N)
  #optimal.design.data = load_design_data(optimal.design, min.matrix, N) 
  #print('Generate density plot: data loaded')
  # get index values for 1 and 10 in  optimal.design
  loc_1 = which(optimal.design %in% 1)
  loc_10 = which(optimal.design %in% 10)
  
  # get corresponding rows for 1 and 10 from optimal.design.data
  data_1 = as.vector(optimal.design.data[loc_1,])
  data_10 = as.vector(optimal.design.data[loc_10,])
  
  df1 <- data.frame(x = data_1, label=rep('1 Data', length(data_1)))
  df2 <- data.frame(x = data_10, label=rep('10 Data', length(data_10)))
  df=rbind(df1, df2)
  print('Generate density plot: data filtered, now creating plot!')
  myplot <- ggplot(df, aes(x=x, color=label, fill = label, ..scaled..)) + geom_density(alpha = 0.3)
  ggsave(file.name, device = "png")
  cat("Saved plot to", file.name, "\n")
  #unlink()
  #png(file.name)
  #print(myplot)
  #dev.off()
}


############################
# function 1 to create the covariance matrix and present it in savedata
# x "the design" is the input of the first function
# For each design calcualte Gamma, generate data from weibull regression model, calculate the covaraince matrix for the dataset.
# n is the Number of values in each simulated dataset
# old
find_cov_matrix_OLD<- function(x) {
  N <- 1000 # NUMBER OF SIMULATED DATA
  #n <- 10  # Number of values in each simulated dataset
  Beta<-5 #shape
  theta0 <- 0.1
  theta1 <- -0.2
  
  #k=1
  #numberofsamples = 1000
  #samplesize = 10
  n = length(x)
  savedata = as.data.frame(matrix(0, ncol = N, nrow = 9))
  Gamma <- exp(theta0+theta1*x) # Parameter Gamma depends on the regression model. Whereas Beta is fixed.  
  #print(Gamma)
  #print('Cov_matri for')
  #print(x)
  for (i in 1:N) {
    # seed needs to be set inside this loop to make sure results are consistent
    #set.seed(seed)
    tim<- rweibull(n, shape = Beta, scale = Gamma)
    tim2 <- Surv(time= tim, event=rep(1,n),type="right") #no censoring
    # the following line sometimes throws the following error
    # error happens when using different random seeds (not clear to me why) 
    # Error in if (any(singular)) fit$coeffients[singular] <- NA : 
    #  missing value where TRUE/FALSE needed
    surv_tim<-survreg(tim2 ~ x, dist="weibull", na.action = na.omit)
    var_cov = surv_tim$var
    var_cov = as.data.frame(var_cov)
    k = melt(setDT(var_cov))
    savedata[,i] = k[,2]
  }
  #print('Done') 
  return(savedata)
}


find_cov_matrix2 <- function(df, des, censored = F, Y = 0) {
  #this function receives a n by N DF created from saved simulated data
  #it loops through it column by column and fits a Surv regression model
  #then it uses the resulting values to create a cov matrix
  ncols = ncol(df)
  n = nrow(df)
  savedata = as.data.frame(matrix(0, ncol = n, nrow = 9))
  #write.csv(df, file = "C:/Users/han-s/Downloads/df_.csv", row.names = F)
  
  # create a censoring_indicator if the data has been censored
  if(censored){
    censoring_indicator = as.data.frame(ifelse(as.matrix(df) == Y, 0, 1))
    #print(table(censoring_indicator))
  } else{
    censoring_indicator = 0
  }
  
  for (i in 1:ncols) {
    # seed needs to be set inside this loop to make sure results are consistent
    #set.seed(seed)
    tim  <- df[, i]
    tim2 = 0
    # if censoring has been applied then create an indicator df
    # 0 means the entry was censored, 1 means otherwise
    if(censored){
      #print(length(censoring_indicator[, i]))
      tim2 <- Surv(time = tim, event=censoring_indicator[, i],type="right") #censoring
    } else{
      tim2 <- Surv(time = tim, event=rep(1,n),type="right") #no censoring
    }
    
    # the following line sometimes throws the following error
    # error happens when using different random seeds (not clear to me why) 
    # Error in if (any(singular)) fit$coeffients[singular] <- NA : 
    # missing value where TRUE/FALSE needed
    
    rr <- try({
      surv_tim <- survreg(tim2 ~ des, dist="weibull", na.action = na.omit)
      var_cov = surv_tim$var
      var_cov = as.data.frame(var_cov)
      k = melt(setDT(var_cov))
      savedata[,i] = k[,2]
    }, silent = T) # END try
    # if the following error has happened then fill the ith column with NAs
    # Error in if (any(singular)) fit$coeffients[singular] <- NA
    # These NAs are dealt with later (i.e. in find_trace & find_det functions)
    if(class(rr) == "try-error"){
      savedata[,i] = NA
    }
    
  }
  #cat('No of rows', nrow(savedata), 'No of cols',ncol(savedata),'\n')
  #print('Done') 
  # this code replaces any NA values with the mean of their corresponding values
  savedata <- t(savedata)
  data2 <- savedata     # Duplicate data frame
  for(i in 1:ncol(savedata)) {    # Replace NA in all columns
    data2[ , i][is.na(data2[ , i])] <- mean(data2[ , i], na.rm = TRUE)
  }
  savedata <- t(data2)
  
  return(savedata)
}  


## This function receives an estimate matrix of 3 by N
## it computes its cov matrix and returns the det of this cov matrix
find_det_estimate_cov_matrix <- function(estimate_matrix, compute_cov = F) {
  #print(estimate_matrix)
  if(compute_cov){
    estimate_cov_matrix = cov(estimate_matrix)
  } else {
    estimate_cov_matrix = estimate_matrix
  }
  
  crit = det(estimate_cov_matrix)^(1/3)
  #crit = det(estimate_cov_matrix)^(1/6)
  return(crit)
}

# This function returns a matrix of three estimates
## i.e. beta_hat, theta0_hat and theta1_hat
find_estimate_matrix <- function(df, des, censored = F, Y = 0) {
  #this function receives a n by N DF created from saved simulated data
  #it loops through it column by column and fits a Surv regression model
  #then it uses the resulting values to create a matrix that has 
  #the shape, scale and des extracted from the model generated by flexsurvreg
  ncols = ncol(df)
  n = nrow(df)
  savedata = as.data.frame(matrix(0, ncol = n, nrow = 3))
  #write.csv(df, file = "C:/Users/han-s/Downloads/df_.csv", row.names = F)
  
  # create a censoring_indicator if the data has been censored
  if(censored){
    censoring_indicator = as.data.frame(ifelse(as.matrix(df) == Y, 0, 1))
    #print(table(censoring_indicator))
  } else{
    censoring_indicator = 0
  }
  
  for (i in 1:ncols) {
    
    # seed needs to be set inside this loop to make sure results are consistent
    #set.seed(seed)
    tim  <- df[, i]
    tim2 = 0
    # if censoring has been applied then create an indicator df
    # 0 means the entry was censored, 1 means otherwise
    if(censored){
      #print(length(censoring_indicator[, i]))
      tim2 <- Surv(time = tim, event=censoring_indicator[, i],type="right") #censoring
    } else{
      tim2 <- Surv(time = tim, event=rep(1,n),type="right") #no censoring
    }
    
    # the following line sometimes throws the following error
    # error happens when using different random seeds (not clear to me why) 
    # Error in if (any(singular)) fit$coeffients[singular] <- NA : 
    # missing value where TRUE/FALSE needed
    
    rr <- try({
      #surv_tim <- survreg(tim2 ~ des, dist="weibull", na.action = na.omit)
      #print(tim2)
      
      surv_tim <- flexsurvreg(tim2 ~ des, dist="weibull", na.action = na.omit)
      coefficients = coef(surv_tim)
      #print(coefficients)
      beta_hat = exp(coefficients[1])
      theta0_hat = coefficients[2]
      theta1_hat = coefficients[3]
      
      #shape = coef(surv_tim)['shape']
      #scale = coef(surv_tim)['scale']
      #des_ = coef(surv_tim)['des']
      v = c(beta_hat, theta0_hat, theta1_hat)
      savedata[,i] = v
    }, silent = T) # END try
    # if the following error has happened then fill the ith column with NAs
    # Error in if (any(singular)) fit$coeffients[singular] <- NA
    # These NAs are dealt with later (i.e. in find_trace & find_det functions)
    if(class(rr) == "try-error"){
      #print("Error")
      savedata[,i] = NA
    }
    
  }
  #cat('No of rows', nrow(savedata), 'No of cols',ncol(savedata),'\n')
  #print('Done') 
  # this code replaces any NA values with the mean of their corresponding values
  savedata <- t(savedata)
  data2 <- savedata     # Duplicate data frame
  for(i in 1:ncol(savedata)) {    # Replace NA in all columns
    data2[ , i][is.na(data2[ , i])] <- mean(data2[ , i], na.rm = TRUE)
  }
  #savedata <- t(data2)
  #cat('Mean theta1 = ',mean(savedata[3,], na.rm = T),'\n')
  return(data2)
}  



# this fun computes weights from a cov matrix
# receives the cov matrix of a DOpt design
# returns weights as a vector of 3 values
find_weight_vector <-function(dopt.design.cov_mx){
  # compute corresponding weights
  w1 = 1/dopt.design.cov_mx[1,1]
  w2 = 1/dopt.design.cov_mx[2,2]
  w3 = 1/dopt.design.cov_mx[3,3]
  return(c(w1,w2,w3))
}

# this fun computes weights from a cov matrix
# receives the cov matrix of a DOpt design
# returns weights as a vector of 6 values
# values are in this order: "shape1", "scale1", "design1", "shape2", "scale2", "design2"
find_weight_vector6 <-function(dopt.design.cov_mx){
  # compute corresponding weights
  w1 = 1/dopt.design.cov_mx[1,1]
  w2 = 1/dopt.design.cov_mx[2,2]
  w3 = 1/dopt.design.cov_mx[3,3]
  w4 = 1/dopt.design.cov_mx[4,4]
  w5 = 1/dopt.design.cov_mx[5,5]
  w6 = 1/dopt.design.cov_mx[6,6]
  return(c(w1,w2,w3,w4,w5,w6))
}


## This function computes the crit of weighted a optimal
# receives a vector of weights and the est matrix of the design
# returns the crit value as a weighted avg
find_a_opt_weighted_crit_value <- function(dweights, est_mx){
  #print(dweights)
  # get the cov matrix of the initial design
  design.cov_mx = cov(est_mx)
  v1 = design.cov_mx[1,1] * dweights[1] #beta
  v2 = design.cov_mx[2,2] * dweights[2] #theta0
  v3 = design.cov_mx[3,3] * dweights[3] #theta1
  crit_value = (v1 + v2 + v3)/sum(dweights)
  
  return(crit_value)
}

## This function computes the crit of weighted a optimal
# receives a vector of weights and the cov matrix of the design
# returns the crit value as a weighted avg
find_a_opt_weighted_crit_value6 <- function(dweights, design.cov_mx){
  v1 = design.cov_mx[1,1] * dweights[1] #beta
  v2 = design.cov_mx[2,2] * dweights[2] #theta0
  v3 = design.cov_mx[3,3] * dweights[3] #theta1
  v4 = design.cov_mx[4,4] * dweights[4] #beta
  v5 = design.cov_mx[5,5] * dweights[5] #theta0
  v6 = design.cov_mx[6,6] * dweights[6] #theta1
  crit_value = (v1 + v2 + v3 + v4 + v5 + v6)/sum(dweights)
  
  return(crit_value)
}


# this fun generates weights for the opt.design (des) given the data matrix and
# other params
find_a_optimal_design_weights <- function(data.matrix, N, des, censored = F, Y = 0){
  # grab the design data
  design.data = load_design_data(des, data.matrix, N) 
  
  # we compute the cov matrix of the design (9 rows by N columns)
  cov_mx = find_cov_matrix2(design.data, des, censored = censored, Y = Y)
  #print(cov_mx)
  # remove columns that are fully NA's
  df <- cov_mx[,colSums(is.na(cov_mx))<nrow(cov_mx)]
  N = ncol(df)
  weights <- as.data.frame(matrix(0, ncol = N, nrow = 3))
  for (i in 1:N)   {
    #print(df[,i])
    A = matrix(df[,i], nrow=3,ncol=3)
    #print(A)
    w0 = 1/A[1,1]
    w1 = 1/A[2,2]
    w2 = 1/A[3,3]
    weights[i] = c(w0,w1,w2)
  }
  # makes the weights DF of the shape 9 by [w0,w1,w2]
  weights = t(weights)
  
  # this code replaces any Inf values with the mean of their corresponding values
  weights[is.infinite(weights)] <- NA
  
  for(i in 1:ncol(weights)) {    # Replace NA in all columns
    weights[is.na(weights[,i]), i] <- mean(weights[,i], na.rm = TRUE)
  }
  
  #log_ret[which(!is.finite(log_ret))] <- 0
  
  
  return(weights)
}

# ==================================================================
# Given a design, a candidate and a 3D matrix of pre-generated data
# Return the data for that design according to the candidate
# ==================================================================
load_design_data2 <- function(candidate, design, data.matrix, N) {
  n <- length(design)
  #cat('Min Matrix Dims = ', dim(data.matrix), '\n')
  #print(design)
  # create an empty df to save the corresponding 1000 values of each treatment, subject combination
  design.data = data.frame(matrix(ncol = n, nrow = N))
  
  for(subject in 1:n){
    treatment = design[subject]
    # get the treatment ID from cand
    treatment_id = which(candidate == treatment)
    #cat('Treatment ID:', treatment_id,'\n')
    #cat('Subj:', subject,'\n')
    v = data.matrix[treatment_id, subject, ]
    #cat('V:', v,'\n')
    design.data[,subject] = v
    #print(length(v))
    
  }
  # transpose the matrix so that each design component is a row
  design.data = t(design.data)
  rownames(design.data) <- design
  return(design.data)
}


# ==================================================================
# Given a design and a 3D matrix of pre-generated data
# Return the data for that design
# ==================================================================
load_design_data <- function(initial_design, data.matrix, N) {
  n <- length(initial_design)
  #print(initial_design)
  #cat('Min Matrix Dims = ', dim(data.matrix), '\n')
  # create an empty df to save the corresponding 1000 values of each treatment, subject combination
  #print(n)
  #print(N)
  design.data = data.frame(matrix(ncol = n, nrow = N))
  
  for(subject in 1:n){
    treatment = initial_design[subject]
    #cat('Treatment:', treatment,'\n')
    #cat('Subj:', subject,'\n')
    v = data.matrix[treatment, subject, ]
    #cat('V:', v,'\n')
    design.data[,subject] = v
    #print(length(v))
    
  }
  # transpose the matrix so that each design component is a row
  design.data = t(design.data)
  rownames(design.data) <- initial_design
  return(design.data)
}

# ==================================================================
# Two functions to compute the determinant and Trace of the cov matrix columns.
# ==================================================================
######################################################
find_trace  <- function(saveddata) {
  # remove columns that are fully NA's
  df <- saveddata[,colSums(is.na(saveddata))<nrow(saveddata)]
  N = ncol(df)
  trace_savedata <- as.data.frame(matrix(0, ncol = N, nrow = 1))
  for (i in 1:N)   {
    A = matrix(df[,i], nrow=3,ncol=3)
    trace_savedata[i] = (sum(diag(A)))/3
  }
  # return the computed trace and det values in a DF
  return(t(trace_savedata))
}

find_det  <- function(saveddata) {
  # remove columns that are fully NA's
  df <- saveddata[,colSums(is.na(saveddata))<nrow(saveddata)]
  N = ncol(df)
  
  det_savedata <- as.data.frame(matrix(0, ncol = N, nrow = 1))
  for (i in 1:N)   {
    A = matrix(df[,i], nrow=3,ncol=3)
    det_savedata[i]   = (det(A))^(1/3)
  }
  # return the computed trace and det values in a DF
  return(t(det_savedata))
}

## This fun is for Ds optimal design, for theta_1
# we're interested in a subset of params S (and here we have special case where S = 1)
find_theta1_var_OLD  <- function(saveddata){
  # remove columns that are fully NA's
  df <- saveddata[,colSums(is.na(saveddata))<nrow(saveddata)]
  N = ncol(df)
  
  theta1_var_savedata <- as.data.frame(matrix(0, ncol = N, nrow = 1))
  for (i in 1:N){
    A = matrix(df[,i], nrow=3,ncol=3)
    theta1_var_savedata[i]   = A[2,2]
  }
  # return the computed trace and det values in a DF
  return(t(theta1_var_savedata))
}

# this function receives an estimate matrix (N by 3)
# returns the variance of the 3rd column (theta1)
find_theta1_var  <- function(est_mx){
  # return the variance of the 3rd column
  return(var(est_mx[,3]))
}

## this is to find the crit value for DS for competing risk
# divide the det of the full cov matrix by the det of the cov matrix after
# dropping the rows and columns of the param of interest (theta11 and theta22)
find_competing_ds_crit <- function(cov_mx){
  tmp.cov.mx = cov_mx[c(-3,-6), c(-3,-6)]
  return((det(cov_mx)/det(tmp.cov.mx)) ^ 0.5)
}
# this fun loops thru the cox mx of the design and multiplies
# the diag elements of each 3x3 max by their corresponding weights
# it sums the values and adds them to a vector then returns it
# this vector is the trace_savedata DF (it has 1 column)
find_weighted_avg  <- function(saveddata, weights) {
  # remove columns that are fully NA's
  df <- saveddata[,colSums(is.na(saveddata))<nrow(saveddata)]
  N = ncol(df)
  trace_savedata <- as.data.frame(matrix(0, ncol = N, nrow = 1))
  for (i in 1:N)   {
    A = matrix(df[,i], nrow=3,ncol=3)
    w = weights[i,]
    v = w[1] * A[1,1] + w[2] * A[2,2] + w[3] * A[3,3]
    # divide by 3 to obtain the avg????
    # divide by the sum of weights
    trace_savedata[i] = v/(sum(w))
  }
  # return the computed trace and det values in a DF
  return(t(trace_savedata))
}


# computes the min and indicator matrices which are needed to compute the cov matrix
find_min_indicator_matrices <- function(data.matrix1, data.matrix2, censor = F, Y = 2){
  # create a new matrix that contains the min of corresponding elements
  # from the two matrices
  min.matrix = ifelse(data.matrix1 < data.matrix2, data.matrix1, data.matrix2)
  cat("Min value in Min Matrix =",min(min.matrix),"\n")
  cat("Max value in Min Matrix =",max(min.matrix),"\n")
  indicator.matrix = ifelse(data.matrix1 < data.matrix2, 1, 0)
  
  censor_frac = 0
  #c.min.matrix = 0# censored min matrix
  c.indicator.matrix = indicator.matrix# censored indicator matrix
  # Apply censoring if the user wants to apply it
  if(censor){
    # compute faction of values that are > 1.17 (e.g. which is one year and two months of two years)
    censor_frac = (length(min.matrix[min.matrix > Y])/length(min.matrix))*100
    cat('Fraction of censored values in the Min Matrix:',censor_frac, '\n')
    ## filter data according to Y so that we can calc small y
    # which is the censored data
    min.matrix = ifelse(min.matrix > Y, Y, min.matrix)
    print("Censoring Applied to the Min Matrix!")
  } else {
    print("No Censoring Applied to the Min Matrix!")
  }
  
  if(censor){
    # in case of censoring create an indicator matrix that tells us if 
    # the min element is from censoring or from min.matrix (1 means from censoring)
    c.indicator.matrix = ifelse(min.matrix == Y, 0, 1)
  }
  # list to store all data.frames
  out <- list()
  out$min.matrix <- min.matrix
  out$indicator.matrix <- indicator.matrix
  out$censor_frac = censor_frac
  out$c.indicator.matrix = c.indicator.matrix
  return(out)
}


compute_competing_risk_cov_matrix <- function(min.matrix, indicator.matrix, c.indicator.matrix, design, N, candidate_length, censor = F){
  savedata1 = as.data.frame(matrix(0, ncol = 3, nrow = N))
  savedata2 = as.data.frame(matrix(0, ncol = 3, nrow = N))
  # loop through the N matrices
  for(d in 1:N){
    shape1 = 0
    shape2 = 0
    scale1 = 0
    scale2 = 0
    design1 = 0
    design2 = 0
    # for each matrix of the N matrices, take the matrix row by row
    # and fit flexsurvreg to retrieve the params of the dist
    # then compute the avg of each of those params
    for(r in 1:candidate_length){
      # print first row of min.matrix 
      time = min.matrix[r,,d]
      # print first row of indicator.matrix
      z1 = indicator.matrix[r,,d]
      z2 = ifelse(z1 == 1, 0, 1)
      #cat('No of ones in z1 = ', length(which(z1 %in% 1)), "\n")
      #cat('No of zeros in z1 = ', length(which(z1 %in% 0)), "\n")
      #cat('No of ones in z2 = ', length(which(z2 %in% 1)), "\n")
      #cat('No of zeros in z2 = ', length(which(z2 %in% 0)), "\n")
      #data <- data.frame(design = design, time = time, status1 = z1, status2 = z2)
      #ret <- write.csv(x=data, file='datadata.csv')
      
      #print(data)
      rr <- try({
        if (censor){
          z3 = c.indicator.matrix[r,,d]
          #cat('z3 * z1 =', z1*z3, "\n")
          fitTT <- flexsurvreg(Surv(time, event=z1*z3, type="right") ~ design, dist = "weibull")
        } else{
          #cat('time = ',time,' - design = ',design,' - z1 =', z1, "\n")
          #cat(paste(shQuote(z1, type="cmd"), collapse=", "), "\n")
          fitTT <- flexsurvreg(Surv(time, z1) ~ design, dist = "weibull")
        }
        des1 = coef(fitTT)['design']
        design1 = design1 + des1
        shp1 = exp(coef(fitTT)['shape'])
        shape1 = shape1 + shp1
        scl1 = coef(fitTT)['scale']
        scale1 = scale1 + scl1
        #print('No Fit1 ERROR')
      }, silent = T) # END try
      if(class(rr) == "try-error"){
        print('Fit1 ERROR')
        #write.csv(x=data, file='datadata.csv')
        des1 = 0
        shp1 = 0
        scl1 = 0
      }
      
      rr <- try({
        if (censor){
          z3 = c.indicator.matrix[r,,d]
          #cat('z3 * z2 =', z2*z3, "\n")
          fitTT2 <- flexsurvreg(Surv(time, z2*z3) ~ design, dist = "weibull")
        } else {
          #cat('z2 =', z2, "\n")
          fitTT2 <- flexsurvreg(Surv(time, z2) ~ design, dist = "weibull")
        }
        des2 = coef(fitTT2)['design']
        design2 = design2 + des2
        shp2 = exp(coef(fitTT2)['shape'])
        shape2 = shape2 + shp2
        scl2 = coef(fitTT2)['scale']
        scale2 = scale2 + scl2
        #print('No Fit2 ERROR')
      }, silent = T) # END try
      if(class(rr) == "try-error"){
        #print('Fit2 ERROR')
        des2 = 0
        shp2 = 0
        scl2 = 0
      }
    }
    
    avg_shape1 = shape1/candidate_length
    avg_shape2 = shape2/candidate_length
    
    avg_scale1 = scale1/candidate_length
    avg_scale2 = scale2/candidate_length
    
    avg_design1 = design1/candidate_length
    avg_design2 = design2/candidate_length
    
    v = c(avg_shape1, avg_scale1, avg_design1)
    savedata1[d,] = v
    
    v = c(avg_shape2, avg_scale2, avg_design2)
    savedata2[d,] = v
    
  }
  ## compute the covariance matrix of the two resulting tables
  cov_mx = cov(cbind(savedata1, savedata2))
  colnames(cov_mx) = c("shape1", "scale1", "design1", "shape2", "scale2", "design2")
  rownames(cov_mx) = c("shape1", "scale1", "design1", "shape2", "scale2", "design2")
  return(cov_mx)
}


# to save results at the end of every iterations
# takes all_best_designs which has info about all best designs found so far
# other arguments are the length of a design, method used, where results should be saved,
# where the design data can be found and the value of N (no of Weibull samples),
# whether censoring is applied or not
save_results <- function(all_best_designs, design_length, method, results_path,data_path, N, censor, censor_frac){
  all_best_designs = all_best_designs
  colnames(all_best_designs)[ncol(all_best_designs)] <- "Iter"
  
  # for plotting density ,, the last row in all_best_designs has the best design
  last_row = as.numeric(tail(all_best_designs, n = 1))
  opt.design = last_row[1:design_length]
  n = length(opt.design)
  #opt.design=c(1, 10, 1, 10, 1, 10, 1, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 10, 1, 1  )
  thetas_betas = paste(as.character(theta0),"_",as.character(theta1),"_",as.character(input_beta),"_", as.character(n),"_",as.character(N), sep = "")
  
  if(method == "det"){
    if(censor){ 
      print("Saving to:")
      print(paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
      
      all_best_designs['CensoringFraction'] = censor_frac
      write.csv(all_best_designs, file = paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
      # create and save density plot plot
      #file.name = paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
      #generate_density_plot(opt.design, data_path, file.name, N)
    } else{
      print("Saving to:")
      print(paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
      
      all_best_designs['CensoringFraction'] = 0
      write.csv(all_best_designs, file = paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
      #file.name = paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
      #generate_density_plot(opt.design, data_path, file.name, N)
    }
  } else {
    if(method == "trace"){
      if(censor){ 
        print("Saving to:")
        print(paste(results_path,"/Censoring/WTrace/Result_WTrace_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
        
        all_best_designs['CensoringFraction'] = censor_frac
        write.csv(all_best_designs, file = paste(results_path,"/Censoring/WTrace/Result_WTrace_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        #file.name = paste(results_path,"/Censoring/WTrace/Result_WTrace_Censor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        #generate_density_plot(opt.design, data_path, file.name, N)
      } else{
        print("Saving to:")
        print(paste(results_path,"/NoCensoring/WTrace/Result_WTrace_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
        
        all_best_designs['CensoringFraction'] = 0
        write.csv(all_best_designs, file = paste(results_path,"/NoCensoring/WTrace/Result_WTrace_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        #file.name = paste(results_path,"/NoCensoring/WTrace/Result_WTrace_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        #generate_density_plot(opt.design, data_path, file.name, N)
      }
    } else{
      if(censor){ 
        print("Saving to:")
        print(paste(results_path,"/Censoring/S/Result_S_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
        
        all_best_designs['CensoringFraction'] = censor_frac
        write.csv(all_best_designs, file = paste(results_path,"/Censoring/S/Result_S_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        #file.name = paste(results_path,"/Censoring/S/Result_S_Censor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        #generate_density_plot(opt.design, data_path, file.name, N)
      } else{
        print("Saving to:")
        print(paste(results_path,"/NoCensoring/S/Result_S_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
        
        all_best_designs['CensoringFraction'] = 0
        write.csv(all_best_designs, file = paste(results_path,"/NoCensoring/S/Result_S_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        #file.name = paste(results_path,"/NoCensoring/S/Result_S_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        #generate_density_plot(opt.design, data_path, file.name, N)
      }
      
    }
  }
}

# to save results at the end of every iterations of the competing risk method
# takes all_best_designs which has info about all best designs found so far
# other arguments are the length of a design, method used, where results should be saved,
# where the design data can be found and the value of N (no of Weibull samples),
# whether censoring is applied or not
save_comp_results <- function(candidate, all_best_designs, design_length, method, min.matrix, results_path, N, censor, censor_frac){
  all_best_designs = all_best_designs
  colnames(all_best_designs)[ncol(all_best_designs)] <- "Iter"
  
  # for plotting density ,, the last row in all_best_designs has the best design
  last_row = as.numeric(tail(all_best_designs, n = 1))
  opt.design = last_row[1:design_length]
  n = length(opt.design)
  #opt.design=c(1, 10, 1, 10, 1, 10, 1, 1, 10, 1, 10, 1, 10, 1, 10, 1, 10, 10, 1, 1  )
  thetas_betas = paste(as.character(theta1),"_",as.character(theta2),"_",as.character(input_beta1),"_",as.character(input_beta2),"_",as.character(n),"_",as.character(N), sep = "")
  if(method == "det"){
    if(censor){ 
      all_best_designs['CensoringFraction'] = censor_frac
      cat("Saving results to",paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""),"\n")
      write.csv(all_best_designs, file = paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
      # create and save density plot plot
      file.name = paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
      generate_comp_density_plot(candidate, opt.design, min.matrix, file.name, N)
    } else{
      all_best_designs['CensoringFraction'] = 0
      cat("Saving results to",paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""),"\n")
      write.csv(all_best_designs, file = paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
      file.name = paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
      generate_comp_density_plot(candidate, opt.design, min.matrix, file.name, N)
    }
  } else {
    if(method == "trace"){
      if(censor){ 
        all_best_designs['CensoringFraction'] = censor_frac
        cat("Saving results to",paste(results_path,"/Censoring/WTrace/Result_WTrace_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""),"\n")
        write.csv(all_best_designs, file = paste(results_path,"/Censoring/WTrace/Result_WTrace_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        file.name = paste(results_path,"/Censoring/WTrace/Result_WTrace_Censor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        generate_comp_density_plot(candidate, opt.design, min.matrix, file.name, N)
      } else{
        all_best_designs['CensoringFraction'] = 0
        cat("Saving results to",paste(results_path,"/NoCensoring/WTrace/Result_WTrace_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""),"\n")
        write.csv(all_best_designs, file = paste(results_path,"/NoCensoring/WTrace/Result_WTrace_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        file.name = paste(results_path,"/NoCensoring/WTrace/Result_WTrace_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        generate_comp_density_plot(candidate, opt.design, min.matrix, file.name, N)
      }
    } else{
      if(censor){ 
        all_best_designs['CensoringFraction'] = censor_frac
        cat("Saving results to",paste(results_path,"/Censoring/S/Result_S_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""),"\n")
        write.csv(all_best_designs, file = paste(results_path,"/Censoring/S/Result_S_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        file.name = paste(results_path,"/Censoring/S/Result_S_Censor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        generate_comp_density_plot(candidate, opt.design, min.matrix, file.name, N)
      } else{
        all_best_designs['CensoringFraction'] = 0
        cat("Saving results to",paste(results_path,"/NoCensoring/S/Result_S_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""),"\n")
        write.csv(all_best_designs, file = paste(results_path,"/NoCensoring/S/Result_S_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""), row.names = F)
        file.name = paste(results_path,"/NoCensoring/S/Result_S_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.png", sep = "")
        generate_comp_density_plot(candidate, opt.design, min.matrix, file.name, N)
      }
      
    }
  }
}

