#if(!require(random)){
#  install.packages("random", repos = "http://cran.us.r-project.org")
#}

library(random)
## open the script that contains functions to compute the Det and Trace
## Change the path to where it is.

# to compute the code runtime
start_time <- Sys.time()

# read commandline args
args=commandArgs(trailingOnly = T)

if(length(args) != 12){
	stop("Please pass valid values for theta01, theta02, theta1, theta2, beta1, beta2, C (T or F), (D or T or S), rep number, n, N and Y")
}


#data_path = "C:/Users/han-s/Downloads/SimDataPH/"
theta01 = as.numeric(args[1])
theta02 = as.numeric(args[2])
theta1 = as.numeric(args[3])
theta2 = as.numeric(args[4])
input_beta1 = as.numeric(args[5])
input_beta2 = as.numeric(args[6])
# Here we decide if we wish to apply censoring
censor = as.logical(args[7])
opt.type = args[8]
method = ""
if(opt.type == "D"){
  method = "det"
} else{
  if(opt.type == "T"){
    method = "trace"
  } else {
    if(opt.type == "S"){
      method = "ds"
    }
  }
}

rep_num = as.numeric(args[9])
n = as.numeric(args[10])# design size
N = as.numeric(args[11]) # number of dataset samples
# If we want censoring, what is the censoring time Y
Y = as.numeric(args[12])


# we choose a random initial design with a specific percentage 
# of 1's and 10's
#n = 50
candidate = c(1, 10, 2.5, 4, 5.5, 7, 8.5)
candidate_length = length(candidate)
cat("Candidate length = ", candidate_length, "\n")
# use a vector of probabilities to generate a design with values in candidate
# we give 1 40%, 10 40%, then distribute the remaining 20% on other values equally
# that's why we subtract 2 from candidate length
prob = c(0.4,0.4, rep(0.2/(candidate_length-2), (candidate_length-2)) )

result = rep(candidate, round(n * prob))

# scramble the values in result
initial_design = sample(result)
diff = n - length(initial_design)
if(diff > 0){
  ss = sample(candidate, diff)
  initial_design = c(initial_design, ss)
}
if(diff < 0){
  diff = abs(diff)
  initial_design = initial_design[1:length(initial_design) - diff]
}


# change path on the cluster
#source("C:/Users/han-s/Downloads/my R scripts/Utils_last.R")
source("/users/k1814473/Utils_last.R")
#source("/Users/csstnns/Downloads/Utils_last.R")

# change path on the cluster
data1_dir_name = paste(as.character(theta1),"_",
                       as.character(input_beta1),"_",
                       as.character(censor),"_",
                       method,"_",as.character(rep_num),
                       "_",as.character(n),
                       "_",as.character(N),
                       "/", sep="")
data2_dir_name = paste(as.character(theta2),"_",
                       as.character(input_beta2),"_",
                       as.character(censor),"_",
                       method,"_",as.character(rep_num),
                       "_",as.character(n),
                       "_",as.character(N),
                       "/", sep="")

#data1_path = paste("/Users/csstnns/Downloads/SampleData/data_",data1_dir_name, sep = "")
#data2_path = paste("/Users/csstnns/Downloads/SampleData/data_",data2_dir_name, sep = "")
data1_path = paste("/scratch/users/k1814473/data_",data1_dir_name, sep = "")
data2_path = paste("/scratch/users/k1814473/data_",data2_dir_name, sep = "")
#data1_path = paste("C:/Users/han-s/Downloads/data_",data1_dir_name, sep = "")
#data2_path = paste("C:/Users/han-s/Downloads/data_",data2_dir_name, sep = "")

cat("Data path 1 = ",data1_path,"\n")
cat("Data path 2 = ",data2_path,"\n")

#results_path = "/Users/csstnns/Downloads"
results_path = "/users/k1814473"
#results_path = "C:/Users/han-s/Downloads/my R scripts/"


# if data path is not there, create it!
#if(!dir.exists(data_path)){
#  dir.create(data_path)
#}

#for(theta1 in c(-0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3)){


# check if folder already has csv files, if so, delete all of them
#temp = list.files(data_path, full.names = T, pattern="*.csv")
#rm = file.remove(temp)

# Generate random data from rweibull dist based on the given params
print('Generating random data.')

## generate two matrices full of data from weibull dist
data.matrix1 = generate_rweibull_data_matrix(initial_design, candidate, data1_path, N, theta01, theta1, input_beta1,censor, method, rep_num)
data.matrix2 = generate_rweibull_data_matrix(initial_design, candidate, data2_path, N, theta02, theta2, input_beta2,censor, method, rep_num)


# compute the min and indicator matrices which are needed to compute the cov matrix
out = find_min_indicator_matrices(data.matrix1, data.matrix2, censor = censor, Y = Y)
min.matrix = out$min.matrix# the min matrix
cat('Min Matrix Dims = ', dim(min.matrix), '\n')
indicator.matrix = out$indicator.matrix# the min indicator matrix
c.indicator.matrix = out$c.indicator.matrix# the censoring indicator matrix
censor_frac = out$censor_frac

#cat('Theta1 = ', theta1, '\n')

#initial_design = as.numeric(randomNumbers(n = 20, min = 1, max = 10,col = 1))
#initial_design = c(10,1,10,1,10,1,1,10,1,10,10,1,10,1,10,10,1,2,1,10)
design_length <- length(initial_design)  # Number of values in each simulated dataset
cat('>> Initial Design: ',initial_design,' Length = ',design_length,'\n')


## create a DF to keep track of best designs the their corresponding det/trace values
all_designs = data.frame()
all_crit_values = data.frame(crit_value=numeric(0))
all_rounds = data.frame(round=numeric(0))

dweights = c(0,0,0,0,0,0)
find_best_design <- function(des, method, pre_crit = -1){
  init_start_time <- Sys.time()
  # we assume the initial design to be the best design so far
  best.des = des
  
  if(pre_crit != -1){
    crit_value = pre_crit
  } else {
    # we compute the det/trace of the initial design
    #cov_mx = find_cov_matrix2(design.data, des, censored = censor, Y = Y)
    # est_mx = find_estimate_matrix(design.data, des, censored = censor, Y = Y)
    # to compute the code runtime
    #start_timex <- Sys.time()
    cov_mx = compute_competing_risk_cov_matrix(min.matrix, indicator.matrix, c.indicator.matrix, des, N, candidate_length, censor = censor)
    # to see how long the code took to finish
    #end_time = Sys.time()
    #run_time <- difftime(end_time,start_timex, units = "hours")[[1]]
    
    #cat('Total time taken to compute cov matrix in hours: ', run_time,'\n')
    
    #print(cov_mx)
    
    crit_value = NULL
    if(method == 'det'){
      #old method
      #crit_value = find_det(cov_mx)
      #crit_value = mean(crit_value[,1], na.rm = T)
      crit_value = find_det_estimate_cov_matrix(cov_mx)
    } else {
      if(method == 'trace'){
        # Load the D-Opt design from the corresponding DOpt experiment
        # That experiment must be run first
        # We use this D-Opt design to compute the weights needed to find crit_value 
        thetas_betas =paste(as.character(theta1),"_",as.character(theta2),"_",as.character(input_beta1),"_",as.character(input_beta2),"_",as.character(n),"_",as.character(N), sep = "")
        if(censor){ 
          #print('Load File!')
          #print(paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = ""))
          file = paste(results_path,"/Censoring/Det/Result_Det_Censor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = "")
          dff = read.csv(file)
          last_row = as.numeric(tail(dff, n = 1))
          dopt.design = last_row[1:length(initial_design)]
          
        } else{
          file = paste(results_path,"/NoCensoring/Det/Result_Det_NoCensor_",thetas_betas,"_",as.character(rep_num),"_results.csv", sep = "")
          dff = read.csv(file)
          last_row = as.numeric(tail(dff, n = 1))
          dopt.design = last_row[1:length(initial_design)]
        }
        cat('Corresponding D-Opt Design: ',dopt.design,'\n')
        # Now we compute the cov matrix of dopt.design 
        dopt.design.cov_mx = compute_competing_risk_cov_matrix(min.matrix, indicator.matrix, c.indicator.matrix, dopt.design, N, candidate_length, censor = censor)
        # compute six weights from the cov matrix
        dweights <<- find_weight_vector6(dopt.design.cov_mx)
        ## now we use the cov matrix of the current design (we've already comupted the matrix above)
        # we use it to compute the crit_value
        crit_value = find_a_opt_weighted_crit_value6(dweights, cov_mx)
      } else{
        if(method == 'ds'){
          crit_value = find_competing_ds_crit(cov_mx)
        } else {
          print('Incorrect method to select best design!')
          stop()
        }
      }
    }
  }
  
  init_finish_time <- Sys.time()
  init_crit_time <- difftime(init_finish_time, init_start_time, units = "mins")[[1]]
  cat('Initial Crit value computation step took ', init_crit_time,' minutes!\n')
  
  #det = find_det(cov_mx)
  print(crit_value)
  # we assume the det/trace of the initial design is the smallest
  min_crit_value = crit_value
  #min_crit_value = mean(crit_value[,1], na.rm = T)
  
  # we need the unique values in the initial design to use when when we do the swap
  #candidate_points= unique(des)
  # here we use a fixed range of candidate points
  # instead of extracting them from the design 
  #candidate_points = c(1,2,3,4,5,6,7,8,9,10)
  
  
  ## we store the very first design in all_designs
  all_designs <<- as.data.frame(rbind(all_designs, best.des))
  # and we store its corresponding det/trace in all_crit_values
  all_crit_values <<- as.data.frame(rbind(all_crit_values, min_crit_value))
  ## round 0 means it's the very beginning of the search as we start
  # with the initial design
  all_rounds <<- as.data.frame(rbind(all_rounds, 0))
  #print(2)
  # now we loop thru the design one component at a time
  for(idx in 1:length(des)){
    # loop thru the unique values of the initial design
    round = 1
    for(u in candidate){
      # keep track of the current component because we need to return it to the
      # design in case the design is not better than the one before
      cc = best.des[idx]
      # swap the ith component of the design with the current unique value
      best.des[idx] = u
      #cat('Trying design: ',best.des,' component',idx,' - round = ',round,'\n')
      
      # compute the det/trace and see if it's smaller than the current smallest value
      if(length(unique(best.des)) != 1){
        #cat("Load data in Round: ",round, ' -----\n')
        #design.data = load_design_data(best.des, data.matrix, N)
        #cat("Find Cov Mx in Round: ",round, ' -----\n')
        #print((design.data))
        #print(best.des)
        #tmp.cov_mx = find_cov_matrix2(design.data, best.des, censored = censor, Y = Y)
        tmp.cov_mx = compute_competing_risk_cov_matrix(min.matrix, indicator.matrix, c.indicator.matrix, best.des, N, candidate_length, censor = censor)
        
        #tmp.est_mx = find_estimate_matrix(design.data, best.des, censored = censor, Y = Y)
        #cat('Mean theta1 for Current Design = ',abs(mean(tmp.cov_mx[3,], na.rm = T)),'\n')
        #cat('Mean Exp Beta for Current Design = ',(mean(exp(tmp.cov_mx[1,]), na.rm = T)),'\n')
        #cat('Mean theta0 for Current Design = ',(mean(exp(tmp.cov_mx[2,]), na.rm = T)),'\n')
        #diff = mean(exp(tmp.cov_mx[2,]), na.rm = T) - 5.5 * (mean(tmp.cov_mx[3,], na.rm = T))
        #cat("Diff = ", diff,'\n')
        #cat("Find det in Round: ",round, ' -----\n')
        tmp.crit_value = NULL
        if(method == 'det'){
          #old method
          #tmp.crit_value = find_det(tmp.cov_mx)
          #tmp.crit_value = mean(tmp.crit_value[,1], na.rm = T)
          tmp.crit_value = find_det_estimate_cov_matrix(tmp.cov_mx)
        } else {
          if(method == 'trace'){
            # Now we compute the cov matrix of dopt.design 
            #dopt.design.data = load_design_data(dopt.design, data.matrix, N) 
            #dopt.design.est_mx = find_estimate_matrix(dopt.design.data, dopt.design, censored = censor, Y = Y)
            #dopt.design.cov_mx = cov(dopt.design.est_mx)
            
            tmp.crit_value = find_a_opt_weighted_crit_value6(dweights, tmp.cov_mx)
            
            
            # Use the weight values loaded at the start of the procedure
            #tmp.crit_value = find_weighted_avg(tmp.est_mx, weights)
            #tmp.crit_value = find_trace(tmp.cov_mx)            
          } else{
            if(method == 'ds'){
              tmp.crit_value = find_competing_ds_crit(tmp.cov_mx)
            } else {
              print('Incorrect method to select best design!')
              stop()
            }
          }
        }
        #tmp.crit_value = find_det(tmp.cov_mx)
        #tmp.crit_value = mean(tmp.crit_value[,1], na.rm = T)
        # if it's indeed smaller
        # then save this design into the table of best designs (all_designs)
        # also save the corresponding det/trace
        if(tmp.crit_value < min_crit_value){
          #print('Found smaller')
          #best.des[idx] = u
          min_crit_value = tmp.crit_value
          #cat(best.des, min_det,'\n')
          # if the currently found best design is in our stores for memoization
          # then remove it and replace it with the new computed det/ trace in all_crit_values
          # this is because the current det/trace is less than the previous one!
          #cat('New best design found:',best.des, ' with det value ',min_crit_value,'\n')
          if(length(which(colSums(t(all_designs) == best.des) == ncol(all_designs)))>0){
            rowids = which(colSums(t(all_designs) == best.des) == ncol(all_designs))
            cat('To Drop: ',rowids,'\n')
            #print(rowids)
            all_designs <<- as.data.frame(all_designs[-rowids,])
            all_crit_values <<- as.data.frame(all_crit_values[-rowids,])
            all_rounds <<- as.data.frame(all_rounds[-rowids,])
            cat('-----','Drop design from all_designs: ',best.des,' -----\n')
            # reset row IDs to make sure they are 1,2,3,... 
            if(nrow(all_designs)!=0){
              rownames(all_designs) = 1:nrow(all_designs)
            }
            if(nrow(all_crit_values)!=0){
              rownames(all_crit_values) = 1:nrow(all_crit_values)
            }
            if(nrow(all_rounds)!=0){
              rownames(all_rounds) = 1:nrow(all_rounds)
            }
          }
          ## for memoization
          ## here we store the current best design in all_designs
          ## we store the current best design in all_designs
          all_designs <<- as.data.frame(rbind(all_designs, best.des))
          # and we store its corresponding det/trace in all_crit_values
          all_crit_values <<- as.data.frame(rbind(all_crit_values, min_crit_value))
          # store the round number
          all_rounds <<- as.data.frame(rbind(all_rounds, round))
          
        } else{
          ## here we return the component value in case the desgin is not good
          best.des[idx] = cc
        }
      }
      #else {
      #  print('All values in current Design are the same and this is problematic for coxph.wtest.')
      #}
      #cat("Round: ",round, ' -----\n')
      round = round + 1
    }
    #print(3)
    cat('Best design after finishing component',idx, 'is', best.des,' - Crit Value:',min_crit_value,'\n')
    #print(paste('-----','Best design so far: ',best.des,'\n-----\n'))
    #print(comp)
  }
  
  all_det_designs = cbind(all_designs, all_crit_values, all_rounds)
  varnames <- paste('C', 1:ncol(all_det_designs), sep = "_")
  #print(4)
  colnames(all_det_designs) = varnames
  if(method == 'det'){
    colnames(all_det_designs)[ncol(all_det_designs)-1] <- "Determinant"
  } else {
    colnames(all_det_designs)[ncol(all_det_designs)-1] <- "Trace"
  }
  #print(5)
  colnames(all_det_designs)[ncol(all_det_designs)] <- "Round"
  
  # return the best_trace_design and best_det_design
  return(all_det_designs)
}
# max num of iterations to perform
# plan b exist strategy in case the tolerance/patience are not satisfied!
num_iterations = 30
# tol is the min acceptable diff between two consecutive crit values in results
tol = 0.0000001
# patience is the number of consecutive crit values to consider
patience = 20
#initial_design = c(3,5,7,8,9,6,4,8,7,9,8,10,1,3,6,7)
#initial_design = c(3,5,3,8,9,6,4,7,2)
all_best_designs = data.frame()# this DF saves results from all iterations

## this will be used to keep track of the crit value of a design
# when we use the best design of an iter as the initial design in another iter
# no need to recompute the crit value!
pre_crit = -1
for(iter in 1:num_iterations){
  # to compute the code runtime
  iter_start_time <- Sys.time()
  
  ## get results of current iter
  best_designs = find_best_design(initial_design,method,pre_crit)
  # this vector has the iter number to add it as a column in all_best_designs
  x = rep(iter,nrow(best_designs))
  ## update all_best_designs by adding results of current iter and iter number
  all_best_designs <<- as.data.frame(rbind(all_best_designs, cbind(best_designs,x)))
  # the last row in best_designs has the best design in the current iter
  last_row = as.numeric(tail(best_designs, n = 1))
  cat('Iter',iter, ' - Initial Design: ',initial_design,'\n')
  cat('Best Design: ',last_row[1:length(initial_design)],'\n')
  # take the best design in the current iter as the initial design for the next iter
  initial_design = last_row[1:length(initial_design)]
  # store the corresponding crit value of the best design in the current iter
  pre_crit = last_row[length(initial_design)+1][1]
  #cat('Iter',iter, ' - Crit of Initial Design: ',pre_crit,'\n')
  ## reset these two tables
  all_designs = data.frame()
  all_crit_values = data.frame(crit_value=numeric(0))
  all_rounds = data.frame(round=numeric(0))
  
  # Here we look at whether the crit value has stopped improving
  # get the column that has the crit value (the third column from the end)
  crit_col = ifelse(method == "det", all_best_designs['Determinant'], all_best_designs['Trace'])
  # compute the difference between consec values at the end of this col
  # number of values = patience
  # get the abs value of differences and count how many are < the tolerance
  # if this number equals patience, then there is no improvement for patience
  # number of times
  num_diff = sum(tail(abs(diff(crit_col[[1]])), patience) < tol)
  if(num_diff == patience){
    cat('Exiting because of no Crit value improvement!. Iter = ', iter, '\n')
    # save results before existing
    save_comp_results(candidate, all_best_designs, design_length, method, min.matrix, results_path, N, censor, censor_frac)
    break
  }
  
  # save results at the end of every iteration
  save_comp_results(candidate,all_best_designs, design_length, method, min.matrix,results_path, N, censor, censor_frac)
  
  iter_finish_time <- Sys.time()
  iter_time <- difftime(iter_finish_time, iter_start_time, units = "mins")[[1]]
  cat('Iter ', iter, ' took ', iter_time,' minutes!\n')
}

# to see how long the code took to finish
end_time = Sys.time()
run_time <- difftime(end_time,start_time, units = "hours")[[1]]

cat('Total time taken to run in hours: ', run_time,'\n')


