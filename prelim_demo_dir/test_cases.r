####
### Test Cases for R -> C++ Conversion
### Checks the output of the C++ to the R function
### Since floats, use tolerance to check.
####

library(RcppArmadillo)
library(Rcpp)
library(VGAM)


#Reads in the source file and does compilation if necessary...
Rcpp::sourceCpp("smp_clib.cpp",
                showOutput = TRUE,
                verbose=TRUE,
                rebuild = FALSE)

n <- 4
temp_s <- 0
Id <- diag(1,n)
source("smp_lib.r")
.Call("createIdentityMatrix", 4)
pft <- p * ft(1.0)


#### Test Case 1 ####
orig_ans <- firstpasft(1.0)
c_ans <- .Call('firstpasft', pft)
isEq <- abs(orig_ans - c_ans) < 0.00001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 2 ####
#Note: the same value passed to firstpasFt should
# be the same second arg sent to the C version.

orig_ans <- firstpasFt(1.0)
c_ans <- .Call('firstpasFt', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.00001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 3 ####
#Note: the same value passed to transprobt should
# be the same second arg sent to the C version.
# TODO: Look into accuracy issues... the 0's, not close to 0...

orig_ans <- transprobt(1.0)
c_ans <- .Call('transprobt', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.01
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 4 ####

orig_ans <- SMP.firstpass.cdf(1.0)
c_ans <- .Call('smp_firstpass_cdf', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.0001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 5 ####

orig_ans <- SMP.firstpass.pdf(1.0)
c_ans <- .Call('smp_firstpass_pdf', pft)
isEq <- abs(orig_ans - c_ans) < 0.0001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 6 ####

orig_ans <- SMP.transprob.t(1.0)
c_ans <- .Call('smp_transprob_t', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.01
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 7 ####

orig_ans <- SMP.expected.visits.to.state(1.0)
c_ans <- .Call('smp_expected_visits_to_state', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.0001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 8 ####
#TODO: Look into tolerance issue...

orig_ans <- SMP.expected.time.in.state(1.0)
c_ans <- .Call('smp_expected_time_in_state', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.01
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 9 ####
#TODO: Look into tolerance issue...

orig_ans <- SMP.prob.trans.to.state.0.times(1.0)
c_ans <- .Call('smp_trans_to_state_0_times', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.0001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))


#### Test Case 10 ####
#TODO: Look into tolerance issue...

orig_ans <- SMP.prob.trans.to.state.1.times(1.0)
c_ans <- .Call('smp_trans_to_state_1_times', pft, 1.0)
isEq <- abs(orig_ans - c_ans) < 0.0001
print(ifelse(all(isEq, TRUE), "PASS", "FAIL"))



################## OUTPUT ######################
PASS
PASS
PASS
PASS
PASS
PASS
PASS
PASS
PASS
PASS



