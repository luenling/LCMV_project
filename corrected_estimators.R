#testing bias correction
# number of ind in pool
n = 1000
# read coverage
M = 1000
# minimal read coverage
b = 2

# correction for theta

M= 100

b = 2

# theta_py inner sum
tp_inner <- function(m,M,n) {
  # takes running index of m, coverage M and pool size n
  r=seq(1,n-1,by=1)
  PXm_Yr=dbinom(m,M,r/n)
  res = sum(m*(M-m)*PXm_Yr/r)
  return(res)
  
}

tw_corr <- function(b,M,n){
  # takes running minimal count of b, coverage M and pool size n
  cm = sum(1/seq(1,M-1,by=1))
  r=seq(1,n-1,by=1)
  return(cm/ sum( (pbinom(M-b,M,r/n)-pbinom(b-1,M,r/n))/r ) )
}

n=1000
M=1500
b=10

theta_corr = choose(M,2) / sum(vapply(seq(b,M-b,by=1), tp_inner, 1, M = M, n = n))
tw_corr(b,M,n)




