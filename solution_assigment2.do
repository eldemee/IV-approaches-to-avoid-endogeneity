log using assignment2.log, replace

clear

* Asymptotic distribution of IV estimates

*Here I set seed to 1000 so that the task can be replicated
set seed 1000

*Generate exogenous regressors and instruments
set obs 2000
gen x1 = runiform()
gen x2 = rnormal()
gen z1 = runiform()*2
gen z2 = rnormal()*2
gen z3 = rnormal()*2
gen z4 = rchi2(2)

*Monte Carlo simulation with different sample sizes
foreach i of numlist 50 200 2000 {

	local j = 1
	tempfile mc`i'

	while `j' <= 1000 {
		
		quietly {
			preserve
			
			*Draw a random sample of size i
			sample `i', count
			
			*Print progress indicator
			noisily di "." _continue
			
			*Generate error terms with mean zero
			gen e = rchi2(2)
			sum
			replace e = e - r(mean)
			
			gen v = rchi2(3)
			sum
			replace v = v - r(mean)
			
			*Generate endogenous variable
			gen y2 = 1 + x1 + x2 + z1 + z2 + 0.15*z3 + 0.15*z4 + e + v
			
			*Generate dependent variable
			gen y1 = 1 + y2 + x1 + x2 + e
			
			*IV regression using z1 and z2 as instruments
			ivregress 2sls y1 x1 x2 (y2 = z1 z2)
			matrix ivb1 = e(b)
			
			*IV regression using z3 and z4 as instruments
			ivregress 2sls y1 x1 x2 (y2 = z3 z4)
			matrix ivb2 = e(b)
			
			*IV regression using all instruments
			ivregress 2sls y1 x1 x2 (y2 = z1 z2 z3 z4)
			matrix ivb3 = e(b)
			
			*Convert coefficient matrices into variables
			svmat ivb1 , name(olsz12`i')
			svmat ivb2 , name(olsz34`i')
			svmat ivb3 , name(olsz_all`i')
			
			*Keep the coefficient values
			collapse (mean) ols*
			
			*Append results across simulations
			if `j' > 1{
				append using `mc`i''
			}
			
			save `mc`i'', replace
			restore
			
			local j = `j' + 1
		}
	}
	
	di "obs = `i' done"
}

*Combine results for all sample sizes
use `mc50', clear
append using `mc200'
append using `mc2000'

set more on
summarize

*Example: distribution of the coefficient of x1
twoway (kdensity olsz12501, xline(1)) (kdensity olsz34501) (kdensity olsz_all501)
graph export iv_N50.png, replace
twoway (kdensity olsz122001, xline(1)) (kdensity olsz342001) (kdensity olsz_all2001)
graph export iv_N200.png, replace
twoway (kdensity olsz1220001, xline(1)) (kdensity olsz3420001) (kdensity olsz_all20001)
graph export iv_N2000.png, replace

/*
The graphs show how the distribution of the estimated coefficient changes depending on the instruments used and the sample size. When weak instruments (z3 and z4) are used, the estimator becomes less precise and the standard deviation increases, leading to a wider distribution. Using stronger instruments (z1 and z2) or all instruments produces estimates that are more concentrated around the true value (β = 1). As the sample size increases, the dispersion decreases and the estimates converge to the true parameter.
*/ 

*Dealing with measurment error

clear all

*Here I set seed to 1000 so that the task can be replicated
set seed 1000

set obs 1000

gen z1 = 2*invnorm(uniform())
gen z2 = 2*invnorm(uniform())

matrix omega=(1, 0.7 \ 0.7, 1)
matrix list omega
matrix L=cholesky(omega)
matrix list L
matrix K=L'
mat list K

gen x1=z1*K[1,1]+z2*K[2,1]
gen x2=z1*K[1,2]+z2*K[2,2]

*z1 and z2 are not correlated
corr z1 z2
*x1 and x2 are correlated
corr x1 x2

*We dont observe for x2 but two proxies that have measurment error
gen x2_e1 = x2 + 0.5*invnorm(uniform())
gen x2_e2 = x2 + invnorm(uniform())

*Our true model is
*y = 1 + x1 + x2 + e
gen e = invnorm(uniform())
gen y = 1 + x1 + x2 + e

reg y x1 x2
reg y x1
reg y x1 x2_e1
reg y x1 x2_e2
reg y x1 x2_e1 x2_e2
ivregress 2sls y x1 (x2_e1 = x2_e2)
ivregress 2sls y x1 (x2_e2 = x2_e1)

*Here I set seed to 1000 so that the task can be replicated

*Monte Carlo again
foreach i of numlist 50 1000 {
	tempfile ivmc`i'
	local j = 1
	while `j' <= 1000 {
		quietly{
			preserve
			sample `i', count
			noisily di "." _continue
			
			*Regenerating errors
			replace e = invnorm(uniform())
			replace y = 1 + x1 + x2 + e
			replace x2_e1 = x2 + 0.5*invnorm(uniform())
			replace x2_e2 = x2 + invnorm(uniform())
			
			ivregress 2sls y x1 (x2_e1 = x2_e2)
			matrix b1 = e(b)
			ivregress 2sls y x1 (x2_e2 = x2_e1)
			matrix b2 = e(b)
			reg y x1 x2_e1 x2_e2
			matrix both = e(b)
			
			svmat b1, name(ols_smaller`i')
			svmat b2, name(ols_larger`i')
			svmat both, name(ols_both`i')
			collapse (mean) ols_smaller* ols_larger* ols_both*
			
			if `j' > 1 {
				append using `ivmc`i''
			}
			save `ivmc`i'', replace
			
			restore
			
			local j = `j' + 1
		 }
	}
}

use `ivmc50', clear
append using `ivmc1000'
set more on
sum

twoway (kdensity ols_smaller501 , xline(1)) (kdensity ols_smaller10001)
graph export iv_proxy_small.png, replace
twoway (kdensity ols_larger501 , xline(1)) (kdensity ols_larger10001)
graph export iv_proxy_large.png, replace
twoway (kdensity ols_both501 , xline(1)) (kdensity ols_both10001)
graph export iv_proxy_both.png, replace
*If we use both proxies bias is not solved

*As we can see, the standard errors of the estimated coefficients are much smaller in the large sample than in the small sample. Moreover, the IV estimator removes the bias present in the OLS estimator caused by measurement error. Finally, when the instrument has larger measurement error, the IV estimator has a larger standard deviation, indicating lower precision.

*Running IV regressions "without instruments"
clear all

*Here I set seed to 1000 so that the task can be replicated
set seed 1000

set obs 3000

gen x1 = invnorm(uniform())

foreach i of numlist 100 3000 {
	tempfile  arthurmc`i'
	local j = 1
	while `j' <= 1000 {
		quietly {
			
			preserve
			sample `i', count
			noisily di "." _continue
	
			gen u = invnorm(uniform())
			gen s1 = invnorm(uniform())
			gen s2 = invnorm(uniform())

			*Heterocedasticity errors
			gen e1 = u + exp(x1)*s1
			gen e2 = u + exp(-x1)*s2

			*Triangular model
			gen y2 = 1 + x1 + e2
			gen y1 = 1 + y2 + x1 + e1

			*First stage regression
			/*
			x1 can be used as Z because satisfies both assumptions
			- Relevance cov(x1,y2) != 0
			- Exogeneity con(x1, u) = 0
			*/
			reg y2 x1
			predict e2hat, resid

			*Create the Lewbel instrument
			sum x1
			gen z = (x1 - r(mean))*e2hat

			*Run IV
			ivregress 2sls y1 x1 (y2 = z)
			matrix b1 = e(b)
			
			reg y1 x1 y2 
			matrix b2 = e(b)
			
			svmat b1, name(ols_consistent`i')
			svmat b2, name(ols_inconsistent`i')
			
			collapse (mean) ols_consistent* ols_inconsistent*
			
			if `j' > 1 {
				append using `arthurmc`i''
			}
			
			save `arthurmc`i'', replace
			restore
			
			local j = `j' + 1
		}
	}	
}

use `arthurmc100', clear
append using `arthurmc3000'
set more on
sum

twoway (kdensity ols_consistent1001, xline(1)) (kdensity ols_inconsistent1001)
graph export arthur100.png, replace

twoway (kdensity ols_consistent30001, xline(1)) (kdensity ols_inconsistent30001) 
graph export arthur3000.png, replace


/*
The IV estimator using the Lewbel instrument is centered around the true parameter β = 1,
showing that the estimator is consistent. The OLS estimator is biased because y2 is
endogenous and correlated with the error term. As the sample size increases, the IV
distribution becomes more concentrated around the true value, confirming consistency.
*/


log close


