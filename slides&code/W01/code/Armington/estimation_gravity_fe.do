clear all
set more off

cd "D:\坚果云Lite\Teaching\2025_Teaching\Structural Models and Numerical Methods in Economics\slides_code\w01\code\Armington"

//Data

import excel using "data_gravity\xtr2001.xlsx", firstrow clear
save "temp\xtr2001.dta", replace

import excel using "data_gravity\tariff2001.xlsx", firstrow clear
save "temp\tariff2001.dta", replace

import excel using "data_gravity\distance2001.xlsx", firstrow clear
save "temp\distance2001.dta", replace


//Estimation by fixed-effect estimator
use "temp\xtr2001.dta", clear

merge 1:1 iso_o iso_d using "temp\tariff2001.dta", nogen
merge 1:1 iso_o iso_d using "temp\distance2001.dta", nogen

sum xtr2001

gen log_value = log(xtr2001)
gen log_tar = log(1+tariff)
gen log_dist = log(dist)

reghdfe log_value log_tar log_dist lang sgdp if iso_o != iso_d, absorb(iso_o iso_d) vce(r)

estimates store ma, title(Structural Gravity Equation)
estout ma using "result\estimate_gravity.txt", cells(b(star fmt(3)) se) starlevels( * 0.10 ** 0.05 *** 0.010) stats(r2 N) replace
estimates store clear


//Data for counterfactual

preserve
keep iso_o iso_d xtr2001
reshape wide xtr2001, i(iso_o) j(iso_d)string
export delimited using "data_matlab\xtr2001.csv", replace
restore

preserve
keep iso_o iso_d tariff
sum tariff
reshape wide tariff, i(iso_o) j(iso_d)string
export delimited using "data_matlab\tariff2001.csv", replace
restore




