set more off
set type double
set graphics off
graph drop _all
set scheme lean2

cd K:\workdata\706248\LRT\

* settings
local pop_loc "_p100" 
*foreach pop_loc in "_p1" "_p100" {

	global pop `pop_loc'

	
***********
* 1. load *
***********

use data\sample${pop}.dta, clear

cap drop rdraw
cap drop keepprob
sort pnr year 
rename y_disp Y


***********************
* 2. cohort selection *
***********************

gen drop_this = 0
label define drop_this_lbl 0 "included" 1 "not 30 obs in range" 2 "non-positive y" 3 "self-employed > 3 yrs"
label values drop_this drop_this_lbl

keep if birthyear >= 1950 & birthyear <= 1955 	
by pnr: gen pnr_N = _N
by pnr: gen pnr_n = _n


*********************
* 3. main selection *
*********************

* a. 30 obs in range
gen in_range = age >= 30 & age <= 59
by pnr: egen obs_in_range = sum(cond(in_range == 1,1,0))

* b. non-positive y
by pnr: egen non_pos_y = sum(cond(Y <= 10,1,0))

* c. self-employment
by pnr: egen selfemployed = sum(cond(pstill < 30,1,0))

* d. drop variable
replace drop_this = 1 if obs_in_range != 30
replace drop_this = 2 if non_pos_y > 0 & drop_this == 0
replace drop_this = 3 if selfemployed > 3 & drop_this == 0

* e. table
tab drop_this if pnr_n == 1, matcell(freq) matrow(names)
	
	putexcel set "output\units${pop}.xls", replace
	putexcel A1=matrix(names) B1=matrix(freq) C1=matrix(freq/r(N))
	
tab drop_this if in_range == 1, matcell(freq) matrow(names)

	putexcel set "output\totobs${pop}.xls", replace
	putexcel A1=matrix(names) B1=matrix(freq) C1=matrix(freq/r(N))

* f. total
gen do_drop = in_range == 0 | drop_this > 0


**********************
* 6. cohort-effects * 
**********************
	
* a. log income 
gen logY_raw = log(Y) 
	
* b. subtract cohort-effects 
su birthyear 
global min_c `r(min)' 
global max_c `r(max)' 
reg logY_raw ib(${min_c}).birthyear i.age if do_drop == 0 

gen cohort_effect = . 
replace cohort_effect = 0 if birthyear==${min_c} 
local low = ${min_c}+1
forvalues c=`low' / $max_c { 
	replace cohort_effect = _b[`c'.birthyear] if birthyear==`c' 
} 
gen logY = logY_raw - cohort_effect

* c. net worth
sort year
merge m:1 year using data\deflator_houseprice.dta // adds variable: hp_correction
drop if _merge == 2
drop _merge

gen h_correction = 0
replace h_correction = h/hp_correction - h if !missing(hp_correction)
replace a = (a+h_correction)/exp(cohort_effect)

sort pnr year


*****************
* 7. unbalanced * 
*****************

preserve

* a. drop
drop if (in_range == 0) | (non_pos_y > 0) | (selfemployed > 3)
by pnr: replace pnr_n = _n

* b. table
gen temp = 1
tab temp if pnr_n == 1, matcell(freq) matrow(names)
	
	putexcel set "output\units_unbalanced${pop}.xls", replace
	putexcel A1=matrix(names) B1=matrix(freq) C1=matrix(freq/r(N))
	
tab temp, matcell(freq) matrow(names)

	putexcel set "output\totobs_unbalanced${pop}.xls", replace
	putexcel A1=matrix(names) B1=matrix(freq) C1=matrix(freq/r(N))
	
* c. export
keep pnr age logY 

su age
global min_x `r(min)' 
global max_x `r(max)' 

reshape wide logY, i(pnr) j(age) 
outfile logY${min_x}-logY${max_x} using "data/logY_unbalanced${pop}.csv", comma wide replace  

restore


******************
* 8. zero cut-off * 
*******************

preserve

* a. drop
drop if (in_range == 0) | (selfemployed > 3) | (obs_in_range != 30)
by pnr: replace pnr_n = _n

* b. table
gen temp = 1
tab temp if pnr_n == 1, matcell(freq) matrow(names)
	
	putexcel set "output\units_zerocutoff${pop}.xls", replace
	putexcel A1=matrix(names) B1=matrix(freq) C1=matrix(freq/r(N))
	
tab temp, matcell(freq) matrow(names)

	putexcel set "output\totobs_zerocutoff${pop}.xls", replace
	putexcel A1=matrix(names) B1=matrix(freq) C1=matrix(freq/r(N))
	
* c. export
keep pnr age Y 

su age
global min_x `r(min)' 
global max_x `r(max)' 

reshape wide Y, i(pnr) j(age) 
outfile Y${min_x}-Y${max_x} using "data/Y_zerocutoff${pop}.csv", comma wide replace

restore


****************
* 10. baseline * 
****************

drop if do_drop == 1
by pnr: replace pnr_n = _n

* a. net worth
preserve
collapse (p50) a, by(age)	
outfile a using "data/a_p50${pop}.csv", comma wide replace  
restore
	
preserve
collapse (mean) a, by(age)	
outfile a using "data/a_mean${pop}.csv", comma wide replace  
restore

* b. birthyears
outfile birthyear using "data/birthyear${pop}.csv" if pnr_n == 1, comma wide replace 

* c. income
keep pnr age logY logY_raw Y

su age
global min_x `r(min)' 
global max_x `r(max)' 

preserve 
reshape wide logY logY_raw Y, i(pnr) j(age) 
outfile Y?? using "data/Y${pop}.csv", comma wide replace
outfile logY?? using "data/logY${pop}.csv", comma wide replace
outfile logY_raw?? using "data/logY_raw${pop}.csv", comma wide replace
restore 

* end loops
* }
