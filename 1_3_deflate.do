set more off
set type double
graph drop _all

cd K:\workdata\706248\LRT\

***********************
* 1. load full-sample *
***********************

use data\sample_raw.dta, clear
sort pnr year 

label var age "age"
label var pstill "pstill"

* a. pnr as int
destring pnr, replace
egen newpnr = group(pnr)
drop pnr
rename newpnr pnr
label var pnr "pnr"

* b. year as int
recast int year
sort pnr year

* c. pstill 
destring pstill, replace 

	
**************
* 2. deflate *
**************

summarize year, meanonly
global min_year = r(min)
global max_year = r(max)
	
sort year
merge m:1 year using data\deflator.dta
drop _merge
drop if year < $min_year
drop if year > $max_year

foreach x in y_disp a h {	
	disp "`x'"
	replace `x' = (`x' / 1000) / deflator
}


*************
** 3. save **
*************

sort pnr year
by pnr year: keep if _n == 1
save data\sample_p100.dta, replace

* random sample
by pnr: gen rdraw = runiform() if _n == 1
by pnr: gen keepprob = rdraw[1]

keep if keepprob < 0.01
save data\sample_p1.dta, replace

drop rdraw keepprob
