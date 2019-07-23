* Description: create house price deflator data set.
* Output: data\deflator_housepirce.dta.

set more off
set type double

cd K:\workdata\706248\LRT\

***********
* 1. load *
***********

import excel "data\EJEN6.xlsx", sheet("EJEN6") cellrange(A4:B27) clear
rename A year
rename B hp

destring year, replace


*****************
* 2. regression *
*****************

gen year_sq = year*year
reg hp year if year < 2005 | year > 2010

gen hp_trend = _b[_cons] + _b[year]*year
gen hp_correction = hp/hp_trend


***********
* 3. save *
***********

keep hp_correction year
sort year

save data\deflator_houseprice.dta, replace
