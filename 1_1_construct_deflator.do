* Description: create deflator data set.
* Output: data\deflator.dta.

set more off
set type double

cd K:\workdata\706248\LRT\

*************
** 1. load **
*************

import excel "data\PRIS61.xlsx", sheet("PRIS61") cellrange(A4:B38) clear
rename A year
rename B deflator_base


****************************
** 2. rebase to year 2014 **
****************************

sort year
destring year, replace
recast int year

egen double value_2014 = max(deflator_base)

gen double deflator = deflator_base / value_2014
drop deflator_base value_2014


*************
** 3. Save **
*************

sort year
save data\deflator.dta, replace
