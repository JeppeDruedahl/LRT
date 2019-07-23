libname data "K:/Workdata/706248/LRT/data";

*****************;
* 1. raw sample *;
*****************;

data data.sample_raw;
set data.full;

if missing(age) = 1 then delete;
if age < 18 then delete;
if age > 70 then delete;

run;


**********************;
* 2. export to STATA *;
**********************;

proc export data=data.sample_raw outfile="K:/Workdata/706248/LRT/data/sample_raw.dta" DBMS=Stata REPLACE;
run;


