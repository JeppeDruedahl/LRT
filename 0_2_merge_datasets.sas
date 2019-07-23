libname data "K:/Workdata/706248/LRT/data";

**********;
* 1. BEF *;
**********;

* static verson: all unique pnr (sex and birthyear is constant over year);
proc sql;

	create table data.bef_static
	as select

		/* join-variables */
		distinct pnr as pnr label="pnr",

		/* content-variables */
		year(FOED_DAG) as birthyear label="birthyear",
		KOEN as sex label="sex"

	from

		data.bef

	group by pnr
	order by pnr

;quit;


***********;
* 2. FAIN *;
***********;

* static verson: all unique pnr (sex and birthyear is constant over year);
proc sql;

	create table
		data.fain_static

	as 
	select

		/* join-variables */
		distinct pnr as pnr label="pnr", 

		/* content-variables */
		KOEN as fain_sex label="fain_sex",
		fainSourceYear-ALDER-1 as fain_birthyear label="fain_birthyear"

	from

		data.fain
	
	group by pnr
	order by pnr

;quit;

**********;
* 3. IND *;
**********;

data data.ind_temp(keep = y_disp a h pnr indupdSourceYear);
set data.ind;	
	
	QRENTUD2 = ERHVERVSINDK_GL + OVERFORSINDK + FORMUEINDK_NY + RESUINK_GL - PERSAMLINKNETRENT_NY;
	
	y_disp = DISPON_13 - FORMUEINDK_BRUTTO - LEJEV_EGEN_BOLIG + RENTUDGPR; 
	if missing(y_disp) then y_disp = DISPON_NY - FORMUEINDK_NY + QRENTUD2;
	
	a = FORMREST_NY05;
	if missing(a) then a = FORM;

	h = EJENDOMSVURDERING;
	if missing(h) then h = KOEJD;

run;


****************;	
* 4. full data *;
****************;	

proc sql;

	create table
		data.full_temp

	as 
	select

		ind.pnr,
		ind.y_disp,
		ind.a,
		ind.h,
		ind.indupdSourceYear as year,
		bef_static.birthYear,
		bef_static.sex,
		fain_static.fain_birthyear,
		fain_static.fain_sex,
		idap.PSTILL

	from

		data.ind_temp as ind

		left join
		data.bef_static as bef_static
			on 
			bef_static.pnr 	= ind.pnr			

		left join
		data.fain_static as fain_static
			on 
			fain_static.pnr	= ind.pnr

		left join
		data.idap as idap
			on 
			idap.pnr = ind.pnr and
			idap.idapSourceYear = ind.indupdSourceYear

	order by pnr, year

;quit;

* define: replace age and sex if missing;
data data.full(drop = fain_birthyear fain_sex);
set data.full_temp;	

	if missing(birthyear)
	then birthyear = fain_birthyear;
	age = year-birthyear;

	if missing(sex)
	then sex = fain_sex;

run;
