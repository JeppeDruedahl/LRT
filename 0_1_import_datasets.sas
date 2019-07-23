libname raw "H:/Rawdata/706248/views/";
libname data "K:/Workdata/706248/LRT/data";

**********;
* 1. IND *;
**********;

proc sql;
	create table
		data.ind
	as 		
	select

		DISPON_NY,
		DISPON_13,

		ERHVERVSINDK_13,
		ERHVERVSINDK_GL,
		
		RENTUDGPR,
		FORMUEINDK_NY,
		FORMUEINDK_BRUTTO,
		LEJEV_EGEN_BOLIG,

		PERSAMLINKNETRENT_NY,
		OVERFORSINDK,
		RESUINK_GL,

		FORM,
		FORMREST_NY05,
		EJENDOMSVURDERING,
		KOEJD,

		pnr,
		indupdSourceYear

	from
		raw.indupdv
;quit;


**********;
* 2. BEF *;
**********;

proc sql;
	create table
		data.bef
	as 		
	select
		KOEN,
		FOED_DAG,
		pnr,
		befSourceYear
	from
		raw.befv
;quit;


**********;
* 3. BEF *;
**********;

proc sql;
	create table
		data.fain
	as 		
	select
		KOEN,
		ALDER,
		pnr,
		fainSourceYear
	from
		raw.fainv
;quit;


***********;
* 3. IDAP *;
***********;

proc sql;
	create table
		data.idap
	as 		
	select
		PSTILL,
		pnr,
		idapSourceYear
	from
		raw.idapv
;quit;
