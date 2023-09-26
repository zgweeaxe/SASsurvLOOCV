%let analysis= ;
%let cv= ;
%let dir= ;
%let cns= ;
%let pat= ;
%let censvr= ;
%let timevr= ;
%let k= ;
%let perc1= ; 
%let perc2= ;
%let studynm= ;
%let survds= ;
%let factords= ;
%let ospfs = ;
%let timepoint= ;

libname analysis "&analysis.";
libname cv "&cv.";


proc sort data= analysis.&factords. out=&factords.;
by &pat. ;
run;

proc sort data=analysis.&survds. out=newdat ;
by &pat. ;
run;


data temp;
merge newdat &factords. (in=ina);
by &pat.;
if ina;
run;
data temp2;
set temp;
myid = _N_;
run;
data _null_;
	set temp2;
	call symput ('sampn',_N_);
run;
%put &sampn;


/* generate the cross-validation sample */
proc surveyselect data=temp2 out=xv seed=789556
samprate=1 outall rep=&sampn;
run;


data xv2; set xv;
if replicate=myid then selected=0;
run;


data xv_select;
set xv2;
where selected=1;
if selected=1 then newprogmos=&timevr.;newcens=&censvr.;
drop &timevr. &censvr.;
run;

proc sort data= xv_select out=cv.xv_select;
by &pat. ;
run;


/** creates variable list for cox proportional hazard models**/
proc contents data=&factords. (drop= &pat.) 
 out=vars (keep= varnum name) 
 noprint; 
run; 
proc sql noprint ; 
select name into :vlist separated by ' ' 
from vars 
order by varnum; 
quit; 

%put &vlist.;

data _null_;
call symput('fctrs',strip("&&vlist"));
run;
%put &fctrs;

data ptlist;
set cv.xv_select(keep=&pat.);
by &pat.;
if last.&pat.;
run;

data &factords._anal;
	merge &factords. ptlist(in=inp);
	if inp;
run;

data xv_noselect;
set xv2;
where selected=0;
if selected=0 then newprogmos=&timevr.;newcens=&censvr.;
drop &timevr. &censvr.;
run;

proc sort data= xv_noselect out=cv.xv_noselect;
by &pat.;
run;

proc sort data=cv.xv_select out=xv_select_&timepoint.;
by replicate;
run;
proc sort data=cv.xv_noselect out=xv_noselect_&timepoint.;
by replicate;
run;

/* create datasets for each replicate*/
%macro createDatasets(replicate);
   %do i=1 %to &replicate;

	data xv_select_&timepoint._&i; set xv_select_&timepoint.;
	if replicate=&i;
	run;

	data xv_noselect_&timepoint._&i; set xv_noselect_&timepoint.;
	if replicate=&i;
	run;
%end;
%mend createDatasets;

%createDatasets(&sampn);

options mprint mlogic;
/* analyze and output results of the datasets to sas datasets*/
%MACRO Analysis(replicate); 
	%do i=1 %to &replicate;
		ods output BestSubsets= cv.score_&k.anal_&timepoint._&i  ;
			proc phreg data=xv_select_&timepoint._&i;
				model newprogmos*newcens(&&cns)=  &fctrs. / selection=score start=&k. stop=&k. best=1;
			run;
		ods output close;
	%end;
%MEND Analysis;

%Analysis(&sampn);

data cv.score_&k.anal_&timepoint.;
length VariablesInModel $ 300;
set cv.score_&k.anal_&timepoint._1;

run;


%macro lengthofVIM(replicate);
 %do i=2 %to &replicate;
 data cv.score_&k.anal_&timepoint._&i;
 	length VariablesInModel $ 300;
 	set cv.score_&k.anal_&timepoint._&i;
 run;
 %end;
%mend lengthofVIM;

%lengthofVIM(&sampn);

%MACRO Append(replicate); 
	%do i=2 %to &replicate;
	proc append base=cv.score_&k.anal_&timepoint. data = cv.score_&k.anal_&timepoint._&i;
	run;
	%end;
%MEND Append;

%Append(&sampn);



/* Obtaining the frequency of each variable's inclusion with &k variables */
data cv.score_&k.anal_&timepoint.final score_&k.anal_&timepoint._freqs;
set cv.score_&k.anal_&timepoint.;
array var{&k.} $ var1-var&k.;
do count = 1 to &k.;
var[count]=scan(VariablesInModel, count, ' ');
if count=&k. then output cv.score_&k.anal_&timepoint.final;
output score_&k.anal_&timepoint._freqs;
end;
run;

data cv.score_&k.anal_&timepoint.forfreqs;
	set score_&k.anal_&timepoint._freqs;
	a=1;
	array varn{&k.} $ var1-var&k.;
	do until (a = &k.+1 ) ;
	if count=a then var=varn[a] ;
	a=a+1;
	end;
	keep var;
	rename var=Analyte;
run;


proc freq data=cv.score_&k.anal_&timepoint.forfreqs;
	tables Analyte/ nopercent nocum nopercent out=cv.freq;
	run;


data cv.freq;
	set cv.freq;
	by Analyte;
	percent = put(((count / &sampn) *100),6.2);
run;

proc sort data=cv.freq out=cv.freqsummary;
	by analyte;
run;


data selection;
set cv.freqsummary;
retain sum 0;
if percent > &perc1. then ind=1;
else if percent <= &perc1. then ind=0;
sum=sum+ind;
run;

ods rtf file="&dir\&studynm. &timepoint. &ospfs. analytesover&perc1. &sysdate9..rtf";
proc print data=selection;
run;

ods rtf close;

data _null_;
set selection;
call symput('total',strip(sum));
run;
%put &total;

/* selcting model with &total # of  analytes;*/
%MACRO Analysis1(replicate); 
    %if &total. ne &k. %then %do;
	%do i=1 %to &replicate;
		ods output BestSubsets= cv.score_&total.anal_&timepoint._&i ;
			proc phreg data=xv_select_&timepoint._&i;
				model newprogmos*newcens(&&cns)= &&vlist. / selection=score start=&total stop=&total best=1;
			run;
		ods output close;
	%end;

data cv.score_&total.anal_&timepoint.;
length VariablesInModel $ 300;
set cv.score_&total.anal_&timepoint._1;

run;


 %do i=2 %to &replicate;
 data cv.score_&total.anal_&timepoint._&i;
 	length VariablesInModel $ 300;
 	set cv.score_&total.anal_&timepoint._&i;
 run;
 %end;


	%do i=2 %to &replicate;
	proc append base=cv.score_&total.anal_&timepoint. data = cv.score_&total.anal_&timepoint._&i;
	run;
	%end;




/* Obtaining the frequency of each variable's inclusion with &total variables */
data cv.score_&total.anal_&timepoint.final score_&total.anal_&timepoint._freqs;
length var1-var&total. $15;
set cv.score_&total.anal_&timepoint.;
array var{&total.} $ var1-var&total.;
do count = 1 to &total.;
var[count]=scan(VariablesInModel, count, ' ');
if count=&total. then output cv.score_&total.anal_&timepoint.final;
output score_&total.anal_&timepoint._freqs;
end;
run;

data cv.score_&total.anal_&timepoint.forfreqs;
	set score_&total.anal_&timepoint._freqs;
	a=1;
	array varn{&total.} $ var1-var&total.;
	do until (a = &total.+1 ) ;
	if count=a then var=varn[a] ;
	a=a+1;
	end;
	keep var;
	rename var=Analyte;
run;

proc freq data=cv.score_&total.anal_&timepoint.forfreqs;
	tables Analyte/ nopercent nocum nopercent out=cv.highfreq;
	run;

data cv.highfreq;
	set cv.highfreq;
	by Analyte;
	percent = put(((count / &sampn) *100),5.1);
run;

proc sort data=cv.highfreq out=cv.highfreqsummary;
	by analyte;
run;

data selection2;
set cv.highfreqsummary;
retain sum 0;
if percent > &perc2. then ind_2=1;
else if percent <= &perc2. then ind_2=0;
sum=sum+ind_2;
run;

ods rtf file="&dir\final &studynm. &timepoint. &ospfs. analytesover&perc2. &sysdate9..rtf";
proc print data=selection2;
run;

ods rtf close;
%end;

%mend Analysis1;

%Analysis1(&sampn);

%macro keepcolumn(replicate);
   %do i=1 %to &replicate;
	data cv.score_&total.anal_&timepoint._&i;
	set cv.score_&total.anal_&timepoint._&i(keep= VariablesinModel);
	run;
	%end;
%Mend keepcolumn;
%keepcolumn(&sampn);

*multivariate analysis for each replicate using proc phreg;
%Macro Analysis2(replicate);
	%do i=1 %to &replicate;
	data a;
	set cv.score_&total.anal_&timepoint._&i;
	call symput("mpred&total.",VariablesInModel);
	run;

	proc phreg data=xv_select_&timepoint._&i;
	model newprogmos*newcens(&&cns)= &&mpred&total.;
	baseline covariates=&factords._anal out=cv.Preds&i xbeta=myxbeta  /*nomean*/ ;
	run;
	%end;
	
%MEND Analysis2;

%Analysis2(&sampn);


* to obtain the median cutoff for the prediction dataset for each replicate ;         
%Macro Analysis3(replicate);

	%do i=1 %to &replicate;

	Proc means data=cv.preds&i median;
	var myxbeta;
	class &&pat.;
	output out=cv.medpatpred&i median=median;
	run;

	data cv.temp&i;
	set cv.medpatpred&i;
    tmp = _N_ - 1;
	run;
	
	data cv.predtemp;
	set cv.medpatpred&i;
	if _type_=1 then delete;
	run;
	
	data cv.predtemp;
	set cv.predtemp;
	call symputx('overall',median);
	run;
	
	data cv.finalpred&i;
	set cv.temp&i;
    if tmp=&i then do;
	if median <= &overall then group=1;
	else if median > &overall then group=2;
	end;
	run;

	%end;
%mend Analysis3;

%Analysis3(&sampn);

* appending the group data;
%MACRO Appendgrp(replicate); 
	%do i=2 %to &replicate;
	proc append base=cv.finalpred1 data = cv.finalpred&i;
	run;
	%end;
%MEND Appendgrp;
%Appendgrp(&sampn);

* deleting if group missing;
data cv.group2;
	set cv.finalpred1;
    if group=. then delete;
run;

proc freq data=cv.group2;
table group;
run;

*getting pfstime;
data prediction1;
merge cv.group2 (in=a) analysis.&survds. (keep= &pat. &timevr. &censvr.);
by &pat.;
if a;
run;


data prediction2;
set prediction1;
if &timevr.=. then delete;
run;
 

proc sort data=prediction2 out=cv.&studynm._&ospfs._&timepoint._pred;
by group;
run;
ods graphics on;

ods rtf file="&dir\&studynm. cvloo_&ospfs._&timepoint. &sysdate9..rtf";
*lifetest for the analysis;
proc lifetest data=prediction1 plots=(S) graphics censorsymbol=plus;
time &timevr.*&censvr.(&cns);
strata group;
symbol1 v=none color=black;
symbol2 v=none color=blue line=3;
ods exclude prooductlimitestimates;
title "Kaplan Mier survival curves for predicted &ospfs.time for low and high &ospfs.time the multivariate &timepoint. model loocv";
run;
proc phreg data=prediction1;
model &timevr.*&censvr.(&cns)=group/risklimits;
run;
ods rtf close;

data listmodels;
	set cv.Score_&total.anal_&timepoint.final;
	drop control_var count;
run;

ods rtf file="&dir.\&studynm. training_&timepoint. &ospfs. models.rtf";
proc print data=listmodels;
var VariablesInModel NumberInModel ScoreChiSq var1--var&total.;
run;
ods rtf close;

title;
