
libname spectra "C:\Users\direction to the file"; run;
data final_spectra; set file name; run;

***define the parameters of interest**;

%LET trait= Heat_stability;run; /*trait to analyse*/
%LET export = Animal_No_; run;  /*unique code for each observation*/
%LET obs=25; run;  /*number of observation to keep in most similar spectra analyses*/

/* remove the missing values*/
data spectra; set final_spectra;
if &trait = . then delete; run;

/*pls with all the observations*/
proc pls data=spectra method=pls nfac=20 cv=one;
 	model &trait = col3--col753;
	output out=outpls p=phat_all;
run; quit;

data outpls_allspectra; set outpls; keep &export phat_all; run;

/*calculate the principal components*/
proc princomp data = spectra std out = total_prins;
var col3--col753;
run;

data total_prins_matrix;
set total_prins;
keep &export prin1--prin4;
run;

/*macro that define for each spectra the nearest observations and then compute PLS*/

/*set the correct number of times that is necessary to repeat the macro in the sentence 
%do n=1 %to X where X are the number of observations in the dataset*/

%macro comp;
%local n;
%do n=1 %to X;  /*to = number of observation in the dataset*/

data reference&n; set total_prins_matrix;
if _n_ =&n then output; run;

data base&n; set reference&n;run;

proc append data = total_prins_matrix base = reference&n;
run;

proc fastclus data = reference&n maxc = 1 replace = none maxiter = 0 noprint
out = maha_to_centroid (rename = (distance = maha_distance_to_centroid) drop = cluster);
var prin:;
run;

proc sort data = maha_to_centroid; by maha_distance_to_centroid; run;

data new;
    set maha_to_centroid(obs=&obs);
run;

proc sql;
create table not_important as
   select * from final_spectra where &export not in (select distinct &export from new);
quit;

proc sql;
create table cluster as
   select * from final_spectra where &export not in (select distinct &export from not_important);
quit;


proc pls data=cluster method=pls nfac=20 cv=one;
 	model &trait = col3--col753;
	output out=outpls&n p=phat;
run; quit;

proc sql;
create table not_important1 as
   select * from outpls&n where &export not in (select distinct &export from base&n);
quit;

proc sql;
create table predicted&n as
   select * from outpls&n where &export not in (select distinct &export from not_important1);
quit;

data predicted&n; set predicted&n;
keep &export &trait phat; run;

%end;
%mend comp;
%comp

/*report to predictX the correct number of prediction done; this value corresponds to the value used in the macro
in the sentence from 1 to X*/
data results; set predicted1-predicted389; run;

proc sort data=outpls_allspectra; by &export; run;

proc sort data = results; by &export; run;

/*final dataset with the results from the PLS with all the observations and with PLS with only the most similar spectra*/
data comparison; merge outpls_allspectra (in=a) results (in=b); if a and b; by &export; run;
