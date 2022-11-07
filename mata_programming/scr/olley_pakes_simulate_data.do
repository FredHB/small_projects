/*
	author: F Bennhoff
	description: generate a dataset with Olley-Pakes style firm entry and exit.
*/

qui {
cls
clear all
set seed 100000

global root "d:/git_local/small_projects/mata_programming/"
global data "./data/"
cd "$root"

* size of dataset
loc nfirms 3000
loc horizon 11
loc n_obs = `nfirms' * `horizon'

* slope parameters production function
loc bcap 0.425
loc blab 0.6
loc bage -0.025

* depreciation and investment share
loc delta = 0.2
loc invshare = 0.1

* AR(1) mixing parameter for productivity
loc alpha 0.9

* basic shaping of dataset
set obs `n_obs'
seq firm, from(1) to(`nfirms') block(`horizon')
seq year, from(1) to(`horizon') block(1)
gen age = runiformint(0, 20) if year == 1
replace age = age[_n-1] + 1 if age == .

* initial capital draw
gen capital = rnormal(50,2) if year == 1

* generate productivity process
gen productivity = rnormal(0,1) if year == 1
gen prod_innovation = rnormal(0,1)
gen prod_noise = rnormal(0,2)
gen prod_fmean = 0 //runiform(-1,1) // firm fixed effects
replace productivity = prod_fmean + `alpha'*productivity[_n-1] + sqrt(1-`alpha'^2)*prod_innovation if productivity == .
drop prod_innovation prod_fmean

* time fixed effects: boom-bust-boom pattern
gen prod_time_fe = 0 // abs(year - 5)/10

* labor as a function of productivity
gen lab_noise = rnormal(0,2) // unobserved labor demand shifters
loc labfun 1/(1-`blab')*(productivity + prod_time_fe) + `bcap'/(1-`blab')*capital + `bage'/(1-`blab')*age + log(`blab')/(1-`blab') + lab_noise
gen labor = `labfun' if year == 1

* output
loc outfun `bage'*age + `bcap'*capital + `blab'*labor + productivity + prod_time_fe + prod_noise
noi di "loc outfun `bage'*age + `bcap'*capital + `blab'*labor + productivity + prod_time_fe + prod_noise"
gen output = `outfun' if year == 1

* investment in t=1
loc invfun 2 + `invshare'*((1.5 + productivity/4)*capital - .01*age - capital*age/100)
gen investment = `invfun' if year == 1

* fill in investment, capital, output recursively
forval t = 2/`horizon' {
	* capital flow
	replace capital = capital[_n-1]*(1-`delta') + investment[_n-1] if year == `t'
	* labor demand
	replace labor = `labfun' if year == `t'	
	* output
	replace output = `outfun' if year == `t'
	* investment as a function of productivity, capital and age
	replace investment = `invfun' if year == `t'
}

* generate exit threshold which depends on age and capital in the way OP argues
sum capital
loc cap_mean r(mean)
sum age
loc age_mean r(mean)
gen exit_threshold = -( /*5*(age-`age_mean') +*/ `bcap' * (capital-`cap_mean') - (capital-`cap_mean'+5)^2/10)

* standardize exit thresholds, add slight leftshit, to compare it to productivity
sum exit_threshold
replace exit_threshold = (exit_threshold - `r(mean)')/`r(sd)' - 1.8

* exit indicator
gen exit = productivity < exit_threshold if !mi(exit_threshold)

* some summary stats on the dataset
sum exit 
loc meanexit = r(mean)

noi di as text "%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%" ///
_n "about `:di %6.3f `=`r(mean)'*100''% of firms exit each period " ///
_n "for more summary statistics, see below. " ///
_n "%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%"

sort firm year
by firm: egen leaver = total(exit)
by firm: replace leaver = leaver > 0
noi by leaver, sort: sum age capital labor productivity investment output

* save full dataset
sort firm year
noi save $dirdata/full_OP.dta, replace

* keep only variables that were present in problem set dataset
keep firm year capital labor age output investment exit

xtset firm year
ren  (age capital labor output investment) (A K L Y I)
by firm, sort: gen X = sum(exit)
replace X = 1 if X > 1
drop if X == 1 & X[_n-1] == 1 
replace X = 1-X
replace Y = . if exit == 1
replace L = . if exit == 1
replace I = . if exit == 1

* next 3 lines make sure dataset looks like the data in the problem set.
xtset firm year 
replace X = F.X
drop if year == `horizon'
*drop if Y == . | X == .
*drop exit
gen leaves = exit

* save observed dataset
noi save $data/observed_OP.dta, replace

}
