/*
	Coding Session: Introduction to Mata
	Author: Frederik Bennhoff
*/
preserve
clear all
set obs 100

global root "d:/git_local/small_projects/mata_programming"
cd $root

cap program drop makedata
program define makedata
clear
set obs 100
gen M = rnormal(0,1)
gen insamp = M < .5
gen x1 = rnormal(0,1)
gen x2 = 0.5*x1 + rnormal(0,2)
gen y = 0.5*x1 + 2*x2 - 4
end 

/*
	// very basic mata commands
	mata :

	// can hand-input matrices like this
	x = (2, 2 \ 3,2)

	// row and column concatenations are very easy
	((x \ x), (x \ x))

	// for loops have a different syntax
	for (i=1; i<=10; i=1.5*i+2) {
		i
	}

	end
*/

mata:

	/* Regression function */
	void reg_fred(depvar, xvars, touse)
	{
		st_view(Y=., ., tokens(depvar), touse)
		st_view(X=., ., tokens(xvars),  touse)
		
		nobs = rows(X)
		
		"number observations: "  
		nobs
		
		X = (J(nobs, 1, 1), X)
		XX = X'X
		XY = X'Y
		b = lusolve(XX, XY)	
		st_matrix("r(b)", b)
	}

end

cap program drop reg_fred
program define reg_fred, rclass
	version 17.0
	syntax varlist [if] [in]
	marksample touse
	loc depvar : word 1 of `varlist'
	loc xvars : list varlist - depvar

	di "regress `depvar' on `xvars'"
	
	mata reg_fred("`depvar'", "`xvars'", "`touse'")
	mat li r(b)
	return add
end

program define sim_reg, rclass
	clear
	makedata
	reg_fred y x1 x2 if insamp == 1
	matrix b = r(b)
	return scalar b0 = b[1,1]
	return scalar b1 = b[2,1]
	return scalar b2 = b[3,1]
end

sim_reg
ret li

simulate cons = r(b0) coef1 = r(b1) coef2 = r(b2), reps(100) saving("data/regrssion_sim.dta", replace) : sim_reg
use "data/regrssion_sim.dta", clear