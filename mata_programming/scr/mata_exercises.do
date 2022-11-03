/*
	Coding Session: Introduction to Mata
	Author: Frederik Bennhoff
	
	Script contains some basic examples of:
	
	 - optimization in mata
	 - integrating mata routines into stata programs
	
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

makedata


mata: 

/* Optimize a function in Mata */
	void myeval(todo, x, y, g, H) 
	{
			y = -exp(-x^2 + x - 7)
	}
	S = optimize_init()
	optimize_init_evaluator(S, &myeval())
	optimize_init_which(S,"min")	
	optimize_init_params(S, 0)
	b = optimize(S)
	b

/* Optimize a function using the gradient in Mata */
	void myeval_grad(todo, p, y, g, H) 
	{
		y = exp(-p[1]^2 - p[2]^2 - p[1]*p[2] + p[1] - p[2] - 3)
		if (todo==1) {
			g[1] = (-2*p[1] - p[2] + 1) * y
			g[2] = (-2*p[2] - p[1] - 1) * y				
		}
	}
		
	S = optimize_init()
	optimize_init_evaluator(S, &myeval_grad())
	optimize_init_evaluatortype(S, "d1")	
	optimize_init_which(S,"max")	
	optimize_init_params(S, (0,0))
	p = optimize(S)
	p	

/* Regression directly minimizing ssr */
	// objective function
	void ssr(todo, beta, Y, X, ssr, g, H)
	{
		e = (Y - X*beta')
		// squared residuals
		ssr = e'*e
	
		// gradient (wrt beta)
		if (todo==1) {
			jac = 2*(-X' * Y + X' * X * beta')
			g[1] = jac[1]
			g[2] = jac[2]
			g[3] = jac[3]
		}
	}

	// function wrapper for optimization problem
	void reg_minssr(y, x) 
	{

		st_view(Y=., ., tokens(y))
		st_view(X=., ., tokens(x))	
		X = (J(rows(X), 1, 1), X)	
		S = optimize_init()
		optimize_init_evaluator(S, &ssr())
		optimize_init_evaluatortype(S, "d1") // use the gradient in estimation :)
		optimize_init_which(S,"min")
		optimize_init_params(S, (0,0,0))
		optimize_init_argument(S, 1, Y)
		optimize_init_argument(S, 2, X)
		
		//optimize_evaluate(S)

		beta = optimize(S)
		
		// save to stata matrix
		st_matrix("beta", beta)
		
		// print solution
		"Initial sum of squared residuals"
		optimize_result_value0(S)
		"Final sum of squared residuals "
		optimize_result_value(S)
		"Estimate"
		beta
	}
	
	// call program and see that it works
	reg_minssr("y", "x1 x2")
	
end

mata:

	/* Regression solving FOC */
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