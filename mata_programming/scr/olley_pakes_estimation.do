/*
	date: 		 11/07/2022
	author: 	 F Bennhoff
	description: Estimates an Olley & Pakes model of firm entry and exit using
			     different approaches. 
				 
				 Main approach: Mata program.
				 Failing approach: Stata nl estimation
				 Benchmark: Prodest command by Rovigatti et al.
*/

qui {

clear all
set more off

* load data
global root "d:/git_local/small_projects/mata_programming/"
global data "./data/"
cd "$root"
use "$data/observed_OP.dta", clear 
xtset firm year


/*  Variable Definitions  */
loc freevars L
loc statevars K A
loc proxyvars I
loc depvar Y
loc exit leaves

loc polyvars `statevars' `proxyvars'

loc num_polyvars: word count `polyvars'

forval i = 1/`num_polyvars' {		
	
	loc var`i' : word `i' of `polyvars'
	
	forval j = 1/`i' {

		loc var`j' : word `j' of `polyvars'
		gen var`i'_`j' = `var`i'' * `var`j''
		
		loc polyvars `polyvars' var`i'_`j'
		
	}
   
}

loc Lpolyvars
foreach v in `polyvars' {
	gen L`v' = L.`v'
	loc Lpolyvars `Lpolyvars' L`v'
}

loc Lstatevars
foreach v in `statevars' {
	* capture because some vars already exist
	cap gen L`v' = L.`v'
	loc Lstatevars `Lstatevars' L`v'
}


/*  First Stage  */
reg `depvar' `freevars' `polyvars'
gen phihat = 0 if e(sample)
foreach v in `polyvars' {
	replace phihat = phihat + _b[`v']*`v'
}
gen Lphihat = L.phihat

gen `depvar'_free = `depvar' if e(sample)
foreach v in `freevars' {
	scalar b_`v' = _b[`v']
	replace `depvar'_free = `depvar'_free - b_`v' * `v'
}


/*  propensity scores  */
probit `exit' `Lpolyvars'
predict probhat
gen Lprobhat = L.probhat


* ---------------- OP SECOND STAGE WITH SLOOOW STATA ROUTINES ---------------- *
// estimation of the second stage with nl-command is likely to fail because of
// non-convergence. 

/*  Second Stage  */
* unrestricted regression equation

loc prregvars Lphihat `Lstatevars' Lprobhat
loc num_prregvars: word count `prregvars'
forval i = 1/`num_prregvars' {		
	
	loc var`i' : word `i' of `prregvars'
	
	forval j = 1/`i' {

		loc var`j' : word `j' of `prregvars'
		gen ppregv`i'_`j' = `var`i'' * `var`j''
		
		loc prregvars `prregvars' ppregv`i'_`j'
		
	}
   
}
di "`prregvars'"

reg `depvar'_free `statevars'

loc nl_Lstatevars
loc nl_statevars
foreach v in `statevars' {
		loc nl_Lstatevars `nl_Lstatevars' {b_`v'=`=_b[`v']'} * L`v' + 
		loc nl_statevars `nl_statevars'   {b_`v'=`=_b[`v']'} * `v'  + 
}

loc op_command nl (`depvar'_free = {b_0 = `=_b[_cons]'} + `nl_statevars' ///
{b_h}       * ( {b_Lphihat} * Lphihat - {b_0 = `=_b[_cons]'} - (`nl_Lstatevars' 0) ) + ///
{b_h_sq}    * ( {b_Lphihat} * Lphihat - {b_0 = `=_b[_cons]'} - (`nl_Lstatevars' 0) )^2 + ///
{b_pr}      * probhat + ///
{b_pr_sq}   * probhat^2 + ///
{b_pr_h}    * ( {b_Lphihat} * Lphihat - {b_0 = `=_b[_cons]'} - (`nl_Lstatevars' 0) ) * probhat ///
) ///
if !mi(Lphihat) & !mi(probhat), iter(1000) 


noi di "`op_command'"
cap noi `op_command'

}

* ---------------- OP SECOND STAGE WITH FAAAAST MATA ROUTINES ---------------- *
// This works really well and gives us basically the same results as PRODEST.
preserve
clear mata

keep if !mi(Lphihat) &  !mi(phihat)

marksample touse					
sum `touse'
cap gen ONES = 1
			
mata : 
	
	void fun_op(todo, beta, X, LX, PHI, LPHI, PR, PR2, ssr, g, H)
	{
		/* 
		Description of estimation and objective function:
			X = statevars, including a constant
			Y-beta_L'L = phi(t) + eta
			phi(t) + eta = beta_X'X(t) 
							+ poly_2[ phi(t-1) - beta_X'X(t-1) ; theta ] 
							+ poly_2( P(t) ; alpha ) + error
			hence, define Psi(t) = phi(t) - beta_X'X(t).
			given beta_X, choose theta and alpha to minimize:
			SSR = [ Psi(t) - poly_2[ Psi(t-1) ; theta ] - poly_2( P(t) ; alpha ) ]^2
			or using the two-variable polynomial of degree 2:
			SSR = [ Psi(t) - poly_2[ Psi(t-1), P(t) ; theta ] ]^2
			we do the former.
		*/
		
		N = rows(X)
		
		// beta is a column vector
		PSI = PHI - X*beta'
		LPSI = LPHI - LX*beta'
		LPSI2 = LPSI:*LPSI
		
		// PMAT := matrix with RHS variables including the polynomial	
		PMAT = (LPSI, LPSI2, PR, PR2)
		PMATPMAT_inv = invsym(PMAT'*PMAT)
			
		// obtain coefficient on polynomial
		theta = PMATPMAT_inv*PMAT'*PSI
		RES = PSI - PMAT*theta 
		ssr =  - RES'*RES/N
		
		// calculate the derivative ssr wrt beta
		if (todo==1) { // there may be some error remaining in here because
		// gradient different from optimizer numerically calculated gradient...
		// I suspect it is a numerical issue with the gradient calculation...
			g = (10^6,10^6,10^6)
			for (k=1; k<=length(beta); k++) {
				// N x 4, derivative of PMAT wrt parameter #k
				PMAT_wrt_betak =(-LX[1..., k], 
								-2*LX[1..., k] :* (LPHI - LX[1..., k] :* beta[k]), 
								J(rows(X),2,0)
								)
				
				// 4 x 1, derivative of first stage esitmate wrt parameter #k
				theta_wrt_betak = -PMATPMAT_inv * ///
				( PMAT_wrt_betak'*PMAT + PMAT'*PMAT_wrt_betak  ) * ///
				PMATPMAT_inv*(PMAT'*PSI) + ///
				PMATPMAT_inv * (PMAT_wrt_betak' * PSI)
				
				// N x 1, derivative of outcome PSI wrt param #k
				PSI_wrt_betak = -X[1...,k]
					
				// scalar, sum of squares wrt param #k
				ssr_wrt_betak = 2*PSI'*PSI_wrt_betak + ///
				(PSI + 2*PMAT*theta)' * ( PMAT * theta_wrt_betak + PMAT_wrt_betak * theta ) - ///
				sum(2* PMAT * theta :* PSI_wrt_betak)
				
				g[k] = - ssr_wrt_betak/N
					
			}
			
			//"gradient"
			//g
		}
		
		
	}
	
	void optimize_op(phi, Lphi, statevars, Lstatevars, prob)
	{
		// read variables into mata
		st_view(PHI=.,.,tokens(phi), "`touse'")
		st_view(LPHI=.,.,tokens(Lphi), "`touse'")
		st_view(X=.,.,tokens(statevars), "`touse'")
		X = (J(rows(X),1,1), X)
		st_view(LX=.,.,tokens(Lstatevars), "`touse'")
		LX = (J(rows(LX),1,1), LX)		
		st_view(PR=.,.,tokens(prob), "`touse'")
		PR2 = PR :* PR	
		
		// inital value
		beta0 = ( 0, 0.5, 0)
		
		// set stage for optimizer
		S = optimize_init()
		optimize_init_argument(S, 1, X)
		optimize_init_argument(S, 2, LX)
		optimize_init_argument(S, 3, PHI)
		optimize_init_argument(S, 4, LPHI)
		optimize_init_argument(S, 5, PR)
		optimize_init_argument(S, 6, PR2)
		optimize_init_evaluator(S, &fun_op())
		optimize_init_params(S, beta0)
		optimize_init_evaluatortype(S, "d0")	
		optimize_init_technique(S,"nr")
		optimize_init_which(S,"max")	
		optimize_init_trace_params(S, "off") 
		optimize_init_trace_dots(S, "on")
		optimize_init_trace_value(S, "off") 
		optimize_init_trace_gradient(S, "off") 
		
		// run optimizer and make output
		b = optimize(S)	

		printf("\n\nEstimates: [beta_0, beta_K, beta_A]\n")
		optimize_result_params(S) 
		printf("\n\n(Sum of squared residuals)/N:\n")
		optimize_result_value(S)
		printf("\n\nGradient at optimum (numerically):\n")
		optimize_result_gradient(S)
		
		// check analytical gradient		
		beta = optimize_result_params(S)
		N = rows(X)
		PSI = PHI - X*beta'
		LPSI = LPHI - LX*beta'
		LPSI2 = LPSI:*LPSI
		PMAT = (LPSI, LPSI2, PR, PR2)
		PMATPMAT_inv = invsym(PMAT'*PMAT)
		theta = PMATPMAT_inv*PMAT'*PSI
		g = (10^6,10^6,10^6)
		for (k=1; k<=length(beta); k++) {
			// N x 4, derivative of PMAT wrt parameter #k
			PMAT_wrt_betak =(-LX[1..., k], 
							-2*LX[1..., k] :* (LPHI - LX[1..., k] :* beta[k]), 
							J(rows(X),2,0)
							)
			
			// 4 x 1, derivative of first stage esitmate wrt parameter #k
			theta_wrt_betak = -PMATPMAT_inv * ///
			( PMAT_wrt_betak'*PMAT + PMAT'*PMAT_wrt_betak  ) * ///
			PMATPMAT_inv*(PMAT'*PSI) + ///
			PMATPMAT_inv * (PMAT_wrt_betak' * PSI)
			
			// N x 1, derivative of outcome PSI wrt param #k
			PSI_wrt_betak = -X[1...,k]
				
			// scalar, sum of squares wrt param #k
			ssr_wrt_betak = 2*PSI'*PSI_wrt_betak + ///
			(PSI + 2*PMAT*theta)' * ( PMAT * theta_wrt_betak + PMAT_wrt_betak * theta ) - ///
			sum(2* PMAT * theta :* PSI_wrt_betak)
			
			g[k] = - ssr_wrt_betak/N
				
		}
		
		printf("\n\nAnalytical gradient\n")
		g
		printf("\n\nAnalytical gradient length\n")
		(g[1]^2 + g[2]^2 + g[3]^2)^0.5
		
	}
	
	optimize_op("phihat", "Lphihat", "K A", "LK LA", "probhat")

end


* ---------------------------------- PRODEST --------------------------------- *
// unsurprisingly, this works.
* compare to PRODEST command (BFGS). (NB: Nelder-Mead yields rubbish here...)
prodest Y, state(A K) proxy(I) free(L) met(op) att poly(2) reps(3) tol(0.000001) opt(bfgs)
