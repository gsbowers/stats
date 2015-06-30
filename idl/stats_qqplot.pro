;+
; Procedure: stats_qqplot
;
; Purpose:
;   creates a normal probability plot for input, data
;
; Inputs: 
;   data: assumed gaussian distribution
; 
; Output: 
;
; Example:
;   stats_qqplot, randomn(0.01, 10000)
;
; Reference: 
;   http://www.itl.nist.gov/div898/handbook/eda/section3/normprpl.htm
;
; Author:
;   Gregory Bowers
;   gsbowers@ucsc.edu
;
; Notes:
;   written for IDL Version 8.0
;+

pro stats_qqplot, data

	n = n_elements(data)

	;define uniform order statistic medians
	U = dblarr(n)

	U[0] = 1-0.5^(1.0d/n)	
	U[-1] = 0.5^(1.0d/n)
	for i=1,n-2 do $ 
		U[i] = (i+1-0.3175)/(n+0.365d)

	;define normal order statistic medians
	GU = dblarr(n)
	for i=0,n-1 do $
	 GU[i] = -GAUSS_CVF(U[i])

	;create ordered response values 
	Z = (data-mean(data))/stddev(data)
	Z = Z[sort(Z)]

	;create plot of ordered response values vs ordered response values
	plot, GU, Z, xtitle='Normal N(0,1) Order Statistic Medians', $
		ytitle='Ordered Response', title='Normal Probability Plot', $
		psym=1, charsize=1.8
	oplot, GU, GU
	xyouts, 0.3,0.8, string(format='(%"corr(x,y)=%6.4g")',correlate(Z, GU)), /normal, charsize=1.5 

	print, "compare corr(x,y) to critical values of Normal PCC Distribution:"
	print, "http://www.itl.nist.gov/div898/handbook/eda/section3/eda3676.htm"

end
