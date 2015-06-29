function zerofunc, z

	common share, prob
	return, prob - gauss_pdf(z)

end

;+
; Function: stats_inv_normalcdf
;
; Purpose:
;   Calculate the inverse of the percent point function (gauss_pdf) 
;   This function returns the standard normal deviate, Z
;   corresponding to the probability p, returned by the cumulative
;   normal distribution function p = CDF(Z)
;
;   This functions determines Z by finding the root to the equation
; 
;                          p - CDF(Z) = 0
;
;   using IDLs fx_root method over the interval Z = [-5,0,5]
;
; Inputs: 
;   P:  Cumulative gaussian probabilty.  Scaler or Array
; 
; Output: 
;   Z:  Standard normal deviate associated with probability P 
;   that all deviates from distribution are less than or equal to Z 
;
; Example: 
;   IDL> print, inv_normalcdf(gauss_pdf(1.0))
;       1.0000001
;   IDL> print, inv_normalcdf(gauss_pdf(2.0))
;       1.9999999
;   IDL> print, inv_normalcdf(gauss_pdf(4.5))
;       4.5000130
;   IDL> print, inv_normalcdf(gauss_pdf(4.8))
;       6.9999892


function stats_inv_normalcdf, p

	COMMON share, prob
	x = [-6.0d,0.0d,6.0d]

	if isa(p,/array) then begin
		root = dblarr(n_elements(p))
		for i=ULONG64(0),n_elements(p)-1 do begin
			prob = p[i]
			root[i] = FX_ROOT(x, 'zerofunc', /double, tol=1e-5, /stop)
		endfor
	endif else begin
		prob = p
		root = FX_ROOT(x, 'zerofunc', /double, tol=1e-5, /stop)
	endelse

	return, root

end
