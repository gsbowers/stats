;+
; Procedure: stats_qqplot
;
; Purpose:
;   finds outliers for given data set using generalized (extreme 
;   studentized deviate) ESD test with hampel identifier 
;
; Inputs: 
;   data: data consisting of gaussian background and outliers
;   k: max number of outliers to look for
; 
; Output: 
;   w: indices into outliers  
;
; Keywords:
;   alpha:  confidence level (default = 0.05)
;
;
; Example:
;    IDL> a=[-0.25,0.68,0.94,1.15,1.20,1.26,1.26,$
;    IDL> 1.34,1.38,1.43,1.49,1.49,1.55,1.56,$
;    IDL> 1.58,1.65,1.69,1.70,1.76,1.77,1.81,$
;    IDL> 1.91,1.94,1.96,1.99,2.06,2.09,2.10,$
;    IDL> 2.14,2.15,2.23,2.24,2.26,2.35,2.37,$
;    IDL> 2.40,2.47,2.54,2.62,2.64,2.90,2.92,$
;    IDL> 2.92,2.93,3.21,3.26,3.30,3.59,3.68,$
;    IDL> 4.30,4.64,5.34,5.42,6.01]
;    IDL> w = stats_outliers(a, 10)           
;    IDL> print, w
;              53          52          51
;    IDL> print, a[w]
;          6.01000      5.42000      5.34000
;    IDL> b = [a[26:-1], 100,1e9,1e256, a[0:25]]
;    IDL> w = stats_outliers(b, 10)                              
;    IDL> print, w   
;          30          29          28          27          26          25
;    IDL> print, b[w]
;      Inf  1.00000e+09      100.000      6.01000      5.42000      5.34000
;
; Reference: 
;   http://www.itl.nist.gov/div898/handbook/eda/section3/eda35h3.htm
;
; Author:
;   Gregory Bowers
;   gsbowers@ucsc.edu
;
; Notes:
;   written for IDL Version 8.0
;+

function stats_outliers, data, k, alpha=alpha

	;significance level
	if ~keyword_set(alpha) then alpha = 0.05

	n = double(n_elements(data))	

	;use hampel identifier 
	mad = 1.4826 * median(abs(data-median(data)))

	;compute test statistics		
	r = abs(data-median(data))/mad
	rmax = reverse((r(sort(r)))[-k:-1]) ;get k values that maximize r
	datai = reverse((sort(r))[-k:-1]) ;associated data indices

	;compute critical values
	lambda = dblarr(k)
	for i=0,k-1 do begin 
		p = 1.0d - alpha/(2.0d*(n-i+1.0d))
		t = - T_CVF(double(p), double(n-i-1.0d))
		lambda[i] = (n-i)*t/sqrt((n - i - 1.0d + t^2.0d)*(n - i + 1.0d))
	endfor

	;compare test statistics to critical values
	test = rmax - lambda
	w = where(test gt 0, count) ;find outliers

	if count ne 0 then begin
		return, datai[w] 
	endif else begin
		return, -1
	endelse

end
