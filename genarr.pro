;+
; NAME: 	GENARR(N, MIN, MAX, LOG=LOG)
;
; PURPOSE: 	Creates an evenly spaced array of numbers.  Input can
; 		be INTEGER, LONG, FLOAT, or DOUBLE.  Output is the
; 		same type as input.  Note that for INTEGER and LONG
; 		types, this may mean that the output is not spaced
; 		exactly evenly, since there is a minimum increment for
; 		the output array.
;
; INPUTS:	N: Number of points
;		MIN, MAX: The minimum and maximum of the output array.
;
; KEYWORDS: 	LOG: Set this keyword TRUE for logarithmically spaced
; 		numbers.  If this keyword is set and the number type
; 		is INTEGER/LONG, a FLOAT/DOUBLE array is returned.
;
; OUTPUTS:	An evenly spaced array of N numbers from MIN to MAX.
;
; EXAMPLE:	x = GENARR(7, 4., 10.)
;		print, x
;
;		Output:
;		4.00000      5.00000      6.00000      7.00000      8.00000      9.00000
;
; MODIFICATION HISTORY: Created Feb 2008 by Marc Pulupa
;
;-

function genarr, n, min, max, log = log

  n_dim = n_elements(n)

  dim_check = 1.

  if n_elements(min) NE n_dim then dim_check = 0
  if n_elements(max) NE n_dim then dim_check = 0
  if keyword_set(log) then if n_elements(log) NE n_dim then dim_check = 0

  if not(dim_check) then begin
     print, 'All inputs must have the same dimension!'
     return, 0
  end

  if n_dim EQ 1 then begin
     
     arr = -1
     
     smin = size(min[0])
     smax = size(max[0])
     
     type = max([smin[1], smax[1]])

     if n[0] GT 1 then begin
        
        if not(keyword_set(log)) then begin
           case type of
              2:arr = fix((max[0]-min[0])*indgen(n[0])/(n[0]-1)+min[0])
              3:arr = fix((max[0]-min[0])*lindgen(n[0])/(n[0]-1)+min[0])
              4:arr = (max[0]-min[0])*findgen(n[0])/(n[0]-1)+min[0]
              5:arr = (max[0]-min[0])*dindgen(n[0])/(n[0]-1)+min[0]
           end
        endif else begin
           if type eq 2 then type = 4
           if type eq 3 then type = 5
           case type of
              4:arr = min[0]*(max[0]/min[0])^(findgen(n[0])/(n[0]-1)) 
              5:arr = min[0]*(max[0]/min[0])^(dindgen(n[0])/(n[0]-1))
           end
        end

     endif else begin

        arr = [min[0]]

     end

     return, arr
     
  endif else if n_dim EQ 2 then begin

     n0 = n[0]
     min0 = min[0]
     max0 = max[0]

     n1 = n[1]
     min1 = min[1]
     max1 = max[1]

     if keyword_set(log) then begin
        log0 = log[0]
        log1 = log[1]
     endif else begin
        log0 = 0
        log1 = 0
     end
     
     arr0 = genarr(n0, min0, max0, log = log0)

     arr1 = genarr(n1, min1, max1, log = log1)

     arr = [[[reform(rebin(arr0, n0, n1), n0, n1, 1)]], $
            [[reform(transpose(rebin(arr1, n1, n0)), n0, n1, 1)]]]

     return, arr

  end

end
