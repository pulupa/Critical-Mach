function crit_mf_u1x, g, b1, th1, th2, u2x, u2z

  return, (2.*(TAN(th2)*COS(th1)^2.-SIN(th1)*COS(th1))*$
           (u2x/u2z*TAN(th2)-1.)/ $
           (b1*TAN(th1)))^(1./2)
  
end

function crit_mf_u2z, g, b1, th1, th2, u2x

  s1 = sin(th1) & c1 = cos(th1) & t1 = tan(th1)
  s2 = sin(th2) & c2 = cos(th2) & t2 = tan(th2)

  return, (u2x*t2-(u2x^2+1.)*t1/u2x)/ $
          (1+(t1*(c1^2*t2^2-s1^2-b1))/(2*(t2*c1^2-s1*c1)))
  
end

function crit_mf_cf, g, b1, th1

  return, SQRT((2+g*b1)/(2*b1) $
               +(((2+g*b1)/(2*b1))^2 $
                 -2*g/b1*COS(th1)^2)^(1./2))

end

function crit_mf_rh4, th2
  
  COMMON crit_mf_com, g, b1, th1, u2x, teti2
  
  s1 = sin(th1) & c1 = cos(th1) & t1 = tan(th1)
  s2 = sin(th2) & c2 = cos(th2) & t2 = tan(th2)
  
  if u2x le 0 then begin
     u2x_old = u2x
     if u2x EQ -1 then u2x = 3.*c2*sqrt(1./(teti2+1)) else $
        if u2x EQ -2 then u2x = 3.*c2*sqrt(g/(teti2+1))
  endif else u2x_old = u2x

  u2z = crit_mf_u2z(g, b1, th1, th2, u2x)
  
  rh4 = b1*g/(g-1) $
        +(t2*c1^2-s1*c1)*(u2x*t2/u2z-1.)/t1 $
        +2*s1^2 $
        -2*t1*(t2*c1^2-s1*c1)/(u2z*(u2x*t2-u2z))*(g/(g-1.)+u2x^2/2.+u2z^2/2.) $
        -2*s1*c1*t2

  u2x = u2x_old

  return, rh4

end

function crit_mf_thz1
  
  COMMON crit_mf_com, g, b1, th1, u2x, teti2

  s1 = sin(th1) & c1 = cos(th1) & t1 = tan(th1)

  if u2x gt 0 then begin
     thz = atan(tan(th1)*(u2x^2+1)/u2x^2)     
     return, thz
  endif else if u2x EQ -1 then begin
     coeffs = [t1*(1+(teti2+1)/9), $
               -1, $
               (teti2+1)/9*t1]
  endif else if u2x EQ -2 then begin
     coeffs = [t1*(1+(teti2+1)/(9*g)), $
               -1, $
               (teti2+1)/(9*g)*t1]
  end
  
  roots = fz_roots(coeffs)

  real_roots = where(imaginary(roots) LT (Machar()).eps, real_count)

  if real_count EQ 0 then return, 0.

  thz = atan(real_part(roots[real_roots]))

  return, thz

end

function crit_mf_thz2

  COMMON crit_mf_com, g, b1, th1, u2x, teti2

  s1 = sin(th1) & c1 = cos(th1) & t1 = tan(th1)

  if u2x gt 0 then begin
     coeffs = [-2*(u2x^2+1.)/u2x^2*t1, $
               2*(u2x^2+1.)/u2x^2-b1/c1^2-t1^2, $
               0., $
               1.]
  endif else if u2x EQ -1 then begin
     coeffs = [-2.*(1.+(teti2+1.)/9)*t1, $
               2.*(1.+(teti2+1.)/9)-b1/c1^2-t1^2, $
               -2*((teti2+1.)/9*t1), $
               1.+2.*(teti2+1.)/9]
  endif else if u2x EQ -2 then begin
     coeffs = [-2.*(1.+(teti2+1.)/(9*g))*t1, $
               2.*(1.+(teti2+1.)/(9*g))-b1/c1^2-t1^2, $
               -2*((teti2+1.)/(9*g)*t1), $
               1.+2.*(teti2+1.)/(9*g)]
  end

  roots = fz_roots(coeffs)

  real_roots = where(imaginary(roots) LT (Machar()).eps, real_count)

  if real_count EQ 0 then return, 0.

  thz = atan(real_part(roots[real_roots]))

  return, thz

end

function crit_mf, gamma, beta1, theta1, $
                  silent = silent, verbose = verbose, $
                  type = type, trial = trial, $
                  teti2_in = teti2_in, u2x_in = u2x_in, $
                  n_rat = n_rat, b_rat = b_rat

  if not keyword_set(trial) then trial = 0

  COMMON crit_mf_com, g, b1, th1, u2x, teti2

  if trial EQ 2 then begin
     db = 0.05
     dt = 0.5/!RADEG
     mf0 = crit_mf(gamma, beta1+db, theta1, $
                   silent = keyword_set(silent), verbose = keyword_set(verbose), $
                   type = type, trial = 3, $
                   teti2_in = teti2, $
                   n_rat = n_rat0, b_rat = b_rat0)
     mf1 = crit_mf(gamma, beta1-db, theta1, $
                   silent = keyword_set(silent), verbose = keyword_set(verbose), $
                   type = type, trial = 3, $
                   teti2_in = teti2, $
                   n_rat = n_rat1, b_rat = b_rat1)
     mf2 = crit_mf(gamma, beta1, theta1+dt, $
                   silent = keyword_set(silent), verbose = keyword_set(verbose), $
                   type = type, trial = 3, $
                   teti2_in = teti2, $
                   n_rat = n_rat2, b_rat = b_rat2)
     mf3 = crit_mf(gamma, beta1, theta1-dt, $
                   silent = keyword_set(silent), verbose = keyword_set(verbose), $
                   type = type, trial = 3, $
                   teti2_in = teti2, $
                   n_rat = n_rat3, b_rat = b_rat3)
     if not keyword_set(silent) then print, 'Using neighbor method'
     if finite(mf0) and finite(mf1) and finite(mf2) and finite(mf3) then begin
        mf = mean([mf0, mf1, mf2, mf3])
        b_rat = mean([b_rat0, b_rat1, b_rat2, b_rat3])
        n_rat = mean([n_rat0, n_rat1, n_rat2, n_rat3])
        return, mf
     end
  end

  g = gamma
  b1 = beta1
  th1 = theta1
     
  if not keyword_set(u2x_in) then begin

     if not keyword_set(teti2_in) then teti2 = 0. else teti2 = teti2_in
     if not keyword_set(type) then type = '1'

     type = strupcase(strcompress(string(type), /remove_all))

     case type of
        '1':u2x = sqrt(g)
        '2':u2x = sqrt(1/(teti2+1))
        '2G':u2x = sqrt(g)*sqrt(1/(teti2+1))
        'F':u2x = -1
        'FG':u2x = -2
        ELSE:u2x = !VALUES.F_NAN
     end

     case type of
        '1':longtype = 'first'
        '2':longtype = 'second'
        '2G':longtype = '(K85 style) second'
        'F':longtype = 'foreshock'
        'W':longtype = '(low beta) whistler'
        'WN':longtype = '(low beta) nonlinear whistler'
        'WG':longtype = '(low beta) group velocity whistler'
        'WB':longtype = 'whistler'
        ELSE:longtype = ''
     end
     
     if not keyword_set(silent) then begin
        if longtype NE '' then $
           print, 'Solving for ' + longtype + ' critical fast Mach number...'
     end

  endif else begin
     
     teti2 = 0.
     u2x = u2x_in
     if not keyword_set(silent) then $
        print, 'Solving for user input u2x = ' + $
               strcompress(string(u2x), /remove_all) + '...'

  end

  if finite(u2x) then begin

     epsilon= (Machar()).eps    ; Make sure to use SINGLE precision here!
     
     if g LE 1. then begin
        if not keyword_set(silent) then $
           print, 'Gamma must be greater than 1!'
        return, !VALUES.D_NAN
     end
     
     if b1 LE 1e-5 then begin
        if not keyword_set(silent) then $
           print, 'Plasma beta must be greater than 1e-5!'
        b1 = 1e-5
     end
     
     if th1 GE 0.999*!DPI/2 then th1 = 0.999*!DPI/2
     if th1 LE 0.001*!DPI/2 then th1 = 0.001*!DPI/2
     
     th1 = double(th1)
     
     angles_00 = [th1, $
                  crit_mf_thz1(), $
                  crit_mf_thz2(), $
                  !DPI/2]
 
     angles_01 = angles_00[where(angles_00 GE th1 and angles_00 LE !DPI/2)]
     
     angles_02 = angles_01[sort(angles_01)]
     
     angles = angles_02[uniq(angles_02)]
     
     n_intervals = n_elements(angles)-1
     
     if type EQ 'F' or type EQ 'FG' and trial gt 0 then begin

        tmp_angles = angles
        for i = 0, n_intervals-1 do begin
           start = tmp_angles[i]+10*epsilon
           stop = tmp_angles[i+1]-10*epsilon

           n_test = 500
           th2_test = genarr(n_test, start, stop)
           rh4_test = crit_mf_rh4(th2_test)
           sign_change = (rh4_test GE 0) * shift(rh4_test LT 0, 1) + $
                         (rh4_test GE 0) * shift(rh4_test LT 0,-1)
           sign_change[0] = 0
           sign_change[n_test-1] = 0
           sc_ind = where(sign_change EQ 1)
           if sc_ind[0] GE 0 then begin
              angles = [angles, $
                        th2_test[sc_ind-1], $
                        th2_test[sc_ind]]
           end
        end

        angles = angles[sort(angles)]
        angles = angles[uniq(angles)]
        n_intervals = n_elements(angles)-1

     end

     th2_trials = dblarr(n_intervals) + !VALUES.D_NAN

     for i = 0, n_intervals-1 do begin
        
        start = angles[i]+10*epsilon
        stop = angles[i+1]-10*epsilon
        
        th2_trial = zbrent(start, stop, func = 'CRIT_MF_RH4', $
                           max_iter = 1e9, tol = 1.e-10)
        
        if th2_trial NE start then th2_trials[i] = th2_trial
        
     end

     finite_th2 = where(finite(th2_trials), th2_fincount)
     
     if th2_fincount GE 1 then th2s = th2_trials[finite_th2] else th2s = 0.
     
     if not keyword_set(silent) then begin
        th2_test = genarr(1000, th1, !DPI/2)
        
        plot, th2_test*!RADEG, crit_mf_rh4(th2_test), yrange = [-4, 4]

        oplot, [angles*!RADEG], [fltarr(n_elements(angles))], psym = 2
        oplot, [th2s*!RADEG], [fltarr(n_elements(th2s))], psym = 4
     end

     if (type EQ 'F' or type EQ 'FG') and trial EQ 0 and beta1 GT 0.7 then begin
        mf = crit_mf(gamma, beta1, theta1, $
                     silent = keyword_set(silent), verbose = keyword_set(verbose), $
                     type = type, trial = 1, $
                     teti2_in = teti2, $
                     n_rat = n_rat, b_rat = b_rat)
        return, mf
     end

     if th2_fincount EQ 0 and (type EQ 'F' or type EQ 'FG') and $
        trial EQ 1 then begin
        mf = crit_mf(gamma, beta1, theta1, $
                     silent = keyword_set(silent), verbose = keyword_set(verbose), $
                     type = type, trial = 2, $
                     teti2_in = teti2, $
                     n_rat = n_rat, b_rat = b_rat)
        return, mf
     end

     if th2_fincount EQ 0 then begin
        if not keyword_set(silent) then print, 'Cannot find solution!'
        n_rat = !VALUES.F_NAN
        b_rat = !VALUES.F_NAN
        return, !VALUES.F_NAN
     end

     if type EQ 'F' then u2x = 3.*cos(th2s)*sqrt(1./(teti2+1))
     if type EQ 'FG' then u2x = 3.*cos(th2s)*sqrt(g/(teti2+1))

     u2zs = crit_mf_u2z(g, b1, th1, th2s, u2x)
     
     u1xs = crit_mf_u1x(g, b1, th1, th2s, u2x, u2zs)
     
     cf = crit_mf_cf(g, b1, th1)
     
     mfs = u1xs/cf

     if type NE 'F' and type NE 'FG' then begin
        if b1 lt 0.1 then begin 
           
           th2 = max(th2s, sol_ind, /NAN)
           
           mf = mfs[sol_ind]
           
        endif else begin
           
           mf = max(mfs, sol_ind, /NAN)
           
           th2 = th2s[sol_ind]
           
        end
     endif else begin

        th2 = max(th2s, sol_ind, /NAN)

        mf = mfs[sol_ind]

     end

     u2z = u2zs[sol_ind]
     
     u1x = u1xs[sol_ind]
     
     n_rat = (tan(th2)-u2z/sqrt(g))/tan(th1)
     
     b_rat = cos(th1)/cos(th2)

     if not finite(u1x) then begin
        if not keyword_set(silent) then print, 'Cannot find solution!'
        n_rat = !VALUES.F_NAN
        b_rat = !VALUES.F_NAN
        return, !VALUES.F_NAN
     end
     
     if mf LT 1.001 then begin
        
        mf = 1.
        
        n_rat = 1.
        
        b_rat = 1.
        
     end
     
     if n_rat lt 1. and n_rat gt 0.5 then n_rat = 1.
     
     if b_rat lt 1. and b_rat gt 0.5 then b_rat = 1.
     
     if not keyword_set(silent) then begin

        if keyword_set(verbose) then begin

           outf = '(3(A10,F8.3))'
           
           print, ''
           print, 'Input Parameters:'
           
           print, 'gamma:', g, 'beta:', b1, 'theta1:', th1*!RADEG, format = outf
           if type NE '1' then print, 'te2/ti2:', teti2, format = outf

           print, ''
           print, 'Results:'
           
           print, 'theta2:', th2*!RADEG, 'u2z:', u2z, format = outf
           print, 'n2/n1:', n_rat, 'b2/b1:', b_rat, format = outf
           print, 'u1x:', u1x, 'cf:', cf, 'Mf*:', mf, format = outf
           
        endif else begin
           print, 'Mf* = ' + strcompress(string(mf), /remove_all)
        end
     end
     
     return, mf
     
  endif else begin

     mu = 1/42.85                 ; Square root of mass ratio

     csca2 = gamma*beta1/2.

     cfca2 = (1.+csca2)/2. + (((1.+csca2)/2.)^2.-csca2*COS(theta1)^2.)^(1./2)

     case type of
        'W':mf = 1./(2.*mu)*abs(cos(theta1))
        'WG':mf = 1./(mu)*sqrt(27./64)*abs(cos(theta1))
        'WN':mf = 1/sqrt(2.*mu)*abs(cos(theta1))
        'WB':mf = 1./(2.*mu)*abs(cos(theta1))/sqrt(cfca2)
        ELSE:mf = -1
     end

     n_rat = !VALUES.D_NAN
     b_rat = !VALUES.D_NAN

     if mf LT 0 then begin
        if not keyword_set(silent) then begin
           print, 'Mach number type unrecognized.'
           print, 'Available types: "1", "2", "2G", "W", "WG", "WN", "WB", "F", "FG"'
        end
        return, !VALUES.D_NAN
     endif else begin
        if not keyword_set(silent) then $
           print, 'Mf* = ' + strcompress(string(mf), /remove_all)
        return, mf
     end

  end

end
