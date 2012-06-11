pro crit_mf_plot, gamma = gamma, beta = beta, $
                  ngamma = ngamma, g_max = gamma_max, g_min = gamma_min, $
                  nbeta = nbeta, b_max = beta_max, b_min = beta_min, $
                  ntheta = ntheta, t_max = theta_max, t_min = theta_min, $
                  nxtheta = nxtheta, tx_max = thetax_max, tx_min = thetax_min, $
                  cpoints = cpoints, surface = surface, $
                  n_rat = n_rat, b_rat = b_rat, type = type, $
                  teti2_in = teti2, u2x_in = u2x_in
  
  if not keyword_set(gamma) then gamma = 5./3
  if not keyword_set(beta) then beta = 0.

  if not keyword_set(nbeta) then nbeta = 50
  if not keyword_set(beta_max) then beta_max = 4.
  if not keyword_set(beta_min) then beta_min = 1.e-5
  
  if not keyword_set(ntheta) then ntheta = 50
  if not keyword_set(theta_max) then theta_max = !DPI/2
  if not keyword_set(theta_min) then theta_min = 0
  
  if not keyword_set(ngamma) then ngamma = 1
  if not keyword_set(gamma_max) then gamma_max = 5./3
  if not keyword_set(gamma_min) then gamma_min = 4./3

  if not keyword_set(u2x_in) then u2x_in = 0.0
  if not keyword_set(teti2) then teti2 = 0.0

  if not keyword_set(type) then type = '1' else begin
     type = strcompress(string(type), /remove_all)
  end

  if ngamma EQ 1 then begin
     
     x=genarr([nbeta, ntheta], [beta_min, theta_min], [beta_max, theta_max])
     
     if keyword_set(nxtheta) then begin
        xx=genarr([nbeta, nxtheta], [beta_min, thetax_min], [beta_max, thetax_max])
        ntheta = ntheta + nxtheta
        x = [[x], [xx]]
        all_thetas = x[0, *, 1]
        sort_thetas = sort(all_thetas)
        for i = 0, nbeta-1 do begin
           x[i, *, 1] = (x[i, *, 1])[sort_thetas]
           x[i, *, 0] = (x[i, *, 0])[sort_thetas]
        end
     end

     betas = x[*, *, 0]
     thetas = x[*, *, 1]
     mfs = fltarr(nbeta, ntheta)
     nrats = fltarr(nbeta, ntheta)
     brats = fltarr(nbeta, ntheta)

;     stop

     for i = 0l, n_elements(betas)-1 do begin
        mfs[i] = crit_mf(gamma, betas[i], thetas[i], /silent, $
                         n_rat = nrati, b_rat = brati, $
                         type = type, teti2 = teti2, u2x_in = u2x_in)
        nrats[i] = nrati
        brats[i] = brati
     end
     
     if keyword_set(cpoints) then nodata = 0 else nodata = 1
     
     if !D.NAME EQ 'PS' then begin
        xthick = 3 & ythick = 3 & thick = 3 & charthick = 3 & charsize = 2
     endif else begin
        xthick = 0 & ythick = 0 & thick = 0 & charthick = 0 & charsize = 0
     end

     case type of
        '1':title = 'M_{f1}^{*}'
        '2':title = 'M_{f2}^{*}'
        '2G':title = 'M_{f2}^{*}(K85)'
        'W':title = 'M_{w} (low \beta)'
        'WN':title = 'M_{nw} (low \beta)'
        'WG':title = 'M_{g} (low \beta)'
        'WB':title = 'M_{w}'
        'F':title = 'M_{fs}'
        'FG':title = 'M_{fs}(K85)'
        ELSE:title = ''
     endcase 

     title = title + ' ('
     title = title + '\gamma = ' + $
             string(gamma, format = '(F5.3)')
     if teti2 NE 0.0 or type EQ '2' or type EQ '2OLD' then $
        title = title + ', T_{e2}/T_{i2} = ' + $
                string(teti2, format = '(F4.2)')
     title = title + ')'

     if not (keyword_set(surface) or $
        keyword_set(n_rat) or $
        keyword_set(b_rat)) then begin
        
        plot, thetas*!RADEG, betas, psym = 3, nodata = nodata, $
              xrange = minmax(thetas)*!RADEG, /xsty, $
              yrange = minmax(betas), /ysty, $
              xtitle = textoidl('\theta_{bn1}'), $
              ytitle = textoidl('\beta_{1}'), $
              title = textoidl(title), $
              xthick = xthick, ythick = ythick, $
              charsize = charsize/2.5, $
              charthick = charthick

        if max(mfs) lt 10 then begin
           clevels = [1.001, findgen(100)/10+1.1]
           clabels = [1, reform(rebin([1, 0], 2, 50), 100)]
           cannots = strarr(n_elements(clevels))
           cannots[where(clabels EQ 1)] = string(clevels[where(clabels EQ 1)], $
                                                 format='(F3.1)')
        endif else if type EQ 'F' or type EQ 'FG' then begin
           clevels = [1.01, 2., 3., 4., 5., 10.]
           clabels = [1, 1, 1, 1, 1, 1]
           cannots = strarr(n_elements(clevels))
           cannots[where(clabels EQ 1)] = string(clevels[where(clabels EQ 1)], $
                                                 format='(I2)')
        endif else begin
           clevels = [findgen(31)+1]
           clabels = [1, reform(rebin([1, 0], 2, 15), 30)]
           cannots = strarr(n_elements(clevels))
           cannots[where(clabels EQ 1)] = string(clevels[where(clabels EQ 1)], $
                                                 format='(I2)')
        end

        if type EQ 'WB' then begin
           clevels = clevels[1:*]
           clabels = clabels[1:*]
           cannots = cannots[1:*]
        end

        contour, mfs, thetas*!RADEG, betas, /over, $
                 levels = clevels, $
                 c_labels = [1, reform(rebin([1, 0], 2, 25), 50)], $
                 c_annot = cannots, c_thick = thick, $
                 c_charsize = charsize/3., c_charthick = charthick
        
     endif else if not (keyword_set(n_rat) or keyword_set(b_rat)) then begin
        
        surface, mfs, thetas*!RADEG, betas, $
                 max_value = max(mfs), min_value = 1.0, $
                 az = 210., ax = 25., xrange = [0., 90.], $
                 yrange = [0., max(betas)], $
                 zrange = [1.0, max(mfs)], $
                 xstyle = 5, ystyle = 5, zstyle = 5, $
                 position = [0.2, 0.25, 1.0, 1.0, 0., 1.], /save
        
        down = 0.05*max(mfs)
        
        axis, 90., 4., 1.0, xaxis = 1, xrange = [0., 90.], $
              /xsty, /ysty, /zsty, $
              xticks = 2, xtickv = [0., 45., 90.], $
              xtickname = [' ', ' ', ' '], charsize = charsize/3, $
              /t3d
        
        xyouts, 0., max(betas)*1.05, z = 1-down, '0!Uo!N', align = 0.5, /t3d
        xyouts, 45., max(betas)*1.05, z = 1-down, '45!Uo!N', align = 0.5, /t3d
        xyouts, 90., max(betas)*1.05, z = 1-down, '90!Uo!N', align = 0.5, /t3d
        xyouts, 45., max(betas)*1.1, z = 1-2*down, textoidl('\theta_{bn}'), /t3d
        
        axis, 0., 0., 1.0, xaxis = 1, xrange = [0., 90.], $
              /xsty, /ysty, /zsty, $
              xticks = 2, xtickv = [0., 45., 90.], $
              xtickname = [' ', ' ', ' '], charsize = charsize/3, $
              /t3d
        
        axis, 90., 0., 1.0, yaxis = 1, yrange = [0., max(betas)], $
              /xsty, /ysty, /zsty, $
              yticks = 2, ytickv = [0., max(betas)/2., max(betas)], $
              ytickname = [' ', ' ', ' '], charsize = charsize/3, /t3d
        
        xyouts, 95., 0., z = 1-down, string('0.', format = '(F3.1)'), align = 0.5, /t3d
        xyouts, 95., max(betas)/2, z = 1-down, string(max(betas)/2, format = '(F3.1)'), align = 0.5, /t3d
        xyouts, 95., max(betas), z = 1-down, string(max(betas), format = '(F3.1)'), align = 0.5, /t3d
        xyouts, 100., max(betas)/2, z = 1-2*down, textoidl('\beta'), /t3d
        
        axis, 00., 0., 1.0, yaxis = 1, yrange = [0., max(betas)], $
              /xsty, /ysty, /zsty, $
              yticks = 2, ytickv = [0., max(betas)/2., max(betas)], $
              ytickname = [' ', ' ', ' '], charsize = charsize/3, /t3d
        
        axis, 90., 0., 1.0, zaxis = 1, zrange = [1., max(mfs)], $
              /xsty, /ysty, /zsty, zticks = 1, ztickv = [1.0, max(mfs)], /t3d, $
              charsize = charsize/3, ztickname = [' ', ' ']
        
        xyouts, 100., 0., z = 1.0, '1.0', /t3d
        xyouts, 100., 0., z = max(mfs), string(max(mfs), format = '(F4.2)'), /t3d
        
        axis, 00., 0., 1.0, zaxis = 1, zrange = [1., max(mfs)], $
              /xsty, /ysty, /zsty, zticks = 1, ztickv = [1.0, max(mfs[*, 0])], /t3d, $
              charsize = charsize/3, ztickname = [' ', ' ']
        
        xyouts, -3., 0., z = max(mfs[*, 0]), string(max(mfs[*, 0]), $
                                                    format = '(F4.2)'), /t3d
        
        
        title = textoidl('M_{f1}^{*} (\gamma = ') + $
                string(gamma, format = '(F5.3)') + ')'
        
        xyouts, 0.5, 0.9, /normal, title, align = 0.5

     endif else if not keyword_set(n_rat) then begin

        plot, thetas*!RADEG, betas, psym = 3, nodata = nodata, $
              xrange = minmax(thetas)*!RADEG, /xsty, $
              yrange = minmax(betas), /ysty, $
              xtitle = textoidl('\theta_{bn}'), $
              ytitle = textoidl('\beta'), $
              title = textoidl('B_{2}/B_{1} (\gamma = ') + $
              string(gamma, format = '(F5.3)') + ')', $
              xthick = xthick, ythick = ythick, charthick = charthick
        
        clevels = [1.001, findgen(50)/10+1.1]
        clabels = [1, reform(rebin([1, 0], 2, 25), 50)]
        cannots = strarr(n_elements(clevels))
        cannots[where(clabels EQ 1)] = string(clevels[where(clabels EQ 1)], $
                                              format='(F3.1)')
        
        contour, brats, thetas*!RADEG, betas, /over, $
                 levels = clevels, $
                 c_labels = [1, reform(rebin([1, 0], 2, 25), 50)], $
                 c_annot = cannots, c_thick = thick, c_charthick = charthick

     endif else begin

        plot, thetas*!RADEG, betas, psym = 3, nodata = nodata, $
              xrange = minmax(thetas)*!RADEG, /xsty, $
              yrange = minmax(betas), /ysty, $
              xtitle = textoidl('\theta_{bn}'), $
              ytitle = textoidl('\beta'), $
              title = textoidl('N_{2}/N_{1} (\gamma = ') + $
              string(gamma, format = '(F5.3)') + ')', $
              xthick = xthick, ythick = ythick, charthick = charthick
        
        clevels = [1.001, findgen(50)/10+1.1]
        clabels = [1, reform(rebin([1, 0], 2, 25), 50)]
        cannots = strarr(n_elements(clevels))
        cannots[where(clabels EQ 1)] = string(clevels[where(clabels EQ 1)], $
                                              format='(F3.1)')
        
        contour, nrats, thetas*!RADEG, betas, /over, $
                 levels = clevels, $
                 c_labels = [1, reform(rebin([1, 0], 2, 25), 50)], $
                 c_annot = cannots, c_thick = thick, c_charthick = charthick

     endelse

  endif else begin
     
     x=genarr([ngamma, ntheta], [gamma_min, theta_min], [gamma_max, theta_max])
     
     gammas = x[*, *, 0]
     thetas = x[*, *, 1]
     mfs = fltarr(ngamma, ntheta)
     
     for i = 0l, n_elements(gammas)-1 do begin
        mfs[i] = crit_mf(gammas[i], beta, thetas[i], /silent)
     end
     
     if keyword_set(cpoints) then nodata = 0 else nodata = 1
     
     if !D.NAME EQ 'PS' then begin
        xthick = 3
        ythick = 3
        thick = 3
        charthick = 3
        charsize = 3
     endif else begin
        xthick = 0
        ythick = 0
        thick = 0
        charthick = 0
        charsize = 0
     end
     
     if not keyword_set(surface) then begin
        
        plot, thetas*!RADEG, gammas, psym = 3, nodata = nodata, $
              xrange = minmax(thetas)*!RADEG, /xsty, $
              yrange = minmax(gammas), /ysty, $
              xtitle = textoidl('\theta_{bn}'), $
              ytitle = textoidl('\gamma'), $
              title = textoidl('M_{f1}^{*} (\beta = ') + $
              string(beta, format = '(F5.3)') + ')', $
              xthick = xthick, ythick = ythick, charthick = charthick
        
        clevels = [1.001, findgen(50)/10+1.1]
        clabels = [1, reform(rebin([1, 0], 2, 25), 50)]
        cannots = strarr(n_elements(clevels))
        cannots[where(clabels EQ 1)] = string(clevels[where(clabels EQ 1)], $
                                              format='(F3.1)')
        
        contour, mfs, thetas*!RADEG, gammas, /over, $
                 levels = clevels, $
                 c_labels = [1, reform(rebin([1, 0], 2, 25), 50)], $
                 c_annot = cannots, c_thick = thick, c_charthick = charthick
        
     endif else begin

        surface, mfs, thetas*!RADEG, gammas, $
                 max_value = max(mfs), min_value = 1.0, $
                 az = 210., ax = 25., xrange = [0., 90.], $
                 yrange = [1., max(gammas)], $
                 zrange = [1.0, max(mfs)], $
                 xstyle = 5, ystyle = 5, zstyle = 5, $
                 position = [0.2, 0.25, 1.0, 1.0, 0., 1.], /save
        
        down = 0.05*max(mfs)
        
        axis, 90., max(gammas), 1.0, xaxis = 1, xrange = [0., 90.], $
              /xsty, /ysty, /zsty, $
              xticks = 2, xtickv = [0., 45., 90.], $
              xtickname = [' ', ' ', ' '], charsize = charsize/3, $
              /t3d
        
        xyouts, 0., max(gammas)*1.05, z = 1-down, '0!Uo!N', align = 0.5, /t3d
        xyouts, 45., max(gammas)*1.05, z = 1-down, '45!Uo!N', align = 0.5, /t3d
        xyouts, 90., max(gammas)*1.05, z = 1-down, '90!Uo!N', align = 0.5, /t3d
        xyouts, 45., max(gammas)*1.1, z = 1-2*down, textoidl('\theta_{bn}'), /t3d
        
        axis, 0., 1., 1.0, xaxis = 1, xrange = [0., 90.], $
              /xsty, /ysty, /zsty, $
              xticks = 2, xtickv = [0., 45., 90.], $
              xtickname = [' ', ' ', ' '], charsize = charsize/3, $
              /t3d
        
        axis, 90., 1., 1.0, yaxis = 1, yrange = [1., max(gammas)], $
              /xsty, /ysty, /zsty, $
              yticks = 2, ytickv = [1.0, (max(gammas)+1.)/2, max(gammas)], $
              ytickname = [' ', ' ', ' '], charsize = charsize/3, /t3d
        
        xyouts, 95., 1., z = 1-down, string(1., format = '(F3.1)'), align = 0.5, /t3d
        xyouts, 95., (max(gammas)+1.)/2, z = 1-down, string((max(gammas)+1.)/2, format = '(F3.1)'), align = 0.5, /t3d
        xyouts, 95., max(gammas), z = 1-down, string(max(gammas), format = '(F3.1)'), align = 0.5, /t3d
        xyouts, 100., (max(gammas)+1.)/2, z = 1-2*down, textoidl('\gamma'), /t3d

        axis, 00., 1., 1.0, yaxis = 1, yrange = [1., max(gammas)], $
              /xsty, /ysty, /zsty, $
              yticks = 2, ytickv = [1.0, (max(gammas)+1.)/2, max(gammas)], $
              ytickname = [' ', ' ', ' '], charsize = charsize/3, /t3d
        
        axis, 90., 1., 1.0, zaxis = 1, zrange = [1., max(mfs)], $
              /xsty, /ysty, /zsty, zticks = 1, ztickv = [1.0, max(mfs)], /t3d, $
              charsize = charsize/3, ztickname = [' ', ' ']
        
        xyouts, 100., 1., z = 1.0, '1.0', /t3d
        xyouts, 100., 1., z = max(mfs), string(max(mfs), format = '(F5.2)'), /t3d
        
        axis, 00., 1., 1.0, zaxis = 1, zrange = [1., max(mfs)], $
              /xsty, /ysty, /zsty, zticks = 1, ztickv = [1.0, max(mfs[*, 0])], /t3d, $
              charsize = charsize/3, ztickname = [' ', ' ']
        
        xyouts, -3., 1., z = max(mfs[*, 0]), string(max(mfs[*, 0]), $
                                                    format = '(F5.2)'), /t3d
        
        
        title = textoidl('M_{f1}^{*} (\beta = ') + $
                string(beta, format = '(F5.3)') + ')'
        
        xyouts, 0.5, 0.9, /normal, title, align = 0.5
        

     end

  end

end
