;+
;
; NAME: f_thick_warm_kappa
;
; PURPOSE: This function returns photon flux from a warm thick-target model
;  Uses simplified form of Equations 24 and 25 from Kontar et al, ApJ 2015
;  http://adsabs.harvard.edu/abs/2015arXiv150503733K
;
; CATEGORY:
;       SPECTRA, XRAYS
;
; CALLING SEQUENCE:
;       Flux = f_thick_warm_kappa(E, params )
;
; CALLS:
;       f_vth, BRM_GAULEG, BRM_BREMCROSS
;
; INPUTS:
; E -energy vector in keV, 2XN edges or N mean energies
; params -model parameters defined below
;   params(0) - Total integrated electron flux, in units of 10^35 electrons sec^-1.
;   params(1) - Power-law index of the electron flux distribution function below
;          eebrk.
;   params(2) - kappa temperature (in keV)
;   params(3) - High energy cutoff in the electron flux distribution function (in keV)
;
;   *********************(above parameters are for non-thermal e-distribution)
;
;   params(4) - Warm plasma density [10^10 cm ^-3]
;   params(5) - Warm plasma temperature [keV]
;   params(6) - Warm plasma length [Mm], 1Mm =10^8 cm
;   params(7)  - relative abundance of elements, default = 1 relative to CHIANTI coronal
;           /basically done as in f_vth/
;
; NOTES:
; the routine is based on f_thick_warm.pro
; and uses elements from f_thin_kappa.pro
;
; MODIFICATION HISTORY:
; eduard(at)astro.gla.ac.uk, April 9, 2018 written
; Kim Tolbert, 25-Apr-2018, Set any NaN and Inf values to small number before returning.
;


function f_thick_warm_kappa, energy_edges, params, _extra=_extra

  ;FTT= f_vth_thick2(e, a[0:5], /no_vth, _extra=_extra)

  print, 'params=', N_elements(params)

  if ((size(energy_edges))[0] eq 2) then edge_products,energy_edges, mean = en else en = energy_edges

  ; constants in CGS units
  mc2 = 510.98d+00 ; keV
  clight = 2.9979d+10 ; speed of light cm/s
  au = 1.496d+13 ;cm
  z = 1.2d0 ; average atomic number
  ; conversion from keV to K
  keV2K = 1.16d7 ;
  KK=2.6d-18 ;

  Ndot = params[0]*1d35
  kappa_index = params[1]
  Tk = params[2]
  ehigh = params[3]

  np =params[4]*1e10
  Tloop=params[5]
  L= params[6]*1e8
  abun=params[7]


  maxfcn = 2048 ; maximum number of points for the Gaussian quadrature integration

  fcoeff = 1.D0/(4.d0*!dpi*(au^2)) ;numerical coefficient for the photon flux
  ; 1/mc2 is due to cross-section

  rerr = 1d-3               ;relative error
  nen = n_elements(en)
  intsum = dblarr(nen)          ; integral sum
  nlim = 2 ; maximum possible order of the Gaussian quadrature integration is 2^nlim.
  conv = 0b ; convergence flag

  for ires = 1, nlim do begin
    npoint = 2L^ires
    ; abscissas and weights
    x = dblarr(nen,npoint)
    w = dblarr(nen,npoint)
    ; brm_gauleg
    brm_gauleg,en,ehigh,npoint,x,w
    lastsum = intsum

    ph_en_array = rebin(en,nen,npoint)
    ; cross section is 2D array according to x,w
    brm_bremcross,x,ph_en_array,z,cross

    nVF_E= Ndot*(x/KK)*($
      (kappa_index-1.)/(1.+x/(kappa_index*Tk))^(kappa_index-2.) $
      -(kappa_index-2.)/(1.+x/(kappa_index*Tk))^(kappa_index-1.))
    ; Equation from Kontar et al paper

    integrand = nVF_E* cross/mc2

    ; integration
    intsum = total(w * integrand,2,/double)
    ; compute relative error
    diff = abs(intsum - lastsum) / intsum
    ind = where(diff gt rerr)
    if (ind[0] eq -1) then begin
      conv = 1b
      break
    endif
  endfor

  if (conv ne 1b) then begin
    print,'f_thick_warm_kappa.pro: Not converged'
  endif

  flux1 =  intsum * fcoeff

  ll=Tloop^2/(2.D0*KK*np) ; collisional stoping distance for electrons of Tloop energy
  Emin=Tloop*3.*(5.d0*ll/L)^4
  EM_add=3.*!PI/2./KK/clight*SQRT(mc2/8.)*Tloop^2/sqrt(Emin)*Ndot
  EM49  =EM_add*1d-49 ; EM in units of 10^49 cm^(-3)
  ;print,'Emin =', Emin, 'Additional EM =', EM49
  ;mean free path of the lowest energy in the injected spectrum
  L_min = Tk^2/(2.*KK*np)/3.
  IF (L_min GT L) THEN message,/info, 'WARNING: min mean free path > loop length: Emin/kT ='+trim(Emin/Tloop)
  ; For low density /or large low energy cut off in the injected spectrum/
  ; thermalization is less important

  IF (Emin GT 0.1 ) THEN BEGIN
    message,/info, 'WARNING: Emin > 0.1 keV,  Setting Emin to 0.1'
    Emin = 0.1
    ; approximation breaks for large Emin
  ENDIF


  Fvth_add=f_vth(energy_edges, [EM49,Tloop,abun])

  flux = flux1 + Fvth_add
  bad = where(~finite(flux), nbad)
  if nbad gt 0 then flux[bad] = 1.0e-27
  return, flux

end
