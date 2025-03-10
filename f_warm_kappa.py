pre_path1='C:\\Users\\xixif\\Dropbox\\' #in working laptop
pre_path2='E:\\Dropbox\\' # in own laptop

import sys

if sys.prefix == 'C:\\Users\\xixif\\anaconda3\\envs\\envte':
    pre_path=pre_path2
if sys.prefix == 'C:\\Users\\xixif\\anaconda3\\envs\\wkenv':
#if sys.prefix == 'C:\\Users\\xixif\\anaconda3\\envs\\envte':
    pre_path=pre_path2
elif sys.prefix == 'D:\\software\\envs\\envte':
    pre_path=pre_path1
import sys
sys.path.append(pre_path+'Kappa study\\Following up\\Programming\\sim_package\\sunkit-spex')

def kap_get_integrand_log(x_log, *, Ndot, kappa_index, Tk, photon_energy, z):
    """
    Return the value of the integrand for the thick-target bremsstrahlung models with kappa distribution
    Wrong now!!! need to correct

    Parameters
    ----------
    x_log : `numpy.array`
        Log of the electron energies
    Ndot, kappa_index, Tk : kappa distribution
        Electron distribution as function of energy
    photon_energy : `numpy.array`
        Photon energies
    z : `float`
        Mean atomic number of plasma

    Returns
    -------
    `numpy.array`
        The values of the integrand at the given electron_energies

    References
    ----------see f_thick_warm_kappa.pro

    """
    import sunkit_spex
    from sunkit_spex.legacy.emission import collisional_loss,bremsstrahlung_cross_section
    import numpy as np
    mc2 = 510.98 # keV
    clight = 2.9979e10 # speed of light cm/s
    au = 1.496e13 #cm
    keV2K = 1.16e7 #
    KK=2.6e-18 # cm^2 KeV^2, take Coulomg logarithm to be 20
    r0 = 2.8179e-13 #classical electron radius, in cm

    # L=log10 (E), E=l0L and dE=10L ln(10) dL hence the electron_energy * np.log(10) below
    ele_band = 10**x_log

    pho_band = photon_energy
    dedt=collisional_loss(ele_band)
    decoeff = 4*np.pi*r0**2*clight*mc2 #coefficient apply to the dedt
    dedt = dedt*decoeff
    v_c=np.sqrt(2 * ele_band / mc2 + ele_band**2 / mc2**2)/(1 + ele_band/mc2) * clight #electron speed
    dedx = dedt / v_c # dE/dx=(dE/dt)v
    ##bremsstrahlung cross section
    cross=bremsstrahlung_cross_section(ele_band, pho_band, z=z)
    nVF_E= Ndot*(1/dedx)*((kappa_index-1)/(1+ele_band/(kappa_index*Tk))**(kappa_index-2)\
    -(kappa_index-2)/(1.+ele_band/(kappa_index*Tk))**(kappa_index-1.))
    integrand = nVF_E* cross/mc2
    return (integrand)

def kap_get_integrand(xt, *, Ndot, kappa_index, Tk, photon_energy, z):
    """
    Return the value of the integrand for the thick-target bremsstrahlung models with kappa distribution

    Parameters
    ----------
    x_log : `numpy.array`
        Log of the electron energies
    Ndot, kappa_index, Tk : kappa distribution
        Electron distribution as function of energy
    photon_energy : `numpy.array`
        Photon energies
    z : `float`
        Mean atomic number of plasma

    Returns
    -------
    `numpy.array`
        The values of the integrand at the given electron_energies

    References
    ----------see f_thick_warm_kappa.pro

    """
    import sunkit_spex
    from sunkit_spex.legacy.emission import collisional_loss,bremsstrahlung_cross_section
    import numpy as np
    mc2 = 510.98 # keV
    clight = 2.9979e10 # speed of light cm/s
    au = 1.496e13 #cm
    keV2K = 1.16e7 #
    KK=2.6e-18 # cm^2 KeV^2, take Coulomg logarithm to be 20
    r0 = 2.8179e-13 #classical electron radius, in cm

    # L=log10 (E), E=l0L and dE=10L ln(10) dL hence the electron_energy * np.log(10) below
    ele_band = xt

    pho_band = photon_energy
    dedt=collisional_loss(ele_band)
    decoeff = 4*np.pi*r0**2*clight*mc2 #coefficient apply to the dedt
    dedt = dedt*decoeff
    v_c=np.sqrt(2 * ele_band / mc2 + ele_band**2 / mc2**2)/(1 + ele_band/mc2) * clight #electron speed
    dedx = dedt / v_c # dE/dx=(dE/dt)v
    ##bremsstrahlung cross section
    cross=bremsstrahlung_cross_section(ele_band, pho_band, z=z)
    nVF_E= Ndot*(1/dedx)*((kappa_index-1)/(1+ele_band/(kappa_index*Tk))**(kappa_index-2)\
    -(kappa_index-2)/(1.+ele_band/(kappa_index*Tk))**(kappa_index-1.))
    integrand = nVF_E* cross/mc2
    return (integrand)

def f_warm_kappa(total_eflux, kappa_index, kappa_temperature , plasma_d, loop_temp, length, high_cut = 1e4, energies=None):
    '''this is the python version of the warm-target model in kappa distribution
    input:
    energies: in kev
    parameters:
;   total_eflux - Total integrated electron flux, in units of 10^35 electrons sec^-1.
;   kappa_index - Kappa index
;   kappa_temperature - kappa temperature (in keV)
;   high_cut - High energy cutoff in the electron flux distribution function (in keV),default to be 1e4 keV
;
;   *********************(above parameters are for non-thermal e-distribution)
;
;   plasma_d - Warm plasma density [10^10 cm ^-3]
;   loop_temp - Warm plasma temperature [keV]
;   length - Warm plasma length [Mm], 1Mm =10^8 cm
    '''
    import sunkit_spex
    from sunkit_spex import fitting_legacy
    from sunkit_spex.legacy.integrate import gauss_legendre
    from sunkit_spex.legacy.emission import collisional_loss,bremsstrahlung_cross_section
    import numpy as np
    
    if type(energies) == 'list':
        energies = np.array(energies)
    
    nen = energies.shape[0]
    if energies.shape[1]==2:
        pho_band=(energies[:,0]+energies[:,1])/2.
    else:
        pho_band = energies
    ##constant
    mc2 = 510.98 # keV
    clight = 2.9979e10 # speed of light cm/s
    au = 1.496e13 #cm
    z = 1.2 # average atomic number
    keV2K = 1.16e7 #
    KK=2.6e-18 # cm^2 KeV^2, take Coulomg logarithm to be 20
    r0 = 2.8179e-13 #classical electron radius, in cm
    
    ##read parameters
    Ndot = total_eflux*1e35
    kappa_index = kappa_index
    Tk = kappa_temperature
    ehigh = high_cut

    np_ =plasma_d*1e10
    Tloop=loop_temp
    L= length*1e8
    
    
    ######
    maxfcn = 2048 # maximum number of points for the Gaussian quadrature integration

    fcoeff = 1/(4 * np.pi * (au**2)) #numerical coefficient for the photon flux
    # 1/mc2 is due to cross-section

    rerr = 1e-3               #relative error
#    intsum = dblarr(nen)          ; integral sum
    nlim = 12 # maximum possible order of the Gaussian quadrature integration is 2^nlim. [Updated]
    conv = True # convergence flag

    intsum = np.zeros(nen)
    for ires in range(1,nlim+1):
        lastsum = np.copy(intsum)
        npoint = int(2**ires)
        if npoint > maxfcn:
            conv = False ##Not converge
            intsum = intsum
        # electron and photon energy
        #ele_band = np.linspace(pho_band[0],ehigh,num=nen)
        
        ##Energy loss
        #dedt=collisional_loss(ele_band)
        #decoeff = 4*np.pi*r0**2*clight*mc2 #coefficient apply to the dedt
        #v_c=sqrt(2 * ele_band / mc2 + ele_band**2 / mc2**2)/(1 + ele_band/mc2) * clight #electron speed
        #dedx = dedt / v_c # dE/dx=(dE/dt)v
        ##bremsstrahlung cross section
        #cross=bremsstrahlung_cross_section(ele_band, pho_band, z=1.2)
        #nVF_E= Ndot*(1/dedx)*((kappa_index-1)/(1+ele_band/(kappa_index*Tk))**(kappa_index-2)\
        #-(kappa_index-2)/(1.+ele_band/(kappa_index*Tk))**(kappa_index-1.))
        
        #integrand = nVF_E* cross/mc2
        #a_lg = np.log10(pho_band) #Lower limit of integration
        #b_lg = np.log10(np.full_like(a_lg, ehigh)) #Upper limit of integration
        a = pho_band #Lower limit of integration
        b = np.full_like(a, ehigh) #Upper limit of integration
        intsum = gauss_legendre(
            kap_get_integrand,
            a,
            b,
            n=npoint,
            func_kwargs={
                'Ndot':Ndot,
                'kappa_index':kappa_index,
                'Tk':Tk,
                'photon_energy': pho_band,
                'z': z,
            },
        )
        ## Convergence criterion
        diff = np.abs(intsum - lastsum) / intsum
        ind = np.where (diff > rerr)[0]
        
        # If all point have reached criterion return value and flags
        if ind.size == 0:
            intsum = intsum
            conv = True

    if not conv:
        print('f_thick_warm_kappa_f.pro: Not converged')

    flux1 =  intsum * fcoeff #nonthermal-conponent

    #thermal component
    ll=Tloop**2/(2*KK*np_) # collisional stoping distance for electrons of Tloop energy
    Emin=Tloop*3*(5*ll/L)**4
    EM_add=3*np.pi/2/KK/clight*np.sqrt(mc2/8)*Tloop**2/np.sqrt(Emin)*Ndot
    EM49  =EM_add*1e-49 # EM in units of 10^49 cm^(-3)
    #print,'Emin =', Emin, 'Additional EM =', EM49
    #mean free path of the lowest energy in the injected spectrum
    L_min = Tk**2/(2*KK*np_)/3. 
    if L_min > L:
        print('WARNING: min mean free path > loop length: Emin/kT = {:.1}'.format(Emin/Tloop))
    # For low density /or large low energy cut off in the injected spectrum/
    # thermalization is less important

    Tloop_mk = Tloop * keV2K /1e6 ## MK
    EM46 = EM49 * 1e3
    ###The range should be within [1,200] keV for thermal
    Fvth_add = np.zeros(nen)
    if energies.shape[1]==2:
        locmin = np.min(np.where (energies[:,0] > 1))
        locmax = np.max(np.where (energies[:,1] < 200.1))
        energies_the = energies [locmin:locmax,:]
    else:
        locmin = np.min(np.where (energies > 1))
        locmax = np.max(np.where (energies < 200.1))
        energies_the = energies [locmin:locmax]

    if Emin > 0.1: ##Means not exact work, just a warning
        print('WARNING: Emin > 0.1 keV,  Setting Emin to 0.1')
        Emin = 0.1
    
    
    Fvth_add[locmin:locmax] = fitting_legacy.photon_models_for_fitting.f_vth(temperature = Tloop_mk,emission_measure46=EM46,energies=energies_the)

    flux = flux1 + Fvth_add ##Total photon flux
    return flux

