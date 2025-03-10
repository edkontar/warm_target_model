# warm_target_model for Hard X-ray solar flare spectroscopy
Warm target model to fit Hard X-ray spectra (available in SSW/OSPEX https://hesperia.gsfc.nasa.gov/ssw/packages/xray/idl/)

**f_thick_warm_kappa.pro** and **f_thick_warm_kappa.py**: Kappa distribution of accelerated electrons (see https://doi.org/10.3847/1538-4357/aafad3 and https://doi.org/10.3847/1538-4357/ad6a59 ) for application to RHESSI and SolO/STIX Hard X-ray spectroscopy.
The model assumes injected (accelerated) electron distribution in the form of kappa-distribution:

$$f_\kappa(v)=\frac{n_\kappa}{\pi^{3 / 2} v_{\mathrm{te}}^3 \kappa^{3 / 2}} \frac{\Gamma(\kappa)}{\Gamma\left(\kappa-\frac{3}{2}\right)}\left(1+\frac{v^2}{\kappa v_{\mathrm{te}}^2}\right)^{-\kappa},$$

where $v_{\mathrm{te}}=\sqrt{2 k_{\mathrm{B}} T_\kappa / m_e}$ is the thermal speed, and $\kappa=$ $\Gamma_{\text {coll }} / 2 D_0$ is $\frac{1}{2} \times$ the ratio of the collisional and diffusional coefficients in the Fokker-Planck equation describing electron acceleration (http://dx.doi.org/10.1088/0004-637X/796/2/142). So the accelerated electron rate spectrum $\dot{N}_\kappa(E)$ is

$$\dot{N} _ \kappa (E)=A v f_\kappa(v) 4 \pi v^2 \frac{d v}{d E}=A n_\kappa\left(\frac{8}{\pi m_e k_{\mathrm{B}} T_\kappa}\right)^{1 / 2} \times \frac{\Gamma(\kappa)}{\kappa^{3/2}\Gamma\left(\kappa-\frac{3}{2}\right)} \frac{E / k_{\mathrm{B}} T_\kappa}{\left(1+E / \kappa k_{\mathrm{B}} T_\kappa\right)^\kappa}\,$$

where $A$ is the cross-sectional area of the loop. The total rate of electron injection is

$$\dot{N}_ 0= \int_0^{\infty} \dot{N}_ \kappa(E) d E=2 A n_\kappa\left(\frac{2 k_{\mathrm{B}} T_\kappa}{m_e}\right)^{1 / 2} \times \frac{\kappa^{1 / 2}}{(\kappa-2) B\left(\kappa-\frac{3}{2}, \frac{1}{2}\right)}\,,$$

where $B(\cdot, \cdot)$ is the beta function and the function to determine

$$\dot{N}_ \kappa(E)=\left(\frac{N_0}{k_{\mathrm{B}} T_\kappa}\right) \frac{(\kappa-1)(\kappa-2)}{\kappa^2} \frac{E / k_{\mathrm{B}} T_\kappa}{\left(1+E / \kappa k_{\mathrm{B}} T_\kappa\right)^\kappa} .$$







