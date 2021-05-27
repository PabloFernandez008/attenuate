# Función de Cardelli para meter atenuación al logaritmo de la luminosidad,
# necesita usar Cosmology.py y datos de las galaxias de Millenium
# Pablo Fernández

import Cosmology as cosmo
from numpy import log10, pi, sqrt, exp

########################## FUNCIÓN DE CARDELLI ################################


Rv = 3.1
Zsun = 0.0134       # Metalicidad solar, Asplund et al. 2009
s = 1.6             # Exponente s, si lambda>2000A
a_disc = 1.68       # En la fórmula de <NH>
mp = 1.67e-27       # Masa del protón en kg
costheta = 0.30          # scattering forward oriented
sectheta = 1./costheta
albedo = 0.56
Mpc_to_cm = 3.086e+24    # pasar de Mpc a cm
    #se usa porque rhalf_mass_dis viene dado en Mpc
Msun_to_kg = 1.989e+30   # masa del Sol en kg, 
    #se usa porque Mcold_disc viene dado en unidades de masas solares

cosmo.set_cosmology(0.307, 0.048, 0.693, 0.823, 0.96)   # definir cosmología de Planck


def cardelliIR(waveA):  # IR: 9090A < λ < 33300A
    wave = waveA/10000.     # la pasa a micrómetros
    x = 1./wave             # x es el inverso
    ax = 0.574*(x**1.61) 
    bx = -0.527*(x**1.61)
    return ax+bx/Rv

def cardelliOPT(waveA):  # Optical/NIR: 3030A < λ < 9090A
    wave = waveA/10000.   # la pasa a micrometros
    x = 1./wave           
    y = x-1.82
    ax = 1.+0.17699*y-0.50447*y**2-0.02427*y**3+0.72085*y**4+0.01979*y**5-0.77530*y**6+0.32999*y**7 
    bx = 1.41338*y+2.28305*y**2+1.07233*y**3-5.38434*y**4-0.62251*y**5+5.30260*y**6-2.09002*y**7
    return ax+bx/Rv


def att_Car_1(waveA, Mcold_disc, rhalf_mass_disc, Z_disc):
    
    if waveA < 9090:            # Dependiendo de λ, calcula Alambda/A_V
        Al_Av = cardelliOPT(waveA)
    else:
        Al_Av = cardelliIR(waveA)
        
    
    mean_col_dens_disc_log = log10(Mcold_disc)+log10(Msun_to_kg)-log10(1.4*mp*pi) \
        -2.*log10(a_disc*rhalf_mass_disc*Mpc_to_cm)
    
    # calcula el log10 de tau_lambda
    tau_disc = log10(Al_Av) + s*log10(Z_disc/Zsun) + mean_col_dens_disc_log - log10(2.1e+21)
    tau_disc = 10.**tau_disc    # deshacer el logaritmo
    
    al_disc = (sqrt(1.-albedo))*tau_disc     # calcula el a_lambda
    
    # calcula la atenuación
    attenuation_disc = -2.5*log10((1.-exp(-al_disc*sectheta))/(al_disc*sectheta))
    
    return attenuation_disc

def att_Car_2(logL, attenuation_disc, line_factor):
    
    # calcula el log10 de la L atenuada
    logL_attenuated = logL - 0.4*attenuation_disc*(line_factor)
    
    return logL_attenuated