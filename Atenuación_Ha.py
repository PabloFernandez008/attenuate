# Este código toma datos de un archivo del Millenium para H_alpha
# y dibuja, para un redshift z dado, la función de luminosidad sin atenuar,
# la atenuada por Millenium, la atenuada por Calzetti y la de por Cardelli
# Pablo Fernández, 14/05/2021

import read_jc_obs as jc

from Att_Calzetti import att_Cal_1, att_Cal_2
from Att_Cardelli import att_Car_1, att_Car_2

from matplotlib import pyplot as plt  # para usar legend
from pylab import plot, errorbar, xlabel, ylabel, title, xlim, ylim
from numpy import zeros, size, histogram, log10, append, arange, loadtxt

###############################################################################

unitsh = True  # True: todo con unidades de h
hu = 0.6777     # parámetro Hubble

###############################################################################

# Toma un redshift y el archivo correspondiente
redshift = 0.69
file = 'LHa_45_200.txt'
waveA = 6563        # Longitud de onda de H_alpha en Angstroms

###############################################################################

# PARÁMETROS DEL PLOT DEL HISTOGRAMA
lbox = 200     # lado de la caja en Mpc/h
vol = lbox**3  # volumen en Mpc^3/h^3

if not unitsh:  # volumen sin h
    vol = (lbox/hu)**3
    
lmin = 39.
lmax = 43.
dl = 0.1
lbins = arange(lmin,lmax,dl)
lhist = lbins
vbins = append(lbins,lmax) 
    
types = ['sin atenuar','Millenium', 'Calzetti', 'Cardelli']
colors = ['blue','green','red','cyan']

    
# DATOS EXPERIMENTALES DEL Ha: no disponemos

# DATOS MILLENIUM: extracción de luminosidades y parámetros de las galaxias

data = loadtxt(file, float, delimiter=',')

if unitsh:
    L_Ha_noext = log10(data[:, 0])+40   # para que esté en log10(L[erg/s h^-2])
    L_Ha_ext = log10(data[:, 1])+40
else: # sin la h
    L_Ha_noext = log10(data[:, 0])+40 -2*log10(hu)
    L_Ha_ext = log10(data[:, 1])+40 -2*log10(hu)

Mcold_disc = data[:, 2]
metals_coldgas = data[:, 3]
rdisk = data[:, 4]

Z_disc = metals_coldgas/Mcold_disc
if unitsh:
    Mcold_disc = Mcold_disc*1e10       # (en masas solares/h)
    rhalf_mass_disc = rdisk            # (en Mpc/h)
else:  # sin la h
    Mcold_disc = Mcold_disc*1e10/hu
    rhalf_mass_disc = rdisk/hu

# PLOT HISTOGRAMA

# v_datos es una matriz que guardará todos los tipos de luminosidades
v_datos = zeros((size(types),size(L_Ha_noext)))

v_datos[0,:] = L_Ha_noext  # sin atenuar
v_datos[1,:] = L_Ha_ext    # atenuado por Millenium
    
    
for j in range(size(L_Ha_noext)):
    
    attenuation_disc_Cal = att_Cal_1(waveA, Mcold_disc[j], rhalf_mass_disc[j], Z_disc[j])
    attenuation_disc_Car = att_Car_1(waveA, Mcold_disc[j], rhalf_mass_disc[j], Z_disc[j])
    
    # atenuación en el continuo
    v_datos[2,j] = att_Cal_2(v_datos[0,j], attenuation_disc_Cal, 1)
    v_datos[3,j] = att_Car_2(v_datos[0,j], attenuation_disc_Car, 1)
        
for i in range(size(types)):
    H, bins_edges = histogram(v_datos[i],bins=vbins)
    y = log10(H/(vol*dl))
         
    plot(lhist,y,color=colors[i],label=types[i])    
    
leg = plt.legend(loc=0)

# PLOT DATOS EXPERIMENTALES: no hay

title('z = %.2f' %redshift) 

if unitsh:
    xlabel("log$_{10}$(LH$\u03B1$[erg s$^{-1}$ h$^{-2}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)
else:
    xlabel("log$_{10}$(LH$\u03B1$[erg s$^{-1}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$])",size=15)
    
xlim(40.5,43)
ylim(-6,-1.5)

plt.savefig('Ha_45_200.pdf')