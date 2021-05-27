# Este código toma datos de un archivo del Millenium para OII
# y dibuja, para un redshift z dado, la función de luminosidad sin atenuar,
# la atenuada con y sin el factor de Saito. Se crean 3 gráficas según el tipo
# de atenuación: Cardelli, Calzetti y Millenium
# Pablo Fernández, 14/05/2021

import read_jc_obs as jc

from Att_Calzetti import att_Cal_1, att_Cal_2
from Att_Cardelli import att_Car_1, att_Car_2

from matplotlib import pyplot as plt  # para usar legend
from pylab import figure, plot, errorbar, xlabel, ylabel, title, xlim, ylim
from numpy import zeros, size, histogram, log10, append, arange, loadtxt

###############################################################################

unitsh = True  # True: todo con unidades de h
hu = 0.6777     # parámetro Hubble

###############################################################################

# Toma un redshift y el archivo correspondiente
#z = [4.89, 2.24, 0.99, 0.51]
#files = ['LOII_21_200.txt', 'LOII_31_200.txt', 'LOII_41_200.txt', 'LOII_48_200.txt']
z = [2.24, 0.99, 0.51]
files = ['LOII_31_200.txt', 'LOII_41_200.txt', 'LOII_48_200.txt']
waveA = 3727        # Longitud de onda del OII en Angstroms

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
    
types = ['sin atenuar','Millenium', 'Calzetti', 'Calzetti con Saito', 
         'Cardelli', 'Cardelli con Saito']
colors = ['blue','green','red','pink']  # cada z se pinta con un color distinto


i = -1 #parámetro de control de archivo

for redshift in z:
    
    i = i + 1
    
    # DATOS MILLENIUM: extracción de luminosidades y parámetros de las galaxias
    
    data = loadtxt(files[i], float, delimiter=',')  # con esa z, toma el archivo correspondiente
    
    if unitsh:
        L_OII_noext = log10(data[:, 0])+40   # para que esté en log10(L[erg/s h^-2])
        L_OII_ext = log10(data[:, 1])+40
    else: # sin la h
        L_OII_noext = log10(data[:, 0])+40 -2*log10(hu)
        L_OII_ext = log10(data[:, 1])+40 -2*log10(hu)
        
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
    
    
    # v_datos es una matriz que guardará todos los tipos de luminosidades
    v_datos = zeros((size(types),size(L_OII_noext)))
    
    v_datos[0,:] = L_OII_noext  # sin atenuar
    v_datos[1,:] = L_OII_ext    # atenuado por Millenium
 
###############################################################################
       
    figure(1)   # Atenuación por Millenium
        
    H, bins_edges = histogram(v_datos[1],bins=vbins)
    y = log10(H/(vol*dl))
             
    plot(lhist,y,color=colors[i],linestyle='solid',label='z = %.2f' %redshift)    
        
    leg = plt.legend(loc=0)
    
###############################################################################    
    
    figure(2)   # Atenuación por Calzetti
     
    for j in range(size(L_OII_noext)):
        
        attenuation_disc_Cal = att_Cal_1(waveA, Mcold_disc[j], rhalf_mass_disc[j], Z_disc[j])
        
        # atenuación en el continuo
        v_datos[2,j] = att_Cal_2(v_datos[0,j], attenuation_disc_Cal, 1)
        v_datos[3,j] = att_Cal_2(v_datos[0,j], attenuation_disc_Cal, 5/(redshift+2.2))  # con Saito
            
    H, bins_edges = histogram(v_datos[2],bins=vbins)
    y = log10(H/(vol*dl))
    plot(lhist,y,color=colors[i],linestyle='solid',label='z = %.2f' %redshift)

    H, bins_edges = histogram(v_datos[3],bins=vbins)
    y = log10(H/(vol*dl))
    plot(lhist,y,color=colors[i],linestyle='dashed',label='z = %.2f con Saito' %redshift)     
        
    leg = plt.legend(loc=0)
    
#############################################################################3#
    
    figure(3)   # Atenuación por Cardelli
     
    for j in range(size(L_OII_noext)):
        
        attenuation_disc_Car = att_Car_1(waveA, Mcold_disc[j], rhalf_mass_disc[j], Z_disc[j])
        
        # atenuación en el continuo
        v_datos[4,j] = att_Car_2(v_datos[0,j], attenuation_disc_Car, 1)
        v_datos[5,j] = att_Car_2(v_datos[0,j], attenuation_disc_Car, 5/(redshift+2.2))
            
    H, bins_edges = histogram(v_datos[4],bins=vbins)
    y = log10(H/(vol*dl))
    plot(lhist,y,color=colors[i],linestyle='solid',label='z = %.2f' %redshift)
    
    H, bins_edges = histogram(v_datos[5],bins=vbins)
    y = log10(H/(vol*dl))
    plot(lhist,y,color=colors[i],linestyle='dashed',label='z = %.2f con Saito' %redshift)  
        
    leg = plt.legend(loc=0)
    
################################################################################    
# Por último, ponemos títulos y labels    
    
figure(1)
title('Atenuación Millenium') 
    
if unitsh:
    xlabel("log$_{10}$(LOII[erg s$^{-1}$ h$^{-2}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)
else:
    xlabel("log$_{10}$(LOII[erg s$^{-1}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$])",size=15)
    
xlim(40.5,43)
ylim(-6,-1.5)

plt.savefig('Millenium_z.pdf')

figure(2)
title('Atenuación por Calzetti') 
    
if unitsh:
    xlabel("log$_{10}$(LOII[erg s$^{-1}$ h$^{-2}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)
else:
    xlabel("log$_{10}$(LOII[erg s$^{-1}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$])",size=15)
    
xlim(40.5,43)
ylim(-6,-1.5)

plt.savefig('Calzetti_z.pdf')

figure(3)
title('Atenuación por Cardelli') 
    
if unitsh:
    xlabel("log$_{10}$(LOII[erg s$^{-1}$ h$^{-2}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)
else:
    xlabel("log$_{10}$(LOII[erg s$^{-1}$]])",size=15)
    ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$])",size=15)
    
xlim(40.5,43)
ylim(-6,-1.5)

plt.savefig('Cardelli_z.pdf')
