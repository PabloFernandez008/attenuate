# Este código dibuja la función de luminosidad a 4 z distintos del OII con 
# atenuación de Calzetti tanto con el factor extra de Saito como sin él
# Pablo Fernández, 14/05/2021


import read_jc_obs as jc

from Att_Calzetti import att_Cal_1, att_Cal_2

from matplotlib import pyplot as plt  # para usar legend
from pylab import figure, plot, errorbar, xlabel, ylabel, xlim, ylim, title
from numpy import zeros, size, histogram, log10, pi, sqrt, exp, append, arange, loadtxt



###############################################################################

unitsh = True  # True: todo con unidades de h
hu = 0.6777 # parámetro Hubble

###############################################################################

# Toma un redshift y el archivo correspondiente
#z = [0.99, 0.69, 0.51, 0.17]
#files = ['LOII_41_200.txt', 'LOII_45_200.txt', 'LOII_48_200.txt', 'LOII_56_200.txt']
z = [0.17]
files = ['LOII_56_200.txt']
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
    
types = ['Sin atenuar','Calzetti', 'Calzetti X Saito']
lines = ['-','-.','--']

fig = 0 # parámetro de control de las 4 figuras

for redshift in z:
    
    # DATOS EXPERIMENTALES DEL OII
    obs_path = '/Users/pablofernandez/Documents/5º curso/TFG2/Datos/Comparat_2016/'

    ox, oy, el, eh = jc.read_jc_lf(obs_path,redshift,
                               infile='O2_3728-data-summary-Planck15.txt',h2unitL=True)
    
    if not unitsh:   # ahora, sin h^-2 ni h^-3
        ox, oy, el, eh = jc.read_jc_lf(obs_path,redshift,
                               infile='O2_3728-data-summary-Planck15.txt',h2unitL=False)

    Merror = zeros((2,size(el)))  # guarda las barras de error
    Merror[0,:]=el
    Merror[1,:]=eh
    
    # DATOS MILLENIUM: 4 ARCHIVOS DIFERENTES
    data = loadtxt(files[fig], float, delimiter=',')

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
    
    # PLOT HISTOGRAMA
    
    fig = fig + 1
    figure(fig)

    v_datos = zeros((size(types),size(L_OII_noext)))
    # sin atenuar
    v_datos[0,:] = L_OII_noext
        
        
    for j in range(size(L_OII_noext)):
        
        attenuation_disc = att_Cal_1(waveA, Mcold_disc[j], rhalf_mass_disc[j], Z_disc[j])
        
        # atenuación en el continuo
        v_datos[1,j] = att_Cal_2(v_datos[0,j], attenuation_disc, 1)
        # atenuación con factor de Saito
        v_datos[2,j] = att_Cal_2(v_datos[0,j], attenuation_disc, 5/(redshift+2.2))
            
    for i in range(size(types)):
        H, bins_edges = histogram(v_datos[i],bins=vbins)
        y = log10(H/(vol*dl))
             
        plot(lhist,y,linestyle=lines[i],label=types[i])    
        
    leg = plt.legend(loc=0)
    
    # PLOT DATOS EXPERIMENTALES
    
    errorbar(ox, oy, Merror, fmt='o', color='grey', capsize=2)
    
    title('z = %.2f' %redshift) 
    
    if unitsh:
        xlabel("log$_{10}$(LOII[erg s$^{-1}$ h$^{-2}$]])",size=15)
        ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$ h$^{-3}$])",size=15)
    else:
        xlabel("log$_{10}$(LOII[erg s$^{-1}$]])",size=15)
        ylabel("log$_{10}$($\Phi$[Mpc$^{-3}$ dex$^{-1}$])",size=15)
        
    xlim(40,43)
    ylim(-6,-1.)
    
    plt.savefig('Saito_56_200.pdf')