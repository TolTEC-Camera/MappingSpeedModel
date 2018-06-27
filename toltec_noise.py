import numpy as nm
import pylab as pl

verbose = True

# load loading data
tmp = nm.loadtxt('LMT_am_runs/DJF_50.dat')
f_GHz  = tmp[1:,0]
atm_tx = tmp[1:,2]
T_atm  = tmp[1:,3]

def calculate_net_and_nep(f_GHz,system_efficiency,T_det,atmosphere_transmission, dish_diameter, P0_override=False):
    # constants [SI]
    h = 6.626068e-34
    k = 1.3806503e-23
    c = 299792458.0

    # convert from GHz
    f = f_GHz*1.0e9

    df = nm.mean(nm.diff(f))

    # calculate power vs frequency for RJ source
    P0 = k*T_det*df

    shot = 2*k*T_det*h*f*df
    wave = 2*(P0)**2.0 / df
    # this nep will be "low" because bad optical efficiency has reduced the as-seen loading
    nep = nm.sqrt(shot+wave)# watts root second

    # do the NET_CMB integral (Lueker's Thesis, 2.23)
    # in this scaling, bad optical efficiency makes the NET high again
    # since NEP_shot drops as sqrt(OE), but this equation embiggens as OE
    # that means NET does scale up as sqrt(OE) as expected
    net_integrand = ((h*f)/(2.725))**2.0 * (1.0/k) * ( nm.exp( (h*f)/(k*2.725) ) / ((nm.exp( (h*f)/(k*2.725) ) - 1)**2.0) )
    net = (nep / (nm.sqrt(2.0) * system_efficiency * net_integrand * df) )*1.0e6 # uK root second 

    # scale the NEP 
    nep = nep*1.0e18 #attowatts root hz

    # multiply both by the fact that detector noise = 0.3 * blip
    nep = nep*nm.sqrt(1.3)
    net = net*nm.sqrt(1.3)

    # calculate neft using the dish size
    nefd = (nep*1e-18*2.0e26*1e3) / ((nm.pi*(dish_diameter/2.0)**2)*df) # mJrtHz
    # and system efficieicny since nep is the at-detector number
    nefd = nefd/system_efficiency
    # and convert from rtHz to rts
    nefd = nefd/nm.sqrt(2.0)

    # scale NET and NEFD up by atmospheric transmission
    net = net/atmosphere_transmission
    nefd = nefd/atmosphere_transmission

    # net, nefd, and nep are all arrays now
    # https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Dealing_with_variance
    # sigma^2 = 1 / sum(sigma^(-2))
    net = nm.sqrt(1.0 / nm.nansum(net**(-2.0)))
    nefd = nm.sqrt(1.0 / nm.nansum(nefd**(-2.0)))
    # nep is sum of squares
    nep = nm.sqrt(nm.nansum(nep**2.0))
    # power just adds
    P0 = nm.nansum(P0)

    return P0,net,nefd,nep

# detector efficiency and passband
detector_optical_efficiency = 0.7
fcent_nom = 273

# load the passband
tmp = nm.load('model_passbands.npz')
f_GHz = tmp['f_GHz']
if fcent_nom==143:
    passband = tmp['band_150'] + 1.0e-10
    fmin = 128.0
    fmax = 170.0
elif fcent_nom==214:
    passband = tmp['band_220'] + 1.0e-10
    fmin = 195.0
    fmax = 245.0
elif fcent_nom==273:
    passband = tmp['band_280'] + 1.0e-10
    fmin = 245.0
    fmax = 310.0

pl.ion()
pl.figure(100)
pl.plot(f_GHz,passband)
pl.xlabel('[GHz]')
pl.ylabel('Passband')

# aperture efficiency for 1flambda horns
aperture_efficiency = 0.35

# ruze scattering formula
primary_mirror_optical_efficiency = nm.exp(-((4.0*nm.pi*(76e-6))/(299792458/(f_GHz*1.0e9)))**2)

## add a factor for the primary mirror emission
## noting that not all of the sky loading is transmitted through the primary
#T_tot = T_atm*primary_mirror_optical_efficiency + 273.0*0.06

# THIS approach says that both the "transmitted" light by the primary, and the scattered light
# all make it out to the sky anyway, so it's all T_atm if the mirror is involved
T_tot = T_atm + 273.0*0.06

# add a factor for 3 mirror emissions, assuming they have perfect transmission
T_tot += 3*290.0*0.01
# this is now the temperature incident onto the cryostat window
# scale this down by detector optical efficiency to make the temperature at the detector
T_det = T_tot*detector_optical_efficiency*passband*aperture_efficiency

telescope_diameter = 48.0

# calculate system efficiency
system_efficiency = detector_optical_efficiency*passband*aperture_efficiency*primary_mirror_optical_efficiency


# single-detector noise
# this calculates the loading assuming the T_det is the at-detector loading
# and then calculates NET and NEFD by scaling that noise by the system efficiency and the atmospheric transmission
P0,net,nefd,nep_at_detector = calculate_net_and_nep(f_GHz,
                                                    system_efficiency,
                                                    T_det,
                                                    atm_tx,
                                                    telescope_diameter)
if verbose:
    print('Center Frequency = '+str(fcent_nom)+' GHz')
    print('Loading = '+str(P0*1.0e12)+' pW')
    print('NEP at Detector = '+str(nep_at_detector)+' aWrts')
    print('NET_CMB = '+str(net)+' uKrts')
    print('NEFD = '+str(nefd)+' mJrts')

# convert nefd into mapping speed of TolTEC instrument
if fcent_nom==143:
    fwhm_arcsec = 9.5
    Ndet = 450*2
elif fcent_nom==214:
    fwhm_arcsec = 6.3
    Ndet = 900*2
elif fcent_nom==273:
    fwhm_arcsec = 5.0
    Ndet = 1800*2
else:
    print('FAIL: Choose 143 or 214 or 273 GHz as a center frequency')

# print noise scaled by Ndet
if verbose:
    print(' ')
    print('Total NET_CMB = '+str(net/nm.sqrt(Ndet))+' uKrts')
    print('Total NEFD = '+str((1000.0*nefd)/nm.sqrt(Ndet))+' uJrts')

# central construct is:
#     if you have a nefd of X in mJrts
#     and you have a small beam such that N beams tile 1 square degree
#     if you spend one second on each beam, across the square degree
#     for a total of N seconds, you will have a pseudo-mapping-speed of
#         (1/X)^2 deg^2/mJ^2/(N seconds) per det
#     say that (N seconds) is less than an hour, you could cover more sky in that hour
#     so to convert to /hour you multiply, yielding
#         (1/X)^2 * (3600 seconds)/(N seconds) deg^2/mJ^2/hour per det
#     the mapping speed also "goes as" the detector count, since it's mJ^2 units
#         (1/X)^2 * Ndet * (3600 seconds)/(N seconds) deg^2/mJ^2/hour per instrument
# putting that into practice
N_beams_in_sqdeg = (60.0*60.0)**2 / ((fwhm_arcsec**2)*(nm.pi/(4.0*nm.log(2.0)))) # do the 2 pi sigma^2 over 2.355 (the log(2) thing) and see if it matches the NIKA2 thing
mapping_speed = (1.0/nefd)**2 * Ndet * (3600.0/N_beams_in_sqdeg)

if verbose:
    print(' ')
    print('Mapping Speed, Raw = '+str(mapping_speed)+' deg^2/mJ^2/hr')
    print('Mapping Speed, Downscaled = '+str(mapping_speed*(26.0/184.0))+' deg^2/mJ^2/hr')

# check what the map RMS would be for a 100 square degree, 100 hour survey
map_rms = ((mapping_speed*(26.0/184.0)) * 100.0 / 100.0)**(-0.5)

# print everything
if verbose:
    print('Map rms 90: '+str(((mapping_speed*(26.0/184.0)) * 100.0 / 90.0)**(-0.5)))
    print('Map rms 2: '+str(((mapping_speed*(26.0/184.0)) * 100.0 / 2.0)**(-0.5)))
    print('Map rms 0.9: '+str(((mapping_speed*(26.0/184.0)) * 100. / 0.9)**(-0.5)))
    print('Map rms 100: '+str(((mapping_speed*(26.0/184.0)) * 100.0 / 100.0)**(-0.5)))