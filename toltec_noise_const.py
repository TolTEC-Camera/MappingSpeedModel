import numpy as nm
import pylab as pl

# always skip first data point, because f=0 causes problem with CMB integral

# load loading data
tmp = nm.loadtxt('LMT_am_runs/DJF_50.dat')
f_GHz  = tmp[1:,0]
atm_tx = tmp[1:,2]
T_atm  = tmp[1:,3]

# overwrite with constant to test integral
T_atm = 20.0*nm.ones_like(f_GHz)
atm_tx = 0.92*nm.ones_like(f_GHz)

def calculate_net_and_nep(f_GHz,system_efficiency,T_det,atmosphere_transmission, dish_diameter, P0_override=False):
    # constants [SI]
    h = 6.626068e-34
    k = 1.3806503e-23
    c = 299792458.0

    # convert from GHz
    f = f_GHz*1.0e9

    df = nm.mean(nm.diff(f))

    print('T_det')
    print(nm.max(T_det))

    # calculate power vs frequency for RJ source
    P0 = k*T_det*df

    shot = 2*k*T_det*h*f*df
    wave = 2*(P0)**2.0 / df
    print('Shot noise = '+str(nm.sqrt(nm.sum(shot))*1.0e18)+' aWrtHz')
    print('Wave noise = '+str(nm.sqrt(nm.sum(wave))*1.0e18)+' aWrtHz')
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

    ## multiply both by the fact that detector noise = 0.3 * blip
    #nep = nep*nm.sqrt(1.3)
    #net = net*nm.sqrt(1.3)

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
    net = nm.sqrt(1.0 / nm.sum(net**(-2.0)))
    nefd = nm.sqrt(1.0 / nm.sum(nefd**(-2.0)))
    # nep is sum of squares
    nep = nm.sqrt(nm.sum(nep**2.0))
    # power just adds
    P0 = nm.sum(P0)

    return P0,net,nefd,nep

# detector efficiency and passband
detector_optical_efficiency = 0.7
fcent_nom = 273
fmin = 245.0
fmax = 310.0

passband = nm.zeros_like(f_GHz) + 1e-10 # literally zero will make noise infinite, which may mess up the sum later
passband[(f_GHz>fmin) & (f_GHz<fmax)] = 1.0

pl.ion()
pl.figure(100)
pl.plot(f_GHz,passband)
pl.xlabel('[GHz]')
pl.ylabel('Passband')

# aperture efficiency for 1flambda horns
aperture_efficiency = 0.35

primary_mirror_optical_efficiency = 0.47

# scale this down by detector optical efficiency to make the temperature at the detector
T_det = T_atm*detector_optical_efficiency*passband*aperture_efficiency

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
print('Center Frequency = '+str(fcent_nom)+' GHz')
print('Loading = '+str(P0*1.0e12)+' pW')
print('NEP at Detector = '+str(nep_at_detector)+' aWrts')
print('NET_CMB = '+str(net)+' uKrts')
print('NEFD = '+str(nefd)+' mJrts')