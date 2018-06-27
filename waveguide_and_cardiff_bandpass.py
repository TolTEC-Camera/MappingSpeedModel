import numpy as nm
import matplotlib.pyplot as pl
import scipy
import scipy.interpolate

def calc_band_tophat(f_GHz_in,f_lo,f_hi):
    band = nm.zeros_like(f_GHz_in)+1.0e-20

    band[(f_GHz_in > f_lo) & (f_GHz_in < f_hi)] = 1

    return band

def calc_scattering_matrix_for_waveguide(frequencies,length,a):

    # define the single-frequency model
    def single_freq_model_matmultiply(freq,length,a):
        mu = 4.0*nm.pi*1e-7
        c = 299792458.0
        epsilon =  1.0 / (mu*c**2)

        # convert freq from GHz
        freq = freq*1.0e9

        # find the impedance of the start guide for TE11 mode in circular waveguide
        Z_space = 119.9169832*nm.pi
        omega = 2.0*nm.pi*freq
        k = omega/c
        kc = 1.8411837813406593/(a*1.0e-3)
        Z = Z_space*(k / nm.sqrt(k**2.0 - kc**2.0 + 0.0j))
        bl=nm.sqrt(k**2.0 - kc**2.0 + 0.0j)*(length*1.0e-3)

        #print 'lambda_g '+str((2.0*nm.pi)/nm.sqrt(k**2.0 - kc**2.0 + 0.0j))
        #print 'f_c '+str( 1.0e-9*(c / ((2.0*nm.pi)/kc)))

        # calculate transmission and reflection
        A = nm.cos(bl); B = 1j*Z*nm.sin(bl); C = (1j/Z)*nm.sin(bl); D = nm.cos(bl)
        # from Pozar, for identical start and stop impedances
        Z = -Z
        t = 2.0 / (A + B/Z + C*Z + D)
        r = (A + B/Z - C*Z - D)/(A + B/Z + C*Z + D)

        # calculate amplitude and phase
        T_amp = t*nm.conj(t)
        T_phase = nm.angle(t)
        R_amp = r*nm.conj(r)
        R_phase = nm.angle(r)

        # return
        return T_amp,T_phase,R_amp,R_phase

    T_amp = nm.zeros_like(frequencies)
    T_phase = nm.zeros_like(frequencies)
    R_amp = nm.zeros_like(frequencies)
    R_phase = nm.zeros_like(frequencies)

    for i in range(len(frequencies)):
        T_amp[i],T_phase[i],R_amp[i],R_phase[i] = single_freq_model_matmultiply(frequencies[i],length,a)

    return T_amp,T_phase,R_amp,R_phase

def calc_waveguide_cuton(f_GHz_in,f_lo):
    # lambda_g at 148 GHz (0.00378334077804+0j) mm
    # f_c 125.0 GHz
    # a = 0.7027938657899518 mm
    # always 2.5 mm of waveguide
    T_amp,T_phase,R_amp,R_phase = calc_scattering_matrix_for_waveguide(f_GHz_in,2.5,0.7027938657899518*(125.0/f_lo))

    print(0.7027938657899518*(125.0/f_lo))

    return T_amp

def load_cardiff_dichroic(f_GHz_in,isLowpass,desired_cutoff):
    import pandas as pd
    if isLowpass:
        df = pd.read_excel('TolTEC Detector Array Testing Filter Feb 2017.xlsx',sheetname='T0710R21',names=['transmission'],parse_cols=1)
    else:
        df = pd.read_excel('TolTEC Detector Array Testing Filter Feb 2017.xlsx',sheetname='T0715R1',names=['transmission'],parse_cols=1)
    nu_twiddle = nm.array(df.index)
    f_GHz = nu_twiddle*1e2*299792458.0*1.0e-9
    # scale the frequency
    f_GHz *= desired_cutoff/176.9
    cardiff_filter_transmission = nm.array(df)[:,0]
    
    # this isn't at the desired frequencies
    func = scipy.interpolate.interp1d(f_GHz,cardiff_filter_transmission,kind='linear',fill_value = 'extrapolate',assume_sorted=True)

    return func(f_GHz_in)+1.0e-20

def load_cardiff_scaled(f_GHz_in,f_hi_base,f_hi):
    d = {175:'T1098R25',186:'T1088R28',255:'T1625R10',285:'T1394R29',325:'S3276R5'}
    import pandas as pd
    df = pd.read_excel('TolTEC Detector Array Testing Filter Feb 2017.xlsx',sheetname=d[f_hi_base],names=['transmission'],parse_cols=1)
    nu_twiddle = nm.array(df.index)
    f_GHz = nu_twiddle*1e2*299792458.0*1.0e-9
    cardiff_filter_transmission = nm.array(df)[:,0]

    # scale the 300 GHz cutoff, to whatever is desired
    f_GHz_scaled = f_GHz * (float(f_hi)/float(f_hi_base))

    # this isn't at the desired frequencies
    func = scipy.interpolate.interp1d(f_GHz_scaled,cardiff_filter_transmission,kind='linear',fill_value = 'extrapolate',assume_sorted=True)

    return func(f_GHz_in)+1.0e-20
