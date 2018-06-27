import numpy as nm
import matplotlib.pyplot as pl
import waveguide_and_cardiff_bandpass as w

# always skip first data point, because f=0 causes problem with CMB integral

# load atmosphere DJF
#              Zenith               LOS
# 25        = (980.799 um_pwv)     (1132.53 um_pwv)
# 25 * 1.02 = (1000.38 um_pwv)     (1155.14 um_pwv)
# 50        = (1942.02 um_pwv)     (2242.45 um_pwv) 
# 50 * 1.02 = (1980.84 um_pwv)     (2287.28 um_pwv)
# 75 =        (3680.74 um_pwv)     (4250.15 um_pwv)
# 75 * 1.02 = (3754.38 um_pwv)     (4335.19 um_pwv)
tmp = nm.loadtxt('LMT_am_runs/DJF_25.dat')
f_GHz  = tmp[1:,0]
atm_tx_25 = tmp[1:,2]
T_atm_25  = tmp[1:,3]
tmp = nm.loadtxt('LMT_am_runs/DJF_50.dat')
f_GHz  = tmp[1:,0]
atm_tx_50 = tmp[1:,2]
T_atm_50  = tmp[1:,3]
tmp = nm.loadtxt('LMT_am_runs/DJF_75.dat')
f_GHz  = tmp[1:,0]
atm_tx_75 = tmp[1:,2]
T_atm_75  = tmp[1:,3]

tmp = nm.loadtxt('LMT_am_runs/DJF_25_1p02.dat')
f_GHz  = tmp[1:,0]
atm_tx_25_hi = tmp[1:,2]
T_atm_25_hi  = tmp[1:,3]
tmp = nm.loadtxt('LMT_am_runs/DJF_50_1p02.dat')
f_GHz  = tmp[1:,0]
atm_tx_50_hi = tmp[1:,2]
T_atm_50_hi  = tmp[1:,3]
tmp = nm.loadtxt('LMT_am_runs/DJF_75_1p02.dat')
f_GHz  = tmp[1:,0]
atm_tx_75_hi = tmp[1:,2]
T_atm_75_hi  = tmp[1:,3]

delta_T_50um_25 = 50.0 * ((T_atm_25_hi - T_atm_25) / (1155.14 - 1132.53))
delta_T_50um_50 = 50.0 * ((T_atm_50_hi - T_atm_50) / (2287.28 - 2242.45))
delta_T_50um_75 = 50.0 * ((T_atm_75_hi - T_atm_75) / (4335.19 - 4250.15))

pl.ion()
pl.figure(1,figsize=[18,8])
pl.clf()

# 280 band
# 245 GHz dichroic
lowpass_245_dichroic = w.load_cardiff_dichroic(f_GHz, True, 245.0)
highpass_245_dichroic = w.load_cardiff_dichroic(f_GHz, False, 245.0)
# waveguide cuton
wg_233 = w.calc_waveguide_cuton(f_GHz, 233.0)
# 315 GHz lowpass
lp_315 = w.load_cardiff_scaled(f_GHz, 285, 315.0)

# multiply the total band
total_280_band = wg_233 * highpass_245_dichroic * lp_315 * lp_315

pl.subplot(233)
pl.plot(f_GHz,wg_233,'k:',label='233 GHz Waveguide')
pl.plot(f_GHz,highpass_245_dichroic,'b',label='245 GHz Dichroic\nHighpass')
pl.plot(f_GHz,lp_315*lp_315,'r:',label='315 GHz Lowpass\n(Two of them)\n(Yields 310 GHz rolloff)')
pl.plot(f_GHz,total_280_band,'k',label='Passband')
pl.legend(loc='center left')
#pl.xlabel('[GHz]')
pl.ylabel('Tranmission')
pl.title('280 GHz Band Definition')
pl.xlim([100,550])
pl.ylim([-0.05,1.05])
pl.grid()




###### 220 band
# waveguide cuton
wg_195 = w.calc_waveguide_cuton(f_GHz,195.0)
# 255 GHz lowpass
lp_255 = w.load_cardiff_scaled(f_GHz,255,255.0)
# 185 GHz dichroic
lowpass_185_dichroic = w.load_cardiff_dichroic(f_GHz,True,185.0)
highpass_185_dichroic = w.load_cardiff_dichroic(f_GHz,False,185.0)

# multiply the total band
total_220_band = highpass_185_dichroic*wg_195*lowpass_245_dichroic*lp_255

# illustrate this
pl.subplot(232)
pl.plot(f_GHz,highpass_185_dichroic,'b:',label='185 GHz\nDichroic\nHighpass')
pl.plot(f_GHz,wg_195,'k:',label='195 GHz\nWaveguide')
pl.plot(f_GHz,lowpass_245_dichroic,'b',label='245 GHz\nDichroic\nLowpass')
pl.plot(f_GHz,lp_255,'r:',label='255 GHz\nLowpass')
pl.plot(f_GHz,total_220_band,'k',label='Passband')
pl.legend(loc='center right')
#pl.xlabel('[GHz]')
pl.ylabel('Tranmission')
pl.title('220 GHz Band Definition')
pl.xlim([100,550])
pl.ylim([-0.05,1.05])
pl.grid()




###### 150 band
# waveguide cuton
wg_128 = w.calc_waveguide_cuton(f_GHz,128.0)
# 175 GHz lowpass
lp_175 = w.load_cardiff_scaled(f_GHz,175,170.0)

# multiply the total band
total_150_band = lowpass_245_dichroic*wg_128*lp_175*lowpass_185_dichroic

# illustrate this
pl.subplot(231)
pl.plot(f_GHz,wg_128,'k:',label='128 GHz Waveguide')
pl.plot(f_GHz,lp_175,'r:',label='170 GHz Lowpass')
pl.plot(f_GHz,lowpass_185_dichroic,'b:',label='185 GHz Dichroic\nLowpass')
pl.plot(f_GHz,lowpass_245_dichroic,'b',label='245 GHz Dichroic\nLowpass')
pl.plot(f_GHz,total_150_band,'k',label='Passband')
pl.legend(loc='center right')
#pl.xlabel('[GHz]')
pl.ylabel('Tranmission')
pl.title('150 GHz Band Definition')
pl.xlim([100,550])
pl.ylim([-0.05,1.05])
pl.grid()










###### logscale and atmosphere
pl.subplot(234)
pl.plot(f_GHz,total_150_band,'cornflowerblue')
pl.plot(f_GHz,total_220_band,'goldenrod')
pl.plot(f_GHz,total_280_band,'seagreen')
pl.plot(f_GHz,total_150_band+total_220_band+total_280_band,':k')
pl.xlabel('[GHz]')
pl.ylabel('Tranmission')
pl.xlim([100,550])
pl.ylim([-0.05,1.05])
pl.title('All Three Bands, Linear Scale')
pl.grid()


pl.subplot(235)
pl.plot(f_GHz,10*nm.log10(total_150_band),'cornflowerblue')
pl.plot(f_GHz,10*nm.log10(total_220_band),'goldenrod')
pl.plot(f_GHz,10*nm.log10(total_280_band),'seagreen')
pl.xlabel('[GHz]')
pl.ylabel('[dB]')
pl.xlim([100,550])
pl.ylim([-60,5])
pl.title('All Three Bands, Log Scale')
pl.grid()

pl.subplot(236)
pl.plot(f_GHz,total_150_band*T_atm_25,linestyle=':',color='cornflowerblue')
pl.plot(f_GHz,total_220_band*T_atm_25,linestyle=':',color='goldenrod')
pl.plot(f_GHz,total_280_band*T_atm_25,linestyle=':',color='seagreen')
pl.plot(f_GHz,T_atm_25,'0.7')
pl.plot(f_GHz,T_atm_50,'0.5')
pl.plot(f_GHz,T_atm_75,'0.3')
pl.xlabel('[GHz]')
pl.ylabel('[K]')
pl.xlim([100,550])
pl.ylim([-10,100])
pl.title('25-50-75 Percentile Atmosphere in Dec-Feb')
pl.grid()

pl.tight_layout()

pl.figure(2,figsize=[16,8])
pl.clf()
pl.subplot(311)
pl.plot(f_GHz,T_atm_25,'0.7')
pl.plot(f_GHz,T_atm_50,'0.7')
pl.plot(f_GHz,T_atm_75,'0.7')
pl.plot(f_GHz,total_150_band*T_atm_25,color='cornflowerblue',linestyle=':')
pl.plot(f_GHz,total_220_band*T_atm_25,color='goldenrod',linestyle=':')
pl.plot(f_GHz,total_280_band*T_atm_25,color='seagreen',linestyle=':')
pl.plot(f_GHz,total_150_band*T_atm_50,color='cornflowerblue',linestyle='--')
pl.plot(f_GHz,total_220_band*T_atm_50,color='goldenrod',linestyle='--')
pl.plot(f_GHz,total_280_band*T_atm_50,color='seagreen',linestyle='--')
pl.plot(f_GHz,total_150_band*T_atm_75,color='cornflowerblue',linestyle='-')
pl.plot(f_GHz,total_220_band*T_atm_75,color='goldenrod',linestyle='-')
pl.plot(f_GHz,total_280_band*T_atm_75,color='seagreen',linestyle='-')
pl.xlim([00,380])#550])
pl.ylim([-3,80])
pl.ylabel('Loading [K]')
pl.title('25-50-75 Percentile Atmosphere Loading in Dec-Feb at LMT Site')
pl.grid()

pl.subplot(312)
pl.plot(f_GHz,delta_T_50um_25,'0.7')
pl.plot(f_GHz,delta_T_50um_50,'0.7')
pl.plot(f_GHz,delta_T_50um_75,'0.7')
pl.plot(f_GHz,total_150_band*delta_T_50um_25,color='cornflowerblue',linestyle=':')
pl.plot(f_GHz,total_220_band*delta_T_50um_25,color='goldenrod',linestyle=':')
pl.plot(f_GHz,total_280_band*delta_T_50um_25,color='seagreen',linestyle=':')
pl.plot(f_GHz,total_150_band*delta_T_50um_50,color='cornflowerblue',linestyle='--')
pl.plot(f_GHz,total_220_band*delta_T_50um_50,color='goldenrod',linestyle='--')
pl.plot(f_GHz,total_280_band*delta_T_50um_50,color='seagreen',linestyle='--')
pl.plot(f_GHz,total_150_band*delta_T_50um_75,color='cornflowerblue',linestyle='-')
pl.plot(f_GHz,total_220_band*delta_T_50um_75,color='goldenrod',linestyle='-')
pl.plot(f_GHz,total_280_band*delta_T_50um_75,color='seagreen',linestyle='-')
pl.xlim([00,380])#550])
pl.ylim([-0.05,0.85])
pl.ylabel('dT [K]')
pl.title('2% Fluctuation on 25-50-75 Percentile Atmosphere')
pl.grid()

pl.subplot(313)
pl.plot(f_GHz,total_150_band,color='cornflowerblue',linestyle='-')
pl.plot(f_GHz,total_220_band,color='goldenrod',linestyle='-')
pl.plot(f_GHz,total_280_band,color='seagreen',linestyle='-')
pl.plot(f_GHz,total_150_band*atm_tx_75,color='cornflowerblue',linestyle=':')
pl.plot(f_GHz,total_220_band*atm_tx_75,color='goldenrod',linestyle=':')
pl.plot(f_GHz,total_280_band*atm_tx_75,color='seagreen',linestyle=':')
pl.plot(f_GHz,atm_tx_25,'0.7')
pl.plot(f_GHz,atm_tx_50,'0.5')
pl.plot(f_GHz,atm_tx_75,'0.3')
pl.xlabel('[GHz]')
pl.ylabel('Tranmission')
pl.xlim([00,380])#550])
pl.ylim([-0.05,1.05])
pl.title('25-50-75 Percentile Transmission')
pl.grid()

# save the resulting passbands to a file
nm.savez('model_passbands.npz', f_GHz=f_GHz,
	                            band_150=total_150_band,
	                            band_220=total_220_band,
	                            band_280=total_280_band)

pl.tight_layout()




pl.figure(100)
pl.clf()
pl.subplot(311)
pl.plot(f_GHz,T_atm_25,'0.7')
pl.plot(f_GHz,T_atm_50,'0.5')
pl.plot(f_GHz,T_atm_75,'0.3')
pl.ylabel('Antenna\nTemperature [K]')
pl.grid('on')
pl.xlim([0,400])
pl.subplot(312)
pl.plot(f_GHz,delta_T_50um_25,'0.7',label='1.0 mm PWV')
pl.plot(f_GHz,delta_T_50um_50,'0.5',label='2.0 mm PWV')
pl.plot(f_GHz,delta_T_50um_75,'0.3',label='3.5 mm PWV')
pl.legend()
pl.ylabel('2% PWV\nChange [K]')
pl.grid('on')
pl.xlim([0,400])
pl.subplot(313)
pl.plot(f_GHz,atm_tx_25,'0.7')
pl.plot(f_GHz,atm_tx_50,'0.5')
pl.plot(f_GHz,atm_tx_75,'0.3')
pl.ylabel('    \nTransmission')
pl.xlabel('[GHz]')
pl.grid('on')
pl.xlim([0,400])
pl.tight_layout()