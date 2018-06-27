import os

# delete the existing files
os.system('rm LMT_am_runs/*')


# run each am configuration
os.system('am ../LMT_am_model/am_models/LMT_DJF_25.amc 0 GHz 550 GHz 10 MHz 30 deg 1.0 > LMT_am_runs/DJF_25.dat')
os.system('am ../LMT_am_model/am_models/LMT_DJF_25.amc 0 GHz 550 GHz 10 MHz 30 deg 1.02 > LMT_am_runs/DJF_25_1p02.dat')
os.system('am ../LMT_am_model/am_models/LMT_DJF_50.amc 0 GHz 550 GHz 10 MHz 30 deg 1.0 > LMT_am_runs/DJF_50.dat')
os.system('am ../LMT_am_model/am_models/LMT_DJF_50.amc 0 GHz 550 GHz 10 MHz 30 deg 1.02 > LMT_am_runs/DJF_50_1p02.dat')
os.system('am ../LMT_am_model/am_models/LMT_DJF_75.amc 0 GHz 550 GHz 10 MHz 30 deg 1.0 > LMT_am_runs/DJF_75.dat')
os.system('am ../LMT_am_model/am_models/LMT_DJF_75.amc 0 GHz 550 GHz 10 MHz 30 deg 1.02 > LMT_am_runs/DJF_75_1p02.dat')