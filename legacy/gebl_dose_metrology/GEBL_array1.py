import os
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt 
import array as arr
from scipy.signal import argrelextrema
from collections import defaultdict
from uncertainties import ufloat,unumpy as unp
import scipy as sp
from scipy.odr import Model, RealData, ODR
from AFM_GEBL_functions import *

plt.close("all"); plt.ion()
plt.rcParams.update({"text.usetex": True, "font.family": "sans-serif", "font.sans-serif": ["Computer Modern Serif"]})
home_dir = os.getcwd()
main_data_dir = home_dir+'/longwrite_meas'

#AFM data exported as txt from Gwyddion
dose_list =[]; clear_space = []; clear_space_err = []
scansize = 1.5 #micron along dispersion direction
aspect_rat = 8.0
default_coords = 0.0,0.0,1.0,0.0 #lev_x1,lev_y1,lev_x2,lev_y2 for flattening
data_split = 8; cutoff1,ind_cut = 40.0,20
nm_cut = 1.0; data_split = 4

#dose_list.append(44.)
#space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose44.0_00020.txt',skip_header=4),
#    main_data_dir+"/GEBL_nodelay_44",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=False,aspect_rat=8.0)
#clear_space.append(space); clear_space_err.append(err)

dose_list.append(45.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose45.0_00019.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_45",scansize,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(46.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose46.0_00018.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_46",scansize,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(47.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose47.0_00017.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_47",scansize,(0.2,0.0,1.5,1.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(48.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose48.0_00016.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_48",scansize,(0.015,0.0,1.26,0.9),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(49.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose49.0_00015.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_49",scansize,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(50.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose50.0_00014.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_50",scansize,(0.0,0.0,1.23,4.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(51.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose51.0_00013.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_51",scansize,(0.03,1.4,1.32,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(52.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/nodelaysample_dose52.0_00012.txt',skip_header=4),
    main_data_dir+"/GEBL_nodelay_52",scansize,(0.015,6.4,1.27,0.4),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space_array1.txt', [dose_list,clear_space,clear_space_err]); 
meas_space_nodelay = np.genfromtxt('meas_space_array1.txt')

dose_list =[]; clear_space = []; clear_space_err = []

dose_list.append(44.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose44.0_00010.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_44",scansize,(0.165,0.0,1.44,7.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(45.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose45.0_00000.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_45",scansize,(0.1,8.1,1.39,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(46.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose46.0_00001.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_46",scansize,(0.093,13.5,1.4,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(47.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose47.0_00002.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_47",scansize,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(48.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose48.0_00003.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_48",scansize,(0.04,14.1,1.354,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(49.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose49.0_00004.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_49",scansize,(0.2,0.0,1.495,13.8),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(50.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose50.0_00006.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_50",scansize,(0.15,0.0,1.468,8.8),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(51.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose51.0_00008.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_51",scansize,(0.013,12.3,1.355,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(52.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-07-26_GEBLarrays_redo/delaysample_dose52.0_00009.txt',skip_header=4),
    main_data_dir+"/GEBL_delay_52",scansize,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space_array2.txt', [dose_list,clear_space,clear_space_err]); 
meas_space_delay = np.genfromtxt('meas_space_array2.txt')


dose_list =[]; clear_space = []; clear_space_err = []
scansize = 2.0 #micron along dispersion direction
aspect_rat = 4.0

dose_list.append(40.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2022-12-21_GEBLref_wafer/dose40_grating.0_00007.txt',skip_header=4),
    main_data_dir+"/oldGEBL1_40",scansize,(0.088,5.3,1.996,0.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(41.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2022-12-21_GEBLref_wafer/dose41_grating.0_00004.txt',skip_header=4),
    main_data_dir+"/oldGEBL1_41",scansize,(0.075,4.3,1.987,0.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space_arrayold1.txt', [dose_list,clear_space,clear_space_err]); 
meas_space_old1 = np.genfromtxt('meas_space_arrayold1.txt')

dose_list =[]; clear_space = []; clear_space_err = []
scansize = 2.0 #micron along dispersion direction
aspect_rat = 4.0

#(0.1,8.0,1.99,0.0)

dose_list.append(44.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-01-09_GEBLref_wafer2/GEBL_dose44.0_00001.txt',skip_header=4),
    main_data_dir+"/oldGEBL2_44",scansize,(0.175,8.2,1.808,0.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)



dose_list.append(45.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-01-09_GEBLref_wafer2/GEBL_dose45.0_00002.txt',skip_header=4),
    main_data_dir+"/oldGEBL2_45",scansize,(0.1,8.0,1.99,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(46.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-01-09_GEBLref_wafer2/GEBL_dose46.0_00003.txt',skip_header=4),
    main_data_dir+"/oldGEBL2_46",scansize,(0.25,0.0,1.85,5.8),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(49.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-01-09_GEBLref_wafer2/GEBL_dose49.0_00007.txt',skip_header=4),
    main_data_dir+"/oldGEBL2_49",scansize,(0.245,0.0,1.853,3.7),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(50.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-01-09_GEBLref_wafer2/GEBL_dose50.0_00006.txt',skip_header=4),
    main_data_dir+"/oldGEBL2_50",scansize,(0.14,3.3,1.85,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

dose_list.append(51.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2023-01-09_GEBLref_wafer2/GEBL_dose51.0_00005.txt',skip_header=4),
    main_data_dir+"/oldGEBL2_51",scansize,(0.4,0.0,1.960,2.4),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space_arrayold2.txt', [dose_list,clear_space,clear_space_err]); 
meas_space_old2 = np.genfromtxt('meas_space_arrayold2.txt')



dose_list =[]; clear_space = []; clear_space_err = []
scansize = 1.5 #micron along dispersion direction
aspect_rat = 4.0

#(0.1,8.0,1.99,0.0)
#dose_list.append(33.)
#space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3300_40nA_lowdose.0_00008.txt',skip_header=4),
#    main_data_dir+"/nom3300_33",scansize,(0.060,10.4,1.36,0.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
#clear_space.append(space); clear_space_err.append(err)

#dose_list.append(35.)
#space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3300_40nA_nomdose.0_00003.txt',skip_header=4),
#    main_data_dir+"/nom3300_35",scansize,(0.116,2.5,1.4,0.2),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
#clear_space.append(space); clear_space_err.append(err)

dose_list.append(35.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3300_40nA_nomdose.0_00004.txt',skip_header=4),
    main_data_dir+"/nom3300_35",scansize,(0.225,0.0,1.216,8.6),data_split,nm_cut,cutoff1,ind_cut,close_fig=False,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)



dose_list.append(37.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3300_40nA_highdose.0_00010.txt',skip_header=4),
    main_data_dir+"/nom3300_37",scansize,(0.042,10.6,1.37,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)


np.savetxt('meas_space_beamtest1.txt', [dose_list,clear_space,clear_space_err]); 
meas_space_beamtest1 = np.genfromtxt('meas_space_beamtest1.txt')

dose_list =[]; clear_space = []; clear_space_err = []
scansize = 1.5 #micron along dispersion direction
aspect_rat = 4.0


#dose_list.append(33.)
#space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3400_40nA_nomdose.0_00006.txt',skip_header=4),
#    main_data_dir+"/nom3400_35",scansize,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=False,aspect_rat=aspect_rat)
#clear_space.append(space); clear_space_err.append(err)


#3400_40nA_nomdose.0_00006.txt'

dose_list.append(35.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3400_40nA_nomdose.0_00006.txt',skip_header=4),
    main_data_dir+"/nom3400_35",scansize,(0.26,0.4,1.232,5.9),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)


dose_list.append(37.)
space,err = space_meas(dose_list[-1],np.genfromtxt('2022-10-17_GEBLbeamtest/3400_40nA_highdose.0_00005.txt',skip_header=4),
    main_data_dir+"/nom3400_37",scansize,(0.076,9.6,1.384,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=aspect_rat)
clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space_beamtest2.txt', [dose_list,clear_space,clear_space_err]); 
meas_space_beamtest2 = np.genfromtxt('meas_space_beamtest2.txt')

def udata_dose(file_str,base_dose):
    data_import = np.genfromtxt(file_str,skip_header=0)
    dat_float = ufloat(np.mean(data_import),np.std(data_import))
    return [ufloat(base_dose,0),dat_float]

refdose37 = udata_dose('ref37.dat',37.0)
refdose38 = udata_dose('ref38.dat',38.0)
refdose39 = udata_dose('ref39.dat',39.0)
refdose40 = udata_dose('ref40.dat',40.0)
refdose41 = udata_dose('ref41.dat',41.0)
refdose42 = udata_dose('ref42.dat',42.0)
refdose43_alt = udata_dose('ref43_wafer1.dat',43.0)

refdose43 = udata_dose('ref43.dat',43.0)
refdose44 = udata_dose('ref44.dat',44.0)
refdose45 = udata_dose('ref45.dat',45.0)
refdose46 = udata_dose('ref46.dat',46.0)
refdose47 = udata_dose('ref47.dat',47.0)
refdose48 = udata_dose('ref48.dat',48.0)
refdose49 = udata_dose('ref49.dat',49.0)
refdose50 = udata_dose('ref50.dat',50.0)
refdose51 = udata_dose('ref51.dat',51.0)

refdose_arr1 = np.array([refdose37,refdose38,refdose39,refdose40,refdose41,refdose42,refdose43_alt]).T #,refdose39
refdose_arr2 = np.array([refdose43,refdose44,refdose45,refdose46,refdose47,refdose48,refdose49,refdose50,refdose51]).T



plt.ion(); plt.figure()
spin1 = 1
plt.errorbar(meas_space_nodelay[0],meas_space_nodelay[1]/spin1,
    yerr=meas_space_nodelay[2]/spin1,xerr=0*meas_space_nodelay[2],fmt='.-',label='no delay')
spin2 = 1
plt.errorbar(meas_space_delay[0],meas_space_delay[1]/spin2,
    yerr=meas_space_delay[2]/spin2,xerr=0*meas_space_delay[2],fmt='.-',label='delay')
spin3 = 1
#plt.errorbar(meas_space_old1[0],meas_space_old1[1]/spin3,
 #   yerr=meas_space_old1[2]/spin3,xerr=0*meas_space_old1[2],fmt='.-',label='old wafer 1')
spin4 = 1.
plt.errorbar(meas_space_old2[0],meas_space_old2[1]/spin4,
    yerr=meas_space_old2[2]/spin4,xerr=0*meas_space_old2[2],fmt='.-',label='old wafer 2')
plt.errorbar(meas_space_beamtest1[0],meas_space_beamtest1[1],
    yerr=meas_space_beamtest1[2],xerr=0*meas_space_beamtest1[2],fmt='.-',label='old beamtest 3300')
plt.errorbar(meas_space_beamtest2[0],meas_space_beamtest2[1],
    yerr=meas_space_beamtest2[2],xerr=0*meas_space_beamtest2[2],fmt='.-',label='old beamtest 3400')
plt.errorbar(unp.nominal_values(refdose_arr1)[0],unp.nominal_values(refdose_arr1)[1],
    yerr=unp.std_devs(refdose_arr1)[1],xerr=unp.std_devs(refdose_arr1)[0],fmt='.-',label='wafer 1')
plt.errorbar(unp.nominal_values(refdose_arr2)[0],unp.nominal_values(refdose_arr2)[1],
    yerr=unp.std_devs(refdose_arr2)[1],xerr=unp.std_devs(refdose_arr2)[0],fmt='.-',label='wafer 2')
plt.xlabel(r'base GEBL dose ($\mu$C/cm$^2$)'); plt.ylabel('cleared feaure width (nm)')
plt.xlim(32,55.); #plt.ylim(0,2)#
plt.ylim(0.,175.)
plt.axhline(y=315.151515/4,color='k',linestyle='--',alpha=0.5,label='nominal target width')
plt.axhline(y=120.,color='r',linestyle='--',alpha=0.5,label='approx space on long write')
plt.title('PMMA: 2700 rpm on 6" wafer, 1.5 min develop')
plt.legend()
plt.show()

os.chdir(home_dir)