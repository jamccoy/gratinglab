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
pos_list =[]; clear_space = []; clear_space_err = []
scansize = 2.0 #micron along dispersion direction
aspect_rat = 4.0
default_coords = 0.0,0.0,1.0,0.0 #lev_x1,lev_y1,lev_x2,lev_y2 for flattening
data_split = 8; cutoff1,ind_cut = 40.0,20
nm_cut = 1.5; data_split = 4


pos_list.append(0.5)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_500um_from_bot.0_00000.txt',skip_header=4),
    main_data_dir+"/center_0p5mm_from_bot",2.0,(0.0175,0.0,1.9,3.5),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(1.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_1mm_from_bot.0_00001.txt',skip_header=4),
    main_data_dir+"/center_1mm_from_bot",2.0,(0.01,0.1,1.9,4.7),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(2.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_2mm_from_bot.0_00002.txt',skip_header=4),
    main_data_dir+"/center_2mm_from_bot",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(3.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_3mm_from_bot.0_00003.txt',skip_header=4),
    main_data_dir+"/center_3mm_from_bot",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(4.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_4mm_from_bot.0_00004.txt',skip_header=4),
    main_data_dir+"/center_4mm_from_bot",2.0,(0.18,0.1,1.77,1.5),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(5.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_5mmfrombot.0_00004.txt',skip_header=4),
    main_data_dir+"/center_4mm_from_bot",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(10.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_10mmfrombot.0_00005.txt',skip_header=4),
    main_data_dir+"/center_10mm_from_bot",2.0,(0.263,0.04,1.847,2.3),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(15.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_15mmfrombot.0_00007.txt',skip_header=4),
    main_data_dir+"/center_15mm_from_bot",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(20.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_20mmfrombot.0_00008.txt',skip_header=4),
    main_data_dir+"/center_20mm_from_bot",2.0,(0.17,1.3,1.75,0.05),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(25.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2023-01-24-morelongwrite/center_25mm_frombot.0_00003.txt',skip_header=4),
    main_data_dir+"/center_25mm_from_bot",1.5,(0.17,1.3,1.75,0.05),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(30.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2023-01-24-morelongwrite/center_30mm_frombot.0_00004.txt',skip_header=4),
    main_data_dir+"/center_30mm_from_bot",1.5,(0.04,0.05,1.34,1.66),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(35.)
space1,err1 = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_35mm_from_top.0_00007.txt',skip_header=4),
    main_data_dir+"/center_35mm_from_bot/file1",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
space_err1 = ufloat(space1,err1)
space2,err2 = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_verycenter.0_00006.txt',skip_header=4),
    main_data_dir+"/center_35mm_from_bot/file2",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=False,aspect_rat=4.0)
space_err2 = ufloat(space2,err2)
space_err = np.mean([space_err1,space_err2])
clear_space.append(space_err.nominal_value); clear_space_err.append(space_err.std_dev)

pos_list.append(40.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_30mm_from_top.0_00006.txt',skip_header=4),
    main_data_dir+"/center_30mm_from_top",2.0,(0.228,0.3,1.84,1.1),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(45.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_25mm_from_top.0_00005.txt',skip_header=4),
    main_data_dir+"/center_25mm_from_top",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(50.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_20mm_from_top.0_00004.txt',skip_header=4),
    main_data_dir+"/center_20mm_from_top",2.0,(0.0315,0.5,1.9740,0.8),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(55.)
space,err = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_15mm_from_top.0_00003.txt',skip_header=4),
    main_data_dir+"/center_15mm_from_top",2.0,(0.175,1.35,1.78,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(60.)
space1,err1 = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_10mm_from_top.0_00002.txt',skip_header=4),
    main_data_dir+"/center_10mm_from_top/file1",2.0,(0.015,3.37,1.95,0.08),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
space_err1 = ufloat(space1,err1)
space2,err2 = space_meas(pos_list[-1],np.genfromtxt('2022-12-19_morelongwrite/longwrite_center_10mm_from+top.0_00011.txt',skip_header=4),
    main_data_dir+"/center_10mm_from_top/file2",2.0,(0.013,1.3,1.926,0.0),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
space_err2 = ufloat(space2,err2)
space3,err3 = space_meas(pos_list[-1],np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_10mmfromtop.0_00001.txt',skip_header=4),
    main_data_dir+"/center_10mm_from_top/file3",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=False,aspect_rat=4.0)
space_err3 = ufloat(space3,err3)
space_err = np.mean([space_err1,space_err2,space_err3])
clear_space.append(space_err.nominal_value); clear_space_err.append(space_err.std_dev)



pos_list.append(65.)
space1,err1 = space_meas(pos_list[-1],np.genfromtxt('2023-01-10_morelongwrite/center_5mm_from_top.0_00000.txt',skip_header=4),
    main_data_dir+"/center_5mm_from_top/file1",2.0,(0.1,8.2,2.0,0.5),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
space_err1 = ufloat(space1,err1)
space2,err2 = space_meas(pos_list[-1],np.genfromtxt('2022-12-19_morelongwrite/longwrite_center_5mm_from_top.0_00004.txt',skip_header=4),
    main_data_dir+"/center_5mm_from_top/file2",2.0,(0.188,0.0,1.811,0.9),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
space_err2 = ufloat(space2,err2)
space_err = np.mean([space_err1,space_err2])
clear_space.append(space_err.nominal_value); clear_space_err.append(space_err.std_dev)




pos_list.append(66.0)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-19_morelongwrite/longwrite_center_4mm_from_top.0_00003.txt',skip_header=4),
    main_data_dir+"/center_4mm_from_top",2.0,(0.25,1.7,1.86,0.6),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(67.0)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-19_morelongwrite/longwrite_center_3mm_from_top.0_00002.txt',skip_header=4),
    main_data_dir+"/center_3mm_from_top",2.0,(0.178,2.2,1.8,0.8),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(68.0)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-19_morelongwrite/longwrite_center_2mm_from_top.0_00001.txt',skip_header=4),
    main_data_dir+"/center_2mm_from_top",2.0,(0.25,1.7,1.97,0.5),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(69.0)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_1mm_from_top.0_00008.txt',skip_header=4),
    main_data_dir+"/center_1mm_from_top",2.0,default_coords,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)

pos_list.append(69.5)
space,err = space_meas(pos_list[-1],np.genfromtxt('2022-12-14_morelongwrite/longwrite_center_500um_from_top.0_00007.txt',skip_header=4),
    main_data_dir+"/center_0p5mm_from_top",2.0,(0.275,0.0,1.9,2.7),data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0)
clear_space.append(space); clear_space_err.append(err)


np.savetxt('meas_space_longwrite.txt', [pos_list,clear_space,clear_space_err]); meas_space = np.genfromtxt('meas_space_longwrite.txt')

plt.ion(); plt.figure()
plt.errorbar(meas_space[0],meas_space[1],
    yerr=meas_space[2],xerr=0*meas_space[2],fmt='.-',label='long write')
plt.xlabel('central groove direction (mm)'); plt.ylabel('cleared feaure width (nm)')
plt.xlim(0,70.); plt.ylim(0.,175.)
plt.axhline(y=315.151515/4,color='k',linestyle='--',alpha=0.5,label='nominal target width')
plt.title(r'35 $\mu$C/cm$^2$, 1.5 min develop time')
plt.legend()
plt.show()

os.chdir(home_dir)

