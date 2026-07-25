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
pos_list =[]; data_list = []; fig_dir_list = []; flatten_list = []; scansize_list =[]; aspect_list =[] 
clear_space = []; clear_space_err = []; pos_from_bot = []; close_fig_list =[]
scansize = 2.0 #micron along dispersion direction
aspect_rat = 4.0
default_coords = 0.0,0.0,1.0,0.0 #lev_x1,lev_y1,lev_x2,lev_y2 for flattening
data_split = 8; cutoff1,ind_cut = 40.0,20
nm_cut = 1.5

pos_list.append(5.)
data_list.append(np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_5mmfrombot.0_00004.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_5mm_from_bot")
flatten_list.append(default_coords); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)
#space,err = space_meas(pos_list[i],data_list[i],fig_dir_list[i],scansize_list[i],flatten_list[i],data_split,nm_cut,cutoff1,ind_cut,close_fig=close_fig_list[i],aspect_rat=4.0)

pos_list.append(10.)
data_list.append(np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_10mmfrombot.0_00005.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_10mm_from_bot")
flatten_list.append((0.263,0.04,1.847,2.3)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)
'''
pos_list.append(15.)
data_list.append(np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_15mmfrombot.0_00007.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_15mm_from_bot")
flatten_list.append(default_coords); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(20.)
data_list.append(np.genfromtxt('2022-12-06_GEBL_lonmgwrite_cont/longwrite_20mmfrombot.0_00008.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_20mm_from_bot")
flatten_list.append((0.17,1.3,1.75,0.05)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(25.)
data_list.append(np.genfromtxt('2023-01-24-morelongwrite/center_25mm_frombot.0_00003.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_25mm_from_bot")
flatten_list.append((0.05,0.7,1.34,0.2)); scansize_list.append(1.5); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(30.)
data_list.append(np.genfromtxt('2023-01-24-morelongwrite/center_30mm_frombot.0_00004.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_30mm_from_bot")
flatten_list.append((0.04,0.05,1.34,1.66)); scansize_list.append(1.5); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(35.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_35mm_from_top.0_00007.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_35mm_from_top")
flatten_list.append(default_coords); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(40.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_30mm_from_top.0_00006.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_30mm_from_bot")
flatten_list.append((0.228,0.3,1.84,1.1)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(45.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_25mm_from_top.0_00005.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_25mm_from_top")
flatten_list.append(default_coords); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(50.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_20mm_from_top.0_00004.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_20mm_from_top")
flatten_list.append((0.0315,0.5,1.9740,0.8)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(55.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_15mm_from_top.0_00003.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_15mm_from_top")
flatten_list.append((0.175,1.35,1.78,0.0)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(60.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_10mm_from_top.0_00002.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_10mm_from_top")
flatten_list.append((0.015,3.37,1.95,0.08)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

pos_list.append(65.)
data_list.append(np.genfromtxt('2023-01-10_morelongwrite/center_5mm_from_top.0_00000.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_5mm_from_top")
flatten_list.append((0.1,8.2,2.0,0.5)); scansize_list.append(2.0); aspect_list.append(aspect_rat); close_fig_list.append(True)

for i in range(0,len(pos_list)):
    space,err = space_meas(pos_list[i],data_list[i],fig_dir_list[i],scansize_list[i],flatten_list[i],data_split,nm_cut,cutoff1,ind_cut,close_fig=close_fig_list[i],aspect_rat=4.0)
    clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space.txt', [pos_list,clear_space,clear_space_err]); meas_space = np.genfromtxt('meas_space.txt')

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
'''
'''
test = np.array(test_).T
aspect = 4.
scan = 2.0
#test2 = np.linspace(0,scan/aspect,len(test[0]))

i=0
plt.figure()
plt.errorbar(test[i],np.linspace(0,scan/aspect,len(test[i])),xerr=test_err[i],yerr=test_err[i],label='stylus contrast',color='g',linestyle='None')
plt.xlim(0,2)
plt.show()

fit_model = sp.odr.Model(line_fit); datafit = sp.odr.RealData(test[i],np.linspace(0,scan/aspect,len(test[i])),sx=test_err[i],sy=test_err[i])

deg_guess = 89.; slope_guess = np.tan(deg_guess*np.pi/180.); yint_guess = -15.

odr = sp.odr.ODR(datafit,fit_model,beta0=[slope_guess,yint_guess]); out = odr.run()
slope,yint = ufloat(out.beta[0],out.sd_beta[0]),ufloat(out.beta[1],out.sd_beta[1])
deg_fit = unp.arctan(slope)*(180./np.pi)
rot_test = unp.sin(unp.arctan(slope))
#print(rot_test)


plt.figure()
for i in range(0,len(test)):
	#plt.plot(test[i],np.linspace(0,scan/aspect,len(test[i])))
	plt.errorbar(test[i],np.linspace(0,scan/aspect,len(test[i])),xerr=test_err[i],yerr=test_err[i],label='stylus contrast',color='g',linestyle='None')
plt.show()
'''
