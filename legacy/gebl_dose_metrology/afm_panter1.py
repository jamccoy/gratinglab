import os
import numpy as np
import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pyplot as plt 
import array as arr
from scipy.signal import argrelextrema
from collections import defaultdict
from uncertainties import ufloat,unumpy as unp

plt.close("all"); plt.ion()
plt.rcParams.update({"text.usetex": True, "font.family": "sans-serif", "font.sans-serif": ["Computer Modern Serif"]})
cmap = mpl.cm.rainbow
home_dir = os.getcwd()
main_data_dir = home_dir+'/panter1_meas'

def raw_data(data,scansize):
    disp_um = np.arange(0,len(data[0]))*scansize/(len(data[0])-1)
    groov_nm = (10**9)*(np.mean(data,axis=0)-min(np.mean(data,axis=0)))
    return disp_um,groov_nm

def linear_flat_corr(disp,groov,lev_x1,lev_y1,lev_x2,lev_y2):
    flat_line = ((lev_y2-lev_y1)/(lev_x2-lev_x1))*(disp-lev_x1)+lev_y1
    return groov-flat_line

def find_maxima(disp,height,cutoff1,ind_cut):
    disp_um_cut,groov_nm_corr_cut = disp[np.where(height>cutoff1)],height[np.where(height>cutoff1)]
    ind1 = argrelextrema(groov_nm_corr_cut, np.greater)[0] #find relative maxima of this subset
    ind_cond = np.diff(ind1) - ind_cut #make sure found maxima aren't too close together 
    #print(ind_cond)
    if len(np.where(ind_cond<0)[0]) > 0:
        print(np.where(ind_cond<0)[0])
        del_ind = np.where(ind_cond<0)[0] + 1
        ind1 = np.delete(ind1,del_ind)
    max_disp,max_groov = disp_um_cut[ind1],groov_nm_corr_cut[ind1]; max_list=[]
    for i in range(0,len(ind1)):
        max_list.append(np.argmin(abs(disp-max_disp[i])))
    return max_disp,max_groov,max_list

def groov_meas(index,groov_dict,maxlist,nm_cut,scan_size,fig_dir):
    meas_list = []; plt.figure(); index = int(index)
    for i in range(0,len(groov_dict)):
        cutoff2 = nm_cut+min(groov_dict[i][0][1]) #1.5 nm above lowest point...
        botrange = np.where(groov_dict[i][0][1]<=cutoff2)[0]
        plt.plot(groov_dict[i][0][0],groov_dict[i][0][1],color=cmap(i/(len(maxlist)-1)))
        plt.axhline(y=cutoff2,color=cmap(i/(len(maxlist)-1)),alpha=0.5)
        plt.axvline(x=groov_dict[i][0][0][min(botrange)],linestyle='--',color=cmap(i/(len(maxlist)-1)),alpha=0.5)
        plt.axvline(x=groov_dict[i][0][0][max(botrange)],linestyle='--',color=cmap(i/(len(maxlist)-1)),alpha=0.5)
        space = (groov_dict[i][0][0][max(botrange)] - groov_dict[i][0][0][min(botrange)])*1000
        if space>40.:
            meas_list.append(space)
    plt.xlim(0,scan_size)
    plt.xlabel('dispersion direction ($\mu$m)'); plt.ylabel('groove depth (nm)'); plt.title(r'AFM of TASTE grating')
    os.chdir(fig_dir); plt.savefig('groovemeas'+str(index)+'.png'); os.chdir(home_dir); plt.close()
    return ufloat(np.round(np.mean(meas_list),2),np.round(np.std(meas_list),2))

def space_meas(pos,data,fig_dir,scansize,flatten,data_split):
    isExist = os.path.exists(fig_dir) # Check whether the specified path exists or not
    if not isExist:
        os.makedirs(fig_dir) # Create a new directory because it does not exist
    disp_um,groov_nm = raw_data(data,scansize); groov_nm_corr = linear_flat_corr(disp_um,groov_nm,flatten[0],flatten[1],flatten[2],flatten[3])
    plt.figure() #have a look at how the data are flattened
    plt.plot(disp_um,groov_nm,color='k',linestyle='--',alpha=0.5)
    plt.plot(disp_um,groov_nm_corr,color='k')
    plt.axhline(y=0,linestyle='--',color='k',alpha=0.2)
    plt.ylim(-5,120); plt.xlim(0,scansize)
    plt.xlabel('dispersion direction ($\mu$m)'); plt.ylabel('groove depth (nm)'); plt.title(r'full scan average at '+str(pos)+' mm')
    os.chdir(fig_dir); plt.savefig('fullscanavg'+str(pos)+'.png'); os.chdir(home_dir)
    #if close_fig==True:
    #    plt.close()
    #split up the image
    data_dict = defaultdict(list); groov_nm_dict = defaultdict(list); meas_space = []; max_disp_list = []
    for i in range (0,data_split):
        entry = data[int(i*len(data)/int(data_split)):int((i+1)*len(data)/int(data_split))]
        data_dict[i].append(entry)
        groov_nm = (10**9)*(np.mean(data_dict[i][0],axis=0)-min(np.mean(data_dict[i][0],axis=0)))
        groov_nm_dict[i].append(linear_flat_corr(disp_um,groov_nm,flatten[0],flatten[1],flatten[2],flatten[3])) 
    for k in range(0,len(groov_nm_dict)):
        groov_nm_corr = groov_nm_dict[k][0]; plt.figure()
        max_disp,max_groov,max_ind_list = find_maxima(disp_um,groov_nm_corr,40.0,20); single_groov_dict = defaultdict(list)
        for i in range(0,len(max_ind_list)-1):
            single_groov_dict[i].append([disp_um[max_ind_list[i]:max_ind_list[i+1]],groov_nm_corr[max_ind_list[i]:max_ind_list[i+1]]])
        plt.plot(disp_um,groov_nm_corr,color='k',alpha=0.5)
        max_disp_list.append(max_disp)
        plt.plot(max_disp,max_groov,'ro')
        for i in range(0,len(max_ind_list)-1):
            plt.plot(single_groov_dict[i][0][0],single_groov_dict[i][0][1],color=cmap(i/(len(max_ind_list)-1)),linestyle='--')
        plt.xlim(0,scansize)
        meas_space.append(groov_meas(k+1,single_groov_dict,max_ind_list,1.0,scansize,fig_dir))
        plt.xlabel('dispersion direction ($\mu$m)'); plt.ylabel('groove depth (nm)')
        plt.title(r'single grooves from partial average '+str(k+1)+' at '+str(pos)+' mm')
        os.chdir(fig_dir); plt.savefig('scanavg'+str(k+1)+'.png'); os.chdir(home_dir); plt.close()
    delta_x = abs(max_disp_list[1] - max_disp_list[int(data_split-1)]); 
    term = (2*delta_x)/(3*scansize)
    cos_term = (1-(term**2))/((term**2)+1); rot = ufloat(np.mean(cos_term),np.std(cos_term)) 
    avg_space = np.mean(meas_space)/rot
    return np.round(avg_space.nominal_value,2),np.round(avg_space.std_dev,2)

clear_space = []; clear_space_err = []; pos_from_bot = []
scansize = 2.0 #micron along dispersion direction
default_coords = 0.0,0.0,1.0,0.0 #lev_x1,lev_y1,lev_x2,lev_y2 for flattening
data_split = 7

#close_fig =[False,False,False,False,False,False,False]
#AFM data exported as txt from Gwyddion
pos_list =[]; data_list = []; fig_dir_list = []; flatten_list = []; close_fig =[]


pos_list.append(65.)
data_list.append(np.genfromtxt('2023-09-21_longwrite/5mm_from_top.0_00000.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_5mm_from_top")
flatten_list.append((0.019,4.7,1.724,0.3))
close_fig.append(True)

#pos_list.append(60.)
#data_list.append(np.genfromtxt('2023-09-21_longwrite/10mm_from_top.0_00001.txt',skip_header=4,skip_footer=100)) #skip_header=104
#fig_dir_list.append(main_data_dir+"/center_10mm_from_top")
#flatten_list.append((0.0,0.0,1.0,0.0))
#close_fig.append(True)

pos_list.append(55.)
data_list.append(np.genfromtxt('2023-09-21_longwrite/15mm_from_top.0_00002.txt',skip_header=4))
fig_dir_list.append(main_data_dir+"/center_15mm_from_top")
flatten_list.append((0.0,0.0,1.0,0.0))
close_fig.append(True)

pos_list.append(50.)
data_list.append(np.genfromtxt('2023-09-21_longwrite/20mm_from_top.0_00003.txt',skip_header=64))
fig_dir_list.append(main_data_dir+"/center_20mm_from_top")
flatten_list.append((0.028,10.7,1.742,0.6))
close_fig.append(True)

pos_list.append(45.)
data_list.append(np.genfromtxt('2023-09-21_longwrite/25mm_from_top.0_00004.txt',skip_header=64))
fig_dir_list.append(main_data_dir+"/center_25mm_from_top")
flatten_list.append((0.006,3.6,1.690,0.5))
close_fig.append(True)



for i in range(0,len(pos_list)):
    space,err = space_meas(pos_list[i],data_list[i],fig_dir_list[i],scansize,flatten_list[i],data_split)
    clear_space.append(space); clear_space_err.append(err)

np.savetxt('meas_space_panter1.txt', [pos_list,clear_space,clear_space_err]) 

meas_space = np.genfromtxt('meas_space_panter1.txt')

plt.ion()
plt.figure()
plt.errorbar(meas_space[0],meas_space[1],
    yerr=meas_space[2],xerr=0*meas_space[2],fmt='.-',label='long write')
plt.ylabel('cleared feaure width (nm)')
plt.xlabel('central groove direction (mm)')
plt.xlim(0,70.)
plt.ylim(0.,175.)
plt.axhline(y=315.151515/4,color='k',linestyle='--',alpha=0.5,label='nominal target width')
plt.title(r'29 $\mu$C/cm$^2$, 2.0 min develop time')
plt.legend()
plt.show()

