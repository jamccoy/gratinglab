import os
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 
import array as arr
from scipy.signal import argrelextrema
from collections import defaultdict
from uncertainties import ufloat,unumpy as unp
import scipy as sp
from scipy.odr import Model, RealData, ODR

plt.close("all"); plt.ion()
plt.rcParams.update({"text.usetex": True, "font.family": "sans-serif", "font.sans-serif": ["Computer Modern Serif"]})
cmap = mpl.cm.rainbow
home_dir = os.getcwd()

def line_fit(C,x):
    return C[0]*x+C[1]

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
        del_ind = np.where(ind_cond<0)[0] + 1
        ind1 = np.delete(ind1,del_ind)
    max_disp,max_groov = np.array(disp_um_cut[ind1]),np.array(groov_nm_corr_cut[ind1]); max_list=[]
    #max_disp_err = np.mean(abs(max_disp-max_disp_offp))
    max_disp_err = 0.5*(np.array(disp_um_cut[ind1+1])-np.array(disp_um_cut[ind1-1]))
    #max_disp_offp = np.array(disp_um_cut[ind1+1])
    for i in range(0,len(ind1)):
        max_list.append(np.argmin(abs(disp-max_disp[i])))
    #try adding uncertainty of one index in mx pos
    #print(max_list)
    return max_disp,max_groov,max_list,max_disp_err

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

def space_meas(pos,data,fig_dir,scansize,flatten,data_split,nm_cut,cutoff1,ind_cut,close_fig=True,aspect_rat=4.0,rot_ang=90.):
    isExist = os.path.exists(fig_dir) # Check whether the specified path exists or not
    if not isExist:
        os.makedirs(fig_dir) # Create a new directory because it does not exist
    disp_um,groov_nm = raw_data(data,scansize); groov_nm_corr = linear_flat_corr(disp_um,groov_nm,flatten[0],flatten[1],flatten[2],flatten[3])
    plt.figure() #have a look at how the data are flattened
    plt.plot(disp_um,groov_nm,color='k',linestyle='--',alpha=0.5)
    plt.plot(disp_um,groov_nm_corr,color='k')
    plt.axhline(y=0,linestyle='--',color='k',alpha=0.2)
    plt.ylim(-5,140); plt.xlim(0,scansize)
    plt.xlabel('dispersion direction ($\mu$m)'); plt.ylabel('groove depth (nm)'); plt.title(r'full scan average at '+str(pos)+' mm')
    os.chdir(fig_dir); plt.savefig('fullscanavg'+str(pos)+'.png'); os.chdir(home_dir)
    if close_fig==True:
        plt.close()
    #split up the image
    data_dict = defaultdict(list); groov_nm_dict = defaultdict(list); meas_space = []; max_disp_list = []
    for i in range (0,data_split):
        entry = data[int(i*len(data)/int(data_split)):int((i+1)*len(data)/int(data_split))]
        data_dict[i].append(entry)
        groov_nm = (10**9)*(np.mean(data_dict[i][0],axis=0)-min(np.mean(data_dict[i][0],axis=0)))
        groov_nm_dict[i].append(linear_flat_corr(disp_um,groov_nm,flatten[0],flatten[1],flatten[2],flatten[3])) 
    for k in range(0,len(groov_nm_dict)):
        groov_nm_corr = groov_nm_dict[k][0]; plt.figure()
        max_disp,max_groov,max_ind_list,max_disp_err = find_maxima(disp_um,groov_nm_corr,cutoff1,ind_cut); single_groov_dict = defaultdict(list)
        for i in range(0,len(max_ind_list)-1):
            single_groov_dict[i].append([disp_um[max_ind_list[i]:max_ind_list[i+1]],groov_nm_corr[max_ind_list[i]:max_ind_list[i+1]]])
        plt.plot(disp_um,groov_nm_corr,color='k',alpha=0.5)
        #print(max_disp)
        max_disp_list.append(max_disp)
        plt.plot(max_disp,max_groov,'ro')
        for i in range(0,len(max_ind_list)-1):
            plt.plot(single_groov_dict[i][0][0],single_groov_dict[i][0][1],color=cmap(i/(len(max_ind_list)-1)),linestyle='--')
        plt.xlim(0,scansize)
        meas_space.append(groov_meas(k+1,single_groov_dict,max_ind_list,nm_cut,scansize,fig_dir))
        plt.xlabel('dispersion direction ($\mu$m)'); plt.ylabel('groove depth (nm)')
        plt.title(r'single grooves from partial average '+str(k+1)+' at '+str(pos)+' mm')
        os.chdir(fig_dir); plt.savefig('scanavg'+str(k+1)+'.png'); os.chdir(home_dir); plt.close() #now, correct for rotation
    fit_arr_x = np.array([max_disp_list],dtype=object).T[0] #to be indexed
    fit_arr_y = np.linspace(0,scansize/aspect_rat,len(max_disp_list)); pix_err = (3.9/1000.)*np.ones(len(fit_arr_y))#not to be indexed
    deg_guess = 89.; slope_guess = np.tan(deg_guess*np.pi/180.); yint_guess = -15.; rot_test_list = []
    #plt.figure(); #print(slope_guess)
    for i in range(0,len(fit_arr_x)):
        #plt.plot(fit_arr_x[i],fit_arr_y)
        #print(len(fit_arr_x[i])), print(len(fit_arr_y))
        if len(fit_arr_x[i])==len(fit_arr_y):
            fit_model = sp.odr.Model(line_fit); datafit = sp.odr.RealData(fit_arr_x[i],fit_arr_y,sx=pix_err,sy=pix_err)
            odr = sp.odr.ODR(datafit,fit_model,beta0=[slope_guess,yint_guess]); out = odr.run()
            slope,yint = ufloat(out.beta[0],out.sd_beta[0]),ufloat(out.beta[1],out.sd_beta[1])
            #print(slope)
            if abs(slope)>40.:
            #deg_fit = unp.arctan(slope)*(180./np.pi)
                rot_test_list.append(abs(unp.sin(unp.arctan(slope))))
            else:
                rot_test_list.append(ufloat(1.0,0.0))
    #delta_x = abs(max_disp_list[0] - max_disp_list[int(data_split-1)]); 
    #term = (2*delta_x)/(3*scansize)
    #cos_term = (1-(term**2))/((term**2)+1); rot = ufloat(np.mean(cos_term),np.std(cos_term)) 
    if rot_test_list and rot_ang!=90.:
        rot_test = np.mean(rot_test_list)
    else:
        rot_test = ufloat(np.sin(rot_ang*np.pi/180.),0.05)
    #print(rot_test)
    avg_space = np.mean(meas_space)*rot_test
    return np.round(avg_space.nominal_value,2),np.round(avg_space.std_dev,2)
