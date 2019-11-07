#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:02:19 2019

@author: carlarossell
"""

#Extract properties from watershed
import glob
import h5py
import numpy as np
from skimage.measure import label, regionprops
import cv2
import scipy

#Define functions to sort datasets
def custom_sorth5(name_fov): 
    [other,value]=name_fov.split('sub')
    [value,other,other]=value.split('_')    
    return value
def custom_sortdapi(name_dapi): 
    [other,value]=name_dapi.split('sub')
    [value,other]=value.split('.')    
    return value



#Read all h5 files in a directory. 
h5_files=glob.glob('*.h5')
h5_files.sort(key=custom_sorth5)
#print(h5_files)

#Read all dapi files
dapi_files=glob.glob('*.tif')
dapi_files.sort(key=custom_sortdapi)
#Get all the FOVs and arrange them in the position
#To get the range we have to add one
number_fovs=np.arange(1,len(h5_files)+1)
#fov_matrix=np.array(number_fovs).reshape(rows,columns)
fov_matrix=np.array(number_fovs).reshape(23,25)
dapi_fov_counter=0
whole_fov_sparses=[]
whole_fov_location=[]
fov_location=[]
for current_fov in h5_files:    
    coord_list=[]
    area_list=[]
    fov_list=[]
    #The h5 extension  is 'exported_data'
    #name_fov=h5_files[0]
    name_fov=current_fov
    fov = h5py.File(name_fov, 'r')
    mask=fov['/exported_watershed_masks'][:]
    mask_reduced=np.squeeze(mask, axis=2)

    #Get DAPI fov 
    dapi_fov= cv2.imread(dapi_files[dapi_fov_counter],cv2.IMREAD_GRAYSCALE)
    
    #Check which position the FOV occupies within the big scan
    #Position of FOV ilastik_masks_watershed1.h5
    [other,value]=name_fov.split('sub')
    [value,other,other]=value.split('_')
    (x,y)=np.where(fov_matrix==int(value))
    
    
   
    
    #Label the  mask
    mask_label=label(mask_reduced)
    #plt.imshow(mask_reduced)
    type(mask)
    [height,width]=mask_reduced.shape
    

    #Set FOV counter
    dapi_fov_counter=0
    #Create a list of the sparses for each FOV
    all_fov_sparses=[]
    z_dimension=0
    for i in range(0,np.amax(mask_label)+1):
        xmask,ymask=np.where(mask_label==i)
        single_cell_mask=np.zeros([1024,1024])
        single_cell_mask[xmask,ymask]=dapi_fov[xmask,ymask]
	#Transform to csr_matrix file
        sparse=scipy.sparse.csr_matrix(single_cell_mask)
        #Append the sparses on the list
        all_fov_sparses.append(sparse)
        fov_location.append((x,y))
        z_dimension= z_dimension+1

    #Save all sparse array in the fov
    whole_fov_sparses.append(all_fov_sparses)
    whole_fov_location.append(fov_location)
    
    #Move to the next dapi fov
    dapi_fov_counter=dapi_fov_counter+1
    
#all_properties={"Fov_Position":fov_list,"Centroids":coord_list,"Area":area_list}
          
np.savez('fov_data.npz',all_fov_sparses=all_fov_sparses)
