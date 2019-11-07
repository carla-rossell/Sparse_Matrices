#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:02:19 2019

@author: carlarossell
"""

#Extract properties from watershed
import glob
import h5py
from skimage.measure import label, regionprops
from matplotlib import pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
import cv2


#Sort .h5 and .tif datasets
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

#Read all tif files in directory
dapi_files=glob.glob('*.tif')
dapi_files.sort(key=custom_sortdapi)
#Get all the FOVs and arrange them in the position
#To get the range we have to add one
number_fovs=np.arange(1,len(h5_files)+1)
#fov_matrix=np.array(number_fovs).reshape(rows,columns)
fov_matrix=np.array(number_fovs).reshape(23,25)

#Set counter to check current fov
dapi_fov_counter=0
#Create a list in which to save the sparses
all_fov_sparses=[]

for current_fov in h5_files:    
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
    

    #Create a 3D sparse array where x,y are FOV size and z is the amount of nuclei in the FOV 
    #sparse_in_fov=np.zeros([1024,1024,np.amax(mask_label)])
    sparse_in_fov=csr((1024,1024,np.amax(mask_label),dtype='uint8')	
    #Set a counter to check current stack
    z_dimension=0
	# 0 is background so it doesn't get included in the range
    for i in range(1,np.amax(mask_label)+1):
        xmask,ymask=np.where(mask_label==i)
	single_cell_mask=csr((1024,1024),dtype='uint8')	        
	single_cell_mask[xmask,ymask]=dapi_fov[xmask,ymask]
	#Add current nuclei sparse on to the FOV array
        sparse_in_fov= np.insert(sparse_in_fov, z_dimension, single_cell_mask, axis=2)
	#Move to the next stack (next nuclei label)
        z_dimension= z_dimension+1

    #Regionprops. Extract properties from  each nuclei
    #Save all sparse array in the fov
    
    all_fov_sparses.append(sparse_in_fov)
    
    
    #Move to the next dapi fov
    dapi_fov_counter=dapi_fov_counter+1
    
#Save properties
          
np.savez('fov_data.npz',all_fov_sparses=all_fov_sparses)
