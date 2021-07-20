# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 10:05:30 2020

@author: MZampini


Reconstruction of mrd passed whose name is passed as input

"""

def MRDtoNIfTI(ExpFolder, MRDname, *argv):
    import sys
    #print(sys.version)
    import os
#    import subprocess
    import nibabel as nb
#    import scipy.ndimage
    import numpy as np
    
#    subprocess.call('python -m pip install --upgrade pip')
#    subprocess.call('pip install nibabel')
    sys.path.insert(1, 'C:/Users/MZampini/OneDrive - MR Solutions/Desktop/Deliverable/CodePython/')
    
#    from sur2nifti.S2N.SUR2NIFTI_MAZ import SUR2NIFTI_MAZ as s2n_MAZ
    import getMRD.get_mrd_3d as MRD
    
    
    ''' data reconstruction '''
    import win32com.client
    import time
#    import shutil
#    import distutils.dir_util
    
    # pip install distutils-pytest
    
    reconfolder = ExpFolder
    reconApp = win32com.client.Dispatch("recon.Application")
#    reconApp.Visible = True
    
    reconfolder = reconfolder.replace('\\', '/')
    
    os.chdir(reconfolder)
    
    now_folder = reconfolder
    MRDfile_path = now_folder + MRDname
    #print(MRDfile_path)
    reconApp.DataFile = MRDfile_path
    reconApp.ImageFile = now_folder + 'temp.sur'
    
#    if argv:
    reconApp.FirstEcho = str(argv[0])
    reconApp.LastEcho = str(argv[0])
    #reconApp. = str(2) # how to use reconmax?
    
    reconApp.Run
    while not reconApp.StatusID==4:
        time.sleep(0.1) 
#    reconApp.Visible = False
#    print('reconstruction FINISHED')






                    
    if not ExpFolder.endswith('/'):
        ExpFolder =  ExpFolder + '/'
    
    img_dir = ExpFolder
    os.chdir(img_dir)
#    odirname = img_dir
                    
    i = os.listdir()[0]        
    dirname = img_dir 
    os.chdir(dirname)
    dirlist = os.listdir()
    
    sur_file = [name for name in dirlist if name.lower().endswith('.sur') ]

    
#    mrd_file = MRDname
    im, dim, par, scaling = MRD.get_mrd_3d(ExpFolder + sur_file[0])
    
    try:
        FOV = par['FOV']
        FOV = FOV.split(", ")
        FOV = float(FOV[-1])
        nx = par['NO_SAMPLES']
        nx = nx.split(", ")
        nx = int(nx[-1])
        dx = FOV/nx
        
        ny = nx
        dy = dx
        
        st = par['SLICE_THICKNESS']
        st = st.split(", ")
        nv2 = par['NO_VIEWS_2']
        nv2 = nv2.split(", ")
        nz = int(nv2[-1])
        dz = float(st[-1])/float(nz)
    except:
        dx = 0.1
        dy = 0.1
        dz = 0.1
        
    
#    im, dim, par, scaling = MRD.get_mrd_3d(dirname + '/' + mrd_file)

    
    #print(scalingSUR)
    
    imSUR, dimSUR, parSUR, scalingSUR = MRD.get_mrd_3d(dirname + '/' + sur_file[0])
    im_tot = np.zeros( np.squeeze(imSUR).shape + (len(sur_file),))
    
    for i in range(0,len(sur_file)):
        imSUR, dimSUR, parSUR, scalingSUR = MRD.get_mrd_3d(dirname + '/' + sur_file[i])
        # nb.viewers.OrthoSlicer3D(np.squeeze(np.abs(np.squeeze(imSUR))/scalingSUR))
        # import matplotlib.pyplot as plt
        # imgplot = plt.imshow(np.squeeze(np.abs(np.squeeze(imSUR))/scalingSUR))
        # plt.show()
        # imSUR_abs = np.squeeze(np.abs(np.squeeze(imSUR)))
        imSUR_scaled = np.squeeze(np.abs(np.squeeze(imSUR))/scalingSUR)
        im_tot[:,:,i] = imSUR_scaled

    time.sleep(1)
    for i in sur_file:
        os.remove(ExpFolder+i)
        
    affmat = [[dx, 0, 0, 0],[0, dy, 0, 0],[0, 0, dz, 0], [0, 0, 0, 1]]
    nb.save(nb.Nifti1Image(im_tot, affmat), ExpFolder +'/' + argv[1] +'.nii') 

    #print('Conversion to NIfTI terminated! :)')
    
    return im_tot, par
    
#MRDtoNIfTI('C:/smis/FLASH_test/test_phasecycling/', '000.MRD', 1)
