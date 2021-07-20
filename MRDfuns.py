# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 11:28:45 2021

@author: MZampini

MRD file and RAMSES image reconstruction
"""
import numpy as np, sys, os, nibabel as nb
from numpy import pi as pi

def img2phase(dataMRI):
    '''
    This function returns the phase image img_phase provided real and imaginary part are provided in dataMRI.
    Unwarping is still needed afterwards
    '''
    img_phase = np.empty(dataMRI.shape) + 1j*np.empty(dataMRI.shape)
    img1a = np.real(dataMRI)*np.sign(np.imag(dataMRI))+1j*np.imag(dataMRI)*np.sign(np.real(dataMRI)) # seems to be working - needs unwrapping
    img1a_phase = np.angle(img1a)*np.sign(np.real(dataMRI)*np.sign(np.imag(dataMRI))) # seems to be working - needs unwrapping
    img_phase = (img1a_phase-(0))/(pi-(0))*(pi-(-pi))+(-pi) 
    return img_phase





def reg_circshift(imgFIX, imgMOV, maxs, metric):
    '''
    maxs = maximum shifts in X and Y directions to check for registration
    '''
    imgMOVonFIX = np.zeros(np.hstack([imgMOV.shape, 2*maxs+1, 2*maxs+1]))
    res = np.zeros(( 2*maxs+1, 2*maxs+1 ))
    xx = 0
    for ii in range(-maxs,maxs+1):
        yy = 0
        for jj in range(-maxs,maxs+1):
            imgMOVonFIX[:,:,:,xx,yy] = np.roll(imgMOV,jj, axis=0) # axis 0 = y, axis 1 = x
            if ii !=0:
                imgMOVonFIX[:,:,:,xx,yy] = np.roll(imgMOVonFIX[:,:,:,xx,yy],ii, axis=1) # axis 0 = y, axis 1 = x
            if metric == 'MSE':
                res[xx,yy] = (abs(imgMOVonFIX[:,:,:,xx,yy]-imgFIX)).sum()/ np.dot(imgFIX.shape,imgFIX.shape)
            yy += 1
        xx += 1
    return imgMOVonFIX, res




def kspace2data2D(kdata):
    if len(kdata.shape)>3:
        n3 = kdata.shape[3]
    else:
        n3 = 1
    img_K= np.zeros(np.hstack([kdata.shape[0:3],n3])) +1j*np.zeros(np.hstack([kdata.shape[0:3],n3]))
    for kk in range(n3):
        for ii in range(kdata.shape[2]):
            img_K[:,:,ii,kk] = np.fft.ifftshift(np.fft.ifftn(kdata[:,:,ii,kk], s = [kdata.shape[0],kdata.shape[1]]))
        if np.mod(kk,2):
            img_K[:,:,:,kk] = np.rot90(np.flip(np.squeeze(img_K[:,:,:,kk])), 2, axes=(0,2))
        #else:
        #    img_K[:,:,:,kk] = img_K[:,:,:,kk]# np.rot90(np.flip(np.squeeze(img_K[:,:,:,kk])), 2, axes=(0,2))
            
        if kk == 0: # if RAMSES_14            
            img_K[:,:,:,0] = np.flip(np.flip(np.squeeze(img_K[:,:,:,0])),2)
    return img_K



def kspace2data2D_LR(kdata,Nr):
    if len(kdata.shape)>3:
        n3 = kdata.shape[3]
    else:
        n3 = 1
    img_K= np.zeros(np.hstack([kdata.shape[0:3],n3])) +1j*np.zeros(np.hstack([kdata.shape[0:3],n3]))
    for kk in range(n3):
        for ii in range(kdata.shape[2]):
            kdatanow = np.squeeze( np.pad(kdata[:,int(n0/2)-Nr:int(n0/2)+Nr,ii,kk], ((0,0),(int((int(n0/2)-Nr)),int((int(n0/2)-Nr)))), 'constant', constant_values=(0)) )
            img_K[:,:,ii,kk] = np.fft.ifftshift(np.fft.ifftn(kdatanow, s = [kdata.shape[0],kdata.shape[1]]))
        
        if np.mod(kk,2):
            img_K[:,:,:,kk] = np.rot90(np.flip(np.squeeze(img_K[:,:,:,kk])), 2, axes=(0,2))
        if kk == 0: # if RAMSES_14            
            img_K[:,:,:,0] = np.flip(np.flip(np.squeeze(img_K[:,:,:,0])),2)
    return img_K





def MRDparameters(parMRD,matsize): # for RAMSES # TODO: MGE
    ''' parameters retrieval '''
    FOV_exp = int(parMRD['FOV'].split()[0])
    FOVf = int(parMRD['VAR FOVf'].split()[0])/8
    slice_thickness = float(parMRD['SLICE_THICKNESS'].split(', ')[2])
    no_echoes = int(parMRD['NO_ECHOES'].split()[1])-1 # ramses
    no_samples = int(parMRD['NO_SAMPLES'].split()[1])
    tacq =  1/ (int(parMRD['SAMPLE_PERIOD sample_period,'].split(',')[0]) * 1e3)
    tacq = tacq*no_samples
    tsel90 = 1332*1e-6 # for rfnum = 1
    tramp = int(parMRD['tramp,'].split()[0])*1e-6
    tdp =  int(parMRD['tdp,'].split()[0])*1e-6                      # not in MGE
    tref = tdp+tramp                                                # not in MGE
    TEmin = (tsel90/2) + 3*tramp + tref + (tacq/2)
    TE =  int(parMRD['te,'].split()[0])/1000
    if TE == 0:
        TE = TEmin
    RO_gr = tacq + 2*tramp
    TEarray = TE + RO_gr*np.array(range(no_echoes))
    n0,n1,n2,no_echoes = matsize
    acqsize = [FOV_exp/n0, FOV_exp*FOVf/n1, slice_thickness]
    
    return acqsize, TEarray




def RAMSESregistrationcomplex(path_RAMSES):
    
    
    os.chdir(path_RAMSES)
    dirlist = os.listdir()
    
    import getMRD.get_mrd_3d as MRD
    mrd_file = [name for name in dirlist if name.lower().endswith('.mrd') ]

    
    ''' Retrieving k-space data from MASTER_FOLDER + magnitude image reconstruction from k-space'''
    imgMRD, dimMRD, parMRD, _ = MRD.get_mrd_3d(path_RAMSES + mrd_file[0])
    n0 = int(parMRD['NO_SAMPLES'].split(", ")[-1])
    n1 = int(parMRD['NO_VIEWS'].split(", ")[-1])
    n2 = max([int(parMRD['NO_SLICES'].split(", ")[-1]), int(parMRD['NO_VIEWS_2'].split(", ")[-1]), ])
    no_echoes = int(parMRD['NO_ECHOES'].split(", ")[-1])
    sequencename = parMRD['PPL C'].split(os.sep)[-1]
    matsize = [n0,n1,n2,no_echoes]
    kdata0 = np.squeeze(imgMRD)
    kdata = np.zeros((n0,n1,n2,no_echoes)) + 1j*np.zeros((n0,n1,n2,no_echoes))
    img   = np.zeros((n0,n1,n2,no_echoes)) + 1j*np.zeros((n0,n1,n2,no_echoes))
    
    if n2>1:
        for jj in range(no_echoes):
            for ii in range(n2):
                if int(parMRD['NO_SLICES'].split(", ")[-1]) > 1:    # 2D multislice case
                    kdata[:,:,ii,jj] = kdata0[jj,ii,:,:] 
                    img[:,:,ii,jj]      = np.fft.ifftshift(np.fft.ifftn(np.squeeze(kdata[:,:,ii,jj])))
                else:                                               # 3D case
                    kdata[:,:,ii,jj] = kdata0[jj,:,ii,:] 
            if int(parMRD['NO_SLICES'].split(", ")[-1]) == 1:       # 3D case
                img[:,:,:,jj] = np.fft.ifftshift(np.fft.ifftshift(np.fft.ifftn(np.squeeze(kdata[:,:,:,jj]))), axes=2)
    else:
        for jj in range(no_echoes):
            kdata[:,:,0,jj] = kdata0[jj,:,:]
            img[:,:,0,jj]   = np.fft.ifftshift(np.fft.ifftn(np.squeeze(kdata[:,:,0,jj])))    

    dirlist = os.listdir()
    regechoes_file = [name for name in dirlist if name.lower().startswith('reg') and name.lower().endswith('.nii') ]
    
    if not regechoes_file:
        
        startecho = 1
        
        sys.path.insert(1, 'C:/smis/FLASH_test/Codes/')
        from MRDtoNIfTI import MRDtoNIfTI as M2N
        im_echo, par = M2N(path_RAMSES, mrd_file[0], 1, 'echo_0_forechocount')
        endecho = int(par['NO_ECHOES'].split()[-1])
        os.remove(os.path.join(os.getcwd(), 'echo_0_forechocount.nii'))
        
        echoes_file = [name for name in dirlist if name.lower().startswith('echo') and name.lower().endswith('.nii') ]
        
        if not echoes_file:
            for ii in range(startecho,endecho+1):
                if ii<10:
                    prename = '0'
                else: 
                    prename = ''
                im_echo, par = M2N(path_RAMSES, mrd_file[0], ii, 'echo_'+prename+str(ii))
        im_echo, par = M2N(path_RAMSES, mrd_file[0], startecho, 'echo_0'+str(startecho))        # doing again for echo_01
        affmat = nb.load('echo_0'+str(startecho)+'.nii').affine
        
        affmat = [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
        
        dirlist = os.listdir()
        echoes_file = [name for name in dirlist if name.lower().startswith('echo') and name.lower().endswith('.nii') ]
        
        ''' if there are no registered files:'''
        data = np.zeros( np.hstack( [np.array(nb.load(echoes_file[1]).dataobj).shape, len(echoes_file)] ) )
        dataR= np.zeros( np.hstack( [np.array(nb.load(echoes_file[1]).dataobj).shape, len(echoes_file)] ) )
        dataI= np.zeros( np.hstack( [np.array(nb.load(echoes_file[1]).dataobj).shape, len(echoes_file)] ) ) 
        kk = 0
        kkk = 0
        for i in echoes_file:
            if i == echoes_file[0] and mrd_file[0].lower().startswith('ram'): # ignore first signal (AFI1) in RAMSES
                if 'RAMSES_13' in sequencename:
                    data[:,:,:,kk]= np.flip(np.array(abs(img[:,:,:,kkk])))
                    dataR[:,:,:,kk]= np.flip(np.array(np.real(img[:,:,:,kkk]))) 
                    dataI[:,:,:,kk]= np.flip(np.array(np.imag(img[:,:,:,kkk]))) 
                else:
                    data[:,:,:,kk]= np.flip(np.flip(np.array(abs(img[:,:,:,kkk])),1),0)
                    dataR[:,:,:,kk]= np.flip(np.flip(np.array(np.real(img[:,:,:,kkk])),1),0)
                    dataI[:,:,:,kk]= np.flip(np.flip(np.array(np.imag(img[:,:,:,kkk])),1),0)
            else:
                if mrd_file[0].lower().startswith('ram'): # or mrd_file[0].lower().startswith('mge')):
                    if not np.mod(kk,2):
                        data[:,:,:,kk]= np.array(abs(img[:,:,:,kkk])) 
                        dataR[:,:,:,kk]= np.array(np.real(img[:,:,:,kkk])) 
                        dataI[:,:,:,kk]= np.array(np.imag(img[:,:,:,kkk])) 
                    else:
                        data[:,:,:,kk]= np.flip(abs(img[:,:,:,kkk]),1)
                        dataR[:,:,:,kk]= np.flip(np.real(img[:,:,:,kkk]),1)
                        dataI[:,:,:,kk]= np.flip(np.imag(img[:,:,:,kkk]),1)                    
                else:
                    if not np.mod(kk,2):
                        data[:,:,:,kk]= np.flip(abs(img[:,:,:,kkk]),1)       
                        dataR[:,:,:,kk]= np.flip(np.real(img[:,:,:,kkk]),1)
                        dataI[:,:,:,kk]= np.flip(np.imag(img[:,:,:,kkk]),1)                    
                    else:
                        data[:,:,:,kk]= np.array(abs(img[:,:,:,kkk]))
                        dataR[:,:,:,kk]= np.array(np.real(img[:,:,:,kkk])) 
                        dataI[:,:,:,kk]= np.array(np.imag(img[:,:,:,kkk])) 
                nb.save(nb.Nifti1Image(data[:,:,:,kk], affmat), path_RAMSES +'rotM' + echoes_file[kkk]) 
                nb.save(nb.Nifti1Image(dataR[:,:,:,kk], affmat), path_RAMSES +'rotR' + echoes_file[kkk]) 
                nb.save(nb.Nifti1Image(dataI[:,:,:,kk], affmat), path_RAMSES +'rotI' + echoes_file[kkk]) 
            kk += 1
            kkk += 1  
                    
            
            
            ## THE PROBLEM IS THE REGISTRATION OF REAL AND IMAGINARY PARTS
            
        ''' registration to overcome jumps '''
        dirlist = os.listdir()
        maxs_move = 5
        nb.save(nb.Nifti1Image(data[:,:,:,0], affmat), path_RAMSES +'regM' + echoes_file[0]) 
        nb.save(nb.Nifti1Image(dataR[:,:,:,0], affmat), path_RAMSES +'regR' + echoes_file[0]) 
        nb.save(nb.Nifti1Image(dataI[:,:,:,0], affmat), path_RAMSES +'regI' + echoes_file[0]) 
        for jj in range(0,data.shape[-1]-1):
            _, res = reg_circshift(data[:,:,:,jj]**2, data[:,:,:,jj+1]**2, maxs_move, 'MSE') #seems to be working!
            xx,yy = np.where(res == np.min(res))
            data[:,:,:,jj+1] = np.roll(data[:,:,:,jj+1],yy-maxs_move, axis=0)
            data[:,:,:,jj+1] = np.roll(data[:,:,:,jj+1],xx-maxs_move, axis=1)
            dataR[:,:,:,jj+1] = np.roll(dataR[:,:,:,jj+1],yy-maxs_move, axis=0)
            dataR[:,:,:,jj+1] = np.roll(dataR[:,:,:,jj+1],xx-maxs_move, axis=1)
            dataI[:,:,:,jj+1] = np.roll(dataI[:,:,:,jj+1],yy-maxs_move, axis=0)
            dataI[:,:,:,jj+1] = np.roll(dataI[:,:,:,jj+1],xx-maxs_move, axis=1)
            nb.save(nb.Nifti1Image(data[:,:,:,jj+1], affmat), path_RAMSES +'regM' + echoes_file[jj+1]) 
            nb.save(nb.Nifti1Image(dataR[:,:,:,jj+1], affmat), path_RAMSES +'regR' + echoes_file[jj+1]) 
            nb.save(nb.Nifti1Image(dataI[:,:,:,jj+1], affmat), path_RAMSES +'regI' + echoes_file[jj+1]) 
            
    ''' registered data upload'''
    dirlist = os.listdir()
    reg_echoes_file = [name for name in dirlist if name.lower().startswith('reg') and name.lower().endswith('.nii') ]        
    data = np.zeros( np.hstack( [np.array(nb.load(reg_echoes_file[1]).dataobj).shape, no_echoes] ) )
    dataR = np.zeros( np.hstack( [np.array(nb.load(reg_echoes_file[1]).dataobj).shape, no_echoes] ) )
    dataI = np.zeros( np.hstack( [np.array(nb.load(reg_echoes_file[1]).dataobj).shape, no_echoes] ) )

    rotsuffix = 'MRI'        
    kkk = 0
    jjj = 0
    iii = 0
    for kk in range(0,3):
        reg_echoes_file = [name for name in dirlist if name.startswith('reg'+ rotsuffix[kk])]               
        for i in reg_echoes_file:
            if kk == 0:
                data[:,:,:,kkk]= np.array(nb.load(i).dataobj)
                kkk +=1
            if kk == 1:
                dataR[:,:,:,jjj]= np.array(nb.load(i).dataobj) 
                jjj +=1
            if kk == 2:
                dataI[:,:,:,iii]= np.array(nb.load(i).dataobj)
                iii +=1
        


    ''' deleting original echoes (non registered) and rotated ones '''
    dirlist = os.listdir()
    deletefiles = [name for name in dirlist if not name.lower().startswith('reg') and name.lower().endswith('.nii') ]        
    for item in deletefiles[1:]:
        os.remove(os.path.join(os.getcwd(), item))

    return data, dataR, dataI, kdata, matsize, parMRD
