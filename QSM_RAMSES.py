# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 12:29:24 2021

@author: MZampini

Even-odd phase correction and QSM map calculation

- Phase correction: (based on Li2015 eq1?)
According to Li2015 eq1, the constant phase offset phi0(r) could be obtained by
a linear regression as the intercept of phi(r,TE), while detheta(r)
Problems: unwrapping

- QSM map calculation: (based on Li2015 eq6?)
"""


def kspace2data2D(kdata):
    import numpy as np
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

        if kk == 0: # if RAMSES_14            
            img_K[:,:,:,0] = np.flip(np.flip(np.squeeze(img_K[:,:,:,0])),2)
    return img_K


def kspace2data2D_LR(kdata,Nr):
    import numpy as np
    if len(kdata.shape)>3:
        n3 = kdata.shape[3]
    else:
        n3 = 1
    img_K= np.zeros(np.hstack([kdata.shape[0:3],n3])) +1j*np.zeros(np.hstack([kdata.shape[0:3],n3]))
    for kk in range(n3):
        for ii in range(kdata.shape[2]):
            kdatanow = np.squeeze( np.pad(kdata[:,int(n1/2)-Nr:int(n1/2)+Nr,ii,kk], ((0,0),(int((int(n1/2)-Nr)),int((int(n1/2)-Nr)))), 'constant', constant_values=(0)) )
            img_K[:,:,ii,kk] = np.fft.ifftshift(np.fft.ifftn(kdatanow, s = [kdata.shape[0],kdata.shape[1]]))
        
        if np.mod(kk,2):
            img_K[:,:,:,kk] = np.rot90(np.flip(np.squeeze(img_K[:,:,:,kk])), 2, axes=(0,2))
        if kk == 0: # if RAMSES_14            
            img_K[:,:,:,0] = np.flip(np.flip(np.squeeze(img_K[:,:,:,0])),2)
    return img_K

def kspace2data3D(kdata):
    import numpy as np
    n3 = kdata.shape[3]
    img_K= np.zeros(kdata.shape) +1j*np.zeros(kdata.shape)
    for kk in range(n3):
        img_K[:,:,:,kk] = np.fft.ifftshift(np.fft.ifftshift(np.fft.ifftn(kdata[:,:,:,kk], s = [kdata.shape[0],kdata.shape[1],kdata.shape[2]])), axes=2)
        if np.mod(kk,2):
            img_K[:,:,:,kk] = np.rot90(np.flip(np.squeeze(img_K[:,:,:,kk])), 2, axes=(0,2))
        if kk == 0: # if RAMSES_14            
            img_K[:,:,:,0] = np.flip(np.flip(np.squeeze(img_K[:,:,:,0])),2)
    return img_K 


def kspace2data3D_LR(kdata,Nr1,Nr2):
    import numpy as np
    n3 = kdata.shape[3]
    img_K= np.zeros(kdata.shape) +1j*np.zeros(kdata.shape)
    for kk in range(n3):
        kdatanow = np.squeeze( np.pad(kdata[:,int(n1/2)-Nr1:int(n1/2)+Nr1,int(n2/2)-Nr2:int(n2/2)+Nr2,kk], ((0,0),(int((int(n1/2)-Nr1)),int((int(n1/2)-Nr1))),(int((int(n2/2)-Nr2)),int((int(n2/2)-Nr2))) ), 'constant', constant_values=(0)) )
        img_K[:,:,:,kk] = np.fft.ifftshift(np.fft.ifftshift(np.fft.ifftn(kdatanow, s = [kdata.shape[0],kdata.shape[1],kdata.shape[2]])), axes=2)
        
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


'''
def kspace2data3D(kdata):
    import numpy as np
    img_K= np.zeros(kdata.shape[0:3]) +1j*np.zeros(kdata.shape[0:3]) 
    for ii in range(kdata.shape[3]):
        img_K[:,:,:,ii] = np.fft.ifftshift(np.fft.ifftn(kdata[:,:,:,ii], s = [kdata.shape[0],kdata.shape[1],kdata.shape[2]]))
        if np.mod(ii,2):
            img_K[:,:,:,ii] = np.rot90(np.flip(img_K[:,:,:,ii]), 2, axes=(0,2))         # to be checked
        else:
            img_K[:,:,:,ii] = np.rot90(img_K[:,:,:,ii], 2, axes=(0,1))                  # to be checked
    return img_K
'''



if __name__ == '__main__':
    
    import numpy as np, sys, os, nibabel as nb
    from tkinter import Tk
    from tkinter.filedialog import askdirectory
    #from numpy import pi as pi
    #import getMRD.get_mrd_3d as MRD
    
    
    sys.path.insert(1,'C://Users//MZampini//OneDrive - MR Solutions//Desktop//Deliverable//CodePython')
    #from masking3Dpy import masking3D
    sys.path.insert(1,'C://smis//FLASH_test//Codes')
    #from Susceptibility_QSM_complexv01 import reg_circshift    
    from saveniiMRIdata import RAMSESregistrationcomplex   
    from Susceptibility_QSM_complexv02 import img2phase    
    from masking3Dpy import masking3D


    root = Tk()
    root.wm_attributes('-topmost', 1)
    path_RAMSES = askdirectory(title='Select Folder') 
    path_RAMSES = path_RAMSES + '/'
    root.destroy()
    
    data, dataR, dataI, kdata, matsize, parMRD = RAMSESregistrationcomplex(path_RAMSES)
    acqsize, TEarray = MRDparameters(parMRD,matsize)
    n0,n1,n2,no_echoes = matsize
    no_echoes = no_echoes-1
    
    dataRI = dataR+1j*dataI
    img_phase = img2phase(dataRI)
    img_phasea = img_phase[:,:,:,0::2]
    img_phaseb = img_phase[:,:,:,1::2]
    dataRI = np.abs(dataRI)*np.exp(1j*img_phase)

    
    if int(parMRD['NO_SLICES'].split()[1]) == 1:                # 3D
        datak = kspace2data3D(kdata)
        img_LR = kspace2data3D_LR(kdata,int(n1/8),int(n2/4))
    else:                                                       # 2D
        datak = kspace2data2D(kdata) # needs alignment
        img_LR = kspace2data2D_LR(kdata,int(n0/16))

    dataRI_LR = np.real(img_LR)+1j*np.imag(img_LR)
    img_phase_LR = img2phase(dataRI_LR)
    img_phasea_LR = img_phase_LR[:,:,:,0::2]
    img_phaseb_LR = img_phase_LR[:,:,:,1::2]
    dataRI_LR = np.abs(dataRI_LR)*np.exp(1j*img_phase_LR)
    
    # nb.viewers.OrthoSlicer3D(img_phase_LR)

    TEsa = TEarray[0::2]
    TEsb = TEarray[1::2]
    
    ''' masking based on Otsu's method '''
    from skimage.filters import gaussian, threshold_otsu
    from scipy.ndimage.morphology import binary_closing as binaryclose    
    mask        = np.zeros((n0,n1,n2))
    
    mask00 = masking3D(np.abs(np.squeeze(dataRI[:,:,:,0])))
    mask00 = mask00[0]
    
    for ii in range(n2):
        original = abs(np.squeeze(dataRI[:,:,ii,0])*np.squeeze(mask00[:,:,ii]))
        blurred = gaussian(original, sigma=.8)
        binary = blurred > threshold_otsu(blurred)
        
        binary_closed = binaryclose(binary, iterations=3)
        mask[:,:,ii] = binary_closed
    # nb.viewers.OrthoSlicer3D(mask)
        
     
    '''
    # nb.viewers.OrthoSlicer3D(np.abs(dataRI))  
    mask0 = np.array(np.abs(dataRI[:,:,:,0]))
    maskthr = 0.2
    mask0[mask0<maskthr] = 0
    mask0[mask0>=maskthr] = 1
    # nb.viewers.OrthoSlicer3D(mask0)          
    '''    

    ''' phase unwrapping '''
    img1a_phaseuw =  np.unwrap(img_phase[:,:,:,0::2])
    img1b_phaseuw =   np.unwrap(img_phase[:,:,:,1::2])
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw)
    # nb.viewers.OrthoSlicer3D(img1b_phaseuw)

    
    from skimage.restoration import unwrap_phase
    TEss = max(len(TEsa),len(TEsb))
    img1a_phaseuw0 =  np.zeros((n0,n1,n2,TEss))
    img1b_phaseuw0 =  np.zeros((n0,n1,n2,TEss))
    img1a_phaseuw0_LR =  np.zeros((n0,n1,n2,TEss))
    img1b_phaseuw0_LR =  np.zeros((n0,n1,n2,TEss))
    jj = 0
    while jj<max(len(TEsa),len(TEsb)):
        try:
            img1a_phaseuw0[:,:,:,jj] = unwrap_phase(img_phasea[:,:,:,jj]*mask)
            img1a_phaseuw0_LR[:,:,:,jj] = unwrap_phase(img_phasea_LR[:,:,:,jj]*mask)
        except:
            pass
        try:
            img1b_phaseuw0[:,:,:,jj] = unwrap_phase(img_phaseb[:,:,:,jj]*mask)
            img1b_phaseuw0_LR[:,:,:,jj] = unwrap_phase(img_phaseb_LR[:,:,:,jj]*mask)
        except:
            pass
        jj +=1
    # nb.viewers.OrthoSlicer3D(img1b_phaseuw0)
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw0)
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw0_LR)
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw0[:,:,:,0]/img1b_phaseuw0[:,:,:,0])
        
    ''' phase correction - Li2015 eq.2 '''
    dphimap = img_phase[:,:,:,2] - (img_phase[:,:,:,1] + img_phase[:,:,:,3])/2
    dphimap = img_phase[:,:,:,1] - (img_phase[:,:,:,0] + img_phase[:,:,:,2])/2
    dphimap_uw =  img1a_phaseuw[:,:,:,1] - (img1b_phaseuw[:,:,:,0] + img1b_phaseuw[:,:,:,1])/2
    # nb.viewers.OrthoSlicer3D(dphimap)
    # nb.viewers.OrthoSlicer3D(dphimap_uw)

    ''' mixing with Lu2008 eq.5 - could this work for RAMSES?'''
    theta_map = np.angle(dataRI[:,:,:,0]*dataRI[:,:,:,1].conj())/2
    dataRI_corr = np.zeros(dataRI.shape) +1j*np.zeros(dataRI.shape)
    for ii in range(no_echoes):
        dataRI_corr[:,:,:,ii] = np.abs(dataRI[:,:,:,ii])*np.exp(1j*(img_phase[:,:,:,ii]+(np.mod(ii,2)*2-1)*theta_map))
    # nb.viewers.OrthoSlicer3D(np.abs(dataRI))
    # nb.viewers.OrthoSlicer3D(np.angle(dataRI))
    # nb.viewers.OrthoSlicer3D(theta_map)
    # nb.viewers.OrthoSlicer3D(np.angle(dataRI_corr))
   
    ''' honestly, I would use theta_map from the Lu2008 approach,
    as this does take into account the only the first 2 echoes
    - this can be then used as dtheta of Li2015 eq.2.
    Assuming eddy currents reach a steady state, dtheta can be expressed 
    as Li2015 eq.3.
    This needs to be fitted --> fit linearly a relatively homogeneous phase
    discrepancy map along all 3 spatial directions
    '''
    
    rmap0, rmap1, rmap2 = np.mgrid[0:n0, 0:n1, 0:n2]
    rmap = np.zeros(np.hstack([rmap0.shape,3]))
    rmap[:,:,:,0] = rmap0
    rmap[:,:,:,1] = rmap1
    rmap[:,:,:,2] = rmap2
    # nb.viewers.OrthoSlicer3D(np.abs(rmap))


    ''' and then what?! how to fit Li2015 eq.3??
    datanow = theta_map
    map_m = np.zeros(datanow.shape[0:-1])
    from sklearn.linear_model import LinearRegression
    xs0 = TEsa.reshape((-1, 1))
    for i0 in range(n0):
        for i1 in range(n1):
            for i2 in range(n2):
                if mask[i0,i1,i2]:
                    for ir in range(3):
                        ys = datanow[i0,i1,i2,:]
                        xs = xs0[ys != 0]
                        ys = ys[ys != 0]
                        if not (len(xs)== 0 or sum(ys)==0):                
                            fitlog = LinearRegression().fit(xs, ys)
                            #map_b[i0,i1,i2] = 
                            map_m[i0,i1,i2] = fitlog.coef_    
    '''





    img1a_phaseuw =   np.unwrap(img_phase[:,:,:,1::2])
    img1b_phaseuw =   np.unwrap(img_phase[:,:,:,2::2])
    img1a_phaseuw_fewpts =   np.unwrap(img_phase[:,:,:,1:5:2])
    # nb.viewers.OrthoSlicer3D(img_phase)
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw)
    # nb.viewers.OrthoSlicer3D(img1b_phaseuw)
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw_fewpts)
    
    img_phaseLR = img2phase(img_LR)
    img1a_phaseuw_LR =   np.unwrap(img_phaseLR[:,:,:,1::2])
    img1b_phaseuw_LR =   np.unwrap(img_phaseLR[:,:,:,2::2])
    
    img1a_dephaseuw_LR = np.diff(img1a_phaseuw_LR)
    img1b_dephaseuw_LR = np.diff(img1b_phaseuw_LR)
    # nb.viewers.OrthoSlicer3D(img_phaseLR)
    # nb.viewers.OrthoSlicer3D(img1b_phaseuw_LR)
    # nb.viewers.OrthoSlicer3D(img1a_dephaseuw_LR)
    
            
    ''' weight definition for linear regression - Kressler2010 eq.12'''
    stdRea = np.std(np.real(dataRI[:,:,:,1::2]), axis=3)
    stdReb = np.std(np.real(dataRI[:,:,:,2::2]), axis=3)
    wa = np.repeat(stdRea[:,:,:,np.newaxis],len(TEsa),axis=3)/np.abs(dataRI[:,:,:,1::2])
    wb = np.repeat(stdReb[:,:,:,np.newaxis],len(TEsb),axis=3)/np.abs(dataRI[:,:,:,2::2])
    # nb.viewers.OrthoSlicer3D(1/wa)

    datanowa = img1b_phaseuw0_LR #img1a_phaseuw_LR
    datanowb = img1a_phaseuw0_LR[:,:,:,1:]# img1b_phaseuw_LR
    map_ma = np.zeros(datanowa.shape[0:-1])
    map_Ba = np.zeros(datanowa.shape[0:-1])
    map_mb = np.zeros(datanowb.shape[0:-1])
    map_Bb = np.zeros(datanowb.shape[0:-1])
    from sklearn.linear_model import LinearRegression
    xs0a = TEsa.reshape((-1, 1))
    xs0b = TEsb.reshape((-1, 1))
    for i0 in range(n0):
        for i1 in range(n1):
            for i2 in range(n2):
                if mask[i0,i1,i2]:
                    ysa = datanowa[i0,i1,i2,:]
                    ysa = ysa[ysa != 0]
                    xsa = xs0a[ysa != 0]
                    weighta  = 1/wa[i0,i1,i2]
                    ysb = datanowb[i0,i1,i2,:]
                    ysb = ysb[ysb != 0]
                    xsb = xs0b[ysb != 0]
                    weightb  = 1/wb[i0,i1,i2]
                    if not (len(xsa)== 0 or sum(xsb)==0 or sum(ysa)==0 or sum(ysb)==0):                
                        fitloga = LinearRegression().fit(xsa, ysa, weighta)
                        map_ma[i0,i1,i2] = fitloga.coef_    
                        map_Ba[i0,i1,i2] = fitloga.intercept_
                        fitlogb = LinearRegression().fit(xsb, ysb, weightb)
                        map_mb[i0,i1,i2] = fitlogb.coef_    
                        map_Bb[i0,i1,i2] = fitlogb.intercept_    
    # nb.viewers.OrthoSlicer3D(map_m)
    # nb.viewers.OrthoSlicer3D(map_Ba)
    # nb.viewers.OrthoSlicer3D(map_Bb)
    # nb.viewers.OrthoSlicer3D(map_Bb-map_Ba)
    # nb.save(nb.Nifti1Image(map_B, [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]), path_RAMSES + 'maptheta0.nii')   
      
    dethetaLR = map_Ba-map_Bb        
    # nb.viewers.OrthoSlicer3D(dethetaLR)
    img_phaseb_corr_LR = img1a_phaseuw0_LR + np.repeat(dethetaLR[:,:,:,np.newaxis],len(TEsa),axis=3)
    img_phaseb_corr_LR[:,:,:,0] = img1a_phaseuw0_LR[:,:,:,0]
    # nb.viewers.OrthoSlicer3D(img_phaseb_corr_LR)
    img_phase_corr_LR = np.zeros((n0,n1,n2,len(TEsa)+len(TEsb)))    
    img_phase_corr_LR[:,:,:,0::2] = img_phaseb_corr_LR
    img_phase_corr_LR[:,:,:,1::2] = img1a_phaseuw0_LR[:,:,:,1:]
    # nb.viewers.OrthoSlicer3D(img_phase_corr_LR)


    ''' phase jumps between odd and even echoes have now been corrected.
    Next step: SWI and QSM'''
    
    ''' SQM: calculation of background gradients '''
    
    
    

    ''' removal of background field? '''
    imga_phaseuwHR = img1a_phaseuw0 - img1a_phaseuw0_LR
    imgb_phaseuwHR = img1b_phaseuw0 - img1b_phaseuw0_LR
    # nb.viewers.OrthoSlicer3D(imga_phaseuwHR)

    dphimapLR_uw =  img1a_phaseuw[:,:,:,1] - (img1b_phaseuw[:,:,:,0] + img1b_phaseuw[:,:,:,1])/2
    
    
    
    ''' img_LR0 and img_LR1 are the s_hat_i of eq5 '''
    ''' High order phase error correction: Yu2010 eq.5 '''
    theta_hat_i = np.angle(img_LR0corr*img_LR1corr.conj())/2
    corr_theta_hat = np.exp(1j*theta_hat_i)
    corr_mtheta_hat = np.exp(-1j*theta_hat_i)

    # nb.viewers.OrthoSlicer3D(theta_hat_i)
    imga_corr = dataRI[:,:,:,1::2]*np.repeat(corr_theta_hat[:, :, :, np.newaxis],len(TEsa),axis=3)
    imgb_corr = dataRI[:,:,:,2::2]*np.repeat(corr_mtheta_hat[:, :, :, np.newaxis],len(TEsb),axis=3)
    # nb.viewers.OrthoSlicer3D(np.angle(imga_corr))
    # nb.viewers.OrthoSlicer3D(np.angle(imgb_corr))

    img_corr = np.zeros((n0,n1,n2,no_echoes-1)) + 1j*np.zeros((n0,n1,n2,no_echoes-1))
    img_corr[:,:,:,0::2] = imga_corr
    img_corr[:,:,:,1::2] = imgb_corr
    # nb.viewers.OrthoSlicer3D(np.angle(img_corr))
    # nb.viewers.OrthoSlicer3D(np.angle(dataRI))


    #from skimage.measure import label, regionprops, regionprops_table
    #props = regionprops_table(label_img, properties=('centroid','major_axis_length','minor_axis_length'))
    #import scipy
    #mask_fill = scipy.ndimage.morphology.binary_fill_holes(mask)
    # nb.viewers.OrthoSlicer3D(mask_fill)      


    
    
    from skimage.restoration import unwrap_phase      # phase unwrapping
    img1a_phaseuw_ = np.zeros((n0,n1,n2,len(TEsa)))
    img1b_phaseuw_ = np.zeros((n0,n1,n2,len(TEsb)+1))
    for ii in range(n2):
        for jj in range(len(TEsb)+1):
            img1b_phaseuw_[:,:,ii,jj] = unwrap_phase(img_phaseb[:,:,ii,jj]*mask0[:,:,ii])
        for jj in range(len(TEsa)):
            img1a_phaseuw_[:,:,ii,jj] = unwrap_phase(img_phasea[:,:,ii,jj]*mask0[:,:,ii])
    
    # nb.viewers.OrthoSlicer3D(img_phase) 
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw_) 
    # nb.viewers.OrthoSlicer3D(imgaphase) 

    img1a_phaseuw__ = np.zeros((n0,n1,n2,len(TEsa)))
    img1b_phaseuw__ = np.zeros((n0,n1,n2,len(TEsb)+1))
    for jj in range(len(TEsb)+1):
        img1b_phaseuw__[:,:,:,jj] = unwrap_phase(img_phaseb[:,:,:,jj]*mask0)
    for jj in range(len(TEsa)):
        img1a_phaseuw__[:,:,:,jj] = unwrap_phase(img_phasea[:,:,:,jj]*mask0)    
    
        
    # nb.viewers.OrthoSlicer3D(img_phase) 
    # nb.viewers.OrthoSlicer3D(img1a_phaseuw__) 
    # nb.viewers.OrthoSlicer3D(img1b_phaseuw__) 
    # nb.viewers.OrthoSlicer3D(imgaphase) 
    
    
    
    
    ''' F) Baudrexel2009 eq 2 and OrosPeusquens2019 eq 5'''


    imgaphase = img1a_phaseuw
    imgbphase = img1b_phaseuw
    
    gamma = 42.577478518*1e6 # Hz/T
    Gsusca = (imgaphase[:,:,:,1]-imgaphase[:,:,:,0])/(gamma*(TEsa[1]-TEsa[0])*matsize[2])
    
    # nb.viewers.OrthoSlicer3D(Gsusca)
    
    
    Gsuscatot = np.zeros((n0,n1,n2,len(TEsa)-1,3))
    Gsusca    = np.zeros((n0,n1,n2,3))
    Gsuscbtot = np.zeros((n0,n1,n2,len(TEsb)-1,3))
    Gsuscb    = np.zeros((n0,n1,n2,3))
    
    for idir in range(0,3):
        for iTE in range(len(TEsa)-1):
            imgaphase1 = np.roll(imgaphase[:,:,:,iTE+1], 1, axis=idir)
            Gsuscatot[:,:,:,iTE,idir] = ( imgaphase1 - imgaphase[:,:,:,0] ) / (gamma * matsize[idir] * TEsa[iTE+1])
        for iTE in range(len(TEsb)-1):
            imgbphase1 = np.roll(imgbphase[:,:,:,iTE+1], 1, axis=idir)
            Gsuscbtot[:,:,:,iTE,idir] = ( imgbphase1 - imgbphase[:,:,:,0] ) / (gamma * matsize[idir] * TEsb[iTE+1])
            
        Gsusca[:,:,:,idir] = np.mean(Gsuscatot[:,:,:,iTE,:],axis=3) # assuming linearity in phase variation
        Gsuscb[:,:,:,idir] = np.mean(Gsuscbtot[:,:,:,iTE,:],axis=3) # assuming linearity in phase variation
    
    # nb.viewers.OrthoSlicer3D(Gsusca)
    # nb.viewers.OrthoSlicer3D(Gsuscb)
    # nb.viewers.OrthoSlicer3D(imgaphase1)
        
    import scipy    
    Gsusca_med = scipy.ndimage.median_filter(Gsusca[:,:,:,2], size=(5,5,1))    
    Gsusca_rep = np.repeat(Gsusca_med[:, :, :, np.newaxis], len(TEsa), axis=3)
    imga_corr = dataRI[:,:,:,0::2]
    sincterm = np.sinc(gamma/2*np.squeeze(Gsusca_rep)*matsize[2]*TEsa)
    imga_corr = imga_corr/sincterm
    imga_corr = abs(imga_corr)
    imga_corr[imga_corr>0.5] = .5
    
    # nb.viewers.OrthoSlicer3D(np.abs(sincterm))
    # nb.save(nb.Nifti1Image(abs(imga_corr[:,:,:,0]), [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]), path_RAMSES + 'imga_corr.nii')   


    import time
    for i in range(3):
        print('\007')
        time.sleep(0.5)

    