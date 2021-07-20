# Description: Port of Get_MRD.m to Python
# Author(s): MR Solutions Ltd (Tom Wilkinson, William Scott-Jackson)
# Copyright: MR solutions Ltd 2014-2018
# Date: 03/03/2015 - created
# Description: Function to open multidimensional mrd images given a filename
# with PPR-parsing
# Original Author: MR Solutions Ltd (Ruslan Garipov)
# Date: 01/03/2014 - swapped views and views2 dimension - now correct
# Date: 05/08/2016 - PEP8 Compliant, added copyright notice
# Date: 24/08/2017 - Fixed Switch-Case block, change condition values to hexidecimal
# Date: 04/09/2017 - Fixed Centric ordering bug (convereted to zero index arthmetic)
# Date: 04/09/2017 - Added Support for 3D data with second degree of centric ordering
# Date: 04/09/2017 - Removed some redudant computations
# Data: 20/09/2017 - Made fixes to PPR parsing, added functions and improved structure, added automatic centric ordering
# Data: 14/11/2019 - MAZ Scaled im with the scaling factor -- removed! Does not work

# Usage example: 
# import get_MRD as MRD 
#
# filename = "example_0.MRD"
#
# k_space_data, k_space_dims, PPR_Params = MRD.get_mrd_3d(filename)
#
# Where k_space_data is a 6 dimensional matrix where the data are ordered by experiments, echoes, slices, views, views2 and samples
# k_space_dims is a 1x6 array which tells you the dimensions, ordered as [nExperiments, nEchoes, nSlices, nViews, nViews2, nSamples]
# PPR_Params is a python dictionary with various parameters pulled from the footer at the bottom of the MRD. You can see the parameters by 
# using print(PPR_Params)
# PPR_Keywords should be in the same folder get_mrd_3d is found

# Import the required modules
from scipy import complex_
import mmap
import numpy as np
import array
from struct import unpack
try:
    from . import PPR_Keywords as PPR
except:#ImportError
    import PPR_Keywords as PPR

def get_dataformat(dt):
    if dt is 0:
        return 'B'
        # datasize = 1
    elif dt is 1:
        return 'b'
        # datasize = 1
    elif dt is 2:
        return 'h'
        # datasize = 2
    elif dt is 3:
        return 'h'
        # datasize = 2
    elif dt is 4:
        return 'l'
        # datasize = 4
    elif dt is 5:
        return 'f'
        # datasize = 4
    elif dt is 6:
        return 'd'
        # datasize = 8
    else:
        return 'i'
        # datasize = 4
    
def get_mrd_3d(filename_):
        
    # Open file as read-binary only
    fidf = open(filename_, 'r+b')
    fid = mmap.mmap(fidf.fileno(), 0)
    
    # Read first 4 values from header
    val = unpack('llll', fid.read(16))
    
    # Get dimensions from this
    no_samples, no_views, no_views_2, no_slices = val[0], val[1], val[2], val[3]
    
    # Work out datatype of data
    fid.seek(18)
    datatype_ = unpack('h', fid.read(2))[0]
    fid.seek(48)
    scaling = unpack('f', fid.read(4))
    bitsperpixel = unpack('f', fid.read(4))
    
    fid.seek(152)
    val = unpack('ii', fid.read(8))
    no_echoes = val[0]
    no_expts = val[1]

    fid.seek(256)
    text = fid.read(256)

    dim = [no_expts, no_echoes, no_slices, no_views_2, no_views, no_samples]

    dt = datatype_
   
    if dt >= 0x10:
        iscomplex = 2
        dt = dt-0x10
    else:
        iscomplex = 1
    
    dataformat = get_dataformat(dt)

    # Compute the number of values expected to be read from the dimensions
    num2read = no_expts*no_echoes*no_slices*no_views*no_views_2 * no_samples*iscomplex
   
    fid.seek(512)
    
    m_total = array.array(dataformat)
    m_total.fromfile(fid, num2read)
    if len(m_total) != num2read:
        print("We have a problem...(file length/read mismatch)")
        return 0
    
    par = PPR.ParseKeywords(fid)
    fid.close()

    if iscomplex is 2:
        m_real = m_total[::2]
        m_imag = m_total[1::2]
        m_C = np.vectorize(complex)(m_real, m_imag)
        m_real = None
        m_imag = None
        m_total = None
    else:
        m_C = m_total
        m_total = None

    n = 0
    
    ord_ = list(range(no_views))

    if 'VAR_centric_on' in par:
        if int(par['VAR centric_on']) == 1:
            val = int(no_views * 0.5)
            ord_ = list(range(no_views))
            for x in range(val):
                ord_[2*x] = val + x
                ord_[2*x+1] = val - x - 1
    elif 'VAR pe1_order' in par:
        if int(par['VAR pe1_order']) == 1:
            val = int(no_views * 0.5)
            ord_ = list(range(no_views))
            for x in range(val):
                ord_[2*x] = val + x
                ord_[2*x+1] = val - x - 1
    
    ord2 = list(range(no_views_2))
    
    if 'pe2_centric_on,' in par:

        if int(par['pe2_centric_on,']) == 1:
            val = int(no_views_2 * 0.5)
            ord2 = list(range(no_views_2))
            for x in range(val):
                ord2[2*x] = val + x
                ord2[2*x + 1] = val - x - 1
            
    im = np.zeros(shape=(no_expts, no_echoes, no_slices, no_views, no_views_2,
                         no_samples), dtype=complex_)

    for a in range(no_expts):
        for b in range(no_echoes):
            for c in range(no_slices):
                for d in range(no_views):
                    for e in range(no_views_2):
                        for f in range(no_samples):
                            im[a][b][c][ord_[d]][ord2[e]][f] = m_C[f+n]
                        n += no_samples

    ord_ = None
    
             
    #print(scaling)
    #im = im/scaling  # maybe works only for SUR
    #print(float(str(scaling)))
    
    return im, dim, par, scaling
