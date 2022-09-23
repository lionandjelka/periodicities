
import numpy as np
import json
from scipy.signal import find_peaks
from traitlets.traitlets import Integer
from scipy import interpolate, optimize
from scipy.stats.mstats import mquantiles
from sklearn.utils import shuffle
from periodicity.utils.correlation import correlation_nd
from periodicity.algorithms.wavelets import *

def get_full_width(x: np.ndarray, y: np.ndarray, peak:np.ndarray,height: float = 0.5) -> float:
    er1=[]
    er3=[]
    from scipy.signal import chirp, find_peaks, peak_widths
    results_half = peak_widths(y, peak)
  
    
    for i in range(len(peak)):
        height_half_max=y[peak[i]]*height
        index_max = peak[i]
        
        x_low = np.interp(height_half_max, y[:index_max-3:index_max],  x[:index_max-3:index_max])
        x_high = np.interp(height_half_max, np.flip(y[index_max::index_max+3]), np.flip(x[index_max::index_max+3]))
        
        if( index_max - 5 > 0 ):
            inversefunction = interpolate.interp1d(y[index_max-5:index_max],  x[index_max-5:index_max], kind='cubic',fill_value="extrapolate")
            inversefunction2 = interpolate.interp1d(y[index_max:index_max+5],  x[index_max:index_max+5], kind='cubic',fill_value="extrapolate")

            arr=y[(x>=x_low)&(x<=x_high)]   
            
            q1,q3 = mquantiles(arr, [0.25,0.75])
        
            xer1=inversefunction(q1)
            xer3=inversefunction2(q3)
        else:
            xer1 = 0
            xer3 = 0
        er1.append(xer1)
        er3.append(xer3)
        
    return er1,er3



def periods (data,ngrid, plot = False, peakHeight = 0.6,  minfq = 500, maxfq = 10): 
    """Perform period determination for the output of hybrid2d data

        Parameters
        ----------
        data :auto correlation matrix
        ngrid : values for controling wwz execution (see inp_param function)
        minfq : minimum frequency
        maxfq : maximum fraquency
        peakHeight: max peak height
        plot: True of Folse

    """
    
    hh1=np.rot90(data).T/np.rot90(data).T.max()
    hh1arr=np.rot90(hh1.T)
    hh1arr1=np.abs(hh1arr).sum(1)/np.abs(hh1arr).sum(1).max()
    
    
    
    fmin =1/minfq
    fmax = 1/maxfq
    df = (fmax - fmin) / ngrid
    
    
    
    osax=np.arange(start=fmin,stop=fmax+df,step=df)

    xax =np.arange(start=fmin, stop = fmax + df, step = df/2)
    from scipy import interpolate
    f = interpolate.interp1d(osax, np.abs(hh1arr1), fill_value="extrapolate")
        
        
    yax = []
    for v in xax:
        yax.append(float(f(v)))
    yax = np.array(yax)
    
    
    peaks,_ = find_peaks(yax,peakHeight, prominence = 0.7)
    
    npeak = len(peaks)
    
    
  
    if( plot == True ):
        plt.plot(xax,np.abs(yax))
        plt.axvline(xax[peaks[0]],ymin=0,ymax=1, linestyle='--', color='k')
        plt.xlabel(r'Freqeuncy [day$^{-1}$]')
        plt.ylabel(r"correlation")
        
    error_upper, error_lower=get_full_width(xax,yax,peaks)
    
    r_peaks = []
    r_peaks_err_upper = []
    r_peaks_err_lower = []
    for i in range(npeak):
        r_peaks.append(1/xax[peaks[i]])
        if( error_upper[i] == 0):
            r_peaks_err_upper.append(-1)
        else:
            r_peaks_err_upper.append(np.abs(1/xax[peaks[i]]-(1/error_upper[i])))
                
        if(error_lower[i] == 0):
            r_peaks_err_lower.append(-1)
        else:
            r_peaks_err_lower.append(np.abs(1/xax[peaks[i]]-( 1/error_lower[i])))
    
    
        
    return r_peaks, r_peaks_err_upper, r_peaks_err_lower
    
    
    
def signif_johnoson(numlc, peak, corr,  tt, yy, ntau,ngrid, f = 2, peakHeight = 0.6, minfq = 500, maxfq = 10,  method = 'linear'):
    """Determination of significance usign Johnson method

        Parameters
        ----------
        numlc : int, number of lc for determination
        peak : determined periodicity peak
        corr : hybrid 2d output
        tt : time
        yy: magnitude
        plot: True of Folse
        ntau, ngrid, f : values for controling wwz execution (see inp_param function)
        minfq : minimum frequency
        maxfq : maximum fraquency
        peakHeight: max peak height
    """
    hh1=np.rot90(corr).T/np.rot90(corr).T.max()
    hh1arr=np.rot90(hh1.T)
    hh1arr1=np.abs(hh1arr).sum(1)/np.abs(hh1arr).sum(1).max()
    peaks,_ = find_peaks(hh1arr1,peakHeight, prominence = 0.8)
    
    if peak > len(peaks):
        return None
    
    idxrep=peaks[peak]
    #peak power larger than red noise peak power
    count=0.
    #peak power of red noise larger than observed peak power
    count11=0.
    bins11=[]
    bins=[]
   
    for i in range(numlc):
        y = shuffle(yy)
        ntau,params,decay_constant, parallel=inp_param(ntau, ngrid, f, minfq, maxfq)
        wwt_removedx = libwwz.wwt(timestamps=tt,
                             magnitudes=y,
                             time_divisions=ntau,
                             freq_params=params,
                             decay_constant=decay_constant,
                             method='linear',
                             parallel=parallel)
        corr1x=correlation_nd(np.rot90(wwt_removedx[2]),np.rot90(wwt_removedx[2]))
        hhx=np.rot90(corr1x).T/corr1x.max()
        hh1x=np.rot90(hhx.T)
        hh1xarr=np.abs(hh1x).sum(1)/np.abs(hh1x).sum(1).max()

        bins.append(hh1xarr[idxrep])
        if ((hh1arr1[idxrep]/hh1xarr[idxrep])>1.):
            count=count+1.
        else:
            count11=count11+1.  
            bins11.append(hh1xarr[idxrep])
            
    
    return bins, bins11, count/numlc, count11/numlc