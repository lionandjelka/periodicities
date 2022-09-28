
import numpy as np
import json
from scipy.signal import chirp, find_peaks, peak_widths
from traitlets.traitlets import Integer
from scipy import interpolate, optimize
from scipy.stats.mstats import mquantiles
from sklearn.utils import shuffle
from periodicity.utils.correlation import correlation_nd
from periodicity.algorithms.wavelets import *

def get_full_width(x: np.ndarray, y: np.ndarray, peak:np.ndarray,height: float = 0.5) -> float:
    """Function is used to calculate the error of the determined period using FWHM method.
        The period uncertainty method of Schwarzenberg-Czerny requires the so called Mean Noise Power Level (MNPL)
        in the vicinity of P. The 1-sigma confidence interval on P then is equal to the width of the line at the P – MNPL level. 
        This method is a so-called ‘post-mortem analysis’. 
        To find the MNPL we detect FWHM of the peak and then use mquantile method to detect points between 25th and 75th quantile 
        Reference: Schwarzenberg-Czerny, A., 1991, Mon. Not. R. astr. Soc., 253, 198-206
        
        Parameters
        ----------
        x,y: numpy array, arrays with data
        peak : numpy array, array with determined peaks
        height: size of peak
        
        Returns:
        ----------
        Arrays with results, and quantiles and half peak size, and low/high x values

    """

    er1=[]
    er3=[]
    height_half_max = 0
    quantiles = []
    phmax =  []
    x_lows =[]
    x_highs = []
    # Iterate over all peaks
    for i in range(len(peak)):
        
        # Half max of a peak as its y axis value multiplied by height parameter (by default 0.5)
        height_half_max=y[peak[i]]*height
        index_max = peak[i]
        
        # Determine the x-axis border of a peak but on its half max height    
        tmp = peak[i]
        tmp2 = peak[i]
                                          
        x_low = 0
        x_high = 0
        while True:
            tmp = tmp - 1;
            if( y[tmp] - height_half_max) < 0:
                x_low = x[tmp+1]
                break
        
        while True:
            tmp = tmp + 1;
            if( y[tmp] - height_half_max) < 0:
                x_high = x[tmp-1]
                break
        
            
        # q1 and q2 represents quantiles    
        q25 = 0
        q75 = 0
        
        # Check if we have sufficient number of points (5 points)
        
        if( index_max - 5 > 0 ):

            arr=y[(x>=x_low)&(x<=x_high)]  # array of data between x_low and x_hight where we determine quantiles
            # For a reference see  Schwarzenberg-Czerny, A., 1991, Mon. Not. R. astr. Soc., 253, 198-206
            
            q25,q75 = mquantiles(arr, [0.25,0.75])
 
            # Calculate the value at specific quantiles (q1 and q3)
            inversefunction = interpolate.interp1d(y[index_max-5:index_max],  x[index_max-5:index_max], kind='cubic',fill_value="extrapolate")
            inversefunction2 = interpolate.interp1d(y[index_max:index_max+5],  x[index_max:index_max+5], kind='cubic',fill_value="extrapolate")
            
            xer1=inversefunction(q25)
            xer3=inversefunction2(q75)
            
        else:
            xer1 = 0
            xer3 = 0
            
        er1.append(xer1)
        er3.append(xer3)
        quantiles.append([q25, q75])
        phmax.append(height_half_max)
        x_lows.append(x_low)
        x_highs.append(x_high)
        
        
    return er1,er3, quantiles, phmax, x_lows, x_highs



def periods (lcID, data, ngrid, plot = False, save = False, peakHeight = 0.6, prominence = 0.7, minfq = 500, maxfq = 10, xlim = None): 
    """Perform period determination for the output of hybrid2d data.

        Parameters
        ----------
        lcId : id of a curve
        data :auto correlation matrix
        ngrid : values for controling wwz execution (see inp_param function)
        minfq : minimum frequency
        maxfq : maximum fraquency
        peakHeight: max peak height
        prominence: prominence for peak determination
        plot: True of Folse if we want a plot
        save: determine to save a plot
        xlim: set xlim for a plot
        
        Returns:
        ---------
        r_peaks: arrays with determined periods
        r_peaks_err_upper: arrays with upper errors of corresponding periods
        r_peaks_err_lower: arrays with lower errors of corresponding periods

    """
    
    hh1=np.rot90(data).T/np.rot90(data).T.max()
    hh1arr=np.rot90(hh1.T)
    hh1arr1=np.abs(hh1arr).sum(1)/np.abs(hh1arr).sum(1).max()
    
    
    
    fmin =1/minfq
    fmax = 1/maxfq
    df = (fmax - fmin) / ngrid
    
    
    
    # osax interpolation (stacked data of h2d along one axis) to obtain more points
    osax=np.arange(start=fmin,stop=fmax+df,step=df)
    xax =np.arange(start=fmin, stop = fmax + df, step = df/2)
    from scipy import interpolate
    f = interpolate.interp1d(osax, np.abs(hh1arr1), fill_value="extrapolate")
    yax = []
    for v in xax:
        yax.append(float(f(v)))
    yax = np.array(yax)
    
    
    # Finding peaks
    peaks,_ = find_peaks(yax, peakHeight, prominence = prominence)

    
    # Polotting if needed
    if( plot == True ):       
        if xlim != None:
            plt.xlim(xlim)
        plt.plot(xax,np.abs(yax))
        plt.axvline(xax[peaks[0]],ymin=0,ymax=1, linestyle='--', color='k')
        plt.title(str(lcID))
        plt.xlabel(r'Freqeuncy [day$^{-1}$]')
        plt.ylabel(r"correlation")
        if( save == True) :
            plt.savefig(str(lcID) + 'stackd_h2d.png')  
        
        
    
    #Get error estimates for each peak (period)
    error_upper, error_lower, quantiles, halfmax, x_lows, x_highs = get_full_width (xax,yax,peaks)
    
    if plot == True:
        plt.plot(xax,np.abs(yax))
        if xlim != None:
            plt.xlim(xlim)
        plt.title(str(lcID))
        plt.xlabel(r'Freqeuncy [day$^{-1}$]')
        plt.ylabel(r"correlation")
        
        for i in range(len(peaks)):
            plt.axvline(xax[peaks[i]],ymin=0,ymax=1, linestyle='--', color='black')
            plt.axhline(quantiles[i][0],linestyle='--', color='green')
            plt.axhline(quantiles[i][1],linestyle='--', color='red')
            plt.axvline(x_lows[i], ymin =0, ymax = 1, linestyle ='--', color='blue')
            plt.axvline(x_highs[i], ymin =0, ymax = 1, linestyle ='--', color='blue')
            plt.axhline(halfmax[i],linestyle='--', color='purple')
        
        
        if( save == True) :
            plt.savefig(str(lcID) + 'stackd_h2d_peaks.png')  
        
        
    
    
    # Prepare the output
    r_peaks = []
    r_peaks_err_upper = []
    r_peaks_err_lower = []
    for i in range(len(peaks)):
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
    
    
    
def signif_johnoson(numlc, peak, corr,  tt, yy, ntau,ngrid, f = 2, peakHeight = 0.6, minfq = 500, maxfq = 10, algorithm ='wwz', method = 'linear'):
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
        peakHeight: max peak height,
        algorithm: wwz or superlets (for now)
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
        if( algorithm == 'wwz'):
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