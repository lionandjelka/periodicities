from superlet import superlet
from superlet import gen_superlet_testdata, 



def superlets_methods(tt, mag,ntau, f = 2, minfq = 500, maxfq = 10,  method = 'linear'):
    """Perform hybrid2d method using superlets

        Parameters
        ----------
        tt : list of time data
        mag : list of magnitude values
        ntau, ngrid : values for controling wwz execution (see inp_param function)
        minfq : minimum frequency
        maxfq : maximum fraquency

    """
    flist=np.arange(minfq,maxfq,ntau)
    scales1=scale_from_period(1/flist)
    gg=superlet(
        m,
        samplerate=1/sampleratess,
        scales=scales1,
        order_max=100,
        #order_min=1, default
        order_min=10,
        #c_1=3 default
        c_1=3,
        adaptive=True,
    )
    gg1=np.abs(gg)
    corr = correlation_nd(gg1, gg1)
    
    extentmin=np.min(corr)
    extentmax=np.max(corr)

    extent=[extentmin,extentmax,extentmin,extentmax]
    
    

    return  corr,  extent
