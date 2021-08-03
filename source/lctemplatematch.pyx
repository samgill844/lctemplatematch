import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sin, asin, cos, acos, sqrt, M_PI, fmod, abs, tan, atan, atan2, floor, pow, log
from cython.parallel import prange
from scipy.stats import  sem




############################################################################################
#                               Time immediatly prior to periastron passage                #
############################################################################################
cdef _t_ecl_to_peri(double t_ecl, double e, double w, double incl, double p_sid, double radius_1) :
    
    # Define variables used
    cdef double efac  = 1.0 - e*2       # check if 2 or power
    cdef double sin2i = sin(incl)**2.

    # Value of theta for i=90 degrees
    cdef double theta_0 = (M_PI/2.) - w       # True anomaly at superior conjunction
    cdef double theta = theta_0 # BODGE!

    cdef double ee
    if (theta == M_PI) :  ee = M_PI
    else :  ee =  2.0 * atan(sqrt((1.-e)/(1.0+e)) * tan(theta/2.0))

    cdef double eta = ee - e*sin(ee)
    cdef double delta_t = eta*p_sid/(2*M_PI)
    return t_ecl  - delta_t 

############################################################################################
#                               Anomalies                                                  #
############################################################################################
cdef _get_Mean_Anomaly(double time, double t_ecl_to_peri, double period, double e=0.) : 
    if e < 1e-5 : return ((time - t_ecl_to_peri)/period - floor(((time - t_ecl_to_peri)/period)))*2.*M_PI
    else        : return 2*M_PI*((time -  t_ecl_to_peri  )/period % 1.)


cdef _get_Eccentric_Anomaly(double Mean_Anomaly, double eccentricity) : 
    cdef double M = fmod(Mean_Anomaly , 2*M_PI)
    cdef int flip = 0
    if eccentricity == 0 : return M
    if M > M_PI:
        M = 2*M_PI - M
        flip = 1
    cdef double alpha = (3*M_PI + 1.6*(M_PI-abs(M))/(1+eccentricity) )/(M_PI - 6/M_PI)
    cdef double d = 3*(1 - eccentricity) + alpha*eccentricity
    cdef double r = 3*alpha*d * (d-1+eccentricity)*M + M**3
    cdef double q = 2*alpha*d*(1-eccentricity) - M**2
    cdef double w = (abs(r) + sqrt(q**3 + r**2))**(2/3)
    cdef double E = (2*r*w/(w**2 + w*q + q**2) + M) / d
    cdef double f_0 = E - eccentricity*sin(E) - M
    cdef double f_1 = 1 - eccentricity*cos(E)
    cdef double f_2 = eccentricity*sin(E)
    cdef double f_3 = 1-f_1
    cdef double d_3 = -f_0/(f_1 - 0.5*f_0*f_2/f_1)
    cdef double d_4 = -f_0/(f_1 + 0.5*d_3*f_2 + (d_3**2)*f_3/6)
    E = E -f_0/(f_1 + 0.5*d_4*f_2 + d_4**2*f_3/6 - d_4**3*f_2/24)
    if flip==1 : E =  2*M_PI - E
    return E

cdef _get_True_anamoly(double Eccentric_Anomaly, double eccentricity) : return 2.*atan(sqrt((1.+eccentricity)/(1.-eccentricity))*tan(Eccentric_Anomaly/2.))

############################################################################################
#                               Projected seperation                                       #
############################################################################################

cdef _getProjectedSeperation(double nu, double e, double incl, double w, double radius_1) : 
    return (1-e**2) * sqrt( 1.0 - sin(incl)**2  *  sin(nu + w)**2) / (1 + e*cos(nu)) /radius_1

@cython.boundscheck(False) # turn off bounds-checking for entire function
@cython.wraparound(False)  # turn off negative index wrapping for entire function
@cython.nonecheck(False)  # turn off negative index wrapping for entire function
cdef __getProjectedSeperation(double * nu, double e, double incl, double w, double radius_1, int n, double * projected_seperations) : 
    cdef int i
    for i in range(n):
        projected_seperations[i] = _getProjectedSeperation(nu[i], e, incl, w, radius_1)

def getProjectedSeperation(np.ndarray[double, ndim=1] nu, double e, double incl, double w, double radius_1) :
    cdef double[::1] projected_seperations= np.empty(nu.shape[0], dtype=nu.dtype)
    __getProjectedSeperation(&nu[0], e, incl, w, radius_1, nu.shape[0], &projected_seperations[0]) 
    return projected_seperations.base

cdef getProjectedPosition(double nu, double w, double incl) : 
    return sin(nu + w)*sin(incl) 


############################################################################################
#                               Utilities                                                  #
############################################################################################
cdef _find_secondary_phase(double fs, double fc):
    cdef double w = atan2(fs,fc)
    cdef double e = fs**2 + fc**2
    cdef double Phi = M_PI + 2*atan((e*cos(w)) / sqrt(1 - e**2))
    return (Phi - sin(Phi)) / (2*M_PI)

def find_secondary_phase(double fs, double fc): return _find_secondary_phase(fs, fc)


############################################################################################
#                               Radial velocity                                            #
############################################################################################

cdef _rv(double * time, double t_zero, double period, double K1, double fs, double fc, double V0, double incl, int n, double * rv_arr ):
    cdef int i
    cdef double w, e, nu
    for i in range(n):
        # Sort out decorr parameters
        w = atan2(fs, fc) 
        e = fs**2 + fc**2

        # Calculate the true anomaly
        nu = _t_ecl_to_peri(t_zero, e, w, incl, period, 0.2) # time prior to periastron passage
        nu = _get_Mean_Anomaly(time[i], nu, period, e) # Get the mean anomaly
        nu =  _get_Eccentric_Anomaly(nu, e) # Get the eccentric anomaly via analytical form
        nu = _get_True_anamoly(nu, e)

        # Calculate te RV model
        rv_arr[i] = K1*(e*cos(w) + cos(nu + w)) + V0



def rv(np.ndarray[double, ndim=1] time, double t_zero = 0., double period = 1., double K1 = 10., double fs=0., double fc=0., double V0=0., double incl = 90.):
    cdef double[::1] rv_arr= np.empty(time.shape[0], dtype=time.dtype)
    incl = M_PI*incl/180.
    _rv(&time[0], t_zero, period, K1, fs,fc, V0, incl, time.shape[0], &rv_arr[0])
    return rv_arr.base


############################################################################################
#                               Overlap geometry of circles                                #
############################################################################################
cdef area(double d, double x, double R):
    if (x <= (R - d)) :  return M_PI*x*x						
    elif (x >= (R + d)) : return M_PI*R*R	
    elif (d > (x+R)) : return 0.			
    else :  
        arg1 = (d*d + x*x - R*R)/(2.*d*x)
        arg2 = (d*d + R*R - x*x)/(2.*d*R)
        arg3 = max((-d + x + R)*(d + x - R)*(d - x + R)*(d + x + R), 0.)
        return x*x*acos(arg1) + R*R*acos(arg2) - 0.5*sqrt(arg3)


############################################################################################
#                               Secodnary eclipse                                          #
############################################################################################

cdef frac_secondary(double d, double x, double R):
    return 1 - area(d, x,R)/(M_PI*x**2)

############################################################################################
#                               qpower2                                                    #
############################################################################################
cdef clip(double a, double b, double c):
    if (a < b)      :  return b
    elif (a > c)    :  return c
    else            :  return a

cdef q1(double z, double p, double c, double a, double g, double I_0):
    cdef double zt = clip(abs(z), 0,1-p)
    cdef double s = 1-zt*zt
    cdef double c0 = (1-c+c*pow(s,g))
    cdef double c2 = 0.5*a*c*pow(s,(g-2))*((a-1)*zt*zt-1)
    return -I_0*M_PI*p*p*(c0 + 0.25*p*p*c2 - 0.125*a*c*p*p*pow(s,(g-1)))

cdef q2(double z, double p, double c, double a, double g, double I_0, double eps):
    cdef double zt = clip(abs(z), 1-p,1+p)
    cdef double d = clip((zt*zt - p*p + 1)/(2*zt),0,1)
    cdef double ra = 0.5*(zt-p+d)
    cdef double rb = 0.5*(1+d)
    cdef double sa = clip(1-ra*ra,eps,1)
    cdef double sb = clip(1-rb*rb,eps,1)
    cdef double q = clip((zt-d)/p,-1,1)
    cdef double w2 = p*p-(d-zt)*(d-zt)
    cdef double w = sqrt(clip(w2,eps,1))
    cdef double c0 = 1 - c + c*pow(sa,g)
    cdef double c1 = -a*c*ra*pow(sa,(g-1))
    cdef double c2 = 0.5*a*c*pow(sa,(g-2))*((a-1)*ra*ra-1)
    cdef double a0 = c0 + c1*(zt-ra) + c2*(zt-ra)*(zt-ra)
    cdef double a1 = c1+2*c2*(zt-ra)
    cdef double aq = acos(q)
    cdef double J1 =  (a0*(d-zt)-(2./3.)*a1*w2 + 0.25*c2*(d-zt)*(2.0*(d-zt)*(d-zt)-p*p))*w + (a0*p*p + 0.25*c2*pow(p,4))*aq 
    cdef double J2 = a*c*pow(sa,(g-1))*pow(p,4)*(0.125*aq + (1./12.)*q*(q*q-2.5)*sqrt(clip(1-q*q,0.0,1.0)) )
    cdef double d0 = 1 - c + c*pow(sb,g)
    cdef double d1 = -a*c*rb*pow(sb,(g-1))
    cdef double K1 = (d0-rb*d1)*acos(d) + ((rb*d+(2./3.)*(1-d*d))*d1 - d*d0)*sqrt(clip(1-d*d,0.0,1.0))
    cdef double K2 = (1/3)*c*a*pow(sb,(g+0.5))*(1-d)
    if (J1 > 1)   : J1 = 0
    cdef double FF =  I_0*(J1 - J2 + K1 - K2)
    #if (FF < 0.9) : FF=1.0
    return -FF


cdef Flux_drop_analytical_power_2(double z, double k, double c, double a, double eps):
    '''
    Calculate the analytical flux drop por the power-2 law.

    Parameters
    z : double
        Projected seperation of centers in units of stellar radii.
    k : double
        Ratio of the radii.
    c : double
        The first power-2 coefficient.
    a : double
        The second power-2 coefficient.
    f : double
        The flux from which to drop light from.
    eps : double
        Factor (1e-9)
    '''
    cdef double I_0 = (a+2)/(M_PI*(a-c*a+2))
    cdef double g = 0.5*a

    if (z < 1-k)           : return q1(z, k, c, a, g, I_0)
    elif (abs(z-1) < k)    : return q2(z, k, c, a, g, I_0, eps)
    else                   : return 0. 


############################################################################################
#                               lc                                                         #
############################################################################################



cdef __lc(double time, double flux, double flux_err, double t_zero, double period, double radius_1, double k, double fs, double fc, double incl, double ldc_1, double ldc_2, double SBR, double light_3, double jitter, double zp,  int switch):
    cdef double w, e, nu, f, F_TRANSIT, loglike=0.
    # Sort out decorr parameters
    w = atan2(fs, fc) 
    e = fs**2 + fc**2

    # Allocate the zeropoint 
    F_TRANSIT = 1.

    # Calculate the true anomaly
    nu = _t_ecl_to_peri(t_zero, e, w, incl, period, 0.2) # time prior to periastron passage
    nu = _get_Mean_Anomaly(time, nu, period, e) # Get the mean anomaly
    nu =  _get_Eccentric_Anomaly(nu, e) # Get the eccentric anomaly via analytical form
    nu = _get_True_anamoly(nu, e)

    # Calculate the projected position to see if primary or secondary
    f = getProjectedPosition(nu, w, incl) 

    # Calculate the projected seperation
    nu = _getProjectedSeperation(nu, e, incl, w, radius_1)

    # Now check to see if in-transit
    if (nu < (1.0+ k)):
        if f > 0:
            # If we are here, this is a primary eclipse
            F_TRANSIT += Flux_drop_analytical_power_2(nu, k, ldc_1, ldc_2, 1E-8) 

            # Now account for dilution of primary caused by SBR
            if SBR > 0 : F_TRANSIT = (F_TRANSIT - 1)*(1. - k*k*SBR) + 1.
        elif SBR>0.:
            F_TRANSIT -= SBR*k**2*(1-frac_secondary(nu, k ,1.))

    # Now account for third light
    if (light_3 > 0.0) :  F_TRANSIT = F_TRANSIT/(1. + light_3) + (1.-1.0/(1. + light_3)) # third light

    # Calculate the loglike
    if switch : return -0.5*( (zp*flux - F_TRANSIT)**2 / (flux_err**2 + jitter**2))
    else      : return F_TRANSIT

cdef _lc(double * time, double * flux, double * flux_err, double t_zero, double period, double radius_1, double k, double fs, double fc, double incl, double ldc_1, double ldc_2, double SBR, double light_3, double jitter, double zp, int n, double * lc_arr, int switch):
    cdef int i
    cdef double loglike = 0.
    if switch==1:
        for i in range(n):
            loglike += __lc(time[i], flux[i], flux_err[i], t_zero, period, radius_1, k, fs,fc, incl, ldc_1, ldc_2, SBR, light_3, jitter, zp, switch)
        return loglike
    else:
        for i in range(n):
            lc_arr[i] = __lc(time[i], flux[i], flux_err[i], t_zero, period, radius_1, k, fs,fc, incl, ldc_1, ldc_2, SBR, light_3, jitter, zp, switch)
    return loglike

def lc(np.ndarray[double, ndim=1] time = np.linspace(0,1,100),  np.ndarray[double, ndim=1] flux = np.array([-1], dtype = np.float64), np.ndarray[double, ndim=1] flux_err = np.array([-1], dtype = np.float64),
        double t_zero = 0., double period = 1., double radius_1 = 0.2, double k=0.2, 
        double fs=0., double fc=0., double incl = 90., 
        double ldc_1 = 0.8 , double ldc_2 = 0.8, double SBR=0., double light_3 = 0., double jitter=0., double zp=1.0):
    cdef double[::1] lc_arr= np.empty(time.shape[0], dtype=time.dtype)
    cdef double loglike
    cdef int switch = 1
    if flux[0]==-1 : switch = 0
    incl = M_PI*incl/180.
    loglike = _lc(&time[0], &flux[0], &flux_err[0], t_zero, period, radius_1, k, fs,fc, incl, ldc_1, ldc_2, SBR, light_3, jitter, zp, time.shape[0], &lc_arr[0], switch)
    if switch==0 : return lc_arr.base
    else         : return loglike


############################################################################################
#                               TLS                                                        #
############################################################################################
cdef transit_width(double r, double k, double b, double P):
    return P*asin(r*sqrt( ((1+k)**2-b**2) / (1-b**2*r**2) ))/M_PI



cdef median_loglike(double * flux, double * flux_err, double median, double jitter,  int n, double * loglike_line):
    for i in range(n):
        loglike_line[i] = -0.5*( (flux[i] - median)**2 / (flux_err[i]**2 + jitter**2))






cdef clipint(int a, int b, int c):
    if (a < b)      :  return b
    elif (a > c)    :  return c
    else            :  return a

cdef _tls(double * time, double * flux, double * flux_err, double * loglike_line, double * loglike_transit,
        double * t_zero, int * idxs, double period = 1., double radius_1 = 0.2, double k=0.2, 
        double fs=0., double fc=0., double incl = 90., 
        double ldc_1 = 0.8 , double ldc_2 = 0.8, double SBR=0., double light_3 = 0., double jitter=0., double zp=1.0,
        int n=1, int n2=1, int nhalfwidth=1, double halfwidth=1.):
        '''
        cdef i, j, idx
        for i in range(n):
            for j in range(clipint(idxs[i] - nhalfwidth, 0, n2), clipint(idxs[i] + nhalfwidth, 0, n2)):
                loglike_transit[i] -= loglike_line[j]
                loglike_transit[i] += __lc(time[j], flux[j], flux_err[j], t_zero[i], period, radius_1, k, fs,fc, incl, ldc_1, ldc_2, SBR, light_3, jitter, 1)
        '''

        cdef double bic_line, bic_transit 
        cdef int i,j
        for i in range(n):
            bic_line = 0.
            bic_transit = 0.
            for j in range(clipint(idxs[i] - nhalfwidth, 0, n2), clipint(idxs[i] + nhalfwidth, 0, n2)):
                if abs(time[j] - t_zero[i]) < halfwidth : 
                    bic_line += loglike_line[j]
                    bic_transit += __lc(time[j], flux[j], flux_err[j], t_zero[i], period, radius_1, k, fs,fc, incl, ldc_1, ldc_2, SBR, light_3, jitter, zp, 1)
                loglike_transit[i] = bic_transit - bic_line

def template_match(np.ndarray[double, ndim=1] time,  np.ndarray[double, ndim=1] flux, np.ndarray[double, ndim=1] flux_err,
        double period = 1., double radius_1 = 0.2, double k=0.2, 
        double fs=0., double fc=0., double incl = 90., 
        double ldc_1 = 0.8 , double ldc_2 = 0.8, double SBR=0., double light_3 = 0., double jitter=0., double zp=1.0, double nsplit = 30):
    '''
    Function to template match data to search for transit events.
    '''
    # Convert inclination
    incl = M_PI*incl/180.

    # Allocate the chi array
    cdef double[::1] loglike_line = np.empty(time.shape[0], dtype=np.double)

    # Calculate the median
    cdef median = zp #np.median((flux))

    # Now populate the chi arrray with the median
    median_loglike(&flux[0], &flux_err[0], median, jitter,  time.shape[0], &loglike_line[0])

    # get the transit width
    cdef b = cos(incl) / radius_1 
    cdef double width = transit_width(radius_1, k, b, period)
    cdef int nhalfwidth = int(width / np.median(np.gradient(time)) / 2) + 10

    # Create the time array to check 
    # To d : reduce this to only points we can use
    cdef double[::1] time_supersampled = np.arange(np.min(time) - width, np.max(time) + width, width / nsplit, dtype=np.double)
    #cdef double[::1] time_supersampled = np.arange(np.min(time) - width, np.max(time) + width, np.median(np.gradient(time)), dtype=np.double)

    cdef double[::1] loglike_transit = np.empty(time_supersampled.shape[0], dtype=np.double)
    cdef int[::1] idxs = np.interp(time_supersampled, time , np.arange(time.shape[0]) ).astype(np.int32)

    print('Transit width = {:.2f} days [{:.2f} hrs]'.format(width, 24*width))
    print('Gradiant is {:.4f} so idx_halfwidth is {:}'.format(np.median(np.gradient(time)), nhalfwidth))

    # Now iterate
    
    _tls(&time[0], &flux[0], &flux_err[0], &loglike_line[0], &loglike_transit[0],
        &time_supersampled[0], &idxs[0], period, radius_1, k, 
        fs, fc, incl, 
        ldc_1, ldc_2, SBR, light_3, jitter, zp,
        time_supersampled.shape[0], time.shape[0], nhalfwidth, width/2.)
    
    return time_supersampled, loglike_transit, loglike_line, width, nhalfwidth
    # Get the index half width
    #cdef int idx_halfwidth = int(np.ceil(width / 2. / np.median(np.gradient(time)) )) + 3





############################################################################################
#                               Utitlity functions                                         #
############################################################################################
def lc_bin(time, flux, bin_width):
    '''
    Function to bin the data into bins of a given width. time and bin_width 
    must have the same units
    '''

    edges = np.arange(np.min(time), np.max(time), bin_width)
    dig = np.digitize(time, edges)
    time_binned = (edges[1:] + edges[:-1]) / 2
    flux_binned = np.array([np.nan if len(flux[dig == i]) == 0 else flux[dig == i].mean() for i in range(1, len(edges))])
    err_binned = np.array([np.nan if len(flux[dig == i]) == 0 else sem(flux[dig == i]) for i in range(1, len(edges))])
    time_bin = time_binned[~np.isnan(err_binned)]
    err_bin = err_binned[~np.isnan(err_binned)]
    flux_bin = flux_binned[~np.isnan(err_binned)]   

    return time_bin, flux_bin, err_bin

def grad_split_array(x, dx_lim):
    '''
    Split array buy time gaps. Can be used to get individual nights in datasets.
    '''
    dx = np.gradient(x) 
    dx_thresh = np.sort(np.where(dx >= dx_lim)[0] +1)

    # Now check for consecutive integers and delte them
    delete_idxs = []
    i = 0 
    while i < (dx_thresh.shape[0]-1):
        if (dx_thresh[i]+1) == dx_thresh[i+1] :
            delete_idxs.append(i+1)
            i +=2
        else : i +=1
    dx_thresh = np.delete(dx_thresh, delete_idxs)

    # create the idx to split
    idx = np.arange(x.shape[0])
    return np.split(idx, dx_thresh)

def lc_mags_to_flux(mags, mags_err):
    """
    Take 2 arrays, light curve and errors
    and convert them from differential magnitudes
    back to relative fluxes
    Rolling back these equations:
        mags = - 2.5 * log10(flux)
        mag_err = (2.5/log(10))*(flux_err/flux)
    """
    flux = 10.**(mags / -2.5)
    flux_err = (mags_err/(2.5/np.log(10))) * flux
    return flux, flux_err

def lc_flux_to_mags(flux, flux_err):
    """
    Take 2 arrays, light curve and errors
    and convert them from differential magnitudes
    back to relative fluxes
    Applying these equations:
        mags = - 2.5 * log10(flux)
        mag_err = (2.5/log(10))*(flux_err/flux)
    """
    mags = -2.5*np.log10(flux)
    mags_err = (2.5/np.log(10))*(flux_err/flux)
    return mags, mags_err

def phase_times(times, epoch, period, phase_offset=0.0):
    """
    Take a list of times, an epoch, a period and
    convert the times into phase space.
    An optional offset can be supplied to shift phase 0
    """
    return (((times - epoch)/period)+phase_offset)%1