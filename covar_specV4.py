import os
import sys # Debug
import warnings
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt

warnings.simplefilter('always', UserWarning)

def excess_variance(cnt_rate:np.ndarray, error:np.ndarray) -> float:
    """Calculates the excess covariance of a light curve

    Args:
        cnt_rate (numpy.ndarray): 1D Array of count rate
        error (numpy.ndarray): 1D Array of errors with length = len(cnt_rate)

    Returns:
        float: excess covariance
    """
    
    #Calculate excesss covariance
    excess_var:float = (np.sum((cnt_rate-np.nanmean(cnt_rate))*(cnt_rate-np.nanmean(cnt_rate)))/(len(cnt_rate)-1)) - np.nanmean(error*error)
    """if excess_var < 0:
        excess_var = 0"""
            
    return excess_var

def excess_variance_2d(cnt_rate_2d, error_2d):
    excess_var = np.zeros(len(cnt_rate_2d))
    
    for i,_ in enumerate(cnt_rate_2d):
        excess_var[i] = excess_variance(cnt_rate_2d[i], error_2d[i])
    
    return excess_var

def covariance(cnt_rate:np.ndarray, ref_cnt_rate:np.ndarray) -> float:
    """Calculates covariance with respect to reference energy band

    Args:
        cnt_rate (numpy.ndarray): 1D array of countrate
        ref_cnt_rate (numpy.ndarray): 1D array of reference countrate

    Returns:
        float: Covariance
    """
    
    #Calculate covariance with respect to reference energy band
    cov:float = np.sum((cnt_rate-np.nanmean(cnt_rate))*(ref_cnt_rate-np.nanmean(ref_cnt_rate))) / (len(cnt_rate)-1)
    
    return cov

def covariance_2d(cnt_rate_2d, ref_cnt_rate_2d):
    cov_2d = np.zeros(len(cnt_rate_2d))
    for i,_ in enumerate(cnt_rate_2d):
        cov_2d[i] = covariance(cnt_rate_2d[i], ref_cnt_rate_2d[i])
        
    return cov_2d

def normalized_covariance(cov:float, ref_cnt_rate: np.ndarray, ref_error:np.ndarray) -> float:
    """Calculates normalized covariance

    Args:
        cov (float): Covariance of lightvurve with respect to reference enery band
        ref_cnt_rate (numpy.ndarray): 1D array of reference countrate
        ref_error (numpy.ndarray): 1D array of reference lightcurve errors

    Returns:
        float: Normalized covariance
    """
    
    #Calculate excess covariance of reference band
    var_excess = excess_variance(ref_cnt_rate, ref_error)
    #Calculate normalized covariance
    cov_norm = cov/np.sqrt(var_excess)
    
    return cov_norm

def normalized_covariance_2d(cov_2d, ref_cnt_rate_2d, ref_err_2d):
    var_excess = excess_variance_2d(ref_cnt_rate_2d, ref_err_2d)
    cov_norm = cov_2d/np.sqrt(var_excess)
    
    return cov_norm, var_excess

def avg_cov(cnt_rate_2d:np.ndarray, ref_cnt_rate_2d:np.ndarray, ref_err_2d:np.ndarray) -> float:
    """Calculates average covariance from all the segments

    Args:
        cnt_rate_2d (np.ndarray): 2D array of count rate divided into segments
        ref_cnt_rate_2d (np.ndarray): 2D array of reference count rate divided into M segments

    Returns:
        float: Average covariance
    """
    #vec_cov = np.vectorize(covariance)
    covar = covariance_2d(cnt_rate_2d, ref_cnt_rate_2d)
    return np.nanmean(normalized_covariance_2d(covar, ref_cnt_rate_2d, ref_err_2d))

def avg_excess_var(cnt_rate_2d:np.ndarray, err_2d:np.ndarray) -> float:
    """Calculates average excess covariance

    Args:
        cnt_rate_2d (np.ndarray): 2D array of count rates
        err_2d (np.ndarray): 2D array of errors

    Returns:
        float: Average excess covariance
    """
    evar = excess_variance_2d(cnt_rate_2d, err_2d)
    evar = evar[evar>0]

    return np.nanmean(evar)

def normalized_covariance_error(cnt_rate_2d:np.ndarray, ref_cnt_rate_2d:np.ndarray, err_2d:np.ndarray, ref_err_2d:np.ndarray) -> float:
    """Calculates normalized covariance

    Args:
        cnt_rate_2d (np.ndarray): 2D array of count rate bronken into M segments
        ref_cnt_rate_2d (np.ndarray): 2D array of reference count rates broken into M segments
        err_2d (np.ndarray): 2D array pf errors containing M segments
        ref_err_2d (np.ndarray): 2D aray of reference band errors containing M segments

    Returns:
        float: Normalized covariance error
    """
    
    evar_x_avg = avg_excess_var(cnt_rate_2d, err_2d) # Average excess covariance of countrate
    evar_y_avg = avg_excess_var(ref_cnt_rate_2d, ref_err_2d) # Average p.excess covariance of reference countrate
    err:np.ndarray = err_2d.flatten()
    ref_err:np.ndarray = ref_err_2d.flatten()
    err2:np.ndarray = err*err
    ref_err2:np.ndarray = ref_err*ref_err
    
    segments:int = len(cnt_rate_2d)
    bins_ps:int = len(cnt_rate_2d[0])
    
    nc_err2:float = (evar_x_avg*np.nanmean(ref_err2) +
               evar_y_avg*np.nanmean(err2) +
               np.nanmean(err2)*np.nanmean(ref_err2)) / (bins_ps*segments*evar_y_avg)
                
    nc_err:float = np.sqrt(nc_err2)
                
    return nc_err
    

def rms_err(n_bins, err, cnt_rate, excess_var):
    #print(err)
    #err2 = np.nanmean((err-np.nanmean(err))*(err-np.nanmean(err)))
    err2 = np.nanmean(err*err)
    cnt_rate_avg2 = np.nanmean(cnt_rate)*np.nanmean(cnt_rate)
    F_var = np.sqrt(excess_var/cnt_rate_avg2)
    p1 = np.sqrt(2/n_bins)*(err2/cnt_rate_avg2)
    p2 = np.sqrt(err2/n_bins)*2*(F_var/np.nanmean(cnt_rate))
    #p1 = np.sqrt(2/n_bins)*(err2)
    #p2 = np.sqrt(err2/n_bins)*2*(F_var)
    
    err = np.sqrt(p1*p1 + p2*p2)
    return err

def rms_err_2d(n_bins, err_2d, cnt_rate_2d, excess_var):
    rms_error = np.zeros(len(err_2d))
    for i,_ in enumerate(err_2d):
        rms_error[i] = rms_err(n_bins, err_2d[i], cnt_rate_2d[i], excess_var[i])
        
    return rms_error

def read_fits(fits_file_path, clip_arr=True):
    hdul = fits.open(fits_file_path)
    data = hdul[1].data
    data = data.view(np.recarray)  # Structured numpy record array
    data = np.array(data.tolist())
    if clip_arr:
        data = data[2:len(data)-2]
    hdul.close()
    
    time, count_rate, err = data.T
    time-=time[0]
    
    #plt.plot(time[:35], count_rate[:35], "k.")
    #plt.ylim(0,130)
    #plt.show()
    
    return np.array([count_rate, err, time]).T

def user_def_input():
    path = input("Enter the absolute path to the directory that contains your data without quotes (Default: Current dir): ")
    if not path.strip():
        path = os.path.dirname(os.path.abspath(__file__))
    
    if not os.path.isdir(path):
        raise Exception("PathNotFound: Check directory path")
    
    ref_lcurve = input("Enter the name of the reference light curve with extension (Default: en-sub-ref.fits): ")
    if not ref_lcurve.strip():
        ref_lcurve = "en-sub-ref.fits"
    if not os.path.exists(f"{path}/{ref_lcurve}"):
        raise Exception("FileNotFound: Check reference light curve file name")
    
    _,_,time = read_fits(f"{path}/{ref_lcurve}").T
    delta_t, total_t = time[1], time[len(time)-1]
    f_max, f_min = 0.5/delta_t, 1/total_t
    print(f"Bin size is {delta_t} seconds")
    #print(f"Maximum allowed frequency: {f_max}")
    print(f"Minimum allowed frequency: {np.format_float_scientific(f_min, precision=4)}")
    
    while(True):
        f_low = float(input("Enter low frequency (Lower bound): "))
        if not isinstance(f_low, float):
            raise Exception("IncorrectDataType: Frequency must be float or int")
        if f_low < f_min:
            warning = f"Lower frequency must be >= {f_min}"
            warnings.warn(warning)
        else:
            print("Accepted")
            print(f"The RMS and covariance will be calculated for the Fourier frequency range {f_low} Hz and {f_max} Hz")
            break
        
    return f_low, total_t, ref_lcurve, path

def bin_time(f_max:float) -> float:
    binning_time = 0.5/f_max
    
    return binning_time

def segment_time(f_min:float) -> float:
    seg_time = 1/f_min
    
    return seg_time

def calc_segments_nbins(seg_time:float, total_time:float, arr:np.ndarray) -> np.ndarray:
    length = len(arr)
    segments = int(np.floor(total_time/seg_time))
    n_bins = int(np.floor(length/segments)) # Number of bins per segment
    new_length = int(n_bins*segments)
    
    return segments, n_bins, new_length

def reshape_array(segments, n_bins, arr):
    new_length = int(n_bins*segments)
    arr = arr[:new_length]
    arr = arr.reshape((int(segments), int(n_bins)))
    
    return arr

def calc_covar_spectra():
    f_low, total_t, ref_lcurve, path = user_def_input()
    seg_time = segment_time(f_low)
    reference_lc = read_fits(f"{path}/{ref_lcurve}", clip_arr=True)
    _,_,ref_time = reference_lc.T
    #print(total_t, seg_time)
    segments, n_bins, new_length = calc_segments_nbins(seg_time, total_t, ref_time)
    print("Total Number of Segments:", segments)
    ref_cnt_rate, ref_err, ref_time = reference_lc[:new_length].T
    ref_cnt_rate = reshape_array(segments, n_bins, ref_cnt_rate)
    
    ref_err = reshape_array(segments, n_bins, ref_err)
    ref_time = reshape_array(segments, n_bins, ref_time)
    covar_arr = np.zeros((13,segments))
    cov_err_arr = np.zeros((13,segments))
    excess_var_arr = np.zeros((13,segments))
    rms_err_arr = np.zeros((13,segments))
    excess_var = excess_variance_2d(ref_cnt_rate, ref_err)
    evar_check = np.sum([excess_var<0])
    _ = print(f"WARNING: {evar_check} negative excess variances detected for reference energy band") if evar_check != 0 else None
    for enr in range(13):
        light_curve = read_fits(f"{path}/en-sub{enr+1}.fits")
        cnt_rate, err, time = light_curve[:new_length].T
        #print("length", len(cnt_rate))
        cnt_rate = reshape_array(segments, n_bins, cnt_rate)
        err = reshape_array(segments, n_bins, err)
        time = reshape_array(segments, n_bins, time)
        cov = covariance_2d(cnt_rate, ref_cnt_rate)
        norm_cov,_ = normalized_covariance_2d(cov, ref_cnt_rate, ref_err)
        
        #norm_cov,_ = normalized_covariance_2d(cov, cnt_rate, err)
        excess_var = excess_variance_2d(cnt_rate, err)
        evar_check = np.sum([excess_var<0])
        _ = print(f"WARNING: {evar_check} negative excess variances detected for energy band number {enr+1}") if evar_check != 0 else None
        excess_var_arr[enr] = excess_var
        covar_err = normalized_covariance_error(cnt_rate, ref_cnt_rate, err, ref_err)
        covar_arr[enr] = norm_cov
        cov_err_arr[enr] = covar_err
        rms_err_arr[enr] = rms_err_2d(n_bins, err, cnt_rate, excess_var)
    
    return covar_arr, excess_var_arr, cov_err_arr, rms_err_arr

cov, excess_var, cov_err, rms_error = calc_covar_spectra()
#print(rms_error)
cov = np.nanmean(cov, axis=1)
excess_var = np.nanmean(excess_var, axis=1)
cov_err = np.nanmean(cov_err, axis=1)
rms_error = np.nanmean(rms_error, axis=1)

        
###### PLOTTING RMS AND COVARIANCE SPECTRA ########################     
    

energy = np.array([0.4, 0.6, 0.8, 1.0, 1.25, 1.55, 1.85, 2.5, 2.75, 3.5, 4.5, 6.0, 8.0])
energy_err = np.array([0.1, 0.1, 0.1, 0.1, 0.15, 0.15, 0.15, 0.25, 0.25, 0.5, 0.5, 1.0, 1.0])
plt.errorbar(energy, cov/energy, xerr=energy_err, yerr=cov_err, linestyle="", fmt=".", label="Covariance")
#plt.errorbar(energy, np.sqrt(excess_var)/energy, xerr=energy_err, yerr=rms_error, linestyle="", fmt=".", label="rms")
#plt.errorbar(energy, rms_error, linestyle="", fmt=".", label="rms")
plt.yscale("log")
plt.xscale("log")
#plt.ylim(8e-2, 1e1)
plt.legend()
plt.xlabel("Energy (keV)")
plt.ylabel(r"Normalized counts s$^{-1}$ keV$^{-1}$")
#plt.show()
plt.savefig('covar-rms.png')
np.savetxt("norm_cov.txt", np.array([energy, energy_err, cov/energy, cov_err]).T, fmt=('%2.2f', '%1.2f', '%2.4f', '%1.4f'))
np.savetxt("rms.txt", np.array([energy, energy_err, np.sqrt(excess_var)/energy, rms_error]).T, fmt=('%2.2f', '%1.2f', '%2.4f', '%1.4f'))
plt.show()
