import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sgn

from raspi_import import raspi_import
import data

def prepros (data, maximum = 4096):
    #normalize to 0 - 1
    d = data / maximum
    try:
        #remove the average, so that the signal goes from -1 to 1
        d = sgn.detrend(d) 
    except ZeroDivisionError:
        d = d - np.mean(d)
    #d = d - np.mean(d)
    return d

def windowedcorr(x,y, maxlag):
    rxy = np.zeros(2*maxlag)
    for l in range(0, 2*maxlag):
        rxy[l] = np.sum(np.roll(x, maxlag - l)*y) 
    
    return rxy

def findCorrPeaks(data, plot=False):
    #interpolate
    resample_coeff = 1
    #data = filter(data) #np.fft.fft(signal)
    data = sgn.resample(data, len(data)*resample_coeff)
    
    maxlag = min(1000, len(data))
    number_of_peaks = 5
    corr = windowedcorr(data,data,maxlag)
    
    #corr = sgn.hilbert(np.abs(corr), len(corr)//2)

    abs_corr = np.abs(corr)[((len(corr)//2)):len(corr)]
    fft = np.abs(np.fft.fft(abs_corr))[0:maxlag//2]
    fft[0] = 0 #remove 0 hz

    max_fft_peak = np.where((fft == np.amax(fft)))[0][0]

    max_peak = 1000 / max_fft_peak

    max_peaks = sgn.find_peaks(abs_corr, prominence=0.1)[0]

    max_peaks = max_peaks[0: min(number_of_peaks, len(max_peaks) ) ]
    if(plot):
        print(max_peaks)
        print(max_fft_peak, max_peak)
        plt.plot(abs_corr, "-d", markevery=max_peaks.tolist())
        plt.show()
        plt.plot(fft)
        plt.show()
    
    return max_peaks/resample_coeff, max_peak/resample_coeff


def findVelocity(data, plot=False):

    q = prepros(data[0][:,0])
    i =  prepros(data[0][:,1])

    fft_i = np.abs(np.fft.fft(i)[:len(i)//2])
    fft_q = np.abs(np.fft.fft(q)[:len(i)//2])
    
    direction = 1 if ( np.mean(np.arctan2(fft_q, fft_i)) > 0 ) else -1

    peaks, peak  = findCorrPeaks(i, plot)
    if(len(peaks) > 0):
        arr = np.zeros(len(peaks) - 1)
        for n in range(len(peaks)-1):
            arr[n] = peaks[n + 1] - peaks[n]
        meanlag = arr.mean()
        fn =  31250/(meanlag*2)

        f_0 = 24.13e9
        c = 3e8
        v_corr = ( c*(fn) )/ (2*(f_0))
    else:
        v_corr = 0
        print("Could not compute velocity by correlation")

    v_corr_fft = ( c*(31250/(peak*2)) )/ (2*(f_0))


    m = np.amax(fft_i)
    i = np.where(fft_i == m)[0][0]
    
    fn = ((i / (len(fft_i)*2))) * 31250

    f_0 = 24.13e9
    c = 3e8
    v = ( c*(fn) )/ (2*(f_0))
    print("Measurement ", data[1], " Actual velocity ", data[2])
    print("Estimated velocity: ", v,",", v_corr, ",", v_corr_fft)
    print("Error ", v-data[2], ",", v_corr-data[2], ",", v_corr_fft-data[2])

    return str(direction*v) + "," + str(direction*v_corr) + "," + str(direction*v_corr_fft)

def findVelocities(data):
    vs = [0]*len(data)
    for i in range(0, len(data)):
        vs[i] = findVelocity(data[i])
    return vs

if __name__ == "__main__":
    for i in findVelocities(data.getRadarDatas(21,30)):
        print(i)

