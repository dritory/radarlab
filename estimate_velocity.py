import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sgn

from raspi_import import raspi_import

import data

def prepros (data, maximum = 4096):
    d = data / maximum
    #d = d - 0.5
    try:
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
    number_of_peaks = 10
    corr = windowedcorr(data,data,maxlag)

    abs_corr = np.abs(corr)[((len(corr)//2)):len(corr)]
    
    max_peaks = sgn.find_peaks(abs_corr, prominence=0.05)[0]

    max_peaks = max_peaks[0: min(number_of_peaks, len(max_peaks) ) ]
    if(plot):
        plt.plot(abs_corr, "-d", markevery=max_peaks.tolist())
        plt.show()
    
    return max_peaks/resample_coeff


def findVelocity(data, start= 0, stop = -1):
    q = prepros(data[start:stop,0])
    i =  prepros(data[start:stop,1])

    fft_i = np.abs(np.fft.fft(i)[:len(i)//2])
    fft_q = np.abs(np.fft.fft(q)[:len(i)//2])
    
    direction = 1 if ( np.mean(np.arctan2(fft_q, fft_i)) > 0 ) else -1

    peaks  = findCorrPeaks(i)
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
    m = np.amax(fft_i)
    i = np.where(fft_i == m)[0][0]
    
    fn = ((i / (len(fft_i)*2))) * 31250

    f_0 = 24.13e9
    c = 3e8
    v = ( c*(fn) )/ (2*(f_0))

    return direction*v, direction*v_corr


data2 = data.data2
print (findVelocity(data2[0][0]))


if(False):
    sample_period, data1 = raspi_import("./raw_data/adcData_1.bin")
    sample_period, data2 = raspi_import("./raw_data/adcData_2.bin")
    sample_period, data3 = raspi_import("./raw_data/adcData_3.bin")
    sample_period, data4 = raspi_import("./radarData/radarData31.bin")

    data = data4
    win_len = 10000
    length = data.shape[0]//win_len
    v1 = np.zeros(length)
    v2 = np.zeros(length)
    num = 200
    for i in np.linspace(0,length - 1,num):
        #print(i*win_len, (i + 1)*win_len)
        v1[int(i)], v2[int(i)] = findVelocity(data, int(i*win_len), int((i + 1)*win_len))
    plt.plot(v1)
    plt.plot(v2)
    plt.show()