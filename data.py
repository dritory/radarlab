import numpy as np
from raspi_import import raspi_import
from matplotlib import pyplot as plt
from scipy import signal as sgn
#v=1,08 for data 61-65
#v=0,77 for data 51-55
#v=1,50 for data 41-45
#v=0,26 for data 31-35
#v=0,69 for data 24-28

roiDict = {31:(170000,192000),32:(163800,170000),33:(140000,180000),34:(130000,150000),
           41:(30000,42700),42:(30000,45000),43:(70000,80000),44:(40000,50000),
           51:(58700,80000),52:(60000,70800),53:(70000,80000),54:(60000,70000),
           61:(50000,60000),62:(60000,70000),63:(60000,70000),64:(60000,70000)
}


sample_period = 32

def prepros (data, maximum = 4096):
    d = sgn.detrend(data,axis=0)

    return d
    
        

def notch_filter (data):
    b1 = [0.99,-1.831,0.99,]
    a1 = [1,-1.831,0.9801]
    data = sgn.lfilter(b1,a1, data, axis=0)

    #ba = sgn.butter(3,[0.01312*2, 0.01440*2],btype='bandstop')
    #data = sgn.lfilter(ba[0],ba[1], data, axis=0)

    return data


def plotData(data):
    for d in data:
        plt.title("Measurement " + str(d[1]) + " Velocity " + str(d[2]) + " m/s")
        plt.plot(d[0][:,1])
        plt.show()
def plotFFTData(data):
    f_0 = 24.13e9
    c = 3e8
    for d in data:
        
        plt.title("Measurement " + str(d[1]) + " Velocity " + str(d[2]) + " m/s")
        f = d[2]*((4*f_0))/(c)
        plt.axvline(x=f, ymin=0, ymax=50, label='Doppler frequency')
        fft = (np.abs(np.fft.fft(d[0][:,1])[:-2]))
        plt.xlim(-3000,3000)
        plt.plot(np.fft.fftshift(np.fft.fftfreq(len(fft),sample_period/(1e6))),np.fft.fftshift(fft))
        plt.show()

def getRadarDatas(start, end, velocity):
    datas = [0]*(end - start)
    for i in range(start,end):
        s = roiDict[i][0]
        e = roiDict[i][1]
        data1 = (raspi_import("./out/radarData" + str(i) + ".bin"))[1][s:e,3:5]
        data1 = prepros(data1)
        data1 = notch_filter(data1)
        #data1 = filter(data1)
        #print(data1)

        datas[i - start] = [data1, i, velocity]
    return datas

def getRadarData(num, velocity):
    s = roiDict[i][0]
    e = roiDict[i][1]
    data1 = (raspi_import("./out/radarData" + str(i) + ".bin"))[1][s:e,3:5]
    data1 = prepros(data1)
    #data1 = notch_filter(data1)
    #data1 = filter(data1)
    #print(data1)

    return [data1, i, velocity]
    


#data1 = [getRadarDatas(24,28), 0.69]
data2 = getRadarData(1, 0.53)
data3 = getRadarData(2, 0.71)
data4 = getRadarData(3, 0.88)
data5 = getRadarData(4, 1.04)

#plotData(data2[0])
#plotData(data3)
plotFFTData(data3)
#plotData(data4)
#plotData(data5)

#print(getRadarDatas(31,35))