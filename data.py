import numpy as np
from raspi_import import raspi_import
from matplotlib import pyplot as plt
from scipy import signal as sgn
#v=1,08 for data 61-65
#v=0,77 for data 51-55
#v=1,50 for data 41-45
#v=0,26 for data 31-35
#v=0,69 for data 24-28

roiDict = {
        1:(70000,190000),2:(40000,190000),3:(90000,190000),4:(60000,150000),
        5:(90000,140000),6:(60000,130000),7:(40000,100000),8:(110000,240000),
        9:(90000,120000),10:(90000,130000),


        31:(170000,192000),32:(163800,170000),33:(140000,180000),34:(130000,150000),
           41:(30000,42700),42:(30000,45000),43:(70000,80000),44:(40000,50000),
           51:(58700,80000),52:(60000,70800),53:(70000,80000),54:(60000,70000),
           61:(50000,60000),62:(60000,70000),63:(60000,70000),64:(60000,70000)
}

velocities = {
    1:0.53,2:0.71,3:0.88,4:1.04,5:0.45,6:1.23,7:1.66,8:0.77,9:0.32,10:1.35,
    21:0.76,22:0.83,23:0.77,24:0.85,25:0.89,26:1.92,27:1.83,28:1.85,29:0.54
    
}


sample_period = 32

def prepros (data, maximum = 4096):
    d = sgn.detrend(data,axis=0)
    ba = sgn.butter(2, [0.0016], btype='high')
    d = sgn.lfilter(ba[0], ba[1], d, axis=0)
    return d
    
        

def notch_filter (data):
    #b1 = [0.99,-1.831,0.99,]
    #a1 = [1,-1.831,0.9801]
    #data = sgn.lfilter(b1,a1, data, axis=0)

    #ba = sgn.butter(3,[0.01312*2, 0.01440*2],btype='bandstop')
    #data = sgn.lfilter(ba[0],ba[1], data, axis=0)

    return data


def plotDatas(data):
    for d in data:
        plt.title("Measurement " + str(d[1]) + " Velocity " + str(d[2]) + " m/s")
        plt.plot(d[0][:,1])
        plt.show()

def plotData(data):
    d = data
    plt.title("Measurement " + str(d[1]) + " Velocity " + str(d[2]) + " m/s")
    plt.plot(d[0][:,1])
    plt.show()
def plotFFTDatas(data):
    f_0 = 24.13e9
    c = 3e8
    for d in data:
        
        plt.title("Measurement " + str(d[1]) + " Velocity " + str(d[2]) + " m/s")
        f = d[2]*((2*f_0))/(c)
        plt.axvline(x=f, ymin=0, ymax=50, label='Doppler frequency', color='g')
        fft = (np.abs(np.fft.fft(d[0][:,1])[:-2]))
        plt.xlim(-3000,3000)
        plt.plot(np.fft.fftshift(np.fft.fftfreq(len(fft),sample_period/(1e6))),np.fft.fftshift(fft))
        plt.show()
def plotFFTData(data):
    f_0 = 24.13e9
    c = 3e8
    d = data
    plt.title("Measurement " + str(d[1]) + " Velocity " + str(d[2]) + " m/s")
    f = d[2]*((2*f_0))/(c)
    plt.axvline(x=f, ymin=0, ymax=50, label='Doppler frequency', color='g')
    fft = (np.abs(np.fft.fft(d[0][:,0])[:-2]))
    plt.xlim(-3000,3000)
    plt.plot(np.fft.fftshift(np.fft.fftfreq(len(fft),sample_period/(1e6))),np.fft.fftshift(fft))
    plt.show()

def getRadarDatas(start, end):
    datas = [0]*(end - start)
    for i in range(start,end):
        try:
            s = roiDict[i][0]
            e = roiDict[i][1]
        except KeyError:
            print("Roi not found")
            s = 0
            e = -1
        try:
            velocity = velocities[i]
        except KeyError:
            raise "Velocity not found"

        data1 = (raspi_import("./out/radar" + str(i) + ".bin"))[1][s:e,3:5]
        data1 = prepros(data1)
        data1 = notch_filter(data1)
        #data1 = filter(data1)
        #print(data1)

        datas[i - start] = [data1, i, velocity]
    return datas

def getRadarData(num):
    s = 0
    e = -1
    try:
        s = roiDict[num][0]
        e = roiDict[num][1]
    except KeyError:
        print("Roi not found")
        s = 0
        e = -1
    velocity = 0
    try:
        velocity = velocities[num]
    except KeyError:
        raise Exception("Velocity not found")
    
    
    data1 = (raspi_import("./out/radar" + str(num) + ".bin"))[1][s:e,3:5]
    data1 = prepros(data1)
    #data1 = notch_filter(data1)
    #data1 = filter(data1)

    return [data1, num, velocity]
    


if __name__ == "__main__":

    #data1 = [getRadarDatas(24,28), 0.69]
    data2 = getRadarData(1)
    data3 = getRadarData(2)
    data4 = getRadarData(3)
    data5 = getRadarData(4)

    #plotData(data2[0])
    #plotData(data3)
    plotFFTDatas(getRadarDatas(21,30))
    #plotData(data4)
    #plotData(data5)

    #print(getRadarDatas(31,35))