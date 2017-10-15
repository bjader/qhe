import numpy as np
import matplotlib.pyplot as plt
import pylab
import os

data_dir = os.path.join(os.path.dirname(__file__),'..','data','density_profile')


x1,x2,x3,x4,x5,x6,x7 = [],[],[],[],[],[],[]
y1,y2,y3,y4,y5,y6,y7 = [],[],[],[],[],[],[]

#Import analytical results
with open(os.path.join(data_dir, 'density_n2_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x1.append(xval)
        y1.append(yval)


with open(os.path.join(data_dir, 'density_n5_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x2.append(xval)
        y2.append(yval)

with open(os.path.join(data_dir, 'density_n10_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x3.append(xval)
        y3.append(yval)

with open(os.path.join(data_dir, 'density_n15_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x4.append(xval)
        y4.append(yval)

with open(os.path.join(data_dir, 'density_n20_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x5.append(xval)
        y5.append(yval)

with open(os.path.join(data_dir, 'density_n25_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x6.append(xval)
        y6.append(yval)

with open(os.path.join(data_dir, 'density_n30_width0.100000_L10.000000_10.000000m_m1.txt') ) as f:
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[1])
        yval = float(currentline[0])
        x7.append(xval)
        y7.append(yval)


plt.title("2-D particle density profile for the m=1 Laughlin state")
plt.xlabel("Radius")
plt.ylabel("2-D particle density")

plt.plot(x1,y1, c='b', label='2')
plt.plot(x2,y2, c='r', label='5')
plt.plot(x3,y3, c='g', label='10')
plt.plot(x4,y4, c='c', label='15')
plt.plot(x5,y5, c='m', label='20')
plt.plot(x6,y6, c='y', label='25')
plt.plot(x7,y7, c='k', label='30')


plt.legend(title = 'n = ', loc = 'upper right')

plt.xlim(0,7)
plt.ylim(0)

plt.show()
