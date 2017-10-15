import numpy as np
import matplotlib.pyplot as plt
import os

data_dir = os.path.join(os.path.dirname(__file__),'..','data','s2')

#Set up empty arrays to fill with data
x1,x2,x3,x4,x5 = [],[],[],[],[]
y1,y2,y3,y4,y5 = [],[],[],[],[]
y_err1,y_err2,y_err3,y_err4,y_err5 = [],[],[],[],[]

#m3 data
with open(os.path.join(data_dir, 'MC_n30_10m_L10_m1.txt') ) as f:
    
    #Skip first n lines of data
    n=0
        
    for x in range(0,n):
        next(f)

    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[2])
        yval = float(currentline[0])
        y_err = float(currentline[1])
        
        x1.append(np.sqrt(xval))
        y1.append(yval)
        y_err1.append(y_err)

plt.title("$S_2$ particle scaling for $m=1$ Laughlin states")
plt.xlabel(r"Number of particles ($\sqrt{n}$)")
plt.ylabel("Entanglement entropy ($S_{2}$)")

#Manually enforce axes not to be negative

#Plot scatter graphs
plt.errorbar(x1,y1, yerr=y_err1, c='b', ls='none', marker='x', label="m=1 10m iterations")

plt.legend(loc = 'upper left')

plt.show()