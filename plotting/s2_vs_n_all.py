import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import statsmodels.api as sm
import os

data_dir = os.path.join(os.path.dirname(__file__),'..','data','s2')

#Set up empty arrays to fill with data
x1,x2,x3,x4,x5 = [],[],[],[],[]
y1,y2,y3,y4,y5 = [],[],[],[],[]
y_err1,y_err2,y_err3,y_err4,y_err5 = [],[],[],[],[]

#m1 data
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

#m3 data
with open(os.path.join(data_dir, 'MC_n30_10m_L10_m3.txt') ) as f:
    
    #Skip first n lines of data
    n=0
        
    for x in range(0,n):
        next(f)

    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[2])
        yval = float(currentline[0])
        y_err = float(currentline[1])
        
        x2.append(np.sqrt(xval))
        y2.append(yval)
        y_err2.append(y_err)

with open(os.path.join(data_dir, 'MC_n23_100m_L10_m3.txt') ) as f:
    
    #Skip first n lines of data
    n=0
    
    for x in range(0,n):
        next(f)
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[2])
        yval = float(currentline[0])
        y_err = float(currentline[1])
        
        x3.append(np.sqrt(xval))
        y3.append(yval)
        y_err3.append(y_err)

#m5 data
with open(os.path.join(data_dir, 'MC_n28_100m_L10_m5.txt') ) as f:
    
    #Skip first n lines of data
    n=0
    
    for x in range(0,n):
        next(f)
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[2])
        yval = float(currentline[0])
        y_err = float(currentline[1])
        
        x4.append(np.sqrt(xval))
        y4.append(yval)
        y_err4.append(y_err)

with open(os.path.join(data_dir, 'MC_n30_10m_L10_m5.txt') ) as f:
    
    #Skip first n lines of data
    n=0
    
    for x in range(0,n):
        next(f)
    
    for line in f:
        
        currentline = line.split(",")
        xval = float(currentline[2])
        yval = float(currentline[0])
        y_err = float(currentline[1])
        
        x5.append(np.sqrt(xval))
        y5.append(yval)
        y_err5.append(y_err)



plt.title("$S_2$ particle scaling for Laughlin states")
plt.xlabel(r"Number of particles ($\sqrt{n}$)")
plt.ylabel("Entanglement entropy ($S_{2}$)")

#Manually enforce axes not to be negative

#Plot scatter graphs
plt.errorbar(x1,y1, yerr=y_err1, c='b', ls='none', marker='x', label="m=1 10m")
plt.errorbar(x2,y2, yerr=y_err2, c='r', ls='none', marker='x', label="m=3 10m")
#plt.errorbar(x3,y3, yerr=y_err3, c='g', ls='none', marker='x', label="m=3 100m")
#plt.errorbar(x4,y4, yerr=y_err4, c='c', ls='none', marker='x', label="m=5 10m")
plt.errorbar(x5,y5, yerr=y_err5, c='m', ls='none', marker='x', label="m=5 100m")


'''#Plot least square fit for each scatter
slope,intercept,r_value,p_value,std_error = linregress(x1_scaled,y1_scaled)

model = sm.OLS(y1_scaled,sm.add_constant(x1_scaled))
result = model.fit()
print(result.summary())


#Create an x-array for the regression line
x_max = x1_scaled[len(x1_scaled)-1]
x_min = x1_scaled[0]
xplot = np.linspace(x_min,x_max,num=100)

#Apply the straight line equation with the parameters from linregress
y1_line = (slope * xplot + intercept)

#Plot the linear regression lines
plt.plot(xplot,y1_line, c='black', ls='dotted', label= 'y=' + "%.3f" % slope + 'x + ' + "%.3f" % intercept)
'''
plt.legend(loc = 'upper left')

plt.show()
