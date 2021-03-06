import numpy as np
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
from scipy.stats import linregress

m=3
gamma=0.5*np.log(m)

data_dir = os.path.join(os.path.dirname(__file__),'..','data','s2')

#Set up empty arrays to fill with data
x1,x2,x3,x4,x5 = [],[],[],[],[]
y1,y2,y3,y4,y5 = [],[],[],[],[]
y_err1,y_err2,y_err3,y_err4,y_err5 = [],[],[],[],[]

#m3 data
with open(os.path.join(data_dir, 'MC_m3_10m_L20_n30.txt') ) as f:
#with open(os.path.join(data_dir, 'MC_n23_100m_L10_m3.txt') ) as f:
    
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


#Normalise the intercept of our data to the analytic true value
y2 = np.array(y2) + gamma

#Find least square fit for each scatter
model = sm.OLS(y2,x2)
result = model.fit()
slope = result.params[0]

#Create an x-array for the regression line
x_max = x2[len(x2)-1]
x_min = x2[0]
xplot = np.linspace(x_min,x_max,num=100)

#Apply the straight line equation with the parameters from linregress
y1_line = (slope * xplot)

#Plot the linear regression lines
plt.plot(xplot,y1_line, c='black', ls='dotted', label= 'y=' + "%.3f" % slope + 'x')

plt.title("$S_2$ particle scaling for $m=3$ Laughlin states")
plt.xlabel(r"Number of particles ($\sqrt{n}$)")
plt.ylabel("Entanglement entropy ($S_{2}$ + " + "%.3f" % gamma + ")")

#Plot scatter graphs
plt.errorbar(x2,y2, yerr=y_err2, c='r', ls='none', marker='x', label="L=20")

plt.legend(loc = 'upper left')


plt.show()
