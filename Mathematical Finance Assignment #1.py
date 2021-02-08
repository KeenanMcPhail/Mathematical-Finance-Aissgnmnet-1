#!/usr/bin/env python
# coding: utf-8

# In[208]:


import numpy as np
import math
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.interpolate import interp1d


# In[209]:


#Bond Data
bd=[["Maturity Date", "Issue Date", "ISIN", "Coupon", "Jan 18", "Jan 19", "Jan 20", "Jan 21", "Jan 22"," Jan 25", "Jan 26", "Jan 27", "Jan 28", "Jan 29"],
    [3/1/2021,10/19/2015,"CA135087F254",0.750,100.08,100.08,100.08,100.07,100.07,100.06,100.06,100.06,100.06,100.06],
          [9/1/2021,4/11/2016,"CA135087F585",0.750,100.40,100.40,100.39,100.38,100.38,100.37,100.37,100.38,100.37,100.37],
          [3/1/2022,10/11/2016,"CA135087G328",0.500,100.44,100.43,100.41,100.40,100.41,100.39,100.39,100.41,100.41,100.41],
          [8/1/2022,5/4/2020,"CA135087L286",0.250,100.20,100.19,100.16,100.16,100.16,100.14,100.15,100.17,100.18,100.18],
          [3/1/2023,10/6/2017,"CA135087H490",1.750,103.29,103.29,103.25,103.24,103.25,103.22,103.22,103.25,103.26,103.24],
          [3/1/2024,3/1/2021,"CA135087J546",2.250,106.21,106.21,106.17,106.15,106.17,106.12,106.18,106.22,106.19,106.18],
          [9/1/2024,4/5/2019,"CA135087J967",1.500,104.27,104.28,104.25,104.19,104.24,104.29,104.26,104.32,104.24,104.24],
          [3/1/2025,10/11/2019,"CA135087K528",1.250,103.60,103.60,103.59,103.51,103.52,103.62,103.59,103.66,103.61,103.57],
           [9/1/2025,4/3/2020,"CA135087K940",0.500,103.33,103.32,103.30,103.22,103.27,103.37,103.33,103.42,103.35,103.33]]

#ytd at maturity in days approximated assuming each month has thirty days
md=[61,241,421,571,781,1141,1321,1501,1681]
# closing prices ytd
cpd=[18,19,20,21,22,25,26,27,28,29]
           


# In[221]:


x_0=1




ytm=[[],[],[],[],[],[],[],[],[],[]]
for j in range(10):
    for i in range(1,10):
        def f(x):
            expytm=[(100+(bd[i][3]/2))*math.exp(-x*((md[i-1]-cpd[j])/365))]
            for v in range(1,50):
                if md[i-1]-180*v>cpd[j]:
                    expytm.append((bd[i][3]/2)*math.exp(-x*((md[i-1]-180*v-cpd[j])/365)))
            y=sum(expytm)
            return bd[i][j+4]-y
        sol=optimize.root(f,x_0,method='lm')
        ytm[j].append(sol.x[0])

print(ytm)


# In[222]:


for i in range(len(ytm)):
             plt.plot(md,ytm[i],label=str(cpd[i]))

plt.title("Yield Curve of aaa Candian Government Bonds")
plt.xlabel("days since January 1st 2021")
plt.ylabel("yield ")

plt.legend()
plt.show()


# In[213]:


spot=[[],[],[],[],[],[],[],[],[],[]]
for j in range(10):
    for i in range(1,10):
        expytm=[math.exp(-spot[j][k]*((md[k]-cpd[j])/365))for k in range(0,len(spot[j]))]
        y=sum(expytm)
        x= bd[i][4+j]-(bd[i][3]/2)*y
        spot[j].append(-math.log(x/(100+(bd[i][3]/2))/((md[i-1]-cpd[j]/365))))




print(spot)


# In[145]:


for i in range(len(spot)):
             plt.plot(md,spot[i],label=str(cpd[i]))
plt.plot(md, spot[0],"r.")


plt.title("Spot Curve of aaa Candian Government Bonds")
plt.xlabel("days since January 2st 2021")
plt.ylabel("yield ")

plt.legend()
plt.show()


# In[151]:


g=interp1d(md,spot[0])
print(g(100))


# In[233]:


fwd=[[],[],[],[],[],[],[],[],[],[]]
for j in range(10):
    g=interp1d(md,spot[j])
    for v in range(0,8):
        fwd[j].append(g(365+cpd[j]+180*v)*(365+cpd[j]+180*v)/365-g(365+cpd[j]))

fwdx=[]


print(fwd)


# In[234]:


for j in range(len(fwd)):
    x=[365+cpd[j]+180*v for v in range(0,8)]
    plt.plot(x,fwd[j],label=str(cpd[j]))
    

plt.title("1yr-tyr Foward curve of aaa Candian Government Bonds")
plt.xlabel("days since January 2st 2021")
plt.ylabel("Forward rate")

plt.legend()
plt.show()


# In[225]:


MAT1=[[],[],[],[]]
#b = np.linspace(0, 2000, num=2000, endpoint=True)
for i in range(0,4):
    for j in range(1,9):
        g=interp1d(md,ytm[j])
        h=interp1d(md,ytm[j+1])
        MAT1[i].append(math.log(h((i+1)*365+cpd[j])/g((i+1)*365+cpd[j])))
print(MAT1)
cov1=np.cov(MAT1,bias=True)
print(cov1)


# In[238]:


MAT2=[[],[],[],[]]
#b = np.linspace(0, 2000, num=2000, endpoint=True)
for i in range(0,4):
    for j in range(1,9):
        x=[365+cpd[j]+180*v for v in range(0,8)]
        g=interp1d(x,fwd[j])
        h=interp1d(x,fwd[j+1])
        MAT2[i].append(math.log(h((i+1)*365+cpd[j])/g((i+1)*365+cpd[j])))
print(MAT2)
cov2=np.cov(MAT2,bias=True)
print(cov2)


# In[244]:


np.linalg.eigh(cov1)


# In[242]:


np.linalg.eigh(cov2)


# In[ ]:




