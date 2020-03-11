import matplotlib.pyplot as plt
import numpy as np

a=[-30.3314216559651,19.0004016698518,-17.150793787408,9.49499574218739,-2.54768404538229,0.265382965410969]
b=[-180.992524120965,168.471004362887,-67.499549702687,13.5075841245848,-1.31983368963974,0.0500087685129987]
c=3.0
d=21.2968837223113


temps=np.logspace(1,6,50)
heats=[]
for temp in temps:
        heating=0.0
        if (temp >1.0e5):
            heating=(c*np.log10(temp))-d            
        elif temp> 891.0:
            for i,B in enumerate(b):
                heating=heating+(B*(np.log10(temp)**(i)))
        else:
            for i,A in enumerate(a):
                heating=heating+(A*(np.log10(temp)**(i)))
        heats.append(heating)

fig,ax=plt.subplots()
ax.plot(temps,heat
	s)
ax.set(xscale='log')
plt.show()