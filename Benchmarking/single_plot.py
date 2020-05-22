import pandas as pd
import matplotlib.pyplot as plt

fig,ax=plt.subplots()
df=pd.read_csv("test/low_rad_full.dat")
df=df.rename(lambda x: x.strip(),axis=1)
species=["CO","H2O","#CO","#H2O"]

ax.set(xscale="log",yscale="log")

for spec in species:
	ax.plot(df["Time"],df[spec],label=spec)
ax.legend()
ax.set(xscale="log",yscale="log",ylim=(1e-12,1e-3))

plt.show()
