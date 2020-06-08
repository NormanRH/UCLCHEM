import pandas as pd
import matplotlib.pyplot as plt
plt.style.use("thesis")

pca_df=pd.read_csv("Dimensionality/principle_components.csv")

fig,ax=plt.subplots(figsize=(4.27,4.27))
ax.plot(pca_df["Component"],pca_df["Cumulative"])
ax.set(xlabel="Number of Components",ylabel="Explained Variance")
fig.savefig("paperplots/pca.pdf",type="PDF",bbox_inches="tight")