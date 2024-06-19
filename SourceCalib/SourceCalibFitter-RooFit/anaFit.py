import uproot
import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
file = uproot.open("paraFile.root")
Pions = file["newTree;1"]
df = Pions.pandas.df(flatten=False)

sns.pairplot(df, x_vars = ["FullPeak","FullWidth","ChiSquared"], y_vars = ["FullPeak","FullWidth","ChiSquared"])
plt.savefig("pairs.eps")
plt.show()

resolution = []
for i,j in enumerate(df["FullPeak"]):
    resolution.append(100*df["FullWidth"][i]/df["FullPeak"][i])

#plt.hist2d(resolution, df["FullPeak"], cmap='Purples',bins=100)

fig, ax = plt.subplots(1,1)
n, bins, patches = ax.hist(resolution,
                           bins=40,
                           range=(0,12),histtype='step',
                           label="Full Fraction")
fig, ax = plt.subplots(1,1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
plt.errorbar(bin_centers, n, yerr=np.sqrt(n), fmt='b.')
emean = round(np.mean(resolution),2)
estd = round(np.std(resolution),2)

plt.text(2,50,'Mean:'+str(emean))
plt.text(2,45, 'Std dev: '+str(estd))
ax.set_ylabel('Crystals per bin')
ax.set_xlabel("Resolution [%]")
fig.savefig('Fullresolution.pdf')


resolution = []
for i,j in enumerate(df["FirstPeak"]):
    resolution.append(100*df["FirstWidth"][i]/df["FirstPeak"][i])

#plt.hist2d(resolution, df["FullPeak"], cmap='Purples',bins=100)

fig, ax = plt.subplots(1,1)
n, bins, patches = ax.hist(resolution,
                           bins=40,
                           range=(0,12),histtype='step',
                           label="Full Fraction")
fig, ax = plt.subplots(1,1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
plt.errorbar(bin_centers, n, yerr=np.sqrt(n), fmt='b.')
emean = round(np.mean(resolution),2)
estd = round(np.std(resolution),2)

plt.text(2,50,'Mean:'+str(emean))
plt.text(2,45, 'Std dev: '+str(estd))
ax.set_ylabel('Crystals per bin')
ax.set_xlabel("Resolution [%]")
fig.savefig('Firstresolution.pdf')


resolution = []
for i,j in enumerate(df["SecondPeak"]):
    resolution.append(100*df["SecondWidth"][i]/df["SecondPeak"][i])

#plt.hist2d(resolution, df["FullPeak"], cmap='Purples',bins=100)

fig, ax = plt.subplots(1,1)
n, bins, patches = ax.hist(resolution,
                           bins=40,
                           range=(0,20),histtype='step',
                           label="Full Fraction")
fig, ax = plt.subplots(1,1)
bin_centers = 0.5 * (bins[:-1] + bins[1:])
plt.errorbar(bin_centers, n, yerr=np.sqrt(n), fmt='b.')
emean = round(np.mean(resolution),2)
estd = round(np.std(resolution),2)

plt.text(2,80,'Mean:'+str(emean))
plt.text(2,75, 'Std dev: '+str(estd))
ax.set_ylabel('Crystals per bin')
ax.set_xlabel("Resolution [%]")
fig.savefig('Secondresolution.pdf')


fig, ax = plt.subplots(1,1)
n, bins, patches = ax.hist(df["vfrFull"],
                           bins=100,
                           range=(0,1),histtype='step',
                           label="Full Fraction")
n, bins, patches = ax.hist(df["vfrFirst"],
                           bins=100,
                           range=(0,1),histtype='step',
                           label="First Escape Fraction")
n, bins, patches = ax.hist(df["vfrSecond"],
                           bins=100,
                           range=(0,1),histtype='step',
                           label="Second Escape Fraction")
ax.set_ylabel('N')
ax.legend()
ax.set_xlabel("Fraction of Events")
fig.savefig('fraction.pdf')

