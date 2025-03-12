#!/usr/bin/env python
# coding: utf-8

# In[2]:


import uproot
import pandas
import numpy
import matplotlib.pyplot as plt


# #### Opening the 10M file

# In[ ]:


input_file = uproot.open("10M.root")
input_tree = input_file["CaloExample"]["Calo"]
# Works only in uproot 3
#df = tree.pandas.df(flatten=False)

# Works only un uproot 4
np_array = input_tree.arrays(library="np")
np_title = np_array.keys()
df = pandas.DataFrame()


# In[ ]:


input_tree.keys()


# #### Opening the paraFile that contains zero and neg chi sq which indicate problem crystals
file = uproot.open("paraFile.root")
feature = file["covar"]
np_array_paraFile = feature.arrays(library='pd')
np_title_paraFile = np_array_paraFile.keys()
#df_paraFile = pandas.DataFrame()
#print(np_title_paraFile)
#file = uproot.open("paraFile.root")
#feature = file["covar"]
#df = feature.arrays(library='pd')
### Using for loop on 10M file info to create a crystal map
# #### Opening the arvivTable for chi2 that contains zero peak error which indicate problem crystals for chi2 model
file_chi2 = uproot.open("arXivTable_chi2.root")
feature = file_chi2["covar"]
np_array_chi2 = feature.arrays(library='pd')
np_title_chi2 = np_array_chi2.keys()
# #### Opening the arvivTable for chi2 that contains zero peak error which indicate problem crystals for nll model
file_nll = uproot.open("arXivTable_nll.root")
feature = file_nll["covar"]
np_array_nll = feature.arrays(library='pd')
np_title_nll = np_array_nll.keys()
# #### Opening the arvivTable for chi2 that contains zero peak error which indicate problem crystals for nll model

# In[ ]:


file_nll = uproot.open("arXivTable_disk0_migrad.root")
feature = file_nll["covar"]
np_array_nll = feature.arrays(library='pd')
np_title_nll = np_array_nll.keys()


# In[7]:


for title in np_title:
    dt = np_array[title].dtype.descr
    for idx in range(len(dt)):
        subtitle = dt[idx][0]
        if subtitle == "":
            df_title = title
            np_data = np_array[title]
            # convert 32 bytes data into 64 bytes for minuit migrad to work
            if type(np_data[0]) == numpy.float32:
                np_data = np_data.astype(numpy.float64)
            if type(np_data[0]) == numpy.int32:
                np_data = np_data.astype(numpy.int64)
            df[df_title] = np_data
        else:
            dd_title = title + " " + subtitle
            np_data = np_array[title][subtitle]
            # convert 32 bytes data into 64 bytes for minuit migrad to work
            if type(np_data[0]) == numpy.float32:
                np_data = np_data.astype(numpy.float64)
            if type(np_data[0]) == numpy.int32:
                np_data = np_data.astype(numpy.int64)


print(np_title)
for i in range(5):
    #print("nCluster: ", df["nCluster"][i], "\n")
    #print("cluList: ", df["cluList"][i], "\n")
    #print("length of cluList: ", len(df["cluList"][i]), "\n")
    print("cryId: ", df["cryId"][i])
    print("cryPosX: ", df["cryPosX"][i])
    print("cryPosY: ", df["cryPosY"][i])
    print("cryPosZ: ", df["cryPosZ"][i])
    print("crySimIdx: ", df["crySimIdx"][i])


max_cryId = 1347
# Initialize cryPos with an extra column for cryId (-1000 for placeholders)
cryPos = -1000 * numpy.ones((max_cryId, 4))

i = 0
# Loop until all placeholders in cryPos are updated
while numpy.any(cryPos[:, 1:] == [-1000, -1000, -1000]):
    for j in range(len(df["cryId"][i])):
        cry_id = df["cryId"][i][j]  # Get the cryId
        cryPos[cry_id - 1] = [
            cry_id,  # First column stores the cryId
            df["cryPosX"][i][j],
            df["cryPosY"][i][j],
            df["cryPosZ"][i][j],
        ]
    i += 1
numpy.savetxt("cryPos_wID.txt", cryPos, fmt = ['%d', '%.4e', '%.4e', '%.4e'], delimiter=',')#fmt='%.4e'fmt = ['%d', '%.4e', '%.4e', '%.4e']


# In[8]:


for title in np_title:
    dt = np_array[title].dtype.descr
    for idx in range(len(dt)):
        subtitle = dt[idx][0]
        if subtitle == "":
            df_title = title
            np_data = np_array[title]
            # convert 32 bytes data into 64 bytes for minuit migrad to work
            if type(np_data[0]) == numpy.float32:
                np_data = np_data.astype(numpy.float64)
            if type(np_data[0]) == numpy.int32:
                np_data = np_data.astype(numpy.int64)
            df[df_title] = np_data
        else:
            dd_title = title + " " + subtitle
            np_data = np_array[title][subtitle]
            # convert 32 bytes data into 64 bytes for minuit migrad to work
            if type(np_data[0]) == numpy.float32:
                np_data = np_data.astype(numpy.float64)
            if type(np_data[0]) == numpy.int32:
                np_data = np_data.astype(numpy.int64)


print(np_title)
for i in range(5):
    #print("nCluster: ", df["nCluster"][i], "\n")
    #print("cluList: ", df["cluList"][i], "\n")
    #print("length of cluList: ", len(df["cluList"][i]), "\n")
    print("cryId: ", df["cryId"][i])
    print("cryPosX: ", df["cryPosX"][i])
    print("cryPosY: ", df["cryPosY"][i])
    print("cryPosZ: ", df["cryPosZ"][i])
    print("crySimIdx: ", df["crySimIdx"][i])


max_cryId = 1347
# Initialize cryPos with an extra column for cryId (-1000 for placeholders)
cryPos = -1000 * numpy.ones((max_cryId, 7))

i = 0
# Loop until all placeholders in cryPos are updated
while numpy.any(cryPos[:, 1:] == [-1000, -1000, -1000]):
    for j in range(len(df["cryId"][i])):
        cry_id = df["cryId"][i][j]  # Get the cryId
        cryPos[cry_id - 1] = [
            cry_id,  # First column stores the cryId
            df["cryPosX"][i][j],
            df["cryPosY"][i][j],
            df["cryPosZ"][i][j],
            df["simStartX"][i][j],
            df["simStartY"][i][j],
            df["simStartZ"][i][j],
        ]
    i += 1
numpy.savetxt("testfile.txt", cryPos, fmt = ['%d', '%.4e', '%.4e', '%.4e','%.4e', '%.4e', '%.4e'], delimiter=',')#fmt='%.4e'fmt = ['%d', '%.4e', '%.4e', '%.4e']


# In[9]:


import numpy as np

max_cryId = 1347
# Initialize pipePos with placeholders
pipePos = -1000 * np.ones((max_cryId, 4))

i = 0
# Loop until all placeholders in pipePos are updated
while np.any(pipePos[:, 1:] == [-1000, -1000, -1000]):
    if i >= len(df["cryId"]):  # Prevent index out of range
        print("No more data to process. Check your input files.")
        break

    for j in range(len(df["cryId"][i])):
        try:
            cry_id = df["cryId"][i][j]  # Get the cryId
            pipePos[cry_id - 1] = [
                cry_id,  # First column stores the cryId
                df["simStartX"][i][j],
                df["simStartY"][i][j],
                df["simStartZ"][i][j],
            ]
        except IndexError:
            print(f"Index error at i={i}, j={j}. Check input data consistency.")
            continue

    i += 1

# Save the array
np.savetxt("pipePosition.txt", pipePos, fmt=['%d', '%.4e', '%.4e', '%.4e'], delimiter=',')


# #### Using a for loop to find locations of problem crystals

# In[15]:


zero_chisq_with_crystal = []
zero_chisq = []
# Iterate over each value in df["ChiSq"] and the corresponding crystalNo in df["crystalNo"]
for chi_sq, crystal_no in zip(np_array_paraFile["ChiSq"], np_array_paraFile["crystalNo"]):
    # Check if the ChiSq value is equal to 0
    if chi_sq == 0:
        # Append a tuple of (crystal_no, chi_sq) to zero_chisq_with_crystal
        zero_chisq_with_crystal.append((crystal_no, chi_sq))
        
for crystal_no, _ in zero_chisq_with_crystal:
    zero_chisq.append(crystal_no)
print(zero_chisq)


# In[16]:


neg_chisq_with_crystal = []
neg_chisq = []
for chi_sq, crystal_no in zip(np_array_paraFile["ChiSq"], np_array_paraFile["crystalNo"]):
    if chi_sq < 0:
        neg_chisq_with_crystal.append((crystal_no, chi_sq))
        
for crystal_no, _ in neg_chisq_with_crystal:
    neg_chisq.append(crystal_no)
print(neg_chisq)


# #### Using a for loop to find locations of zero peak errors in chi2 model

# In[1]:


zero_error_with_crystal = []
zero_error = []
# Iterate over each value in df["ChiSq"] and the corresponding crystalNo in df["crystalNo"]
for peakerr, crystal_no in zip(np_array_chi2["PeakErr"], np_array_chi2["crystalNo"]):
    # Check if the ChiSq value is equal to 0
    if peakerr == 0:
        # Append a tuple of (crystal_no, chi_sq) to zero_chisq_with_crystal
        zero_error_with_crystal.append((crystal_no, peakerr))
        
for crystal_no, _ in zero_error_with_crystal:
    zero_error.append(crystal_no)


# #### Using a for loop to find locations of zero peak errors in nll model
zero_error_nll_with_crystal = []
zero_error_nll = []
# Iterate over each value in df["ChiSq"] and the corresponding crystalNo in df["crystalNo"]
for peakerr_nll, crystal_no in zip(np_array_nll["PeakErr"], np_array_nll["crystalNo"]):
    # Check if the ChiSq value is equal to 0
    if peakerr_nll == 0:
        # Append a tuple of (crystal_no, chi_sq) to zero_chisq_with_crystal
        zero_error_nll_with_crystal.append((crystal_no, peakerr_nll))
        
for crystal_no, _ in zero_error_nll_with_crystal:
    zero_error_nll.append(crystal_no)
# In[17]:





# In[25]:


import numpy as np
import matplotlib.pyplot as plt
import uproot

# Load cryPos data
cryPos = np.loadtxt("cryPos_wID.txt", delimiter=',')
cryPos = cryPos[673:]  # Use only rows from the 674th row onward

# Number of crystals
num_cry = len(cryPos)

# Create figure and axes
fig = plt.figure(figsize=(30,30))
ax = fig.add_subplot(1, 1, 1)
ax.set_aspect('equal')

# Draw rectangles with colors
for i in range(num_cry):
    x_val = cryPos[i, 1]  # x-coordinate
    y_val = cryPos[i, 2]  # y-coordinate
    cry_id = int(cryPos[i, 0])  # Cry ID from the first column

    # Check if cry_id is in zero_chisq
    if cry_id in zero_chisq:
        face_color = 'red'  # Color for zero_chisq crystals
    elif cry_id in neg_chisq:
            face_color = 'yellow'  # Color for zero_chisq crystals
    else:
        face_color = 'white'  # Default color

    # Draw the rectangle
    ax.fill_between(
        [x_val - 17, x_val + 17],
        y_val - 17, y_val + 17,
        color=face_color,
        edgecolor='black'
    )

    # Add text for cry_id
    ax.text(x_val, y_val, str(cry_id), va='center', ha='center', fontsize=16, color='blue',fontweight='bold')
# Scatter overlay with pipePosition data
pipe_x = pipePos[:, 1]  # x-coordinates from pipePosition
pipe_y = pipePos[:, 2]  # y-coordinates from pipePosition

# Scatter plot for pipe positions
ax.scatter(pipe_x, pipe_y, color='green', s=100, label="Pipe Positions", alpha=0.7)

# Set limits, grid, and ticks
min_x, max_x = np.min(cryPos[:, 1]) - 20, np.max(cryPos[:, 1]) + 20
min_y, max_y = np.min(cryPos[:, 2]) - 20, np.max(cryPos[:, 2]) + 20
ax.set_xlim(min_x, max_x)
ax.set_ylim(min_y, max_y)
ax.set_xticks(np.arange(min_x, max_x, 38))
ax.set_yticks(np.arange(min_y, max_y, 38))
ax.grid(color='gray', linestyle='--', linewidth=0.5)
# Adjust axis number font sizes
ax.tick_params(axis='both', which='major', labelsize=20)  # Major ticks
ax.tick_params(axis='both', which='minor', labelsize=10)  # Minor ticks (optional)
# Display the plot
plt.savefig("Problem_crys.png",dpi=150, bbox_inches='tight')  # Save as PNG with high resolution
plt.show()


# ##### making crystal map to indicate which crystals have zero peak

# In[24]:


import numpy as np
import matplotlib.pyplot as plt

# Load cryPos data
cryPos = np.loadtxt("cryPos_wID.txt", delimiter=',')
cryPos = cryPos[673:]  # Use only rows from the 674th row onward

# Number of crystals
num_cry = len(cryPos)

# Create figure and axes
fig = plt.figure(figsize=(21, 13))
ax = fig.add_subplot(1, 1, 1)
ax.set_aspect('equal')

# Draw rectangles with colors
for i in range(num_cry):
    x_val = cryPos[i, 1]  # x-coordinate
    y_val = cryPos[i, 2]  # y-coordinate
    cry_id = int(cryPos[i, 0])  # Cry ID from the first column

    # Determine face color based on dataset membership
    if cry_id in zero_error and cry_id in zero_error_nll:
        face_color = 'green'  # Overlapping crystals are green
    elif cry_id in zero_error:
        face_color = 'red'  # Crystals in zero_error only
    elif cry_id in zero_error_nll:
        face_color = 'yellow'  # Crystals in zero_error_nll only
    else:
        face_color = 'white'  # Default color

    # Draw the rectangle
    ax.fill_between(
        [x_val - 17, x_val + 17],
        y_val - 17, y_val + 17,
        color=face_color,
        edgecolor='black'
    )

    # Add text for cry_id
    ax.text(x_val, y_val, str(cry_id), va='center', ha='center', fontsize=8, color='blue')

# Set limits, grid, and ticks
min_x, max_x = np.min(cryPos[:, 1]) - 20, np.max(cryPos[:, 1]) + 20
min_y, max_y = np.min(cryPos[:, 2]) - 20, np.max(cryPos[:, 2]) + 20
ax.set_xlim(min_x, max_x)
ax.set_ylim(min_y, max_y)
ax.set_xticks(np.arange(min_x, max_x, 34))
ax.set_yticks(np.arange(min_y, max_y, 34))
ax.grid(color='gray', linestyle='--', linewidth=0.5)

# Display the plot
plt.savefig("Problem_crys.png", dpi=150, bbox_inches='tight')  # Save as PNG with high resolution
plt.show()


# In[ ]:




