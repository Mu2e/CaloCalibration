#!/usr/bin/env python
# coding: utf-8

import uproot
import pandas
import numpy as np
import matplotlib.pyplot as plt


# Opening the 10M file
input_file = uproot.open("10M.root")
input_tree = input_file["CaloExample"]["Calo"]
np_array = input_tree.arrays(library="np")
np_title = np_array.keys()
df = pandas.DataFrame()

# Opening the paraFile that contains zero and neg chi sq which indicate problem crystals
file = uproot.open("paraFile.root")
feature = file["covar"]
np_array_paraFile = feature.arrays(library='pd')
np_title_paraFile = np_array_paraFile.keys()

#Opening the arvivTable for chi2 that contains zero peak error which indicate problem crystals for nll model
file_dik0 = uproot.open("arXivTable_disk0_migrad.root")
feature_disk0 = file_disk0["covar"]
np_array_disk0 = feature_disk0.arrays(library='pd')
np_title_disk0 = np_array_disk0.keys()


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

#create pipe positions as simstart points
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


# Use a for loop to find locations of problem crystals

fullfrac_with_crystal = []
fullfrac_zero = []
# Iterate over each value in df["-param-"] and the corresponding crystalNo in df["crystalNo"]
for frfull, crystal_no in zip(np_array_paraFile["frFull"], np_array_paraFile["crystalNo"]):
    # Check if the ChiSq value is equal to 0
    if frfull == 0:
        # Append a tuple of (crystal_no, -param'-) to -param-_with_crystal
        fullfrac_with_crystal.append((crystal_no, frfull))
        
for crystal_no, _ in fullfrac_with_crystal:
    fullfrac_zero.append(crystal_no)


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

    # Add crystals of interest to look over
    if cry_id in fullfrac_zero:
        face_color = 'red'  # Color for problem param
    #elif cry_id in neg_chisq:
     #       face_color = 'yellow'  # Color for other probme param
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

