#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# # glob, os, and requests are used to check paths and retrieve remote data in this notebook. you do 
# # not need to import them for most uses of pdr.
# import glob
# import os 
# import numpy as np

# import requests
# import pandas as pd
# import pdr

# url = "https://search-pdsppi.igpp.ucla.edu/ditdos/download?id=pds://PPI/mess-mag-calibrated/data/mso-avg/2015"
# filename = url.split('/')[-1]

# print(filename)
# print("Checking for file")

# if not os.path.exists(filename):
#     print("Getting file")
#     req = requests.get(url)
#     #Doesn't work for folder
#     open(filename, 'wb').write(req.content)
#     print(len(req.content))
# else:
#     print("File found")


# Test_Data = pd.read_table('2015/001_031_JAN/MAGMSOSCIAVG15001_01_V08.TAB')
# print(np.shape(Test_Data))


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Create a figure
fig = plt.figure(figsize=(10, 15))

# Create a GridSpec with 5 rows and 3 columns
gs = gridspec.GridSpec(5, 3)

# Create the 4 vertically aligned subplots (spanning all 3 columns each)
ax1 = fig.add_subplot(gs[0, :])
ax2 = fig.add_subplot(gs[1, :])
ax3 = fig.add_subplot(gs[2, :])
ax4 = fig.add_subplot(gs[3, :])

# Create the 3 horizontally aligned subplots (spanning only 1 column each)
ax5 = fig.add_subplot(gs[4, 0])
ax6 = fig.add_subplot(gs[4, 1])
ax7 = fig.add_subplot(gs[4, 2])

# Example data to plot
x = range(10)
y = [i**2 for i in x]

# Plot example data in the subplots
for ax in [ax1, ax2, ax3, ax4, ax5, ax6, ax7]:
    ax.plot(x, y)

# Set titles to identify the subplots
ax1.set_title('Vertical Subplot 1')
ax2.set_title('Vertical Subplot 2')
ax3.set_title('Vertical Subplot 3')
ax4.set_title('Vertical Subplot 4')
ax5.set_title('Horizontal Subplot 1')
ax6.set_title('Horizontal Subplot 2')
ax7.set_title('Horizontal Subplot 3')

# Adjust layout
plt.tight_layout()
plt.show()