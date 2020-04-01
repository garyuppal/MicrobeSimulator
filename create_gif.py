# -*- coding: utf-8 -*-
"""
Created on Mon 22 Jan 2020

@author: guppa
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import sys
from celluloid import Camera

dataDir = sys.argv[1] 
dataDir = dataDir.rstrip('/')

path = os.path.dirname(os.path.realpath(__file__)) + "/" + dataDir + "/"

if len(sys.argv) == 3:
    video_name = sys.argv[2]
    if not video_name.endswith(".gif"):
        video_name = video_name + ".gif"
else:
    video_name = "video_" + dataDir + "_v.gif"

# plot_folder =
folder = os.fsencode(path)

# filenames = []

fig, ax = plt.subplots(1)
camera = Camera(fig)

# get boundary data:
with open(path + "boundary.dat", "r") as f:
    boundary_data = np.loadtxt(f)

# get rectangles:
with open(path + "rectangles.dat", "r") as f:
    rectangles = np.loadtxt(f)
# rectangle[i,:] = [x1, y1, x2, y2]


def rescaleColor(sc, smax):
    temp = sc/smax
    if temp > 1.0:
        return 1.0
    else:
        return temp


def plotAndSnap(data, boundary, rect):
    if data.size > 0:
        x = data[:, 0]
        y = data[:, 1]
        # plt.figure()
        srates = data[:, 2]
        myColors = [[0, rescaleColor(sr, 100.0), 0] for sr in srates]
        plt.scatter(x, y, s=10, c=myColors)
        plt.xlim(boundary[0, 0], boundary[1, 0])
        plt.ylim(boundary[0, 1], boundary[1, 1])
     
    ax.set_aspect('equal')

    # print(rectangles.shape)

    if rect.size !=0 and rect.ndim == 1:
        rect = np.reshape(rect,(1,4))

    for i in range(rect.shape[0]):
        rectpatch = patches.Rectangle((rect[i, 0], rect[i, 1]), width=(rect[i, 2]-rect[i, 0]),
                                 height=(rect[i, 3]-rect[i, 1]), linewidth=1, edgecolor='black', facecolor='grey')
        ax.add_patch(rectpatch)
    # plt.show()
    camera.snap()


# get bacteria data and plot:
for file in os.listdir(folder):
    filename = os.fsdecode(file)
    if filename.startswith("bacteria"):
        with open(path + filename, "r") as f:
            dataB = np.loadtxt(f, ndmin=2)
            plotAndSnap(dataB, boundary_data, rectangles)

# save animation:
print("saving animation...")
animation = camera.animate()
animation.save(video_name, writer='imagemagick')
print("done")
