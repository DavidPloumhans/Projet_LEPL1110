import numpy as np
import cv2
import sys
import os

files = os.listdir('images')  # liste des fichiers images


case = 0.2 # facteur entre sigma_fluid et sigma_bottom
meshname = str(case)  # nom du mesh dont on fait la vidéo
# Create a VideoWriter object
size = 650 * 2
fps = 24
output = cv2.VideoWriter(f'output_{case}.avi', cv2.VideoWriter_fourcc(*'XVID'), fps, (size, size))

# find the sigmas
sigmas = []
for file in files:
    if file.endswith('.png'):
        sigma_bottom = file.split('_')[1] # la contrainte en bas
        if sigma_bottom not in sigmas:
            sigmas.append(float(sigma_bottom))

sigmas.sort()
for sigma in sigmas:
    images = [f"{meshname}_{sigma}_def.png", f"{meshname}_{sigma}_VM.png", f"{meshname}_{sigma}_ZZ.png", f"{meshname}_{sigma}_CM.png"]
    # Concatenate the images
    images = [cv2.imread(os.path.join('images', image)) for image in images]
    # concatenate the images 2x2
    imgTop = np.concatenate((images[0], images[1]), axis=1)
    imgBottom = np.concatenate((images[2], images[3]), axis=1)
    img = np.concatenate((imgTop, imgBottom), axis=0)
    # Write the sigma
    sigmaMPa = float(sigma) / 1e6
    cv2.putText(img, f"Cas ou Sigma fluid = {case} * Sigma Bottom", (30, 30), cv2.FONT_HERSHEY_SIMPLEX, 0.7, (255, 0, 0), 3)
    cv2.putText(img, f"Sigma bottom: {int(sigmaMPa)} MPa", (int(1.3 * size/4), int(size/4)), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 4)
    cv2.putText(img, f"Sigma fluid: {int(sigmaMPa * case)} MPa", (int(1.3 * size/4), int(1.1 * size/4)), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 255, 0), 4)
    # Display the image
    # cv2.imshow('image', img)
    # cv2.waitKey(1000)
    # Write the image to the video
    output.write(img)

output.release()