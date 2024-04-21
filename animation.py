import numpy as np
import cv2
import sys
import os
import keyboard as kb
import threading
import time
import pyautogui

def take_screenshot(img_name):
    image = pyautogui.screenshot(region = (0, 0, 650, 650))
    image.save(img_name)

def capture_images(mesh, sigma_bottom, initial_sleep, interval=2):
    time.sleep(initial_sleep)
    path = "C:/Users/Boss/Desktop/LEPL1110/Projet/images"
    take_screenshot(path + f"/{mesh}_{sigma_bottom}_def.png")  # deformation
    time.sleep(interval)
    take_screenshot(path + f"/{mesh}_{sigma_bottom}_VM.png")  # Von Mises
    time.sleep(interval)
    take_screenshot(path + f"/{mesh}_{sigma_bottom}_ZZ.png")  # Sigma ZZ
    time.sleep(interval)
    take_screenshot(path + f"/{mesh}_{sigma_bottom}_CM.png")  # Coulomb Mohr YY

def run(sigma_bottom):
    # Run les 3 projets
    # Le current directory
    os.chdir("c:/Users/Boss/Desktop/LEPL1110/Projet/")
    original_dir = os.getcwd()

    # faut changer les valeurs dans le fichier problem.txt
    lines = []
    with open("data/problem.txt", "r") as file:
        lines = file.readlines()
        file.close()
    case = 0.1
    sigma_fluid = case * sigma_bottom
    # je change les valeurs
    lines[9] = f"Boundary condition :  Neumann-N          = -{sigma_fluid:.7e},            nan : Upper_curvature\n"
    lines[10] = f"Boundary condition :  Neumann-N          = -{sigma_fluid:.7e},            nan : Concave_curvature\n"
    lines[11] = f"Boundary condition :  Neumann-X          = -{sigma_fluid:.7e},            nan : Purple_line\n"
    lines[13] = f"Boundary condition :  Neumann-N          = -{sigma_bottom:.7e},            nan : Bottom_curve\n"
    lines[14] = f"Boundary condition :  Neumann-Y          =  {sigma_bottom:.7e},            nan : Bottom\n"
    with open("data/problem.txt", "w") as file:
        file.writelines(lines)
        file.close()
    # PreProcessor  => Il n'y a besoin de rien faire

    # Project
    os.chdir(original_dir)
    project_dir = "Project"
    os.chdir(project_dir)
    # os.system("make")
    os.chdir("build")
    command = '.\myFem.exe'
    os.system(command)

    # PostProcessor
    os.chdir(original_dir)
    # il va falloir prendre des screens avec la fonction capture_images
    thread = threading.Thread(target=capture_images, args=(str(case), sigma_bottom, 1.75, 1.5))
    command = '.\myFem.exe'
    thread.start()
    post_dir = "ProjectPostProcessor"
    os.chdir(post_dir)
    os.chdir("build")
    os.system(command)
    thread.join()



def main():
    nbr_image = 50
    sigma_bottoms = np.linspace(10**8, 8 * 10**9, nbr_image)
    for sigma_bottom in sigma_bottoms:
        run(sigma_bottom)

main()
