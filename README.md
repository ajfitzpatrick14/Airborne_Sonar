# Photoacoustic Airborne Sonar System (PASS): 
**_Multi-modal sensor fusion for 3D airborne sonar imaging in hydrodynamic conditions_**


This repo contains code for image reconstruction, simulations, and system analysis for our photoacoustic airborne sonar imaging system. More information on airborne sonar imaging can be found [here](https://ieeexplore.ieee.org/document/9228880).

All example experimental and simulation data can be found [here](https://drive.google.com/file/d/1dF0gsPam3Wmj6cTYlOa-9wAmoAFxPWH3/view?usp=sharing). Place the data in a folder on the path called 'Data'.

**In the 'Image_Reconstruction' folder:**

  - ImageReconMain.m: This script is the main image reconstruction script which can be used to reconstruct images from the simulated and experimental datasets.
  - Functions: This folder contains custom MATLAB functions which perform the image reconstruction (GPW-SAR). Note: some functions have been written explicitly for specific data files and therefore may not always generalize without prior edits.
   

**In the 'Airborne_Sonar_Simulator' folder:**

  - kWave_PASS_Simulator.m: This script uses the k-Wave MATLAB toolbox to simulate photoacoustic airborne sonar imaging in user-defined conditions. The k-Wave package must be downloaded from [here](http://www.k-wave.org/index.php) and added to the path.
  - Surface_Generator.m: This function is used within the simulator to define surface waves based on the user-defined wind conditions. 

**In the 'System_Scripts' folder:**

  - PASS_Equation.m: This script can be used to analyze the maximum imaging depth of PASS as a function of user-defined inputs. 
  - Surface_Processing_Pipeline.m: This script converts raw SR305 point clouds to post-processed surface maps defined over the a region of interest.


Scripts and functions are well-annotated to provide further details on the user-defined inputs/outputs.





