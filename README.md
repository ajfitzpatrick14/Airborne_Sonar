# Photoacoustic Airborne Sonar System (PASS): 
**_Multi-modal sensor fusion for 3D airborne sonar imaging in hydrodynamic conditions_**


This repo contains code for image reconstruction, simulations, and system analysis for photoacoustic airborne sonar imaging. More information on airborne sonar imaging can be found here.

All example experimental and simulation data can be found here. Place the data in a folder on the path called 'Data'.

**In the 'Image Reconstruction' folder:**

  - ImageReconMain.m: This script is the main image reconstruction script which can be used to reconstruct images from the simulated and experimental datasets.
  - Functions: This folder contains custom MATLAB functions which perform data pre-processing and image reconstruction (GPW-SAR).
   

**In the 'Airborne Sonar Simulator' folder:**

  - kWave_Airborne_Sonar_Simulator.m: This script uses the k-Wave MATLAB toolbox to simulate photoacoustic airborne sonar imaging in user-defined conditions. The k-Wave package must be downloaded from here and added to the path.
  - Surface_Generator.m: This function is used within the simulator to define surface waves based on the user-defined wind conditions. 

**In the 'System Scripts' folder:**

  - PASS_Equation.m: This script can be used to analyze the maximum imaging depth of PASS as a function of user-defined inputs. 
  - Surface_Processing_Pipeline.m: This script converts raw SR305 point clouds to post-processed surface maps defined over the a region of interest.


Scripts and functions are well-annotated to provide further details on the user-defined inputs/outputs.





