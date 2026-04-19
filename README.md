# Capra

We present Capra, a HEALPix-native pipeline for reconstructing all-sky energetic neutral atom (ENA) intensity maps from the NASA Interstellar Boundary Explorer (IBEX) and  Interstellar Mapping and Acceleration (IMAP) missions' event-counting data. As one of the possibilities, Capra distributes events according to the normal distribution over the HEALPix grid. This produces smooth, physically interpretable maps with consistent resolution changes and reference-frame transformations. 

The reconstructed maps capture the large-scale morphology of the large-scale intensity while filtering out small-scale mottling and strip artifacts.  

We verify invariance (preservation) of counts, signals, rates, and intensities across tessellations under coordinate transformations, and demonstrate numerical convergence by N_side=64 for a representative dataset. 

We show that thread-based parallelism performance exhibit approximately linear weak scaling, with bounded memory usage, enabling efficient high-resolution reconstruction and making Capra well-suited for high-resolution IMAP reconstruction.

In the "Run the code" notebook:
1. Set up paths:
    a) baseDir - full path to your Capra folder
    b) healpixDir - in the example healpix folder there are a few tesselations generated. You can use healpix C library to generate files of the desired tesselation.
    c) ibexDataDir - in this folder you should have folders with your files, say "Map2018A_1deg_binned" and corresponding folders with the instruments spin axis say "goodTimes_Axis"
    d) outputDir - full path to your output folder, where all the maps data will be generated

2. Import the packages

3. Load the data and create the ibexHi data structure

4. Load the instrument properties. For IbexHi: 
scanRadius = 90 
colRadius = 10
energyStep = 1,2,3,4,5,6

5. Export the data. You will get two types of data:
a) One without the "Null" in the dataname - those serve for post-processing of the map (when loading to Mathematica from files), that will have pixels with 0 exposure time (unobserved) replaced with the word "bad" (for "bad pixel").
b) One with the "Null" in the dataname - those serve for straight up visualization with Python healpy library, that requires "Nulls", so the pixels with 0 exposure time will be replaced with Null.


In the "Postprocess the code" notebook:
You will find functionality such as:
1. Creating relative/flux maps

2. Changing the tessellationg

    a) uniformly 
    
    b) at a chosen point

3. Changing Pixel Shapes, ex: give them nxm resolution

4. Transforming the map to Ribbon Centered Coordinates

5. Transforming maps to one standard layout for comparison:

    a) HEALPix maps
    
    b) HEALPix -> Ribbon centered maps
    
    c) Raw 1 deg binned maps
    
    d) 6 deg binned maps (Re-binned from 1 deg binned maps)

Maciej Bzowski from Polish Space Research Centre, Polish Academy of Sciences contributed with geometricPackage.wl, and with general design and science advice to the project.
Marzena Kubiak from Polish Space Research Centre, Polish Academy of Sciences contributed with the collimator function. 
This study was supported by the Polish Ministry for Education and Science under contract MEiN/2021/2/DIR. 

