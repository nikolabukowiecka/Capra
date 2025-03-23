# Capra
HEALPix based system to create all sky maps of probability density distributions from event counting sky surveys. The algorithm was developed based on the IBEX mission characteristics, with the applications to IMAP mission, however could be applied to any event counting sky surveying missions. 

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

Maciej Bzowski from Polish Space Research Centre, Polish Academy of Sciences contributed with geometricPackage.wl, and with general design and science advice to the project.
Marzena Kubiak from Polish Space Research Centre, Polish Academy of Sciences contributed with the collimator function. 
This study was supported by the Polish Ministry for Education and Science under contract MEiN/2021/2/DIR. 
