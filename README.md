# Cell-Cycle
Matlab program to stage cells in respective cell cycle phases using fluorescence microscope images

This is a brief guide to using the Matlab routine developed by Shivnarayan Dhuppar, TIFR Hyderabad, India.
if you use or modify it to your purpose, please cite the source: Dhuppar and Mazumder, Cell Cycle, 2018.

Copyright (c) 2018, Shivnarayan Dhuppar


Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the conditions in License file
are met.

GUIDE TO THE MATLAB ROUTINE:

First of all you need to add the *functions* folder to your search path in Matlab.
You can do this by selecting Home tab in your Matlab window and selecting Set Path
option from the tool bar. After this you can use the functions defined in the folder
from anywhere in your computer.

For background correction, you need to take a blankimage before every experiment or you can
add another way to background correct your images. I have provided a blank-image with each example
folder.

The function which segments the nuclei is NMask (you can find it in *functions* folder).
The arguments it takes are in the following order:

Image,DiskRadius,minArea,maxArea,minCircularity,maxCircularity,minRoundness,maxRoundness

You can change disk radius, minimum and maximum areas and minimum and maximum circularity and
roundness depending on the cell type. Generally, circularity and roundness are set equivalent but can
in priniciple can be different.

Deciding on the disk radius to smoothen out the nuclear intensity obtained from DNA channel is very
important -- you kept it lower and nuclei will be over-segmented; you kept it higher, nearby nuclei
will fuse. You have to play around with a few values -- you can first estimate the size using ImageJ
or some such software.

The other important thing to take in account is the order in which the data are printed.
The Intensity.m script which is to quantify immunofluorescence staining in a cell prints 
the output in the following order:
Column 1: DNA Content (DAPI) (Channel 1)
Column 2: p-H2A.X Intensity (Channel 3)
Column 3: Number of p-H2A.X Foci
Column 4: p53 Nuclear Intensity (CHannel 2)
Column 5L Nuclear to Cytoplasmic ratio of mean levels of p53

For Auto_mRNACount.m which combines cell cycle staging with the counting of mRNA (based on 
manual segmentation; for more info read Dhuppar et al. Cell Cycle 2018) the order is as following:
Column 1: DNA Content (DAPI) (Channel 1)
Column 2: Nuclear intensity of p53 (Channel 2)
Column 3: Nuclear to Cytoplasmic ratio of mean levels of p53
Column 4: mRNA Count called ManG (Channel 3)
Column 5: Cellular footprint

EXAMPLE IMAGES

An example set of images is given for both the cases in the respective folders. After adding the
*functions* folder in the search path you just have to run the main script, which is Intensity.m or 
Auto_mRNACount.m as the case be.
