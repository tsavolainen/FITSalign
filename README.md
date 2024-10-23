# FITSalign
Software to align VLBI images based on 2D cross-correlation.
Copyright (C) 2017 Tuomas Savolainen


usage: FITSalign.pl [-aco] <image1> <image2> [<x1> <y1> <x2> <y2>
                     <xtemp1> <ytemp1> <xtemp2> <ytemp2>] 

The purpose of this program is to align two FITS images 
by normalized 2D cross-correlation. The input images need to 
have the same pixel size and their PSFs (or restoring beams) 
should be identical. 

The program first plots the input images and asks the user
to select an appropriate area for the subsequent analysis.
Then the user is prompted to select a feature in the image to be 
matched. In the case of aligning radio images of AGN jets taken at 
different frequencies, it is important to select a feature 
that is optically thin and does not show strong spectral gradients. 
 
The program calculates a normalized cross-correlation between the 
first image and the selected feature of the second image using a 
fast method described by J.P. Lewis (1995). The correlation plot is 
shown along with useful information about the shift between the 
input images (calculated from the maximum of the correlation plot).  

The program also calculates two spectral index images: one 
with the calculated shift applied and one without. Check the 
results! The program gives a warning if the feature selected 
for matching turns out to contain optically thick parts. Finally, it 
is possible to write out the spectral index image (with the shift 
applied) as well as the correlation plot as FITS files.

Switches:

  -o=<outfile>   If <outfile> is given, the output is written to 
                 <outfile>.out and the spectral index image to 
                 <outfile>.fits.  

  -c=<icut>      Parameters <icut> and is the cut-off level used 
                 in calculating spectral index maps. They are 
                 optional: if not specified by the user, a value 
                 of 0.001 Jy/beam is used. 

  -a             Automatic mode: The user still needs to provide the 
                 coordinates of the image region and template in 
                 [<x1> <y1> <x2> <y2> <xtemp1> <ytemp1> <xtemp2> <ytemp2>]. 
                 No output graphics.

NOTE THAT BAD VALUE SUPPORT IN PDL IS NEEDED TO RUN THIS 
PROGRAM.

