This document explains how the Line Segment Detector code works. 

LSD is called from Python, but requires most files in this directory, and the Python code must point to this directory. The Matlab function called is call_lsd.m, and it is called from the Python function [vp_analysis.py] in the Python branch ../../project/src.

HOW TO MAKE SURE CODE IS WORKING
Copy an image into this directory, or use an image already existing in this directory. Open compile.m in a Matlab editor and set the filename which is read into the variable 'chairs' to the filename of the image you want to detect lines in. Then in the Matlab Command Window, navigate the Matlab path to this directory and call compile.m. The image will pop up in a figure, and if working, lines will overlay on top of the image after 10 seconds or so, depending on image size.