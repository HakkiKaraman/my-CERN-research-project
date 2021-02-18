# my-CERN-research-project
This repository contains some of my code, mathematical derivations, result figures as well as a reference paper for my research project at CERN "Novel investigations on vertical two-phase CO2 flow to automatically identify the flow patterns and produce flow regime maps by using pattern recognition algorithms on high speed camera images for the new generation CO2 cooling systems of the ATLAS Experiment". My first implementation to improve the results of existing gas void fraction estimation pipeline is 'Rotation Invariant, Geometric Ellipse Fitting'.

On preprocessed 2D bubbly flow images, detecting bubbles by using ellipse detection algorithms are crucial part of gas void fraction estimation pipeline. In 2D, the classical least square algorithm finds a curve that minimizes the sum of the squares of the vertical distances from the data points to the curve. In contrast, the geometric fitting finds the curve that minimizes the sum of the squares of the orthogonal distances from the data points to the curve (which is equivalent to minimize the Geometric Root Mean Square Error - Geometric RMSE) in stead of vertical distances. So, geometric fitting is coordinate system invariant. 

The references are the original paper from Taubin 1991 - 'Estimation of Planar Curves, Surfaces, and Nonplanar Space Curves Defined by Implicit Equations with Applications to Edge and Range Image Segmentation' that described geometric fitting method in more general terms. Furthermore, there are my derivations for 3D case with the name 'Hakki Karaman's Derivations for Trivariate Implicit Full Quadric Surface Fitting (Taubin's Method)'. Here, my implementation is a 2D reduced version of 3D equations.

coordinateInvariantEllipseFit.m is a function that is ready to use for Rotation Invariant Ellipse Fitting with real data. This function uses ellipseParametersCalculation.m function to calculate the ellipse parameters from the general quadric equation of ellipse.

mainSyntheticDataCoordinateInvariantEllipseFit.m can be used use for Rotation Invariant Ellipse Fitting with synthetic data. First, it creates the data by using ellipticDataGenerator.m then fits the ellipse by using coordinateInvariantEllipseFit.m.

Also, there are 3 result images.
