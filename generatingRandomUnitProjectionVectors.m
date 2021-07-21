function [randomUnitProjectionVectors] = generatingRandomUnitProjectionVectors(numberOfRandomUnitProjectionVectors)
%--- DESCRIPTION ----------------------------------------------------------
%generatingRandomUnitProjectionVectors is a function that generates random
%unit projection vectors in spherical coordinates. The projection vectors 
%can be expressed as points on a unit sphere. The probability density
%function is increaing with the increasing polar angle to maintain uniform 
%probability per unit spherical portion.
%
%--- INPUT ----------------------------------------------------------------
%-numberOfRandomUnitProjectionVectors: the number of random unit projection
%                                      vectors.
%
%--- OUTPUT ---------------------------------------------------------------
%-randomUnitProjectionVectors: numberOfRandomUnitProjectionVectors×3 matrix
%                              & each raw represent a 3D projection vector 
%
%--- DEPENDENCIES ---------------------------------------------------------
%
%--- REFERENCES -----------------------------------------------------------
%-W.M.Yan et al.2017, 'Inferring 3D particle size and shape characteristics 
% from projected 2D images: Lessons learned from ellipsoids'
% https://www.sciencedirect.com/science/article/abs/pii/S0266352X17303178
%
%--- DEVELOPPER & PROJECT -------------------------------------------------
%This code developped by Hakki Karaman (hakki.karaman@cern.ch) in 
%2021 for CERN research project 'Novel investigations on vertical two-phase 
%CO2 flow to automatically identify the flow patterns and produce flow 
%regime maps by using pattern recognition algorithms on high speed camera 
%images for the new generation CO2 cooling systems of the ATLAS Experiment'
%--------------------------------------------------------------------------

%Assigning azimuthal angle randomly
theta = 2*pi*rand([numberOfRandomUnitProjectionVectors 1]); 

%Assigning polar angle(adjusted to maintain the uniform probability per 
%unit spherical portion)
phi = acos(1-2*rand([numberOfRandomUnitProjectionVectors 1]));

%phi = pi*rand([numberOfRandomUnitProjectionVectors 1]); %to demonstrate if phi is uniformly choosen non-uniform distribution 


%From Spherical (on a Unit Sphere) to Cartesian Coordinate Transformation 
nx = cos(theta).*sin(phi);
ny = sin(theta).*sin(phi);
nz = cos(phi);

%Saving the Vectors in numberOfRandomUnitProjectionVectors×3 matrix
randomUnitProjectionVectors = [nx ny nz];

%Plotting to check the randomness. In normal cases no need to plot 
%plotUnitProjectionVectors can be set to 0
plotUnitProjectionVectors = 0;
if plotUnitProjectionVectors
    plot3(nx,ny,nz,'.')
    grid on
    axis equal
    xlabel('x') 
    ylabel('y')
    zlabel('z')
end
end
