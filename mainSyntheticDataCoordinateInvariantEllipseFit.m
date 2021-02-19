%--- DESCRIPTION -------------------------------------------------------------------------------------------------------------------
%mainSyntheticDataCoordinateInvariantEllipseFit is a main code that 
%generates synthetic data by using ellipticDataGenerator and fits an 
%ellipse by using coordinateInvariantEllipseFit. The user can chance the
%input struct:
%-EllipticDataPointGeneratorInput = 
%   struct with fields:
%               xCoordinateCenter: a real number (example 2)
%               yCoordinateCenter: a real number (example -1)
%             lengthSemiMajorAxis: a positive real number (example 5.3).
%             lengthSemiMinorAxis: a positive real number (example 3). Note that minor axis must not be longer than major axis
%             obliqueAngleRadians: a real number in an interval [0, pi] (example pi/3)
%     dataPointsPolarAngleRadians: a list of real numbers in an interval [0, 2*pi] (example [pi:pi/80:3*pi/2])
%         syntheticNoiseAmplitude: a positive real number (example 0.0500).
%
%--- DEPENDENCIES ---------------------------------------------------------
%This main code is using ellipticDataGenerator and
%coordinateInvariantEllipseFit functions
%
%--- REFERENCES -----------------------------------------------------------
%Osmund Francis (https://math.stackexchange.com/users/149569/osmund-francis), 
%Determining the major/minor axes of an ellipse from general form, URL 
%(version: 2017-02-11): https://math.stackexchange.com/q/820896
%
%--- DEVELOPPER & PROJECT -------------------------------------------------
%This code developped by Hakki Karaman (hakki.karaman@cern.ch) in February 
%2021 for CERN research project 'Novel investigations on vertical two-phase 
%CO2 flow to automatically identify the flow patterns and produce flow 
%regime maps by using pattern recognition algorithms on high speed camera 
%images for the new generation CO2 cooling systems of the ATLAS Experiment'
%--------------------------------------------------------------------------

%The user can change following parameters
EllipticDataPointGeneratorInput.xCoordinateCenter = 1;
EllipticDataPointGeneratorInput.yCoordinateCenter = -7;
EllipticDataPointGeneratorInput.lengthSemiMajorAxis = 6;
EllipticDataPointGeneratorInput.lengthSemiMinorAxis = 3;
EllipticDataPointGeneratorInput.obliqueAngleRadians = -pi/3;
EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians = [pi:pi/80:3*pi/2];
EllipticDataPointGeneratorInput.syntheticNoiseAmplitude = .05;
%--------------------------------------------------------------------------

%Creating synthetic elliptic data points
[dataPoints] = ellipticDataGenerator(EllipticDataPointGeneratorInput);

%It is fitting and plotting the optimal ellipse to the data points. 
%The ellipse is based on geometric fitting i.e. minimizing the sum of the 
%squares of the orthogonal distances from the data points to the curve. 
%It is rotation invariant.
[Ellipse] = coordinateInvariantEllipseFit(dataPoints)