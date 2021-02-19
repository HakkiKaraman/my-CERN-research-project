function [dataPoints] = ellipticDataGenerator(EllipticDataPointGeneratorInput)
%--- DESCRIPTION ----------------------------------------------------------
%ellipticDataGenerator is a function that takes struct 
%EllipticDataPointGeneratorInput and creates synthetic dataPoints matrix.
%After creation of the points artificial noise is added for each point.
%Adding Gaussian noise to x-y coordinate is not realistic, so x-y noise are
%calculated from the Gaussian noise amplitude and random direction.
%
%--- INPUT ----------------------------------------------------------------
%EllipticDataPointGeneratorInput = 
%  struct with fields:
%               xCoordinateCenter: a real number (example 2)
%               yCoordinateCenter: a real number (example -1)
%             lengthSemiMajorAxis: a positive real number (example 5.3)
%             lengthSemiMinorAxis: a positive real number (example 3). Note that minor axis must not be longer than major axis
%             obliqueAngleRadians: a real number in an interval [0, pi] (example pi/3)
%     dataPointsPolarAngleRadians: a list of real numbers in an interval [0, 2*pi] (example [pi:pi/80:3*pi/2])
%         syntheticNoiseAmplitude: a positive real number (example 0.0500).
%
%--- OUTPUT ---------------------------------------------------------------
%-dataPoints: numberOfDataPoints×2 matrix, the first coloumn is x coordinate,
%the second coloumn is y coordinate of the data points.
%
%--- DEPENDENCIES ---------------------------------------------------------
%This function is used by mainSyntheticDataCoordinateInvariantEllipseFit
%
%--- REFERENCES -----------------------------------------------------------
%Gil Epshtain (https://math.stackexchange.com/users/529131/gil-epshtain), 
%What is the parametric equation of a rotated Ellipse (given the angle of rotation), 
%URL (version: 2018-05-06): https://math.stackexchange.com/q/2647450
%
%--- DEVELOPPER & PROJECT -------------------------------------------------
%This code developped by Hakki Karaman (hakki.karaman@cern.ch) in February 
%2021 for CERN research project 'Novel investigations on vertical two-phase 
%CO2 flow to automatically identify the flow patterns and produce flow 
%regime maps by using pattern recognition algorithms on high speed camera 
%images for the new generation CO2 cooling systems of the ATLAS Experiment'
%--------------------------------------------------------------------------

%Creating 2D noise in polar coordinates to have Gaussian distribution in
%radial direction and randomness in azhimutal

%Gaussian noise amplitude, the length is numberOfDataPoints 
noise_amplitude = EllipticDataPointGeneratorInput.syntheticNoiseAmplitude*randn(1,length(EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians)); 

%random angle in an interval[0,2*pi], the length is numberOfDataPoints
noise_angle = 2*pi*rand(1,length(EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians)); 

%for mathematical derivation see the reference 
xCoordinateDataPoints = EllipticDataPointGeneratorInput.lengthSemiMajorAxis*cos(EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians)*cos(EllipticDataPointGeneratorInput.obliqueAngleRadians) - ...
                        EllipticDataPointGeneratorInput.lengthSemiMinorAxis*sin(EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians)*sin(EllipticDataPointGeneratorInput.obliqueAngleRadians) + ...
                        EllipticDataPointGeneratorInput.xCoordinateCenter + noise_amplitude.*cos(noise_angle);

%for mathematical derivation see the reference                   
yCoordinateDataPoints = EllipticDataPointGeneratorInput.lengthSemiMajorAxis*cos(EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians)*sin(EllipticDataPointGeneratorInput.obliqueAngleRadians) + ...
                        EllipticDataPointGeneratorInput.lengthSemiMinorAxis*sin(EllipticDataPointGeneratorInput.dataPointsPolarAngleRadians)*cos(EllipticDataPointGeneratorInput.obliqueAngleRadians) + ...
                        EllipticDataPointGeneratorInput.yCoordinateCenter + noise_amplitude.*sin(noise_angle); 
                                        
dataPoints = zeros(length(xCoordinateDataPoints),2); %initialization, numberOfDataPoints×2 matrix
dataPoints(:,1) = xCoordinateDataPoints; %the first coloumn corresponds x coordinates the data points
dataPoints(:,2) = yCoordinateDataPoints; %the second coloumn corresponds y coordinates the data points

%plotting the data points
figure(2)
plot(xCoordinateDataPoints,yCoordinateDataPoints,'.b');
grid on
axis equal
end