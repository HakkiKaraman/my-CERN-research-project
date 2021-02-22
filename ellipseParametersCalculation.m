function [Ellipse] = ellipseParametersCalculation(coefficientsOfQuadric)
%--- DESCRIPTION ----------------------------------------------------------
%ellipseParametersCalculation is a function that takes 6 coefficients of 
%full quadric equations and calculates 5 independent ellipse parameters 
%(lengthSemiMajorAxis, lengthSemiMinorAxis, xCoordinateCenter, yCoordinateCenter, obliqueAngleRadians)
%and calculates additional dependent parameters by using 5 parameters
%
%--- INPUT ----------------------------------------------------------------
%-coefficientsOfQuadric == [c1 c2 c3 c4 c5 c6] 
%       where c1 + c2*x + c3*y + c4*x^2 + c5*x*y + c6*y^2 == 0
%
%--- OUTPUT ---------------------------------------------------------------
%-Ellipse = 
%   struct with fields:
%               lengthSemiMajorAxis 
%               lengthSemiMinorAxis 
%                 xCoordinateCenter 
%                 yCoordinateCenter 
%               obliqueAngleRadians
%     xCoordinateMajorAxisEndPoints 
%     yCoordinateMajorAxisEndPoints 
%     xCoordinateMinorAxisEndPoints
%     yCoordinateMinorAxisEndPoints 
%                       ellipseArea 
%
%--- DEPENDENCIES ---------------------------------------------------------
%This function is used by coordinateInvariantEllipseFit function
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

[f, d, e, a, b, c] = feval(@(x) x{:}, num2cell(coefficientsOfQuadric)); 
%assigning the coefficients [f, d, e, a, b, c] = [c1 c2 c3 c4 c5 c6] 
%initial form c1 + c2*x + c3*y + c4*x^2 + c5*x*y + c6*y^2 = 0
%final form a*x^2 + b*x*y + c*y^2 + d*x + e*y + f = 0

delta = b^2 - 4*a*c; %delta < 0 is the indication of ellipse
q = 64*(-f*delta - a*e^2 + b*d*e - c*d^2)/delta^2; %q is the coefficient normalising factor
s = (1/4)*(q^2*(b^2 + (a - c)^2))^(1/4); %distance between center and focal point

Ellipse.lengthSemiMajorAxis = 1/8*sqrt(2*abs(q)*sqrt(b^2 + (a - c)^2) - 2*q*(a + c)); %Ellipse semi major axis lenght
Ellipse.lengthSemiMinorAxis = sqrt(Ellipse.lengthSemiMajorAxis^2 - s^2); %Ellipse semi minor axis lenght
Ellipse.xCoordinateCenter = (2*c*d - b*e)/delta; %Center of the ellipse - x coordinate
Ellipse.yCoordinateCenter = (2*a*e - b*d)/delta; %Center of the ellipse - y coordinate

%calculating the oblique angle of the ellipse 'theta' for different cases
%For the mathematical derivations see the references
if (q*a - q*c) == 0
    if q*b == 0
        theta = 0;
    elseif q*b > 0
        theta = pi/4;
    else
        theta = 3/4*pi;
    end
elseif (q*a - q*c) > 0
    if q*b >= 0
        theta = 1/2*atan(b/(a-c));
    else
        theta = 1/2*atan(b/(a-c)) + pi;
    end
else
    theta = 1/2*atan(b/(a-c)) + pi/2;
end  

Ellipse.obliqueAngleRadians = theta;

%For the mathematical derivations see the references
Ellipse.xCoordinateMajorAxisEndPoints = [(Ellipse.xCoordinateCenter - Ellipse.lengthSemiMajorAxis*cos(theta))...
                                         (Ellipse.xCoordinateCenter + Ellipse.lengthSemiMajorAxis*cos(theta))];
                                     
Ellipse.yCoordinateMajorAxisEndPoints = [(Ellipse.yCoordinateCenter - Ellipse.lengthSemiMajorAxis*sin(theta))...
                                         (Ellipse.yCoordinateCenter + Ellipse.lengthSemiMajorAxis*sin(theta))];
                                     
Ellipse.xCoordinateMinorAxisEndPoints = [(Ellipse.xCoordinateCenter + Ellipse.lengthSemiMinorAxis*sin(theta))...
                                         (Ellipse.xCoordinateCenter - Ellipse.lengthSemiMinorAxis*sin(theta))];
                                     
Ellipse.yCoordinateMinorAxisEndPoints = [(Ellipse.yCoordinateCenter - Ellipse.lengthSemiMinorAxis*cos(theta))...
                                         (Ellipse.yCoordinateCenter + Ellipse.lengthSemiMinorAxis*cos(theta))];
                                     
Ellipse.ellipseArea = pi*Ellipse.lengthSemiMajorAxis*Ellipse.lengthSemiMinorAxis;
end
