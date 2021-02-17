function [Ellipse] = coordinateInvariantEllipseFit(dataPoints)
%--- DESCRIPTION -------------------------------------------------------------------------------------------------------------------
%coordinateInvariantEllipseFit is a function, fits and plots the optimal 
%ellipse to the real world or synthetic elliptic data points. 
%The ellipse is based on geometric fitting i.e. minimizes the sum of the 
%squares of the orthogonal distances from the data points to the curve. It
%is rotation invariant.
%--- INPUT ----------------------------------------------------------------
%-dataPoints: numberOfDataPoints×2 matrix, the first coloumn is x 
%coordinate, the second coloumn is y coordinate of the data points.
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
%             coefficientsOfQuadric
%                taubinsErrorMetric
%                     geometricRMSE
%        taubinsErrorMetricVersion2
%        taubinsErrorMetricVersion3
%        orderOfEigenvalue
%--- DEPENDENCIES ---------------------------------------------------------
%This function is used by mainSyntheticDataCoordinateInvariantEllipseFit
%and it is using ellipseParametersCalculation.
%--- REFERENCES -----------------------------------------------------------
%-Taubin 1991, 'Estimation of planar curves, surfaces, and nonplanar space curves 
%defined by implicit equations with applications to edge and range image segmentation'
%https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=103273
%-Hakki's Derivations for Trivariate Implicit Full Quadric Surface Fitting
%(Taubin's Method) - uploaded on GitHub
%--- DEVELOPPER & PROJECT -------------------------------------------------
%This code developped by Hakki Karaman (hakki.karaman@cern.ch) in February 
%2021 for CERN research project 'Novel investigations on vertical two-phase 
%CO2 flow to automatically identify the flow patterns and produce flow 
%regime maps by using pattern recognition algorithms on high speed camera 
%images for the new generation CO2 cooling systems of the ATLAS Experiment'
%--------------------------------------------------------------------------
n = size(dataPoints,1); %number of points in the dataset
xCoordinateDataPoints = (dataPoints(:,1))'; %the first coloumn is x coordinates of the data
yCoordinateDataPoints = (dataPoints(:,2))'; %the second coloumn is y coordinates of the data
% General ellipse equation c1 + c2*x + c3*y + c4*x^2 + c5*x*y + c6*y^2 == 0
% OR [c1 c2 c3 c4 c5 c6]*[1 x y x^2 x*y y^2]' == 0 or c*X ==0
%matrixM = [1; x; y; x^2; x*y; y^2]*[1; x; y; x^2; x*y; y^2]' 
%size of matrixM is (6,n)×(n,6)=(6,6)
matrixM = 1/n*([ones(1,n); xCoordinateDataPoints; yCoordinateDataPoints; xCoordinateDataPoints.^2; xCoordinateDataPoints.*yCoordinateDataPoints; yCoordinateDataPoints.^2]*...
          [ones(1,n); xCoordinateDataPoints; yCoordinateDataPoints; xCoordinateDataPoints.^2; xCoordinateDataPoints.*yCoordinateDataPoints; yCoordinateDataPoints.^2]');
%matrixN = Jacobian([1; x; y; x^2; x*y; y^2])*(Jacobian([1; x; y; x^2; x*y; y^2]))'
%size of matrixN is (6,2*n)×(2*n,6)=(6,6)
matrixN = 1/n*([zeros(1,2*n); ones(1,n) zeros(1,n); zeros(1,n) ones(1,n); 2*xCoordinateDataPoints zeros(1,n); yCoordinateDataPoints xCoordinateDataPoints; zeros(1,n)  2*yCoordinateDataPoints]*...
          [zeros(1,2*n); ones(1,n) zeros(1,n); zeros(1,n) ones(1,n); 2*xCoordinateDataPoints zeros(1,n); yCoordinateDataPoints xCoordinateDataPoints; zeros(1,n)  2*yCoordinateDataPoints]');
%Solving Generalized Eigenvalue Problem M*c == lambda*N*c
[matrixEigenVectors,matrixEigenValues] = eig(matrixM, matrixN);
%matrixEigenValues is a (6,6) diagonal matrix. There are 6 eigenvalues in descending order,
%the last one is the smallest and represents the minimum lambda i.e. the
%"Taubin's error metric" or "Geometric Mean Square Error" (Geometric MSE) 
%matrixEigenVectors is a (6,6) matrix. Each coloumn corresponds the
%eigenvector [c1 c2 c3 c4 c5 c6] of kth eigenvalue.  
k = 6; %the get the smallest eigenvalue and corrsponding eigenvector
%Checks whether the best fitting quadric to the data is an ellipse
while matrixEigenVectors(5,k)^2-4*matrixEigenVectors(4,k)*matrixEigenVectors(6,k) >= 0 
    k = k - 1; %if it was not an ellipse checks the eigenvalues in order till it finds the ellipse
    %if the smallest eigenvalue was not corresponding to an ellipse most 
    %probably the Geometric RMSE will not be small enough
end
%M*c == lambda*N*c if and only if lambda == (c'*M*c)/(c'*N*c). If
%c1=c/sqrt(c'*N*c) then lambda == (c1'*M*c1). Lambda and 2 equations of
%lambda they are all Taubin's error metric and they need to be the same number.
%Normalization of the coefficientsOfQuadric to achive lambda == (c1'*M*c1)
coefficientsOfQuadric = matrixEigenVectors(:,k)/sqrt((matrixEigenVectors(:,k))'*matrixN*matrixEigenVectors(:,k));
%From 6 coefficients of full quadric equations calculating 5 independent ellipse parameters
[Ellipse] = ellipseParametersCalculation(coefficientsOfQuadric); 
Ellipse.coefficientsOfQuadric = coefficientsOfQuadric;
Ellipse.taubinsErrorMetric = coefficientsOfQuadric'*matrixM*coefficientsOfQuadric;
Ellipse.geometricRMSE = sqrt(Ellipse.taubinsErrorMetric);
Ellipse.taubinsErrorMetricVersion2 = matrixEigenValues(k,k);
Ellipse.taubinsErrorMetricVersion3 = ((matrixEigenVectors(:,k))'*matrixM*matrixEigenVectors(:,k))/((matrixEigenVectors(:,k))'*matrixN*matrixEigenVectors(:,k));
Ellipse.orderOfEigenvalue = k; 
%If k is not 6, it means there was another quadric fitting better to the data
%3 versions of Taubins Error Must be the same, after testing and verifying 
%the code with big amount of real data two of them can be removed 
%plotting
figure(1)
subplot(1,2,1); %plotting the data points
plot(xCoordinateDataPoints,yCoordinateDataPoints,'.b'); 
title('Data Points')

subplot(1,2,2); 
plot(xCoordinateDataPoints,yCoordinateDataPoints,'.b',Ellipse.xCoordinateCenter,Ellipse.yCoordinateCenter,'+r'); %plotting the data points
hold on
text(Ellipse.xCoordinateCenter,Ellipse.yCoordinateCenter, sprintf('(%.2f,%.2f)', Ellipse.xCoordinateCenter, Ellipse.yCoordinateCenter)) %ploting the center of the ellipse
plot(Ellipse.xCoordinateMajorAxisEndPoints,Ellipse.yCoordinateMajorAxisEndPoints,'-or'); %ploting the major axis of the ellipse
plot(Ellipse.xCoordinateMinorAxisEndPoints,Ellipse.yCoordinateMinorAxisEndPoints,'-or'); %ploting the minor axis of the ellipse
fimplicit(@(X,Y) coefficientsOfQuadric(1) + coefficientsOfQuadric(2)*X + coefficientsOfQuadric(3)*Y + ...
    coefficientsOfQuadric(4)*X.^2 + coefficientsOfQuadric(5)*X.*Y + coefficientsOfQuadric(6)*Y.^2, 'r'); %plotting the fitting ellipse
text((Ellipse.xCoordinateCenter + Ellipse.xCoordinateMajorAxisEndPoints(2))/2,(Ellipse.yCoordinateCenter + Ellipse.yCoordinateMajorAxisEndPoints(2))/2,sprintf('a = %.2f', Ellipse.lengthSemiMajorAxis'))
text((Ellipse.xCoordinateCenter + Ellipse.xCoordinateMinorAxisEndPoints(2))/2,(Ellipse.yCoordinateCenter + Ellipse.yCoordinateMinorAxisEndPoints(2))/2,sprintf('b = %.2f', Ellipse.lengthSemiMinorAxis'))
text(Ellipse.xCoordinateMajorAxisEndPoints(1), Ellipse.yCoordinateMajorAxisEndPoints(1),sprintf(' \\angle \\theta = %.1f\\circ',(Ellipse.obliqueAngleRadians*180/pi)));

 xlim_reference = xlim;
 ylim_reference = ylim;
text(xlim_reference(1), (0.9*ylim_reference(2)+0.1*ylim_reference(1)),sprintf('geometric RMSE = %.3f',Ellipse.geometricRMSE));
title('Coordinate System Invariant Ellipse Fitting to the Data')
grid on
axis equal
xlabel('x (mm)') 
ylabel('y (mm)')
%using the same x-y limits for both graphs
xlim([xlim_reference(1), xlim_reference(2)]);
ylim([ylim_reference(1), ylim_reference(2)]);
hold off 
subplot(1,2,1);
grid on
axis equal
xlim([xlim_reference(1), xlim_reference(2)]);
ylim([ylim_reference(1), ylim_reference(2)]);
xlabel('x (mm)') 
ylabel('y (mm)') 
end