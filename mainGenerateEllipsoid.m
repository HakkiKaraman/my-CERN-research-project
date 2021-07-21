%--- DESCRIPTION ----------------------------------------------------------
%mainGenerateEllipsoid is a function that creates Monte Carlo Simulations
%with different models. It creates ellipsoids and their projections in random
%directions and fits models to solve an inverse problem (from projections to 
%find out the mean ellipsoids) specifically for bubbly two phase flow. 
%
%--- DEPENDENCIES ---------------------------------------------------------
%uses: generatingRandomUnitProjectionVectors
%used by: ANNmyNeuralNetworkFunctionMatrixOnly(majorMinorAxisProjectedEllipse)
%--- REFERENCES -----------------------------------------------------------
%-W.M.Yan et al.2017, 'Inferring 3D particle size and shape characteristics 
% from projected 2D images: Lessons learned from ellipsoids'
% https://www.sciencedirect.com/science/article/abs/pii/S0266352X17303178
%
%--- DEVELOPPER & PROJECT -------------------------------------------------
%This code developped by Hakki Karaman (hakki.karaman@cern.ch) in March 
%2021 for CERN research project 'Novel investigations on vertical two-phase 
%CO2 flow to automatically identify the flow patterns and produce flow 
%regime maps by using pattern recognition algorithms on high speed camera 
%images for the new generation CO2 cooling systems of the ATLAS Experiment'
%--------------------------------------------------------------------------

close all
clear all

%can be modified
numberOfEllipsoids = 10^5;
numberOfBins = 20;
numberOfRandomUnitProjectionVectors = 1;
semiAxisEllipsoidProjectionEllipse = zeros(numberOfEllipsoids,11);
maximumBubbleRadius = 4; %the inner diameter of the tube is 8 mm
minimumDetectableBubbleRadius = 0.05; %in mm
aspectRatio = 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:numberOfEllipsoids
    characteristicR = minimumDetectableBubbleRadius + (maximumBubbleRadius-minimumDetectableBubbleRadius)*rand(1);
    supR = min(characteristicR*sqrt(aspectRatio),maximumBubbleRadius);
    infR = max(characteristicR/sqrt(aspectRatio),minimumDetectableBubbleRadius);
    R = sort(infR + (supR-infR)*rand([1 3]),'descend');    
    A = diag(R.^-2);
    %Check the reference
    [randomUnitProjectionVectors] = generatingRandomUnitProjectionVectors(numberOfRandomUnitProjectionVectors);
    v = randomUnitProjectionVectors';
    B = A - (A*v*v'*A)/(v'*A*v); 
    r = [1/median(sqrt(eig(B))) 1/max(sqrt(eig(B)))]; %note that eig(B) is a vector with 3 dimensions and min(eig(B))== 0 is meaningless
    volumeEllipsoid.groundTruth = 4/3*pi*R(1)*R(2)*R(3);
    areaProjectedEllipse = pi*r(1)*r(2); 
    
    volumeEllipsoid.arithmetricAverageApproximation = 4/3*pi*r(1)*r(2)*(r(1)+r(2))/2; %Arithmetric average based approximation
    volumeEllipsoid.geometricAverageApproximation = 4/3*pi*(r(1)*r(2))^1.5; %Geometric average based approximation
    volumeEllipsoid.harmonicAverageApproximation = 8/3*pi*(r(1)*r(2))^2/(r(1)+r(2)); %Harmonic average based approximation
    volumeEllipsoid.l2normApproximation = 4/3*pi*r(1)*r(2)*norm([r(1) r(2)],2)/sqrt(2); %L2 norm based approximation
    volumeEllipsoid.arithmetricAverageApproximationCoefficient = ...
        volumeEllipsoid.groundTruth/volumeEllipsoid.arithmetricAverageApproximation;
    volumeEllipsoid.geometricAverageApproximationCoefficient = ...
        volumeEllipsoid.groundTruth/volumeEllipsoid.geometricAverageApproximation;
    volumeEllipsoid.harmonicAverageApproximationCoefficient = ...
        volumeEllipsoid.groundTruth/volumeEllipsoid.harmonicAverageApproximation;
    volumeEllipsoid.l2normApproximationCoefficient = ...
        volumeEllipsoid.groundTruth/volumeEllipsoid.l2normApproximation;
    
    semiAxisEllipsoidProjectionEllipse(i,:) = [areaProjectedEllipse volumeEllipsoid.groundTruth ...
        volumeEllipsoid.arithmetricAverageApproximationCoefficient volumeEllipsoid.geometricAverageApproximationCoefficient ...
        volumeEllipsoid.harmonicAverageApproximationCoefficient volumeEllipsoid.l2normApproximationCoefficient R r];
end

majorMinorAxisProjectedEllipse = semiAxisEllipsoidProjectionEllipse(:,10:11);
volumeEllipsoid = semiAxisEllipsoidProjectionEllipse(:,2); 

xlimMax = 50;
figure(1)
plot(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,2),'.')
grid on

 f1 = fit(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,2),'poly3');
 hold on
 %plot(f1)
[N, edges, bins] = histcounts(semiAxisEllipsoidProjectionEllipse(:,1),numberOfBins); %put the numbers into 36 bins
midPoints = (edges(1:end-1))'+(edges(2)-edges(1))/2;
meanY = splitapply(@mean, semiAxisEllipsoidProjectionEllipse(:,2), bins);
stairs(edges,[meanY; meanY(end)]);
size(meanY)
size(edges)
title('Ellipsoid Volume vs. Projected Ellipsoid Area');
xlabel('Projected Ellipsoid Area (mm^2)') 
ylabel('Ellipsoid Volume (mm^3)')
xlim([0 xlimMax]);
hold off

figure(2)
plot(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,3),'.')
grid on
f2 = fit(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,3),'poly3');
hold on
%plot(f2)
[N, edges, bins] = histcounts(semiAxisEllipsoidProjectionEllipse(:,1),numberOfBins); %put the numbers into bins
midPoints = (edges(1:end-1))'+(edges(2)-edges(1))/2;
meanY = splitapply(@mean, semiAxisEllipsoidProjectionEllipse(:,3), bins);
%plot(midPoints,meanY,'.-r');
stairs(edges,[meanY; meanY(end)]);
title('Volume Approx. Arithmetic Mean Coefficient vs Projected Ellipsoid Area');
xlabel('Projected Ellipsoid Area (mm^2)') 
ylabel('Particular Arithmetic Mean Coefficient')
ylim([0 4])
xlim([0 xlimMax])
xlim_reference = xlim;
ylim_reference = ylim;
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.9*ylim_reference(1)+0.9*ylim_reference(2),...
    sprintf('Arithmetic Mean Coefficient = %.3f',mean(semiAxisEllipsoidProjectionEllipse(:,3))));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.95*ylim_reference(1)+0.85*ylim_reference(2),...
    sprintf('Bins Fluctuation Range = %.3f',(max(meanY)-min(meanY)))); 
hold off

l2normApproximationErrorPercentage = 100*(mean(semiAxisEllipsoidProjectionEllipse(:,6)) - 0.9142) %percentage

figure(3)
plot(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,4),'.')
grid on
f3 = fit(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,4),'poly3');
hold on
[N, edges, bins] = histcounts(semiAxisEllipsoidProjectionEllipse(:,1),numberOfBins); %put the numbers into 36 bins
midPoints = (edges(1:end-1))'+(edges(2)-edges(1))/2;
meanY = splitapply(@mean, semiAxisEllipsoidProjectionEllipse(:,4), bins);
%plot(midPoints,meanY,'.-r');
%plot(f3)
stairs(edges,[meanY; meanY(end)]);
title('Volume Approx. Geometric Mean Coefficient vs Projected Ellipsoid Area');
xlabel('Projected Ellipsoid Area (mm^2)') 
ylabel('Particular Geometric Mean Coefficient')
ylim([0 4])
xlim([0 xlimMax])
xlim_reference = xlim;
ylim_reference = ylim;
%text(xlim_reference(1), (0.9*ylim_reference(2)+0.1*ylim_reference(1)),...
%     sprintf('Geometric RMSE = %.3f',Ellipse.geometricRMSE));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.9*ylim_reference(1)+0.9*ylim_reference(2),...
    sprintf('Geometric Mean Coefficient = %.3f',mean(semiAxisEllipsoidProjectionEllipse(:,4))));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.95*ylim_reference(1)+0.85*ylim_reference(2),...
    sprintf('Bins Fluctuation Range = %.3f',(max(meanY)-min(meanY)))); 
hold off

figure(4)
plot(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,5),'.')
grid on
f4 = fit(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,5),'poly3');
hold on
%plot(f4)
 [N, edges, bins] = histcounts(semiAxisEllipsoidProjectionEllipse(:,1),numberOfBins); %put the numbers into 36 bins
 %midPoints = (edges(1:end-1))'+(edges(2)-edges(1))/2;
meanY = splitapply(@mean, semiAxisEllipsoidProjectionEllipse(:,5), bins);
%plot(midPoints,meanY,'.-r');
stairs(edges,[meanY; meanY(end)]);
title('Volume Approx. Harmonic Mean Coefficient vs Projected Ellipsoid Area');
xlabel('Projected Ellipsoid Area (mm^2)') 
ylabel('Particular Harmonic Mean Coefficient')
ylim([0 4])
xlim([0 xlimMax])
xlim_reference = xlim;
ylim_reference = ylim;
%text(xlim_reference(1), (0.9*ylim_reference(2)+0.1*ylim_reference(1)),...
%     sprintf('Geometric RMSE = %.3f',Ellipse.geometricRMSE));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.9*ylim_reference(1)+0.9*ylim_reference(2),...
    sprintf('Harmonic Mean Coefficient = %.3f',mean(semiAxisEllipsoidProjectionEllipse(:,5))));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.95*ylim_reference(1)+0.85*ylim_reference(2),...
    sprintf('Bins Fluctuation Range = %.3f',(max(meanY)-min(meanY)))); 
hold off
figure(5)
plot(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,6),'.')
grid on
hold on
%f5 = fit(semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,6),'poly3');
%data = [semiAxisEllipsoidProjectionEllipse(:,1),semiAxisEllipsoidProjectionEllipse(:,6)];
[N, edges, bins] = histcounts(semiAxisEllipsoidProjectionEllipse(:,1),numberOfBins); %put the numbers into 36 bins
%midPoints = (edges(1:end-1))'+(edges(2)-edges(1))/2;
meanY = splitapply(@mean, semiAxisEllipsoidProjectionEllipse(:,6), bins);
%plot(midPoints,meanY,'.-r');
stairs(edges,[meanY; meanY(end)]);
hold on
%plot(f5)
title('Volume Approx. L2-norm Mean Coefficient vs Projected Ellipsoid Area');
xlabel('Projected Ellipsoid Area (mm^2)') 
ylabel('Particular L2-norm Mean Coefficient')
ylim([0 4])
xlim([0 xlimMax])
xlim_reference = xlim;
ylim_reference = ylim;
%text(xlim_reference(1), (0.9*ylim_reference(2)+0.1*ylim_reference(1)),...
%     sprintf('Geometric RMSE = %.3f',Ellipse.geometricRMSE));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.9*ylim_reference(1)+0.9*ylim_reference(2),... 
    sprintf('L2-norm Mean Coefficient = %.3f',mean(semiAxisEllipsoidProjectionEllipse(:,6))));
text(0.5*xlim_reference(1)+0.5*xlim_reference(2), 0.95*ylim_reference(1)+0.85*ylim_reference(2),...
    sprintf('Bins Fluctuation Range = %.3f',(max(meanY)-min(meanY)))); 
hold off

figure(6) 
subplot(1,2,1)
ellipsoid(0,0,0,R(1),R(2),R(3))
scaling = 1.1*norm([R(1) R(2) R(3)]);
hold on
quiver3(0,0,0,scaling*v(1),scaling*v(2),scaling*v(3));
axis equal
title('Ellipsoid with projection direction');

subplot(1,2,2)
ellipsoid(0,0,0,r(1),r(2),0)
axis equal
title('Projection Ellipse');
view(0,90)
hold off