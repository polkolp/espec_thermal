function [NoOutlier] = Outlier(phasecube,info)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalibrateStage v .1 Sep/18/2016            
% Bill Green, James Corsetti          
% Description: Look at the way that the distribution evolves with time. Notice how there
% are lobes that develop at regular intervals, like an interference
% pattern. This script looks for the data points that are more than one
% standard deviation from the mean measured phase for the region, and then
% replaces that with an interpolated value. 

if nargin < 2, info = []; end
% supply default parameters
%if isempty(numsteps), numsteps = 20; end

if info ==1
figure;
clear F
k = 1;
for i = 1:10:size(phasecube,3)
histogram(phasecube(:,:,k),100);
%plot(OPL(50*i,:))
ylim([0,size(phasecube,1)*size(phasecube,2)])
xlim([-6*pi,6*pi])
xlabel('Measured Drift \mum')
ylabel('Pixel Count')
F(k) = getframe
k = k+1
end
end

% If, in each frame, the measurement for a pixel is further away from the
% mean than one standard deviation, flag it.

% This creates a grid on which the pixels sit
[X,Y] = meshgrid(1:size(phasecube,2),1:size(phasecube,1));

% This is the loop that corrects outliers
for i = 1:size(phasecube,3);

phaseframe = phasecube(:,:,i);

%threshold = std2(phaseframe);
threshold = pi;
    
%For each frame, create a logical matrix that flags values more than the threshold from the mean    
isOutlier = abs(mean2(phaseframe) - phaseframe) > abs(mean2(phaseframe) - threshold);

%quality(i) = length(find(isOutlier>0)); %For the 5-1-SNR test, this shows that about 10% of the data points are outliers
% Creates an interpolating function using only the data points that are not
% outliers
F = scatteredInterpolant(X(~isOutlier & ~mod(X,5)),Y(~isOutlier & ~mod(X,5)),phaseframe(~isOutlier & ~mod(X,5)));
% Queries the interpolant to populate the entire X,Y grid of pixels.
% Outlier values are replaced with the interpolated value.

% This is a temporary matrix that handles the index assignment between
% logicals 
phaseframe(isOutlier) = F(X(isOutlier),Y(isOutlier));
NoOutlier(:,:,i) = phaseframe;
%before(i) = mean2(phasetemp);
%after(i) = mean2(

% Alternate method not working
% Need to update 2-D interpolation to use the ~argument in the X and Y
% data. Otherwise, the outlier values are used to interpolate as well
%all_idx(~outlier_idx), x(~outlier_idx), all_idx(outlier_idx))
%phasetemp(isOutlier) = interp2(X(~isOutlier),Y(~isOutlier),phasetemp(~isOutlier),X(isOutlier),Y(isOutlier));
i
end


