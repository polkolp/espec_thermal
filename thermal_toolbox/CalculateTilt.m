function [xtilt,ytilt] = CalculateTilt(phasetemp,coordsTilt)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalculateTilt v .1 Jan/6/2016            
% Bill Green, James Corsetti          
% Description: Uses an interferogram to figure out what the best fit tilt
% is in the x and y directions of the image. Outputs the xtilt and ytilt in
% units of microns

% tiltplane - formally "b" - the 2D spatially unwrapped phase of tiltregion
% columns - the x pixel values of the interferogram
% rows - the y pixel values of the interferogram

columns = coordsTilt(2): coordsTilt(2)+coordsTilt(4);
rows = coordsTilt(1): coordsTilt(1)+coordsTilt(3);

regionTilt = phasetemp(columns,rows);


%Spatially unwraps the region designated for tilt
tiltplane = puma_ho(regionTilt,2);
 
% Fit 2-D tilt to background data
[x1,y1] = meshgrid(columns,rows);
model_function0 = @(vars,r) vars(3) + vars(1)*r(:,1) + vars(2)*r(:,2); % fit background tilt to a plane 
guess2 = [(tiltplane(round(size(tiltplane,1)/2),end)-tiltplane(round(size(tiltplane,1)/2),1))/(columns(end)-columns(1)) (tiltplane(end,round(size(tiltplane,2)/2))-tiltplane(1,round(size(tiltplane,2)/2)))/(rows(end)-rows(1)) mean(mean(tiltplane))]; % guess is taking endpoints and finding slope
[varfit2,resid,J,Sigma] = nlinfit([x1(:) y1(:)],tiltplane(:),model_function0,guess2,statset('TolFun',1e-30,'TolX',1e-20,'MaxIter',1000));
        
xtilt = varfit2(1);        
ytilt = varfit2(2);

