%% Tilt Extraction
fileloc = 'I:\2017_05_18_Saph_p8_to18ramp_4hrs\PhaseMaps';
cd(fileloc);
clear; clc; close all;
path(path, 'I:\thermal_toolbox'); 
load('060917wkspace.mat');


for frame = f1:10:f2;%f2; %4:10;%26493%length(idx)
    frame
    PM = D(frame).name;
    load(PM);
    
    b = puma_ho(phasetemp(backgroundrows,backgroundcols),2);
    % Fit 2-D tilt to background data
    [x1,y1] = meshgrid(backgroundcols,backgroundrows);
    model_function0 = @(vars,r) vars(3) + vars(1)*r(:,1) + vars(2)*r(:,2); % fit background tilt to a plane 
    guess = [(b(round(size(b,1)/2),end)-b(round(size(b,1)/2),1))/(backgroundcols(end)-backgroundcols(1)) (b(end,round(size(b,2)/2))-b(1,round(size(b,2)/2)))/(backgroundrows(end)-backgroundrows(1)) mean(mean(b))]; % guess is taking endpoints and finding slope
    [varfit1,resid,J,Sigma] = nlinfit([x1(:) y1(:)],b(:),model_function0,guess,statset('TolFun',1e-30,'TolX',1e-20,'MaxIter',1000));
    xtilt(frame-f1+1) = varfit1(1);%/(2*nair(frame-f1+1)*k); % putting the tilt into units of waves
    ytilt(frame-f1+1) = varfit1(2);%/(2*nair(frame-f1+1)*k);
end
cd ..
save('wkspace_from_2014.mat');

