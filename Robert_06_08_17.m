%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ThermalInterferometer Control Room v .1 Jan/15/2016            
% Bill Green, James Corsetti          
% 
% Description: This script arranges, for the first time user, all the
% functions necessary to operate the thermal interferometer for a variety
% of samples. All functions are commented, please read these while
% troubleshooting, and have fun!

% Features still wanted: 
% -Calculate OPD considering the spatial index variation of the air
% -Error analysis for many of the measurables
% -Interpolation of the puma-ho algorithm to increase processing speed
% -Better implementation of temperature data integration.
% -Integrate James' low pass fourier filter
% -Implement a function that will allow more than one sample to be measured
% at the same time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revised by Tianyi Yang 051517
% purpose: Adapt to be used in espec
% changes: 
% - packaged data acquisition into separate file. This file only works on
% analysis
%
%%%%%%%%%%%%%%051617%%%%%%%%%%%%%%%%%%%%%%%%
% The folder containing interferogram should contain ONLY interferograms,
% not any other subfolders for the code to work.
%
%
%
%%%%%%%%%%%%%%%052517%%%%%%%%%%%%%%%%%%%%%%
% -put the meshsize part into initializtion 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% path = 'E:\Documents\Research\Ellis\MGRIN\Data\';


% %% Sets up the directories you will use for data
% clear bgcnt D datadir i idx j loop_progress numPhaseMaps OPLBackground phaseMetaData phasetemp t Time vv
% %Where do you want to store your data?
% 
% nameOfTest = 'Zerodur_try3';
% 
% when = clock;
% datadir = strcat(path,num2str(when(1)),'_',num2str(when(2)),'_',num2str(when(3)),'_',nameOfTest)
% mkdir(datadir)
% cd(datadir)
% mkdir('phasemaps')
% phasedir = strcat(datadir,'\phasemaps')

% %% Crop the region that is being watched by the camera! 
% % A figure will appear, drag and right click to select the region you would
% % like to record interferograms for. Cropping is important because it
% % allows you to store smaller size phasemaps and reduce the size of the raw
% % data.
% 
% vidobj.ROIPosition = DefineRegion(vidobj,2);

%% This is where all the inititialization variables are defined
clear; clc; close all;
path(path,'C:\Users\Chopin\OneDrive\espec_thermal\thermal_toolbox') %add analysis functions in path
temploc = 'F:\TMB01033.XLS'; %temperature file location of Omega measurements
if (exist('TempLog.txt')) 
    delete('TempLog.txt');
end

% This setting determines the kind of test. A envalue of 1 will ask for
% information about the corrsponding region. Ex. Set FlagTransmit = 0 if
% you are only interested in the CTE of a sample.


FlagReflect = 1; FlagTransmit = 1; FlagGRIN = 0;
% Do you want to look at the wedge in the sample?
FlagWedge = 0;

%Important information about the test and sample
waveLength = .6328; %microns
pressure = 1000.*98.9259776400000; %pascals
relativeHumidity = 26.5; %percent
thickness = 5000; %microns
inDexSample = 1.7659; %Sample index at ambient temperature
wavenumber = 2*pi/waveLength;
% The meshsize will determine what pixel density you want to collect from
% the data. For example, setting a meshsize of 3 will collect the following
% pixels:

% o x x o x x o
% x x x x x x x
% x x x x x x x
% o x x o x x o
% x x x x x x x
% x x x x x x x
% o x x o x x o

% Each o represents a pixel that will be used for analysis. The meshsize of
% 3 means that 1/9 of the total pixels are used.
Bmesh = 1;
Rmesh = 1;
Tmesh = 1;


Ptab = 1000.*98.9259776400000; %pressure
RHtab = 26.5; %relative humidity  

%% If you are collecting and analyzing the data in different sessions, use
% this section to browse for the phasemaps you are interested in analyzing.
% You need to go into the "Phasemap folder" and click OK.
disp('Select Phasemap Folder ...')
phasedir = uigetdir 

%% Parses the selected directory to make phasemaps available for analysis

%Information from the files
clear D vv idx phasecube
cd(phasedir)
Ph = dir(phasedir); %Ph stands for raw phasemap data structure.
%load('test.mat');% This is only for testing, when the time information is needed to be imported from other file,you need to go to the folder which has this file


[vv,idx] = sort([Ph.datenum]); %just in case, sort files in order of acquisition time
phaseMetaData = dir(phasedir);
numPhaseMaps = length(find(~[Ph.isdir])); %%why minus 1? try not to

%Brings up the data selection
load(Ph(3).name) % why starting from 3? because the first two are name, '.' and '..')
%phasetemp= PhaseMapfromsnap(snap);

DataDisplay = figure(49);
figure(DataDisplay)
imagesc(phasetemp);
colormap summer
axis equal
hold on
title('Data Selection Regions')

% Pick up the region you want to analyze: Background , Reflection, Transmission
% Number of the files which are not the data files. May come from the
% system 
other = 0;
% Grab and display the background region
disp('Currently Selecting Background Region...')
coordsB = DefineRegion(phasetemp,1)
disp('Background Selection Complete.');pause(.4);

figure(DataDisplay);rectangle('Position',coordsB,'EdgeColor','b','LineWidth',3)
text(coordsB(1),(coordsB(2)+coordsB(4)/2),['Background' num2str(coordsB(3)) '*' num2str(coordsB(4))],'Color','blue','FontSize',14)

% Grab and display the reflected region
if FlagReflect == 1
    disp('Currently Selecting Reflected Region...')
    coordsR = DefineRegion(phasetemp,1)
    disp('Reflect Selection Complete.')
    
    figure(DataDisplay);rectangle('Position',coordsR,'EdgeColor','r','LineWidth',3);
    text(coordsR(1),(coordsR(2)+coordsR(4)/2),['Reflection' num2str(coordsR(3)) '*' num2str(coordsR(4))],'Color','red','FontSize',14)
end

% Grab and display the transmitted region
if FlagTransmit == 1
    disp('Currently Selecting Transmitted Region...')
    coordsT = DefineRegion(phasetemp,1)
    disp('Transmission Selection Complete.')
    
    figure(DataDisplay);rectangle('Position',coordsT,'EdgeColor','k','LineWidth',3);
    text(coordsT(1),(coordsT(2)+coordsT(4)/2),['Transmission' num2str(coordsT(3)) '*' num2str(coordsT(4))],'FontSize',14)
end

% Records the way that you cropped the data
cd ..
saveas(gcf,'datasetup','png')
cd(phasedir)

% Now we go and grab the pixels that we're interested in
clear Bcube Rcube Tcube

% Can we have     Bcube = phasetemp(coordsB(2):Bmesh:coordsB(2)+coordsB(4)-1 ...
%        ,   coordsB(1):Bmesh:coordsB(1)+coordsB(3)-1,:);



for i = 1:numPhaseMaps-other
    load(Ph(i+2).name); %load snap
    i %display progress
    %phasetemp = PhaseMapfromsnap(snap); %calculate phasemap

    Bcube(:,:,i) = phasetemp(coordsB(2):Bmesh:coordsB(2)+coordsB(4)-1 ...
        ,   coordsB(1):Bmesh:coordsB(1)+coordsB(3)-1);
    
    if FlagReflect ==1
        Rcube(:,:,i) = phasetemp(coordsR(2):Rmesh:coordsR(2)+coordsR(4)-1 ...
            ,   coordsR(1):Rmesh:coordsR(1)+coordsR(3)-1);
    end
    
    if FlagTransmit==1
        Tcube(:,:,i) = phasetemp(coordsT(2):Tmesh:coordsT(2)+coordsT(4)-1 ...
            ,   coordsT(1):Tmesh:coordsT(1)+coordsT(3)-1);
    end
    
end

%% Grab the time information from the phase

% tstart is the time in seconds for which you want to pin all your phase
% measurements to zero. For the CGA V1 test process, this can safely be set to a number
% between 0 and 400 seconds.

tstart = 0;

% This block retrieve the time of aquisition and calculate the elapsed time
% in seconds

Time = [Ph.datenum];
DS.RefTime = Time(3);
Time = (Time(3:end)-Time(3))*86400;


%Designate tstart as the beginning of meaningful data, and set the phase in
%all the data to be 0 at this point. 
startIndex = find(abs(Time-tstart) < 1);
startIndex = startIndex(1);
offset = Time(startIndex);
Time = Time - Time(startIndex);

% % The reference time records the beginning of the test
DS.Time = Time';

%% Unwrap all of the data (this should be replaced with the pumaho algorithm ?)
% and extract tilt
%please enter the mesh size which you want to use in fit
fit_mesh=10; %spacing between tilt data points

clear B_X B_Y B_fit R_X R_Y R_fit T_X T_Y T_fit Bcube_tilt Rcube_tilt Tcube_tilt 
clear fB fR fT
disp('processing Background')
% make the phasemap continuousm,same as following
BcubeUnwrapped = unwrap(unwrap(unwrap(Bcube),[],2),[],3);
% extract the tilt, we care about the phase difference,same as following
%Bcubezero = bsxfun(@minus,BcubeUnwrapped,BcubeUnwrapped(:,:,startIndex));
[Bn,Bm,~] = size(BcubeUnwrapped);
%set up the fit data
B_fit.data=BcubeUnwrapped(1:fit_mesh:end,1:fit_mesh:end,1);
[B_X,B_Y] = meshgrid(1:fit_mesh:Bm,1:fit_mesh:Bn);
B_fit.plot=[B_Y(:),B_X(:),B_fit.data(:)]
%fB: function after fitting;
%fotfB: structure include fitting stats
[fB,fotfB] = fit([B_fit.plot(:,1),B_fit.plot(:,2)],B_fit.plot(:,3),'poly11');  
for i=1:Bm
    Bcube_tilt(1:Bn,i)=fB(1:Bn,i);
end
%extract tilt!
Bcubezero = bsxfun(@minus,BcubeUnwrapped,Bcube_tilt);
DS.Background = Bcubezero;
disp('Background rsquare fit:')
disp(fotfB.rsquare)
disp('Background Done')

disp('processing Reflection zone')
if FlagReflect == 1
    RcubeUnwrapped = unwrap(unwrap(unwrap(Rcube),[],2),[],3);
    [Rn,Rm,~] = size(RcubeUnwrapped);
    R_fit.data=RcubeUnwrapped(1:fit_mesh:end,1:fit_mesh:end,1);
    [R_X,R_Y] = meshgrid(1:fit_mesh:Rm,1:fit_mesh:Rn);
    R_fit.plot=[R_Y(:),R_X(:),R_fit.data(:)]
    [fR,fotfR] = fit([R_fit.plot(:,1),R_fit.plot(:,2)],R_fit.plot(:,3),'poly11');
    for i=1:Rm
        Rcube_tilt(1:Rn,i)=fR(1:Rn,i);
    end
    % Sets the first entry in the cube to be 0
    Rcubezero = bsxfun(@minus,RcubeUnwrapped,Rcube_tilt);
    DS.Reflect = Rcubezero;
    
end
disp('Reflection zone rsquare fit:')
disp(fotfR.rsquare)
disp('Reflection zone done')

disp('processing Transmission zone')
if FlagTransmit == 1
    TcubeUnwrapped = unwrap(unwrap(unwrap(Tcube),[],2),[],3);
    [Tn,Tm,~] = size(TcubeUnwrapped);
    T_fit.data=TcubeUnwrapped(1:fit_mesh:end,1:fit_mesh:end,1);
    [T_X,T_Y] = meshgrid(1:fit_mesh:Tm,1:fit_mesh:Tn);
    T_fit.plot=[T_Y(:),T_X(:),T_fit.data(:)]
    [fT,fotfT] = fit([T_fit.plot(:,1),T_fit.plot(:,2)],T_fit.plot(:,3),'poly11');
    for i=1:Tm
        Tcube_tilt(1:Tn,i)=fT(1:Tn,i);
    end
    % Sets the first entry in the cube to be 0
    %Tcubezero = TcubeUnwrapped(:,:,startIndex)- Tcube_tilt;
    Tcubezero = bsxfun(@minus,TcubeUnwrapped,Tcube_tilt);
    DS.Transmit = Tcubezero;  
end
disp('Transmission zone rsquare fit:')
disp(fotfT.rsquare)
disp('Transmission zone done')
% The DataSet structure is useful for storing information about the test,
% along with the processed data from it. At the end of the session, you
% will only need to save this variable

DS.location = phasedir;
DS.numphasemaps = numPhaseMaps;

DS.Thickness = thickness;
DS.Wavelength = waveLength;
DS.Wavenumber = 2*pi/(DS.Wavelength);
disp('unwrap done')
%% Grabs the temperature thermocouple measurements

% This section retrives record time (in datenum form) and the temperature
% readings of corresponding probes used in measurement. The section plots
% the mentioned data for a sanity check. Interpolation will be done in
% later sections.

% WARNING: Ensure that the only data in the file is from the test you just
% performed. The section can only be used in 2 modes: 4 probes or 12 probes

xlsdata = GrabTemp2(temploc,4); %temploc is the xls file location of Omega measurements

Th.Time = xlsdata(1:41,1); %for testing purpose
Th.Temps = xlsdata(1:41,2:end);

% Th.Time = xlsdata(:,1); %Th stands for raw thermal data structure.
% Th.Temps = xlsdata(:,2:end);
disp(['end time difference [s]: ',num2str((Th.Time(end,1)-Th.Time(1,1))*86400-DS.Time(end))]);
rawtime = (Th.Time(:,1)-Th.Time(1,1))*24;

figure(44)
plot(rawtime,Th.Temps,'linewidth',2)
xlabel('Time, [h]')
ylabel('Temp, [\circC]')
set(gca,'fontsize', 18)
legend('Probe1','Probe2','Probe3','Probe4','Location','northwest')
grid on

cd ..
saveas(gcf,'temps','png')

cd(phasedir)

%% Interpolate temperature data
% In this section we interpolate the temperature data to match the rate of
% phase map acquisition. 

timetab = (Th.Time(:,1)-Th.Time(1,1))*86400; % The elapsed time in thermal measurements in second

% Based on the way the .xls file is written, this will arrange
% temperature measurements from each thermocouple into a matrix. The
% nth row of temp matches T(n) where T is the time that the measurement
% was taken.

temp_interp = interp1(timetab,Th.Temps,DS.Time);

% figure;
% plot(DS.Time,'linewidth',2)
% xlabel('Time, [h]')
% ylabel('Temp, [\circC]')
% set(gca,'fontsize', 18)
% legend('Probe1','Probe2','Probe3','Probe4','Location','northwest')
% grid on

TA_interped = temp_interp(:,1); % Ambient
TS_interped = (temp_interp(:,2)+temp_interp(:,3)+temp_interp(:,4))./3; % Sample
nairtab = air_index_calc(TA_interped,Ptab,RHtab,.6328);

SecTime=DS.Time/3600; % in hour
DS.Temps = TS_interped;

figure;
plot(SecTime,TA_interped,SecTime,TS_interped,'linewidth',2)
xlabel('Time, [h]')
ylabel('Temp, [\circC]')
set(gca,'fontsize', 18)
legend('Ambient','Sample','Location','northwest')
grid on




%% deviation of the background pixels over time, the "width" of the curve is the deviation of the background

% Figure 43 plots the deviation of the background pixels over time. If you ran a full
% test, do not use all of the pixels, because this will take a very long
% time to plot!

%you may now input the pixel pitch
pixel=10;

figure(43)
hold on
for i = 1:pixel:coordsB(3)
    for j = 1:pixel:coordsB(4)
        plot(squeeze(DS.Background(j,i,:))/DS.Wavenumber)
    end
end
xlabel('Time[s]')
ylabel('Total Optical Path Change [\mu m]')
grid on

%%

% Look at the random walk of the individual pixel measurements for optical
% path change. Remember that when you average the pixels to obtain one mean
% value for the optical path change, you are now dealing with the standard
% deviation of the mean. For a homogeneous sample, you are averaging over
% all the pixels. For a GRIN sample, you need to compute the average based
% on the size of the regions you are using.

% Calculate standard deviation of the mean instead of standard deviation
clear i
for i = 1:DS.numphasemaps-other
    Berror(i) = std2(DS.Background(:,:,i)/DS.Wavenumber)/sqrt(Bn*Bm);
    Rerror(i) = std2(DS.Reflect(:,:,i)/DS.Wavenumber)/sqrt(Rn*Rm);
    Terror(i) = std2(DS.Transmit(:,:,i)/DS.Wavenumber)/sqrt(Tn*Tm);
end




figure(42)
hold on
plot(1000*Berror,'linewidth',2)
plot(1000*Rerror,'linewidth',2)
plot(1000*Terror,'linewidth',2)
%plot(Terror,'linewidth',2)

grid on
legend('Background','Reflection','Transmission')

ylabel('Error in Region Measurement [nm]')
xlabel('Phasemap Number')
grid on

%% Averages the pixels together to form the path lengths that you will use

for i = 1:DS.numphasemaps-other
    DS.BOPD(i) = mean2(DS.Background(:,:,i))./DS.Wavenumber;
    if FlagReflect == 1;
        DS.ROPD(i) = mean2(DS.Reflect(:,:,i))./DS.Wavenumber;
    end
    if FlagTransmit == 1;
        DS.TOPD(i) = mean2(DS.Transmit(:,:,i))./DS.Wavenumber;
    end
end

%% 
% Calculate OPD_R - OPD_B = 2 * nair * L(T)
%This is where we subtract the optical paths to perform linear fitting
%to get the CTE. We need to plot L vs Temperature.

%tstart = 200;
%tend = 1500;


DS.dL = -.5*(DS.ROPD-DS.BOPD)'; %air?

%This is the thickness change of the sample in um.

figure(40)
plot(DS.Temps,DS.dL,'o','linewidth',2)
xlabel('Temperature [\circC]')
ylabel('Sample Thickness Change \mum')
grid on
%xlim([24.8,26.3])
%ylim([-.1,-.05])

[DS.ctefit,DS.ctefitstats] = polyfit(DS.Temps,DS.dL,1);
DS.CTE = DS.ctefit(1)/DS.Thickness; %Where is this factor of 6 hiding?
hold on
plot(DS.Temps,DS.Temps.*DS.ctefit(1)+ DS.ctefit(2),'r--','linewidth',2)

rsquared_CTE = 1 - DS.ctefitstats.normr^2 / norm(DS.dL-mean(DS.dL))^2;

legend('Measured Thickness Change',['Linear Fit, r^2 = ',num2str(rsquared_CTE)])
set(gca,'FontSize',14)


%% Error estimate on the CTE

%uncertainty: obtained from figure(42)
deltaL = 1e-6;

%diag of the covariance matrix, variance of the fit results
%b_err = sqrt(diag((DS.ctefitstats.R)\inv(DS.ctefitstats.R'))./DS.ctefitstats.normr.^2./DS.ctefitstats.df);
b_err = sqrt(diag((DS.ctefitstats.R)\inv(DS.ctefitstats.R')).*DS.ctefitstats.normr.^2./DS.ctefitstats.df);


deltadLdT = b_err(1);

%propagation of the uncertainty
%deltaalpha = sqrt( (1/thickness^2) *DS.CTE(1)* deltaL)^2 + ( (1/thickness) * deltadLdT)^2;
deltaaphaa = sqrt((1/DS.Thickness * DS.CTE * deltaL)^2 + ((1/thickness) * deltadLdT)^2);

%% This will plot all three paths as a function of time


figure;plot(DS.Time/3600,DS.BOPD,DS.Time/3600,DS.ROPD,DS.Time/3600,DS.TOPD,'linewidth',2);legend('OPD_1','OPD_2','OPD_3')
xlabel('Time (hours)')
grid on
ylabel('Optical Path Length Change (\mum)')
%DataSet.Info = ; %Test length, temperatures, comments, sample type, thermocouples, directories, settings
set(gca,'FontSize',14)
figure;plot(DS.Time/3600,[-(DS.ROPD - DS.BOPD)./2]);
hold on
plot(DS.Time/3600,[-(DS.ROPD - DS.BOPD)./2] + Berror + Rerror)
plot(DS.Time/3600,[-(DS.ROPD - DS.BOPD)./2] - Berror - Rerror)
legend('R(t)','R + \sigma','R - \sigma')
xlabel('Time (hours)')
ylabel('Expansion of Sample (\mum)')
figure;plot(DS.Time/3600,(DS.TOPD- DS.ROPD)./(DS.ROPD- DS.BOPD));legend('T(t)')
xlabel('Time (hours)')
ylabel('Index Change of Sample')


%% Temperature Stability Analysis - This is an example script for looking at the steady state temperature variation with time.
DriftTemps = DS.Temps(7500:8010,:);
MeanTemps = mean(DriftTemps,1)
for i = 1:12
    DriftTemps(:,i) = DriftTemps(:,i) - MeanTemps(i);
end
figure;
plot(DriftTemps(:,1:4))
xlabel('time (s)')
ylabel('Temperature Variation \circ C')
legend('Mirror1','Mirror2','Mirror3','Air at Mirror','Inlet','Outlet','Inner Chamber','Outer Chamber','Bottom Window','Top of Window','Bottom Plate','Ambient')

SampleTemp = mean(DriftTemps(:,1:4),2)
figure;
plot(SampleTemp)
eTemp = 3*std(SampleTemp)

%% Calculate the dn/dT
% calculate the OPD_T - OPD_R=2*n(T)*l(T)

DS.l=DS.Thickness*(1+ DS.CTE*(DS.Temps-DS.Temps(startIndex)))

DS.dn= .5* (DS.TOPD-DS.ROPD)./DS.l
figure(39)
plot(DS.Temps,DS.dn,'o','linewidth',2)
xlabel('Temperature [\circC]')
ylabel('Sample index Change')
grid on

[DS.dndTfit,DS.dndTfitstats] = polyfit(DS.Temps,DS.dn,1);
DS.dndT = DS.dndTfit; %Where is this factor of 6 hiding?
hold on
plot(DS.Temps,DS.Temps.*DS.dndTfit(1)+ DS.dndTfit(2),'r--','linewidth',2)

rsquared_dndT = 1 - DS.dndTfitstats.normr^2 / norm(DS.dn-mean(DS.dn))^2;

legend('Measured Index Change',['Linear Fit, r^2 = ',num2str(rsquared_dndT)])
set(gca,'FontSize',14)