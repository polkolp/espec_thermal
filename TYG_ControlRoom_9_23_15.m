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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
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
temploc = 'F:\TMA29\TMA29030.xls'; %temperature file location of Omega measurements

% This setting determines the kind of test. A value of 1 will ask for
% information about the corrsponding region. Ex. Set FlagTransmit = 0 if
% you are only interested in the CTE of a sample.



FlagReflect = 1; FlagTransmit = 1; FlagGRIN = 0;
% Do you want to look at the wedge in the sample?
FlagWedge = 0;

%Important information about the test and sample
waveLength = .6328; %microns
pressure = 1000.*98.9259776400000; %pascals
relativeHumidity = 26.5; %percent
thickness = 14080; %microns
inDexSample = 1.4935; %Sample index at ambient temperature
wavenumber = 2*pi/waveLength;

%% If you are collecting and analyzing the data in different sessions, use
% this section to browse for the phasemaps you are interested in analyzing.
% You need to go into the "Phasemap folder" and click OK.
disp('Select Phasemap Folder ...')
phasedir = uigetdir 

%% Parses the selected directory to make phasemaps available for analysis

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

%Information from the files
clear D vv idx phasecube
Ph = dir(phasedir); %Ph stands for raw phasemap data structure.
cd(phasedir)
[vv,idx] = sort([Ph.datenum]); %just in case, sort files in order of acquisition time
phaseMetaData = dir(phasedir);
numPhaseMaps = length(find(~[Ph.isdir]))-1 %%why minus 1

%Brings up the data selection
load(Ph(4).name) % why starting from 4? because the first threes are name, '.' and '..')
phasetemp = PhaseMapfromsnap(snap);

DataDisplay = figure(45);
figure(DataDisplay)
imagesc(phasetemp);
colormap summer
axis equal
hold on
title('Data Selection Regions')

%%
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

%% Now we go and grab the pixels that we're interested in
clear Bcube Rcube Tcube

for i = 1:numPhaseMaps
    load(Ph(i+2).name); %load snap
    i %display progress
    phasetemp = PhaseMapfromsnap(snap); %calculate phasemap
    Bcube(:,:,i) = phasetemp(coordsB(2):Bmesh:coordsB(2)+coordsB(4)...
        ,   coordsB(1):Bmesh:coordsB(1)+coordsB(3));
    
    if FlagReflect ==1
        Rcube(:,:,i) = phasetemp(coordsR(2):Rmesh:coordsR(2)+coordsR(4)...
            ,   coordsR(1):Rmesh:coordsR(1)+coordsR(3));
    end
    
    if FlagTransmit==1
        Tcube(:,:,i) = phasetemp(coordsT(2):Tmesh:coordsT(2)+coordsT(4)...
            ,   coordsT(1):Tmesh:coordsT(1)+coordsT(3));
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

% The reference time records the beginning of the test
DS.Time = Time';

%% Unwrap all of the data (this should be replaced with the pumaho algorithm ?)

BcubeUnwrapped = unwrap(Bcube,[],3);
Bcubezero = bsxfun(@minus,BcubeUnwrapped,BcubeUnwrapped(:,:,startIndex));
DS.Background = Bcubezero;
[Bn,Bm,~] = size(DS.Background);

if FlagReflect == 1
    RcubeUnwrapped = unwrap(Rcube,[],3);
    % Sets the first entry in the cube to be 0
    Rcubezero = bsxfun(@minus,RcubeUnwrapped,RcubeUnwrapped(:,:,startIndex));
    DS.Reflect = Rcubezero;
    [Rn,Rm,~] = size(DS.Background);
end
if FlagTransmit == 1
    TcubeUnwrapped = unwrap(Tcube,[],3);
    Tcubezero = bsxfun(@minus,TcubeUnwrapped,TcubeUnwrapped(:,:,startIndex));
    DS.Transmit = Tcubezero;
    [Tn,Tm,~] = size(DS.Background);
end
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
Th.Time = xlsdata(:,1); %Th stands for raw thermal data structure.
Th.Temps = xlsdata(:,2:end);

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

Thtimeelapse = (Th.Time(:,1)-Th.Time(1,1))*86400; % The elapsed time in thermal measurements in second
% p1temp = Th.Temps(:,1);
% p2temp
% p3temp
% p4temp

% 
% offset = etime(TempBegin,tref)
% timeDiff = offset:pollrate:offset+pollrate*(numTempPts-1);

% Based on the way the .xls file is written, this will arrange
% temperature measurements from each thermocouple into a matrix. The
% nth row of temp matches T(n) where T is the time that the measurement
% was taken.

ptemp_interp = interp1(Thtimeelapse,Th.Temps,DS.Time);

figure;
plot(DS.Time,,'linewidth',2)
xlabel('Time, [h]')
ylabel('Temp, [\circC]')
set(gca,'fontsize', 18)
legend('Probe1','Probe2','Probe3','Probe4','Location','northwest')
grid on

% Bcube Rcube Tcube

f1 =3;%%3;%15201%13502%11501%09501%7501%4992%3201%0003;%%%%3001%0001
f2 =13343;%45003;%65263;%78248;%63;%size(time_PMS,2);%length(time_PMS);%16942%15200%13501%11500%9500%7500%4991%3200;%%%%6100%3000
tic
for frame = f1:1:f2;%f2; %4:10;%26493%length(idx)
    frame
    PM = Ph(frame).name;
%     load(PM);
    clear phasetemp_unfil phasetemp_fft phasetemp_combo
% %     time(frame-f1+1) = D(frame).datenum;
% % %         time(frame-f1+1) = time_PMS(frame);
        TAmeas(frame-f1+1) = interp1(timetab,TAtab,time(frame-f1+1));%,'cubic','extrap'); %interpolate temperature data (every 5s) to match phasemap acquisition rate(2 per s) 
        TSmeas(frame-f1+1) = interp1(timetab,TStab,time(frame-f1+1));%,'cubic','extrap');
%         TS_20mm_tab(frame-f1+1) = interp1(timetab,TS_20mm_tab,time(frame-f1+1));
        RHmeas(frame-f1+1) = RHtab;
        nair(frame-f1+1) = interp1(timetab,nairtab,time(frame-f1+1));%,'cubic','extrap');
end




% DS.Time;
% DS.Temps;
%% Calculations to get the OPD

% Figure 43 plots the deviation of the background pixels over time. If you ran a full
% test, do not use all of the pixels, because this will take a very long
% time to plot!

figure(43)
hold on
for i = 1:100
    for j = 50
plot(squeeze(DS.Background(j,i,:))*DS.Wavelength/(2*pi))
    end
end
xlabel('Time')
ylabel('Total Optical Path Change [\mu m]')
grid on

%%

% Look at the random walk of the individual pixel measurements for optical
% path change. Remember that when you average the pixels to obtain one mean
% value for the optical path change, you are now dealing with the standard
% deviation of the mean. For a homogeneous sample, you are averaging over
% all the pixels. For a GRIN sample, you need to compute the average based
% on the size of the regions you are using.


for i = 1:DS.numphasemaps
Berror(i) = 2*std2(DS.Background(:,:,i)*DS.Wavelength/(4*pi))/sqrt(Bn*Bm);
Rerror(i) = 2*std2(DS.Reflect(:,:,i)*DS.Wavelength/(4*pi))/sqrt(Rn*Rm);
%Terror(i) = 2*std2(DS.Transmit(:,:,i)*DS.Wavelength/(4*pi))/sqrt(Tn*Tm);
end


figure(42)
hold on
plot(1000*Berror,'linewidth',2)
plot(1000*Rerror,'linewidth',2)
%plot(Terror,'linewidth',2)
xlabel('Time [s]')
grid on
legend('Background','Reflection','Transmission')

ylabel('Error in Region Measurement [nm]')
xlabel('Time [s]')
grid on

%% Averages the pixels together to form the path lengths that you will use

for i = 1:DS.numphasemaps
    
DS.BOPD(i) = mean2(DS.Background(:,:,i))./DS.Wavenumber;
if FlagReflect == 1;DS.ROPD(i) = mean2(DS.Reflect(:,:,i))./DS.Wavenumber;end
if FlagTransmit == 1;DS.TOPD(i) = mean2(DS.Transmit(:,:,i))./DS.Wavenumber;end
end

%% This is where we subtract the optical paths to perform linear fitting to get the CTE. We need to plot vs Temperature.
tstart = 200;
tend = 1500;


DS.dL = -.5*(DS.ROPD(tstart:tend)-DS.BOPD(tstart:tend));%This is the thickness change of the sample in um.

figure(40)
plot(DS.Temps(tstart:tend,3),DS.dL,'linewidth',2)
xlabel('Temperature [\circC]')
ylabel('Sample Thickness Change \mum')
grid on
%xlim([24.8,26.3])
%ylim([-.1,-.05])

[DS.ctefit,DS.ctefitstats] = polyfit(DS.Temps(tstart:tend,3),DS.dL',1);
DS.CTE = DS.ctefit(1)/DS.Thickness %Where is this factor of 6 hiding?
hold on
plot(DS.Temps(tstart:tend,3),DS.Temps(tstart:tend,3).*DS.ctefit(1)+ DS.ctefit(2),'linewidth',2)

rsquared = 1 - DS.ctefitstats.normr^2 / norm(DS.dL-mean(DS.dL))^2;

legend('Measured Thickness Change',['Linear Fit, r^2 = ',num2str(rsquared)])
set(gca,'FontSize',14)


%% Error estimate on the CTE
thickness = 3.26e-3;
deltaL = 1e-6; %10 microns uncertainty
b_err = sqrt(diag((DS.ctefitstats.R)\inv(DS.ctefitstats.R'))./DS.ctefitstats.normr.^2./DS.ctefitstats.df);
deltadLdT = b_err(1); %uncertainty of the slope

deltaalpha = sqrt( (1/thickness^2) *DS.CTE(1)* deltaL)^2 + ( (1/thickness) * deltadLdT)^2);

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
legend('R(t)','R + 3\sigma','R - 3\sigma')
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


