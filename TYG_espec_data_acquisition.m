%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Espec Main data acquisition created on 051517 by Tianyi Yang
%===================================================================
% modified on absis of thermalstream control room v .1 Jan/15/2016 by Bill Green, James Corsetti           
%           
% 
% Description: This script arranges, for the first time user, all the
% functions necessary to operate the thermal interferometer for a variety
% of samples. All functions are commented, please read these while
% troubleshooting, and have fun!

%=====================================================
% Features still wanted: 




% ==============notes from original control room================
% -Calculate OPD considering the spatial index variation of the air
% -Error analysis for many of the measurables
% -Interpolation of the puma-ho algorithm to increase processing speed
% -Better implementation of temperature data integration.
% -Integrate James' low pass fourier filter
% -Implement a function that will allow more than one sample to be measured
% at the same time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Revised by Tianyi Yang 051117
% Adapt to be used in espec
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear; clc; close all; format long e
path = 'E:\Documents\Research\Ellis\MGRIN\Data\';

%% Initial Setup

clear all
% This is the location of the toolbox on the computer

%Start up the camera! This starts up the camera with a target framerate of
%60 fps, and an exposure time of 16 ms.
[src vidobj] = InitCamera(100);

%%

%Start up the Stage! This creates an object which sends commands to the
%piezocontroller.
htrans = InitStage;

%% Sets up the directories you will use for data
clear bgcnt D datadir i idx j loop_progress numPhaseMaps OPLBackground phaseMetaData phasetemp t Time vv
%Where do you want to store your data?

nameOfTest = 'Zerodur_try3';

when = clock;
datadir = strcat(path,num2str(when(1)),'_',num2str(when(2)),'_',num2str(when(3)),'_',nameOfTest)
mkdir(datadir)
cd(datadir)
mkdir('phasemaps')
phasedir = strcat(datadir,'\phasemaps')

%% Crop the region that is being watched by the camera! 
% A figure will appear, drag and right click to select the region you would
% like to record interferograms for. Cropping is important because it
% allows you to store smaller size phasemaps and reduce the size of the raw
% data.

vidobj.ROIPosition = DefineRegion(vidobj,2);
%% Contrast Check. Adjust the Polarizer if you see the brightness clipping.
check = getsnapshot(vidobj);
figure(21); subplot(2,1,1); histogram(check); xlim([0,255])
subplot(2,1,2); title('linescan'); plot(check(floor(end/2),:)); ylim([0,255]);
saveas(gcf,'ContrastStart','png');
%% Calibrate the stage!
% Ensure that the selected pixel for calibration is representative of the background
PSramp = CalibrateStage(htrans,vidobj,20);
saveas(gcf,'Calibration','png');

%%

%Edit this loop to set a timer before phasemap recording begins
% for i = 1:1800
%     pause(1)
%     1800-i
% end

% If you are wary of vibrations in the system, it is good practice to run a
% short test ~5 minutes, and check the standard deviation of the drift on
% the background. Use - std2(DataSet.Background(:,:,end))

RecordPhaseMap(PSramp,htrans,vidobj,phasedir,1)
%RecordInterferograms(PSramp,htrans,vidobj,phasedir,1)

%% This is the end of the data acquisition section of the script
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%
%%%%%%