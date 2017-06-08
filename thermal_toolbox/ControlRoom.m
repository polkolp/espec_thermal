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
%% Initial Setup

clear all

% This is the location of the toolbox on the computer
cd('C:\Users\Chopin\OneDrive\espec_thermal\thermal_toolbox')
path(path,'C:\Users\Chopin\OneDrive\espec_thermal\thermal_toolbox') %add analysis functions in path

% This setting determines the kind of test. A value of 1 will ask for
% information about the corrsponding region. Ex. Set FlagTransmit = 0 if
% you are only interested in the CTE of a sample.


FlagReflect = 0;
FlagTransmit = 0;
FlagGRIN = 0;

% Do you want to look at the wedge in the sample?
FlagWedge = 0;

%Important information about the test and sample
waveLength = .6328; %micronsh
pressure = 1000.*98.9259776400000; %pascals
relativeHumidity = 26.5; %percent
thickNess = 5280; %micronsc
inDexSample = 1.4935; %Sample index at ambient temperature
wavenumber = 2*pi/waveLength;

%Start up the camera! This starts up the camera with a target framerate of
%60 fps, and an exposure time of 16 ms.
[src vidobj] = InitCamera(60);

%Start up the Stage! This creates an object which sends commands to the
%piezocontroller.
htrans = InitStage;

%Calibrate the stage!
%% Ensure that the selected pixel for calibration is representative of the background
PSramp = CalibrateStage(htrans,vidobj,[50,50],20);


%% Crop the region that is being watched by the camera! 
% A figure will appear, drag and right click to select the region you would
% like to record interferograms for. Cropping is important because it
% allows you to store smaller size phasemaps and reduce the size of the raw
% data.

vidobj.ROIPosition = DefineRegion(vidobj,2);

%% Runs the loop to save phase maps
clear bgcnt D DataBackgroundNormalized DataBackgroundPhase DataBackgroundUnwrapped datadir dotsize i idx j loop_progress numPhaseMaps OPLBackground phaseMetaData phasetemp t Time vv
%Where do you want to store your data?

datadir = 'F:\02_09_2016_testFourier\PhaseMaps';
mkdir(datadir)

RecordPhaseMap(PSramp,htrans,vidobj,datadir,2)

%% This is the end of the data acquisition section of the script
%%%%%%
%%%%%%
%%%%%%
%%%%%%

%% Parses the selected directory to make phasemaps available for analysis

clear D vv idx
datadir = 'F:\02_09_2016_testFourier\PhaseMaps';
D = dir(datadir);
[vv,idx] = sort([D.datenum]);

%Information from the files
phaseMetaData = dir(datadir);
numPhaseMaps = length(find(~[D.isdir]))

%Loads the first phase map so that its image can be displayed
cd(datadir)
load(D(3).name)

% The DataSet structure is useful for storing information about the test,
% along with the processed data from it

DataSet = {};
DataSet.name = '5_01_16_SNR';
DataSet.location = datadir;
DataSet.numphasemaps = numPhaseMaps;
DataSet.Reflect = FlagReflect;
DataSet.Transmit = FlagTransmit;
DataSet.Thickness = thickness;
DataSet.Wavelength = waveLength;




%% Select the distinct measurement regions

disp('Currently Selecting Background Region...')
coordsBackground = DefineRegion(phasetemp,1)
disp('Background Selection Complete.')

coordsTilt = [round(1.2*coordsBackground(1)),round(1.2*coordsBackground(2)),...
    100,100];

if FlagReflect == 1
disp('Currently Selecting Reflected Region...');
coordsReflect = DefineRegion(phasetemp,1);
disp('Reflectance Selection Complete.');
end
 
if FlagTransmit == 1
disp('Currently Selecting Transmitted Region...')
coordsTransmit = DefineRegion(phasetemp,1)
disp('Transmittance Selection Complete.')
end

%Every tenth pixel will be used to determine fringe movement
meshsize = 10;

DataDisplay = figure(45);
figure(DataDisplay)
imagesc(phasetemp);
colormap summer
axis equal
hold on

%Plotting the regions for tilt and wedge
rectangle('Position',coordsTilt,'LineWidth',2);

if FlagWedge == 1
    rectangle('Position',coordsWedge,'LineWidth',2);
end

%Graphically displays the cropped region in the figure
figure(DataDisplay)
dotsize = 4;

regionBackground = SubSample(coordsBackground,meshsize);
for j = 1:size(regionBackground,1)
    hold on
    plot(regionBackground(j,1),regionBackground(j,2),'ro','MarkerSize',dotsize,'MarkerEdgeColor','k','MarkerFaceColor',[ j./size(regionBackground,1) 0 0])
end

% Homogenous reflection region if sample is not a GRIN
if FlagReflect == 1
    regionReflect = SubSample(coordsReflect,meshsize);
    for j = 1:size(regionReflect,1)
        hold on
        plot(regionReflect(j,1),regionReflect(j,2),'ro','MarkerSize',dotsize,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 j./size(regionReflect,1)])
    end
end

% Transmission for a homogenous sample
if FlagTransmit == 1
    regionTransmit = SubSample(coordsTransmit,meshsize);
    for j = 1:size(regionTransmit,1)
        hold on
        plot(regionTransmit(j,1),regionTransmit(j,2),'ro','MarkerSize',dotsize,'MarkerEdgeColor','k','MarkerFaceColor',[0 j./size(regionTransmit,1) 0])
    end
end

saveas(gcf,'datasetup','png')



%% Code for loading all temperature and phase data from saved files, and creating measurables from it 

tic
cd(datadir)
DataBackgroundPhase = zeros(numPhaseMaps,size(regionBackground,1));
if FlagReflect ==1; DataReflectPhase = zeros(numPhaseMaps,size(regionReflect,1));end
if FlagTransmit ==1; DataTransmitPhase = zeros(numPhaseMaps,size(regionTransmit,1));end



for i = 1:numPhaseMaps
    load(D(i+2).name)
    
    %Pulls the datapoints specified by the subsampling of the background
    for bgcnt = 1:size(regionBackground,1);
        DataBackgroundPhase(i,bgcnt)=phasetemp(regionBackground(bgcnt,2)+1,regionBackground(bgcnt,1)+1);
    end
    
    %Pulls the datapoints for homogenous CTE
    if FlagReflect == 1
        for rfcnt = 1:size(regionReflect,1);
            DataReflectPhase(i,rfcnt)=phasetemp(regionReflect(rfcnt,2)+1,regionReflect(rfcnt,1)+1);
        end
    end
    
    %Pulls the datapoints for homogeneous Transmission
    if FlagTransmit == 1
        for tmcnt = 1:size(regionTransmit,1);
            DataTransmitPhase(i,tmcnt)=phasetemp(regionTransmit(tmcnt,2)+1,regionTransmit(tmcnt,1)+1);
        end
    end
    
    %Runs the puma-ho algorithm to determine the background tilt changes
    %over time - only every tenth data point is required.
%     if ~mod(i,10)
%     [xtilt(i/10),ytilt(i/10)] = CalculateTilt(phasetemp,coordsTilt); %microns/pixel
%     end
%     
     %Runs the puma-ho algorithm to determine wedge in the sample
    if FlagWedge == 1
        [xwedge(i),ywedge(i)] = CalculateTilt(phasetemp,coordsWedge); %microns/pixel
    end
    
    
    %Arranges into a matrix the timestamp for each phase measumement
    t(i,:) = datevec(phaseMetaData(i+2).datenum);
    Time(i) = etime(t(i,:),t(1,:));
    loop_progress = 100*i/numPhaseMaps
end
toc

DataSet.Time = Time';

%% Grabs the temperature thermocouple measurements
temploc = 'E:\Documents\Research\Ellis\MGRIN\TestBed\CGA_V1\4_28_16_CaF2\CaF2ColdTemps.xls';

tref = t(1,:); %The time at which the test began
DataSet.Temps = GrabTemp(Time,tref,temploc);

% Specify which thermocouple channels represent what quantities

% Warning: This section is sketchy! You want to end up with a n x 12
% matrix, where n is the same as the number of phasemaps. This should
% output the interpolated temperature data.


%% Vetting all the data in order to calculate OPD

% tstart is the time in seconds for which you want to pin all your phase
% measurements to zero. For the CGA V1, this can safely be set to a number
% between 0 and 400 seconds.
tstart = 0;


DataBackgroundUnwrapped = unwrap(DataBackgroundPhase,[],1);
if FlagReflect == 1;DataReflectUnwrapped = unwrap(DataReflectPhase,[],1);end
if FlagTransmit == 1;DataTransmitUnwrapped = unwrap(DataTransmitPhase,[],1); end

%Designate tstart as the beginning of meaningful data, and set the phase in
%all the data to be 0 at this point. 


startIndex = find(abs(Time-tstart) < 2);
startIndex = startIndex(1);

% These loops pin the optical path difference to be zero at the time
% selected above.
for j=1:size(DataBackgroundUnwrapped,2)
DataBackgroundNormalized(:,j) = DataBackgroundUnwrapped(:,j) - DataBackgroundUnwrapped(startIndex,j);
end

if FlagReflect == 1 %Homogenous Reflection
    for k=1:size(DataReflectUnwrapped,2)
        DataReflectNormalized(:,k) = DataReflectUnwrapped(:,k) - DataReflectUnwrapped(startIndex,k);
    end
end

if FlagTransmit == 1 %Homogenous Trasmission
    for k=1:size(DataTransmitUnwrapped,2)
        DataTransmitNormalized(:,k) = DataTransmitUnwrapped(:,k) - DataTransmitUnwrapped(startIndex,k);
    end
end

%% Calculations to get the OPD

OPLBackground = mean(DataBackgroundNormalized,2)./wavenumber;
if FlagReflect == 1;OPLReflect = mean(DataReflectNormalized,2)./wavenumber;end
if FlagTransmit == 1;OPLTransmit = mean(DataTransmitNormalized,2)./wavenumber;end

DataSet.BOPL = OPLBackground
DataSet.ROPL = OPLReflect
DataSet.TOPL = OPLTransmit
DataSet.Temps = Temperatures;

%% Plotting of time series data

figure;plot(DataSet.Time/3600,[DataSet.BOPL,DataSet.ROPL,DataSet.TOPL]);legend('B','R','T')
xlabel('Time (hours)')
ylabel('Physical Length Change (\mum)')
%DataSet.Info = ; %Test length, temperatures, comments, sample type, thermocouples, directories, settings
figure;plot(DataSet.Time/3600,[-(DataSet.ROPL - DataSet.BOPL)./2]); legend('R(t)')
xlabel('Time (hours)')
ylabel('Expansion of Sample (\mum)')
figure;plot(DataSet.Time/3600,(DataSet.TOPL- DataSet.ROPL)./(DataSet.BOPL - DataSet.ROPL))
xlabel('Time (hours)')
ylabel('Index Change of Sample')
%%


%inDexAir = air_index_calc(DataSet.Temps(:,12),pressure,relativeHumidity,waveLength);

figure;
%plot(Time/3600,OPLBackground)
[ax,p1,p2] = plotyy(DataSet.Time/3600,DataSet.BOPL,DataSet.Time/3600,DataSet.Temps(:,11),'plot','plot')
legend('Background Path','Aluminum Plate Temperature')
xlabel(ax(1),'Time (hours)') % label x-axis
ylabel(ax(1),'Optical Path Change (\mum)') % label left y-axis
ylabel(ax(2),'Temperature (C)') % label right y-axis
grid on

% figure;
% [ax,p1,p2] = plotyy(Time/3600,OPLBackground,...
%     [Time/3600,Time/3600],[DataSet.Temps(:,1),DataSet.Temps(:,6)],'plot','plot');

figure
set(gca,'fontsize',18)
hold on
h = plot(DataSet.Time/3600,DataSet.Temps(:,1),'color',[1 0 0],'linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,2),'color',[0 1 0],'linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,3),'color',[0 0 1],'linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,4),'color',[1 1 0],'linestyle','--','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,5),'color',[0 1 1],'linestyle','--','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,6),'color',[1 0 1],'linestyle','--','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,7),'color',[0 0 0],'linestyle','--','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,8),'color',[0 1 0],'linestyle','-','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,9),'color',[.5 .5 1],'linestyle','-','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,10),'color',[.2 .4 .6],'linestyle',':','linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,11),'color',[.6 .3 .2],'linewidth',1)
plot(DataSet.Time/3600,DataSet.Temps(:,12),'color',[.1 .8 .5],'linewidth',1)
grid on

title('COLD CTE, ThermoTron Convection Test')
xlabel('time (h)')
ylabel('Temperature')
legend('Mirror1','Mirror2','Mirror3','Air at Mirror','Inlet','Outlet','Inner Chamber','Outer Chamber','Bottom Window','Top of Window','Bottom Plate','Ambient')
set(h,'xlim',[0,3])
set(h,'ylim',[25,-30])
grid on

%% Temperature Stability Analysis - This is an example script for looking at the steady state temperature variation with time.
DriftTemps = DataSet.Temps(7500:8010,:);
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


%% Getting the CTE out from the Temp vs. OPL curve
figure
plot(DataSet.Temps(:,2),-(DataSet.ROPL - DataSet.BOPL)./2)
coefficients = polyfit(DataSet.Temps(:,2),-(DataSet.ROPL - DataSet.BOPL)./2,2)
hold on
plot(DataSet.Temps(:,2),-1.88e-4.*DataSet.Temps(:,2) +  .1324.*  DataSet.Temps(:,2) - 2.6858)
xlabel('Temperature \circ C')
ylabel('Physical Length Change \mum')
legend('Measured Data','Polynomial fit')
CTE = coefficients(2)/.00528
grid on




