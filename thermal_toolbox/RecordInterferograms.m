function RecordPhaseMap(PSramp,htrans,vidobj,fileloc,pollrate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RecordPhaseMap v .1 Jan/11/2016            
% Bill Green, James Corsetti          
% Description: Continuously creates and saves phasemaps
% until the dialog box is terminated.

% PSramp - the voltage vector that was correlated to a 2pi phase shift
% htrans - the handle for the piezo stage
% vidobj - the handle for the camera
% fileloc - where the phasemaps will be save
% pollrate - the amount of time in seconds that a phase measurement will
% take. WARNING: the maximum value is the framerate divided by the number
% of voltage steps used in the calibration.

if nargin < 5, pollrate = []; end
% supply default parameters
if isempty(pollrate), pollrate = 1; end

snapshotinterval = 1/(pollrate*length(PSramp));
originalDirectory = pwd;
cd(fileloc)


% Creates a file which records the datenum of each phasemap
timefile = fopen('00000_timestamps.txt','a+');
fmt = '%f \n';




FS = stoploop;
PhaseMap = 0;
while ~FS.Stop()
    startburst = clock
    %_1 Name the count for the phase map
    PhaseMap = PhaseMap + 1;
    if PhaseMap < 10;
        tag = ['0000' num2str(PhaseMap)];
    elseif (PhaseMap > 9) && (PhaseMap < 100)
        tag = ['000' num2str(PhaseMap)];
    elseif (PhaseMap > 99) && (PhaseMap < 1000)
        tag = ['00' num2str(PhaseMap)];
    elseif (PhaseMap > 999) && (PhaseMap < 10000)
        tag = ['0' num2str(PhaseMap)];
    else
        tag = num2str(PhaseMap);
    end
    
    PhaseMap
	
for i = 1:length(PSramp)
htrans.SetVoltOutput(0,PSramp(i));
phasetemp(:,:,i) = getsnapshot(vidobj);
end

fprintf(timefile,fmt,now);


save(tag,'phasetemp');
stopburst = clock
Rtime = 1 - etime(stopburst,startburst)
pause(Rtime)
end

fclose(timefile);
cd(originalDirectory)