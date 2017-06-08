function [DataTemperature] = GrabTemp(Time,tref,temploc,pollrate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GrabTemp v .1 Feb/12/2016
% Bill Green, James Corsetti
% Description: Acquires thermocouple data from the SD card, and
% interpolates it in order to line up with the phase timestamps that are
% provided.

% NOTE: In order for the function to perform properly, the data file must
% be manually coaxed a bit. Open the file in Notepad++, hit Ctrl+F, and
% replace all '/' with commas. Also replace all ':' with commas. This code
% will return an error if this has not been performed. Then, Shift-select
% to delete all the data that is not relevant to your test.



% We need to replace all the :'s and /'s with commas (,), so that tdfread
% can properly parse the file. You have no idea how long it took me to get
% MATLAB to do this.

TempLogFile = 'TempLog.txt';

%Performs the string replacement, and then creates a new text file
%containing all the information. The file is placed in the folder above the
%phasemaps.
data = fileread(temploc);
data = strrep(data, ':', ',');
data = strrep(data, '/', ',');
fid = fopen(TempLogFile, 'w');
fwrite(fid, data, 'char');
fclose(fid)

%Brings all the information into a structure for reading
table = tdfread(TempLogFile);

% Check to see where the relevant data begins
numTempPts = length(table.Place)

%Creates date vectors for the data
if length(fieldnames(table)) == 11
    channels = 4;
elseif length(fieldnames(table)) == 27
    channels = 12;
end

%% Check to see which thermocouples are active

if channels == 4
    isActive(1) = isnumeric(table.Value);
    isActive(2) = isnumeric(table.Value1);
    isActive(3) = isnumeric(table.Value2);
    isActive(4) = isnumeric(table.Value3);
    Valuenames = {'table.Value','table.Value1','table.Value2','table.Value3'}
elseif channels == 12
    isActive(1) = isnumeric(table.Value);
    isActive(2) = isnumeric(table.Value1);
    isActive(3) = isnumeric(table.Value2);
    isActive(4) = isnumeric(table.Value3);
    isActive(5) = isnumeric(table.Value4);
    isActive(6) = isnumeric(table.Value5);
    isActive(7) = isnumeric(table.Value6);
    isActive(8) = isnumeric(table.Value7);
    isActive(9) = isnumeric(table.Value8);
    isActive(10) = isnumeric(table.Value9);
    isActive(11) = isnumeric(table.Value10);
    isActive(12) = isnumeric(table.Value11);
end
Valuenames = {'table.Value','table.Value1','table.Value2','table.Value3','table.Value4','table.Value5','table.Value6',...
    'table.Value7','table.Value8','table.Value9','table.Value10','table.Value11'};

% Takes the relevant thermocouple data, checking for whether a
% thermocouple was hooked into the channel
j = 1;
for k = 1:channels;
    if isActive(k)
        ActiveValue = eval(Valuenames{k});
        temperatures(:,j) = ActiveValue;
        j = j+1
    end
end

% TempBegin is when you hit the button of the Omega Unit. We want to figure
% out when the temperature measurements occured relative to the phase
% measurements
TempBegin = [table.Date(1,:),table.Time(1,:)];
offset = etime(TempBegin,tref)
timeDiff = offset:pollrate:offset+pollrate*(numTempPts-1);

% Based on the way the .xls file is written, this will arrange
% temperature measurements from each thermocouple into a matrix. The
% nth row of temp matches T(n) where T is the time that the measurement
% was taken.

DataTemperature= interp1(timeDiff,temperatures,Time','spline');
end