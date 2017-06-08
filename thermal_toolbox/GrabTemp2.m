function thermaltable = GrabTemp2(temploc,channels)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Created by Tianyi Yang 051617
% read temperature and time data from Omega measurements
% Modified from GrabTemp by Bill Green and James Corsetti
% numprobe = number of probes used
%=============original readme ====================
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

%% Check to see which thermocouples are active
thermaldate = [table.Date];
thermaltime = [table.Time];

datetimecol = datenum([thermaldate thermaltime]);
if channels == 4
      pb1 = [table.Value];%probe readings
      pb2 = [table.Value1];
      pb3 = [table.Value2];
      pb4 = [table.Value3]; 
      thermaltable = [datetimecol, pb1, pb2, pb3, pb4];
elseif channels == 12
      pb1 = [table.Value0];%probe readings
      pb2 = [table.Value1];
      pb3 = [table.Value2];
      pb4 = [table.Value3]; 
      pb5 = [table.Value4];
      pb6 = [table.Value5];
      pb7 = [table.Value6];
      pb8 = [table.Value7]; 
      pb9 = [table.Value8];
      pb10 = [table.Value9];
      pb11 = [table.Value10];
      pb12 = [table.Value11]; 
      thermaltable = [datetimecol, pb1, pb2, pb3, pb4, pb5, pb6, pb7, pb8, pb9, pb10, pb11, pb12];
end
end