
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GrabTemp v .1 Jan/4/2016            
% Bill Green, James Corsetti          
% Description: Acquires thermocouple data from the SD card, and
% interpolates it in order to line up with the phase timestamps that are
% provided.

% Features to add: instead of Startplace, have this function find the
% correct temp data automatically based on the timestamps.

startPlace = 1;  %Designate the row for which temp data begins
 
logFile = fopen('E:\TMB01\TMB01028.xls');
dateFormat = 'yyyy/mm/ddHH:MM:SS'

%This is the time that will be used to indicate the start of the test (T=0)
tref = [2015,12,23,09,45,30];
phaseTimes = linspace(0,5000);

for i=1:startPlace
    advanceset = fgets(logFile);
end
    %Loginfo parses each line of the thermocouple data, delimited by tabs, into the time of the
    %measurement and the value for each thermocouple. 
    
    % Each '%f' is a thermocouple reading. 
    % Each '%s' records the string that contains the time. 
    % A formatspec such as '%*s' ignores the string that it parses.
    
    logInfo = textscan(logFile,'%*f %s %s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f %*s %*s %f')
    numTempPts = length(logInfo{1,1}); %Number of data points collected
    
    %Finds the time difference between the start of phase tracking and each
    %temperature point
    for j = 1:numTempPts;
        logTime(j,:) = datevec(strcat(logInfo{1,1}{j,1},logInfo{1,2}{j,1}),dateFormat);
        
        % Relative to T=0, this is the time that each temperature
        % measurement occured
        timeDiff(j) = etime(logTime(j,:),tref);
    end
   
    numLogEntries = size(logInfo);
    
    % Based on the way the .xls file is written, this will arrange
    % temperature measurements from each thermocouple into a matrix. The
    % nth row of temp matches T(n) where T is the time that the measurement
    % was taken. 
   
    for k = 3:numLogEntries(2)-1
    temp(:,k-2)= interp1(timeDiff,logInfo{1,k},phaseTimes,'pchip');
    end
    fclose(logFile)