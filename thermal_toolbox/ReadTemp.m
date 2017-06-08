function [temperatures,temptime] = ReadTemp(thermocouplefile)

logFile = fopen(thermocouplefile);
table = tdfread(thermocouplefile);
    
% %% Check to see where the relevant data begins
% 
numTempPts = length(table.Place)

%Creates date vectors for the data
%temptime = [table.Date,table.Time];
temptime = 0
if length(fieldnames(table)) == 11
    channels = 4;
elseif length(fieldnames(table)) == 27
    channels = 12;
end

% measurementCompare = -100
% i = 0
% while measurementCompare < -10
%     measurementCompare = etime(temptime(i+1,:),tref)
%     i = i + 1
% end


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
   % Takes the relevant thermocouple data, starting from the time
   % determined in the above loop
   j = 1;
    for k = 1:channels;
    if isActive(k)
       ActiveValue = eval(Valuenames{k});
       temperatures(:,j) = ActiveValue;
       j = j+1
    end
    end
    fclose('all');
end    