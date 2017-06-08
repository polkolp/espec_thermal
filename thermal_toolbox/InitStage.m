function htrans = InitStage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InitStage v .1 Jan/11/2016            
% Bill Green, James Corsetti          
% Description: Starts up the Thorlabs piezometric stage

% htrans - a handle for the stage, used to set the voltage



disp('Initializing Piezometric stage...')
close(figure(2))
fig = figure(2); % Define figure for stage control activex GUI
% set(fig,'Position',[200 200 1100 400]);

htrans = actxcontrol('MGPIEZO.MGPiezoCtrl.1',[0 0 549 400],fig);% Define control for translation stage
%set(htrans,'HWSerialNum',29500337);% Serial Number for the CGA Stage
set(htrans,'HWSerialNum',81834010);% Serial Number for the eSpec Stage
htrans.StartCtrl;% Start control
disp('Phase Shifting Stage Initialized Successfully.')

