function [PSramp] = CalibrateStage(htrans,vidobj,numsteps)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalibrateStage v .1 May/09/2016            
% Bill Green, James Corsetti          
% Description: Steps the stage through a few volts, corresponding to a few
% microns of reference mirror displacement, then, by looking at a 2pi
% displacement, correlates the voltage in the piezo-stage to the
% displacement of the reference mirror. Default value for the number of
% steps is 20 per phasemap.

% Last Updated Sep/20/2016
% Added lines that allow the user to select a pixel for use in calibration

% htrans - handle for the stage
% vidobj - handle for the video feed
% pxl    - the row and column of the pixel being used to perform the
%          calibration

% PSramp - a vector that contains a voltage range over which the mirror
% will move by a distance corresponding to a 2pi phase shift


if nargin < 3, numsteps = []; end
% supply default parameters
if isempty(numsteps), numsteps = 20; end

snap = getsnapshot(vidobj);
figure(22);imagesc(snap); axis equal
title('Please Select a point on the background for calibration')
pxl = floor(ginput(1));

close(22)

pause(1)
tic
clear vlt_stps V_app V_ramp clbrtn_frms
disp('Moving Stage to Obtain Calibration Sinusoid...')
for vlt_stps = 1:201; % Number of voltage 
    V_app = 4+.05.*vlt_stps./2; % applied voltage per step
    V_ramp(vlt_stps) = V_app; % generate voltage ramp
    htrans.SetVoltOutput(0,V_app);
    pause(.07)
    clbrtn_frms(:,:,vlt_stps) = getsnapshot(vidobj);
end
toc


%% Plot acquired linescan of a single pixel intensity vs. applied voltage
clear abc pxl_R pxl_C pxl_int amp1 shftd_pxl_int;


pxl_R = pxl(2); pxl_C = pxl(1);
sqr_sz = 2;
for abc = 1:length(V_ramp);
    pxl_int(abc) = mean(mean(clbrtn_frms(pxl_R:pxl_R+sqr_sz,pxl_C:pxl_C+sqr_sz,abc)));
end;

amp1 = 0.5.*(max(pxl_int) - min(pxl_int));
shftd_pxl_int = pxl_int - amp1 - min(pxl_int);

close(figure(3))
% figure(3)
% plot(V_ramp,shftd_pxl_int,'r.')
grid on

% close(figure(4));figure(4);
hold on; box on;
xx_data = V_ramp;
yy_data=shftd_pxl_int;
% plot(xx_data, yy_data, '.b')

d = xx_data;
Irrad = yy_data;
dmax = d(find(Irrad == max(Irrad),1,'first'));

Ifft=ifft(Irrad-mean(Irrad));
fftmax = find(abs(Ifft(1:fix(length(d)/2))) == max(abs(Ifft(1:fix(length(d)/2)))),1,'first');
dfx=(1/abs(d(2)-d(1)))/length(d);
fxguess = dfx*(fftmax-1)*2*pi;
phiguess = -fxguess*dmax;

AA = amp1;
BB = fxguess;
CC = phiguess;
DD = 0;
plot(xx_data, AA .* cos(BB .* xx_data + CC) + DD, '-g') % <-- Adjust these
title('Manual Fit')

% Now feed the starting point to Matlab
myfit = fittype('a*cos(b*x+c)+d', 'independent', 'x', 'dependent', 'y');
myopt = fitoptions('Method', 'NonlinearLeastSquares');
myopt.StartPoint = [AA BB CC DD]; % <-- Numbers from Line 9
[f, g] = fit(xx_data', yy_data', myfit, myopt);
yfit = f.a .* cos(f.b .* xx_data + f.c)  + f.d;

% Plot James' data and the regression result
close(figure(12));figure(12)
hold on; box on
plot(xx_data, yy_data, '.b')
plot(xx_data, yfit, '-r')
title('Calibration Fit')
xlabel('Voltage on Stage')
ylabel('Intensity')
legend('Data Points','Fitted Sinusoid','Stored Voltage Points')

V_fit = linspace(xx_data(1),xx_data(length(xx_data)),10001);
yy_fit = f.a .* cos(f.b .* V_fit + f.c)  + f.d;
[ee V_max] = max(yy_fit(1:round(length(yy_fit)./2.5)));
V1 = V_fit(V_max); V2 = V1 + 2.*pi./f.b;
hold on; box on;



% plot(linspace(V1,V2,num_st),linspace(0,0,num_st),'go')
% yyfit = f.a .* cos(f.b .* linspace(V1,V2,num_st) + f.c)  + f.d;
PSramp = linspace(V1,V2,numsteps+1);
PSramp = PSramp(1:numsteps);
plot(PSramp,f.a .* cos(f.b .* PSramp + f.c)  + f.d,'ko')
VoltageRange = PSramp(numsteps)-PSramp(1);
title(['Calibration Voltage: ',num2str(VoltageRange),' V.'])

disp(['The range that corresponds to' char(10) 'a 2pi phase shift is: ' ,num2str(VoltageRange), ' Volts.'])