%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GeneratePhaseMap v .1 Jan/11/2016            
% Bill Green, James Corsetti          
% Description: Continuously creates and saves phasemaps
% until the dialog box is terminated.

% Last Updated: 9/20/2016
% Changed the timestamp to a datevector, which the Recordphasemap function
% will write to a text file




function [phasemap,timestamp] = GeneratePhaseMap(PSramp,htrans,vidobj,tag,pollrate)


% This section sets up the time intervals that take place inside the
% phasemap generator. We want to capture the phasemap as quickly as
% possible, then wait the remainder of the time to get an equally spaced
% sampling rate.

%This is the time that the stage halts, to try to attenuate any transient
%motion before the snapshot is taken. It should never be shorter than the
%camera framerate.
steppause = .02;

 %Take a series of snapshots used to create a phasemap
 tic
 pause(.1)
    for ijk = 1:length(PSramp) % Number of voltage steps
        htrans.SetVoltOutput(0,PSramp(ijk));
        pause(steppause)
%         for lmn = 1 % number of images you are averaging
%             onefrm_stp = getsnapshot(vidobj); % store a single frame
%             allfrm_stp(:,:,lmn)=onefrm_stp; % Frames to be averaged
%             clear onefrm_stp  
%         end
%         avfrm_stp=mean(allfrm_stp,3);
        snap(:,:,ijk) = getsnapshot(vidobj);
%         clear avfrm_stp
    end
    
    
    %Create a timestamp for the phasemap
    clk = clock; 
    timestamp = [tag '_' num2str(clk(4)) 'h' num2str(clk(5)) 'm_' num2str(clk(2)) '_' num2str(clk(3)) '_' num2str(clk(1))];
    
    %Malacara Algorithm to generate a phasmap
    del = (1:length(PSramp)).*2.*pi./(length(PSramp)-1);
    I = double(snap);%(:,:,1:PSramp);
    N = 0;
    D = 0;
    for qrs = 1:(length(del))
        N = N - I(:,:,qrs).*sin(del(qrs));
        N=double(N);
        D = D + I(:,:,qrs).*cos(del(qrs));
        D=double(D);
    end
    phasemap = atan2(N,D);
    
    %This will hold the program until the full second is over.
    duration = toc
    pause(pollrate-duration);  