function [src vidobj] = InitCamera(framerate,exposuretime)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% InitCamera v .1 Apr/29/2016            
% Bill Green, James Corsetti          
% Description: Starts up the CMOS sensor used for capturing phasemaps

% src - a structure that contains video settings such as frame rate,
% exposure time, etc.

% vidobj - a structure that references the video feed

% framerate - provided in FPS. This script will provide a default if not
% specified

% exposuretime - prodivded in seconds. A default is provided if not specified. 
% This script will terminate if the exposure time is longer than the time between frames.
imaqreset

if nargin < 2, exposuretime = []; end
if nargin < 1, framerate = []; end


% supply default parameters
if isempty(framerate), framerate = 30; end
if isempty(exposuretime), exposuretime = .03000; end


if exist('vidobj')
delete(vidobj)
clear vidobj
end

vidobj=videoinput('pointgrey', 1, 'F7_Raw8_640x512_Mode1');
viewport = preview(vidobj);
vidobj.ReturnedColorspace='grayscale';
src = getselectedsource(vidobj);
vidobj.ROIPosition = [0 0 640 512];


disp('Please wait 10 seconds...')

pause(2)
src.ExposureMode = 'Manual';
pause(2)
maxfps = 120;
src.FrameRatePercentageMode = 'Manual';
src.FrameRatePercentage = 100*framerate/maxfps;
pause(2)
src.GainMode = 'Manual';
pause(2)
src.ShutterMode = 'Manual';
if exposuretime > 1/framerate
    exposuretime = 1/(1.1*framerate);
    disp('Exposure time exceeds time inbetween frames, shutterspeed changed accordingly')
end
%src.Shutter = exposuretime;
pause(2)
disp('Camera Initialized Successfully')