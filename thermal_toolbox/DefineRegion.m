function [cropCoords] = DefineRegion(image,imagetype)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DefineRegion v .1 Jan/15/2016            
% Bill Green, James Corsetti          
% Description: Grabs a screenshot from the camera, and allows you
% to interactively select a region to subsample.

% image - The image that will be looked at in order to crop a region. This can
% either be the handle for the video feed, or a static image of the phase
% map. Designate this with the imagetype flag.

% imagetype - a flag that gets a screenshot differently depending on
% whether you are pulling from a static image or the video feed. 
% imagetype = 1, use a static screenshot
% imagetype = 2, use the video feed

% cropCoords - vector containing the information about the cropped
% region: [xorigin yorigin width height]

switch imagetype
    case 1
        imagetocrop = image;
        
    case 2
        imagetocrop = getsnapshot(image);
        
end

fig = figure;
[sampleRegion,cropCoords] = imcrop(imagetocrop);

%imshow(sampleRegion);
%crop = 'This is the region you have selected';
cropCoords = floor(cropCoords)+[1,1,0,0]; %The +1 is the offset for indexing the matrix
close(fig)
%j = text(.1*cropCoords(3),.1*cropCoords(4), crop, 'FontSize',14,'Color','r');




        
    