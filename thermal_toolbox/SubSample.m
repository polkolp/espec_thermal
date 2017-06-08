function bg_pts = SubSample(coords,meshsize)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subsample v .1 Jan/4/2016            
% Bill Green, James Corsetti          
% Description: Used to create a grid of points that will be taken
% off the phasemaps. OPD will be calculated at each of these points.

% coords - a vector that describes the desired cropped region (Get this
% from DefineRegion)

% meshsize - A meshsize of 10 indicates that every 10th pixel will be used

% bg_pts - An nx2 matrix where n is a point at which the phase will
% be taken from the phasemaps in order to calculate the OPD.

rows_bg=[coords(1):meshsize:coords(1)+coords(3)];
cols_bg=[coords(2):meshsize:coords(2)+coords(4)];
[X_bg, Y_bg] = meshgrid(rows_bg,cols_bg);
X_bg_R = reshape(X_bg,length(rows_bg)*length(cols_bg),1);
Y_bg_R = reshape(Y_bg,length(rows_bg)*length(cols_bg),1);

for bg_cnt = 1:length(rows_bg)*length(cols_bg);
    bg_pts(bg_cnt,1) = X_bg_R(bg_cnt);
    bg_pts(bg_cnt,2) = Y_bg_R(bg_cnt);
end