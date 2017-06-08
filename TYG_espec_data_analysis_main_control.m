%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main control of espec data analysis code created by Tianyi Yang on 051517. Modified on the base of
% Analyze_CTE_dndt_.m by James Corsetti.
% 
% 
% 
%
%=================================================
%% defining path and other basic parameters
clear; clc; close all;

path(path,'C:\Users\Chopin\OneDrive\espec_thermal\thermal_toolbox') %add analysis functions in path

path = 'E:\Documents\Research\Ellis\MGRIN\Data\';
fileloc = 'F:\02_09_2016_testFourier\PhaseMaps'; %phasemap folder directory
D = dir(fileloc); %lists all files in directory
[vv,idx] = sort([D.datenum]); %sort files by creation time, idx is index after sorting
%%
figure;
load([fileloc '\' D(3).name]);
phasetemp = im2double(snap);
imshow(snap,[])

colorbar
colormap summer
axis equal
hold on
MKs = 6;

% %%
% % Background region (calculate dn/dT)
% areasize = 80; %in pixels 
% areagap = 5;
% 
% rows_bg=[70:5:150]; %80*80 area
% cols_bg=[30:5:110]';
% [X_bg, Y_bg] = meshgrid(rows_bg,cols_bg);
% X_bg_R = reshape(X_bg,length(rows_bg)*length(cols_bg),1);
% Y_bg_R = reshape(Y_bg,length(rows_bg)*length(cols_bg),1);
% 
% for bg_cnt = 1:length(rows_bg)*length(cols_bg); %combine two backgroud vectors into one matrix
%     bg_pts(bg_cnt,1) = X_bg_R(bg_cnt);
%     bg_pts(bg_cnt,2) = Y_bg_R(bg_cnt);
% end
% 
% for bg_cnt2 = 1:size(bg_pts,1)
%     plot(bg_pts(bg_cnt2,2),bg_pts(bg_cnt2,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[0 bg_cnt2./size(bg_pts,1) 0])
% end
% 
% % for bg_ww_cnt = bg_ww;
%     plot(bg_pts(bg_ww,2),bg_pts(bg_ww,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% % end
% bg = [mean(bg_pts(bg_ww,1)) mean(bg_pts(bg_ww,2))];
% % bg = [500 190];
% text(bg(2)+5,bg(1),'BG')
% plot(bg(2),bg(1),'kx')
% plot(bg(2),bg(1),'ko')