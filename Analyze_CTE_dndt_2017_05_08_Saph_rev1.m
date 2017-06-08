% Code for analyzing TG
clear timetab TAtab Ptab RHtab nairtab TStab
timetab = excel_5_8_2017(:,1);
% TA_Left_tab = excel_5_8_2017(:,2); 
% TA_Right_tab = excel_5_8_2017(:,3);
Ptab = 1000.*98.9259776400000; %pressure
RHtab = 26.5; %relative humidity  
% nair_Left_tab = air_index_calc(TA_AW_tab,Ptab,RHtab,.6328);
% nair_Right_tab = air_index_calc(TA_MM_tab,Ptab,RHtab,.6328);
% % TS_Left_tab = excel_5_8_2017(:,4);
% % TS_Right_tab = excel_5_8_2017(:,5);
% TS_C_tab = excel_01_27_2016_HERO_ramp(:,6);
TAtab = excel_5_8_2017(:,2);%(TA_Left_tab + TA_Right_tab)./2; %TA_AW_tab; %air temperature probe
TStab = (excel_5_8_2017(:,3) + excel_5_8_2017(:,4) + excel_5_8_2017(:,5))./3; % TS_C_tab;% %average of three surface temperature probes
nairtab = air_index_calc(TAtab,Ptab,RHtab,.6328); %calculate index of air
% TS_20mm_tab = excel_10_30_2015_KF_148_8_minus30_run1(:,5);

k = 2*pi/.6328;
%%
close(figure(154));figure(154)
plot(excel_5_8_2017(:,1),TAtab,'m--')
grid on; hold on
plot(excel_5_8_2017(:,1),excel_5_8_2017(:,3),'r')
plot(excel_5_8_2017(:,1),excel_5_8_2017(:,4),'g')
plot(excel_5_8_2017(:,1),excel_5_8_2017(:,5),'b')
plot(excel_5_8_2017(:,1),TStab,'k')

%% starts here
clear D vv 
fileloc = 'I:\2017_05_08_Saph_p20_to32ramp_3hrs\PhaseMaps';
D = dir(fileloc); %lists all files in directory
[vv,idx] = sort([D.datenum]); %sort files by creation time
%%
close(figure(32))
figure(32);
imagesc(phasetemp(:,:));
% imagesc(phasetemp_unfil);
colorbar
colormap summer
axis equal
hold on
MKs = 6;

bg_ww_sft = 200; %314; 1:323;
bg_ww = what_bg_pxls; %[74 82 147];%(305:307);%-19+(305:323);%[1+bg_ww_sft 2+bg_ww_sft 3+bg_ww_sft 20+bg_ww_sft 21+bg_ww_sft 22+bg_ww_sft 39+bg_ww_sft 40+bg_ww_sft 41+bg_ww_sft];%[275 276 277 294 295 296 313 314 315];
% ref_ww_sft = 1; %1:323;
ref_ww = what_ref_pxls; %[137 179 280];%19+(1:19);%[1+ref_ww_sft 2+ref_ww_sft 3+ref_ww_sft 20+ref_ww_sft 21+ref_ww_sft 22+ref_ww_sft 39+ref_ww_sft 40+ref_ww_sft 41+ref_ww_sft];
trans_ww = what_trans_pxls;
regi4_ww = 1;%[1 2 3 4 6 7 8]; % not 5 9
reg55_ww =1;

% Background region (calculate dn/dT)
clear rows_bg cols_bg X_bg Y_bg X_bg_R Y_bg_R bg_cnt bg_pts

% rows_bg=[5:10:75];
% cols_bg=[5:10:75]';
rows_bg=[70:5:150];
cols_bg=[30:5:110]';
[X_bg, Y_bg] = meshgrid(rows_bg,cols_bg);
X_bg_R = reshape(X_bg,length(rows_bg)*length(cols_bg),1);
Y_bg_R = reshape(Y_bg,length(rows_bg)*length(cols_bg),1);

for bg_cnt = 1:length(rows_bg)*length(cols_bg); %combine two backgroud vectors into one matrix
    bg_pts(bg_cnt,1) = X_bg_R(bg_cnt);
    bg_pts(bg_cnt,2) = Y_bg_R(bg_cnt);
end

for bg_cnt2 = 1:size(bg_pts,1)
    plot(bg_pts(bg_cnt2,2),bg_pts(bg_cnt2,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[0 bg_cnt2./size(bg_pts,1) 0])
end

% for bg_ww_cnt = bg_ww;
    plot(bg_pts(bg_ww,2),bg_pts(bg_ww,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% end
bg = [mean(bg_pts(bg_ww,1)) mean(bg_pts(bg_ww,2))];
% bg = [500 190];
text(bg(2)+5,bg(1),'BG')
plot(bg(2),bg(1),'kx')
plot(bg(2),bg(1),'ko')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % % % % % % % % Reflection region (calculate CTE)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear ref rows_ref cols_ref X_ref Y_ref X_ref_R Y_ref_R ref_cnt ref_pts
% ref = [160 510];
% plot(ref(2),ref(1),'kx')

hold on
% rows_ref=[270:10:370];
% cols_ref=[280:10:440]';
% rows_ref=[5:10:75];
% cols_ref=[105:10:175]';
rows_ref=[20:5:100];
cols_ref=[200:5:280]';
[X_ref, Y_ref] = meshgrid(rows_ref,cols_ref);
X_ref_R = reshape(X_ref,length(rows_ref)*length(cols_ref),1);
Y_ref_R = reshape(Y_ref,length(rows_ref)*length(cols_ref),1);

for ref_cnt = 1:length(rows_ref)*length(cols_ref);
    ref_pts(ref_cnt,1) = X_ref_R(ref_cnt);
    ref_pts(ref_cnt,2) = Y_ref_R(ref_cnt);
end


for ref_cnt2 = 1:size(ref_pts,1)
    plot(ref_pts(ref_cnt2,2),ref_pts(ref_cnt2,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[ref_cnt2./size(ref_pts,1) 0 1])
end

plot(ref_pts(ref_ww,2),ref_pts(ref_ww,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
ref = [mean(ref_pts(ref_ww,1)) mean(ref_pts(ref_ww,2))];
% ref = [150 500];
text(ref(2)+5,ref(1),'ref')
plot(ref(2),ref(1),'kx')
plot(ref(2),ref(1),'ko')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Transmission region (calculate dn/dT)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
clear trans rows_trans cols_trans X_trans Y_trans X_trans_R Y_trans_R trans_cnt trans_pts


hold on
rows_trans=[130:5:210];
% cols_trans=[2:5:52]';
cols_trans=[200:5:280]';
[X_trans, Y_trans] = meshgrid(rows_trans,cols_trans);
X_trans_R = reshape(X_trans,length(rows_trans)*length(cols_trans),1);
Y_trans_R = reshape(Y_trans,length(rows_trans)*length(cols_trans),1);

for trans_cnt = 1:length(rows_trans)*length(cols_trans);
    trans_pts(trans_cnt,1) = X_trans_R(trans_cnt);
    trans_pts(trans_cnt,2) = Y_trans_R(trans_cnt);
end


for trans_cnt2 = 1:size(trans_pts,1)
    plot(trans_pts(trans_cnt2,2),trans_pts(trans_cnt2,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1-trans_cnt2./size(trans_pts,1) 1-trans_cnt2./size(trans_pts,1) 1-trans_cnt2./size(trans_pts,1)])
end


plot(trans_pts(trans_ww,2),trans_pts(trans_ww,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
trans = [mean(trans_pts(trans_ww,1)) mean(trans_pts(trans_ww,2))];
% trans = [400 500];
text(trans(2)+5,trans(1),'trans')
plot(trans(2),trans(1),'kx')
plot(trans(2),trans(1),'ko')

xlabel('Pixels')
ylabel('Pixels')

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % region 4 region (CTE of stainless steel)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% clear regi4 rows_regi4 cols_regi4 X_regi4 Y_regi4 X_regi4_R Y_regi4_R regi4_cnt regi4_pts
% % regi4 = [140 190];
% % plot(regi4(2),regi4(1),'kx')
% 
% hold on
% rows_regi4=[20:5:70];
% cols_regi4=[230:5:300]';
% [X_regi4, Y_regi4] = meshgrid(rows_regi4,cols_regi4);
% X_regi4_R = reshape(X_regi4,length(rows_regi4)*length(cols_regi4),1);
% Y_regi4_R = reshape(Y_regi4,length(rows_regi4)*length(cols_regi4),1);
% 
% for regi4_cnt = 1:length(rows_regi4)*length(cols_regi4);
%     regi4_pts(regi4_cnt,1) = X_regi4_R(regi4_cnt);
%     regi4_pts(regi4_cnt,2) = Y_regi4_R(regi4_cnt);
% end
% 
% for regi4_cnt2 = 1:size(regi4_pts,1)
%     plot(regi4_pts(regi4_cnt2,2),regi4_pts(regi4_cnt2,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 regi4_cnt2./size(regi4_pts,1)])
% end
% 
% 
% plot(regi4_pts(regi4_ww,2),regi4_pts(regi4_ww,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% regi4 = [regi4_pts(regi4_ww,1) regi4_pts(regi4_ww,2)];
% % regi4 = [50 190];
% text(regi4(2)+5,regi4(1),'SS')
% plot(regi4(2),regi4(1),'kx')
% plot(regi4(2),regi4(1),'ko')
% 
% xlabel('Pixels')
% ylabel('Pixels')
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% % reg55mission region (calculate dn/dT)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% clear reg55 rows_reg55 cols_reg55 X_reg55 Y_reg55 X_reg55_R Y_reg55_R reg55_cnt reg55_pts
% 
% 
% hold on
% rows_reg55=[92:5:142];
% cols_reg55=[230:5:300]';
% [X_reg55, Y_reg55] = meshgrid(rows_reg55,cols_reg55);
% X_reg55_R = reshape(X_reg55,length(rows_reg55)*length(cols_reg55),1);
% Y_reg55_R = reshape(Y_reg55,length(rows_reg55)*length(cols_reg55),1);
% 
% for reg55_cnt = 1:length(rows_reg55)*length(cols_reg55);
%     reg55_pts(reg55_cnt,1) = X_reg55_R(reg55_cnt);
%     reg55_pts(reg55_cnt,2) = Y_reg55_R(reg55_cnt);
% end
% 
% 
% for reg55_cnt2 = 1:size(reg55_pts,1)
%     plot(reg55_pts(reg55_cnt2,2),reg55_pts(reg55_cnt2,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1-reg55_cnt2./size(reg55_pts,1) 1-reg55_cnt2./size(reg55_pts,1) 1-reg55_cnt2./size(reg55_pts,1)])
% end
% 
% 
% plot(reg55_pts(reg55_ww,2),reg55_pts(reg55_ww,1),'ro','MarkerSize',MKs,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0])
% reg55 = [reg55_pts(reg55_ww,1) reg55_pts(reg55_ww,2)];
% % reg55 = [400 500];
% text(reg55(2)+5,reg55(1),'reg55')
% plot(reg55(2),reg55(1),'kx')
% plot(reg55(2),reg55(1),'ko')
% 
% xlabel('Pixels')
% ylabel('Pixels')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Background for fitting tilt
clear backgroundrows backgroundcols
backgroundrows = 70:150;%597:796;
backgroundcols = 30:110;%180:380;%;978:1339;
line([backgroundcols(1) backgroundcols(length(backgroundcols))],[backgroundrows(1) backgroundrows(1)],'Color',[1 1 1],'LineWidth',2.5)
line([backgroundcols(1) backgroundcols(length(backgroundcols))],[backgroundrows(length(backgroundrows)) backgroundrows(length(backgroundrows))],'Color',[1 1 1],'LineWidth',2.5)
line([backgroundcols(1) backgroundcols(1)],[backgroundrows(1) backgroundrows(length(backgroundrows))],'Color',[1 1 1],'LineWidth',2.5)
line([backgroundcols(length(backgroundcols)) backgroundcols(length(backgroundcols))],[backgroundrows(1) backgroundrows(length(backgroundrows))],'Color',[1 1 1],'LineWidth',2.5)

%%
clear time TAmeas TSmeas RHmeas nair refphase bgphase transphase regi4phase xtilt ytilt
%%
clear bgphase_unfil bgphase_fft refphase_unfil refphase_fft transphase_unfil transphase_fft
%%
f1 =3;%%3;%15201%13502%11501%09501%7501%4992%3201%0003;%%%%3001%0001
f2 =13343;%45003;%65263;%78248;%63;%size(time_PMS,2);%length(time_PMS);%16942%15200%13501%11500%9500%7500%4991%3200;%%%%6100%3000
tic
for frame = f1:1:f2;%f2; %4:10;%26493%length(idx)
    frame
    PM = D(frame).name;
%     load(PM);
    clear phasetemp_unfil phasetemp_fft phasetemp_combo
% %     time(frame-f1+1) = D(frame).datenum;
% % %         time(frame-f1+1) = time_PMS(frame);
        TAmeas(frame-f1+1) = interp1(timetab,TAtab,time(frame-f1+1));%,'cubic','extrap'); %interpolate temperature data (every 5s) to match phasemap acquisition rate(2 per s) 
        TSmeas(frame-f1+1) = interp1(timetab,TStab,time(frame-f1+1));%,'cubic','extrap');
%         TS_20mm_tab(frame-f1+1) = interp1(timetab,TS_20mm_tab,time(frame-f1+1));
        RHmeas(frame-f1+1) = RHtab;
        nair(frame-f1+1) = interp1(timetab,nairtab,time(frame-f1+1));%,'cubic','extrap');
% 
% %         b = puma_ho(phasetemp(backgroundrows,backgroundcols),2);
% %         % Fit 2-D tilt to background data
% %         [x1,y1] = meshgrid(backgroundcols,backgroundrows);
% %         model_function0 = @(vars,r) vars(3) + vars(1)*r(:,1) + vars(2)*r(:,2); % fit background tilt to a plane 
% %         guess = [(b(round(size(b,1)/2),end)-b(round(size(b,1)/2),1))/(backgroundcols(end)-backgroundcols(1)) (b(end,round(size(b,2)/2))-b(1,round(size(b,2)/2)))/(backgroundrows(end)-backgroundrows(1)) mean(mean(b))]; % guess is taking endpoints and finding slope
% %         [varfit1,resid,J,Sigma] = nlinfit([x1(:) y1(:)],b(:),model_function0,guess,statset('TolFun',1e-30,'TolX',1e-20,'MaxIter',1000));
% %         xtilt(frame-f1+1) = varfit1(1);%/(2*nair(frame-f1+1)*k); % putting the tilt into units of waves
% %         ytilt(frame-f1+1) = varfit1(2);%/(2*nair(frame-f1+1)*k);
      
% % % [phasetemp_unfil,phasetemp_fft]=Make_Phasemap_FFT_4(snap,PSramp,0); 
% % % % % % [phasetemp_unfil] = Make_Phasemap_NO_FFT(snap,PSramp,0);

% % % phasetemp_combo = zeros(size(phasetemp_unfil,1)+size(phasetemp_fft,1),size(phasetemp_unfil,2));
% % % phasetemp_combo(1:size(phasetemp_unfil,1),:)=phasetemp_unfil;
% % % phasetemp_combo(165:328,1:176)=phasetemp_fft;

% % % % % % % % % % % % % 
% %         b = puma_ho(phasetemp(backgroundrows,backgroundcols),2);
% %         % Fit 2-D tilt to background data
% %         [x1,y1] = meshgrid(backgroundcols,backgroundrows);
% %         model_function0 = @(vars,r) vars(3) + vars(1)*r(:,1) + vars(2)*r(:,2); % fit background tilt to a plane 
% %         guess = [(b(round(size(b,1)/2),end)-b(round(size(b,1)/2),1))/(backgroundcols(end)-backgroundcols(1)) (b(end,round(size(b,2)/2))-b(1,round(size(b,2)/2)))/(backgroundrows(end)-backgroundrows(1)) mean(mean(b))]; % guess is taking endpoints and finding slope
% %         [varfit1,resid,J,Sigma] = nlinfit([x1(:) y1(:)],b(:),model_function0,guess,statset('TolFun',1e-30,'TolX',1e-20,'MaxIter',1000));
% %         xtilt(frame-f1+1) = varfit1(1);%/(2*nair(frame-f1+1)*k); % putting the tilt into units of waves
% %         ytilt(frame-f1+1) = varfit1(2);%/(2*nair(frame-f1+1)*k);
        
% % %         b = puma_ho(phasetemp_fft(backgroundrows,backgroundcols),2);
% % %         % Fit 2-D tilt to background data
% % %         [x1,y1] = meshgrid(backgroundcols,backgroundrows);
% % %         model_function0 = @(vars,r) vars(3) + vars(1)*r(:,1) + vars(2)*r(:,2); % fit background tilt to a plane 
% % %         guess = [(b(round(size(b,1)/2),end)-b(round(size(b,1)/2),1))/(backgroundcols(end)-backgroundcols(1)) (b(end,round(size(b,2)/2))-b(1,round(size(b,2)/2)))/(backgroundrows(end)-backgroundrows(1)) mean(mean(b))]; % guess is taking endpoints and finding slope
% % %         [varfit1,resid,J,Sigma] = nlinfit([x1(:) y1(:)],b(:),model_function0,guess,statset('TolFun',1e-30,'TolX',1e-20,'MaxIter',1000));
% % %         xtilt3_fft(frame-f1+1) = varfit1(1);%/(2*nair(frame-f1+1)*k); % putting the tilt into units of waves
% % %         ytilt3_fft(frame-f1+1) = varfit1(2);%/(2*nair(frame-f1+1)*k);

% % % % % % % % 


% %         for bgcnt = 1:size(bg_pts,1);
% %             bgphase(frame-f1+1,bgcnt)=phasetemp(bg_pts(bgcnt,1),bg_pts(bgcnt,2));
% % %             bgphase_fft(frame-f1+1,bgcnt)=phasetemp_fft(bg_pts(bgcnt,1),bg_pts(bgcnt,2));
% %         end
% % %         
% %         for refcnt = 1:size(ref_pts,1);
% %             refphase(frame-f1+1,refcnt)=phasetemp(ref_pts(refcnt,1),ref_pts(refcnt,2));
% % %             refphase_fft(frame-f1+1,refcnt)=phasetemp_fft(ref_pts(refcnt,1),ref_pts(refcnt,2));
% %         end
% %         
% %          for transcnt = 1:size(trans_pts,1);
% %             transphase(frame-f1+1,transcnt)=phasetemp(trans_pts(transcnt,1),trans_pts(transcnt,2));
% % %             transphase_fft(frame-f1+1,transcnt)=phasetemp_fft(trans_pts(transcnt,1),trans_pts(transcnt,2));
% %          end
%         
%          for regi4cnt = 1:size(regi4_pts,1);
%             regi4phase(frame-f1+1,regi4cnt)=phasetemp(regi4_pts(regi4cnt,1),regi4_pts(regi4cnt,2));
%          end
%          
%          for reg55cnt = 1:size(reg55_pts,1);
%             reg55phase(frame-f1+1,reg55cnt)=phasetemp(reg55_pts(reg55cnt,1),reg55_pts(reg55cnt,2));
%          end
        
         
%         000001__02_10_2016
%         for transcnt = 1:size(trans_pts,1);
%             transphase(frame-f1+1,transcnt)=phasetemp(trans_pts(transcnt,1),trans_pts(transcnt,2));
%         end
        
%         for regi4cnt = 1:size(regi4_pts,1);
%             regi4phase(frame-f1+1,regi4cnt)=phasetemp(regi4_pts(regi4cnt,1),regi4_pts(regi4cnt,2));
%         end
% % % PM_name = strcat('PM',D(frame).name(1:6),'__02_22_2016');
% % % save(PM_name,'p_combo','-v7.3')
end
toc
%%
xtilt=(xtilt_interp./(2.*nair.*k));
ytilt=(ytilt_interp./(2.*nair.*k));
% toc
disp('done')
%%
xtilt_interp = interp1(1:10:13341,xtilt(1:10:13341),1:1:13341,'spline');
ytilt_interp = interp1(1:10:13341,ytilt(1:10:13341),1:1:13341,'spline');
figure;plot(xtilt_interp,'b.');
hold on;
plot(1:10:13341,xtilt,'ro')
figure;plot(ytilt_interp,'b.');
hold on;
plot(1:10:13341,ytilt,'ro')

%%
% Unwrap background
clear bgphaseuw bguw_cnt
close(figure(34))
figure(34)
% figure
subplot(2,2,1)
bgphaseuw = zeros(size(bgphase,1),size(bgphase,2));
for bguw_cnt = 1:size(bgphase,2);
    bgphaseuw(:,bguw_cnt) = unwrap(bgphase(:,bguw_cnt));
    plot(bgphaseuw(:,bguw_cnt),'Color',[0 1 0]);%%[0 bguw_cnt./size(bg_pts,1) 0])
    hold on
end
grid on
plot(mean(bgphaseuw,2),'k')
xlim([0 150])
xlabel('Count,time')
ylabel('Fringes')
title('Background')

% Unwrap reflection
clear refphaseuw refuw_cnt
subplot(2,2,2)
refphaseuw = zeros(size(refphase,1),size(refphase,2));
for refuw_cnt = 1:size(refphase,2);
    refphaseuw(:,refuw_cnt) = unwrap(refphase(:,refuw_cnt));
    plot(refphaseuw(:,refuw_cnt),'Color',[1 0 1])%[refuw_cnt./size(ref_pts,1) 0 1])
    hold on
end
grid on
xlim([0 150])

xlabel('Count,time')
ylabel('Fringes')
title('Reflection')
% % % % % %%
% Unwrap transmission
clear transphaseuw transuw_cnt
subplot(2,2,3)
transphaseuw = zeros(size(transphase,1),size(transphase,2));
for transuw_cnt = 1:size(transphase,2);
    transphaseuw(:,transuw_cnt) = unwrap(transphase(:,transuw_cnt));
    plot(transphaseuw(:,transuw_cnt),'Color',[0 0 0]);%[1-transuw_cnt./size(trans_pts,1) 1-transuw_cnt./size(trans_pts,1) 1-transuw_cnt./size(trans_pts,1)])
    hold on
end
grid on
xlabel('Count,time')
ylabel('Fringes')
title('Transmission')
% 
% % Unwrap region 4
% clear regi4phaseuw regi4uw_cnt
% subplot(2,2,4)
% regi4phaseuw = zeros(size(regi4phase,1),size(regi4phase,2));
% for regi4uw_cnt = 1:size(regi4phase,2);
%     regi4phaseuw(:,regi4uw_cnt) = unwrap(regi4phase(:,regi4uw_cnt));
%     plot(regi4phaseuw(:,regi4uw_cnt),'Color',[0 0 regi4uw_cnt./size(regi4_pts,1)])
%     hold on
% end
% grid on
% xlabel('Count,time')
% ylabel('Fringes')
% title('Background')
% figure;plot(mean(bgphaseuw_shft(:,62:66),2));hold on;plot(mean(refphaseuw_shft(:,4:9),2));plot(mean(transphaseuw_shft,2));xlim([10240 12556])

% Shifted background
clear bgphaseuw_shft
close(figure(40))
figure(40)
% figure
subplot(2,2,1)
shft = 1;
aaa = shft;
bbb = 14000;%101300;%length(timetime);

clear bgphaseuw_shft
bgphaseuw_shft = zeros(size(bgphaseuw,1),size(bgphaseuw,2));
for bguw_shft_cnt = 1:size(bgphaseuw,2);
    bgphaseuw_shft(:,bguw_shft_cnt) = bgphaseuw(:,bguw_shft_cnt)-bgphaseuw(shft,bguw_shft_cnt);
    plot(bgphaseuw_shft(:,bguw_shft_cnt),'Color',[0 bguw_shft_cnt./size(bgphaseuw,2) 0])
    hold on
end
plot(mean(bgphaseuw_shft,2),'k','LineWidth',2.5)
grid on
% title(['Background, Shift = ' num2str(shft) ' counts, ' num2str(time_PMS(shft)) ' hrs'])
xlabel('Count,time')
ylabel('Fringes')
xlim([aaa bbb])
% ylim([-60 60])

% Shifted reference
clear refphaseuw_shft
subplot(2,2,2)
% shft = 3000;
clear refphaseuw_shft
refphaseuw_shft = zeros(size(refphaseuw,1),size(refphaseuw,2));
for refuw_shft_cnt = 1:size(refphaseuw,2);
    refphaseuw_shft(:,refuw_shft_cnt) = refphaseuw(:,refuw_shft_cnt)-refphaseuw(shft,refuw_shft_cnt);
    plot(refphaseuw_shft(:,refuw_shft_cnt),'Color',[refuw_shft_cnt./size(refphaseuw,2) 0 1])
    hold on
end
grid on
plot(mean(refphaseuw_shft,2),'k','LineWidth',2.5)
% title(['Reflection, Shift = ' num2str(shft) ' counts, ' num2str(time_PMS(shft)) ' hrs'])
xlabel('Count,time')
ylabel('Fringes')
xlim([aaa bbb])
% ylim([-60 60])

% Shifted Transmission
clear transphaseuw_shft
subplot(2,2,3)
% shft = 3000;
clear transphaseuw_shft
transphaseuw_shft = zeros(size(transphaseuw,1),size(transphaseuw,2));
for transuw_shft_cnt = 1:size(transphaseuw,2);
    transphaseuw_shft(:,transuw_shft_cnt) = transphaseuw(:,transuw_shft_cnt)-transphaseuw(shft,transuw_shft_cnt);
    plot(transphaseuw_shft(:,transuw_shft_cnt),'Color',[1-transuw_shft_cnt./size(transphaseuw,2) 1-transuw_shft_cnt./size(transphaseuw,2) 1-transuw_shft_cnt./size(transphaseuw,2)])
    hold on
end
grid on
% title(['Transmission, Shift = ' num2str(shft) ' counts, ' num2str(time_PMS(shft)) ' hrs'])
xlabel('Count,time')
ylabel('Fringes')
xlim([aaa bbb])
% 
% % Shifted regi4mission
% clear regi4phaseuw_shft
% subplot(2,2,4)
% % shft = 3000;
% clear regi4phaseuw_shft
% regi4phaseuw_shft = zeros(size(regi4phaseuw,1),size(regi4phaseuw,2));
% for regi4uw_shft_cnt = 1:size(regi4phaseuw,2);
%     regi4phaseuw_shft(:,regi4uw_shft_cnt) = regi4phaseuw(:,regi4uw_shft_cnt)-regi4phaseuw(shft,regi4uw_shft_cnt);
%     plot(regi4phaseuw_shft(:,regi4uw_shft_cnt),'Color',[0 0 regi4uw_shft_cnt./size(regi4_pts,1)])
%     hold on
% end
% grid on
% title(['regi4mission, Shift = ' num2str(shft) ' counts, ' num2str(time_PMS(shft)) ' hrs'])
% xlabel('Count,time')
% ylabel('Fringes')


