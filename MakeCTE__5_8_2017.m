clear startframe endframe
%
bgphaseuw1 = mean(bgphaseuw_shft(:,bg_ww),2);%bguw_shft_cnt);
refphaseuw1 = mean(refphaseuw_shft(:,ref_ww),2);
% refphaseuw1 = mean(regi4phaseuw_shft(:,regi4_ww),2);
transphaseuw1 = mean(transphaseuw_shft(:,trans_ww),2);
% figure;plot(bgphaseuw1);hold on;plot(refphaseuw1);plot(transphaseuw1)

dOPD1 = (bgphaseuw1)./k;
dOPD2 = (refphaseuw1)./k;
dOPD3 = (transphaseuw1)./k;

% dOPD1 = dOPD1(1:78241);
% dOPD2 = dOPD2(1:78241);
% dOPD3 = dOPD3(1:78241);


close(figure(54));figure(54);plot(dOPD1,'g');hold on;plot(dOPD2,'m');plot(dOPD3,'k');grid on

close(figure(55));figure(55);plot(TSmeas,dOPD1,'g');hold on;plot(TSmeas,dOPD2,'m');grid on;plot(TSmeas,dOPD3,'k');grid on
xlabel('Temp')
ylabel('Change in OPD ({\mu}m)');
legend('Background','Sample','Transmission')


startframe=shft;%shft;
endframe = length(TSmeas);%35500;%16990;
dOPD1_bnd = dOPD1(startframe:endframe)';
dOPD2_bnd = dOPD2(startframe:endframe)';
dOPD3_bnd = dOPD3(startframe:endframe)';
% ytilt(4797:length(ytilt))=ytilt(4797:length(ytilt))-(ytilt(4797)-ytilt(4782));

nair_bnd = nair(startframe:endframe);
xtilt_bnd = xtilt(startframe:endframe);
ytilt_bnd = ytilt(startframe:endframe);
TAmeas_bnd = TAmeas(startframe:endframe);
TSmeas_bnd = TSmeas(startframe:endframe);
time_bnd = time(startframe:endframe);

% xtilt_bnd = zeros(1,length(nair_bnd));
% ytilt_bnd = zeros(1,length(nair_bnd));

% t0=3479.8;%3000;
t0=5144;
n0=1.4329;


% dOPD2_bnd=dOPD2_bnd.*.987;
% dOPD3_bnd=dOPD3_bnd.*1;


tsamp = (nair_bnd(1).*t0 + 0.5.*(dOPD1_bnd-dOPD2_bnd) - 1.*((nair_bnd(1).*xtilt_bnd(1) - nair_bnd.*xtilt_bnd).*(ref(2)-bg(2)) + (nair_bnd(1).*ytilt_bnd(1) - nair_bnd.*ytilt_bnd).*(ref(1)-bg(1))))./nair_bnd;
nsamp = (n0.*t0 + .5.*(dOPD3_bnd-dOPD2_bnd) + 1.*((nair_bnd.*xtilt_bnd - nair_bnd(1).*xtilt_bnd(1)).*(ref(2) - trans(2)) + (nair_bnd.*ytilt_bnd - nair_bnd(1).*ytilt_bnd(1)).*(ref(1) - trans(1))))./tsamp;

dt = tsamp - tsamp(1);
% % dtA=2584;dtB=2578;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=2865;dtB=2863;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=2944;dtB=2942;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=7209;dtB=7206;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=8992;dtB=8984;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=11292;dtB=11273;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=42507;dtB=42505;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=42721;dtB=42719;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=42767;dtB=42762;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=42786;dtB=42783;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=43015;dtB=42969;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=44221;dtB=44214;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));
% % dtA=51192;dtB=51188;dt(dtA:length(dt))=dt(dtA:length(dt))-(dt(dtA)-dt(dtB));


dn = nsamp - nsamp(1);
% % % dnA=2584;dnB=2578;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % % dnA=2865;dnB=2863;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % % dnA=2944;dnB=2942;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=5747;dnB=5742;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=7209;dnB=7201;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=8992;dnB=8984;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=11287;dnB=11276;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=42510;dnB=42503;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=42721;dnB=42718;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=42767;dnB=42763;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=42786;dnB=42783;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=42974;dnB=42965;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=43004;dnB=43000;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=43015;dnB=43014;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=44221;dnB=44218;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));
% % dnA=51191;dnB=51187;dn(dnA:length(dn))=dn(dnA:length(dn))-(dn(dnA)-dn(dnB));


close(figure(78))
figure(78)

% testCTE_calc([1 19423 36333 53900],TSmeas_bnd,dt,t0)
% testdndT_calc([1 19423 36333 53900],TSmeas_bnd,dn,n0,10)

subplot(1,2,1)
start2 = 2000;%
end2 = 14000;
xtilt_bnd=xtilt_bnd(start2:end2);
ytilt_bnd=ytilt_bnd(start2:end2);
dOPD1_bnd=dOPD1_bnd(start2:end2);
dOPD2_bnd=dOPD2_bnd(start2:end2);
dOPD3_bnd=dOPD3_bnd(start2:end2);
TSmeas_bnd=TSmeas_bnd(start2:end2);
TAmeas_bnd=TAmeas_bnd(start2:end2);
dt=dt(start2:end2);
dn=dn(start2:end2);
nair_bnd=nair_bnd(start2:end2);

disp(['Temp range --> from ' num2str(TSmeas_bnd(1)) '',char(176) 'C to ' num2str(TSmeas_bnd(length(TSmeas_bnd))) '',char(176) ' C'])

cnfd_val=.999;
dt_FIT_CTE = nlinfit((TSmeas_bnd),dt,@(b1,TSmeas_bnd)(b1(1).*TSmeas_bnd+b1(2)),[.06 -1.7]);
% hold on
CTE_slp = dt_FIT_CTE(1)./(t0);
plot(TSmeas_bnd,dt,'bo','MarkerSize',2);xlabel(['Temp (',char(176),'C)'],'FontSize',12);ylabel('Change in thickness ({\mu}m)','FontSize',12);title(['Sample: JC022-7 (t = ' sprintf('%.3f',t0./1000) 'mm), CTE = ' sprintf('%.2f',CTE_slp*(1e6))  'x10^{-6} [1/',char(176) 'C]'],'FontSize',12)
hold on; %plot(TAmeas(eva_end),dt(eva_end),'ro');
xpy = linspace(min(TSmeas_bnd),max(TSmeas_bnd),201); plot(xpy,dt_FIT_CTE(2)+dt_FIT_CTE(1).*xpy,'r','LineWidth',2)
pub_CTE=t0*(7.5e-6);
fit_diff = mean(dt_FIT_CTE(2)+dt_FIT_CTE(1).*xpy)-mean(dt_FIT_CTE(2)+pub_CTE.*xpy);
% % plot(xpy,fit_diff+dt_FIT_CTE(2)+pub_CTE.*xpy,'g','LineWidth',2,'Color',[0 .75 0]) % SS .375"
pub_CTE_up = t0*(8e-6);
pub_CTE_low = t0*(7e-6);
fit_diff1 = mean(dt_FIT_CTE(2)+dt_FIT_CTE(1).*xpy)-mean(dt_FIT_CTE(2)+pub_CTE_up.*xpy);
% % plot(xpy,fit_diff1+dt_FIT_CTE(2)+pub_CTE_up.*xpy,'g--','LineWidth',2,'Color',[0 .75 0]) % SS .375"
fit_diff1 = mean(dt_FIT_CTE(2)+dt_FIT_CTE(1).*xpy)-mean(dt_FIT_CTE(2)+pub_CTE_low.*xpy);
% % plot(xpy,fit_diff1+dt_FIT_CTE(2)+pub_CTE_low.*xpy,'g--','LineWidth',2,'Color',[0 .75 0]) % SS .375"
fitresult = fit(TSmeas_bnd',dt','poly1');
ci = confint(fitresult,cnfd_val);
disp(['dL/dT CTE of JC022-7 = ' sprintf('%.3f',CTE_slp*(1e6)) 'x10^{-6} [1/',char(176) 'C]'])%  ''])
% disp(['Error of CTE of JC022-7 = ' sprintf('%.2f',(CTE_slp-ci(1,1)./t0)*(1e6)) 'x10^{-6} [1/',char(176) 'C]'])
legend('Measured','Linear Fit');%,'Published')
% xlim([-40 20])
grid on
% ylim([-1 3])
CTE_delta=(1./t0).*((dt(end)-dt(1))./(TSmeas_bnd(end)-TSmeas_bnd(1)));
disp(['DL/DT CTE of JC022-7 = ' sprintf('%.3f',CTE_delta*(1e6)) 'x10^{-6} [1/',char(176) 'C]'])%  ''])

nsamp=n0;
% close(figure(100))
% figure(100)
subplot(1,2,2)
nsamp_FIT_dndT = nlinfit((TSmeas_bnd),dn,@(b2,TSmeas_bnd)(b2(1).*TSmeas_bnd+b2(2)),[.06 -1.7]);
dndT_slp = nsamp_FIT_dndT(1);
plot(TSmeas_bnd,dn,'bo','MarkerSize',2);xlabel(['Temp (',char(176),'C)'],'FontSize',12);ylabel('Change in index of refraction','FontSize',12);title(['Sample: JC022-7 (n =' sprintf('%.3f',n0) '), dn/dT = ' sprintf('%.2f',dndT_slp*(1e6))  'x10^{-6} [1/',char(176) 'C]'],'FontSize',12)
hold on; %plot(TAmeas(eva_end),dt(eva_end),'ro');
xpy = linspace(min(TSmeas_bnd),max(TSmeas_bnd),201); plot(xpy,nsamp_FIT_dndT(2)+nsamp_FIT_dndT(1).*xpy,'r','LineWidth',2)
pub_dndT=14.3e-6;
fit_diff = mean(nsamp_FIT_dndT(2)+nsamp_FIT_dndT(1).*xpy)-mean(nsamp_FIT_dndT(2)+pub_dndT.*xpy);
% % plot(xpy,fit_diff+nsamp_FIT_dndT(2)+pub_dndT.*xpy,'g','LineWidth',2,'Color',[0 .75 0]) % SS .375"
pub_dndT_up = 14.8e-6;
pub_dndT_low = 13.8e-6;
fit_diff1 = mean(nsamp_FIT_dndT(2)+nsamp_FIT_dndT(1).*xpy)-mean(nsamp_FIT_dndT(2)+pub_dndT_up.*xpy);
% % plot(xpy,fit_diff1+nsamp_FIT_dndT(2)+pub_dndT_up.*xpy,'g--','LineWidth',2,'Color',[0 .75 0]) % SS .375"
fit_diff1 = mean(nsamp_FIT_dndT(2)+nsamp_FIT_dndT(1).*xpy)-mean(nsamp_FIT_dndT(2)+pub_dndT_low.*xpy);
% % plot(xpy,fit_diff1+nsamp_FIT_dndT(2)+pub_dndT_low.*xpy,'g--','LineWidth',2,'Color',[0 .75 0]) % SS .375"
fitresult = fit(TSmeas_bnd',dn','poly1');
ci = confint(fitresult,cnfd_val);
% % % % plot(xpy, ci(1,2)+ci(1,1).*xpy,'r--','LineWidth',2)
% % % % plot(xpy, ci(2,2)+ci(2,1).*xpy,'r--','LineWidth',2)
disp(['dn/dT of JC022-7 = ' sprintf('%.2f',dndT_slp*(1e6)) 'x10^{-6} [1/',char(176) 'C]'])%  ''])
disp(['Error of dn/dT of JC022-7 = ' sprintf('%.2f',(dndT_slp-ci(1,1))*(1e6)) 'x10^{-6} [1/',char(176) 'C]'])
legend('Measured','Linear Fit');%,'Published')
% xlim([19 33])
grid on

% close all;figure;plot((1:length(dt))-20000,dt);xlim([0 10000])

