function UnwrapPhase

close(figure(3))
figure(3)
subplot(2,4,[1:2 5:6])
surf(phasetemp)
% imagesc(phasetemp)
% ylim([0 550])
view(2)
shading interp
colorbar
hold on
vcut = 220;
hcut =80;
line([vcut vcut],[0 size(phasetemp,1)],[4 4],'LineWidth',1,'Color','k')
line([0 size(phasetemp,2)],[hcut hcut],[4 4],'LineWidth',1,'Color','k')
axis equal

subplot(2,4,3)
plot((phasetemp(:,vcut)),'k')
title('Vertical')
grid on
xlim([0 size(phasetemp,1)])
subplot(2,4,7)
plot((phasetemp(hcut,:)),'k')
title('Horizontal')
grid on
xlim([0 size(phasetemp,2)])


subplot(2,4,4)
v_unwrap = unwrap(phasetemp(:,vcut)');
vymin = 1; vymax = size(phasetemp,1); vxdata = vymin:vymax;
vydata = v_unwrap(vymin:vymax)./(2.*pi);
myfit = fittype('a*x+b', 'independent', 'x', 'dependent', 'y');
myopt = fitoptions('Method', 'NonlinearLeastSquares');
myopt.StartPoint = [-1 20]; %
[f, g] = fit(vxdata', vydata', myfit, myopt);
plot(vxdata,vydata,'r.');hold on;plot(vxdata,f.a.*vxdata+f.b,'LineWidth',2)
xlabel('Pixel')
ylabel('Fringes')
grid on
title('Vertical')
xlim([0 512])

subplot(2,4,8)
h_unwrap = unwrap(phasetemp(hcut,:));
hymin = 1; hymax = size(phasetemp,2); hxdata = hymin:hymax;
hydata = h_unwrap(hymin:hymax)./(2.*pi);
myfit = fittype('a*x+b', 'independent', 'x', 'dependent', 'y');
myopt = fitoptions('Method', 'NonlinearLeastSquares');
myopt.StartPoint = [-1 20]; %
[f, g] = fit(hxdata', hydata', myfit, myopt);
plot(hxdata,hydata,'r.');hold on;plot(hxdata,f.a.*hxdata+f.b,'LineWidth',2)
xlabel('Pixel')
ylabel('Fringes')
grid on
title('Horizontal')
xlim([0 640])