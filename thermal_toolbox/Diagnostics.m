%% Diagnostic Tools
load('E:\Documents\Research\Ellis\MGRIN\Data\5_01_16_SNR\9_19_16_DriftData.mat')

%This code films the evolution of the data with time.
i = 1
k = 1
figure
for i = 1:2:size(unwrapD,3)
    imagesc(unwrapDzero(:,:,i))
    colorbar
    title('5/1/16 Drift Test with Outliers')
    xlabel('Unwrapped Phase [rad]')
    drift(k) = mean2(unwrapD(:,:,i));
    G(k) = getframe;
    k = k+1
end

unwrapDzero = bsxfun(@minus,unwrapD,unwrapD(:,:,1));

for i = 1:size(unwrapDzero,3)
    drift(i) = mean2(unwrapDzero(:,:,i));
end

figure;plot(drift)

figure;hold on;plot(DriftOut);plot(DriftNoOut);plot(DriftNoOutTight)
ylabel('Mean Drift [rad]')
legend('With Outliers','Without Outliers','pi/3 Threshold')

