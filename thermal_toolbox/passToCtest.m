%% Test code to pass a phasecube to the .C unwrapping algorithm
clear Z

[X,Y] = meshgrid(0:.005:1,0:.005:1);


figure
% Creates a spatial frequency
for i = 1:10
   Z(:,:,i) = i.*X - (10-i).*Y;
   %Znoise = Z + noise
   Zwrap(:,:,i) = wrapToPi(Z(:,:,i));
   imagesc(Zwrap(:,:,i))
   colorbar
   pause(.3)
end
%% METHOD A: Unwrap each pixel in time

tic
unwrapA = unwrap(Zwrap,[],3);
timeA = toc;


%% METHOD B: Unwrap each sheet in 2D

tic
for i = 1:10
unwrapB(:,:,i) = unwrap2D(Zwrap(:,:,i));
end
timeB = toc;

%% METHOD C: Export to a 2D algorithm written in C++
tic
ZwrapC = Zwrap + pi; %This unwrapper takes input from 0 to 2Pi.
ZwrapC = single(Zwrap);
for i = 1:10
unwrapC(:,:,i) = Miguel_2D_unwrapper(Zwrap(:,:,i));
i
end
timeC = toc;


%% METHOD D: Export to a 3D algorithm written in C++

fileID = fopen('phasecube.bin','w')
fwrite(fileID,Zwrap,'single')
fclose(fileID)

%Open up Developer Command prompt from VStudio2015 and type the following
%command:

% 3DBPPUASL phasecube.bin uphasecube.bin ROW COL FRAME

fileID = fopen('uphasecube.bin')
phasedata = fread(fileID,'single');
phasedata(end+1) = 0;
unwrapD = reshape(phasedata,[201 201 10]);

%% Trying this on the 5_1_16 SNR data

fileID = fopen('phasecube.bin','w')
fwrite(fileID,phasecube(:,:,1:200),'single')
fclose(fileID)

%Open up Developer Command prompt from VStudio2015 and type the following
%command:

% 3DBPPUASL phasecube.bin uphasecube.bin ROW COL FRAME

fileID = fopen('husseinwrap.bin')
phasedata = fread(fileID,'single');
%phasedata(end+1) = 0;
unwrapD = reshape(phasedata,[256 160 200]);

unwrapDzero = bsxfun(@minus,unwrapD,unwrapD(:,:,1));
for i = 1:200
    drift(i) = mean2(unwrapDzero(:,:,i));
end

figure;plot(drift)


%% METHOD E: Load the entire phasecube in, and then make lengthXtime 
%  slices. Stitch the cube back together from these slices.

tic
[dimL,dimB,dimT] = size(phasecube);

unwrapSheet = zeros(dimL,dimB,dimT);

for i = 1:dimB
unwrapSheet(:,i,:) = Miguel_2D_unwrapper(squeeze(phasecube(:,i,:)));
% The averaged time-series drift according to each unwrapped sheet
drift(:,i) = mean(unwrapSheet(:,i,:),2);
i
end


unwrapEzero = bsxfun(@minus,unwrapSheet,unwrapSheet(:,:,1));

timeE = toc;



%% Visually Compare all the methods
figure(3)
for i = 1:10
  subplot(3,1,1)
  imagesc(unwrapC(:,:,i))
  axis equal
  subplot(3,1,2)
  imagesc(unwrapD(:,:,i))
  axis equal
  subplot(3,1,3)
  imagesc(Z(:,:,i))
  axis equal
  pause(.3)
  i
end

%% Call the mex function to unwrap the data
loadlibrary

Zunwrap = calllib('trythis','