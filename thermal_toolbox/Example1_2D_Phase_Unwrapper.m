clear; clc; close all
%This Matlab program calls a 2D Phase unwrapper written in C language
%The wrapped phase image is floating point data type. 
%Also, the unwrapped phase image is floating data type.

%read the wrapped phase image from a file. 
image_width = 512;
image_height = 512;
fid = fopen('wrapped phase map float 512X512.dat');
WrappedPhase = fread(fid, image_width * image_height, 'float');
fclose(fid);
WrappedPhase = reshape(WrappedPhase, image_width, image_height);
figure(1)
colormap(gray(256))
imagesc(WrappedPhase);

tic
%call the 2D phase unwrapper from C language
%To compile the C code: in Matlab Command Window type
         mex Miguel_2D_unwrapper.cpp
%The wrapped phase should have the single (float in C) data type
WrappedPhase = single(WrappedPhase);
UnwrappedPhase = Miguel_2D_unwrapper(WrappedPhase);
figure(2)
colormap(gray(256))
imagesc(UnwrappedPhase);
toc
