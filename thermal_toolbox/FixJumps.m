function [DataOut] = FixJumps(DataIn)

phaseDifferences = zeros(size(DataIn));
phaseDifferences(2:end,:) = diff(DataIn,1,1);

%If the phase difference is larger than Pi, flag the data point for
%overwriting
jumpIndex = find(abs(phaseDifferences)>3.14);
jumpFlags = zeros(size(DataIn));
jumpFlags(jumpIndex) = 1;
jumpFlags = logical(jumpFlags);


% This command will store vectors that record the corresponding Time index
% at which that jump occured. For example, an entry of 6100 indicates that
% the fault occured at T(6100)

%jumptimes = [mod(jumpIndexReflect,70734);mod(jumpIndexBackground,70734)];

%Reassigns the value of the Unwrapped phase based on the flaggedpoints
%DataOut(jumpFlags) = ;