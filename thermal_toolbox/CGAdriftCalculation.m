%% The predicted optical path change in the Background

cteFG = 8.3e-6;
cteFS = 0.55e-6;
cteBK7= 7.1e-6;
cteSS = 16.9e-6;
cteUltem=5.5e-5; %10 times HIGHER than I thought this whole time
cteInvar=1.3e-6;

dT = 10; %degC

lFG = 4e-3;
lFS = 12.7e-3;
lSS = 3e-3; % The thickness of the table underneath the chamber
lUltem = 32e-3;
lair = 35e-3

dndt = -.1e-6; %air

CavityChange = cteUltem*lUltem*dT + cteFS*lFS*dT + lSS*cteSS*dT + 3* cteFG*lFG*dT + dndt*1.00029*lair*dT;
OPD1 = 2*CavityChange;

% Answer is 18 micrometers. I saw 60 fringes -> 20 micrometers counting fringes
% yesterday. There's our answer. For my 10 degree test, I saw 45
% micrometers of OPD on the background. Now this number looks reasonable.

CavityChangeFix = cteInvar*lUltem*dT + cteFS*lFS*dT + lSS*cteSS*dT + 3* cteFG*lFG*dT;
OPD1Fix = 2*CavityChangeFix;


