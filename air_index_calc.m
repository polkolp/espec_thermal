function n_air = air_index_calc(varargin)
% Usage is air_index_calc(t,p,RH,lam,xCO2)
% t is temperature in Celcius, p is pressure in Pascals, RH is relative
% humidity in percent (0 to 100) lam is wavelength in microns, and xCO2 is
% carbon dioxide concentration (default 450). All input parameters are
% required except for xCO2, which defaults to 450 if not entered.
% To enter 40degC, 120kPa, 75% relative humidity, 633nm, default CO2:
% air_index_calc(40,120000,75,.633,450) = 1.000299418310159

if size(varargin,2) < 4 || size(varargin,2) > 5
    disp('Error: wrong number of input arguments')
    return
end
t = varargin{1};
p = varargin{2};
RH = varargin{3};
lam = varargin{4};
if size(varargin,2) == 4
    xCO2 = 450;
else
    xCO2 = varargin{5};
end

% A-I. Saturation Vapor Pressure
K1 = 1.16705214528e3;   K2 = -7.24213167032e5;      K3 = -1.70738469401e1;
K4 = 1.20208247025e4;   K5 = -3.23255503223e6;      K6 = 1.49151086135e1;
K7 = -4.82326573616e3;  K8 = 4.05113405421e5;       K9 = -2.38555575678e-1;
K10 = 6.50175348448e2;
T = t + 273.15;     Q = T + K9./(T - K10);      A = Q.^2 + K1.*Q + K2;
B = K3.*(Q.^2) + K4.*Q + K5;    C = K6.*(Q.^2) + K7.*Q + K8;
X = -B + sqrt(B.^2 -4.*A.*C);   psv = (10.^6).*((2.*C./X).^4);

% A-II. Determining Humidity
alpha = 1.00062;    beta = 3.14e-8;     gamma = 5.60e-7;
f = alpha + beta.*p + gamma.*t.*t; % enhancement factor f
xv = (RH./100).*f.*psv./p; % xv is mole fraction

% A-III. Ciddor Calculation of Index of Refraction
% part b
w0 = 295.235;   w1 = 2.6422;    w2 = -0.03238;  w3 = 0.004028;
k0 = 238.0185;  k1 = 5792105;   k2 = 57.362;    k3 = 167917;
a0 = 1.58123e-6;    a1 = -2.9331e-8;    a2 = 1.1043e-10;
b0 = 5.707e-6;      b1 = -2.051e-8;
c0 = 1.9898e-4;     c1 = -2.376e-6;
d = 1.83e-11;       e = -0.765e-8;
pR1 = 101325;       TR1 = 288.15;
Za = 0.9995922115;  pvs = 0.00985938;
R = 8.314472;       Mv = 0.018015;
% part c
S = lam.^-2;
% part d
ras = (10^-8).*((k1./(k0 - S)) + (k3/(k2 - S)));
rvs = (1.022e-8).*(w0 + w1.*S + w2.*S.*S + w3.*S.*S.*S);
% part e
Ma = 0.0289635 + (1.2011e-8).*(xCO2 - 400);
raxs = ras.*(1 + (5.34e-7).*(xCO2 - 450));
% part f
T = t + 273.15;
Zm = 1 - (p./T).*((a0 + a1.*t + a2.*t.*t + (b0 + b1.*t).*xv + (c0 + c1.*t).*xv.*xv)) + (d + e.*xv.*xv).*((p./T).^2);
paxs = pR1.*Ma./(Za.*R.*TR1);
pv = xv.*p.*Mv./(Zm.*R.*T);
pa = (1 - xv).*p.*Ma./(Zm.*R.*T);
% part g (index of refraction)
n_air = 1 + (pa./paxs).*raxs + (pv./pvs).*rvs;