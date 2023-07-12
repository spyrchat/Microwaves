clear all;
N =201; fmin = 0; fmax = 4e9; f0 = 1e9; 
f = fmin:((fmax-fmin)/N):fmax;

L0 = 1/5;
L1 = 1/10;
L2 = 0.13;
beta_L0 = 2 * pi * L0 * f / f0; 
beta_L1 = 2 * pi * L1 * f / f0; 
beta_L2 = 2 * pi * L2 * f / f0; 

ZL = 100 + 1./(1j*2*pi*f*2e-12); Z0 = 50; Z2 = 1./(1j*2*pi*f*2.7e-12);

ZinA = Z0 *(ZL+1j*Z0.*tan(beta_L0))./(Z0 + 1j*ZL.*tan(beta_L0));

ZinSC = 1j *Z0.*tan(beta_L2);

ZL1 = (ZinA.*ZinSC)./(ZinA+ZinSC);

ZinC = Z0 .* (ZL1 + 1j*Z0.*tan(beta_L1)) ./ (Z0 + 1j*ZL1.*tan(beta_L1));


ZL2 = (ZinC.*Z2)./(ZinC + Z2);

Zin = ZL2;

S11 = (Zin-Z0)./(Zin+Z0);

%plot(f/1e9,abs(S11)); 
%xlabel('Frequency (GHz)'); ylabel('Reflection coefficient ');
plot(f/1e9,20*log10(abs(S11))); 
xlabel('Frequency (GHz)'); ylabel('Reflection coefficient (dB)');

