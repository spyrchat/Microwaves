clear all;
Z0 = 50; Z1 = 101.6; Z2 = 101.6; Z3 = 50;
Zs0 = 98.45; Zs1 = 43.6; Zs2 = 98.45;

N = 201; fmin = 0; fmax = 3e9; f0 = 1e9; 
f = fmin:((fmax-fmin)/N):fmax;
L = 1/8; % It means that L = lambda/8 for f=f0;
beta_L = 2*pi*L*f/f0; 

Zin_s0 = -1j*Zs0./tan(beta_L); 
ZL0 = (Zin_s0*Z0)./(Zin_s0 + Z0);
Zin1 = Z1*(ZL0 + 1j*Z1*tan(beta_L))./(Z1+ 1j*ZL0.*tan(beta_L));

Zin_s1 = -1j*Zs1./tan(beta_L); 
ZL1 = (Zin_s1.*Zin1)./(Zin_s1 + Zin1);
Zin2 = Z2*(ZL1 + 1j*Z2*tan(beta_L))./(Z2 + 1j*ZL1.*tan(beta_L));

Zin_s2 = -1j*Zs2./tan(beta_L); 
ZL2 = (Zin_s2.*Zin2)./(Zin_s2 + Zin2);
Zin = Z3*(ZL2 + 1j*Z3*tan(beta_L))./(Z3 + 1j*ZL2.*tan(beta_L));

S11 = (Zin-Z0)./(Zin+Z0);
S11_dB = 20*log10(abs(S11));

SWR = (1+ abs(S11))./(1-abs(S11));
SWR(SWR>=10)=10;
%plot(f/1e9,SWR); 
%xlabel('Frequency (GHz)'); ylabel('SWR');

plot(f/1e9,20*log10(abs(S11))); 
xlabel('Frequency (GHz)'); ylabel('Reflection coefficient (dB)');
