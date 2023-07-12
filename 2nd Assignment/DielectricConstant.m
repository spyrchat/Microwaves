c = 3e8;
m = 1;
n = 0;
f = 9e9;
a = 22.86 * 1e-3;
b = 10.16* 1e-3;
d1 = 1.5e-3;
d2 = 1.517e-3;

%Cutoff Frequency
fc=c/(2*a);


%Measurements Without Dielectric
WD = [48.4, 73.0, 97.5, 122.2] * 1e-3;

lambdaG = 2 * ((WD(2)-WD(1)) + (WD(3)-WD(2)) + (WD(4)-WD(3))) * 1/(length(WD)-1);

%=========================Material A======================%
SC = [46.9, 71.5, 96.0, 120.7] * 1e-3;
OC = [55.6, 80.2, 104.7, 129.4] * 1e-3;


dAmin = 1/3 * (SC(2)-WD(1) + SC(3)-WD(2) + SC(4)-WD(3));


dBmin = 1/4 * (OC(1)-WD(1) + OC(2)-WD(2) + OC(3)-WD(3) + OC(4)-WD(4));

er =  (fc/f)^2 - (1-(fc/f)^2)/(tan(2*pi/(lambdaG)*dAmin) * tan(2*pi/(lambdaG)*dBmin));

%=====================Material B=========================%
SC2 = [46.7, 71.3, 95.8, 120.5] * 1e-3;
OC2 = [52.6, 77.0, 101.7, 126.3] * 1e-3;

dAmin2 = 1/3 * (SC2(2)-WD(1) + SC2(3)-WD(2) + SC2(4)-WD(3));
dBmin2 = 1/4 * (OC2(1)-WD(1) + OC2(2)-WD(2) + OC2(3)-WD(3) + OC2(4)-WD(4));


er2 =  (fc/f)^2 - (1-(fc/f)^2)/(tan(2*pi/(lambdaG)*dAmin2) * tan(2*pi/(lambdaG)*dBmin2));


%=====================Only using the first measurment set==============%


%-------------------------------Material 1-----------------------------%
x= fzero(@(x) (tan(x)/x + c/(2*pi*d1*f)*tan(2*pi/lambdaG*dAmin)/sqrt(1-(fc/f)^2) ), 2);
er_est = (fc/f)^2+(c*x/(2*pi*d1*f))^2;

%-------------------------------Material 2-----------------------------%
y = fzero(@(y) (tan(y)/y + c/(2*pi*d2*f)*tan(2*pi/lambdaG*dAmin2)/sqrt(1-(fc/f)^2) ), 2);
er_est2 = (fc/f)^2+(c*y/(2*pi*d2*f))^2;

%======================================================================%
