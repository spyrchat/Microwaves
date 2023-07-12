% Define the angle in radians
phi = deg2rad(0:0.001:2*180);
theta = deg2rad(0:0.001:180);


%Lambda value
f = 1e9;
lambda = physconst('LightSpeed') / f;

% Dipole length
d = lambda/2;                       % Dipole spacing
k = 2*pi/lambda;                    % Propagation Constant
I = 1;                              % Current that goes through the Dipole

% Distance from the origin
xi = (-7/2:1:7/2);


%=================================Select Modes for the Patterns================================%


n = 1;                        % 1 for alternating and 0 for non alternating currents    
solid = true;                 % true for 3D radiation pattern
vertical = false;              % true for vertical radiation pattern and flase for horizontal


%==============================================================================================%





%-------Depending on the modes above the "If statements" bellow will output the Desired Graph

if n == 0 && vertical == false && solid == false
    theta = pi/2;
    Etotal = 0;
    for i = 1:1:8
        Ei = (cos(pi/2*cos(theta))./sin(theta)).*exp(1j * k * xi(i) * d .* cos(phi) * sin(theta)); 
        Etotal = Etotal + Ei;
    end
    Etotal = Etotal / max(Etotal(:));
    polarplot(phi,abs(Etotal))
    title('Radiation Pattern')
    subtitle(['d = ',num2str(d/lambda),'λ','  All Currents have the same direction'])

%------------------------------------------------------------------------------------------------------%

elseif n == 1 && vertical == false && solid == false
    theta = pi/2;
    Etotal = 0;
    for i = 1:1:8
        Ei = (cos(pi/2*cos(theta))./sin(theta)).*cos((i-1)*pi*n).*exp(1j * k * xi(i) * d .* cos(phi) * sin(theta));
        Etotal = Etotal + Ei;
    end
    Etotal = Etotal / max(Etotal(:));
    polarplot(phi,abs(Etotal))
    title('Radiation Pattern')
    subtitle(['d = ',num2str(d/lambda),'λ','  Alternating Currents'])

%------------------------------------------------------------------------------------------------------%

elseif n == 0 && vertical == true && solid == false
    theta = deg2rad(0.0001:0.001:360);
    phi = 0;
    Etotal = 0;
    for i = 1:1:8
        Ei = (cos(pi/2*cos(pi/2 - theta))./sin(pi/2 - theta)).*exp(1j .* k .* xi(i) .* d .* cos(phi) .* sin(pi/2 - theta));
        Etotal = Etotal + Ei;
    end
    Etotal = Etotal / max(Etotal(:));
    polarplot(theta,abs(Etotal))
    title('Vertical Radiation Pattern')
    subtitle(['d = ',num2str(d/lambda),'λ','  All Currents have the same direction','   phi =', num2str(phi)])

%-----------------------------------------------------------------------------------------------------------%


elseif n == 1 && vertical == true && solid == false
    theta = deg2rad(0.0001:0.001:360);
    phi = 0;
    Etotal = 0;
    for i = 1:1:8
        Ei = (cos(pi/2*cos(pi/2 - theta))./sin(pi/2 - theta)).*cos((i-1)*pi*n).*exp(1j * k * xi(i) .* d .* cos(phi) .* sin(pi/2 - theta));
        Etotal = Etotal + Ei;
    end
    Etotal = Etotal / max(Etotal(:));
    polarplot(theta,abs(Etotal))
    title('Vertical Radiation Pattern')
    subtitle(['d = ',num2str(d/lambda),'λ','  Alternating Currents','  phi =', num2str(phi)])

%----------------------------------------------------------------------------------------------------------%


elseif n == 0 && solid == true
   theta = linspace(0.00001, pi, 500);
   phi = linspace(0.00001, 2*pi, 500);
   [Theta, Phi] = meshgrid(phi, theta);
   
   Etotal = zeros(length(theta), length(phi));
   for i = 1:length(xi) 
        Ei = (exp(1j * k * xi(i) * d .* cos(Phi) ).* sin(Theta)).*cos(pi/2 .*cos(Theta))./sin(Theta);
        Etotal = Ei + Etotal;
   end
    
   Dir = Directivity(Etotal,theta,phi);
   fprintf('The Directivity Calculated using the Riemann Approximation is: %f\n', Dir)

   % E=@(phi,theta) 1/(max(max(abs(Etotal)))) * (cos(pi/2*cos(theta))./sin(theta)).*abs(2*cos(7/2 .* k .*d.*cos(phi).*sin(theta)) + 2*cos(5/2 .*k *d.*cos(phi).*sin(theta)) + 2*cos(3/2 .*k .*d.*cos(phi).*sin(theta)) + 2*cos(1/2 .*k *d.*cos(phi).*sin(theta)));
   % x=@(phi,theta) E(phi,theta).*cos(phi).*sin(theta);
   % y=@(phi,theta) E(phi,theta).*sin(phi).*sin(theta);
   % z=@(phi,theta) E(phi,theta).*cos(theta);
   % fsurf(x,y,z,[0 2*pi 0 pi]);
   % colorbar;
   % title(['3D Radiation pattern for d = ', num2str(d/lambda),'λ']);
   % subtitle(['All Currents have the same direction']);
   % xlabel('X');
   % ylabel('Y');
   % zlabel('Z');


%----------------------------------------------------------------------------------------------------------%

elseif n == 1 && solid == true 
   
   theta = linspace(0.01, pi, 600);
   phi = linspace(0.01, 2*pi, 600);
   [Theta, Phi] = meshgrid(theta, phi);
   Etotal = zeros(length(theta), length(phi));
   
   ctr = 0;
   Etotal = (cos(pi/2*cos(Theta))./sin(Theta)).*abs(-sin(7*pi/2.*cos(Phi).*sin(Theta))+sin(5*pi/2.*cos(Phi).*sin(Theta))-sin(3*pi/2.*cos(Phi).*sin(Theta))+sin(pi/2.*cos(Phi).*sin(Theta)));
   
   %Etotal = (cos(pi/2*cos(Theta))./sin(Theta)).*abs(cos(7*pi/2.*cos(Phi).*sin(Theta))+cos(5*pi/2.*cos(Phi).*sin(Theta))+cos(3*pi/2.*cos(Phi).*sin(Theta))+cos(pi/2.*cos(Phi).*sin(Theta))); 
   Pr = (abs(Etotal))^2 /(2*120*pi);
    Prmax = max(max(Pr));
    Ptotal = 0;
    t = deg2rad(1);
    for phiIndex = 1:(length(phi)-1)
        for thetaIndex = 1:(length(theta)-1)
            Ptotal = (Ptotal + Pr(phiIndex, thetaIndex).*sin(Theta(phiIndex,thetaIndex)).*t^2);
        end
    end
    
    D = 4*pi*Prmax/Ptotal;
   
   fprintf('The Directivity Calculated using the Riemann Approximation is: %f\n', D)
    
   % % E=@(phi,theta) 1/(max(max(abs(Etotal)))) * (cos(pi/2*cos(theta))./sin(theta)).*abs(2j*sin(7/2 .* k .*d.*cos(phi).*sin(theta)) - 2j*sin(5/2 *k *d.*cos(phi).*sin(theta)) + 2j*sin(3/2 .*k .*d.*cos(phi).*sin(theta)) - 2j*sin(1/2 .*k *d.*cos(phi).*sin(theta)));
   % % x=@(phi,theta) E(phi,theta).*cos(phi).*sin(theta);
   % % y=@(phi,theta) E(phi,theta).*sin(phi).*sin(theta);
   % % z=@(phi,theta) E(phi,theta).*cos(theta);
   % % fsurf(x,y,z,[0 2*pi 0 pi]);
   % % colorbar;
   % % title(['3D Radiation pattern for d = ', num2str(d/lambda),'λ']);
   % % subtitle(['Alternating Currents']);
   % % xlabel('X');
   % % ylabel('Y');
   % % zlabel('Z');
end
