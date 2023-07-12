function D = Directivity(Etotal,theta,phi)

    Pr = (abs(Etotal))^2 /(2*120*pi);
    Prmax = max(max(Pr));
    Ptotal = 0;
    
    [Theta, Phi] = meshgrid(theta, phi); 
    
    for phiIndex = 1:(length(phi)-1)
        for thetaIndex = 1:(length(theta)-1)
            dphi = phi(phiIndex+1) - phi(phiIndex);
            dtheta = theta(thetaIndex+1) - theta(thetaIndex);
            Ptotal = (Ptotal + Pr(phiIndex,thetaIndex).*sin(Theta(thetaIndex,phiIndex)).*dphi.*dtheta);
        end
    end
    
    D = 4*pi*Prmax/Ptotal;
end


