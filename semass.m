function [mass, e] = semass(A, Z)
    %all below in MeV
    mp = 938.28; %proton energy equivalent; 1.6726231*10^-27 kg 
    mn = 939.57; %proton energy equivalent; 1.6749286*10^-27 kg  

    av = 15.8; %volume coeff
    as = 18.3; %surface coeff
    ac = 0.714; %coulomb coeff
    aA = 23.2; %asymmetry coeff
    ap = 12.0; %pairing coeff
    
    delta = double(~mod(A,2)).*(double(~mod(Z,2)*2-1));
    
    e = av.*A - as.*A.^(2/3) - ac.*(Z.*(Z-1))./A.^(1/3) - aA.*((A-2.*Z).^2)./A + ap./(A.^(1/2)).*delta; 
    e = e(:,1);
    mass = Z.*mp + (A-Z).*mn - e;
end