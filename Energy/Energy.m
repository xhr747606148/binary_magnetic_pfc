
function r = Energy(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,pp,qq,hx,hy,k2,w0,a1,a2)
    phi1 = real(phi1); 
    fx1 = real(ifft2(pp.*fft2(phi1))); fy1 =  real(ifft2(qq.*fft2(phi1)));
    phi2 = real(phi2); 
    fx2 = real(ifft2(pp.*fft2(phi2))); fy2 =  real(ifft2(qq.*fft2(phi2)));
    norm1 = M1.^2 + M2.^2;
    fM1 = real(ifft2(pp.*fft2(M1))); fM2 =  real(ifft2(qq.*fft2(M1)));
    norm2_1 = (M1 .* fx1 + M2 .* fy1).^2;
    norm2_2 = (M1 .* fx2 + M2 .* fy2).^2;  
    normH = M1.*H1 + M2.*H2;
    N =  -normH + phi1.^4 / 4  - epsilon .* phi1.^2 / 2 + phi2.^4 / 4  - epsilon .* phi2.^2 / 2 + eta .* (abs(phi1).^3 - phi1.^3) / 3 + eta .* (abs(phi2).^3 - phi2.^3) / 3 + (gamma12 .* phi1.^2 .* phi2.^2)./2  + beta .* norm1.^2 / 4  - alpha .* norm1 / 2 - gamma1 .* phi1 .* norm1 - gamma2 .* phi2 .* norm1 + theta .* norm1.^3 / 6 + theta1 .* (fx1.^2 + fy1.^2).^2 / 4 + theta2 .* (fx2.^2 + fy2.^2).^2 / 4 - eta1 .* norm2_1 / 2 - eta2 .* norm2_2 / 2;
    L = 0.5.*phi1.*(lap_diff(lap_diff(phi1,k2),k2)+2.*a1.*lap_diff(phi1,k2)+a1.^2.*phi1)+0.5*phi2.*(lap_diff(lap_diff(phi2,k2),k2)+2*a1.*lap_diff(phi2,k2)+a1.^2*phi2)+0.5*phi1.*(lap_diff(lap_diff(phi2,k2),k2)+2*a2*lap_diff(phi2,k2)+a2.^2.*phi2)+0.5.*w0.*(fM1+fM2).^2;
    r = hx*hy*sum(sum(N+L));

    function lap=lap_diff(phi,k2)
   lap=real(ifft2((k2.*fft2(phi))));
    end
end 