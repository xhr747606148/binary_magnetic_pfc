function [phi1,phi2,M, u, xx, yy, hx, hy] = BMPFC_Spatial_1st(pde,domain,Nx,Ny,time)

epsilon = pde.epsilon;%[0,0.1]
M1 = pde.M1;
M2 = pde.M2;
S1 = pde.S1;
S2 = pde.S2;
alpha = pde.alpha;
beta = pde.beta;
eta = pde.eta;
gamma12 = pde.gamma12;
w0 =pde.w0;
gamma1=pde.gamma1;
gamma2=pde.gamma2;
eta1=pde.eta1;
eta2=pde.eta2;
theta = pde.theta;
theta1=pde.theta1;
theta2=pde.theta2;
B = pde.B;
a1 = pde.a1;
a2 = pde.a2;
a12 = pde.a12;
    
    % time
   T  = time.T;
    t  = time.t0;
    dt = time.dt;
    tsave = time.tsave;
    
    Lx = domain.right - domain.left;
    Ly = domain.top   - domain.bottom;
    hx = Lx/Nx; hy = Ly/Ny;
    x  = domain.left   + hx.*(0:Nx-1);
    y  = domain.bottom + hy.*(0:Ny-1);

    [xx,yy] = ndgrid(x,y);
    
    p = 1i.*[0:Nx/2-1 0 -Nx/2+1:-1].*(2.*pi/Lx);
    q = 1i.*[0:Ny/2-1 0 -Ny/2+1:-1].*(2.*pi/Ly);

    [pp,qq] = ndgrid(p,q);
    
    k_x =1i.*[0:Nx/2 -Nx/2+1:-1].*(2.*pi/Lx);
    k_y =1i.*[0:Ny/2 -Ny/2+1:-1].*(2.*pi/Ly);
%     k_x =[0:Nx/2 -Nx/2+1:-1].*(2.*pi/Lx);
%     k_y =[0:Ny/2 -Ny/2+1:-1].*(2.*pi/Ly);
    k2x = k_x.^2; k2y = k_y.^2;
    [kxx, kyy] = ndgrid(k2x,k2y);
    k2 = kxx + kyy; 
    
    phi1_0 = exact1(xx,yy,0);
    phi2_0 = exact2(xx,yy,0);
    M0_1 = exact3(xx,yy,0);
    M0_2 = exact4(xx,yy,0);
    M0 = [M0_1,M0_2];  
    H_1 = HH1(xx,yy);
    H_2 = HH2(yy);
    E1 = fft2(F(phi1_0,phi2_0,M0_1,M0_2,H_1,H_2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,a1,a2));
    if E1(1,1).*hx.*hy + B < 0
        disp("Root < 0");
        return;
    end
    u0 = sqrt(E1(1,1).*hx.*hy + B);


    % It is important to make the highest frequency N/2 to zero 
    % in odd derivatives due to the symmetry.
    % k_x = 1i.*[0:Nx/2-1 0 -Nx/2+1:-1].*(2.*pi/Lx);
    % k_y = 1i.*[0:Ny/2-1 0 -Ny/2+1:-1].*(2.*pi/Ly);

%     k_x = 1i.*[0:Nx/2 -Nx/2+1:-1].*(2.*pi/Lx);
%     k_y = 1i.*[0:Ny/2 -Ny/2+1:-1].*(2.*pi/Ly);
%     
%     k2x = k_x.^2; k2y = k_y.^2;
%     [kxx, kyy] = ndgrid(k2x,k2y);
%     k2 = kxx + kyy;

    nplot = round((T-t)/dt);
    nsave = round(tsave/dt);
%     t = 0;
    tstart = tic;
    for nt = 1:nplot
        t = t+dt;
        
        H1 = fun_H1(phi1_0,phi2_0,M0_1,M0_2,H_1,H_2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2);
        H1_hat = fft2(H1);
        H1 = H1 - H1_hat(1,1).*hx.*hy/(Lx.*Ly);
        
        H2 = fun_H2(phi1_0,phi2_0,M0_1,M0_2,H_1,H_2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2);
        H2_hat = fft2(H2);
        H2 = H2 - H2_hat(1,1).*hx.*hy/(Lx.*Ly);
        
        fx1 = real(ifft2(pp.*fft2(phi1_0))); fy1 =  real(ifft2(qq.*fft2(phi1_0)));
        fx2 = real(ifft2(pp.*fft2(phi2_0))); fy2 =  real(ifft2(qq.*fft2(phi2_0)));
        
%         H_1 = HH1(xx,yy);
%         H_2 = HH2(yy);
        
        R_1 = fun_R(phi1_0,phi2_0,fx1,fx2,M0_1,M0_1,M0_2,H_1,H_1,H_2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2);
        R_2 = fun_R(phi1_0,phi2_0,fy1,fy2,M0_2,M0_1,M0_2,H_2,H_1,H_2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2);
        R = [R_1,R_2];
        
%         g1_hat = (phi1_0 + phi2_0)/(dt.*M1) + rhs1(xx,yy,t,M1,epsilon,eta,eta1,gamma12,gamma1,theta1,a1,a12,hx,hy)./M1 + rhs2(xx,yy,t,hx,hy)./M1;
        g1_hat = (phi1_0 + phi2_0)/(dt.*M1) ;
%         g1_hat_g = rhs12(xx,yy,t,epsilon,gamma12,gamma1,eta1,theta1,a1,a12,eta) + rhs22(xx,yy,t,epsilon,gamma12,gamma1,eta1,theta1,a1,a12,eta);
        gphia = fft2(g1_hat);
%         gphiag=fft2(g1_hat_g);
%         ga = g1_hat + (a1.^2 + 0.5 .* a2.^2 + S) .* (dt.*M1) .* gphia(1,1).*hx.*hy/(Lx.*Ly) + rhs1(xx,yy,t,M1,epsilon,gamma)./M1 + rhs2(xx,yy,t,M1,epsilon,gamma)./M1;
        ga = g1_hat - 2.*a1.*lap_diff(phi1_0,k2)- 2.*a2.*lap_diff(phi2_0,k2) + (0.5 .* a12.^2 + S1) .* (dt.*M1) .* gphia(1,1).*hx.*hy/(Lx.*Ly) ;
        phia = inv_A(ga,dt,M1,k2,a12,S1); 
        
%         g2_hat = (phi1_0 - phi2_0)/(dt.*M1) + rhs1(xx,yy,t,M1,epsilon,eta,eta1,gamma12,gamma1,theta1,a1,a12,hx,hy)./M1 - rhs2(xx,yy,t,hx,hy)./M1;
        g2_hat = (phi1_0 - phi2_0)/(dt.*M1);
%         g2_hat_g = rhs12(xx,yy,t,epsilon,gamma12,gamma1,eta1,theta1,a2,a12,eta) - rhs22(xx,yy,t,epsilon,gamma12,gamma1,eta1,theta1,a2,a12,eta);
        gphib = fft2(g2_hat);
%         gphibg=fft2(g2_hat_g);
%         gb = g2_hat + (a1.^2 - 0.5 .* a2.^2 + S) .* (dt.*M1) .* gphib(1,1).*hx.*hy/(Lx.*Ly) + rhs1(xx,yy,t,M1,epsilon,gamma)./M1 - rhs2(xx,yy,t,M1,epsilon,gamma)./M1;
        gb = g2_hat - 2.*a1.*lap_diff(phi1_0,k2)+ 2.*a2.*lap_diff(phi2_0,k2) + (- 0.5 .* a12.^2 + S1) .* (dt.*M1) .* gphib(1,1).*hx.*hy/(Lx.*Ly);
        phib = inv_B(gb,dt,M1,k2,a12,S1);
        
        g3_hat = - H1 - H2;
        gphic = fft2(g3_hat);
        gc = g3_hat + (0.5 .* a12.^2 + S1) .* (dt.*M1) .* gphic(1,1).*hx.*hy/(Lx.*Ly);
%         gc = g3_hat ;
        phic = inv_A(gc,dt,M1,k2,a12,S1); 
        
        g4_hat = - H1 + H2;
        gphid = fft2(g4_hat);
        gd = g4_hat + (- 0.5 .* a12.^2 + S1) .* (dt.*M1) .* gphid(1,1).*hx.*hy/(Lx.*Ly);
%         gd = g4_hat ;
        phid = inv_B(gd,dt,M1,k2,a12,S1);
        
        
        f1_1 = M0_1./(M2 .* dt)   ;
        f1_2 = M0_2./(M2 .* dt)   ;
%         f1_1 = M0_1./(M2 .* dt)  ;
%         f1_2 = M0_2./(M2 .* dt)  ;
        M_11_1 = inv_C(f1_1,dt,M2,k2,w0,S2);
        M_11_2 = inv_C(f1_2,dt,M2,k2,w0,S2);
        M_11 = [M_11_1,M_11_2];
        
        f2_1 = - R_1;
        f2_2 = - R_2;
        M_12_1 = inv_C(f2_1,dt,M2,k2,w0,S2);
        M_12_2 = inv_C(f2_2,dt,M2,k2,w0,S2);
        M_12 = [M_12_1,M_12_2];
        
        phi_11 = (phia + phib)./2;
        phi_21 = (phia - phib)./2;
        phi_12 = (phic + phid)./2;
        phi_22 = (phic - phid)./2;
        
        H1phi12 = fft2(H1.*phi_12);
        H2phi22 = fft2(H2.*phi_22);
        RM12 = fft2(R_1.*M_12_1 + R_2.*M_12_2);
        A = 1- 0.5 .* H1phi12(1,1).*hx.*hy - 0.5 .* H2phi22(1,1).*hx.*hy - 0.5.* RM12(1,1).*hx.*hy;
        
        H1phi0 = fft2(H1 .* phi1_0);
        H2phi0 = fft2(H2 .* phi2_0);
        RM0 = fft2(R_1.*M0_1 + R_2.*M0_2);
        g1 = u0 - 0.5 .* H1phi0(1,1).*hx.*hy - 0.5 .* H2phi0(1,1).*hx.*hy - 0.5 .* RM0(1,1).*hx.*hy;
        
        H1phi11 = fft2(H1 .* phi_11);
        H2phi21 = fft2(H2 .* phi_21);
        RM11 = fft2(R_1.*M_11_1 + R_2.*M_11_2);
        G = g1 + 0.5 .* H1phi11(1,1).*hx.*hy + 0.5 .* H2phi21(1,1).*hx.*hy + 0.5 .* RM11(1,1).*hx.*hy;
        
        u = G./A;
        phi1 = phi_11 + u .* phi_12;
        phi2 = phi_21 + u .* phi_22;
        M_1 = M_11_1 + u .* M_12_1;
        M_2 = M_11_2 + u .* M_12_2;
        M = [M_1,M_2];
        
        phi1_0 = phi1;
        phi2_0 = phi2;
        M0_1 = M_1;
        M0_2 = M_2;
        M0 = [M0_1,M0_2];
        u0 = u;
    end
end
function r = fun_H1(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2)
    E1 = fft2(F(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,a1,a2));
    if E1(1,1).*hx.*hy + B < 0
        disp("Root < 0");
        return;
    end
    r = F_derivative_phi1(phi1,phi2,M1,M2,epsilon,eta,eta1,gamma12,gamma1,theta1,S1,pp,qq,a1) ./ sqrt(E1(1,1).*hx.*hy + B);
end
function r = fun_H2(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2)
    E1 = fft2(F(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,a1,a2));
    if E1(1,1).*hx.*hy + B < 0
        disp("Root < 0");
        return;
    end
    r = F_derivative_phi2(phi1,phi2,M1,M2,epsilon,eta,eta2,gamma12,gamma2,theta2,S1,pp,qq,a2) ./ sqrt(E1(1,1).*hx.*hy + B);
end
function r = fun_R(phi1,phi2,nabla1,nabla2,M,M1,M2,H,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,hx,hy,B,a1,a2)
    E1 = fft2(F(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,a1,a2));
    if E1(1,1).*hx.*hy + B < 0
        disp("Root < 0");
        return;
    end
    r = F_derivative_M(phi1,phi2,nabla1,nabla2,M,M1,M2,H,alpha,beta,gamma1,gamma2,theta,eta1,eta2,S2,pp,qq) ./ sqrt(E1(1,1).*hx.*hy + B);
end

function r = inv_A(phi,dt,M1,k2,a12,S1)
% global dt k2 beta 
    r = real(ifft2(fft2(phi)./(1 ./ (dt.*M1) + (k2).^2 + 0.5.*(a12 + k2).^2 + S1)));
end
function r = inv_B(phi,dt,M1,k2,a12,S1)
% global dt k2 beta 
    r = real(ifft2(fft2(phi)./(1 ./ (dt.*M1) + (k2).^2 - 0.5.*(a12 + k2).^2 + S1)));
end
function r = inv_C(phi,dt,M2,k2,w0,S2)
% global dt k2 beta 
    r = real(ifft2(fft2(phi)./(1 ./ (dt.*M2) - w0.*k2 + S2)));
end

function r = F(phi1,phi2,M1,M2,H1,H2,epsilon,alpha,beta,eta,eta1,eta2,gamma12,gamma1,gamma2,theta,theta1,theta2,S1,S2,pp,qq,a1,a2)
    phi1 = real(phi1); 
    fx1 = real(ifft2(pp.*fft2(phi1))); fy1 =  real(ifft2(qq.*fft2(phi1)));
    phi2 = real(phi2); 
    fx2 = real(ifft2(pp.*fft2(phi2))); fy2 =  real(ifft2(qq.*fft2(phi2)));
    norm1 = M1.^2 + M2.^2;
    norm2_1 = (M1 .* fx1 + M2 .* fy1).^2;
    norm2_2 = (M1 .* fx2 + M2 .* fy2).^2;
    normH=(M1.*H1 + M2.*H2);
    r = -normH + 0.5.*a1.*phi1.^2 + 0.5.*a2.*phi2.^2 + phi1.^4 / 4  - epsilon .* phi1.^2 / 2 + phi2.^4 / 4  - epsilon .* phi2.^2 / 2 + eta .* (abs(phi1).^3 - phi1.^3) / 3 + eta .* (abs(phi2).^3 - phi2.^3) / 3 + (gamma12 .* phi1.^2 .* phi2.^2)./2 - 0.5 .* S1 .* phi1.^2 - 0.5 .* S1 .* phi2.^2 - 0.5 .* S2 .*norm1 + beta .* norm1.^2 / 4  - alpha .* norm1 / 2 - gamma1 .* phi1 .* norm1 - gamma2 .* phi2 .* norm1 + theta .* norm1.^3 / 6 + theta1 .* (fx1.^2 + fy1.^2).^2 / 4 + theta2 .* (fx2.^2 + fy2.^2).^2 / 4 - eta1 .* norm2_1 / 2 - eta2 .* norm2_2 / 2;
%     r = phi1.^4 / 4  - epsilon .* phi1.^2 / 2 + phi2.^4 / 4  - epsilon .* phi2.^2 / 2  + (gamma .* phi1.^2 .* phi2.^2)./2 - 0.5 .* S .* phi1.^2 - 0.5 .* S .* phi2.^2;
end
function r = F_derivative_phi1(phi1,phi2,M1,M2,epsilon,eta,eta1,gamma12,gamma1,theta1,S1,pp,qq,a1)
    phi1 = real(phi1); 
    fx = real(ifft2(pp.*fft2(phi1))); fy =  real(ifft2(qq.*fft2(phi1)));
    norm1 = M1.^2 + M2.^2;
    f1 = (fx.^2 + fy.^2).*fx ;  f2 = (fx.^2 + fy.^2).*fy;
    div1 = real(ifft2(pp.*fft2(f1) + qq.*fft2(f2)));

    f3 = (M1 .* fx + M2 .* fy).*M1;  f4 = (M1 .* fx + M2 .* fy).*M2;
    div2 = real(ifft2(pp.*fft2(f3) + qq.*fft2(f4)));
    
    r = a1.*phi1 + phi1 .^3 - epsilon .* phi1 + eta .* (abs(phi1) - phi1) .* phi1 + gamma12 .* phi1 .* phi2.^2 - S1 .* phi1 - gamma1 .* norm1 - theta1 .* div1 + eta1 .* div2;
%     r = phi1 .^3 - epsilon .* phi1  + gamma .* phi1 .* phi2.^2 - S .* phi1;
end
function r = F_derivative_phi2(phi1,phi2,M1,M2,epsilon,eta,eta2,gamma12,gamma2,theta2,S1,pp,qq,a2)
    phi2 = real(phi2); 
    fx = real(ifft2(pp.*fft2(phi2))); fy =  real(ifft2(qq.*fft2(phi2)));
    norm1 = M1.^2 + M2.^2;
    f1 = (fx.^2 + fy.^2).*fx ;  f2 = (fx.^2 + fy.^2).*fy;
    div1 = real(ifft2(pp.*fft2(f1) + qq.*fft2(f2)));

    f3 = (M1 .* fx + M2 .* fy).*M1;  f4 = (M1 .* fx + M2 .* fy).*M2;
    div2 = real(ifft2(pp.*fft2(f3) + qq.*fft2(f4)));
    
    r = a2.*phi2 + phi2 .^3 - epsilon .* phi2 + eta .* (abs(phi2) - phi2) .* phi2 + gamma12 .* phi2 .* phi1.^2 - S1 .* phi2 - gamma2 .* norm1 - theta2 .* div1 + eta2 .* div2;
%      r = phi2 .^3 - epsilon .* phi2  + gamma .* phi2 .* phi1.^2 - S .* phi2;
end
function r = F_derivative_M(phi1,phi2,nabla1,nabla2,M,M1,M2,H,alpha,beta,gamma1,gamma2,theta,eta1,eta2,S2,pp,qq)
    phi1 = real(phi1); 
    fx1 = real(ifft2(pp.*fft2(phi1))); fy1 =  real(ifft2(qq.*fft2(phi1)));
    phi2 = real(phi2); 
    fx2 = real(ifft2(pp.*fft2(phi2))); fy2 =  real(ifft2(qq.*fft2(phi2)));
    norm1 = M1.^2 + M2.^2;
    norm2_1 = (M1 .* fx1 + M2 .* fy1);
    norm2_2 = (M1 .* fx2 + M2 .* fy2);
    r = -S2.*M + beta .* norm1.*M - alpha .* M - H - 2 .* gamma1 .* M .* phi1 - 2 .* gamma2 .* M .* phi2 + theta .* (norm1.^2) .* M - eta1 .* norm2_1 .* nabla1 - eta2 .* norm2_2 .* nabla2;
%     r = beta .* norm1.*M - alpha .* M - 2 .* gamma .* M .* phi + theta1 .* (norm1.^2) .* M ;
end
function z = exact1(x,y,t)  
       z = cos(2.*x).*cos(2.*y).*cos(t);
end
function z = exact2(x,y,t)  
       z = cos(x).*cos(y).*cos(t);
end
function z = exact3(x,y,t)  
       z = sin(4*x).*sin(4.*y).*cos(t);
end
function z = exact4(x,y,t)  
       z = cos(4.*x).*cos(4.*y).*cos(t);
end
function z = HH1(x,y)  
       z = sin(4.*x).*cos(4.*y);
end
function z = HH2(y)  
       z = cos(2.*y);
end
% function z = exact1(x,y,t)  
%        z = cos(2.*pi.*x/32).*cos(2.*pi.*y/32).*cos(t);
% end
% function z = exact2(x,y,t)  
%        z = cos(pi.*x/32).*cos(pi.*y/32).*cos(t);
% end
% function z = exact3(x,y,t)  
%        z = sin(8.*pi.*x/128).*sin(8.*pi.*y/128).*cos(t);
% end
% function z = exact4(x,y,t)  
%        z = cos(8.*pi.*x/128).*cos(8.*pi.*y/128).*cos(t);
% end
% function z = HH1(x,y)  
%        z = sin(8.*pi.*x/128).*cos(8.*pi.*y/256);
% end
% function z = HH2(y)  
%        z = cos(8.*pi.*y/256);
% end
% function z = rhs11(x,y,t,M1) 
% z = -cos((x.*pi)./16).*cos((y.*pi)./16).*(sin(t) - M1.*cos(t));
% end
% function z = rhs12(x,y,t,epsilon,gamma12,gamma1,eta1,theta1,a1,a12,eta)
% z=cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t) - cos((x.*pi)./16).^3.*cos((y.*pi)./16).^3.*cos(t).^3 - a1.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t) - (a12.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t))./2 + epsilon.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t) + eta.*cos((x.*pi)./16).^2.*cos((y.*pi)./16).^2.*cos(t).^2 + gamma1.*cos((x.*pi)./16).^2.*cos((y.*pi)./16).^2.*cos(t).^2 + gamma1.*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^2 - (pi.^4.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t))./16384 - (pi.^4.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t))./524288 + (a1.*pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t))./64 + (a12.*pi.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t))./512 - eta.*cos((x.*pi)./16).*cos((y.*pi)./16).*abs(cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t)).*cos(t) + (eta1.*pi.^2.*cos((x.*pi)./16).^3.*cos((y.*pi)./16).^3.*cos(t).^3)./256 + (eta1.*pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).^3.*sin((x.*pi)./16).^2.*cos(t).^3)./256 - (eta1.*pi.^2.*cos((x.*pi)./16).^3.*cos((y.*pi)./16).*sin((y.*pi)./16).^2.*cos(t).^3)./256 - (theta1.*pi.^4.*cos((x.*pi)./16).*cos((y.*pi)./16).^3.*sin((x.*pi)./16).^2.*cos(t).^3)./16384 - (theta1.*pi.^4.*cos((x.*pi)./16).^3.*cos((y.*pi)./16).*sin((y.*pi)./16).^2.*cos(t).^3)./16384 - gamma12.*cos((x.*pi)./16).*cos((x.*pi)./32).^2.*cos((y.*pi)./16).*cos((y.*pi)./32).^2.*cos(t).^3 - (eta1.*pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^3)./256 + (theta1.*pi.^4.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^3)./16384;
% end
% function z = rhs21(x,y,t,M1)  
% z=-cos((x.*pi)./32).*cos((y.*pi)./32).*(sin(t) - M1.*cos(t));
% end
% function z = rhs22(x,y,t,epsilon,gamma12,gamma2,eta2,theta2,a2,a12,eta)  
% z=gamma2.*cos(t).^2 - (a12.^2.*cos(t))./2 - (pi.^4.*cos(t))./32768 + (a12.*pi.^2.*cos(t))./128 + cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t) - cos((x.*pi)./32).^3.*cos((y.*pi)./32).^3.*cos(t).^3 + (pi.^4.*cos((x.*pi)./32).^2.*cos(t))./16384 + (pi.^4.*cos((y.*pi)./32).^2.*cos(t))./16384 + a12.^2.*cos((x.*pi)./32).^2.*cos(t) + a12.^2.*cos((y.*pi)./32).^2.*cos(t) - 4.*gamma2.*cos((x.*pi)./32).^2.*cos(t).^2 + 4.*gamma2.*cos((x.*pi)./32).^4.*cos(t).^2 - 4.*gamma2.*cos((y.*pi)./32).^2.*cos(t).^2 + 4.*gamma2.*cos((y.*pi)./32).^4.*cos(t).^2 - a2.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t) - gamma12.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t).^3 - (a12.*pi.^2.*cos((x.*pi)./32).^2.*cos(t))./64 - (a12.*pi.^2.*cos((y.*pi)./32).^2.*cos(t))./64 + 4.*gamma12.*cos((x.*pi)./32).*cos((y.*pi)./32).^3.*cos(t).^3 + 4.*gamma12.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).*cos(t).^3 - 4.*gamma12.*cos((x.*pi)./32).*cos((y.*pi)./32).^5.*cos(t).^3 - 4.*gamma12.*cos((x.*pi)./32).^5.*cos((y.*pi)./32).*cos(t).^3 + epsilon.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t) - (pi.^4.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*cos(t))./8192 - 2.*a12.^2.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*cos(t) + eta.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*cos(t).^2 + 32.*gamma2.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*cos(t).^2 - 32.*gamma2.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^4.*cos(t).^2 - 32.*gamma2.*cos((x.*pi)./32).^4.*cos((y.*pi)./32).^2.*cos(t).^2 + 32.*gamma2.*cos((x.*pi)./32).^4.*cos((y.*pi)./32).^4.*cos(t).^2 - 16.*gamma12.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).^3.*cos(t).^3 + 16.*gamma12.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).^5.*cos(t).^3 + 16.*gamma12.*cos((x.*pi)./32).^5.*cos((y.*pi)./32).^3.*cos(t).^3 - 16.*gamma12.*cos((x.*pi)./32).^5.*cos((y.*pi)./32).^5.*cos(t).^3 - (pi.^4.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t))./262144 + (a12.*pi.^2.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*cos(t))./32 + (eta2.*pi.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).^3.*cos(t).^3)./128 + (eta2.*pi.^2.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).*cos(t).^3)./32 - (eta2.*pi.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).^5.*cos(t).^3)./256 - (5.*eta2.*pi.^2.*cos((x.*pi)./32).^5.*cos((y.*pi)./32).*cos(t).^3)./256 - (theta2.*pi.^4.*cos((x.*pi)./32).*cos((y.*pi)./32).^3.*cos(t).^3)./131072 - (theta2.*pi.^4.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).*cos(t).^3)./131072 + (a2.*pi.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t))./256 - eta.*cos((x.*pi)./32).*cos((y.*pi)./32).*abs(cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t)).*cos(t) - (5.*eta2.*pi.^2.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).^3.*cos(t).^3)./128 + (eta2.*pi.^2.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).^5.*cos(t).^3)./64 + (eta2.*pi.^2.*cos((x.*pi)./32).^5.*cos((y.*pi)./32).^3.*cos(t).^3)./64 + (3.*theta2.*pi.^4.*cos((x.*pi)./32).^3.*cos((y.*pi)./32).^3.*cos(t).^3)./262144 - (7.*eta2.*pi.^2.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t).^3)./1024 + (theta2.*pi.^4.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t).^3)./262144;
% end
% function z = rhs3(x,y,t,w0,alpha,beta,gamma1,gamma2,eta1,eta2,theta,M2)  
%   z=- M2.*(cos((y.*pi)./32).*sin((x.*pi)./16) + alpha.*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t) - beta.*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).*(cos((x.*pi)./16).^2.*cos((y.*pi)./16).^2.*cos(t).^2 + sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^2) - (w0.*pi.^2.*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t))./128 - theta.*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).*(cos((x.*pi)./16).^2.*cos((y.*pi)./16).^2.*cos(t).^2 + sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^2).^2 + 2.*gamma1.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^2 + 2.*gamma2.*cos((x.*pi)./32).*cos((y.*pi)./32).*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^2 + (eta1.*pi.*cos((y.*pi)./16).*sin((x.*pi)./16).*cos(t).*((pi.*cos((y.*pi)./16).*sin((x.*pi)./16).^2.*sin((y.*pi)./16).*cos(t).^2)./16 + (pi.*cos((x.*pi)./16).^2.*cos((y.*pi)./16).*sin((y.*pi)./16).*cos(t).^2)./16))./16 + (eta2.*pi.*cos((y.*pi)./32).*sin((x.*pi)./32).*cos(t).*((pi.*cos((x.*pi)./16).*cos((x.*pi)./32).*cos((y.*pi)./16).*sin((y.*pi)./32).*cos(t).^2)./32 + (pi.*cos((y.*pi)./32).*sin((x.*pi)./16).*sin((x.*pi)./32).*sin((y.*pi)./16).*cos(t).^2)./32))./32) - sin((x.*pi)./16).*sin((y.*pi)./16).*sin(t);
% end
% function z = rhs4(x,y,t,w0,alpha,beta,gamma1,gamma2,eta1,eta2,theta,M2)  
%   z=M2.*beta.*cos((x.*pi)./16).^3.*cos((y.*pi)./16).^3.*cos(t).^3 - cos((x.*pi)./16).*cos((y.*pi)./16).*sin(t) - M2.*cos((y.*pi)./32) - 2.*M2.*gamma1.*cos((x.*pi)./16).^2.*cos((y.*pi)./16).^2.*cos(t).^2 + M2.*theta.*cos((x.*pi)./16).^5.*cos((y.*pi)./16).^5.*cos(t).^5 - M2.*alpha.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t) + (M2.*w0.*pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t))./128 + 2.*M2.*theta.*cos((x.*pi)./16).^3.*cos((y.*pi)./16).^3.*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^5 + M2.*beta.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^3 + M2.*theta.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).^4.*sin((y.*pi)./16).^4.*cos(t).^5 - (M2.*eta1.*pi.^2.*cos((x.*pi)./16).^3.*cos((y.*pi)./16).*sin((y.*pi)./16).^2.*cos(t).^3)./256 - 2.*M2.*gamma2.*cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t).^2.*(2.*cos((x.*pi)./32).^2 - 1).*(2.*cos((y.*pi)./32).^2 - 1) - (M2.*eta1.*pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2.*cos(t).^3)./256 - (M2.*eta2.*pi.^2.*cos((x.*pi)./32).^2.*sin((y.*pi)./32).^2.*cos(t).^3.*(2.*cos((x.*pi)./32).^2 - 1).*(2.*cos((y.*pi)./32).^2 - 1))./1024 - (M2.*eta2.*pi.^2.*cos((x.*pi)./32).^2.*cos((y.*pi)./32).^2.*sin((x.*pi)./32).^2.*sin((y.*pi)./32).^2.*cos(t).^3)./256;
% end
function lap=lap_diff(phi,k2)
   lap=real(ifft2((k2.*fft2(phi))));
end

