clear; clc; close all;

% Parameters
para.epsilon = 0.1;%[0,0.1]
para.M1 = 0.01;
para.M2 = 0.01;
para.S1 = 10;
para.S2 = 10;
para.alpha = 1;
para.beta = 10;
para.eta = 10;
para.gamma12 = 0.5;
para.w0 =1;
para.gamma1=1;
para.gamma2=1;
para.eta1=-100;
para.eta2=-100;
para.theta = 1e-9;
para.theta1=1e-9;
para.theta2=1e-9;
para.B = 1e7;
para.a1 = 1;
para.a2 = 1;
para.a12 = 1.2;

% Space: Domain and N
domain.left   = 0; domain.right  = 2.*pi;
domain.bottom = 0; domain.top    = 2.*pi;
Nx = 128; Ny = Nx;

% Time: dt T
T = 0.5;
t0 = 0;
tsave = T;
pde = ex_BMPFC_data(para);
%% time step
% dt_array = 1e-3.*[ 16 8 4 2 1]';
% dt_array = 1./[ 5 10 20 100 200 500 1000]';
% dt_array = 1./[ 1 2 4 8 16 32 64 128 256 512]';
% dt_array = 1e-9.*[ 1 ]';
% maxIt = length(dt_array);
% exact_phi = cell(0,1);
% exact_M = cell(0,1);
%% Run:
% error=zeros(maxIt,1); order=zeros(maxIt,1);
% for k=1

    dt = 1e-4;

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
    E1 = fft2(F(phi1_0,phi2_0,M0_1,M0_2,H_1,H_2,para.epsilon,para.alpha,para.beta,para.gamma12,para.gamma1,para.gamma2,para.eta,para.eta1,para.eta2,para.theta,para.theta1,para.theta2,para.S1,para.S2,pp,qq,para.a1,para.a2));
    if E1(1,1).*hx.*hy + para.B < 0
        disp("Root < 0");
        return;
    end
    u0 = sqrt(E1(1,1).*hx.*hy + para.B);
    

    time = struct('T',T,'t0',t0,'dt',dt,'tsave',tsave);
    time1 = struct('T',dt,'t0',t0,'dt',dt,'tsave',tsave);
    [phi1_1, phi2_1,M11, u1, ~] = BMPFC_1st(pde,domain,Nx,Ny,time1);

    nplot = round((T-t0)/dt);
    nsave = round(tsave/dt);
    t = dt;
    tstart = tic;
    for nt = 2:nplot
        t = t+dt;

        phi1_bar = 2 .* phi1_1 - phi1_0;
        phi2_bar = 2 .* phi2_1 - phi2_0;
        M11_1=M11(:,1:size(M11,2)/2);
        M11_2=M11(:,size(M11,2)/2+1:end);
        M_1_bar = 2.*M11_1-M0_1;
        M_2_bar = 2.*M11_2-M0_2;
        M_bar = [M_1_bar,M_2_bar];   
        H_1 = HH1(xx,yy);
        H_2 = HH2(yy);
        H1 = fun_H1(phi1_bar,phi2_bar,M_1_bar,M_2_bar,H_1,H_2,para.epsilon,para.alpha,para.beta,para.gamma12,para.gamma1,para.gamma2,para.eta,para.eta1,para.eta2,para.theta,para.theta1,para.theta2,para.S1,para.S2,pp,qq,hx,hy,para.B,para.a1,para.a2);
        H1_hat = fft2(H1);
        H1 = H1 - H1_hat(1,1).*hx.*hy/(Lx.*Ly);
        
        H2 = fun_H2(phi1_bar,phi2_bar,M_1_bar,M_2_bar,H_1,H_2,para.epsilon,para.alpha,para.beta,para.gamma12,para.gamma1,para.gamma2,para.eta,para.eta1,para.eta2,para.theta,para.theta1,para.theta2,para.S1,para.S2,pp,qq,hx,hy,para.B,para.a1,para.a2);
        H2_hat = fft2(H2);
        H2 = H2 - H2_hat(1,1).*hx.*hy/(Lx.*Ly);
        
        fx1 = real(ifft2(pp.*fft2(phi1_bar))); fy1 =  real(ifft2(qq.*fft2(phi1_bar)));
        fx2 = real(ifft2(pp.*fft2(phi2_bar))); fy2 =  real(ifft2(qq.*fft2(phi2_bar)));
        
        R_1 = fun_R(phi1_bar,phi2_bar,fx1,fx2,M_1_bar,M_1_bar,M_2_bar,H_1,H_1,H_2,para.epsilon,para.alpha,para.beta,para.gamma12,para.gamma1,para.gamma2,para.eta,para.eta1,para.eta2,para.theta,para.theta1,para.theta2,para.S1,para.S2,pp,qq,hx,hy,para.B,para.a1,para.a2);
        R_2 = fun_R(phi1_bar,phi2_bar,fy1,fy2,M_2_bar,M_1_bar,M_2_bar,H_2,H_1,H_2,para.epsilon,para.alpha,para.beta,para.gamma12,para.gamma1,para.gamma2,para.eta,para.eta1,para.eta2,para.theta,para.theta1,para.theta2,para.S1,para.S2,pp,qq,hx,hy,para.B,para.a1,para.a2);
        R = [R_1,R_2];

        g1_hat = (4 .* phi1_1 - phi1_0 + 4 .* phi2_1 - phi2_0)/(2 .* dt.*para.M1);
        gphia = fft2(g1_hat);
%         ga = g1_hat + (a1.^2 + 0.5 .* a2.^2 + S) .* (dt.*M1) .* gphia(1,1).*hx.*hy/(Lx.*Ly) + rhs1(xx,yy,t,M1,epsilon,gamma)./M1 + rhs2(xx,yy,t,M1,epsilon,gamma)./M1;
        ga = g1_hat - 2.*para.a1.*lap_diff(phi1_bar,k2) - 2.*para.a2.*lap_diff(phi2_bar,k2) + (0.5 .* para.a12.^2 + para.S1) .* (2.*dt.*para.M1) .* gphia(1,1).*hx.*hy/(3 .* (Lx.*Ly));
        phia = inv_A(ga,dt,para.M1,k2,para.a12,para.S1); 
        
        g2_hat = (4 .* phi1_1 - phi1_0 - 4 .* phi2_1 + phi2_0)/(2 .* dt.*para.M1);
        gphib = fft2(g2_hat);
%         gb = g2_hat + (a1.^2 - 0.5 .* a2.^2 + S) .* (dt.*M1) .* gphib(1,1).*hx.*hy/(Lx.*Ly) + rhs1(xx,yy,t,M1,epsilon,gamma)./M1 - rhs2(xx,yy,t,M1,epsilon,gamma)./M1;
        gb = g2_hat - 2.*para.a1.*lap_diff(phi1_bar,k2) + 2.*para.a2.*lap_diff(phi2_bar,k2) + (- 0.5 .* para.a12.^2 + para.S1) .* (2.*dt.*para.M1) .* gphib(1,1).*hx.*hy/(3 .* (Lx.*Ly));
        phib = inv_B(gb,dt,para.M1,k2,para.a12,para.S1);
        
        g3_hat = - H1 - H2;
        gphic = fft2(g3_hat);
        gc = g3_hat + (0.5 .* para.a2.^2 + para.S1) .* (2.*dt.*para.M1) .* gphic(1,1).*hx.*hy/(3 .* (Lx.*Ly));
%         gc = g3_hat ;
        phic = inv_A(gc,dt,para.M1,k2,para.a12,para.S1); 
        
        g4_hat = - H1 + H2;
        gphid = fft2(g4_hat);
        gd = g4_hat + (- 0.5 .* para.a2.^2 + para.S1) .* (2.*dt.*para.M1) .* gphid(1,1).*hx.*hy/(3 .* (Lx.*Ly));
%         gd = g4_hat ;
        phid = inv_B(gd,dt,para.M1,k2,para.a12,para.S1);
        
        
        
        f1_1 = (4.*M11_1 - M0_1)./(2.* para.M2 .* dt)   ;
        f1_2 = (4.*M11_2 - M0_2)./(2.* para.M2 .* dt)   ;
        M_11_1 = inv_C(f1_1,dt,para.M2,k2,para.w0,para.S2);
        M_11_2 = inv_C(f1_2,dt,para.M2,k2,para.w0,para.S2);
        M_11 = [M_11_1,M_11_2];
        
        f2_1 = - R_1;
        f2_2 = - R_2;
        M_12_1 = inv_C(f2_1,dt,para.M2,k2,para.w0,para.S2);
        M_12_2 = inv_C(f2_2,dt,para.M2,k2,para.w0,para.S2);
        M_12 = [M_12_1,M_12_2];
        
        phi_11 = (phia + phib)./2;
        phi_21 = (phia - phib)./2;
        phi_12 = (phic + phid)./2;
        phi_22 = (phic - phid)./2;
        
        H1phi12 = fft2(H1.*3.*phi_12);
        H2phi22 = fft2(H2.*3.*phi_22);
        RM12 = fft2(R_1.*3.*M_12_1 + R_2.*3.*M_12_2);
        A = 3- 0.5 .* H1phi12(1,1).*hx.*hy - 0.5 .* H2phi22(1,1).*hx.*hy - 0.5.* RM12(1,1).*hx.*hy;
        
        H1phi0 = fft2(H1 .* (4.*phi1_1 - phi1_0));
        H2phi0 = fft2(H2 .* (4.*phi2_1 - phi2_0));
        RM0 = fft2(R_1.*(4.*M11_1 - M0_1) + R_2.*(4.*M11_2 - M0_2));
        g1 = 4.*u1 - u0 - 0.5 .* H1phi0(1,1).*hx.*hy - 0.5 .* H2phi0(1,1).*hx.*hy - 0.5 .* RM0(1,1).*hx.*hy;
        
        H1phi11 = fft2(H1 .* 3 .* phi_11);
        H2phi21 = fft2(H2 .* 3 .* phi_21);
        RM11 = fft2(R_1.*3.*M_11_1 + R_2.*3.*M_11_2);
        G = g1 + 0.5 .* H1phi11(1,1).*hx.*hy + 0.5 .* H2phi21(1,1).*hx.*hy + 0.5 .* RM11(1,1).*hx.*hy;
        
        u = G./A;
        phi1 = phi_11 + u .* phi_12;
        phi2 = phi_21 + u .* phi_22;
        M_1 = M_11_1 + u .* M_12_1;
        M_2 = M_11_2 + u .* M_12_2;
        M = [M_1,M_2];
        
        phi1_0 = phi1_1;
        phi1_1 = phi1;
        phi2_0 = phi2_1;
        phi2_1 = phi2;
        M0_1 = M11_1;
        M0_2 = M11_2;
        M0 = [M0_1,M0_2];
        M11_1 = M_1;
        M11_2 = M_2;
        M11 = [M11_1,M11_2];
        u0 = u1;
        u1 = u;
 
        if  0 == mod(nt,100)
            timeElapsed = toc(tstart)
%             fprintf('epsilon=%.3f,t=%.4f/%.3f, dt=%.2e, Nx=%d, Ny=%d, timeElapsed=%f\n',epsilon,t,T,dt,Nx,Ny,timeElapsed);
%             s = pcolor(xx,yy,phi);
%             s.FaceColor='interp';
%             s.EdgeColor='interp';   
%             colorbar;
%             axis square;
%             drawnow;
        end
    end
exact_phi1 = phi1;
exact_phi2 = phi2;
exact_M = M;
save('exact_BMPFC_BDF2(eta=-100).mat','exact_phi1','exact_phi2','exact_M')
%     err1 = fft2((exact1(xx,yy,T) - phi).^2);  % L2
%     error1(k,1) = sqrt(err1(1,1) .* hx .* hy);  % L2
%     err2 = fft2(([exact2(xx,yy,T),exact3(xx,yy,T)] - M).^2);  % L2
%     error2(k,1) = sqrt(err2(1,1) .* hx .* hy);  % L2
% end
% order1(2:maxIt) = log(error1(1:maxIt-1)./error1(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
% order2(2:maxIt) = log(error2(1:maxIt-1)./error2(2:maxIt))./log(dt_array(1:maxIt-1)./dt_array(2:maxIt));
% %% Display error and order
% fprintf('    dt    &  Error_phi  &  Order_phi  &  Error_M  &  Order_M \n');
% for k = 1:maxIt
%     fprintf('%.4e   %.4e     %.4f      %.4e    %.4f \n',dt_array(k),error1(k),order1(k),error2(k),order2(k));
% end

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
    r = real(ifft2(fft2(phi)./(3 ./ (2.*dt.*M1) + k2.^2 + 0.5.*(a12 + k2).^2 + S1)));
end
function r = inv_B(phi,dt,M1,k2,a12,S1)
% global dt k2 beta 
    r = real(ifft2(fft2(phi)./(3 ./ (2.*dt.*M1) + k2.^2 - 0.5.*(a12 + k2).^2 + S1)));
end
function r = inv_C(phi,dt,M2,k2,w0,S2)
% global dt k2 beta 
    r = real(ifft2(fft2(phi)./(3 ./ (2.*dt.*M2) - w0.*k2 +S2)));
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
% function z = rhs1(x,y,t,M1,epsilon,eta,eta1,gamma12,gamma1,theta1,a1,a12,hx,hy) 
% Mat=eta*cos((x*pi)./16)*cos((y*pi)./16)*cos(t)*(abs(cos((x*pi)./16)*cos((y*pi)./16)*cos(t)) - cos((x*pi)./16)*cos((y*pi)./16)*cos(t)) -( 64*gamma1*cos((y*pi)./16).^2*cos(t).^2 - 64*gamma1*sin((y*pi)./16).^2*cos(t).^2 + 32*gamma12*cos((y*pi)./16)*cos((y*pi)./32).^2*cos(t).^3)./128;
%   mu=fft2(Mat);
% z = M1*(theta1*((pi.^2*cos((x*pi)./16)*cos((y*pi)./16)*cos(t)*((pi.^2*cos((x*pi)./16).^2*sin((y*pi)./16).^2*cos(t).^2)./256 + (pi.^2*cos((y*pi)./16).^2*sin((x*pi)./16).^2*cos(t).^2)./256))./128 + (pi*cos((y*pi)./16)*sin((x*pi)./16)*cos(t)*((pi.^3*cos((x*pi)./16)*cos((y*pi)./16).^2*sin((x*pi)./16)*cos(t).^2)./2048 - (pi.^3*cos((x*pi)./16)*sin((x*pi)./16)*sin((y*pi)./16).^2*cos(t).^2)./2048))./16 + (pi*cos((x*pi)./16)*sin((y*pi)./16)*cos(t)*((pi.^3*cos((x*pi)./16).^2*cos((y*pi)./16)*sin((y*pi)./16)*cos(t).^2)./2048 - (pi.^3*cos((y*pi)./16)*sin((x*pi)./16).^2*sin((y*pi)./16)*cos(t).^2)./2048))./16) - gamma1*(cos((x*pi)./16).^2*cos((y*pi)./16).^2*cos(t).^2 + sin((x*pi)./16).^2*sin((y*pi)./16).^2*cos(t).^2) + cos((x*pi)./16).^3*cos((y*pi)./16).^3*cos(t).^3 + a1.^2*cos((x*pi)./16)*cos((y*pi)./16)*cos(t) + (a12.^2*cos((x*pi)./32)*cos((y*pi)./32)*cos(t))./2 - epsilon*cos((x*pi)./16)*cos((y*pi)./16)*cos(t) + (pi.^4*cos((x*pi)./16)*cos((y*pi)./16)*cos(t))./16384 + (pi.^4*cos((x*pi)./32)*cos((y*pi)./32)*cos(t))./524288 + eta*cos((x*pi)./16)*cos((y*pi)./16)*cos(t)*(abs(cos((x*pi)./16)*cos((y*pi)./16)*cos(t)) - cos((x*pi)./16)*cos((y*pi)./16)*cos(t)) - (a1*pi.^2*cos((x*pi)./16)*cos((y*pi)./16)*cos(t))./64 - (a12*pi.^2*cos((x*pi)./32)*cos((y*pi)./32)*cos(t))./512 - eta1*cos((x*pi)./16)*cos((y*pi)./16)*cos(t)*((pi.^2*cos((x*pi)./16).^2*cos((y*pi)./16).^2*cos(t).^2)./256 - (pi.^2*cos((x*pi)./16).^2*sin((y*pi)./16).^2*cos(t).^2)./256 + (pi.^2*cos((y*pi)./16).^2*sin((x*pi)./16).^2*cos(t).^2)./256 - (pi.^2*sin((x*pi)./16).^2*sin((y*pi)./16).^2*cos(t).^2)./256) + gamma12*cos((x*pi)./16)*cos((x*pi)./32).^2*cos((y*pi)./16)*cos((y*pi)./32).^2*cos(t).^3) - (M1*mu(1,1)*hx*hy)./16384 - cos((x*pi)./16)*cos((y*pi)./16)*sin(t);
% end
% function z = rhs2(x,y,t,hx,hy)  
% Mat=cos((x*pi)./32)*cos((y*pi)./32)*cos(t)*(10*abs(cos((x*pi)./32)*cos((y*pi)./32)*cos(t)) - 10*cos((x*pi)./32)*cos((y*pi)./32)*cos(t)) - ((16*sin((y*pi)./16).^2*cos(t).^2)./25 - (16*cos((y*pi)./16).^2*cos(t).^2)./25)/128;
%  mu=fft2(Mat);
% z=(18*cos((x*pi)./16)*cos((y*pi)./16)*cos(t))./25 - mu(1,1)*hx*hy./16384 + (9*cos((x*pi)./32)*cos((y*pi)./32)*cos(t))./10 - cos((x*pi)./32)*cos((y*pi)./32)*sin(t) - (cos((x*pi)./16).^2*cos((y*pi)./16).^2*cos(t).^2)./100 + cos((x*pi)./32).^3*cos((y*pi)./32).^3*cos(t).^3 - (sin((x*pi)./16).^2*sin((y*pi)./16).^2*cos(t).^2)./100 + cos((x*pi)./32)*cos((y*pi)./32)*cos(t)*(10*abs(cos((x*pi)./32)*cos((y*pi)./32)*cos(t)) - 10*cos((x*pi)./32)*cos((y*pi)./32)*cos(t)) - (sin((x*pi)./16)*sin((y*pi)./16)*cos(t)*((pi.^2*cos((x*pi)./16)*cos((y*pi)./16)*sin((x*pi)./32)*sin((y*pi)./32)*cos(t).^2)./1024 - (pi.^2*cos((x*pi)./16)*cos((y*pi)./32)*sin((x*pi)./32)*sin((y*pi)./16)*cos(t).^2)./512 + (pi.^2*cos((x*pi)./32)*cos((y*pi)./16)*sin((x*pi)./16)*sin((y*pi)./32)*cos(t).^2)./512 - (pi.^2*cos((x*pi)./32)*cos((y*pi)./32)*sin((x*pi)./16)*sin((y*pi)./16)*cos(t).^2)./1024))./10 - (3*pi.^2*cos((x*pi)./16)*cos((y*pi)./16)*cos(t))./320 + (pi.^4*cos((x*pi)./16)*cos((y*pi)./16)*cos(t))./32768 - (pi.^2*cos((x*pi)./32)*cos((y*pi)./32)*cos(t))./256 + (pi.^4*cos((x*pi)./32)*cos((y*pi)./32)*cos(t))./262144 + (cos((x*pi)./16)*cos((y*pi)./16)*cos(t)*((pi.^2*cos((x*pi)./16)*cos((x*pi)./32)*cos((y*pi)./16)*cos((y*pi)./32)*cos(t).^2)./1024 - (pi.^2*cos((x*pi)./16)*cos((x*pi)./32)*sin((y*pi)./16)*sin((y*pi)./32)*cos(t).^2)./512 + (pi.^2*cos((y*pi)./16)*cos((y*pi)./32)*sin((x*pi)./16)*sin((x*pi)./32)*cos(t).^2)./512 - (pi.^2*sin((x*pi)./16)*sin((x*pi)./32)*sin((y*pi)./16)*sin((y*pi)./32)*cos(t).^2)./1024))./10 + (cos((x*pi)./16).^2*cos((x*pi)./32)*cos((y*pi)./16).^2*cos((y*pi)./32)*cos(t).^3)./2 + (pi.^2*cos((x*pi)./32)*cos((y*pi)./32)*cos(t)*((pi.^2*cos((x*pi)./32).^2*sin((y*pi)./32).^2*cos(t).^2)./1024 + (pi.^2*cos((y*pi)./32).^2*sin((x*pi)./32).^2*cos(t).^2)./1024))./512000000000 + (pi*cos((y*pi)./32)*sin((x*pi)./32)*cos(t)*((pi.^3*cos((x*pi)./32)*cos((y*pi)./32).^2*sin((x*pi)./32)*cos(t).^2)./16384 - (pi.^3*cos((x*pi)./32)*sin((x*pi)./32)*sin((y*pi)./32).^2*cos(t).^2)./16384))./32000000000 + (pi*cos((x*pi)./32)*sin((y*pi)./32)*cos(t)*((pi.^3*cos((x*pi)./32).^2*cos((y*pi)./32)*sin((y*pi)./32)*cos(t).^2)./16384 - (pi.^3*cos((y*pi)./32)*sin((x*pi)./32).^2*sin((y*pi)./32)*cos(t).^2)./16384))./32000000000;
% end
% function z = rhs3(x,y,t)  
%   z=(pi.^2.*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t))./12800 - (sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t))./100 - sin((x.*pi)./16).*sin((y.*pi)./16).*sin(t) - (cos((y.*pi)./32).*sin((x.*pi)./16))./100 + (sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^3.*(2.*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2 - sin((x.*pi)./16).^2 - sin((y.*pi)./16).^2 + 1))./10 + (sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^5.*(2.*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2 - sin((x.*pi)./16).^2 - sin((y.*pi)./16).^2 + 1).^2)./100000000000 + (pi.*cos((y.*pi)./32).*sin((x.*pi)./32).*cos(t).*((pi.*cos((x.*pi)./32).*sin((y.*pi)./32).*cos(t).^2.*(2.*cos((x.*pi)./32).^2 - 1).*(2.*cos((y.*pi)./32).^2 - 1))./320 + (pi.*cos((x.*pi)./32).*cos((y.*pi)./32).^2.*sin((x.*pi)./32).^2.*sin((y.*pi)./32).*cos(t).^2)./80))./3200 + (pi.^2.*cos((y.*pi)./16).^2.*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^3)./256000 - (cos((x.*pi)./16).*cos((y.*pi)./16).*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^2)./5000 - (cos((x.*pi)./32).*cos((y.*pi)./32).*sin((x.*pi)./16).*sin((y.*pi)./16).*cos(t).^2)./5000;
% end
% function z = rhs4(x,y,t)  
%   z=(cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t).^3.*(2.*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2 - sin((x.*pi)./16).^2 - sin((y.*pi)./16).^2 + 1))./10 - (cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t))./100 - cos((x.*pi)./16).*cos((y.*pi)./16).*sin(t) - (cos((x.*pi)./16).^2.*cos((y.*pi)./16).^2.*cos(t).^2)./5000 - cos((y.*pi)./32)./100 + (cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t).^5.*(2.*sin((x.*pi)./16).^2.*sin((y.*pi)./16).^2 - sin((x.*pi)./16).^2 - sin((y.*pi)./16).^2 + 1).^2)./100000000000 + (pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*cos(t))./12800 - (cos((x.*pi)./32).*cos((y.*pi)./32).*cos(t).^2.*(2.*cos((x.*pi)./32).^2 - 1).*(2.*cos((y.*pi)./32).^2 - 1))./5000 + (pi.^2.*cos((x.*pi)./32).^2.*sin((y.*pi)./32).^2.*cos(t).^3.*(2.*cos((y.*pi)./32).^2 - 2.*cos((x.*pi)./32).^2 + 1))./1024000 + (pi.^2.*cos((x.*pi)./16).*cos((y.*pi)./16).*sin((y.*pi)./16).^2.*cos(t).^3)./256000;
% end
function lap=lap_diff(phi,k2)
   lap=real(ifft2((k2.*fft2(phi))));
end

