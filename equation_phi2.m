%% Phase field crystral Lagrange equation
clear 
clc

syms x y t epsilon M1 M2 S1 S2 alpha beta eta gamma12 gamma1 gamma2 eta1 eta2 w0 theta theta1 theta2 B a1 a2 a12   
mu2  = cos(pi*x/32)*cos(pi*y/32)*cos(t);
phi1 = cos(2*pi*x/32)*cos(2*pi*y/32)*cos(t);
phi2 = cos(pi*x/32)*cos(pi*y/32)*cos(t);
M = [sin(8*pi*x/128)*sin(8*pi*y/128)*cos(t),cos(8*pi*x/128)*cos(8*pi*y/128)*cos(t)];
% H = [sin(8*pi*x/128)*cos(8*pi*y/256),cos(8*pi*y/256)];
k1 = diff(phi1,x,2) + diff(phi1,y,2);
k2 = diff(phi2,x,2) + diff(phi2,y,2);

M_1 = sin(8*pi*x/128)*sin(8*pi*y/128)*cos(t);
M_2 = cos(8*pi*x/128)*cos(8*pi*y/128)*cos(t);
% norm1 = M_1.^2 + M_2.^2;
% 
% norm2 = M_1 .* diff(phi2,x) + M_2 .* diff(phi2,y);
% gradient_M = [norm2.*M_1,norm2.*M_2];
% div1 = diff(norm2.*M_1,x)+diff(norm2.*M_2,y);
% 
% gradient_phi = diff(phi2,x).^2 + diff(phi2,y).^2;
% gradient_phiphi = [gradient_phi.* diff(phi2,x),gradient_phi.*diff(phi2,y)];
% div2 = diff(gradient_phi.* diff(phi2,x),x)+diff(gradient_phi.*diff(phi2,y),y);
% 
% R1 = phi2.^3 - epsilon * phi2 + eta.*(abs(phi2)-phi2).*phi2 + gamma12.*phi2.*phi1.^2  - gamma1 * norm1 - theta1 .* div2 + eta1 .* div1;
% mu = (a2.^2).*phi2 +2.*a2.*k2 + ( diff(k2,x,2) + diff(k2,y,2)) + 0.5.*((a12.^2).*phi1 +2.*a12.*k1 + ( diff(k1,x,2) + diff(k1,y,2))) + R1;

g = simplify( diff(phi2,t) + M1.* mu2 - M1 .* 1./(128*128).*int(int(mu2,x,0,128),y,0,128));
% g =  diff(phi2,t) + M1.* mu - M1 .* (1./(128*128).*int(int(mu,x,0,128),y,0,128));

g=char(g);
g= strrep(g,'*','*');
g= strrep(g,'/','./');
g= strrep(g,'^','.^')

