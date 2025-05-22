clear 
clc

syms x y  t  x y t epsilon M1 M2 S1 S2 alpha beta eta gamma12 gamma1 gamma2 eta1 eta2 w0 theta theta1 theta2 B a1 a2 a12

phi1 = cos(2*pi*x/32)*cos(2*pi*y/32)*cos(t);
phi2 = cos(pi*x/32)*cos(pi*y/32)*cos(t);
M = [sin(8*pi*x/128)*sin(8*pi*y/128)*cos(t),cos(8*pi*x/128)*cos(8*pi*y/128)*cos(t)];
H = [sin(8*pi*x/128)*cos(8*pi*y/256),cos(8*pi*y/256)];
H1 = sin(8*pi*x/128)*cos(8*pi*y/256);
H2 = cos(8*pi*y/256);
k1 = diff(phi1,x,2) + diff(phi1,y,2);
k2 = diff(phi2,x,2) + diff(phi2,y,2);

M_1 = sin(8*pi*x/128)*sin(8*pi*y/128)*cos(t);
M_2 = cos(8*pi*x/128)*cos(8*pi*y/128)*cos(t);
norm1 = M_1.^2 + M_2.^2;
k = [diff(M_1,x,2)+diff(M_1,y,2) ,diff(M_2,x,2) + diff(M_2,y,2)];

norm2_1 = M_1.*diff(phi1,x) + M_2.*diff(phi1,y);
norm2_2 = M_1.*diff(phi2,x) + M_2.*diff(phi2,y);
R2_1 = -alpha .* M_1 +beta .* (M_1.^2+M_2.^2).*M_1 - H1 - 2.* gamma1 * M_1 * phi1 - 2.* gamma2 * M_1 * phi2 + theta * (M_1.^2+M_2.^2).^2 .* M_1 - eta1 * norm2_1 .* diff(phi1,x) - eta2 * norm2_2 .* diff(phi2,x);

R2_2 = -alpha .* M_2 +beta .* (M_1.^2+M_2.^2).*M_2 - H2 - 2 * gamma1 * M_2 * phi1 - 2.* gamma2 * M_2 * phi2 + theta * (M_1.^2+M_2.^2).^2 .* M_2 - eta1 * norm2_1 .* diff(phi1,y) - eta2 * norm2_2 .* diff(phi2,y);

R2 = [R2_1,R2_2];
mu = -w0*k + R2;
M_t=[diff(M_1,t),diff(M_2,t)];
% tmp = phi.^3 +(1-epsilon)*phi;

g = simplify(M_t + M2* mu );
% n = g-[-(sin((x.*pi)./16).*sin((y.*pi)./16).*(128.*sin(t) + 128.*M2.*alpha.*cos(t) - M2.*w0.*pi.^2.*cos(t)))./128, -(cos((x.*pi)./16).*cos((y.*pi)./16).*(128.*sin(t) + 128.*M2.*alpha.*cos(t) - M2.*w0.*pi.^2.*cos(t)))./128];
% simplify(n)
g=char(g);
g= strrep(g,'*','.*');
g= strrep(g,'/','./');
g= strrep(g,'^','.^')
