function z = init_growth1_cutoff2(x,y) 
load('(0,0)gamma12=0.2eta=0.3theta=-45.mat');
%         phi_0 = 0.14;
        d = 10;
%         L = 512;
        
        x0 = 256; y0=256;
        w1 = ((x-x0).^2+(y-y0).^2<=d.^2);
%         w1  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
%         w1 = w1.*((x-x0).^2+(y-y0).^2<=d^2);
        
%         x0 = 306; y0=256;
%         w2 = ((x-x0).^2+(y-y0).^2<=d.^2);
%         w2  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
%         w2 = w2.*((x-x0).^2+(y-y0).^2<=d^2);
        
%         x0 = 600; y0=300;
%         w3 = (abs(x-x0)<d&abs(y-y0)<=d);
%         w3  = (1-((x-x0).^2+(y-y0).^2)/d^2).^2;
%         w3 = w3.*((x-x0).^2+(y-y0).^2<=d^2);
%         theta = -;
%         x2=x;y2=y;
%         x1=x;y1=y;
%         theta = pi/3;
%         x1 =  sin(theta)*x + cos(theta)*y;
%         y1 = -cos(theta)*x + sin(theta)*y;
%         
%         theta = -pi/4;
%         x2 =  cos(theta)*x + sin(theta)*y;
%         y2 = -sin(theta)*x + cos(theta)*y;
        

        
%         eta = 0;
%         s   = 0;
%         
%         kx = sqrt(3/4*(1+eta^2*s^2)/(1+(3/2)*eta^2*s^2));
%         ky = sqrt(3/4*(1+eta)^2/(1+(3/2)*eta^2*s^2));
        phi_bar = 0.285;       
%         R = (1/2)*eta^2*s^2/(1+(3/2)*eta^2*s^2);        
%         A_cof = -4/5*(phi_bar+sqrt(15*(epsilon+R)-36*phi_bar.^2)/3);        
%         z1 = cos(kx*x).*cos(ky*(y-s*x)/sqrt(3))+1/2*cos(2*ky*(y-s*x)/sqrt(3));
%         z2 = cos(kx*x2).*cos(ky*(y2-s*x2)/sqrt(3))+1/2*cos(2*ky*(y2-s*x2)/sqrt(3));
%         z3 = cos(kx*x3).*cos(ky*(y3-s*x3)/sqrt(3))+1/2*cos(2*ky*(y3-s*x3)/sqrt(3));  
%         z1 = C*(cos(q*y/sqrt(3)).*cos(q*x) - 0.5*cos(2*q*y/sqrt(3)) );
%         z2 = C*(cos(q*y2/sqrt(3)).*cos(q*x2) - 0.5*cos(2*q*y2/sqrt(3)) );
%         z3 = C*(cos(q*y3/sqrt(3)).*cos(q*x3) - 0.5*cos(2*q*y3/sqrt(3)) );
 z1 = phi2;
%  z2=0.446*(cos(0.66/sqrt(3)*y2).*cos(0.66*x2)-0.5*cos(2*0.66/sqrt(3)*y2));
%  z3=0.446*(cos(0.66/sqrt(3)*y3).*cos(0.66*x3)-0.5*cos(2*0.66/sqrt(3)*y3));
%  z=0.285+z1+z2+z3;
%         z=phi_bar + w1.*z1 + w2.*z2+ w3.*z3;
          z=(1-w1).*phi_bar + w1.*z1 ;
%         z=phi_0 + w1.*(A_cof*z1) + w2.*(A_cof*z2)+ w3.*(A_cof*z3);
%         save('ex17_1_MPFCdata_Init_u0.mat','z');
%         load('ex17_1_MPFCdata_Init_u0.mat');
    end