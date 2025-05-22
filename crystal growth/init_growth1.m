function z = init_growth1(x,y,theta) 
        d = 10;
%         L = 512;
%         x0 = 128; y0=128;
        x0 = 256; y0=256;
        w1 = ((x-x0).^2+(y-y0).^2<=d.^2);

        x1 =  sin(theta)*x + cos(theta)*y;
        y1 = -cos(theta)*x + sin(theta)*y;
        
%         theta = -pi/4;-
%         x1 =  cos(theta)*x + sin(theta)*y;
%         y1 = -sin(theta)*x + cos(theta)*y;
          phi_bar = 0.285;       

          z1 = 0.446*(cos(0.66/sqrt(3)*y1).*cos(0.66*x1)-0.5*cos(2*0.66/sqrt(3)*y1));

          z = phi_bar + w1.*z1 ;
    end