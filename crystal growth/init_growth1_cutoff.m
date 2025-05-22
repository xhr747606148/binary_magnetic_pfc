function z = init_growth1_cutoff(x,y) 
load('(0,0)gamma12=0.2eta=0.3theta=-45.mat');
%         phi_0 = 0.14;
        d = 10;
%         L = 512;
        
        x0 = 256; y0=256;
        w1 = ((x-x0).^2+(y-y0).^2<=d.^2);

        phi_bar = 0.285;       

         z1 = phi1;

          z=(1-w1).*phi_bar + w1.*z1 ;

    end