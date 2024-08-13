%% Define decision variables within structure
        variable  P11(nx,nx) symmetric  
        variable  P22(nx*nr,nx*nr) symmetric  
        variable  P33(nx*nr,nx*nr) symmetric  
        variables P12(nx,nx*nr) P13(nx,nx*nr) P23(nx*nr,nx*nr)
%         variable  P121(nx,nx) symmetric  
%         variable  P122(nx,nx) symmetric  
%         variable  P123(nx,nx) symmetric  
%         variable  P131(nx,nx) symmetric         
%         variable  P132(nx,nx) symmetric         
%         variable  P133(nx,nx) symmetric         
        variables K11(nu,nx) K121(nu,nx) K122(nu,nx) K123(nu,nx) K131(nu,nx) K132(nu,nx) K133(nu,nx)
        variables K2211(nu,nx) K2212(nu,nx) K2213(nu,nx) K2222(nu,nx) K2223(nu,nx) K2233(nu,nx) 
        variables K2311(nu,nx) K2312(nu,nx) K2313(nu,nx) K2322(nu,nx) K2323(nu,nx) K2333(nu,nx) 
        variables K3311(nu,nx) K3312(nu,nx) K3313(nu,nx) K3322(nu,nx) K3323(nu,nx) K3333(nu,nx)  
        variables R(nx*nr,nx*nr)      S(nx*nr,nx*nr)            T(nx*nr,nx*nr)  
        variables U(nx,nx)            V(nx,nx)                  W(nx,nx)    
        variables X(nx*nr,nx*nr)      Y(nx*nr,nx*nr)            Z(nx*nr,nx*nr)          
        variable M(3*nr*nx, 5*nr*nx)
 
        P = [P11  P12  P13;
             P12' P22  P23;
             P13' P23' P33];
%         P = [P22  P23;
%              P23' P33];
%         P12 = [P121 P122 P123];
%         P13 = [P131 P132 P133];
        K11 = K11;
        K12 = [K121 K122 K123]; 
        K13 = [K131 K132 K133]; 
        K21 = [K121; K122; K123];
        K22 = [K2211 K2212 K2213;
               K2212 K2222 K2223;
               K2213 K2223 K2233];
        K23 = [K2311 K2312 K2313;
               K2312 K2322 K2323;
               K2313 K2323 K2333];   
        K31 = [K131; K132 ; K133];   
        K32 = [K2311 K2312 K2313;
               K2312 K2322 K2323;
               K2313 K2323 K2333]; 
        K33 = [K3311 K3312 K3313;
               K3312 K3322 K3323;
               K3313 K3323 K3333];    
        K = [K11 K12 K13;
             K21 K22 K23;
             K31 K32 K33];
%         K = [K22 K23;
%              K32 K33];


%% Another expression...
%         variable  sP(3*nx*nr,3*nx*nr) symmetric
%         for i = 1:nr
%             for j = 1:nr
%                 P(:,:,i,j) = kron(e(:,:,i),Inr)'*sP*kron(e(:,:,j),Inr)
%             end            
%         end
 











