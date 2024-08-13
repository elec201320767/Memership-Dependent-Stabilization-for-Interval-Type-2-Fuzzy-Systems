clc
close all
clear

cwd = pwd;
addpath(genpath(cwd));
%% Define System Parameters
nr = 3;                                      % # of rules; Parallel Distirbuted Compensation
Interval =  2;                               % Interval to test Stability Region

x0 = [1; 2];                                     % initial state
step = 0.01;
x = -10:step:10;                                 % defines operating domain

[LowerMF, UpperMF, mf_min, mf_max, mfsum_min, mfsum_max, mfsub_min, mfsub_max] = MF3_numeric(x);
%% Plot Figs for Distributions (for journal)
PlotPolytopes(LowerMF,UpperMF,mfsum_min,mfsum_max);
%%
Vertices = Type2FuzzyVertices([LowerMF; UpperMF]);% Find vertex set
%% Plot Membership Function 
PlotMF(x,LowerMF, UpperMF);
[VerticesDL, VerticesDU, rho] = PlotDMF(x,step,LowerMF, UpperMF); 
%% 3dPlot Vertex Coverage
[VerticesL,VerticesU] = Plot3dVertex(LowerMF,UpperMF,mfsum_min,mfsum_max);
%% Define Size of Vertex Sets
SizeVertices   = size(Vertices);   SizeVerticesL  = size(VerticesL);  SizeVerticesU = size(VerticesU);  
SizeVerticesDL = size(VerticesDL); SizeVerticesDU = size(VerticesDU);
%% Define extrema  
xiL    = rho(1,:)';     xiU = rho(2,:)';     
zetaL  = rho(3,:)';   zetaU = rho(4,:)';   
deltaL = mf_min';    deltaU = mf_max';
sigmaL = mfsum_min;  sigmaU = mfsum_max;
kappaL = mfsub_min'; kappaU = mfsub_max';

%% LMI Construction
B1 = zeros(2,1);
dimsys = size(B1);
nx = dimsys(1);   one_nr = ones(nr,1);
nu = dimsys(2);    

    Inr = eye(nr)   ;  Inx = eye(nx);  Onx = zeros(nx);
    Onenr = ones(nr,1);
    Onrx = zeros(nr*nx); Inrx = eye(nr*nx);
    onelr = ones(nr,1); 
    e1 = [Inx; Onx; Onx];  e2 = [Onx; Inx; Onx];  e3 = [Onx; Onx; Inx];
    for i =1:nr
        e(:,:,i) = Inrx(:,nx*(i-1)+1:(nx*(i-1)+1)+1);
    end
    
    E = zeros(5*nx*nr,nx*nr,5);
    E(:,:,1) = [Inrx; Onrx; Onrx; Onrx; Onrx];
    E(:,:,2) = [Onrx; Inrx; Onrx; Onrx; Onrx]; 
    E(:,:,3) = [Onrx; Onrx; Inrx; Onrx; Onrx]; 
    E(:,:,4) = [Onrx; Onrx; Onrx; Inrx; Onrx]; 
    E(:,:,5) = [Onrx; Onrx; Onrx; Onrx; Inrx]; 
    
%     e1 = [Inx; Onrx];  e2 = [Onrx'; Inrx];   
    
% Set of Feasibility region
FeasRegion = zeros(1,2);            FeasRegionInf = zeros(1,2);   
FeasRegionInacUnb  = zeros(1,2);    FeasRegionInacSol  = zeros(1,2);
FeasRegionUnb = zeros(1,2);         FeasRegionInacInf  = zeros(1,2); 
FeasRegionFail = zeros(1,2);        FeasRegionOver = zeros(1,2); 
FeasRegion(1,:) = [];           FeasRegionInacInf(1,:) = [];
FeasRegionInacSol(1,:) = [];    FeasRegionInacUnb(1,:) = [];
FeasRegionUnb(1,:) = [];        FeasRegionInf(1,:) = [];
FeasRegionFail(1,:) = [];       FeasRegionOver(1,:) = [];
close all;
% Define LMI
a_min = 14; a_max =  14;
b_min = 26; b_max = 26;
for avar = a_min:Interval:a_max
    for bvar = b_min:Interval:b_max
        %Example: 2008 Stability Analysis of Interval Type-2 Fuzzy Model-Based Control Systems
        [A, B] = SubSys3_numeric(avar,bvar);    
        A_Col = [A(:,:,1); A(:,:,2); A(:,:,3)];
        B_Col = [B(:,:,1); B(:,:,2); B(:,:,3)];
        
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    P = [P11  P12  P13;  K = [K11 K12 K13; 
%         P12' P22  P23;       K21 K22 K23;
%         P13' P23' P33]       K31 K32 K33] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cvx_clear
 tic       
        cvx_begin sdp
        cvx_precision best                                                                                      %       # of var    cumulat
        Cor1_deflmivar_cvx
       
        
        P >= eps*eye(nx+2*nr*nx); 
        for i = 1:nr
            for j = 1:nr
                e(:,:,i)'*R*e(:,:,j) + (e(:,:,i)'*R*e(:,:,j))' >= eps*Inx;   % R_ij + R_ij' > 0
                e(:,:,i)'*S*e(:,:,j) + (e(:,:,i)'*S*e(:,:,j))' >= eps*Inx;   % S_ij + S_ij' > 0
                e(:,:,i)'*T*e(:,:,j) + (e(:,:,i)'*T*e(:,:,j))' >= eps*Inx;   % T_ij + T_ij' > 0
                e(:,:,i)'*X*e(:,:,j) + (e(:,:,i)'*X*e(:,:,j))' >= eps*Inx;   % X_ij + X_ij' > 0
                e(:,:,i)'*Y*e(:,:,j) + (e(:,:,i)'*Y*e(:,:,j))' >= eps*Inx;   % Y_ij + Y_ij' > 0
                e(:,:,i)'*Z*e(:,:,j) + (e(:,:,i)'*Z*e(:,:,j))' >= eps*Inx;   % Z_ij + Z_ij' > 0
            end
        end
               U + U' >= eps*Inx;   % U + U' > 0
               V + V' >= eps*Inx;   % V + V' > 0
               W + W' >= eps*Inx;   % W + W' > 0
               
         
               
        Rlx_DMF = zeros(5*nx*nr);
         Rlx_MF = zeros(5*nx*nr);
                    %%% Relaxation by Extrema : Diff MF
                    for i = 1:nr %%%%%%%%%% (1,1) %%%%%%%%%% % By shrink Lemma
                        for j = 1:nr
                            for k = 1:nr
                                for l = 1:nr
                                    Rlx_DMF = Rlx_DMF + E(:,:,1)*e(:,:,i)*(-1)*xiL(k)*xiU(l)*e(:,:,k)'*Y*e(:,:,l)*e(:,:,j)'*E(:,:,1)';
                                    Rlx_DMF = Rlx_DMF + E(:,:,1)*e(:,:,i)*(-1)*zetaL(k)*zetaU(l)*e(:,:,k)'*Z*e(:,:,l)*e(:,:,j)'*E(:,:,1)';
                                end
                            end
                        end
                    end 
                      
                    for i = 1:nr %%%%%%%%%% (1,4),(1,5) %%%%%%%%%% % By shrink Lemma
                        Rlx_DMF = Rlx_DMF + E(:,:,1)*e(:,:,i)*kron(xiL,Inx)'*Y*E(:,:,4)';
                        Rlx_DMF = Rlx_DMF + E(:,:,1)*e(:,:,i)*kron(xiU,Inx)'*Y'*E(:,:,4)';
                        Rlx_DMF = Rlx_DMF + E(:,:,1)*e(:,:,i)*kron(zetaL,Inx)'*Z*E(:,:,5)';
                        Rlx_DMF = Rlx_DMF + E(:,:,1)*e(:,:,i)*kron(zetaU,Inx)'*Z'*E(:,:,5)';
                    end
                    %%%%%%%%%% (4,4),(5,5) %%%%%%%%%%
                    Rlx_DMF = Rlx_DMF + E(:,:,4)*(-1)*Y*E(:,:,4)';
                    Rlx_DMF = Rlx_DMF + E(:,:,5)*(-1)*Z*E(:,:,5)';
                    
                    %%% Relaxation by Extrema : MF (1), mathcal{M}_1
                    for i = 1:nr %%%%%%%%%% (1,1)~(1,3) %%%%%%%%%% % By shrink Lemma
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*kron(deltaL,Inx)'*T'*E(:,:,1)';
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*kron(deltaU,Inx)'*R'*E(:,:,1)';
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*kron(deltaL,Inx)'*T'*E(:,:,2)';
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*kron(deltaU,Inx)'*R'*E(:,:,3)';
                    end
                    %%%%%%%%%% (1,1)~(1,3),(2,2)~(2,3),(3,3) %%%%%%%%%%
                    Rlx_MF = Rlx_MF + E(:,:,1)*(-1)*S*E(:,:,1)';
                    Rlx_MF = Rlx_MF + E(:,:,1)*S*E(:,:,2)';
                    Rlx_MF = Rlx_MF + E(:,:,1)*T*E(:,:,2)';
                    Rlx_MF = Rlx_MF + E(:,:,1)*R*E(:,:,3)';
                    Rlx_MF = Rlx_MF + E(:,:,1)*S'*E(:,:,3)';
                    Rlx_MF = Rlx_MF + E(:,:,2)*(-1)*T*E(:,:,2)';
                    Rlx_MF = Rlx_MF + E(:,:,2)*(-1)*S'*E(:,:,3)';
                    Rlx_MF = Rlx_MF + E(:,:,3)*(-1)*R*E(:,:,3)';
                    
                    %%% Relaxation by Extrema : MF (2), mathcal{M}_2
                    for i = 1:nr %%%%%%%%%% (1,1) %%%%%%%%%% % By shrink Lemma
                        for j = 1:nr
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*sigmaL*W*e(:,:,j)'*E(:,:,1)';
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*V*e(:,:,j)'*E(:,:,1)';
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*sigmaU*U*e(:,:,j)'*E(:,:,1)';
                        end
                    end
                    for i = 1:nr %%%%%%%%%% (1,2)~(1,3) %%%%%%%%%% % By shrink Lemma
                        for j =1:nr
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(1+sigmaL)*W*e(:,:,j)'*E(:,:,2)';
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*V*e(:,:,j)'*E(:,:,2)';
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(1+sigmaU)*U*e(:,:,j)'*E(:,:,3)';
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*V*e(:,:,j)'*E(:,:,3)';
                        end
                    end
                    %%%%%%%%%% (2,2),(2,3),(3,3) %%%%%%%%%%
                    for i = 1:nr
                        for j = 1:nr
                            Rlx_MF = Rlx_MF + E(:,:,2)*(-1)*e(:,:,i)*W*e(:,:,j)'*E(:,:,2)';
                            Rlx_MF = Rlx_MF + E(:,:,2)*(-1)*e(:,:,i)*V*e(:,:,j)'*E(:,:,3)';
                            Rlx_MF = Rlx_MF + E(:,:,3)*(-1)*e(:,:,i)*U*e(:,:,j)'*E(:,:,3)';
                        end
                    end
                   
                    %%% Relaxation by Extrema : MF (3), mathcal{M}_3
                    for i = 1:nr %%%%%%%%%% (1,1) %%%%%%%%%% % By shrink Lemma
                        for j = 1:nr
                            Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*kron(kappaL,Inx)'*X*kron(kappaU,Inx)*e(:,:,j)'*E(:,:,1)';
                        end
                    end
                    for i = 1:nr %%%%%%%%%% (1,2),(1,3) %%%%%%%%%% % By shrink Lemma
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*kron(kappaU,Inx)'*X'*E(:,:,2)';
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*(-1)*kron(kappaL,Inx)'*X*E(:,:,2)';
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*kron(kappaU,Inx)'*X'*E(:,:,3)';
                        Rlx_MF = Rlx_MF + E(:,:,1)*e(:,:,i)*kron(kappaL,Inx)'*X*E(:,:,3)';
                    end
                    %%%%%%%%%% (2,2),(2,3),(3,3) %%%%%%%%%%
                    Rlx_MF = Rlx_MF + E(:,:,2)*(-1)*X*E(:,:,2)';
                    Rlx_MF = Rlx_MF + E(:,:,2)*X*E(:,:,3)';
                    Rlx_MF = Rlx_MF + E(:,:,2)*X'*E(:,:,3)';
                    Rlx_MF = Rlx_MF + E(:,:,3)*(-1)*X*E(:,:,3)';
                    
        % System within LMIs
        for h = 1:SizeVertices(2) % SizeVertices = [#Rules, #Vertices]
            Ah = zeros(nx); Bh = zeros(nx,nu);
            for iter = 1:nr
                Ah = Ah + Vertices(iter,h)*A(:,:,iter);
                Bh = Bh + Vertices(iter,h)*B(:,:,iter);
            end
                     Sys = zeros(5*nx*nr);  
                    for i = 1:nr %%%%%%%%%% (1,1) %%%%%%%%%% % By shrink Lemma
                        for j = 1:nr
                            Sys = Sys + E(:,:,1)*e(:,:,i)*Ah*P11*e(:,:,j)'*E(:,:,1)';
                            Sys = Sys + E(:,:,1)*e(:,:,i)*Bh*K11*e(:,:,j)'*E(:,:,1)';
                        end
                    end
                    for i = 1:nr %%%%%%%%%% (1,2) ~ (1,5) %%%%%%%%%% % By shrink Lemma
                        Sys = Sys + E(:,:,1)*e(:,:,i)*Ah*P12*E(:,:,2)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*P12*kron(Inr,Ah)'*E(:,:,2)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*Bh*K12*E(:,:,2)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*K21'*kron(Inr,Bh)'*E(:,:,2)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*Ah*P13*E(:,:,3)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*P13*kron(Inr,Ah)'*E(:,:,3)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*Bh*K13*E(:,:,3)';
                        Sys = Sys + E(:,:,1)*e(:,:,i)*K31'*kron(Inr,Bh)'*E(:,:,3)';
                        Sys = Sys +-E(:,:,1)*e(:,:,i)*P12*E(:,:,4)';
                        Sys = Sys +-E(:,:,1)*e(:,:,i)*P13*E(:,:,5)';
                    end
                    %%%%%%%%%% (2,2) ~ (2,5) %%%%%%%%%%
                    Sys = Sys + E(:,:,2)*kron(Inr,Ah)*P22*E(:,:,2)';
                    Sys = Sys + E(:,:,2)*kron(Inr,Bh)*K22*E(:,:,2)';
                    Sys = Sys + E(:,:,2)*kron(Inr,Ah)*P23*E(:,:,3)';
                    Sys = Sys + E(:,:,2)*P23*kron(Inr,Ah)'*E(:,:,3)';
                    Sys = Sys + E(:,:,2)*kron(Inr,Bh)*K23*E(:,:,3)';
                    Sys = Sys + E(:,:,2)*K32'*kron(Inr,Bh)'*E(:,:,3)';
                    Sys = Sys +-E(:,:,2)*P22*E(:,:,4)';
                    Sys = Sys +-E(:,:,2)*P23*E(:,:,5)';
                    %%%%%%%%%% (3,3) ~ (3,5) %%%%%%%%%%
                    Sys = Sys + E(:,:,3)*kron(Inr,Ah)*P33*E(:,:,3)';
                    Sys = Sys + E(:,:,3)*kron(Inr,Bh)*K33*E(:,:,3)';
                    Sys = Sys +-E(:,:,3)*P23'*E(:,:,4)';
                    Sys = Sys +-E(:,:,3)*P33*E(:,:,5)';
                    
            for hL = 1:SizeVerticesL(2)
                for hU = 1:SizeVerticesU(2)
                    
                    %%% Orthogonal Complement
                    Orth1 = kron(kron(Onenr,Vertices(:,h)'),Inx) ;
                    Orth2 = kron(kron(Onenr,VerticesL(:,hL)'),Inx) ;
                    Orth3 = kron(kron(Onenr,VerticesU(:,hU)'),Inx) ;
                    
                    Orth = [Orth1-Inrx Orth2  Orth3;
                            Onrx       -Inrx  Onrx;
                            Onrx       Onrx   -Inrx;
                            Onrx       Onrx  Onrx;
                            Onrx       Onrx  Onrx];
                    Zero = zeros(5*nx*nr);
                    Zero = Orth*M;
                    
                    Cond = Sys + Zero + Rlx_MF + Rlx_DMF;
%                     Cond = Sys + Rlx_MF + Rlx_DMF;
                    
                    Cond + Cond' <= -eps*eye(5*nx*nr);
                end
            end
        end
disp('Time(s) for setting LMIs:');        
 toc         
 tic
        cvx_end
disp('Time(s) for solving LMIs:');         
 toc
        switch cvx_status
            case 'Solved'
                FeasRegion = [FeasRegion;avar bvar];
            case 'Unbounded'
               FeasRegionUnb = [FeasRegionUnb;avar bvar];
            case 'Infeasible'
                FeasRegionInf = [FeasRegionInf;avar bvar];  
            case 'Failed'
                FeasRegionFail = [FeasRegionFail;avar bvar];  
            case 'Overdetermined'
                FeasRegionOver = [FeasRegionOver;avar bvar];     
            case 'Inaccurate/Solved' 
                FeasRegionInacSol = [FeasRegionInacSol;avar bvar];
            case 'Inaccurate/Unbounded' 
                FeasRegionInacUnb = [FeasRegionInacUnb;avar bvar];
            case 'Inaccurate/Infeasible' 
                FeasRegionInacInf = [FeasRegionInacInf;avar bvar];    
        end
        avar
        bvar                
    end
end
 
    Region2008Lam = [23 23;24 23;25 24];
    Region2020Wang = [17 23;17 24;17 25;17 26;17 27;
                      18 23;18 24;18 25;18 26;18 27;18 28;
                      19 23;19 24;19 25;19 26;19 27;19 28;19 29;
                      20 23;20 24;20 25;20 26;20 27;20 28;20 29;20 30;
                      21 23;21 24;21 25;21 26;21 27;21 28;21 29;21 30;
                      22 23;22 24;22 25;22 26;22 27;22 28;22 29;22 30;22 31];
    

%%
% NumFeas = size(FeasRegion);
% % Plot the Feasible Region
% figure()
% for i = 1:NumFeas(1)
% % plot(FeasRegion(i,1),FeasRegion(i,2),'MarkerColor',[0 0 1],'o')
% plot(FeasRegion(i,1),FeasRegion(i,2),'o','MarkerSize',10,'MarkerEdgeColor','b')
% hold on
% end
% title('Feasibility region')
%%
  mdl = 'Cor1_Ex1';
  open_system(mdl)
  simOut = sim(mdl);
           state = simOut.StateResponse;
           input = simOut.InputResponse;
           
           figure()
           plot(state.Time(:,1),state.Data(:,1),'k'); hold on;
           plot(state.Time(:,1),state.Data(:,2),'b'); hold on;
           grid on
           ylim([-5 3])
           xlabel('Time(s)')
            ylabel('States')
           
            figure()
            plot(input.Time(:,1),input.Data(:,1),'k');
            grid on
            xlabel('Time(s)')
            ylabel('Input')
%            ylim([-5 3])

figure()
subplot(2,1,1);
plot(state.Time(:,1),state.Data(:,1),'k'); hold on;
plot(state.Time(:,1),state.Data(:,2),'b'); hold on;
grid on
ylim([-5 3])
xlabel('Time(s)')
ylabel('States')
legend('$x_1(t)$','$$x_2(t)$$','Interpreter','Latex')

subplot(2,1,2);
plot(input.Time(:,1),input.Data(:,1),'k');
grid on
xlabel('Time(s)')
ylabel('Input')
legend('$u(t)$','Interpreter','Latex')