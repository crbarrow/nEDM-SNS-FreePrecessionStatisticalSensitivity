function CBarrow_FrPrecStatSens_Diff_f_MLmle01
%--- Calculating statistical impact of SNS nEDM experiment for f-factor and f-factor w peak patch loss model instead of single tau
clc; set(0,'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultLegendInterpreter', 'latex'); set(gcf,'units','centimeters'); linecolors = lines(20); set(0,'DefaultAxesFontSize', 14); set(0,'DefaultTextFontSize', 14);
%warning('on','all');
warning('off','all');
global plotsOn; plotsOn = 0;
     
%---- Fitting Method Settings -------------
expFit_s = 1;    % exponential fitting parameter -- 1==single ; 2==double

%---- Generating Variables ----------------
f_walls_s = 2;   % f_walls parameter -- 1==1.0e-5 ; 2==2.0e-5

%---- Fitting Variable Loops ---- !!!! Values Set at Loops !!!! ----
phi0Loops = 1;    % phi0 for loops on/off 1/0 --       % phi0 parameter -- 0==free ; 1==fixed value
T_wallsLoops = 0; % tau_walls for loops on/off 1/0 --  % Tau_walls parameter -- 0==(E)dependent ; 1==single valued
patchLoops = 1;   % weak patch for loops on/off 1/0 -- % weak_patch parameter -- 0==no patch ; 1==patch considered

%--- operating parameters 
tau_n3_scan = 500; %average for unpolarized n
Tm_scan = 1000; %[s] measurement time
Tf_scan = 1000; %[s] cold neutron fill time

%--- fixed parameters
global Emax t_width linecolors
VCell = 40.0*10.1*7.6; %[cm^3]
AWalls = 10.1*7.6*2+40.0*7.6*2+40.0*10.1*2; %[cm^2]
VoptWalls = 160; %[neV]
%VoptWalls = 200; %[neV]
fWalls = (f_walls_s)*1.0e-5;

fpatch = 5e-5; %[unitless]
Upatch = 80; %[neV]
Apatch = 0.2; %[cm^2]
BG = 5; %[s^-1] other background rates

Emax = VoptWalls; %Max energy of stored UCNs
vmax = sqrt(2*Emax*10^-9*1.6e-19/1.67e-27);

dE = 1; %[neV]
E = [dE/2:dE:Emax];
v = sqrt(2*E*1E-9*1.6e-19/1.67e-27); %[m/s]
n_v = v.^2;
n_E = sqrt(E);
n_E = n_E/max(n_E);
 
mubarWalls = 2*fWalls*(VoptWalls./E.*asin(sqrt(E/VoptWalls))-sqrt(VoptWalls./E-1));
mubarWalls(E/VoptWalls>1)=1;
if f_walls_s==1
    tau_walls_scalar = 2000; %%%%---- SEE LOOPS FOR VARIABLE VALUE ----
    tau0 = 2000; % bet0 value of tau_walls before fitted
else
    tau_walls_scalar = 1000; %%%%---- SEE LOOPS FOR VARIABLE VALUE ----
    tau0 = 1000; % bet0 value of tau_walls before fitted
end

mubarPatch = 2*fpatch*(Upatch./E.*asin(sqrt(E/Upatch))-sqrt(Upatch./E-1));
mubarPatch(E/Upatch>1)=1;
tau_patch = 4*VCell*1e-6./(mubarPatch.*v*Apatch*1e-4);

P_UCN = 0.31*(VoptWalls/160)^(3/2); %[UCN/cc/s]
%P_UCN = 0.26*(VoptWalls/160)^(3/2); %[UCN/cc/s] %used in Brad F's calculations
tau_beta = 885;
tau_buildup = (1/tau_beta+1./tau_walls_scalar).^-1; %when accumulating, UCNs are aligned with 3He. ---- note tau_walls dependence
PE = 3/(2*Emax^(3/2))*P_UCN*sqrt(E)*VCell; % P(E) [UCN/s/neV]

Efield = 74 * 1E3 * 100; %[V/m] = [kV/cm] * 1E3 [V/kV] * 100 [cm/m]
Td = 400; %[s] dead time between cycles
P3 = 0.98; % initial 3He pol
Pn = 0.98; %initial UCN pol
Gammap = 1/20000; %[s^-1] 3He and UCN depolarization rate
eps3 = 0.93; % n-3He absorption detection efficency
epsbeta = 0.5; %beta decay detection efficiency

B0 = 30E-3* 1E-4; %[T] = [G]* 1E-4 [T/G]
gamma_n = 1.83247172e8; %s^-1*T^-1 or 18.32472 [Hz/mG] * 1E3 [mG/G] * 1E4 [G/T]
gamma_3 = 2.037894659e8; %s^-1*T^-1  [CODATA]
omega3n = (gamma_3-gamma_n)*B0;
f3n = omega3n/(2*pi);
phi0 = 30*pi/180; % phase angle of UCN precession 

t_width = 20E-3; %[s]

% ---- 'for loops' settings and values ----
%%%% NUMLOOPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLoops = 99;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
loopCount = (phi0Loops+1)*(T_wallsLoops+1)*(patchLoops+1);
doneRes = 100/(numLoops*loopCount);
phi0Range = (1:phi0Loops+1);
TwallsRange = (1:T_wallsLoops+1);
WP_Range = (1:patchLoops+1);

D_f_LS = zeros([1 numLoops]);D_f_RLS = zeros([1 numLoops]);D_f_MLE = zeros([1 numLoops]);
D_phi_LS = zeros([1 numLoops]);D_phi_RLS = zeros([1 numLoops]);D_phi_MLE = zeros([1 numLoops]);
D_tau1_LS = zeros([1 numLoops]);D_tau1_RLS = zeros([1 numLoops]);D_tau1_MLE = zeros([1 numLoops]);
D_tau2_LS = zeros([1 numLoops]);D_tau2_RLS = zeros([1 numLoops]);D_tau2_MLE = zeros([1 numLoops]);
S_f_LS = zeros([1 numLoops]);S_f_RLS = zeros([1 numLoops]);S_f_MLE = zeros([1 numLoops]);
LStime = zeros([1 numLoops]);RLStime = zeros([1 numLoops]);MLEtime = zeros([1 numLoops]);

tic
for i = phi0Range
    phi0_s = (i-1);
    iProg = (patchLoops+1)*(T_wallsLoops+1)*100/loopCount*(i-1);   
    for ii = WP_Range
        patch_s = (ii-1);
        iiProg = (T_wallsLoops+1)*100/loopCount*(ii-1);
        for j = TwallsRange
            tau_walls_s = (j-1);
            jProg = 100/loopCount*(j-1);
            if patch_s==0
                if tau_walls_s==0
                    tau_walls = 4*VCell*1e-6./(mubarWalls.*v*AWalls*1e-4);
                    tau_buildup = (1/tau_beta+1./tau_walls).^-1;
                else
                    tau_walls = tau_walls_scalar*ones(size(v)); 
                    tau_buildup = (1/tau_beta+1./tau_walls).^-1;
                end
            else
                if tau_walls_s==0
                    tau_walls = 4*VCell*1e-6./(mubarWalls.*v*AWalls*1e-4);
                    tau_buildup = (1/tau_beta+1./tau_walls+1./tau_patch).^-1;
                else
                    tau_walls = tau_walls_scalar*ones(size(v)); 
                    tau_buildup = (1/tau_beta+1./tau_walls+1./tau_patch).^-1;
                end
            end
            
 %%%%%%%%%%   Random (Poisson) Generation    %%%%%%%%%%%%           
            for k = 1:numLoops
            
                tau_n3 = tau_n3_scan;
                Tm = Tm_scan;
                Tf = Tf_scan;
                
                N0E = PE.*tau_buildup.*(1-exp(-Tf./tau_buildup)); % N0(E)[UCN/neV]
                tau_tot = (1/tau_beta+1/tau_n3+1./tau_walls).^-1; %average loss rate of UCNs
                Gammatot = 1./tau_tot;
                
                n = round(Tm/t_width);
                t = linspace(0,Tm,n);
                
                N0 = sum(N0E*dE);
                clear NEt;
                NEt = zeros([length(N0E),length(t)]);
                for l = 1:length(N0E)
                    for m = 1:length(t)
                        NEt(l,m) = N0E(l).*exp(-Gammatot(l)*t(m)); %Z(E,t)
                    end
                end
        %%%%%%%%%   Random (Poisson) Light Generation   %%%%%%%%
                Nt = sum(NEt,1)*dE;
                y = (Nt.*(epsbeta/tau_beta + eps3/tau_n3.*(1-P3*Pn*cos(2*pi*f3n*t+phi0)))+BG)*t_width;
                
                data = random('Poisson',y);
                dataErr = sqrt(data); %first guess at error assuming gaussian
                dataErr(dataErr<=0)=1; %have to save the case when error = 0
                
                timeFit = toc;
                
%SS%%%%%%%%%%%%%     SINGLE Exponential Fitting     %%%%%%%%%%%%%%%%%%%%%%
                if expFit_s==1
                    %---- single exponential fit
                    if phi0_s==1
                        FIDsignal = @(p,t) p(1).*exp(-t./p(2)).*(1-p(3)*cos(2*pi*p(4)*t + phi0))+abs(p(5));
                        beta0 = [Nt(1)*t_width*((epsbeta/tau_beta)+(eps3/tau_n3)) 1/(1/tau_n3+1/880+1/tau0) P3*Pn*(eps3/tau_n3)/((epsbeta/tau_beta)+(eps3/tau_n3)) f3n BG*t_width];
                    else
                        FIDsignal = @(p,t) p(1).*exp(-t./p(2)).*(1-p(3)*cos(2*pi*p(4)*t + p(6)))+abs(p(5));
                        beta0 = [Nt(1)*t_width*((epsbeta/tau_beta)+(eps3/tau_n3)) 1/(1/tau_n3+1/880+1/tau0) P3*Pn*(eps3/tau_n3)/((epsbeta/tau_beta)+(eps3/tau_n3)) f3n BG*t_width phi0];
                    end
                    
                %%%%%---- Least Squares fitting to single exponential --%%
            % fitnlm without "Robust" option 
            % fits Least Squares "linear"
                    lsmdl = fitnlm(t,data,FIDsignal,beta0);
            % delta_frequency after Least Squares fit 
            % f3n - 4th fitted coef.
                    D_f_LS(k) = f3n - lsmdl.Coefficients.Estimate(4);
            % sigma (standard error) for freq from LS fitting
                    S_f_LS(k) = lsmdl.Coefficients.SE(4);
            % delta_tau_walls after LS fit 
                    D_tau1_LS(k) = tau0 - 1/(1/lsmdl.Coefficients.Estimate(2)-1/tau_beta-1/tau_n3);
            % delta_phi0 after LS fit
                    if phi0_s==0
                        D_phi_LS(k) = phi0 - lsmdl.Coefficients.Estimate(6);
                    else
                        D_phi_LS(k) = 0;
                    end
                    LStime(k) = toc - timeFit;
                                                          
                    %---- Recursive Least Squares fitting to single exponential
            % function defined below
                    [betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr);
            % saved values from RLS fitting
                    S_f_RLS(k) = errBest2(4);
                    D_f_RLS(k) = f3n - betaBest2(4);
                    D_tau1_RLS(k) = tau0 - 1/(1/betaBest2(2)-1/tau_beta-1/tau_n3);
                    if phi0_s==0
                        D_phi_RLS(k) = phi0 - betaBest2(6);
                    else
                        D_phi_RLS(k) = 0;
                    end
                    
                    RLStime(k) = toc - (timeFit+LStime(k));
                    
                    %---- log likelihood (Poisson Nonlinear Regression)
            % negative log likelihood for FIDsignal
                    mynegloglik = @(beta,data1,cens,freq) -sum(log(poisspdf(data,FIDsignal(beta,t))));
            % mle fitting defined by Matlab mle
                    [phatMLE,pciMLE] = mle(data,'nloglf',mynegloglik,'start',beta0);
                    phatMLE
                    pciMLE
            % extracted values from Matlab mle
                    D_f_MLE(k) = f3n - phatMLE(4);
                    S_f_MLE(k) = phatMLE(4) - pciMLE(1,4);
                    D_tau1_MLE(k) = tau0 - 1/(1/phatMLE(2)-1/tau_beta-1/tau_n3);
                    if phi0_s==0
                        D_phi_MLE(k) = phi0 - phatMLE(6);
                    else
                        D_phi_MLE(k) = 0;
                    end
                    MLEtime(k) = toc - (timeFit+LStime(k)+RLStime(k));
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    mynegloglik = @(beta) -sum(log(poisspdf(data,FIDsignal(beta,t))));
%                    opts = optimset('fminsearch');
%                    opts.MaxFunEvals = Inf;
%                    opts.MaxIter = 10000;
%                    betaHatML = fminsearch(mynegloglik,beta0,opts);
%                    %disp('Poisson Nonlinear regression:')
%                    %disp(['Coeff  =  ' num2str(betaHatML,'%0.8e  ')])
%                    D_f_MLE(k) = f3n - betaHatML(4);
%                    MLEresiduals = FIDsignal(betaHatML,t)-data;
%                    [muMLE,SigmaMLE] = normfit(MLEresiduals);
%                    S_f_MLE(k) = SigmaMLE;   %%%%%%%%%%%% sigma MLE %%%%%%%%%
%                    D_tau1_MLE(k) = 2000 - 1/(1/betaHatML(2)-1/tau_beta-1/tau_n3);
%                    if phi0_s==0
%                        D_phi_MLE(k) = phi0 - betaHatML(6);
%                    else
%                        D_phi_MLE(k) = 0;
%                    end
%                    MLEtime(k) = toc - (timeFit+LStime(k)+RLStime(k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    %MLEnegloglik = @(params,d,cens,freq) -sum(data.*log(FIDsignal(params,t))-FIDsignal(params,t));
                    %MLEnegloglik(beta0,data);
                    %options = statset('MaxIter',10000, 'MaxFunEvals',Inf);
                    %[phat pci] = mle(data,'nloglf',MLEnegloglik,'start',beta0,'alpha',0.37,'options',options);
                    %disp('MLE Poisson Nonlinear regression:')
                    %disp(['Coeff  =  ' num2str(phat,'%0.8e  ')])
                    %disp(['Low    =  ' num2str(pci(1,:),'%0.8e  ')])
                    %disp(['High   =  ' num2str(pci(2,:),'%0.8e  ')])
                    
%DD%%%%%%%%%%%%%     DOUBLE Exponential Fitting     %%%%%%%%%%%%%%%%%%%%%%%%
                else
                    %--- double exponential decay fit
                    if phi0_s==1
                        FIDsignal = @(p,t) (abs(p(1)).*exp(-t./abs(p(2)))+abs(p(3)).*exp(-t./abs(p(4)))).*(1-abs(p(5))*cos(2*pi*abs(p(6))*t + phi0))+abs(p(7));%p(3)).*exp(-t./abs(p(4)))).*(abs(p(5))+abs(p(6))*(1-abs(p(7))*cos(2*pi*abs(p(8))*t + phi0)))+abs(p(9));
                        beta0 = [Nt(1)*t_width*((epsbeta/tau_beta)+(eps3/tau_n3)) 1/(1/tau_n3+1/880+2/tau0) Nt(1)*t_width*((epsbeta/tau_beta)+(eps3/tau_n3)) 1/(1/tau_n3+1/880+2/(3*tau0)) P3*Pn*(eps3/tau_n3)/((epsbeta/tau_beta)+(eps3/tau_n3)) f3n BG*t_width];
                        lb = zeros(1,length(beta0));ub = [50 1000 50 1000 5 10 5]; % bounds for mle function
                    else
                        FIDsignal = @(p,t) (abs(p(1)).*exp(-t./abs(p(2)))+abs(p(3)).*exp(-t./abs(p(4)))).*(1-abs(p(5))*cos(2*pi*abs(p(6))*t + abs(p(8))))+abs(p(7));%p(3)).*exp(-t./abs(p(4)))).*(abs(p(5))+abs(p(6))*(1-abs(p(7))*cos(2*pi*abs(p(8))*t + phi0)))+abs(p(9));
                        beta0 = [Nt(1)*t_width*((epsbeta/tau_beta)+(eps3/tau_n3)) 1/(1/tau_n3+1/880+2/tau0) Nt(1)*t_width*((epsbeta/tau_beta)+(eps3/tau_n3)) 1/(1/tau_n3+1/880+2/(3*tau0)) P3*Pn*(eps3/tau_n3)/((epsbeta/tau_beta)+(eps3/tau_n3)) f3n BG*t_width phi0];
                        lb = zeros(1,length(beta0));ub = [50 1000 50 1000 5 10 5 1];
                    end
                        %warning('off','all');
                    
                    %---- Least Squares fitting to double exponential
                    
            % fitnlm without "Robust" option 
            % fits Least Squares "linear"
                    lsmdl = fitnlm(t,data,FIDsignal,beta0);
            % delta_frequency after Least Squares fit 
            % f3n - 6th fitted coef.
                    D_f_LS(k) = f3n - lsmdl.Coefficients.Estimate(6);
            % sigma (standard error) for freq from LS fitting
                    S_f_LS(k) = lsmdl.Coefficients.SE(6);
            % delta_tau_walls after LS fit 
                    D_tau1_LS(k) = tau0 - 1/(1/lsmdl.Coefficients.Estimate(2)-1/tau_beta-1/tau_n3);
                    D_tau2_LS(k) = tau0 - 1/(1/lsmdl.Coefficients.Estimate(4)-1/tau_beta-1/tau_n3);
            % delta_phi0 after LS fit
                    if phi0_s==0
                        D_phi_LS(k) = phi0 - lsmdl.Coefficients.Estimate(8);
                    else
                        D_phi_LS(k) = 0;
                    end
                    LStime(k) = toc - timeFit;
                    
                    %---- Recursive Least Squares fitting to double exponential
                    [betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr);
                    S_f_RLS(k) = errBest2(6);
                    D_f_RLS(k) = f3n - betaBest2(6);
                    D_tau1_RLS(k) = tau0 - 1/(1/betaBest2(2)-1/tau_beta-1/tau_n3);
                    D_tau2_RLS(k) = tau0 - 1/(1/betaBest2(4)-1/tau_beta-1/tau_n3);
                    if phi0_s==0
                        D_phi_RLS(k) = phi0 - betaBest2(8);
                    else
                        D_phi_RLS(k) = 0;
                    end
                    %disp([num2str(tau_n3) ', ' num2str(Tm) ', ' num2str(Tf) ', ' num2str(sigma_f,'%8.2e') ', ' num2str(diff_f_RLS,'%8.2e') ', ' num2str(chi2_3,'%8.4f') ', ' num2str(finalSensitivity(sigma_f,Tm,Tf,Td,Efield),'%8.2e')]);
                    RLStime(k) = toc - (timeFit+LStime(k));
                    
                    %---- log likelihood (Poisson Nonlinear Regression)
 %%% MLE currently not working for Double Exp model
            % negative log likelihood for FIDsignal
                    mynegloglik = @(beta,data1,cens,freq) -sum(log(poisspdf(data,FIDsignal(beta,t))));
            % mle fitting defined by Matlab mle
                    [phatMLE,pciMLE] = mle(data,'nloglf',mynegloglik,'start',beta0); %,'LowerBound',lb,'UpperBound',ub);
                    phatMLE
                    pciMLE
            % extracted values from Matlab mle
                    D_f_MLE(k) = f3n - phatMLE(6);
                    S_f_MLE(k) = phatMLE(6) - pciMLE(1,6);
                    D_tau1_MLE(k) = tau0 - 1/(1/phatMLE(2)-1/tau_beta-1/tau_n3);
                    D_tau2_MLE(k) = tau0 - 1/(1/phatMLE(4)-1/tau_beta-1/tau_n3);
                    if phi0_s==0
                        D_phi_MLE(k) = phi0 - phatMLE(8);
                    else
                        D_phi_MLE(k) = 0;
                    end
                    MLEtime(k) = toc - (timeFit+LStime(k)+RLStime(k));
                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    mynegloglik = @(beta) -sum(log(poisspdf(data,FIDsignal(beta,t))));
%                    opts = optimset('fminsearch');
%                    opts.MaxFunEvals = Inf;
%                    opts.MaxIter = 10000;
%                    betaHatML = fminsearch(mynegloglik,beta0,opts);
%                    %disp('Poisson Nonlinear regression:')
%                    %disp(['Coeff  =  ' num2str(betaHatML,'%0.8e  ')])
%                    D_f_MLE(k) = f3n - betaHatML(6);
%                    MLEresiduals = FIDsignal(betaHatML,t)-data;
%                    [muMLE,SigmaMLE] = normfit(MLEresiduals);
%                    S_f_MLE(k) = SigmaMLE;   %%%%%%%%%%%% sigma MLE %%%%%%%%%
%                    D_tau1_RLS(k) = 2000 - 1/(1/betaHatML(2)-1/tau_beta-1/tau_n3);
%                    D_tau2_RLS(k) = 2000 - 1/(1/betaHatML(4)-1/tau_beta-1/tau_n3);
%                    if phi0_s==0
%                        D_phi_MLE(k) = phi0 - betaHatML(8);
%                    else
%                        D_phi_MLE(k) = 0;
%                    end
%                    MLEtime(k) = toc - (timeFit+LStime(k)+RLStime(k));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                    
%DWS%DWS%%%%%%%%%   SMOOTHED DOUBLE Exponential fitting  %%%%%%%%%%%%%%%%
                    %--- for double exponential decay fit - smoothing first - fix all 4 double exp parameters
                    %             disp('Double-exponential: smooth first, fixed all 4 double exponential parameters:')
                    %             DE = smoothDataAndDoubleExpFit(t,data,tau_n3);
                    %             FIDsignal = @(p,t) (DE(1).*exp(-t./DE(2))+DE(3).*exp(-t./DE(4))).*(abs(p(1))+abs(p(2))*(1-abs(p(3))*cos(2*pi*p(4)*t + phi0)))+abs(p(5));
                    %             beta0 = [epsbeta/tau_beta eps3/tau_n3 P3*Pn f3n BG*t_width];
                    %             [betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr);
                    %             sigma_f = errBest2(4);
                    %             diff_f_RLS = f3n - betaBest2(4);
                    %             sigma_f_RLS_SmoothFix4ExpParams(i) = sigma_f;
                    %             diff_f_RLS_SmoothFix4ExpParams(i) = diff_f_RLS;
                    %             disp([num2str(tau_n3) ', ' num2str(Tm) ', ' num2str(Tf) ', ' num2str(sigma_f,'%8.2e') ', ' num2str(diff_f_RLS,'%8.2e') ', ' num2str(chi2_3,'%8.4f') ', ' num2str(finalSensitivity(sigma_f,Tm,Tf,Td,Efield),'%8.2e')]);
                    %
                    
                    %--- for double exponential decay fit - smoothing first - fix all 2 double exp parameters
                    %disp('Double-exponential: smooth first, fixed 2 x time constants, leave free 2x amplitudes: \n')
                    %DE = smoothDataAndDoubleExpFit(t,data,tau_n3);
                    %FIDsignal = @(p,t) (abs(p(1)).*exp(-t./DE(2))+abs(p(2)).*exp(-t./DE(4))).*(abs(p(3))+abs(p(4))*(1-abs(p(5))*cos(2*pi*p(6)*t + phi0)))+abs(p(7));
                    %beta0 = [DE(1) DE(3) epsbeta/tau_beta eps3/tau_n3 P3*Pn f3n BG*t_width];
                    %[betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr);
                    %sigma_f = errBest2(6);
                    %diff_f_RLS = f3n - betaBest2(6);
                    %disp([num2str(tau_n3) ', ' num2str(Tm) ', ' num2str(Tf) ', ' num2str(sigma_f,'%8.2e') ', ' num2str(diff_f_RLS,'%8.2e') ', ' num2str(chi2_2,'%8.4f') ', ' num2str(finalSensitivity(sigma_f,Tm,Tf,Td,Efield),'%8.2e')]);
                    
                    
                    %--- from Filipone2009: (90% CL) sigma_d = h * 1.64 * sigma_f / (4* Efield) * sqrt(2)
                end
                
                % percent done
                progress = iProg+iiProg+jProg + doneRes*k
                
            end
%%%%%% DATA to be saved  %%%%%%
            %D_f_'fit' ; & D_phi_'fit' & D_tau1_'fit' & D_tau2_'fit'
            % tau1 & tau2 are fits for 2 coefficients in beta0
            Diffs = [transpose(D_f_LS),transpose(D_f_RLS),transpose(D_f_MLE),transpose(D_phi_LS),transpose(D_phi_RLS),transpose(D_phi_MLE),transpose(D_tau1_LS),transpose(D_tau1_RLS),transpose(D_tau1_MLE),transpose(D_tau2_LS),transpose(D_tau2_RLS),transpose(D_tau2_MLE)];
            times = [mean(LStime),mean(RLStime),mean(MLEtime)];
            Sigmas = [mean(S_f_LS) mean(S_f_RLS) mean(S_f_MLE) std(S_f_LS) std(S_f_RLS) std(S_f_MLE) times];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
            
            % FILENAME %
            % ---- save data - 3 arrays of diffs - LS vs RLS vs MLE - sized by number
            %      of loops, k
            
            filename = 'diffs';
            filenameSig = 'sigmas';
            if expFit_s==1      % exponential fitting parameter -- 1==single ; 2==double
                filename = strcat(filename,'_Ex1');
                filenameSig = strcat(filenameSig,'_Ex1');
            else
                filename = strcat(filename,'_Ex2');
                filenameSig = strcat(filenameSig,'_Ex2');
            end
            
            if phi0_s==1      % phi0 parameter -- 0==free ; 1==fixed value
                filename = strcat(filename,'_Ph1');
                filenameSig = strcat(filenameSig,'_Ph1');
            else
                filename = strcat(filename,'_Ph0');
                filenameSig = strcat(filenameSig,'_Ph0');
            end
            
            if tau_walls_s==0 % Tau_walls parameter -- 0==(E)dependent ; 1==single valued
                filename = strcat(filename,'_Tw0');
                filenameSig = strcat(filenameSig,'_Tw0');
            else
                filename = strcat(filename,'_Tw1');
                filenameSig = strcat(filenameSig,'_Tw1');
            end
            
            if f_walls_s==1   % f_walls parameter -- 1==1.0e-5 ; 2==2.0e-5
                filename = strcat(filename,'_Fw1');
                filenameSig = strcat(filenameSig,'_Fw1');
            else
                filename = strcat(filename,'_Fw2');
                filenameSig = strcat(filenameSig,'_Fw2');
            end
            
            if patch_s==0     % weak_patch parameter -- 0==no patch ; 1==patch considered
                filename = strcat(filename,'_WP0');
                filenameSig = strcat(filenameSig,'_WP0');
            else
                filename = strcat(filename,'_WP1');
                filenameSig = strcat(filenameSig,'_WP1');
            end
            filename = strcat(filename,'_');
            filename = strcat(filename,num2str(numLoops));
            filename = strcat(filename,'.txt');
            
            filenameSig = strcat(filenameSig,'_');
            filenameSig = strcat(filenameSig,num2str(numLoops));
            filenameSig = strcat(filenameSig,'.txt');
            
%%%%%%%%%%%%%  SAVE   %%%%%%%%%%%%%%
            save(filename,'Diffs','-ascii','-tabs')
            save(filenameSig,'Sigmas','-ascii','-tabs')
        end
    end
end

toc

%finalizePlots()
end

%%%%%   RLS    %%%%%
function [betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr)

mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
% disp('First Fit:')
% disp(['beta0  =  ' num2str(beta0,'%0.8e  ')])
% disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
% disp(['StdErr =  ' num2str(mdl.Coefficients.SE','%0.8e  ')])
betaBest1 = mdl.Coefficients.Estimate;
errBest1 = mdl.Coefficients.SE;
chi2_1 = sum(((data-FIDsignal(betaBest1,t))./dataErr).^2)/(length(t)-length(beta0));

dataErrRecurse = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErrRecurse.^2); %weights = 1/variance
%disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
dataErrRecurse = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErrRecurse.^2); %weights = 1/variance
%disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
dataErrRecurse = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErrRecurse.^2); %weights = 1/variance
%disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
dataErrRecurse = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErrRecurse.^2); %weights = 1/variance
%plot(linspace(t(1),t(end),10*n),FIDsignal(mdl.Coefficients.Estimate,linspace(t(1),t(end),10*n)),'r','linewidth',1.5);

betaBest2 = mdl.Coefficients.Estimate;
errBest2 = mdl.Coefficients.SE;
chi2_2 = sum(((data-FIDsignal(betaBest2,t))./dataErrRecurse).^2)/(length(t)-length(beta0));

chi2_3 = sum(((data-FIDsignal(betaBest2,t))./dataErr).^2)/(length(t)-length(beta0));
%CM = mdl.CoefficientCovariance
%SE = diag(sqrt(CM))
%how these errors are calculated in https://www.mathworks.com/help/stats/coefficient-standard-errors-and-confidence-intervals.html

%disp('Final recursive fit:')
%disp(['beta0  =  ' num2str(beta0,'%0.8e  ')])
%disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
%disp(['StdErr =  ' num2str(mdl.Coefficients.SE','%0.8e  ')])

end

function smoothedBestBeta = smoothDataAndDoubleExpFit(t,data,tau_n3)
global plotsOn
averageTime = 2;
numToAverage = floor(averageTime/20E-3);
for i = 1:floor(length(data)/numToAverage)
    startIndex = 1+(i-1)*numToAverage;
    endIndex = i*numToAverage;
    smoothedData(i) = sum(data(startIndex:endIndex))/numToAverage;
    smoothedTime(i) = t(floor((startIndex+endIndex)/2));
end

doubleExpDecay = @(p,t) abs(p(1)).*exp(-t/abs(p(2)))+abs(p(3)).*exp(-t/abs(p(4)));
beta0UCNdecay = [0.5*smoothedData(1)/2 1/(1/tau_n3+1/880+1/2000) 0.5*smoothedData(1)/2 1/(1/tau_n3+1/880+1/500)];
mdl = fitnlm(smoothedTime,smoothedData,doubleExpDecay,beta0UCNdecay);
smoothedBestBeta = abs(mdl.Coefficients.Estimate);
%disp('UCN decay fit:')
%disp(['beta0  =  ' num2str(beta0UCNdecay,'%0.8e  ')])
%disp(['Coeff  =  ' num2str(smoothedBestBeta','%0.8e  ')])
%disp(['StdErr =  ' num2str(mdl.Coefficients.SE','%0.8e  ')])

if plotsOn ==1
    figure(3);
    plot(smoothedTime,smoothedData,'-o','linewidth',2,'Color','r');
    plot(smoothedTime,doubleExpDecay(smoothedBestBeta,smoothedTime),'-','linewidth',2,'Color','g');
    drawnow
end

end

function nEDM90percentCI300days = finalSensitivity(sigma_f,Tm,Tf,Td,Efield)
sigma_dn_90CL = 6.63E-34*1.64*sqrt(2)/(4*Efield)*sigma_f/1.6e-19*100; %[e.cm] = [C.m] /1.6E-19[C/e] *100[cm/m]
numPer300days = 60*60*24*300/(Tm+Tf+Td);
nEDM90percentCI300days = sigma_dn_90CL/sqrt(numPer300days);
end

function initializePlots()
global Emax t_width
figure(1); clf; hold all; box on;
plot([Emax Emax],get(gca,'Ylim'),'k');
set(gca,'Xlim',[0 Emax*1.1],'Ylim',[0 1.1]);
xlabel('$E$ [neV]'); ylabel('production spectrum P(E)')
set(gcf,'units','centimeters'); pos = get(gcf,'position'); set(gcf,'position',[pos(1) 0 10 8]);

figure(2); clf; hold all; box on;
plot([Emax Emax],[0 10000],'k');
set(gca,'Xlim',[0 210],'Ylim',[0 4000]);
title('$f_{\rm walls}$ = 1E-5')
xlabel('$E$ [neV]'); ylabel('accumulated spectrum $N_0(E)$ [neV$^{-1}$]')
set(gcf,'units','centimeters'); pos = get(gcf,'position'); set(gcf,'position',[pos(1) 0 14 8]);

figure(3); clf; hold all; box on;
set(gcf,'units','centimeters'); pos = get(gcf,'position'); set(gcf,'position',[pos(1) 0 14 12]);
xlabel('time [s]')
ylabel(['light signal [counts/timebin] (timebin = ' num2str(t_width*1000,'%4.0f') ' ms )'])

figure(4); clf; hold all;
set(gcf,'units','centimeters'); pos = get(gcf,'position'); set(gcf,'position',[pos(1) 0 10 8]);

end

function finalizePlots(a,b,c)

figure(1);
%pos = get(gcf,'Position'); set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
%print(gcf,'./figure1.pdf','-dpdf','-r0')
%pdfcrop('./figure1.pdf')

figure(2);
legend('dPS cut-off','200s','400s','600s','800s','1000s','1200s','1400s','Location','EastOutside')
pos = get(gcf,'Position'); set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
%print(gcf,'./figure2.pdf','-dpdf','-r0')
%pdfcrop('./figure2.pdf')

figure(3);
set(gca,'Xlim',[0 +Inf],'Ylim',[-1 +Inf])
pos = get(gcf,'Position'); set(gcf,'PaperPositionMode','Auto','PaperSize',[pos(3), pos(4)])
%print(gcf,'./figure3.pdf','-dpdf','-r0')
%pdfcrop('./figure3.pdf')
end
