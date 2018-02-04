function CBarrow_FrPrecStatSens_VarsLoops01
%--- Calculating statistical impact of SNS nEDM experiment for f-factor and f-factor w peak patch loss model instead of single tau
clc; set(0,'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultLegendInterpreter', 'latex'); set(gcf,'units','centimeters'); linecolors = lines(20); set(0,'DefaultAxesFontSize', 14); set(0,'DefaultTextFontSize', 14);
%warning('on','all');
warning('off','all');
global plotsOn; plotsOn = 0;

% Test Variable Tree
 % ---- single exponential vs double exponential
  % ---- phi0 fixed vs phi0 free
   % ---- Tau_walls single valued vs Tau_walls (E) dependent
    % ---- f_walls = 1e-5 vs 2e-5
     % ---- w/ weak_patch vs w/o weak_patch
     
%---- Test Variable Settings
expFit_s = 1;    % exponential fitting parameter -- 1==single ; 2==double
phi0_s = 1;      % phi0 parameter -- 0==free ; 1==fixed value
tau_walls_s = 0; % Tau_walls parameter -- 0==(E)dependent ; 1==single valued
f_walls_s = 1;   % f_walls parameter -- 1==1.0e-5 ; 2==2.0e-5
patch_s = 0;     % weak_patch parameter -- 0==no patch ; 1==patch considered

%--- variable operating parameters
% tau_n3_scan = 400:100:700; %average for unpolarized n
% Tm_scan = 700:100:1300; %[s] measurement time
% Tf_scan = 700:100:1300; %[s] cold neutron fill time
tau_n3_scan = 400:100:500; %average for unpolarized n
Tm_scan = 1000; %[s] measurement time
Tf_scan = 400:200:1400; %[s] cold neutron fill time

%--- fixed parameters
global Emax t_width linecolors
VCell = 40.0*10.1*7.6; %[cm^3]
AWalls = 10.1*7.6*2+40.0*7.6*2+40.0*10.1*2; %[cm^2]
VoptWalls = 160; %[neV]
%VoptWalls = 200; %[neV]
fWalls = f_walls_s*10e-5;

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
if tau_walls_s==0
    if patch_s==0
        tau_walls = 4*VCell*1e-6./(mubarWalls.*v*AWalls*1e-4);
    else
        tau_walls = 1./(1./tau_walls+1./tau_patch);
    end
else
   tau_walls = 2000*ones(size(v)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CHECK VALUE
end

mubarPatch = 2*fpatch*(Upatch./E.*asin(sqrt(E/Upatch))-sqrt(Upatch./E-1));
mubarPatch(E/Upatch>1)=1;
tau_patch = 4*VCell*1e-6./(mubarPatch.*v*Apatch*1e-4);
%tau_walls = 1./(1./tau_walls+1./tau_patch); %comment this line for no weak patch CHANGE HERE
%tau_walls = 2000*ones(size(v)); % [s] - if single-exponential decay, CHANGE HERE

P_UCN = 0.31*(VoptWalls/160)^(3/2); %[UCN/cc/s]
%P_UCN = 0.26*(VoptWalls/160)^(3/2); %[UCN/cc/s] %used in Brad F's calculations
tau_beta = 885;
tau_buildup = (1/tau_beta+1./tau_walls).^-1; %when accumulating, UCNs are aligned with 3He.
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
if phi0_s==1
   phi0 = 30*pi/180;
else
   phi0 = rand*pi/180;
end

t_width = 20E-3; %[s]

%initializePlots()
%figure(1);
%plot(E,n_E,'Color',linecolors(2,:));

tic
for i = 1
    for j = 1:length(Tm_scan)
        for k = 1:length(Tf_scan)
            
            tau_n3 = tau_n3_scan(1);
            Tm = Tm_scan(j);
            Tf = Tf_scan(k);
            
            N0E = PE.*tau_buildup.*(1-exp(-Tf./tau_buildup)); % N0(E)[UCN/neV]
            tau_tot = (1/tau_beta+1/tau_n3+1./tau_walls).^-1; %average loss rate of UCNs
            Gammatot = 1./tau_tot;
            
            n = round(Tm/t_width);
            t = linspace(0,Tm,n);
            
            %figure(2);
            %plot(E,N0E,'linewidth',2);
            
            N0 = sum(N0E*dE);
            clear NEt;
            NEt = zeros([length(N0E),length(t)]);
            for l = 1:length(N0E)
                for m = 1:length(t)
                    NEt(l,m) = N0E(l).*exp(-Gammatot(l)*t(m)); %Z(E,t)
                end
            end
            Nt = sum(NEt,1)*dE;
            %disp(num2str(Nt(1)))
            y = (Nt.*(epsbeta/tau_beta + eps3/tau_n3.*(1-P3*Pn*cos(2*pi*f3n*t+phi0)))+BG)*t_width;
            
            data = random('Poisson',y);
            dataErr = sqrt(data); %first guess at error assuming gaussian
            dataErr(dataErr<=0)=1; %have to save the case when error = 0
            %if plotsOn ==1
            %    figure(3);
             %   h=errorbar(t,data,dataErr,'o','color',linecolors(1,:));
              %  drawnow
               % errorbar_noHorzTab(h,1E10);
                %drawnow
            %end
            
            %---- single exponential fit
            %disp('Single-exponential:')
            FIDsignal = @(p,t) p(1).*exp(-t./p(2)).*(p(3)+p(4)*(1-p(5)*cos(2*pi*p(6)*t + phi0)))+ abs(p(7));
            beta0 = [Nt(1)*t_width 1/(1/tau_n3+1/880+1/2000) epsbeta/tau_beta eps3/tau_n3 P3*Pn f3n BG*t_width];
            FIDsignal = @(p,t) p(1).*exp(-t./p(2)).*(p(3)+p(4)*(1-p(5)*cos(2*pi*p(6)*t + p(8))))+ abs(p(7));
            beta0 = [Nt(1)*t_width 1/(1/tau_n3+1/880+1/2000) epsbeta/tau_beta eps3/tau_n3 P3*Pn f3n BG*t_width 30*pi/180];
            
            %---- Least Squares fitting to single exponential
            LSray = lsqcurvefit(FIDsignal,beta0,t,data);  % returns array of coefficients from FIDsignal single exponential
            diff_f_LS = f3n - LSray(6);
            diff_f_LS_SingleExp(i) = diff_f_LS;
            
            %---- Recursive Least Squares fitting to single exponential
            [betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr);
            sigma_f = errBest2(6);
            diff_f_RLS = f3n - betaBest2(6);
            sigma_f_RLS_SingleExp(i) = sigma_f;
            diff_f_RLS_SingleExp(i) = diff_f_RLS;
            %disp([num2str(tau_n3) ', ' num2str(Tm) ', ' num2str(Tf) ', ' num2str(sigma_f,'%8.2e') ', ' num2str(diff_f_RLS,'%8.2e') ', ' num2str(chi2_3,'%8.4f') ', ' num2str(finalSensitivity(sigma_f,Tm,Tf,Td,Efield),'%8.2e')]);
            
            %---- log likelihood (Poisson Nonlinear Regression)
            mynegloglik = @(beta) -sum(log(poisspdf(data,FIDsignal(beta,t))));
            opts = optimset('fminsearch');
            opts.MaxFunEvals = Inf;
            opts.MaxIter = 10000;
            betaHatML = fminsearch(mynegloglik,beta0,opts);
            %disp('Poisson Nonlinear regression:')
            %disp(['Coeff  =  ' num2str(betaHatML,'%0.8e  ')])
            diff_f_ML(i) = f3n - betaHatML(6);
            
            %MLEnegloglik = @(params,d,cens,freq) -sum(data.*log(FIDsignal(params,t))-FIDsignal(params,t));
            %MLEnegloglik(beta0,data);
            %options = statset('MaxIter',10000, 'MaxFunEvals',Inf);
            %[phat pci] = mle(data,'nloglf',MLEnegloglik,'start',beta0,'alpha',0.37,'options',options);
            %disp('MLE Poisson Nonlinear regression:')
            %disp(['Coeff  =  ' num2str(phat,'%0.8e  ')])
            %disp(['Low    =  ' num2str(pci(1,:),'%0.8e  ')])
            %disp(['High   =  ' num2str(pci(2,:),'%0.8e  ')])
            
            %--- double exponential decay fit
%             disp('Double-exponential: all 4 parameters free')
%             FIDsignal = @(p,t) (abs(p(1)).*exp(-t./abs(p(2)))+abs(p(3)).*exp(-t./abs(p(4)))).*(abs(p(5))+abs(p(6))*(1-abs(p(7))*cos(2*pi*abs(p(8))*t + phi0)))+abs(p(9));
%             beta0 = [0.5*Nt(1)*t_width 1/(1/tau_n3+1/880+1/2000) 0.5*Nt(1)*t_width 1/(1/tau_n3+1/880+1/200) epsbeta/tau_beta eps3/tau_n3 P3*Pn f3n BG*t_width];
%             warning('off','all');
%             [betaBest1,errBest1,chi2_1,betaBest2,errBest2,chi2_2,chi2_3] = recursiveLeastSquaresFitting(t,data,FIDsignal,beta0,dataErr);
%             sigma_f = errBest2(8);
%             diff_f_RLS = f3n - betaBest2(8);
%             sigma_f_RLS_DblExp4ParamsFree(i) = sigma_f;
%             diff_f_RLS_DblExp4ParamsFree(i) = diff_f_RLS;
%             disp([num2str(tau_n3) ', ' num2str(Tm) ', ' num2str(Tf) ', ' num2str(sigma_f,'%8.2e') ', ' num2str(diff_f_RLS,'%8.2e') ', ' num2str(chi2_3,'%8.4f') ', ' num2str(finalSensitivity(sigma_f,Tm,Tf,Td,Efield),'%8.2e')]);
%            
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
    end
    D_f_LS(i) = diff_f_LS_SingleExp(i);  % populates array of frequency difference (Least Squares fitting)
    D_f_RLS(i) = diff_f_RLS_SingleExp(i);  % populates array of frequency difference (Recursive Least Squares)
    D_f_MLE(i) = diff_f_ML(i);  % populates array of frequency difference (MLE with Poisson distribution)
end
D_f_LS,D_f_RLS,D_f_MLE
toc
%save 2016-12-07_SmoothFix4ExpParams_400_1000_1000_1E-5_200neV
%finalizePlots()
end

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

function finalizePlots()

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
