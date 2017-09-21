clc; clear all; set(0,'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultLegendInterpreter', 'latex'); set(gcf,'units','centimeters'); linecolors = lines(20); set(0,'DefaultAxesFontSize', 14); set(0,'DefaultTextFontSize', 14);
figure(1); clf; hold all;
figure(2); clf; hold all;
figure(3); clf; hold all;

%--- variable operating parameters
tau_n3_scan = 300:100:600; %average for unpolarized n
%tau_n3_scan = 500; %average for unpolarized n
%Tm_scan = 200:200:800; %[s] measurement time
Tm_scan = 1000; %[s] measurement time
%Tf_scan = 500:200:1000; %[s] cold neutron fill time
Tf_scan = 1000; %[s] cold neutron fill time

BG = 5; %[s^-1] other background rates

for i = 1:length(tau_n3_scan)
    for j = 1:length(Tm_scan)
        for k = 1:length(Tf_scan)
            tau_n3 = tau_n3_scan(i);
            Tm = Tm_scan(j);
            Tf = Tf_scan(k);
            
            P_UCN = 0.26; %[UCN/cc/s]
            tau_beta = 885;
            tau_walls = 2000;
            tau_tot = (1/tau_beta+1/tau_n3+1/tau_walls)^-1; %average loss rate of UCNs
            tau_buildup = (1/tau_beta+1/tau_walls)^-1; %when accumulating, UCNs are aligned with 3He.
            Gammatot = 1/tau_tot;
            E = 74; %[kV/cm]
            Td = 400; %[s] dead time between cycles
            P3 = 0.98; % initial 3He pol
            Pn = 0.98; %initial UCN pol
            Gammap = 1/20000; %[s^-1] 3He and UCN depolarization rate
            eps3 = 0.93; % n-3He absorption detection efficency
            epsbeta = 0.5; %beta decay detection efficiency
            
            N0 = P_UCN*tau_buildup*3000*(1-exp(-Tf/tau_buildup));
            F = eps3*P3*Pn/(tau_n3*(epsbeta/tau_beta+eps3/tau_n3));
            B0 = 30E-3* 1E-4; %[T] = [G]* 1E-4 [T/G]
            gamma_n = 1.83247172e8; %s^-1*T^-1 or 18.32472 [Hz/mG] * 1E3 [mG/G] * 1E4 [G/T]
            gamma_3 = 2.037894659e8; %s^-1*T^-1  [CODATA]
            omega3n = (gamma_3-gamma_n)*B0;
            f3n = omega3n/(2*pi);
            
            t_width = 10E-3; %[s]
            tend = 1000; %[s]
            n = round(tend/t_width);
            t = linspace(0,tend,n);
            phi0 = 10*pi/180;

            y = ( N0*(epsbeta/tau_beta+eps3/tau_n3)*exp(-Gammatot*t).*(1 - F.*cos(2*pi*f3n*t+phi0))+BG)*t_width;
            
            figure(1); clf; hold all;
            %plot(t,y,'-')
            
            figure(2); clf; hold all; box on;
            %set(gca,'Xlim',[0 1000])
            %xlabel('time [s]')
            %ylabel(['light signal [counts/timebin] (timebin = ' num2str(t_width*1000,'%4.0f') ' ms )'])
            
            %first guess at error assuming gaussian
            data = poissrnd(y);
            dataErr = sqrt(data);
            dataErr(dataErr<=0)=1; %have to save the case when error = 0
            
            figure(2); hold all; box on;
            h=errorbar(t,data,dataErr,'o','color',linecolors(1,:));
            drawnow
            errorbar_noHorzTab(h,1E10);
            drawnow
                
            %--- with unknown phase
            FIDsignal = @(p,t) p(1)*exp(-t./p(2)).*(1-p(3)*cos(2*pi*p(4)*t + p(5)))+ abs(p(6));
            beta0 = [N0*t_width*(epsbeta/tau_beta+eps3/tau_n3) 1/(1/tau_n3+1/880+1/2000) F f3n 10*pi/180 BG*t_width];
            
            mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
            %---recursive least-squares fitting
            %disp('Recursive least squares fitting:')
            dataErr = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
            mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
            %disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
            dataErr = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
            mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
            %disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
            dataErr = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
            mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
            %disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
            dataErr = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
            mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
            %disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
            dataErr = sqrt(FIDsignal(mdl.Coefficients.Estimate,t));
            mdl = fitnlm(t,data,FIDsignal,beta0,'Weights',1./dataErr.^2); %weights = 1/variance
            %disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
            plot(linspace(t(1),t(end),10*n),FIDsignal(mdl.Coefficients.Estimate,linspace(t(1),t(end),10*n)),'r','linewidth',1.5);
            %set(gca,'Xlim',[0 5])
            %disp(' ')
            %disp('Final fit:')
            disp(['beta0  =  ' num2str(beta0,'%0.8e  ')])
            disp(['Coeff  =  ' num2str(mdl.Coefficients.Estimate','%0.8e  ')])
            disp(['StdErr =  ' num2str(mdl.Coefficients.SE','%0.8e  ')])
            
            sigma_f = mdl.Coefficients.SE(4);
            sigma_dn_90CL = 6.63E-34*1.64*sqrt(2)/(4*74*1E3*100)*sigma_f/1.6e-19*100; %
            
            numPer300days = 60*60*24*300/(Tm+Tf+Td);
            sigma_f/sqrt(numPer300days);
            %disp(['Number Cycles per 300 days live time = ' num2str(numPer300days)]);
            %disp(['sigma_dn_90CL_300days = ' num2str(sigma_dn_90CL/sqrt(numPer300days))]);
            %disp(['tau_n3 = ' num2str(tau_n3) ', Tm = ' num2str(Tm) ', Tf = ' num2str(Tf) ', sigma_dn_90CL_1shot = ' num2str(sigma_dn_90CL,'%8.2e') ', sigma_dn_90CL_300days = ' num2str(sigma_dn_90CL/sqrt(numPer300days),'%8.2e')]);
            disp([num2str(tau_n3) ', ' num2str(Tm) ', ' num2str(Tf) ', ' num2str(sigma_dn_90CL,'%8.2e') ', ' num2str(sigma_dn_90CL/sqrt(numPer300days),'%8.2e')]);
        end
    end
end
