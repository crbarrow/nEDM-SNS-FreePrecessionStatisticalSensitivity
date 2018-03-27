clear all; set(0,'defaultTextInterpreter','latex'); set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'DefaultLegendInterpreter', 'latex'); set(gcf,'units','centimeters'); linecolors = lines(20); set(0,'DefaultAxesFontSize', 14); set(0,'DefaultTextFontSize', 14);

% pull files from directory
filenames = dir('C:\Users\Chazz\Documents\Huffman Lab\codes');

% pull files with name containing 'diffs'
for i = 1:numel(filenames)
    if contains(filenames(i).name, 'diffs')

        % WorkFile pulls the name of the file to be processed
        WorkFile = filenames(i).name;
        %%% SigmaFile names the file of information related to the
        % distribution of plotted data and standard error produced by the
        % fitting algorithms for individual data points
        %%% SigmaFile was created when the 'diffs' data was created and is
        % merely called here
        SigmaFile = 'sigmas';
        % adds the appropriate Workfile info to 'sigmas' so that it may be
        % read with 'dlmread'
        SigmaFile = strcat(SigmaFile,WorkFile(6:length(WorkFile)));
        
        M = dlmread(WorkFile,'\t');
        % frequencies produced in the fitting routines
        f_LS = M(:,1);f_RLS = M(:,2);f_MLE = M(:,3);
        
        S = dlmread(SigmaFile,'\t');
        % loads info from 'sigmas_*'
        figure(1); clf; hold all; box on; grid on;set(figure(1), 'Visible', 'off');

        % size the number of elements for each population and prepare
        % histogram binning
        distMax = max(max(M(:,1))*1E6);
        distSpan = 2*distMax;
        bins = 50;
        distRes = distSpan/bins;

        bin_edges = -(distMax):distRes:distMax; %[1E-6 Hz]
        % create histograms representations for each fitting method 
        % population
        hist_diff_f_LS = histogram(f_LS*1E6,bin_edges);
        hist_diff_f_RLS = histogram(f_RLS*1E6,bin_edges);
        hist_diff_f_MLE = histogram(f_MLE*1E6,bin_edges);

        bin_centers = (bin_edges(1:end-1)+bin_edges(2:end))/2;
        % plot histograms
        figure(2); clf; hold all; box on; grid on;set(figure(2), 'Visible', 'off'); % visible off for scan program
        plot(bin_centers,hist_diff_f_LS.Values,'-o','LineWidth',2,'Color',linecolors(1,:))
        plot(bin_centers,hist_diff_f_RLS.Values,'-o','LineWidth',2,'Color',linecolors(2,:))
        plot(bin_centers,hist_diff_f_MLE.Values,'-o','LineWidth',2,'Color',linecolors(3,:))
        xlabel('$f_{\rm true}- f_{\rm extracted}$ [${\rm \mu Hz}$]')
        ylabel(['number per ' num2str(distRes,0.2) '${\rm \mu Hz}$'])
        title('Histogram of difference between true frequency and extracted frequency')
        legend('LS','RLS','MLE')

        %%% STD DEV of population
        % analyze histogram representations
        % produce mean and standard deviation of population
        [mu_LS,sigma_LS,muci,sigmaci] = normfit(hist_diff_f_LS.Data);
        sigma_LS;
        %%% pm_* = confidence interval for Mean of plotted population
        pm_LS = sigma_LS/(sqrt(length(M(:,1))));
        % repeated for LS, RLS, and MLE
        [mu_RLS,sigma_RLS,muci,sigmaci] = normfit(hist_diff_f_RLS.Data);
        sigma_RLS;
        pm_RLS = sigma_RLS/(sqrt(length(M(:,1))));
        [mu_ML,sigma_ML,muci,sigmaci] = normfit(hist_diff_f_MLE.Data);
        sigma_ML;
        pm_ML = sigma_ML/(sqrt(length(M(:,1))));

        % create a multiplier that will scale the Gaussian line
        % to the plotted data points
        mult = distMax*(length(f_LS)*2/bins);

        % plot the data points and the Gaussian line
        figure(3); clf; hold all; box on; grid on;set(figure(3), 'Visible', 'off'); % visible off for scan program
        plot(bin_centers,hist_diff_f_LS.Values,'o','LineWidth',2,'Color',linecolors(1,:));
        h1=plot(bin_centers,normpdf(bin_centers,mu_LS,sigma_LS)*mult,'-','LineWidth',2,'Color',linecolors(1,:));
        plot(bin_centers,hist_diff_f_RLS.Values,'o','LineWidth',2,'Color',linecolors(2,:));
        h2=plot(bin_centers,normpdf(bin_centers,mu_RLS,sigma_RLS)*mult,'-','LineWidth',2,'Color',linecolors(2,:));
        plot(bin_centers,hist_diff_f_MLE.Values,'o','LineWidth',2,'Color',linecolors(3,:));
        h3=plot(bin_centers,normpdf(bin_centers,mu_ML,sigma_ML)*mult,'-','LineWidth',2,'Color',linecolors(3,:));

        xlabel('$f_{\rm true}- f_{\rm extracted}$ [${\rm \mu Hz}$]')
        ylabel(['number per ' num2str(distRes,0.2) '${\rm \mu Hz}$'])
        %title('Histogram of difference between true frequency and extracted frequency')
        legend([h1 h2 h3],'LS','RLS','MLE')
        %set(gca,'xscale','log');    % set x-axis to logarithmic scaling

        % create matrix for extracted MEAN, SE and STD DEV info
        Sigmas = zeros(2,9);
        %%% STD - standard deviation of population  
        Sigmas(1,2) = sigma_LS; Sigmas(1,5) = sigma_RLS; Sigmas(1,8) = sigma_ML;
        %%% SE - standard error from fitting algorithm
        % usually averaged
        % note on MLE: algorithm gives a 95% confidence interval rather
        % than a sigma or std dev, which corresponds to 2 sigmas, so S(3)
        % must be divided by 2 to reflect a "sigma" value
        Sigmas(1,1) = S(1)*1.0e06; Sigmas(1,4) = S(2)*1.0e06; Sigmas(1,7) = (S(3)/2)*1.0e06;
        %%% confidence interval for SE
        Sigmas(2,1) = S(4)*1.0e06; Sigmas(2,4) = S(5)*1.0e06; Sigmas(2,7) = S(6)*1.0e06;
        %%% Mean from plotted population
        Sigmas(1,3) = mu_LS; Sigmas(1,6) = mu_RLS; Sigmas(1,9) = mu_ML;
        %%% Confidence interval of plotted population Mean
        Sigmas(2,3) = pm_LS; Sigmas(2,6) = pm_RLS; Sigmas(2,9) = pm_ML;
        % display Sigmas on Command window
        Sigmas

        % Create filename to save 'sigmas' info
        filenameSig = 'hSig';
        filenameSig = strcat(filenameSig,WorkFile(6:length(WorkFile)));

        save(filenameSig,'Sigmas','-ascii','-tabs')

        % save plotted images
        savdir = 'C:\Users\Chazz\Documents\Huffman Lab\codes\images';
        saveas(figure(3),fullfile(savdir,[WorkFile(1:(length(WorkFile)-4)) '.jpg']));
        clf
    else % if filename does not contain 'diffs'
        continue
    end
end

