function out = LinearModelFit(fullpath, data, CMIdata)
    % LinearModelfit: fit deltaR/R0 using exponential kernel
    %
    % USAGE:
    %       out = LinearModelFit(fullpath, data, CMIdata)
    % 
    % INPUTS:
    %       fullpath: location of the *.mat file from IOSTool GUI
    %       data:     data contains in *.mat file from IOSTool GUI
    %       CMIdata:  data from ComputeMeanIntensity.m
    %
    % OUTPUTS:
    %       out: data structure contains all the fitting results
    %
    % Last Modified: 07-29-2016
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)
    
    % See also:
    %       HRFModel.m
    %       getRGBcomponent.m
    %       bipolar.m
    %       screen2file.m
    %       Fit_HRF_ave_plot.m
    
    %% Get data for fitting
    dRR0filtall = CMIdata.dRR0filtall; 
                                        % Format of dRR0filtall
                                        % cell 1: left hemisphere
                                        % cell 2: right hemisphere
                                        % cell 3: left frontal
                                        % cell 4: right frontal
                                        % cell 5: left parietal
                                        % cell 6: right parietal
                                        % cell 7: left vibrissae
                                        % cell 8: right vibrissae
                                        % cell 9: left visual
                                        % cell 10: right visual
                                        % cell 11: left auditory
                                        % cell 12: right auditory
            
    xloc = CMIdata.xloc; % x coordinate
    yloc = CMIdata.yloc; % y coordinate
    
    s = data.speed.aveBinary; % binary speed event
    t = data.time;            % time, in seconds
    xLH = dRR0filtall{1};     % filted dR/R0 for whole trial, left hemisphere
    xRH = dRR0filtall{2};     % filted dR/R0 for whole trial, right hemisphere

    xLHidx = isempty(xLH);    % left window empty?
    xRHidx = isempty(xRH);    % right window empty?

    %% Fit data using Exponential Model
    if ~xLHidx % left window not empty
        disp('Linear Model Fit: Left Window Exist');
        [a, v, c, HRFa, HRFv, HRF, xFit_HRFa, xFit_HRFv, xFit_HRF, cc] = HRFModel(s,xLH,t);
        fitResults.LH.a = a;                    % arterial fitting weight parameter
        fitResults.LH.v = v;                    % venous fitting weight parameter
        fitResults.LH.c = c;                    % constant parameter
        fitResults.LH.HRFa = HRFa;              % arterial impulse response
        fitResults.LH.HRFv = HRFv;              % venous impulse response
        fitResults.LH.HRF = HRF;                % arterial + venous impulse response
        fitResults.LH.xFit_HRFa = xFit_HRFa;    % arterial fit
        fitResults.LH.xFit_HRFv = xFit_HRFv;    % venous fit
        fitResults.LH.xFit_HRF = xFit_HRF;      % arterial + venous fit
        fitResults.LH.cc = cc;                  % pixelwise correlation coefficient for the fitting
    end

    if ~xRHidx % right window not empty
        disp('Linear Model Fit: Right Window Exist');
        [a, v, c, HRFa, HRFv, HRF, xFit_HRFa, xFit_HRFv, xFit_HRF, cc] = HRFModel(s,xRH,t);
        fitResults.RH.a = a;
        fitResults.RH.v = v;
        fitResults.RH.c = c;
        fitResults.RH.HRFa = HRFa;
        fitResults.RH.HRFv = HRFv;
        fitResults.RH.HRF = HRF;
        fitResults.RH.xFit_HRFa = xFit_HRFa;
        fitResults.RH.xFit_HRFv = xFit_HRFv;
        fitResults.RH.xFit_HRF = xFit_HRF;
        fitResults.RH.cc = cc;
    end
    
    %% Re-organize fitting results for plotting figure
    % Pre-allocate parameters
    ccplot = NaN*zeros(256,256); % correlation coefficient plot
    ccLH = NaN*zeros(256,256); % correlation coefficient for LH
    ccRH = NaN*zeros(256,256); % correlation coefficient for RH
    
    % for plot purpose ONLY
    aImg = zeros(256, 256, 3); % artery image only, for figure purpose ONLY
    vImg = zeros(256, 256, 3); % vein image only, for figure purpose ONLY
    
    % for computational purpose, use the following parameters
    aWindow = NaN*zeros(256,256); % collect all a value in one matrix
    vWindow = NaN*zeros(256,256); % collect all v value in one matrix
    cWindow = NaN*zeros(256,256); % collect all c value in one matrix
    aLH = NaN*zeros(256, 256); % a, left hemisphere
    vLH = NaN*zeros(256, 256); % v, left hemisphere
    aRH = NaN*zeros(256, 256); % a, right hemisphere
    vRH = NaN*zeros(256, 256); % v, right hemisphere


    if ~xLHidx
        xlocLH = xloc{1}; % x coordinate of LH ROI
        ylocLH = yloc{1}; % y coordinate of LH ROI
        N = length(fitResults.LH.cc);
        aImgLH = getRGBcomponent(fitResults.LH.a); % normalize a, for plotting
        vImgLH = getRGBcomponent(fitResults.LH.v); % normalize v, for plotting
        for i = 1:N
            ccplot(xlocLH(i), ylocLH(i)) = fitResults.LH.cc(i);
            ccLH(xlocLH(i), ylocLH(i)) = fitResults.LH.cc(i);
            
            % for plot purpose ONLY
            aImg(xlocLH(i), ylocLH(i),1) = aImgLH(i); % Red color
            vImg(xlocLH(i), ylocLH(i),3) = vImgLH(i); % Blue color
            
            % for computational purpose
            aLH(xlocLH(i), ylocLH(i)) = fitResults.LH.a(i);
            vLH(xlocLH(i), ylocLH(i)) = fitResults.LH.v(i);
            
            aWindow(xlocLH(i), ylocLH(i)) = fitResults.LH.a(i);
            vWindow(xlocLH(i), ylocLH(i)) = fitResults.LH.v(i);            
            cWindow(xlocLH(i), ylocLH(i)) = fitResults.LH.c(i);
        end
        aLH_ave_AP = nanmean(aLH, 2); % mean 'a' for Left hemisphere, from anterior to posterior
        aLH_ave_ML = nanmean(aLH, 1); % mean 'a' for Left hemisphere, from midline to lateral
        aLH_ave_ML = aLH_ave_ML(:);
        
        vLH_ave_AP = nanmean(vLH, 2);
        vLH_ave_ML = nanmean(vLH, 1);
        vLH_ave_ML = vLH_ave_ML(:);
        
        ccLH_ave_AP = nanmean(ccLH, 2);
    else
        aLH_ave_AP = NaN*zeros(256,1);
        aLH_ave_ML = NaN*zeros(256,1);
        
        vLH_ave_AP = NaN*zeros(256,1);
        vLH_ave_ML = NaN*zeros(256,1);
        
        ccLH_ave_AP = NaN*zeros(256,1);
    end


    if ~xRHidx
        xlocRH = xloc{2}; % x coordinate of RH ROI
        ylocRH = yloc{2}; % y coordinate of RH ROI
        N = length(fitResults.RH.cc);
        aImgRH = getRGBcomponent(fitResults.RH.a); % normalize a
        vImgRH = getRGBcomponent(fitResults.RH.v); % normalize v
        for i = 1:N
            ccplot(xlocRH(i), ylocRH(i)) = fitResults.RH.cc(i);
            ccRH(xlocRH(i), ylocRH(i)) = fitResults.RH.cc(i);
            
            % for plot purpose ONLY
            aImg(xlocRH(i), ylocRH(i),1) = aImgRH(i);
            vImg(xlocRH(i), ylocRH(i),3) = vImgRH(i);
            
            % for computational purpose
            aRH(xlocRH(i), ylocRH(i)) = fitResults.RH.a(i);
            vRH(xlocRH(i), ylocRH(i)) = fitResults.RH.v(i);
            
            aWindow(xlocRH(i), ylocRH(i)) = fitResults.RH.a(i);
            vWindow(xlocRH(i), ylocRH(i)) = fitResults.RH.v(i);            
            cWindow(xlocRH(i), ylocRH(i)) = fitResults.RH.c(i);
        end
        aRH_ave_AP = nanmean(aRH, 2); % mean 'a' for Right hemisphere, from anterior to posterior
        aRH_ave_ML = nanmean(aRH, 1); % mean 'a' for Right hemisphere, from midline to lateral
        aRH_ave_ML = aRH_ave_ML(:);
        
        vRH_ave_AP = nanmean(vRH, 2);
        vRH_ave_ML = nanmean(vRH, 1);
        vRH_ave_ML = vRH_ave_ML(:);
        
        ccRH_ave_AP = nanmean(ccRH,2);
    else
        aRH_ave_AP = NaN*zeros(256,1);
        aRH_ave_ML = NaN*zeros(256,1);
        
        vRH_ave_AP = NaN*zeros(256,1);
        vRH_ave_ML = NaN*zeros(256,1);
        
        ccRH_ave_AP = NaN*zeros(256,1);
    end
    
    
    a_ave_ML = nanmean(aWindow, 1);
    a_ave_ML = a_ave_ML(:);
    
    v_ave_ML = nanmean(vWindow, 1);
    v_ave_ML = v_ave_ML(:);
    
    cc_ave_ML = nanmean(ccplot, 1);
    cc_ave_ML = cc_ave_ML(:);
    
    vesselImg = aImg + vImg; % combine a and v together, for plot purpose ONLY
    
    %% Organize output data
    out.fitResults = fitResults;
    
    % for plotting ONLY, artery, vein and vessel image
    out.fitMatrix.aImg = aImg;
    out.fitMatrix.vImg = vImg;
    out.fitMatrix.vesselImg = vesselImg;
    
    % for computational purpose, left, right and whole window
    out.fitMatrix.aLH = aLH;
    out.fitMatrix.vLH = vLH;
    
    out.fitMatrix.aRH = aRH;
    out.fitMatrix.vRH = vRH;
    
    out.fitMatrix.aWindow = aWindow;
    out.fitMatrix.vWindow = vWindow;
    out.fitMatrix.cWindow = cWindow;
    out.fitMatrix.ccplot = ccplot;
    
    % average for left hemisphere, AP and ML
    out.fitMatrix.aLH_ave_AP = aLH_ave_AP;
    out.fitMatrix.aLH_ave_ML = aLH_ave_ML;    
    out.fitMatrix.vLH_ave_AP = vLH_ave_AP;
    out.fitMatrix.vLH_ave_ML = vLH_ave_ML;
    
    out.fitMatrix.ccLH_ave_AP = ccLH_ave_AP;
    
    % average for right hemisphere, AP and ML
    out.fitMatrix.aRH_ave_AP = aRH_ave_AP;
    out.fitMatrix.aRH_ave_ML = aRH_ave_ML;
    out.fitMatrix.vRH_ave_AP = vRH_ave_AP;
    out.fitMatrix.vRH_ave_ML = vRH_ave_ML;

    out.fitMatrix.ccRH_ave_AP = ccRH_ave_AP;
        
    % average for whole window, ML ONLY
    out.fitMatrix.a_ave_ML = a_ave_ML;
    out.fitMatrix.v_ave_ML = v_ave_ML;
    out.fitMatrix.cc_ave_ML = cc_ave_ML;
    
    %% Change arterial map color (in the frontal cortex, make it yellow)
    a_tmp = getRGBcomponent(-aRH); % use the fit matrix data
    idx = a_tmp>0; % pixels with fitting parameters being negative
    tt = aImg(:,:,1);
    tt(idx) = a_tmp(idx);
    aImg(:,:,1) = tt; % set red channel
    
    tt2 = aImg(:,:,2);
    tt2(idx) = a_tmp(idx);
    aImg(:,:,2) = tt2;
    
    vesselImg = aImg + vImg;  
    
    %% 1. Plot a and v values to check whether they are positive or negative [OK]
    monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
    w = 1200; ht = 400; % Define window size: width(w) and height(ht)
    mon = 1; % Only one monitor
    
    % Determine positioning for main figure
    pos = [ ((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
            ((monitorPos(mon,4)-ht)/2)+monitorPos(mon,2)/2,...
            w ht];
    
    figure('name','Check fitting parameters','position',pos,'color','w');
    
    % ---- Plot arterial component image ---- %
    axes('units','normalized',...
        'position',[0.03 0.05 0.3 0.9],...
        'box','off',...
        'xtick',[],'ytick',[]);
    aWindowtmp = aWindow;
    aWindowtmp(isnan(aWindow))=0;
    imagesc(aWindowtmp);
    title('Arterial Component');
    set(gca, 'clim', [-max(abs(aWindowtmp(:))) max(abs(aWindowtmp(:)))]);
    axis xy image off;
    colorbar('location', 'southoutside', 'position', [0.03 0.03 0.3 0.02]);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','k','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','k','linewidth',2);
    end
    hold off;
        
    % ---- Plot venous component image ---- %
    axes('units','normalized',...
        'position',[0.35 0.05 0.3 0.9],...
        'box','off',...
        'xtick',[],'ytick',[]);
    vWindowtmp = vWindow;
    vWindowtmp(isnan(vWindowtmp)) = 0;
    imagesc(vWindowtmp);
    title('Venous Component');
    set(gca, 'clim', [-max(abs(vWindowtmp(:))) max(abs(vWindowtmp(:)))]);
    axis xy image off;
    colorbar('location', 'southoutside', 'position', [0.35 0.03 0.3 0.02]);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','k','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','k','linewidth',2);
    end
    hold off;
    
    % ---- Plot constant image ---- %
    axes('units','normalized',...
        'position',[0.67 0.05 0.3 0.9],...
        'box','off',...
        'xtick',[],'ytick',[]);
    cWindowtmp = cWindow;
    cWindowtmp(isnan(cWindow))=0;
    imagesc(cWindowtmp);
    title('Constant');
    set(gca, 'clim', [-max(abs(cWindowtmp(:))) max(abs(cWindowtmp(:)))]);
    axis xy image off;
    colorbar('location', 'southoutside', 'position', [0.67 0.03 0.3 0.02]);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','k','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','k','linewidth',2);
    end
    hold off;
    
    % set colormap for the whole figure;
    % if you want to set different colormaps for each axis, use
    % 'freezeColor' after each colormap setting
    % you can download 'freezeColor' from matlab fileexchange
    colormap bwr; % blue-white-red colormap
    
    % ---- Save figure to file ---- %
    screen2file(strrep(fullpath, '.mat', '_a_v_check'));   

    
    %% 3. Plot the average fitting results for whole hemisphere, including HRF
    if ~xLHidx % Left Hemisphere exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'LH');
    end
    
    if ~isempty(dRR0filtall{3}) % Left Frontal exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'LF');
    end
    
    if ~isempty(dRR0filtall{5}) % Left Parietal exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'LP');
    end
    
    if ~isempty(dRR0filtall{7}) % Left Vibrissae exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'LV');
    end
    
    if ~isempty(dRR0filtall{9}) % Left Visual exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'LE');
    end
    
    if ~isempty(dRR0filtall{11}) % Left Auditory exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'LA');
    end
    
    if ~xRHidx % Right Hemisphere exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'RH');
    end
    
    if ~isempty(dRR0filtall{4}) % Right Frontal exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'RF');
    end
    
    if ~isempty(dRR0filtall{6}) % Right Parietal exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'RP');
    end
    
    if ~isempty(dRR0filtall{8}) % Right Vibrissae exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'RV');
    end
    
    if ~isempty(dRR0filtall{10}) % Right Visual exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'RE');
    end
    
    if ~isempty(dRR0filtall{12}) % Right Auditory exist
        Fit_HRF_ave_plot(fullpath, data, fitResults, dRR0filtall, xloc, yloc, 'RA');
    end

    
    %% 4. Plot correlation coefficients btwn raw data and fitting data, including averaged cc along different direction
    monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
    w = 800; ht = 600; % Define window size: width(w) and height(ht)
    mon = 1; % Only one monitor

    % Determine positioning for main figure
    pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
        ((monitorPos(mon,4)-ht)/2)+monitorPos(mon,2)/2,...
        w ht];
    
    figure('name','Correlation Coefficient', 'position',pos,'color','w');
    % axes for anterior-posterior average (Left)
    if ~xLHidx
        axes('units','normalized',...
            'position',[0.05 0.37 0.21 0.56],...
            'box','off',...
            'xtick',[]);
        hold on;
        plot(ccLH_ave_AP,'color','k','linewidth',2);
        axis([1 256 0 1]); % 256 pixels, cc is ranging from 0 to 1
        xlim([1 256]);
        set(gca, 'xAxisLocation','top');
        view(-90,90); % switch x and y axis
        hold off;
    end
    
    
    % axes for correlation coefficient
    axes('units','normalized',...
        'position',[0.29 0.37 0.42 0.56],...
        'box','off',...
        'xtick',[],'ytick',[]);
    imagesc(ccplot);
    axis image xy off;
    colormap hot;
    colorbar('location', 'northoutside', 'position', [0.29 0.95 0.42 0.02]);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','w','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','w','linewidth',2);
    end
    hold off;
    
    % axes for midline-lateral average
    axes('units','normalized',...
        'position',[0.29 0.05 0.42 0.28],...
        'box','off',...
        'xtick',[]);
    hold on;
    plot(cc_ave_ML,'color','k','linewidth',2);
    axis([1 256 0 1]);
    xlim([1 256]);
    set(gca, 'ydir', 'reverse');
    hold off;
    
    % axes for anterior-posterior average (Right)
    if ~xRHidx
        axes('units','normalized',...
            'position',[0.74 0.37 0.21 0.56],...
            'box','off',...
            'xtick',[]);
        hold on;
        plot(ccRH_ave_AP,'color','k','linewidth',2);
        axis([1 256 0 1]);
        view(-90,90); % switch x and y axis
        set(gca, 'ydir', 'reverse');
        hold off;
    end
  
    % ---- Save figure to file ---- %
    screen2file(strrep(fullpath, '.mat', '_ccplot'));

    %% 5. Plot normalized a and v, including averaged a and v along different directions
    monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
    w = 1200; ht = 600; % Define window size: width(w) and height(ht)
    mon = 1; % Only one monitor

    % Determine positioning for main figure
    pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
        ((monitorPos(mon,4)-ht)/2)+monitorPos(mon,2)/2,...
        w ht];
    
    figure('name', 'normalized a and v parameters', 'position',pos,'color','w');
    
    % ---- a Image (i.e., arterial component)---- %
    % axes for anterior-posterior average (Left)
    if ~xLHidx
        axes('units','normalized',...
            'position',[0.01 0.32 0.06 0.48],...
            'box','off',...
            'xtick',[]);
        hold on;
        plot(aLH_ave_AP,'color','r','linewidth',2);
        hold on;
        plot(zeros(256,1),'.k');
        xlim([1 256]);
        set(gca, 'xAxisLocation','top','ydir','reverse');
        view(-90,90); % switch x and y axis
        hold off;
    end
    % axes for arterial component image
    axes('units','normalized',...
        'position',[0.07 0.32 0.24 0.48],...
        'box','off',...
        'xtick',[],'ytick',[]);
    imagesc(aImg);
    axis xy image off;
    title('Artery','fontsize',24);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','w','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','w','linewidth',2);
    end
    hold off;   
    % axes for midline-lateral average
    axes('units','normalized',...
        'position',[0.07 0.2 0.24 0.12],...
        'box','off',...
        'xtick',[]);
    hold on;
    plot(a_ave_ML,'color','r','linewidth',2);
    hold on;
    plot(zeros(256,1),'.k');
    xlim([1 256]);   
    hold off;
    % axes for anterior-posterior average (Right)
    if ~xRHidx
        axes('units','normalized',...
            'position',[0.31 0.32 0.06 0.48],...
            'box','off',...
            'xtick',[]);
        hold on;
        plot(aRH_ave_AP,'color','r','linewidth',2);
        hold on;
        plot(zeros(256,1),'.k');
        xlim([1 256]);
        view(-90,90); % switch x and y axis
        hold off;
    end
        
    % ---- v Image (i.e., venous component) ---- %
    % axes for anterior-posterior average (Left)
    if ~xLHidx
        axes('units','normalized',...
            'position',[0.38 0.32 0.06 0.48],...
            'box','off',...
            'xtick',[]);
        hold on;
        plot(vLH_ave_AP,'color','b','linewidth',2);
        hold on;
        plot(zeros(256,1),'.k');
        xlim([1 256]);
        set(gca, 'xAxisLocation','top','ydir','reverse');
        view(-90,90); % switch x and y axis
        hold off;
    end
    % axes for v image
    axes('units','normalized',...
        'position',[0.44 0.32 0.24 0.48],...
        'box','off',...
        'xtick',[],'ytick',[]);
    imagesc(vImg);
    axis xy image off;
    title('Vein','fontsize',24);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','w','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','w','linewidth',2);
    end
    hold off;
    % axes for midline-lateral average
    axes('units','normalized',...
        'position',[0.44 0.2 0.24 0.12],...
        'box','off',...
        'xtick',[]);
    hold on;
    plot(v_ave_ML,'color','b','linewidth',2);
    hold on;
    plot(zeros(256,1),'.k');
    xlim([1 256]);
    hold off;
    % axes for anterior-posterior average (Right)
    if ~xRHidx
        axes('units','normalized',...
            'position',[0.68 0.32 0.06 0.48],...
            'box','off',...
            'xtick',[]);
        hold on;
        plot(vRH_ave_AP,'color','b','linewidth',2);
        hold on;
        plot(zeros(256,1),'.k');
        xlim([1 256]);
        view(-90,90); % switch x and y axis
        hold off;
    end
    
    
    % ---- Vessel Image (a + v image) ---- %
    axes('units','normalized',...
        'position',[0.75 0.32 0.24 0.48],...
        'box','off',...
        'xtick',[],'ytick',[]);
    imagesc(vesselImg);
    axis xy image off;
    title('Artery + Vein','fontsize',24);
    hold on;
    if ~xLHidx
        plot(data.centerline.LH.xdata, data.centerline.LH.ydata,'color','w','linewidth',2);
    end
    hold on;
    if ~xRHidx
        plot(data.centerline.RH.xdata, data.centerline.RH.ydata,'color','w','linewidth',2);
    end
    hold off;
        
    % ---- Save figure to file ---- %
    screen2file(strrep(fullpath, '.mat', '_a_v'));
    
    
% % % % % % % %     
% % % % % % % %         %% 2. Plot the fitting results for single pixel, including HRF
% % % % % % % %     monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
% % % % % % % %     w = 800; ht = 600; % Define window size: width(w) and height(ht)
% % % % % % % %     mon = 1; % Only one monitor
% % % % % % % % 
% % % % % % % %     % Determine positioning for main figure
% % % % % % % %     pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
% % % % % % % %         ((monitorPos(mon,4)-ht)/2)+monitorPos(mon,2)/2,...
% % % % % % % %         w ht];
% % % % % % % % 
% % % % % % % %     figure('name','Linear Fit with Exponential Model',...
% % % % % % % %         'position',pos,'color','w');
% % % % % % % % 
% % % % % % % %     % axes for binarized speed data
% % % % % % % %     h1 = axes('units','normalized',...
% % % % % % % %         'position',[0.07 0.93 0.9 0.01],...
% % % % % % % %         'box','off',...
% % % % % % % %         'xtick',[],'ytick',[]);
% % % % % % % %     hold on;
% % % % % % % %     plot(h1, data.time, data.speed.aveBinary,...
% % % % % % % %         'line','none',...
% % % % % % % %         'marker','o',...
% % % % % % % %         'color',[0.5 0.5 0.5],...
% % % % % % % %         'markerfacecolor',[0.5 0.5 0.5]);
% % % % % % % %     set(gca,'XTickLabel','','YTickLabel','');
% % % % % % % % 
% % % % % % % %     % axes for raw speed data
% % % % % % % %     h2 = axes('units','normalized',...
% % % % % % % %         'position',[0.07 0.65 0.9 0.28],...
% % % % % % % %         'box','on');
% % % % % % % %     hold on;
% % % % % % % %     plot(h2, data.time, data.speed.data,'k');
% % % % % % % %     ylabel(h2,'Walk Speed [m/s]');
% % % % % % % %     legend(h2, 'Walking Speed', 'Baseline','location','southeast');
% % % % % % % %     set(gca,'XTickLabel','');
% % % % % % % %     
% % % % % % % %     pix = 1000;
% % % % % % % %     % axes for fractional reflection change [LEFT]
% % % % % % % %     h3 = axes('units','normalized',...
% % % % % % % %         'position',[0.07 0.35 0.9 0.28],...
% % % % % % % %         'box','on');
% % % % % % % %     hold on;
% % % % % % % %     if ~xLHidx
% % % % % % % %         plot(h3, data.time, dRR0filtall{1}(pix,:)*100, 'k');
% % % % % % % %     elseif ~xRHidx
% % % % % % % %         plot(h3, data.time, dRR0filtall{4}(pix,:)*100, 'k');
% % % % % % % %     else
% % % % % % % %         error('No data of LH or RH available');
% % % % % % % %     end
% % % % % % % %     hold on;
% % % % % % % %     if ~xLHidx
% % % % % % % %         xFit_HRFplot = fitResults.LH.xFit_HRF(pix,:)*100;
% % % % % % % %         xFit_HRFaplot = fitResults.LH.xFit_HRFa(pix,:)*100;
% % % % % % % %         xFit_HRFvplot = fitResults.LH.xFit_HRFv(pix,:)*100;
% % % % % % % %     else
% % % % % % % %         xFit_HRFplot = fitResults.RH.xFit_HRF(pix,:)*100;
% % % % % % % %         xFit_HRFaplot = fitResults.RH.xFit_HRFa(pix,:)*100;
% % % % % % % %         xFit_HRFvplot = fitResults.RH.xFit_HRFv(pix,:)*100;
% % % % % % % %     end
% % % % % % % %     plot(h3, data.time, xFit_HRFplot, 'color', 'm','linewidth',2);
% % % % % % % %     hold off;
% % % % % % % %     ylabel(h3,'\DeltaR/R_0 [%]');
% % % % % % % %     legend(h3, 'Raw Data', 'Exponential Fit','location','southeast');
% % % % % % % %     set(gca,'XTickLabel','');
% % % % % % % %     h4 = axes('units','normalized',...
% % % % % % % %         'position',[0.07 0.05 0.9 0.28],...
% % % % % % % %         'box','on');
% % % % % % % %     hold on;
% % % % % % % %     plot(h4, data.time, xFit_HRFaplot,'color', 'r','linewidth',2);
% % % % % % % %     hold on;
% % % % % % % %     plot(h4, data.time, xFit_HRFvplot, 'color', 'b','linewidth',2);
% % % % % % % %     hold off;
% % % % % % % %     xlabel(h4, 'Time [seconds]')
% % % % % % % %     ylabel(h4,'\DeltaR/R_0 [%]');
% % % % % % % %     legend(h4, 'artery','vein','location','southeast');
% % % % % % % %     linkaxes([h1,h2,h3,h4],'x');
% % % % % % % %     
% % % % % % % %     screen2file(strrep(fullpath, '.mat', '_fit'));
end

