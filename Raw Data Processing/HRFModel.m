function [a, v, c, HRFa, HRFv, HRF, xFit_HRFa, xFit_HRFv, xFit_HRF, cc] = HRFModel(s,x,t)
    % HRFModel -- Construct a linear convolution model to fit the observed
    %       changes in hemodynamic signals based on the assumption that the
    %       hemodyanmic response to locomotion is a linear, time-invariant system.
    %
    % USAGE:
    %       h = HRFModel(s,x,t);
    %
    % INPUTS:
    %       s: binarized locomotion signal, same length as Dalsa image data
    %       x: filtered intensity for each pixel in Regions of Interest (ROI)
    %       t: time, in seconds
    %
    % OUTPUTS:
    %       a: pixel specific weights for arteries
    %       v: pixel specific weights for veins
    %       c: constant offset associated with each pixel
    %       HRFa: impulse response using only arterial component
    %       HRFv: impulse response using only venous component
    %       HRF: impulse response using both arterial and venous component
    %       xFit_HRFa: Fitted data using only arterial component
    %       xFit_HRFv: Fitted data using only venous component
    %       xFit_HRF: Fitted data using both arterial and venous components
    %       cc: Pearson correlation for each pixel, numPixels by 1 matrix
    %
    % Last Modified: 07-23-2015
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)
    %
    % See also:
    %       ModelFitFcn

    % Reference:
    % Bingxing Huo et al., Quantitative spearation of arterial and venous
    % cerebral blood volume increases during voluntary locomotion. NeuroImage
    % 105 (2015): 369-379.
    
    % NOTE:
    %       DeltaR/R0 = s*h+c+e;
    %       h = a*exp(-t/tauA)+v*exp(-t/tauV);
    
    % Make sure there is no NaN data in binary speed data
    idx = find(isnan(s));
    s(idx) = 0;
    clear idx;
    
    if size(s,1) < size(s,2) % make sure s is a column vector
        s = s(:);
    end
    
    % Initialize option parameters for fminsearch
    options = optimset('Display','off',...      % display no output
                       'FunValCheck','off',...  % check objective values
                       'MaxFunEvals', 1000,...  % max number of function evaluations allowed
                       'MaxIter', 1000,...      % max number of iteration allowed
                       'TolFun',1e-8,...        % termination tolerance on the function value
                       'TolX',1e-8,...          % termination tolerance on x
                       'UseParallel','always'); % always use parallel computation
                   
% %     % Initialize starting parameters
% %     numPixels = size(x,1);
% %     c0 = x(:,1); % constant offset associated with each pixel in first frame
% %     x0 = [-0.1*ones(numPixels,2), c0];                   
% %                    
% %     % Pre-allocte memory for fit results
% %     fitParameters=zeros(numPixels,3); % numPixel by 3 matrix, [scaleA, scaleV, dc]
% % 
% %     parfor i = 1:numPixels
% %         fitParameters(i,:)=fminsearch(@(xx) ModelFitFcn(s,x,xx,i,t),x0(i,:),options);
% %     end

    % -- Based on Bing's paper, all the starting parameter are set to 0.1    
    % Initialize starting parameters
        y0 = [0.1 0.1 0.1];                                   
    % Pre-allocte memory for fit results
    numPixels = size(x,1); % number of pixels
    fitParameters=zeros(numPixels,3); % numPixels by 3 matrix, [a, v, c]

    disp('>>> Start fitting data using least mean square method');
    tic
    parfor i = 1:numPixels
        fitParameters(i,:)=fminsearch(@(y) ModelFitFcn(s,x,y,i,t),y0,options);
    end
    toc
    disp('>>> Finish fitting procedure');
    
    % Get all the weight parameters
    a = fitParameters(:,1); % Pixel specific weights for arteries
    v = fitParameters(:,2); % Pixel specific weights for veins
    c = fitParameters(:,3); % Pixel specific constant offsets
    
    disp('>>> Compute HRF and Fitted data');
    tic
    % Compute HRF for each Pixel (HRFa, HRFv, HRF are numPixels by Time matrix)
    HRFa = a*exp(-t/4); % HRF using arterial component
    HRFv = v*exp(-t/40); % HRF using venous component
    HRF = HRFa + HRFv; % HRF using both arterial and venous component
    
    % Compute fitted data (the fitted data are numPixels by length(s) matrix)
    xFit_HRFa = conv2(s', HRFa); % Fitted data using arterial HRF
    xFit_HRFv = conv2(s', HRFv); % Fitted data using venous HRF
    xFit_HRF = bsxfun(@plus, conv2(s',HRF), c); % Fitted data using both arterial and venous HRF
    
    % Truncate matrix to fit the length of x
    N = length(s);
    xFit_HRFa = xFit_HRFa(:,1:N);
    xFit_HRFv = xFit_HRFv(:,1:N);
    xFit_HRF = xFit_HRF(:,1:N);
    
    cc = corr(xFit_HRF', x'); % compute Pearson correlation for each pixel
    cc=diag(cc); % get the values on diagonal line, numPixels by 1 vector
    toc
    disp('>>> Finish Computing HRF');
end


function e = ModelFitFcn(s,x,y,Pixidx,t)
    % ModelFitFcn -- This is the function need to be evaluated in fminsearch
    %
    % INPUTS:
    %       s: binarized locomotion signal
    %       x: IOS signal
    %       y: fitting parameters, [a,v,c]
    %           a: weight parameter for arterial component
    %           v: weight parameter for venous component
    %           c: constant offset
    %       Pixidx: the pixel number
    %
    % OUTPUTS:
    %       e: error term (NOTE: the goal of fminsearch is to minimize e)
    %
    % Last Modified: 07-20-2015 
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)

    a = y(1); % pixel specific weights for arteries
    v = y(2); % pixel specific weights for veins
    c = y(3); % constant offset associated with each pixel
    
    tauA = 4; % fast decay time constant, fixed at 4 s
    tauV = 40; % slow decay time constant, fixed at 40 s
    
    N = length(s);
    
    h = a*exp(-t/tauA)+v*exp(-t/tauV);
    
    xFit = conv(s, h) + c; % The length of fit is 2*N-1
    
    e = sqrt(sum((xFit(1:N)-x(Pixidx,1:N)).^2));
    
end
