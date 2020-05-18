function out = PowerLawFit(freq,power)
% function out = PowerLawFit(freq,power,method)
% PowerLawFit: find the best parameter of slope and intercept, so that
% power = slope*freq+intercept
%
% INPUT:
%       freq: raw frequency time series
%       power: raw power spectral density time series
%       method: fitting method for the dataset

% make sure they are column vector
freq = freq(2:end);
power = power(2:end);

% log transform
freqlog = log10(freq); 
powerlog = log10(power);

LSEout = PowerLawFit_LSE(freqlog,powerlog);
LSERout = PowerLawFit_LSE_resamp(freqlog,powerlog);
%MLEout = PowerLawFit_MLE(freq, power);



% figure; loglog(freq, power, 'k'); hold on; 
% loglog(freq, 10.^(LSEout.predict), 'r'); 
% loglog(LSERout.freq, 10.^(LSERout.predict),'b');



out.LSEout = LSEout;
out.LSERout = LSERout;
%out.MLEout = MLEout;

% % We may need the following structure to make this code script more
% % functional
% method = lower(method);
% switch method
%     case 'lse'
%         out = linearfit_lse(freqlog,powerlog);
%     case 'lse resamp'
%         out = linearfit_lse_resamp(freqlog,powerlog);
%     case 'mle'
%     otherwise
%         error('Not valid method');
% end
end


function out = PowerLawFit_LSE(x,y)
% PowerLawFit_LSE: Linear fit using least square estimation
% INPUT:
%       x: log-transformed frequency
%       y: log-transformed power
%
% OUTPUT:
%       out


b = robustfit(x, y, 'ols'); % ordinary least square fit (No weighing)
slope = b(2); % slope
exponent = -slope; % make sure scaling exponent is positive
intercept = b(1); % intercept

y_predict = slope*x+intercept; % predicted log-transformed power
y_residual = y-y_predict; % residual power = raw power - power law fit

% generate output structure
out.slope = slope;
out.exponent = exponent;
out.intercept = intercept;
out.predict = y_predict;
out.residual = y_residual;
end

function out = PowerLawFit_LSE_resamp(x,y)
% PowerLawFit_LSE_resamp: Linear fit using least square estimation
%       We resampled x and y in the frequency domain to avoid the bias
%       towards higher frequency
% INPUT:
%       x: log-transformed frequency
%       y: log-transformed power
%
% OUTPUT:
%       out

% NOTE:
% when linearly-spaced frequency bins are considered under a logarithmic scale, 
% bins in higher-frequencies become progressively more dense, and thus gain 
% disproportionate weight with respect to lower frequency-bins in a subsequent linear
% regression. To avoid this potential bias, we upsampled the PSD curve with 
% logarithmically spaced frequency bins, resulting in equally-spaced frequency bins 
% under logarithmic scale, required to properly estimate the spectral exponent.

% resample x and y to get logarithmically spaced vector
x1 = logspace(x(1), x(end), length(x)); x1 = x1(:); % x1 <=> frequency
logx1 = log10(x1);                                  % logx1 <=> logfreq
logy1 = interp1(x,y,logx1);                         % logy1 <=> logpower
y1 = 10.^(logy1);                                   % y1 <=> power

% Working on the logarithmically spaced vector to get power law fit using
% least square estimation
b = robustfit(logx1,logy1,'ols'); % ordinary least square fit (No weighing)
slope = b(2); % slope
exponent = -slope; % make sure scaling exponent is positive
intercept = b(1); % intercept

logy1_predict = slope*logx1+intercept; % predicted log-transformed power
logy1_residual = logy1-logy1_predict; % residual power = raw power - power law fit

% generate output structure
out.freq_resamp = x1; % for plot purpose
out.power_resamp = y1; % for plot purpose
out.slope = slope;
out.exponent = exponent;
out.intercept = intercept;
out.predict = logy1_predict;
out.residual = logy1_residual;

end

% function out = PowerLawFit_MLE(x,y)
% % PowerLawFit_MLE: fit a power-law distribution using Maximum-Likelihood
% % Estimation (MLE)
% %
% % INPUT:
% %       x: raw frequency
% %       y: raw power
% 
% % See also:
% %       plfit.m
% 
% % preset parameters for the power law fit using MLE
% range = 1.01:0.01:3.50;
% 
% if length(y) <= 100
%     [alpha, xmin, loglik] = plfit(y, 'range',range, 'finite');
%     % small-size correction to the power-law exponent
% else
%     [alpha, xmin, loglik] = plfit(y,'range',range);
% end
% 
% out.exponent = alpha; % power-law exponent
% out.powercut = xmin; % minumum spectral fit of power law distribution
% out.PLawLikelihood = loglik; % loglikelihood of the power-law distribution
% 
% end