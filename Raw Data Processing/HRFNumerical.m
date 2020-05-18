function out = HRFNumerical(s,x,Fr)
    % HRFNumerical -- Use Analytic Methods to compute HRF functions of the observed
    %       changes in hemodynamic signals based on the assumption that the
    %       hemodyanmic response to locomotion is a linear, time-invariant system.
    %
    % USAGE:
    %       h = HRFNumerical(s,x,t);
    %
    % INPUTS:
    %       s: binarized locomotion signal, same length as Dalsa image
    %           data, sample by 1
    %       x: filtered intensity for each pixel in Regions of Interest
    %           (ROI), pixel by samples
    %       Fr: frame rate, Hz, currently we only analyze 3 Hz data
    %
    % OUTPUTS:
    %       HRF: impulse response, pixels by samples
    %       timeShift: shift input data by timeShift seconds to prevent
    %                  boundary artifacts
    %       IRLength: the length of impulse response, in samples
    %
    % Last Modified: 09-27-2016
    % Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)
    %
    % See also:
    %       ImpResp
    %       HRFfeatures

    % Reference:
    % Bingxing Huo et al., Quantitative spearation of arterial and venous
    % cerebral blood volume increases during voluntary locomotion. NeuroImage
    % 105 (2015): 369-379.
    
    % Mean subtract the Intrinsic data to remove constant trend
    x = bsxfun(@minus, x, mean(x,2));    
    
    % Make sure there is no NaN data in binary speed data
    idx = find(isnan(s));
    s(idx) = 0;
    clear idx;
    % make sure s is a column vector
    if size(s,1) < size(s,2) 
        s = s(:);
    end
    
    timeShift = 5; % shift stimulus backward by 10 seconds to avoid boundary effects
    IRLength = 500; % length of impulse response fuction, 1000 samples
    
    %% Compute HRF directly, sample by pixel
    HRF = ImpResp(s,x',Fr,timeShift,IRLength); % HRF contains timeshift data
    HRF = HRF(1:IRLength,:);
    HRF = HRF'; % pixels by samples
    HRF_time = -timeShift:1/Fr:(IRLength-timeShift*Fr-1)/Fr;
    HRF_time = HRF_time(:);
    
    out.timeShift = timeShift;
    out.IRLength = IRLength;
    out.HRF = HRF;
    out.HRF_time = HRF_time;

    %% Get HRF peak response and onset time for each pixel
    [~, ~, HRFpeak, HRFpeaktime] = ...
        HRFfeatures(HRF_time, HRF, Fr, timeShift, IRLength);
    
    out.HRFpeak = HRFpeak;
    out.HRFpeaktime = HRFpeaktime;        
end



function H = ImpResp(datain,dataout,Fr, timeShift, IRLength)
% ImpResp: compute impulse response for each pixel using numerically
%          analysis, this script works for single-input-single-output (SISO) system
%
% INPUTS:
%       datain: stimulus data (input) for the model
%       dataout: response data (output) for the model, sample by pixels
%       Fr: sample rate for the time series
%       timeShift: time shift for the computing model, in seconds
%       IRLength: the length of the impulse response, in samples
%
% OUTPUTS:
%       H: impulse response, samples by pixels
%
% USUAGE:
%       H = ImpResp(dataStimulus, dataResponse, timeShift, IRLength);
%
% Last Modified: 09-27-2016 
% Modified By: Qingguang Zhang (qingguang.zhang@gmail.com)

if nargin < 2
    error('HRFNumerical -> Not enough inputs');
elseif size(datain,1)~=size(dataout,1)
    error('Input and output should be of same lengths.')
end

if nargin < 3
    Fr = 3; % default frame rate, 3 Hz
    timeShift = 10 * Fr; % default time shift, 10 seconds, convert to samples
    IRLength = round((length(datain) - timeShift)/5); % Length of the impulse response
elseif nargin < 4
    timeShift = 10 * Fr; % shift 10 seconds data, convert to samples
    IRLength = round((length(datain) - timeShift)/5); % Length of the impulse response
else
    timeShift = timeShift * Fr;
end

datain = datain(timeShift+1:end); % shift input time backward, assuming non-causal system
dataout=dataout(1:end-timeShift,:); % truncate output to keep same lengths

N = length(datain); % length of processed input data
pixels = size(dataout, 2); % how many pixels in the output data

% Generate a numerical matrix for the input data
c = zeros(2*N-1,1);
c(1:N) = datain; % this is the first column of the Toeplitz Matrix

r = zeros(IRLength,1);
r(1) = datain(1); % this is the first column of the Toeplitz Matrix

TM = toeplitz(c,r); % generate Toeplitz Matrix (TM)

S = [ones(2*N-1,1), TM]; % the numerical matrix for input data, s

% Generate a numerical matrix for the output data
X = zeros(2*N-1, pixels);
X(1:N,:) = dataout(1:N,:);

H = S\X; % this is more acurate and faster than H = inv(S'*S)*S'*X;
end

function [HRF_interp, HRF_time_interp, HRFpeak, HRFpeaktime] = ...
                    HRFfeatures(HRF_time, HRF, Fr, timeShift, IRLength)
% Interpolate HRF
step = 1/300;
HRF_time_interp = -timeShift:step:(IRLength-timeShift*Fr-1)/Fr;
HRF_time_interp = HRF_time_interp(:);
HRF_interp = NaN*zeros(size(HRF,1),length(HRF_time_interp));

for i = 1:size(HRF,1)
    % turn off warning to prevent MATLAB throwing shit
    warning('off','MATLAB:interp1:EmptyMethod');
    HRF_interp(i,:) = interp1(HRF_time, HRF(i,:), HRF_time_interp,'spline');
end

% find peak value of HRF for each pixel within [0 20] second
st = find(HRF_time_interp >= 0, 1, 'first');
fin = find(HRF_time_interp <= 20, 1, 'last');
[~,HRFpeakIdx] = max(abs(HRF_interp(:,st:fin)), [],2); % peak value for each pixel
HRFpeakIdx = HRFpeakIdx + st-1;
HRFpeaktime = HRF_time_interp(HRFpeakIdx);  

HRFpeak = NaN * zeros(length(HRFpeakIdx),1); % initiate matrix
for i = 1:length(HRFpeakIdx)
    HRFpeak(i) = HRF_interp(i,HRFpeakIdx(i)); % find peak response for each pixel
end
end