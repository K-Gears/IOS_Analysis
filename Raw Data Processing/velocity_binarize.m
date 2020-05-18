   %% velocity_binarize.m
% Bingxing Huo, March 2015
% Based on velocity_proc.m
% This script reads and binarizes the velocity data. Also downsample it to match the
% frequency of Ft.
% #### This code calls the following function:
% ReadInAcquireTDMS.m
%
% Input:
%   velocitydata: read or define raw velocity data
%       - to read: a string containing the file name (example: '
%       140318_11_43_2018') (note: do not use suffix)
%       - to define: directly input the vector containing velocity data
%   Fs: acquisition frequency
%   Ft: Optional. target frequency. Default = Fs.
%   acc_cutoff: Optional. acceleration threshold. Default = 1e-6.
% Output:
%   - imp_bin: binarized locomotion (impulses)
%   - velocity: raw velocity data
% Example:
% [imp_bin, velocity]=velocity_binarize(' 140318_11_43_2018');
%   - it will read the velocity data contained in the tdms file '
%   140318_11_43_2018.tdms'. And gives the raw velocity data as well as the
%   binarized locomotion data. 
%   - the binarized locomotion is at the default 30 kHz.
% [imp_bin, ~]=velocity_binarize(velocity,'',10);
%   - since there is already velocity data, you may skip this output by
%   using '~' in the output argument. 
%   - the binarized locomotion is at 10 Hz
%%
function [imp_bin,velocity]=velocity_binarize(velocitydata,Fs,Ft,acc_cutoff)
% 1. read or define velocity data
if ischar(velocitydata)==1 % if the input is a file name
    % read velocity data
    volfile=[fileid '.tdms'];
    volt=ReadInAcquireTDMS(volfile);
    %convert voltage -> velocity
    maxvel=2*pi*0.06*10;
    maxVol=10; % maximum voltage output (read from USDigital device)
    rate=maxvel/maxVol;
    % take averge velocity for every frame
    velocity=-volt.Data(1,:)*rate; % velocity data is recorded in analog channel 1
elseif isnumeric(velocitydata)==1 % if the input is the velocity data itself
    % define velocity data
        maxvel=2*pi*0.06*10;
    maxVol=10; % maximum voltage output (read from USDigital device)
    rate=maxvel/maxVol;
     velocity=velocitydata*rate;
end
% 2. assume or define sampling frequency
if nargin<2
    Fs=30000; % default = 30 kHz 
elseif isempty(Fs)==1 
    Fs=30000;
end
% 3. low-pass filter velocity below 10 Hz
Flp=10;
% [zeroa,poleb,gain]=butter(5,2/Fs*[Fhp,Flp],'bandpass');
[zeroa,poleb,gain]=butter(5,(2/Fs)*Flp,'low');
[sos,g]=zp2sos(zeroa,poleb,gain);
% h2=dfilt.df2sos(sos,g);
velfilt=filtfilt(sos,g,velocity);
% 3. define target frequency
if nargin<3
    Ft=Fs; % default target frequency is the same as the acquisition frequency
elseif isempty(Ft)==1
    Ft=Fs;
end
% 4. calculate acceleration
acc=diff(velfilt); % acceleration
N=length(acc);
Flight=Fs/Ft;
L=floor(N/Flight);
Index=floor(Flight);
% 5. binarize
if nargin<4 % if user specified the acceleration threshold
    acc_cutoff=1e-6;
    %Use ~1.5e-7 for under p12, 1e-6 under p15, 1e-5 p20, 1e-4 rest of
    %animals 4-30-16 KG
end
accbin=(abs(acc)>acc_cutoff);
% 6. average within each window
if  Ft==Fs
    imp_bin=accbin;
    imp_ave=imp_bin;
else
    imp_bin=zeros(L,1);
    for K=1:L
        imp_bin(K)=(mean(accbin(((K-1)*Index)+1:K*Index))>.1);
    end
    imp_ave=zeros(L,1);
    for K=floor(Ft/2)+1:L-floor(Ft/2)-1
        imp_ave(K)=sum(imp_bin(K-floor(Ft/2):(K+floor(Ft/2)+1)));
    end
end
% % 7. visualize
% figure
% bar([1:L]/Ft,imp_ave)
% axis tight
% set(gca,'fontsize',14)
% xlabel('time (second) ','fontsize',18)
% ylabel('thresholded acceleration (bin=1sec) ','fontsize',14)
