% matchlength.m
% by Bingxing Huo, 2011
% This script matches the high frequency data to the length of the low
% frequency data
% Input:
%   - rawdata: the high frequency data
%   - L: the total length of the low frequency data
%   - rate: the theoretical frequency mutiples Fhigh/Flow
function dataave=matchlength(rawdata,L,rate)
[Ndata,Nchannel]=size(rawdata);
if Ndata<Nchannel
    rawdata=rawdata';
    [Ndata,Nchannel]=size(rawdata);
end
s=floor(Ndata/rate); % theoretical length of L
% datare=zeros(L*rate,Nchannel);
for c=1:Nchannel
datare(:,c)=resample(rawdata(:,c),L,s); % resample to match the real length of L
end
dataave=zeros(L,Nchannel);
for i=1:L
    dataave(i,:)=mean(datare((i-1)*rate+1:i*rate,:)); % average within each frame
end