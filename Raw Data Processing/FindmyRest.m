function [RestInds]=FindmyRest(binarized_running,Whisker_Stim_times,Pad_time,Seg_duration,Fs)
Runtimes=find(binarized_running);
StimTypes=fieldnames(Whisker_Stim_times);
StimEvents(1:length(binarized_running))=0;
for k=1:numel(StimTypes)
    StimEvents(Whisker_Stim_times.(StimTypes{k})==1)=1;
end
Stimtimes=find(StimEvents);
BehavingTimes=union(Runtimes,Stimtimes);
Pad_frames=Pad_time*Fs;
Seg_frames=Seg_duration*Fs;
RestPoints=diff(BehavingTimes,1,2);
RealRest=find(RestPoints>(Pad_frames+Seg_frames));
Restframes=RestPoints(RealRest);
for k=1:numel(Restframes)
    RestInds(1,k)=Restframes(k)+Pad_frames;
    RestInds(2,k)=Restframes(k)+Pad_frames+Seg_frames;
end