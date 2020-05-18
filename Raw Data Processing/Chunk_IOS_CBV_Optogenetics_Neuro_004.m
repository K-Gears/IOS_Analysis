function Chunk_IOS_CBV_Optogenetics_Neuro_004(filename,theLED,bandFilter,cutFilter,animalAge,varargin)
%Function to determine time frame around whiskerpuff for CBV change
%Solenoid identification changed for use w/Neotate Whisker Trial v5+
%Written by Kyle Gheres, Oct 2014
% Identify Files
if nargin==0
    filename = uigetfile('*rawdata.mat','MultiSelect', 'on');
end
% Control for single file instead of a list
if iscell(filename) == 0
    filename = {filename};
end

%% Constants
if ~exist('theLED','var')
    theLED=input('Enter model of LED used for illumination e.g. M530L3. ','s');
end
if ~exist('bandFilter','var')
    bandFilter=input('Enter model of optical bandpass filter in light path. e.g FB570-10. ','s');
end
if ~exist('cutFilter','var')
    cutFilter=input('Enter model of optical cut off filter in light path. e.g FEL0500. ','s');
end
createmptystruct=1;
Run_State={'Still','Run'};
Run_Type={'Volitional','Evoked'};
Run_Length={'Twitch','Short','Medium','Long','Very_Long','Extra_Long'};
HemoResp={'Positive','Negative'};
Run_Cut_Off=[0.5 2 5 10 15];
T_seg=5;
T_beg=10;
Acc_Thresh=5e-5;
Lead_Time=5;%Time to preceed whisker puff
Follow_Time=15;% Time to follow Whisker Puff
Norm_Time=2;%Time after Lead Time to use as normalization window
RunEventThresh=3;%Time in seconds to bin for shortest running events
order=3;
ftype='low';
params.tapers=[3 5];
params.trialave= 1;
params.fpass=[0 100];
movingwin=[1 0.05];
stimcount=1;
convert2uM=1e6;
EMGthresh=0.15;

dateData=char(datetime('now'));
temp_a=strrep(dateData,':','');
temp_b=strrep(temp_a,' ','_');
saveDate=strrep(temp_b,'-','_');


%% Create Empty Locomotion structure
for q=1:length(Run_Type)
    EMG_Data.(Run_Type{q}).All.IOS.forepaw.ChunkRefl=[];
    Loco_Data.(Run_Type{q}).All.IOS.forepaw.ChunkRefl=[];
end



%% Get constants for converting reflectance to Hb concentration
[weightedcoeffHbO,weightedcoeffHbR,weightedcoeffHbT]=getHbcoeffs(theLED,bandFilter,cutFilter);

%% Get resting LFP data 
% eventCount=1;
% for fil=1:length(filename)
%     load(filename{fil});
%     ExpectedLength=RawData.an_fs*300;
%     Ball_Vel=RawData.vBall(1:ExpectedLength);
%     [imp_bin]=velocity_binarize(Ball_Vel,RawData.an_fs,RawData.an_fs,1e-5);
%     [~,T_stand]=motion_cont_Puff(imp_bin,RawData.an_fs,T_seg,T_beg);
%     for eventNum=1:size(T_stand,2)
%         EventLength=(T_stand(2,eventNum)-T_stand(1,eventNum))/RawData.an_fs;
%         if EventLength>=15
%             RefLFP{eventCount}=RawData.Neuro(T_stand(1,eventNum):T_stand(2,eventNum));
%             WorkLFP(:,eventCount)=RawData.Neuro((T_stand(1,eventNum)+2.5*RawData.an_fs):(T_stand(1,eventNum)+12.5*RawData.an_fs));
%             eventCount=eventCount+1;
%         end
%     end
% end
% params.Fs=RawData.an_fs;
% [S,f]=mtspectrumc(diff(WorkLFP,1,1),params);
% RestRefLFP=S;
% RestRefFreq=f;

%% Process Data
for fil = 1:length(filename)
    close all;
    indfile = filename{fil};
    dashLocs=strfind(indfile,'_');
    animal=indfile(1:(dashLocs(2)-1));%(1:9);
    date=indfile((dashLocs(3)+1):(dashLocs(4)-1));%(10:15);(14:19);
    hem=indfile((dashLocs(2)+1):(dashLocs(3)-1));%(7:8);%(11:12);
    filenm= indfile((dashLocs(4)+1):(dashLocs(end)-1));%(10:26);%(14:30);
    load(indfile);
    Loco_Data.filenames{fil}=indfile;
    if fil==1
        [SharVars]=GetSharVars_2(animal,hem);
    end
    %% File Specific Constants
    ROI_name=fieldnames(RawData.IOS);
    reOrderROI=strcmpi(ROI_name,'Pixelwise')==1;
    ROI_name{reOrderROI}=ROI_name{length(ROI_name)};
    ROI_name{length(ROI_name)}='Pixelwise';
    HbCoeffTypes={'HbO','HbR'};
    orderROI=1:length(ROI_name);
    orderROI(reOrderROI)=length(ROI_name);
    orderROI(length(ROI_name))=find(reOrderROI==1);
    Wn= 1/(RawData.dal_fr*0.5); %Bandpass filter properties [low freq, high freq]/nyquist freq
    [zeroa,poleb,gain]=butter(order,Wn,ftype);
    [sos,g]=zp2sos(zeroa,poleb,gain);
    [ball_b,ball_a]=butter(3,30/(RawData.an_fs*0.5));
    params.Fs=RawData.dal_fr;
    
    %% Binarize locomotion events
    if exist('T_run','var')==1
        clear T_run;
    end
    [imp_bin,velocity]=velocity_binarize(RawData.vBall,RawData.an_fs,RawData.dal_fr,Acc_Thresh);
    if max(imp_bin)==0
        RunEvents(fil,(1:size(imp_bin,1)))=0;
        T_run=[];
        new_T_run=[];
        fprintf(1,'Animal did not run during trial\n')
    else
    [T_run,~,new_T_run,run_frac]=motion_cont_Puff(imp_bin,RawData.dal_fr,T_seg,T_beg);
    RunEvents(fil,(1:size(imp_bin,1)))=0;
    if isempty(T_run)==0
        for k=1:size(T_run,2)
            Run_dur=T_run(2,k)-T_run(1,k);
            if gt(Run_dur,(2*RawData.dal_fr))
                RunEvents(fil,T_run(1,k):T_run(2,k))=1;
            end
        end
    else
        fprintf(1,'Animal did not run during trial\n')
    end
    end
    velocity=filtfilt(ball_b,ball_a,velocity);
    RunInds=find(RunEvents(fil,:)==1);
    %Puff_Data.Ball.Velociy(fil,:)=velocity(1:(5*60*RawData.an_fs));
    %Puff_Data.Ball.Run_Frac(fil)=run_frac;
    Puff_Data.dal_fr=RawData.dal_fr;
    Puff_Data.animalAge=animalAge;
    %% Find when animal was puffed
    if ~isfield(RawData,'LED')
        RawData.LED(1:length(RawData.Sol))=0;
    end
    stimchk=union(RawData.Sol,round(RawData.LED,0));
    if max(stimchk)>0
        Sol_sample=downsample(RawData.Sol,floor(RawData.an_fs/RawData.dal_fr));
        Round_Sample=round(RawData.LED,0);% rounding to nearest integer to eliminate any bounce in the signal for binarizing of behavior
        TTL_Change=diff(Round_Sample);%Find where TTL trigger changes state
        if max(TTL_Change)==0
            Laser_State_Change=[];
            Stim_Count=0;
            fprintf('No Optostimuli this trial\n')
        else
            Laser_State_Change=find(TTL_Change==max(TTL_Change));
            Laser_Win=RawData.AcquistionParams.Laser_Duration*RawData.an_fs; %time laser pulses after on set
            Stim_On=Laser_State_Change(1);
            Stim_Off=Stim_On+Laser_Win;
            Stim_Count=numel(Laser_State_Change);
        end
        Counter=1;
        
        while Counter<=Stim_Count % for the purpose of windowing IOS data only keep the first on signal of each stimulus train
            Low_Bound=find(Laser_State_Change>Stim_On);
            Upper_Bound=find(Laser_State_Change<Stim_Off);
            The_Stim_Win=intersect(Low_Bound,Upper_Bound);
            Laser_State_Change(The_Stim_Win)=[];
            if (Counter+1)<=numel(Laser_State_Change)
                Stim_On=Laser_State_Change(Counter+1);
                Stim_Off=Stim_On+Laser_Win;
            end
            Stim_Count=numel(Laser_State_Change);
            Counter=Counter+1;
        end
        Laser_State_Change=round(((Laser_State_Change/RawData.an_fs)*RawData.dal_fr),0);
        Puff.Laser_Stim(1:size(RawData.IOS.(ROI_name{1}).CBVrefl,2))=0;
        Puff.Laser_Stim(Laser_State_Change)=1;
        Puff_Event=Sol_sample(1:size(RawData.IOS.(ROI_name{1}).CBVrefl,2));
        Puff.Contra_Puff=(Puff_Event==1);%Left Whisker Pad when window is on right
        Puff.Ipsi_Puff=(Puff_Event==3);
        Puff.Control_Puff=(Puff_Event==2);
        LaserLogical(1:length(Puff.Laser_Stim))=0;
        for stimnum=1:length(Laser_State_Change)
            LaserLogical(Laser_State_Change(stimnum):(Laser_State_Change(stimnum)+RawData.AcquistionParams.Laser_Duration*RawData.dal_fr))=1;
        end
        LaserLogical=logical(LaserLogical);
        StimTimes=vertcat(Puff.Contra_Puff,Puff.Ipsi_Puff,Puff.Control_Puff);
        StimTimes=max(StimTimes,[],1);
        StimPoints=find(StimTimes==1);
        LaserPoints=find(LaserLogical==1);
    else
        fprintf('No stimulus during trial\n')
        Puff=[];
        LaserPoints=0;
    end
    %% Binarized EMG
    eventDur=5*RawData.dal_fr; %How long an event must be in seconds to be categorized
    restDur=10*RawData.dal_fr; % Minumum time between events in seconds
    
    [binarized_EMG,perc_aboveThresh]=binarizeEMG(filenm,RawData.MUA,RawData.an_fs,EMGthresh,RawData.dal_fr);
    [EMGinds,EventLengths]=getEMGevents(binarized_EMG,eventDur,restDur);
    %% Acquire reflectance data for specific timeframe
     Plot_ROI_CBV_Trial_Optogenetics_005(filenm,RawData,imp_bin,T_run,new_T_run,velocity,weightedcoeffHbO,weightedcoeffHbR,weightedcoeffHbT)
    
    %% Power Spectra of neural activity
    movingwin=[3 0.1];
    params.Fs=RawData.an_fs;
    params.fpass=[1 100];
    params.tapers=[5 9];
    
    [S,t,f]=mtspecgramc_GPU(RawData.LFP,movingwin,params);
    AvgPwr=mean(S,1);
    AvgMat=repmat(AvgPwr,size(S,1),1);
    NormLFPPwr=((S-AvgMat)./AvgMat)*100;
    
        %% EMG Triggered Behavior
        if ~isempty(EMGinds)
            fprintf(1,'chunking reflectance around running\n')
            for k=1:size(EMGinds,2)
                if (EMGinds(2,k)/RawData.dal_fr)<t(length(t)) && ((EMGinds(1,k)+(Follow_Time*RawData.dal_fr))/RawData.dal_fr)<t(length(t))
                    if (EMGinds(1,k)-(Lead_Time*RawData.dal_fr))>=(movingwin(1)*RawData.dal_fr)
                        if (EMGinds(1,k)+(Follow_Time*RawData.dal_fr))<=length(RawData.IOS.barrels.CBVrefl)
                            StimWin=(EMGinds(1,k)-(Lead_Time*RawData.dal_fr)):(EMGinds(1,k)+(Follow_Time*RawData.dal_fr));
                            if max(ismember(LaserPoints,StimWin))==0
                                if isempty(Puff)
                                    TriggerLabel=Run_Type{1};
                                elseif max(ismember(StimPoints,StimWin))==1
                                    TriggerLabel=Run_Type{2};
                                else
                                    TriggerLabel=Run_Type{1};
                                end
                                EventLength=(EMGinds(2,k)-EMGinds(1,k))/RawData.dal_fr;
                                if EventLength< Run_Cut_Off(2)
                                    RunLabel=Run_Length{2};
                                elseif EventLength>Run_Cut_Off(2) && EventLength<Run_Cut_Off(3)
                                    RunLabel=Run_Length{3};
                                elseif EventLength>Run_Cut_Off(3) && EventLength<Run_Cut_Off(4)
                                    RunLabel=Run_Length{4};
                                elseif EventLength>Run_Cut_Off(4) && EventLength<Run_Cut_Off(5)
                                    RunLabel=Run_Length{5};
                                elseif EventLength>Run_Cut_Off(5)
                                    RunLabel=Run_Length{6};
                                end
                                RunInds=(EMGinds(1,k)-(Lead_Time*RawData.dal_fr)):(EMGinds(1,k)+(Follow_Time*RawData.dal_fr));
                                NeuroInds=round((RunInds(1)/RawData.dal_fr)*RawData.an_fs,0):1:round((RunInds(end)/RawData.dal_fr)*RawData.an_fs,0);
                                SpecInds=find(t<(RunInds(1)/RawData.dal_fr),1,'last'):find(t>(RunInds(end)/RawData.dal_fr),1,'first');
                                EvntInds=find(t<(EMGinds(1,k)/RawData.dal_fr),1,'last'):find(t>(EMGinds(2,k)/RawData.dal_fr),1,'first');
                                NormSpec=find(t<(RunInds(1)/RawData.dal_fr),1,'last'):find(t>((EMGinds(1,k)-RawData.dal_fr)/RawData.dal_fr),1,'first');
                                
                                if isempty(EMG_Data.(TriggerLabel).All.IOS.forepaw.ChunkRefl)
                                    AllCount=1;
                                else
                                    AllCount=size(EMG_Data.(TriggerLabel).All.IOS.forepaw.ChunkRefl,1)+1;
                                end
                                for ROI_num=1:length(ROI_name)
                                    ChunkRefl=filtfilt(sos,g,RawData.IOS.(ROI_name{ROI_num}).CBVrefl(:,RunInds));
                                    EventRefl=RawData.IOS.(ROI_name{ROI_num}).CBVrefl(:,(EMGinds(1,k):EMGinds(2,k)));
                                    if strcmpi(ROI_name{ROI_num},'Pixelwise')
%                                         NormalizationHbt=mean(ChunkRefl(:,(1:(Norm_Time*RawData.dal_fr))),2);
%                                         NormChunk=repmat(NormalizationHbt,1,size(ChunkRefl,2));
%                                         NormEvent=repmat(NormalizationHbt,1,size(EventRefl,2));
%                                         ChunkHbT=(weightedcoeffHbT.*log(ChunkRefl./NormChunk)).*convert2uM;
%                                         EventHbT=(weightedcoeffHbT.*log(EventRefl./NormEvent)).*convert2uM;
%                                         EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkRefl(:,:,AllCount)=ChunkRefl;
%                                         EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkHbT(:,:,AllCount)=ChunkHbT;
%                                         EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedChunkHbT(:,:,AllCount)=ChunkHbT-repmat(ChunkHbT(:,(Lead_Time*RawData.dal_fr)),1,size(ChunkHbT,2));
%                                         EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventRefl{AllCount}=EventRefl;
%                                         EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventHbT{AllCount}=EventHbT;
%                                         EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedEventHbT{AllCount}=EventHbT-repmat(EventHbT(:,(1)),1,size(EventHbT,2));
                                    else
                                        NormalizationHbt=mean(ChunkRefl(1:(Norm_Time*RawData.dal_fr)),2);
                                        ChunkHbT=(weightedcoeffHbT.*log(ChunkRefl/NormalizationHbt)).*convert2uM;
                                        EventHbT=(weightedcoeffHbT.*log(EventRefl/NormalizationHbt)).*convert2uM;
                                        EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkRefl(AllCount,:)=ChunkRefl;
                                        EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkHbT(AllCount,:)=ChunkHbT;
                                        EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedChunkHbT(AllCount,:)=ChunkHbT-ChunkHbT(:,(Lead_Time*RawData.dal_fr));
                                        EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventRefl{AllCount}=EventRefl;
                                        EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventHbT{AllCount}=EventHbT;
                                        EMG_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedEventHbT{AllCount}=EventHbT-EventHbT(:,(1));
                                    end
                                end
                                if length(EMG_Data.(TriggerLabel).All.IOS.forepaw.ZeroedEventHbT{AllCount})<Follow_Time*RawData.dal_fr
                                    EMG_Data.(TriggerLabel).All.HemoRespVolTrial(AllCount)=sum(EMG_Data.(TriggerLabel).All.IOS.forepaw.ZeroedChunkHbT(AllCount,((Lead_Time*RawData.dal_fr):end)))/length(EMG_Data.(TriggerLabel).All.IOS.forepaw.ZeroedChunkHbT(AllCount,((Lead_Time*RawData.dal_fr):end)));
                                else
                                    EMG_Data.(TriggerLabel).All.HemoRespVolTrial(AllCount)=sum(EMG_Data.(TriggerLabel).All.IOS.forepaw.ZeroedEventHbT{AllCount})/length(EMG_Data.(TriggerLabel).All.IOS.forepaw.ZeroedEventHbT{AllCount});
                                end
                                EMG_Data.(TriggerLabel).All.Neuro.ChunkRawLFP(AllCount,:)=RawData.LFP(NeuroInds);
                                EMG_Data.(TriggerLabel).All.Neuro.ChunkSpectrogram(:,:,AllCount)=S(SpecInds,:);
                                EMG_Data.(TriggerLabel).All.Neuro.NormChunkSpec(:,:,AllCount)=((S(SpecInds,:)-repmat(mean(S(NormSpec,:),1),size(S(SpecInds,:),1),1))./...
                                    repmat(mean(S(NormSpec,:),1),size(S(SpecInds,:),1),1))*100;
                                EMG_Data.(TriggerLabel).All.EMG.ChunkRawEMG(AllCount,:)=RawData.MUA(NeuroInds);
                                EMG_Data.(TriggerLabel).All.Ball.ChunkRawBall(AllCount,:)=RawData.vBall(NeuroInds);
                                EMG_Data.(TriggerLabel).All.Neuro.EventSpectrogram{AllCount}=S(EvntInds,:);
                                EMG_Data.(TriggerLabel).All.Neuro.NormEventSpec{AllCount}=((S(EvntInds,:)-repmat(mean(S(NormSpec,:),1),size(S(EvntInds,:),1),1))./...
                                    repmat(mean(S(NormSpec,:),1),size(S(EvntInds,:),1),1))*100;
                                if EMGinds(2,k)==9001
                                    EMG_Data.(TriggerLabel).All.Neuro.EventRawLFP{AllCount}=RawData.LFP(round(((EMGinds(1,k)/RawData.dal_fr)*RawData.an_fs),0):end);
                                    EMG_Data.(TriggerLabel).All.EMG.EventRawEMG{AllCount}=RawData.MUA(round(((EMGinds(1,k)/RawData.dal_fr)*RawData.an_fs),0):end);
                                    EMG_Data.(TriggerLabel).All.Ball.EventRawBall{AllCount}=RawData.vBall(round(((EMGinds(1,k)/RawData.dal_fr)*RawData.an_fs),0):end);
                                else
                                    EMG_Data.(TriggerLabel).All.Neuro.EventRawLFP{AllCount}=RawData.LFP(round((((EMGinds(1,k):EMGinds(2,k))/RawData.dal_fr)*RawData.an_fs),0));
                                    EMG_Data.(TriggerLabel).All.EMG.EventRawEMG{AllCount}=RawData.MUA(round((((EMGinds(1,k):EMGinds(2,k))/RawData.dal_fr)*RawData.an_fs),0));
                                    EMG_Data.(TriggerLabel).All.Ball.EventRawBall{AllCount}=RawData.vBall(round((((EMGinds(1,k):EMGinds(2,k))/RawData.dal_fr)*RawData.an_fs),0));
                                end
                                EMG_Data.(TriggerLabel).All.f=f;
                                EMG_Data.(TriggerLabel).All.t=((1:size(EMG_Data.(TriggerLabel).All.Neuro.NormChunkSpec(:,:,AllCount),1))*movingwin(2))-(movingwin(2)+Lead_Time);
                                EMG_Data.(TriggerLabel).All.DurationLabel{AllCount}=RunLabel;
                                EMG_Data.(TriggerLabel).All.fileName{AllCount}=filenm;
                                EMG_Data.(TriggerLabel).All.NegFlag(AllCount)=EMG_Data.(TriggerLabel).All.HemoRespVolTrial(AllCount)>0;
                            else
                                fprintf('Optostim during event. Skipping\n')
                            end
                        end
                    end
                end
            end
        end
    
    %% Run Triggered Behavior
    if ~isempty(new_T_run)
        fprintf(1,'chunking reflectance around running\n')
        for k=1:size(new_T_run,2)
            if (new_T_run(2,k)/RawData.dal_fr)<t(length(t)) && ((new_T_run(1,k)+(Follow_Time*RawData.dal_fr))/RawData.dal_fr)<t(length(t))
                if (new_T_run(1,k)-(Lead_Time*RawData.dal_fr))>=(movingwin(1)*RawData.dal_fr)
                    if (new_T_run(1,k)+(Follow_Time*RawData.dal_fr))<=length(RawData.IOS.barrels.CBVrefl)
                        StimWin=(new_T_run(1,k)-(Lead_Time*RawData.dal_fr)):(new_T_run(1,k)+(Follow_Time*RawData.dal_fr));
                        if max(ismember(LaserPoints,StimWin))==0
                            if isempty(Puff)
                                TriggerLabel=Run_Type{1};
                            elseif max(ismember(StimPoints,StimWin))==1
                                TriggerLabel=Run_Type{2};
                            else
                                TriggerLabel=Run_Type{1};
                            end
                            EventLength=(new_T_run(2,k)-new_T_run(1,k))/RawData.dal_fr;
                            if EventLength< Run_Cut_Off(2)
                                RunLabel=Run_Length{2};
                            elseif EventLength>Run_Cut_Off(2) && EventLength<Run_Cut_Off(3)
                                RunLabel=Run_Length{3};
                            elseif EventLength>Run_Cut_Off(3) && EventLength<Run_Cut_Off(4)
                                RunLabel=Run_Length{4};
                            elseif EventLength>Run_Cut_Off(4) && EventLength<Run_Cut_Off(5)
                                RunLabel=Run_Length{5};
                            elseif EventLength>Run_Cut_Off(5)
                                RunLabel=Run_Length{6};
                            end
                            RunInds=(new_T_run(1,k)-(Lead_Time*RawData.dal_fr)):(new_T_run(1,k)+(Follow_Time*RawData.dal_fr));
                            NeuroInds=round((RunInds(1)/RawData.dal_fr)*RawData.an_fs,0):1:round((RunInds(end)/RawData.dal_fr)*RawData.an_fs,0);
                            SpecInds=find(t<(RunInds(1)/RawData.dal_fr),1,'last'):find(t>(RunInds(end)/RawData.dal_fr),1,'first');
                            EvntInds=find(t<(new_T_run(1,k)/RawData.dal_fr),1,'last'):find(t>(new_T_run(2,k)/RawData.dal_fr),1,'first');
                            NormSpec=find(t<(RunInds(1)/RawData.dal_fr),1,'last'):find(t>((new_T_run(1,k)-RawData.dal_fr)/RawData.dal_fr),1,'first');
                            if isempty(Loco_Data.(TriggerLabel).All.IOS.forepaw.ChunkRefl)
                                AllCount=1;
                            else
                                AllCount=size(Loco_Data.(TriggerLabel).All.IOS.forepaw.ChunkRefl,1)+1;
                            end
                            for ROI_num=1:length(ROI_name)
                                ChunkRefl=filtfilt(sos,g,RawData.IOS.(ROI_name{ROI_num}).CBVrefl(:,RunInds));
                                EventRefl=RawData.IOS.(ROI_name{ROI_num}).CBVrefl(:,(new_T_run(1,k):new_T_run(2,k)));
                                if strcmpi(ROI_name{ROI_num},'Pixelwise')
%                                     NormalizationHbt=mean(ChunkRefl(:,(1:(Norm_Time*RawData.dal_fr))),2);
%                                     NormChunk=repmat(NormalizationHbt,1,size(ChunkRefl,2));
%                                     NormEvent=repmat(NormalizationHbt,1,size(EventRefl,2));
%                                     ChunkHbT=(weightedcoeffHbT.*log(ChunkRefl./NormChunk)).*convert2uM;
%                                     EventHbT=(weightedcoeffHbT.*log(EventRefl./NormEvent)).*convert2uM;
%                                     Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkRefl(:,:,AllCount)=ChunkRefl;
%                                     Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkHbT(:,:,AllCount)=ChunkHbT;
%                                     Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedChunkHbT(:,:,AllCount)=ChunkHbT-repmat(ChunkHbT(:,(Lead_Time*RawData.dal_fr)),1,size(ChunkHbT,2));
%                                     Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventRefl{AllCount}=EventRefl;
%                                     Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventHbT{AllCount}=EventHbT;
%                                     Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedEventHbT{AllCount}=EventHbT-repmat(EventHbT(:,(1)),1,size(EventHbT,2));
                                else
                                    NormalizationHbt=mean(ChunkRefl(1:(Norm_Time*RawData.dal_fr)),2);
                                    ChunkHbT=(weightedcoeffHbT.*log(ChunkRefl/NormalizationHbt)).*convert2uM;
                                    EventHbT=(weightedcoeffHbT.*log(EventRefl/NormalizationHbt)).*convert2uM;
                                    Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkRefl(AllCount,:)=ChunkRefl;
                                    Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ChunkHbT(AllCount,:)=ChunkHbT;
                                    Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedChunkHbT(AllCount,:)=ChunkHbT-ChunkHbT(:,(Lead_Time*RawData.dal_fr));
                                    Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventRefl{AllCount}=EventRefl;
                                    Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).EventHbT{AllCount}=EventHbT;
                                    Loco_Data.(TriggerLabel).All.IOS.(ROI_name{ROI_num}).ZeroedEventHbT{AllCount}=EventHbT-EventHbT(:,(1));
                                end
                            end
                            if length(Loco_Data.(TriggerLabel).All.IOS.forepaw.ZeroedEventHbT{AllCount})<15*RawData.dal_fr
                                Loco_Data.(TriggerLabel).All.HemoRespVolTrial(AllCount)=sum(Loco_Data.(TriggerLabel).All.IOS.forepaw.ZeroedChunkHbT(AllCount,((Lead_Time*RawData.dal_fr):end)))/length(Loco_Data.(TriggerLabel).All.IOS.forepaw.ZeroedChunkHbT(AllCount,((Lead_Time*RawData.dal_fr):end)));
                            else
                                Loco_Data.(TriggerLabel).All.HemoRespVolTrial(AllCount)=sum(Loco_Data.(TriggerLabel).All.IOS.forepaw.ZeroedEventHbT{AllCount})/length(Loco_Data.(TriggerLabel).All.IOS.forepaw.ZeroedEventHbT{AllCount});
                            end
                            Loco_Data.(TriggerLabel).All.Neuro.ChunkRawLFP(AllCount,:)=RawData.LFP(NeuroInds);
                            Loco_Data.(TriggerLabel).All.Neuro.ChunkSpectrogram(:,:,AllCount)=S(SpecInds,:);
                            Loco_Data.(TriggerLabel).All.Neuro.NormChunkSpec(:,:,AllCount)=((S(SpecInds,:)-repmat(mean(S(NormSpec,:),1),size(S(SpecInds,:),1),1))./...
                                repmat(mean(S(NormSpec,:),1),size(S(SpecInds,:),1),1))*100;
                            Loco_Data.(TriggerLabel).All.EMG.ChunkRawEMG(AllCount,:)=RawData.MUA(NeuroInds);
                            Loco_Data.(TriggerLabel).All.Ball.ChunkRawBall(AllCount,:)=RawData.vBall(NeuroInds);
                            Loco_Data.(TriggerLabel).All.Neuro.EventSpectrogram{AllCount}=S(EvntInds,:);
                            Loco_Data.(TriggerLabel).All.Neuro.NormEventSpec{AllCount}=((S(EvntInds,:)-repmat(mean(S(NormSpec,:),1),size(S(EvntInds,:),1),1))./...
                                repmat(mean(S(NormSpec,:),1),size(S(EvntInds,:),1),1))*100;
                            if new_T_run(2,k)==9001
                                Loco_Data.(TriggerLabel).All.Neuro.EventRawLFP{AllCount}=RawData.LFP(round(((new_T_run(1,k)/RawData.dal_fr)*RawData.an_fs),0):end);
                                Loco_Data.(TriggerLabel).All.EMG.EventRawEMG{AllCount}=RawData.MUA(round(((new_T_run(1,k)/RawData.dal_fr)*RawData.an_fs),0):end);
                                Loco_Data.(TriggerLabel).All.Ball.EventRawBall{AllCount}=RawData.vBall(round(((new_T_run(1,k)/RawData.dal_fr)*RawData.an_fs),0):end);
                            else
                                Loco_Data.(TriggerLabel).All.Neuro.EventRawLFP{AllCount}=RawData.LFP(round((((new_T_run(1,k):new_T_run(2,k))/RawData.dal_fr)*RawData.an_fs),0));
                                Loco_Data.(TriggerLabel).All.EMG.EventRawEMG{AllCount}=RawData.MUA(round((((new_T_run(1,k):new_T_run(2,k))/RawData.dal_fr)*RawData.an_fs),0));
                                Loco_Data.(TriggerLabel).All.Ball.EventRawBall{AllCount}=RawData.vBall(round((((new_T_run(1,k):new_T_run(2,k))/RawData.dal_fr)*RawData.an_fs),0));
                            end
                            Loco_Data.(TriggerLabel).All.f=f;
                            Loco_Data.(TriggerLabel).All.t=((1:size(Loco_Data.(TriggerLabel).All.Neuro.NormChunkSpec(:,:,AllCount),1))*movingwin(2))-(movingwin(2)+Lead_Time);
                            Loco_Data.(TriggerLabel).All.DurationLabel{AllCount}=RunLabel;
                            Loco_Data.(TriggerLabel).All.fileName{AllCount}=filenm;
                            Loco_Data.(TriggerLabel).All.NegFlag(AllCount)=Loco_Data.(TriggerLabel).All.HemoRespVolTrial(AllCount)>0;
                        else
                            fprintf('Optostim during event. Skipping\n')
                        end
                    end
                end
            end
        end
    end
                             
    %% Identify Stimulus time point(s) and set observation time frame
    params.tapers=[3 5];
    params.trialave= 1;
    params.fpass=[3 15];
    movingwin=[1 0.05];
    params.Fs=RawData.dal_fr;

    if ~isempty(Puff)
        Stim_Type=fieldnames(Puff);
        for n=1:size(Stim_Type,1)
            if createmptystruct==1
                for k=1:size(ROI_name,1)
                    Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).Run.Refl=[];
                    Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).Still.Refl=[];
                end
            end
            Stim=find(diff(Puff.(Stim_Type{n}))==1);
            %% Chunk reflectace data around stimulus
            for U=1:size(Stim,2)
                Follow_Ind=(Stim(U)+(Follow_Time*RawData.dal_fr));
                Loco_Ind=(Stim(U)+(4*RawData.dal_fr));
                Lead_Ind=(Stim(U)-(Lead_Time*RawData.dal_fr));
                Norm_Ind=(Stim(U)-(Norm_Time*RawData.dal_fr));
                Detrend_Ind=(Stim(U)+(Lead_Time*RawData.dal_fr));
                if max(RawData.LED)>1
                    Laser_Ind=(Lead_Time*RawData.dal_fr):((Lead_Time*RawData.dal_fr)+(RawData.AcquistionParams.Laser_Duration*RawData.dal_fr));
                end
                Stim_Ind=Lead_Ind:Follow_Ind;
                Loco_Win=Lead_Ind:Loco_Ind;
                if gt(Lead_Ind,0)
                    if lt(Follow_Ind,size(RawData.IOS.(ROI_name{1}).CBVrefl,2))
                        Isrunning=ismember(Loco_Win,RunInds);
                        if max(Isrunning)==1 %If locomotion takes place during whisker stim window
                            j=2;
                        else
                            j=1;
                        end
                        for k=1:size(ROI_name,1)
                            if strcmpi(ROI_name{k},'Pixelwise')
                                normMatrix=repmat(mean(RawData.IOS.(ROI_name{k}).CBVrefl(:,(Lead_Ind:Norm_Ind)),2),1,length(Stim_Ind));
                                RawRefl=RawData.IOS.(ROI_name{k}).CBVrefl(:,Stim_Ind);
                                NormRefl=(weightedcoeffHbT.*log((RawData.IOS.(ROI_name{k}).CBVrefl(:,Stim_Ind)./normMatrix))).*convert2uM;
                                NormChunk=repmat(NormRefl(:,(Lead_Time*RawData.dal_fr)-1),1,size(NormRefl,2));
                                NormRefl=NormRefl-NormChunk;
                            else
                                Normalizationconstant=mean(RawData.IOS.(ROI_name{k}).CBVrefl(Lead_Ind:Norm_Ind),2);
                                RawRefl=RawData.IOS.(ROI_name{k}).CBVrefl(Stim_Ind);
                                NormRefl=(weightedcoeffHbT.*log(RawData.IOS.(ROI_name{k}).CBVrefl(Stim_Ind)/Normalizationconstant)).*convert2uM;
                                NormRefl=NormRefl-NormRefl((Lead_Time*RawData.dal_fr)-1);
                            end
                            %% Correct for Optostim artifacts
                            if strcmpi(Stim_Type{n},'Laser_Stim')==1
                                if k==1
                                    FlashCatch=abs(diff(RawData.IOS.Optogenetics.CBVrefl(Stim_Ind)));
                                end
                                if max(FlashCatch)>=20
                                    if strcmpi(ROI_name{k},'Pixelwise')
                                        normMatrix=repmat(mean(RawData.IOS.(ROI_name{k}).CBVrefl(:,(Lead_Ind:Norm_Ind)),2),1,length(Stim_Ind));
                                        %HoldRefl=RawRefl-normMatrix;
                                        HoldRefl=RawRefl;
                                        HoldRefl(:,Flash_Vals)=NaN;
                                        [Interp_Data]=fillmissing(HoldRefl,'spline',2);
                                        NormRefl=(weightedcoeffHbT.*log(Interp_Data./normMatrix)).*convert2uM;
                                        NormChunk=repmat(NormRefl(:,(Lead_Time*RawData.dal_fr)-1),1,size(NormRefl,2));
                                        NormRefl=NormRefl-NormChunk;
                                    else
                                        % HoldRefl=RawRefl-mean(RawData.IOS.(ROI_name{k}).CBVrefl(Lead_Ind:Norm_Ind),2);
                                        HoldRefl=RawRefl;
                                        if k==1
                                            %Flash_Points=HoldRefl(Laser_Ind)>=(4*std(RawData.IOS.(ROI_name{k}).CBVrefl(Lead_Ind:Norm_Ind)));
                                            Flash_Points=find(FlashCatch>=20)+1;
                                            First_Flash=Flash_Points(diff(Flash_Points)==1);
                                            %Flash_Vals=Laser_Ind(Flash_Points);
                                            Flash_Vals=Laser_Ind(ismember(Laser_Ind,First_Flash));%Flash_Points));
                                        end
                                        HoldRefl(Flash_Vals)=NaN;
                                        [Interp_Data]=fillmissing(HoldRefl,'spline');
                                        NormRefl=(weightedcoeffHbT.*log(Interp_Data/mean(RawData.IOS.(ROI_name{k}).CBVrefl(Lead_Ind:Norm_Ind),2))).*convert2uM;
                                        NormRefl=NormRefl-NormRefl((Lead_Time*RawData.dal_fr)-1);
                                    end
                                    
                                else
                                    fprintf('no opto stim correction necessary\n')
                                end
                            end
                            %% Separate Puff response based on behavior
                            if strcmpi(ROI_name{k},'Pixelwise')
                                if isempty(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).Refl)
                                    move=1;
                                else
                                    move=size(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).Refl,3)+1;
                                end
                                Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).Refl(:,:,move)=NormRefl;
                                %Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).RawRefl(:,:,move)=RawRefl;
                            else
                                if isempty(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).Refl)
                                    move=1;
                                else
                                    move=size(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).Refl,1)+1;
                                end
                                Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).Refl(move,:)=NormRefl;
                                Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{j}).RawRefl(move,:)=RawRefl;
                             end
                        end
                    end
                end
                Puff_Data.(Stim_Type{n}).StimInds(stimcount,:)=Stim_Ind;
                Puff_Data.(Stim_Type{n}).StimFile{stimcount}=indfile;
                stimcount=stimcount+1;
            end
        end
        createmptystruct=0;
    else
        Stim_Type=[];
    end
    for k=1:size(ROI_name,1)
        Puff_Data.pixelMaps.(ROI_name{k})=RawData.IOS.(ROI_name{k}).PixelMap;
    end
 end
Puff_Data.leadTime=Lead_Time;
Puff_Data.followTime=Follow_Time;
Puff_Data.normTime=Norm_Time;
%save([animal '_' date '_Puff_Data'],'Puff_Data','-v7.3');

%% Separate EMG triggered data based on hemodynamic response volume
% for trigType=1:size(Run_Type,2)
%     if ~isempty(EMG_Data.(Run_Type{trigType}).All.IOS.forepaw.ChunkRefl)
%         for eventNum=1:2
%             if eventNum==1
%                 EventFlag=find(EMG_Data.(Run_Type{trigType}).All.NegFlag==0);%Negative HbT response ==0
%                 EventLabel='negativeHbT';
%             else
%                 EventFlag=find(EMG_Data.(Run_Type{trigType}).All.NegFlag==1);
%                 EventLabel='positiveHbT';
%             end  
%             subfields=fieldnames(EMG_Data.(Run_Type{trigType}).All);
%             for fieldnum=1:size(subfields,1)
%                 if isstruct(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}))
%                     finFields=fieldnames(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}));
%                     for finnum=1:size(finFields,1)
%                         if isstruct(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}))
%                             termfields=fieldnames(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}));
%                             for termNum=1:size(termfields,1)
%                                 if iscell(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}))
%                                     for cellNum=1:length(EventFlag)
%                                         EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}){cellNum}=EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}){EventFlag(cellNum)};
%                                     end
%                                 else
%                                     if size(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),3)>1
%                                         GroupedLabel=['Avg' termfields{termNum}];
%                                         EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})=EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})(:,:,EventFlag);
%                                         EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(GroupedLabel)=mean(EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),3);
%                                     else
%                                         GroupedLabel=['Avg' termfields{termNum}];
%                                         EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})=EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})(EventFlag,:);
%                                         EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(GroupedLabel)=mean(EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),1);
%                                     end
%                                 end
%                             end
%                         else
%                             if iscell(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}))
%                                 for cellNum=1:length(EventFlag)
%                                     EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}){cellNum}=EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}){EventFlag(cellNum)};
%                                 end
%                             else
%                                 if size(EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}),3)>1
%                                     GroupedLabel=['Avg' finFields{finnum}];
%                                     EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum})=EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum})(:,:,EventFlag);
%                                     EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(GroupedLabel)=median(EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}),3);
%                                 else
%                                     GroupedLabel=['Avg' finFields{finnum}];
%                                     EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum})=EMG_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum})(EventFlag,:);
%                                     EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(GroupedLabel)=mean(EMG_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}),1);
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
save([animal '_' date '_processedData_' saveDate],'EMG_Data','-v7.3');
clear EMG_Data
%% Separate Locomotion triggered data by hemodynamic response
% for trigType=1:size(Run_Type,2)
%     if ~isempty(Loco_Data.(Run_Type{trigType}).All.IOS.forepaw.ChunkRefl)
%         for eventNum=1:2
%             if eventNum==1
%                 EventFlag=find(Loco_Data.(Run_Type{trigType}).All.NegFlag==0);%Negative HbT response ==0
%                 EventLabel='negativeHbT';
%             else
%                 EventFlag=find(Loco_Data.(Run_Type{trigType}).All.NegFlag==1);
%                 EventLabel='positiveHbT';
%             end  
%             subfields=fieldnames(Loco_Data.(Run_Type{trigType}).All);
%             for fieldnum=1:size(subfields,1)
%                 if isstruct(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}))
%                     finFields=fieldnames(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}));
%                     for finnum=1:size(finFields,1)
%                         if isstruct(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}))
%                             termfields=fieldnames(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}));
%                             for termNum=1:size(termfields,1)
%                                 if iscell(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}))
%                                     for cellNum=1:length(EventFlag)
%                                         Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}){cellNum}=Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}){EventFlag(cellNum)};
%                                     end
%                                 else
%                                     if size(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),3)>1
%                                         GroupedLabel=['Avg' termfields{termNum}];
%                                         Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})=Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})(:,:,EventFlag);
%                                         Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(GroupedLabel)=mean(Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),3);
%                                     else
%                                         GroupedLabel=['Avg' termfields{termNum}];
%                                         Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})=Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})(EventFlag,:);
%                                         Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(GroupedLabel)=mean(Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),1);
%                                     end
%                                 end
%                             end
%                         else
%                             if iscell(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}))
%                                 for cellNum=1:length(EventFlag)
%                                     Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}){cellNum}=Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}){EventFlag(cellNum)};
%                                 end
%                             else
%                                 if size(Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}),3)>1
%                                     GroupedLabel=['Avg' finFields{finnum}];
%                                     Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum})=Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum})(:,:,EventFlag);
%                                     Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(GroupedLabel)=median(Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}),3);
%                                 else
%                                     GroupedLabel=['Avg' finFields{finnum}];
%                                     Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum})=Loco_Data.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum})(EventFlag,:);
%                                     Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(GroupedLabel)=mean(Loco_Data.(Run_Type{trigType}).(EventLabel).(subfields{fieldnum}).(finFields{finnum}),1);
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end
%%  Perform Data Averaging
Stim_Type=fieldnames(Puff_Data);
for k=1:size(ROI_name,1)
    for n=1:size(Stim_Type,1)
        if isstruct(Puff_Data.(Stim_Type{n}))
            if ~strcmpi(Stim_Type{n},'pixelMaps')
                for m=1:size(Run_State,2)
                    if isempty(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl)==0
                        if strcmpi(ROI_name{k},'Pixelwise')
                            Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Avg_Refl=mean(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl,3);
                            Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Std_Refl=std(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl,0,3);
                        else
                            Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Avg_Refl=mean(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl,1);
                            Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Std_Refl=std(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl,0,1);
                        end
                    else
                        fprintf(1,['No ' Stim_Type{n} ' ' Run_State{m} ' trial data to average\n']);
                    end
                end
            end
        end
    end
end
for u=1:size(Stim_Type,1)
    if isstruct(Puff_Data.(Stim_Type{n}))
        if ~strcmpi(Stim_Type{n},'pixelMaps')
            Stim_Title{u}=strrep(Stim_Type{u},'_',' ');
        end
    end
end
save([animal '_' date '_processedData_' saveDate],'Puff_Data','Loco_Data','-append');
%% Power spectral analysis for heartrate/breathing
%     for k=1:(size(ROI_name{k})-1)
%         for n=1:size(Stim_Type,1)
%             if isempty(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{1}).Refl)==0
%                 [S,f]= mtspectrumc(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).Normalizedrefl',params);
%                 figure(95);hold on;
%                 semilogy(f,S);
%                 axis tight
%                 title(['Power Spectrum of Whole Whisker Stimulation Session ' ROI_name{k} ' ' Stim_Type{n}],'FontSize',16,'FontWeight','bold','FontName','Arial');
%                 xlabel('Frequency (Hz)','FontSize',10,'FontWeight','bold','FontName','Arial');
%                 ylabel('Power','FontSize',10,'FontWeight','bold','FontName','Arial');
%                 legend(Stim_Title, 'Location','southeast','Orientation','vertical','FontSize',8,'FontWeight','bold','FontName','Arial');
%                 savefig([animal '_' date '_' ROI_name{k} '_' Stim_Type{n} '_Power spectrum of Hemodynamics']);
%             end
%         end
%     end

%% Plot and Save Files
if ~isempty(Stim_Type)
    for m=1:size(Run_State,2)
        lgnd_cnt=1;
        for k=1:(size(ROI_name,1)-1)
            figure(33);
            hold on;
            for n=1:size(Stim_Type,1)
                if isstruct(Puff_Data.(Stim_Type{n}))
                    if ~strcmpi(Stim_Type{n},'pixelMaps')
                        if isempty(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl)==0 %isempty(Puff_Data.(Stim_Type{n}).(Run_State{m}).Image)==0
                            Time=((1:size(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Avg_Refl,2))/RawData.dal_fr)-Lead_Time;
                            
                            plot(Time,filtfilt(sos,g,Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Avg_Refl)); %Plot of avg CBV change in trials when animal was still 3-3-15 KG
                            LgndTxt{lgnd_cnt}=[strrep(Stim_Type{n},'_',' ') ' '  strrep(ROI_name{k},'_',' ')];
                            lgnd_cnt=lgnd_cnt+1;
                        end
                    end
                end
            end
        end
        if exist('LgndTxt','var')
        %if isempty(Puff_Data.(Stim_Type{n}).IOS.(ROI_name{k}).(Run_State{m}).Refl)==0
            Zero_Line(1:(((Lead_Time+Follow_Time)*Puff_Data.dal_fr)+1))=0;
            Time=((1:size(Zero_Line,2))/Puff_Data.dal_fr)-Lead_Time;
            plot(Time,Zero_Line,'k');
            title(['Average stimulus evoked reflectance change ' Run_State{m} ' trials'],'FontSize',14,'FontWeight','bold','FontName','Arial');
            xlabel('Time (s)','FontSize',10,'FontWeight','bold','FontName','Arial');
            ylabel('Delta \muM HbT','FontSize',10,'FontWeight','bold','FontName','Arial');
            legend(LgndTxt, 'Location','southeast','Orientation','vertical','FontSize',8,'FontWeight','bold','FontName','Arial');
            savefig([animal '_' date '_' Run_State{m} '_Normalized Average Stimulus Evoked Reflectance Change Per Session']);
            close;
        end
        %end
    end
%% Visualize pixelwise data
% frameTimes=[(Lead_Time-1), Lead_Time,(Lead_Time+0.5),(Lead_Time+1),(Lead_Time+1.5),(Lead_Time+2),(Lead_Time+2.5),(Lead_Time+3),(Lead_Time+3.5),(Lead_Time+4),(Lead_Time+4.5),(Lead_Time+5)]*Puff_Data.dal_fr;
if ~isempty(Puff_Data.Laser_Stim.IOS.Pixelwise.(Run_State{m}).Refl)
frameTimes=round(((Lead_Time-1):0.5:(Lead_Time+5.5))*Puff_Data.dal_fr,0);
labelTimes=(frameTimes/Puff_Data.dal_fr)-Lead_Time;
EmptyFrame((1:256),(1:256))=NaN;
PixMap=Puff_Data.pixelMaps.Pixelwise;
colnum=3;
rownum=ceil(numel(frameTimes)/colnum);
for n=1:size(Stim_Type,1)
    if isstruct(Puff_Data.(Stim_Type{n}))
        if ~strcmpi(Stim_Type{n},'pixelMaps')
            roiFrame=EmptyFrame;
            if strcmpi(Stim_Type{n},'Laser_Stim')
                %         for pic=1:length(SharVars.ROIs.Optogenetics.yi)
                %             roiFrame(SharVars.ROIs.Optogenetics.yi(pic),SharVars.ROIs.Optogenetics.xi(pic))=1;
                %         end
                roiFrame(round(SharVars.ROIs.Optogenetics.ycenter,0),round(SharVars.ROIs.Optogenetics.xcenter,0))=1;
                rotROI=rot90(roiFrame,3);
                [row,col]=find(rotROI==1);
                TheROI(:,1)=row;
                TheROI(:,2)=col;
            else
                for pic=1:length(SharVars.ROIs.barrels.xi)
                    roiFrame(round(SharVars.ROIs.barrels.yi(pic),0),round(SharVars.ROIs.barrels.xi(pic),0))=1;
                end
                rotROI=rot90(roiFrame,3);
                [row,col]=find(rotROI==1);
                [TheROI(:,1),Inds]=sort(col);
                TheROI(:,2)=row(Inds);
                midRow=median(TheROI(:,2));
                count=1;
                bump=1;
                for vert=1:length(TheROI)
                    if TheROI(vert,2)<midRow
                        UpperInds(count,1)=TheROI(vert,1);
                        UpperInds(count,2)=TheROI(vert,2);
                        count=count+1;
                    else
                        LowerInds(bump,1)=TheROI(vert,1);
                        LowerInds(bump,2)=TheROI(vert,2);
                        bump=bump+1;
                    end
                end
                [newLower(:,1),Inds]=sort(LowerInds(:,1),'descend');
                newLower(:,2)=LowerInds(Inds,2);
                visROI=[UpperInds;newLower];
                visROI(length(visROI)+1,:)=visROI(1,:);
                clear UpperInds LowerInds
            end
            theStim=strrep(Stim_Type{n},'_',' ');
            for m=1:size(Run_State,2)
                if ~isempty(Puff_Data.(Stim_Type{n}).IOS.Pixelwise.(Run_State{m}).Refl)
                    for framenum=1:length(frameTimes)
                        theFrame(:,framenum)=mean(Puff_Data.(Stim_Type{n}).IOS.Pixelwise.(Run_State{m}).Avg_Refl(:,(frameTimes(framenum)-2):(frameTimes(framenum)+2)),2);
                        MakePic(:,:,framenum)=EmptyFrame;
                        for pixnum=1:size(theFrame,1)
                            MakePic(PixMap(1,pixnum),PixMap(2,pixnum),framenum)=theFrame(pixnum,framenum);
                        end
                        plotLabel=['t= ' num2str(labelTimes(framenum))];
                        if framenum==1
                            theImg=EmptyFrame;
                            for frame=1:5
                                for pixnum=1:size(RawData.IOS.Pixelwise.PixelMap,2)
                                    theImg(PixMap(1,pixnum),PixMap(2,pixnum),frame)=RawData.IOS.Pixelwise.CBVrefl(pixnum,frame);
                                end
                                RefImg=mean(theImg,3);
                                [row,col]=find(~isnan(rot90(RefImg,3)));
                                row=sort(row); %find bounds for image display
                                col=sort(col); %find bounds for image display
                                
                            end
                            figure(99);imgAX=subplot(rownum,colnum,framenum);imagesc(rot90(RefImg,3));
                            colormap(imgAX,gray);caxis([0 4095]); axis image; axis off;
                            xlim([col(1) col(end)]);ylim([row(1),row(end)]);
                        end
                        figure(99);theAX=subplot(rownum,colnum,framenum+1);imagesc(rot90(real(MakePic(:,:,framenum)),3));
                        colormap(theAX,parula); axis image; axis off;
                        xlim([col(1) col(end)]);ylim([row(1),row(end)]);
                        title(plotLabel);
                        if strcmpi(Stim_Type{n},'Laser_Stim')
                            caxis([-5 5]);
                            drawcircle('Center',[TheROI(2),TheROI(1)],'Radius',15,'Color',[1 1 1],'LineWidth',0.25,'FaceAlpha',0,'InteractionsAllowed','none');
                        else
                            caxis([-5 5]);
                            drawpolygon('Position',visROI,'Color',[1 1 1],'FaceAlpha',0,'LineWidth',0.25,'InteractionsAllowed','none');
                        end
                        
                    end
                    runTitle={'still','running'};
                    theTitle=['Average response to ' theStim ' while animal was ' runTitle{m}];
                    sgtitle(theTitle);
                    savefig([animal '_' date '_' Stim_Type{n} '_' Run_State{m} '_Pixelwise Average Stimulus Evoked Reflectance Change']);
                    close;
                end
            end
            clear TheROI
        end
    end
end
%save([animal '_' date '_Puff_Data'],'Puff_Data','-v7.3');
close all;
end
end
end