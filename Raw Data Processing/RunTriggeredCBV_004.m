 function RunTriggeredCBV_004(filename,theLED,bandFilter,cutFilter,varargin)
%Function to determine CBV change due to running in two ROI defined in
%ProcessWhiskerTrials.m from data gathered on Intrinsic imaging rig using
%Dalsa 1M60 Camera
%Written by Kyle Gheres, April 2016

%Files Input: RawData.mat generated by ProcessWhiskerTrials.m
%Functions Called: velocity_binarize.m, motion_cont_3.m,
%mtspecgramc.m,Puff_Win.m, Plot_ROI_CBV_Trial.m

%% Identify and load files
if nargin==0
filename = uigetfile('*.mat','MultiSelect', 'on');
end
% Control for single file instead of a list
if iscell(filename) == 0
    filename = {filename};
end

%% Get constants for converting reflectance to Hb concentration
if ~exist('theLED','var')
    theLED=input('Enter model of LED used for illumination e.g. M530L3. ','s');
end
if ~exist('bandFilter','var')
    bandFilter=input('Enter model of optical bandpass filter in light path. e.g FB570-10. ','s');
end
if ~exist('cutFilter','var')
    cutFilter=input('Enter model of optical cut off filter in light path. e.g FEL0500. ','s');
end
[weightedcoeffHbO,weightedcoeffHbR,weightedcoeffHbT]=getHbcoeffs(theLED,bandFilter,cutFilter);
convert2uM=1e6; % Constant to convert \DeltaHbT measurements from M to \muM

%% Constants
T_seg=0.20; %Length of running duration to be considered a single locomotion event
T_beg=5.0; %Length of still time between events to consider animal at rest
Acc_Thresh=1e-5;
Running_Event_Lead=5;%Time in seconds to preceed running onset
Lead_Time=5;%Time to preceed whisker puff
Follow_Time=15;% Time to follow Whisker Puff
Norm_Time=2;%Time after Lead Time to use as normalization window
Run_Type={'Volitional','Evoked'};
Run_Length={'Twitch','Short','Medium','Long','Very_Long','Extra_Long'};
Run_Cut_Off=[0.5 2 5 10 15];
params.tapers=[3 5];
params.trialave= 1;
params.fpass=[0 100];
movingwin=[1 0.05];
stimcount=1;

createmptystruct=1;
Proc_Data.swipefiles=[];
Proc_Data.Volitional.ROI.Run_Events={};
Proc_Data.Volitional.Run_length=[];
Proc_Data.Evoked.ROI.Run_Events={};
Proc_Data.Evoked.Run_length=[];
new_T_run=[];
Whisker_Puff=zeros(1,10);

%% Process RawData.mat files for Running and or Puffing
for fil = 1:length(filename)
    close all;
    indfile = filename{fil};
    dashLocs=strfind(indfile,'_');
    animal=indfile(1:(dashLocs(2)-1));%(1:9);
    date=indfile((dashLocs(3)+1):(dashLocs(4)-1));%(10:15);(14:19);
    hem=indfile((dashLocs(2)+1):(dashLocs(3)-1));%(7:8);%(11:12);
    filenm= indfile((dashLocs(4)+1):(dashLocs(end)-1));%(10:26);%(14:30);
    load(indfile);
    
    ROI_name=fieldnames(RawData.IOS);
    reOrderROI=strcmpi(ROI_name,'Pixelwise')==1;
    ROI_name{reOrderROI}=ROI_name{length(ROI_name)};
    ROI_name{length(ROI_name)}='Pixelwise';
    
    params.Fs=RawData.an_fs;

    Frame_spacing=2;
    order=3;
    Wn= 1/(0.5*RawData.dal_fr); %Bandpass filter properties [low freq, high freq]/nyquist freq
    ftype='low';
    
    [zeroa,poleb,gain]=butter(order,Wn,ftype);
    [CBV_sos,CBV_g]=zp2sos(zeroa,poleb,gain);
    
    %% Solenoid firing times
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
        
    else
        fprintf('No stimulus during trial\n')
        Puff=[];
    end
    
    %%  Run Triggered CBV Analysis Barrels and Fp/Hp ROI

    fprintf(1,'binarizing locomotion\n')
    [imp_bin,velocity]=velocity_binarize(RawData.vBall,RawData.an_fs,RawData.dal_fr,Acc_Thresh);
    [T_run,~,new_T_run]=motion_cont_2(imp_bin,RawData.dal_fr,T_seg,T_beg);
    Isrunning=isempty(new_T_run);
    if Isrunning==0 %If animal is running==0
        RunEvents(1:length(RawData.vBall))=0;
        fprintf(1,'chunking reflectance around running\n')
        for k=1:size(new_T_run,2)
            try
                RunEvents(new_T_run(1,k):new_T_run(2,k))=1;
                if lt((new_T_run(2,k)-new_T_run(1,k)),(15*RawData.dal_fr))
                    RunInds=((new_T_run(1,k)-2*RawData.dal_fr):(new_T_run(2,k)+3*RawData.dal_fr));
                else
                    if lt((new_T_run(2,k)-new_T_run(1,k)),(30*RawData.dal_fr))
                        try
                            RunInds=((new_T_run(1,k)-Running_Event_Lead*RawData.dal_fr):(new_T_run(2,k)+((30*RawData.dal_fr)-(new_T_run(2,k)-new_T_run(1,k)))));
                        catch
                            RunInds=((new_T_run(1,k)-Running_Event_Lead*RawData.dal_fr):(new_T_run(2,k)+(size(RawData.CBVrefl_fp,2)-new_T_run(2,k))));
                        end
                    else
                        RunInds=((new_T_run(1,k)-Running_Event_Lead*RawData.dal_fr):(new_T_run(2,k)+3*RawData.dal_fr));
                    end
                end
                Puff_Triggered_Run=ismember((RunInds(1)+(2*RawData.dal_fr)),Puff_windows);
                if Puff_Triggered_Run==0
                    r=1;
                else
                    r=2;
                end
                try
                    if max(abs(diff(RawData.barrels.CBVrefl_barrels(RunInds))),[],2)>50 %Filters out trials where hand obstructs window 3-23-15 KG
                        Swipe_Trial=RawData.barrels.CBVrefl_barrels(RunInds);
                    else
                        move=size(Proc_Data.(Run_Type{r}).ROI.Run_Events,2)+1;
                        Proc_Data.(Run_Type{r}).Run_length(move)=((new_T_run(2,k)-new_T_run(1,k))/RawData.dal_fr);
                        Proc_Data.(Run_Type{r}).ROI.Run_Events{1,move}=RawData.barrels.CBVrefl_barrels(RunInds);
                        Proc_Data.(Run_Type{r}).ROI.Run_Events{2,move}=RawData.fp.CBVrefl_fp(RunInds);
                        Proc_Data.(Run_Type{r}).ROI.Norm_Run_Events{1,move}=(Proc_Data.(Run_Type{r}).ROI.Run_Events{1,move}-mean(Proc_Data.(Run_Type{r}).ROI.Run_Events{1,move}(1:45)))/mean(Proc_Data.(Run_Type{r}).ROI.Run_Events{1,move}(1:45))*100; % Normalize reflectance values to avg reflectance of 1.5s before running
                        Proc_Data.(Run_Type{r}).ROI.Norm_Run_Events{2,move}=(Proc_Data.(Run_Type{r}).ROI.Run_Events{2,move}-mean(Proc_Data.(Run_Type{r}).ROI.Run_Events{2,move}(1:45)))/mean(Proc_Data.(Run_Type{r}).ROI.Run_Events{2,move}(1:45))*100; % Normalize reflectance values to avg reflectance of 1.5s before running
                    end
                catch
                    fprintf(1,'Error in Chunking reflectance to running events\n')
                end
            catch
                fprintf(1,'Error in setting running event Index\n')
            end
        end
        clear theimage;
    else
        fprintf(1,'No running in this trial\n')
    end
    %% Plot whole trial binarized locomotion with normalized reflectance from FP/HP and barrels ROI
    Plot_ROI_CBV_Trial_Optogenetics_005(filenm,RawData,imp_bin,T_run,new_T_run,velocity,weightedcoeffHbO,weightedcoeffHbR,weightedcoeffHbT)
end
%% Run Event Aggregation/Averaging
for q=1:size(Run_Type,2)
    for d=1:(size(Run_Cut_Off,2)+1)
        if gt(d,size(Run_Cut_Off,2))
            Chunk_Running.(Run_Type{q}).(Run_Length{d}).TF=gt(Proc_Data.(Run_Type{q}).Run_length,Run_Cut_Off(d-1 ));
        else
            Chunk_Running.(Run_Type{q}).(Run_Length{d}).TF=le(Proc_Data.(Run_Type{q}).Run_length,Run_Cut_Off(d));
        end
        Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index=find(Chunk_Running.(Run_Type{q}).(Run_Length{d}).TF==1);
        if gt(d,1)
            p=1;
            while gt((d-p),0)
                Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index=setdiff(Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index,Chunk_Running.(Run_Type{q}).(Run_Length{d-p}).Index);
                p=p+1;
            end
        end
        Chunk_Running.(Run_Type{q}).(Run_Length{d}).Is_Running=isempty(Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index);
    end
    
    %% Visualize Run Evoked CBV Response
     for d=1:(size(Run_Cut_Off,2)+1)
        clear j;
        if Chunk_Running.(Run_Type{q}).(Run_Length{d}).Is_Running==0
            fprintf(1,'Averaging responses\n')
            for j=1:length(Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index)
                try
                    Chunk_Running.(Run_Type{q}).(Run_Length{d}).lengths(j)=size(Proc_Data.(Run_Type{q}).ROI.Run_Events{1, Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index(j)},2);
                catch
                    fprintf(1,['Error in obtaining ' Run_Length{d} ' running lengths\n'])
                end
            end
            
            Index=min(Chunk_Running.(Run_Type{q}).(Run_Length{d}).lengths);
            Frame_num=floor(Index/(Frame_spacing*RawData.dal_fr));
            for j=1:length(Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index)
                try
                    Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).barrels(j,(1:Index))=Proc_Data.(Run_Type{q}).ROI.Norm_Run_Events{1,Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index(j)}(1:Index);
                    Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).fp(j,(1:Index))=Proc_Data.(Run_Type{q}).ROI.Norm_Run_Events{2,Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index(j)}(1:Index);
                    for k=1:(Frame_num+1)
                        if k==1
                            Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Frame{k,j}=Proc_Data.(Run_Type{q}).Image.RunVideo{Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index(j)}(:,:,(1));
                        elseif le(k,Frame_num)
                            Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Frame{k,j}=Proc_Data.(Run_Type{q}).Image.RunVideo{Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index(j)}(:,:,((k*RawData.dal_fr)/5));
                        else
                            Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Frame{k,j}=Proc_Data.(Run_Type{q}).Image.RunVideo{Chunk_Running.(Run_Type{q}).(Run_Length{d}).Index(j)}(:,:,(end));
                        end
                    end
                catch
                    fprintf('Error in Frame assembly\n')
                end
            end
            
            for w=1:size(Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Frame,1)
                Single_Time_Image=[];
                for p=1:size(Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Frame,2)
                    Single_Time_Image=cat(3,Single_Time_Image,Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Frame{w,p});
                end
                Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Avg_Images(:,:,w)=mean(Single_Time_Image,3);           
            end
            Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Avg_Response_barrels=mean(Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).barrels,1);
            Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Avg_Response_fp=mean(Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).fp,1);
            Filtered_Avg_Barrels=filtfilt(sos,g,Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Avg_Response_barrels);
            Filtered_Avg_Forepaw=filtfilt(sos,g,Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Avg_Response_fp);
            Time=((1:length(Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).Avg_Response_barrels))/RawData.dal_fr)-2;
            
            figure(d); hold on; %plot run evoked changes in barrel ROI
            plot(Time,Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).barrels');
            plot(Time,Filtered_Avg_Barrels,'LineWidth',4);
            title([Run_Type{q} ' ' Run_Length{d} ' Triggered Normalized Reflectance Change Barrels ROI'],'FontSize',14,'FontWeight','bold');
            xlabel('Time (s)','FontSize',10,'FontWeight','bold');
            ylabel('Percent Change (%)','FontSize',10,'FontWeight','bold');
            savefig([animal '_' date '_' Run_Type{q} '_' Run_Length{d} '_Triggered Normalized Reflectance Change Barrels']);
            close;
            
            figure(d*10); hold on;
            plot(Time,Proc_Data.(Run_Type{q}).Run_Events.(Run_Length{d}).fp');
            plot(Time,Filtered_Avg_Forepaw,'LineWidth',4);
            title([Run_Type{q} ' ' Run_Length{d} ' Triggered Normalized Reflectance Change Forepaw ROI'],'FontSize',14,'FontWeight','bold');
            xlabel('Time (s)','FontSize',10,'FontWeight','bold');
            ylabel('Percent Change (%)','FontSize',10,'FontWeight','bold');
            savefig([animal '_' date '_' Run_Type{q} '_' Run_Length{d} '_Triggered Normalized Reflectance Change Forepaw']);
            close;
            
            figure((d*10)+1);hold on;
            plot(Time,Filtered_Avg_Barrels,'b');
            plot(Time,Filtered_Avg_Forepaw,'r');
            title([Run_Type{q} ' ' Run_Length{d} 'Average Triggered Normalized Reflectance Change'],'FontSize',14,'FontWeight','bold');
            xlabel('Time (s)','FontSize',10,'FontWeight','bold');
            ylabel('Percent Change (%)','FontSize',10,'FontWeight','bold');
            legend({'Normalized Reflectance Barrels','Normalized Reflectance Forepaw/Hindpaw'}, 'Location','southeast','Orientation','vertical','FontSize',8,'FontWeight','bold');
            savefig([animal '_' date '_' Run_Type{q} '_' Run_Length{d} '_Avg Twitch Triggered Normalized Reflectance Change']);
            close;
        else
            fprintf(1,'No Events\n')
        end
    end
    
    %% Visualize Running Evoked Pixelwise Time Course
    if isfield(Proc_Data.(Run_Type{q}),'Run_Events')==1
    Run_Event=fieldnames(Proc_Data.(Run_Type{q}).Run_Events);
    for d=1:(size(Run_Event,1))
        
        for u=1:size(Run_Event,1)
            Run_Title{u}=strrep(Run_Event{u},'_',' ');
        end
            Num_Rows=ceil((size(Proc_Data.(Run_Type{q}).Run_Events.(Run_Event{d}).Avg_Images,3)+1)/5);
            if Num_Rows==1
                Num_Columns=size(Proc_Data.(Run_Type{q}).Run_Events.(Run_Event{d}).Avg_Images,3)+1;
            else
                Num_Columns=5;
            end
            The_Images=(double(Proc_Data.(Run_Type{q}).Run_Events.(Run_Event{d}).Avg_Images))/10000;
            Seconds_After_Running_Onset=((1:size(The_Images,3))-Running_Event_Lead)*2;
            for c=1:size(Seconds_After_Running_Onset,2)
                Run_Time{c}=num2str(Seconds_After_Running_Onset(c));
            end
            figure(d);hold on;
            for v=1:(size(The_Images,3)+1)
                if v==1
                subplot(Num_Rows,Num_Columns,1);imagesc(Avg_Img);
                axis image
                ax=gca;
                ax.XTick=[];
                ax.XTickLabel=[];
                ax.XLabel.String='Lateral';
                ax.XLabel.FontSize=10;
                ax.XLabel.FontWeight='bold';
                ax.YTick=[];
                ax.YTickLabel=[];
                ax.YLabel.String='Rostral';
                ax.YLabel.FontSize=10;
                ax.YLabel.FontWeight='bold';
                caxis([-0 4500]);
                colormap(hot);
                title('Average Reflectance','FontSize',10,'FontWeight','bold');
                else
                subplot(Num_Rows,Num_Columns,v);imagesc(The_Images(:,:,(v-1)));
                axis image
                axis off
                caxis([-0.05 0.05]);
                colormap(hot);
                title([Run_Time{v-1} 's after onset'],'FontSize',10,'FontWeight','bold');
                end
            end
            ax=axes('Units','Normal','Position',[.075 .1 .85 .85],'Visible','off');
            set(get(ax,'Title'),'Visible','on')
            title([Run_Type{q} ' ' Run_Title{d} ' ' animal ' ' date],'FontSize',14,'FontWeight','bold');
            savefig([animal '_' date '_' Run_Type{q} '_' Run_Event{d} '_Time Series']);    
    end
    end
    close all;
end
close all;
Proc_Data.Volitional.Image.RunVideo=[];
Proc_Data.Evoked.Image.RunVideo=[];
clearvars -except Proc_Data animal date;
fprintf(1,'Saving processed data\n')
save([animal '_' date '_procdata'],'Proc_Data','-v7.3');
end

       