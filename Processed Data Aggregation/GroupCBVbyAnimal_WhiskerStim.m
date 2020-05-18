function GroupCBVbyAnimal_WhiskerStim(~)
%Function to aggregate hemodynamic data by animal and day.
%Written by Kyle Gheres, Nov 2014
%% Identify Files
directoryContents=dir('D:/');
folderNum=4:10;%(size(directoryContents,1)-2);
for q=1:length(folderNum)
    if directoryContents(folderNum(q)).isdir==1
        cd([directoryContents(folderNum(q)).folder  directoryContents(folderNum(q)).name]);
        searchDir=cd;
        fileList=dir(fullfile(searchDir,'**','*processeddata.mat'));
        %% Get Data from each animal
        for filnum=1:size(fileList,1)
            close all;
            load([fileList(filnum).folder '\' fileList(filnum).name]);
            %% Create Filter
            FrameRate=ProcessedData.dal_fr;%Dalsa image acquisition rateorder=3;
            Wn=1/(0.5*FrameRate);
            order=3;
            ftype='low';
            [z,p,k]=butter(order,Wn,ftype);
            [sos,g]=zp2sos(z,p,k);
            Animal_Age=searchDir(4:end);%input(prompt,'s');
            %TempTypes=fieldnames(ProcessedData);
            StimTypes={'Laser_Stim','Contra_Puff','Ipsi_Puff','Control_Puff'};%TempTypes(3:6);
            %% Correct channel order
            for count=2:4
                stimCounts(count-1)=size(ProcessedData.(StimTypes{count}).IOS.barrels.Run.Refl,1)+size(ProcessedData.(StimTypes{count}).IOS.barrels.Still.Refl,1);
            end
            MaxStim=find(stimCounts==max(stimCounts));
            if MaxStim~=1
                HoldStim=ProcessedData.Contra_Puff;
                ProcessedData.Contra_Puff=ProcessedData.(StimTypes{MaxStim+1});
                ProcessedData.(StimTypes{MaxStim+1})=HoldStim;
            end
            animalname=fileList(filnum).name(1:5);
            if strcmpi(animalname,'NC170')
                HoldStim=ProcessedData.Ipsi_Puff;
                ProcessedData.Ipsi_Puff=ProcessedData.Control_Puff;
                ProcessedData.Control_Puff=HoldStim;
            end
            %% Get Stimulus data
            for stimNum=1:length(StimTypes)
                if isfield(ProcessedData.(StimTypes{stimNum}),'IOS')
                    Tempnames=fieldnames(ProcessedData.(StimTypes{stimNum}).IOS);
                    ROInames=Tempnames(1:2);
                    for ROInum=1:length(ROInames)
                        behaviorTypes=fieldnames(ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}));
                        for beNum=1:length(behaviorTypes)
                            if isfield(ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}),'Avg_Refl')
                                if ProcessedData.dal_fr~=30
                                    TempRefl=ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}).Refl;
                                    FlashCatch=diff(TempRefl,1,2);
                                    [FlashRow,FlashCol]=find(abs(FlashCatch)>10);
                                    TempRefl(FlashRow,FlashCol)=NaN;
                                    CorrectedRefl=fillmissing(TempRefl,'spline',2);
                                    SmoothRefl=filtfilt(sos,g,CorrectedRefl')';
                                    OffsetConst=SmoothRefl(:,150);
                                    OffsetMat=repmat(OffsetConst,1,size(SmoothRefl,2));
                                    OffsetRefl=SmoothRefl-OffsetMat;%mean(SmoothRefl(1:(ProcessedData.leadTime*ProcessedData.dal_fr)));
                                    AnimalAvg=mean(OffsetRefl,1);
                                    AnimalStd=std(OffsetRefl,[],1);
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Avg_HbT{filnum}=resample(AnimalAvg,30,ProcessedData.dal_fr);%resample(ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}).Avg_Refl,30,ProcessedData.dal_fr);
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Std_HbT{filnum}=resample(AnimalStd,30,ProcessedData.dal_fr);%resample(ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}).Std_Refl,30,ProcessedData.dal_fr);
                                else
                                    TempRefl=ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}).Refl;
                                    FlashCatch=diff(TempRefl,1,2);
                                    [FlashRow,FlashCol]=find(abs(FlashCatch)>15);
                                    TempRefl(FlashRow,FlashCol)=NaN;
                                    CorrectedRefl=fillmissing(TempRefl,'spline',2);
                                    SmoothRefl=filtfilt(sos,g,CorrectedRefl')';
                                    OffsetConst=SmoothRefl(:,150);
                                    OffsetMat=repmat(OffsetConst,1,size(SmoothRefl,2));
                                    OffsetRefl=SmoothRefl-OffsetMat;%mean(SmoothRefl(1:(ProcessedData.leadTime*ProcessedData.dal_fr)));
                                    AnimalAvg=mean(OffsetRefl,1);
                                    AnimalStd=std(OffsetRefl,[],1);
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Avg_HbT{filnum}=AnimalAvg;%ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}).Avg_Refl;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Std_HbT{filnum}=AnimalStd;%ProcessedData.(StimTypes{stimNum}).IOS.(ROInames{ROInum}).(behaviorTypes{beNum}).Std_Refl;
                                end
                            end
                        end
                    end
                else
                    fprintf('Something has gone horribly wrong\n')
                    keyboard
                end
            end
        end
        %% Constants
        FrameRate=ProcessedData.dal_fr;%Dalsa image acquisition rateorder=3;
        Wn=1/(0.5*FrameRate);
        order=3;
        ftype='low';
        [z,p,k]=butter(order,Wn,ftype);
        [sos,g]=zp2sos(z,p,k);
        leadtime=ProcessedData.leadTime;%time in seconds before puff
        followtime=ProcessedData.followTime; %time in seconds after puff
        peakWin=5;
        earlyResp=0.5;
        midResp=3;
        lateResp=ProcessedData.followTime;
        expectedLength=(ProcessedData.leadTime+ProcessedData.followTime)*ProcessedData.dal_fr;
        StimTypes=fieldnames(GroupedData.AnimalData);
        
        for stimNum=1:length(StimTypes)
            for ROInum=1:length(ROInames)
                for beNum=1:length(behaviorTypes)
                    if isfield(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS,'Avg_HbT')
                        for ind=1:length(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Avg_HbT)
                            if isempty(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Avg_HbT{ind})
                                tempRefl(ind,(1:expectedLength))=NaN;
                                lowRefl(ind,(1:expectedLength))=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinVal=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxVal=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxLoc=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinLoc=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.early(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.mid(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.late(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.core(ind)=NaN;
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.total(ind)=NaN;
                            else
                                tempRefl(ind,:)=GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Avg_HbT{ind}(1:expectedLength);
                                lowRefl(ind,:)=tempRefl(ind,:);
                                %% Find Peak and duration of hemodynamic response
                                %[peaks,locs,width]=findpeaks(lowRefl(ind,(leadtime*FrameRate):((leadtime*FrameRate)+(peakWin*FrameRate))),'MinPeakDistance',(2*FrameRate));
                                
                                PosPeakFind=find(diff((diff(lowRefl(ind,:))>0))==-1)+1; %Sample # from t=-5s at which change in HbT goes from + to -
                                PeaksInWin=PosPeakFind>=(leadtime*FrameRate) & PosPeakFind<=((leadtime*FrameRate)+(peakWin*FrameRate));
                                keepPeaks=PosPeakFind(PeaksInWin);
                                PositivePeaks=lowRefl(ind,keepPeaks)>0;
                                finalPeaks=keepPeaks(PositivePeaks);
                                maxPosPeaks=max(lowRefl(ind,finalPeaks));
                                locPosPeaks=finalPeaks(lowRefl(ind,finalPeaks)==maxPosPeaks);
                                if ~isempty(locPosPeaks)
                                    leftEdge=locPosPeaks-find(flip(lowRefl(ind,((leadtime*FrameRate):locPosPeaks)))<=(0.5*maxPosPeaks),1,'first')+1;
                                    rightEdge=find(lowRefl(ind,(locPosPeaks:end))<=(0.5*maxPosPeaks),1,'first')+locPosPeaks;%((leadtime*FrameRate)+(peakWin*FrameRate))
                                    PosPeakWidth=rightEdge-leftEdge;
                                else
                                    fprintf('No positive peak found\n')
                                end
                                
                                NegPeakFind=find(diff((diff(lowRefl(ind,:))>0))==1)+1; %Sample # from t=-5s at which change in HbT goes from - to +
                                PeaksInWin=NegPeakFind>=(leadtime*FrameRate) & NegPeakFind<=((leadtime*FrameRate)+(peakWin*FrameRate));
                                keepPeaks=NegPeakFind(PeaksInWin);
                                NegativePeaks=lowRefl(ind,keepPeaks)<0;
                                finalPeaks=keepPeaks(NegativePeaks);
                                maxNegPeaks=min(lowRefl(ind,finalPeaks));
                                locNegPeaks=finalPeaks(lowRefl(ind,finalPeaks)==maxNegPeaks);
                                if ~isempty(locNegPeaks)
                                    leftEdge=locNegPeaks-find(flip(lowRefl(ind,(leadtime*FrameRate):locNegPeaks))<=(0.5* maxNegPeaks),1,'first')+1;
                                    rightEdge=find(lowRefl(ind,(locNegPeaks:end))<=(0.5* maxNegPeaks),1,'first')+locNegPeaks;%((leadtime*FrameRate)+(peakWin*FrameRate))
                                    NegPeakWidth=rightEdge-leftEdge;
                                else
                                    fprintf('No negative peak found\n')
                                end
                                
                                if ~isempty(maxPosPeaks)
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak(ind)=maxPosPeaks;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT(ind)=locPosPeaks;
                                    if ~isempty(PosPeakWidth)
                                        GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM(ind)=PosPeakWidth;
                                    else
                                        GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM(ind)=NaN;
                                    end
                                else
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak(ind)=NaN;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT(ind)=NaN;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM(ind)=NaN;
                                end
                                
                                if ~isempty(maxNegPeaks)
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak(ind)=maxNegPeaks;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT(ind)=locNegPeaks;
                                    if ~isempty(NegPeakWidth)
                                        GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM(ind)=NegPeakWidth;
                                    else
                                        GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM(ind)=NaN;
                                    end
                                else
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak(ind)=NaN;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT(ind)=NaN;
                                    GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM(ind)=NaN;
                                end
                                %% Find min and max values in post stim windows ~5s
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinVal(ind)=min(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(peakWin*FrameRate)))));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinLoc(ind)=find(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(peakWin*FrameRate))))==GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinVal(ind));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxVal(ind)=max(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(peakWin*FrameRate)))));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxLoc(ind)=find(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(peakWin*FrameRate))))==GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxVal(ind));
                                %% Calculate volume of hemodynamic response
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.early(ind)=sum(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(earlyResp*FrameRate)))))...
                                    /length(((leadtime*FrameRate):((leadtime*FrameRate)+(earlyResp*FrameRate))));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.mid(ind)=sum(lowRefl(ind,(((leadtime*FrameRate)+(earlyResp*FrameRate)):((leadtime*FrameRate)+(midResp*FrameRate)))))...
                                    /length((((leadtime*FrameRate)+(earlyResp*FrameRate)):((leadtime*FrameRate)+(midResp*FrameRate))));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.late(ind)=sum(lowRefl(ind,(((leadtime*FrameRate)+(peakWin*FrameRate)):((leadtime*FrameRate)+(lateResp*FrameRate)))))...
                                    /length((((leadtime*FrameRate)+(midResp*FrameRate)):((leadtime*FrameRate)+(lateResp*FrameRate))));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.core(ind)=sum(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(midResp*FrameRate)))))...
                                    /length(((leadtime*FrameRate):((leadtime*FrameRate)+(midResp*FrameRate))));
                                GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.total(ind)=sum(lowRefl(ind,((leadtime*FrameRate):((leadtime*FrameRate)+(lateResp*FrameRate)))))...
                                    /length(((leadtime*FrameRate):((leadtime*FrameRate)+(lateResp*FrameRate))));
                            end
                        end
                        %% Get max peaks
                        TestPos=GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak;
                        TestNeg=GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak;
                        TestNeg(isnan(TestNeg))=0;
                        TestPeaks=TestPos>abs(TestNeg);
                        for p=1:size(TestPeaks,2)
                            if TestPeaks(p)==1
                               GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllPeak(p)=TestPos(p);
                               GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllpkPT(p)=GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT(p);
                            else
                              GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllPeak(p)=TestNeg(p);
                              GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllpkPT(p)=GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT(p);
                            end
                        end
                                
                        %% Average Individual Animal Data
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).animalNum=size(tempRefl,1)-size(find(isnan(tempRefl(:,1))),1);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.Avg_HbT=mean(lowRefl,1,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).IOS.StD_HbT=std(lowRefl,0,1,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllPeak=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllPeak,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllpkPT=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.AllpkPT,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosPeak,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PospkPT,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.PosFWHM,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegPeak,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegpkPT,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).PeakParams.NegFWHM,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.early=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.early,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.mid=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.mid,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.late=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.late,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.core=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.core,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.total=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.total,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.early_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.early,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.mid_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.mid,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.late_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.late,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.core_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.core,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.total_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).AUC.total,'omitnan');
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinVal=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinVal);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.Min_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinVal);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinLoc=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinLoc);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinLoc_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MinLoc);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxVal=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxVal);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.Max_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxVal);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxLoc=mean(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxLoc);
                        GroupedData.PopulationData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxLoc_std=std(GroupedData.AnimalData.(StimTypes{stimNum}).(ROInames{ROInum}).(behaviorTypes{beNum}).MinMax.MaxLoc);
                        tempRefl=[];
                        lowRefl=[];
                    end
                end
            end
        end
        save([Animal_Age '_GroupedDataAnimal'],'GroupedData','-v7.3');
        clearvars -except directoryContents folderNum q
    end
end
end