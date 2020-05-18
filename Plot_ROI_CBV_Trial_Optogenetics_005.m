function Plot_ROI_CBV_Trial_Optogenetics_005(filename,RawData,imp_bin,T_run,new_T_run,velocity,weightedcoeffHbO,weightedcoeffHbR,weightedcoeffHbT)
%% Notes
%Output data is plotted as change in concentration of total hemoglobin HbT.
ColorMap=brewermap(22,'RdBu');
ROI_name=fieldnames(RawData.IOS);
reOrderROI=strcmpi(ROI_name,'Pixelwise')==1;
ROI_name{reOrderROI}=ROI_name{length(ROI_name)};
ROI_name{length(ROI_name)}='Pixelwise';
coeffType={'HbO','HbR','HbT'};
orderROI=1:length(ROI_name);
orderROI(reOrderROI)=length(ROI_name);
orderROI(length(ROI_name))=find(reOrderROI==1);
Isrunning=isempty(new_T_run);
Binarized_Running(1:size(imp_bin,1))=0;
Chunked_Running(1:size(imp_bin,1))=0;
BallAcc=diff(velocity);
BallAcc(BallAcc<(1e-5))=0;
if Isrunning==0
    RunStart=new_T_run(1,:); %Find time(s) when animal starts running
    RunStop=new_T_run(2,:); %Find time(s) when animal stops running
    for L=1:length(RunStart)-1 %Determine time between running events
        try
            Stopped(L)=RunStart(L+1)-RunStop(L);
        catch
            fprintf(1,'animal ran through whole trial\n')
        end
    end
    clear Rest;
    try
        Rest=find(Stopped==max(Stopped)); %Find maximum still period
    catch
        fprintf(1,'No rest in trial\n')
    end
    clear Stopped;
    for k=1:size(ROI_name,1)
        try
            Normalizationconstant.(ROI_name{k})=mean(RawData.IOS.(ROI_name{k}).CBVrefl((RunStop(Rest):RunStart(Rest+1))));
        catch
            fprintf(1,'Used trial ave as normalization constant\n')
            Normalizationconstant.(ROI_name{k})= mean(RawData.IOS.(ROI_name{k}).CBVrefl);
        end
    end
    fprintf(1,'filter and plot whole trial reflectance and locomotion\n')
    
    
    for k=1:size(T_run,2)
        Binarized_Running(T_run(1,k):T_run(2,k))=1;
    end
    for q=1:size(new_T_run,2)
        Chunked_Running(new_T_run(1,q):new_T_run(2,q))=1;
    end
else
    for k=1:size(ROI_name,1)
    Normalizationconstant.(ROI_name{k})= mean(RawData.IOS.(ROI_name{k}).CBVrefl);
    end
    Binarized_Running=imp_bin';
end
[Order,lowpass]=butter(3,(1/(RawData.dal_fr*0.5)), 'low');
for HbType=1:length(coeffType)
    if HbType==1
        HbCoeff=weightedcoeffHbO;
    elseif HbType==2
        HbCoeff=weightedcoeffHbR;
    else
        HbCoeff=weightedcoeffHbT;
    end
    %% Correct ROI for Opto Stim
    if isfield(RawData.IOS,'Optogenetics')
        FlashCatch=abs(diff(RawData.IOS.Optogenetics.CBVrefl));
        for k=1:(size(ROI_name,1)-1)
            if max(FlashCatch)>=20
                HoldRefl=RawData.IOS.(ROI_name{k}).CBVrefl;
                %HoldRefl=RawData.IOS.(ROI_name{k}).CBVrefl-Normalizationconstant.(ROI_name{k});
                %Flash_Points=HoldRefl>=(3*std(RawData.IOS.(ROI_name{k}).CBVrefl));
                if k==1
                Flash_Points=find(FlashCatch>=8)+1;
                First_Flash=Flash_Points(diff(Flash_Points)==1);
                end
                HoldRefl(First_Flash)=NaN;
                %Trial_Length=1:length(RawData.IOS.(ROI_name{k}).CBVrefl);
                %Interp_Data=interp1(Trial_Length,HoldRefl,Flash_Points,'spline');
                [Interp_Data]=fillmissing(HoldRefl,'spline');
                %Interp_Refl=HoldRefl;
                Interp_Refl=Interp_Data;
                % for thepoint=1:length(Interp_Data)
                %     Interp_Refl(Flash_Points(thepoint))=Interp_Data(thepoint);
                % end
                Normalizedrefl.(ROI_name{k}).(coeffType{HbType})=(HbCoeff.*log(Interp_Refl/Normalizationconstant.(ROI_name{k})))*1e6;
            else
                Normalizedrefl.(ROI_name{k}).(coeffType{HbType})=(HbCoeff.*log(RawData.IOS.(ROI_name{k}).CBVrefl/Normalizationconstant.(ROI_name{k})))*1e6;
            end
        end
    else
        for k=1:(size(ROI_name,1)-1)
            Normalizedrefl.(ROI_name{k}).(coeffType{HbType})=(HbCoeff.*log(RawData.IOS.(ROI_name{k}).CBVrefl/Normalizationconstant.(ROI_name{k})))*1e6;
        end
    end
    %% Lowpass Filter IOS data
    
    for k=1:(size(ROI_name,1)-1)
        Normalizedrefl.(ROI_name{k}).(coeffType{HbType})=filtfilt(Order,lowpass,Normalizedrefl.(ROI_name{k}).(coeffType{HbType}));
    end
end
%% Resume Normal Plotting

Camera_time(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))=(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))/RawData.dal_fr; %Plots x-axis in seconds 4-23-16 KG
Max_Refl=[max(Normalizedrefl.(ROI_name{1}).(coeffType{1})),max(Normalizedrefl.(ROI_name{2}).(coeffType{1}))];
Binarized_Running(Binarized_Running==1)=max(Max_Refl)+.25*max(Max_Refl); %Binarized locomotion events will be scatter plotted +1 above maximal value of Reflectance changes 4-23-16 KG
Binarized_Running(Binarized_Running==0)=NaN; %Eliminates scatter plotting of 0 values (Still time) of locomotion 4-23-16 KG
Chunked_Running(Chunked_Running==1)=max(Max_Refl)+.25*max(Max_Refl);
Chunked_Running(Chunked_Running==0)=NaN;
if any(RawData.Sol)==1
    Left_Puff=unique(floor(find(RawData.Sol==1)/RawData.an_fs));
    Left_Puff_Time(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))=0;
    Left_Puff_Time((Left_Puff*RawData.dal_fr))=max(Max_Refl)+.1*max(Max_Refl);
    Left_Puff_Time(Left_Puff_Time==0)=NaN;
    Right_Puff=unique(floor(find(RawData.Sol==3)/RawData.an_fs));
    Right_Puff_Time(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))=0;
    Right_Puff_Time((Right_Puff*RawData.dal_fr))=max(Max_Refl)+.1*max(Max_Refl);
    Right_Puff_Time(Right_Puff_Time==0)=NaN;
    Noise_Puff=unique(floor(find(RawData.Sol==2)/RawData.an_fs));
    Noise_Puff_Time(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))=0;
    Noise_Puff_Time((Noise_Puff*RawData.dal_fr))=max(Max_Refl)+.1*max(Max_Refl);
    Noise_Puff_Time(Noise_Puff_Time==0)=NaN;
end
if any(round(RawData.LED,0))==1
    State_Change=diff(round(RawData.LED,0));
    Laser_Switch=unique(find(State_Change==min(State_Change))/RawData.an_fs);
    Laser_Time(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))=0;
    Laser_Time(round(Laser_Switch*RawData.dal_fr,0))=max(Max_Refl)+.1*max(Max_Refl);
    Laser_Time(Laser_Time==0)=NaN;
else
    Laser_Time(1:length(RawData.IOS.(ROI_name{1}).CBVrefl))=NaN;
end
Lgnd_cnt=1;
figure(22)
h1=subplot(4,1,(2:4));
hold on;
if ge(length(Camera_time),length(Binarized_Running))
    scatter(Camera_time(1:length(Binarized_Running)),Binarized_Running,'filled','k');
    LgndTxt{Lgnd_cnt}='Binarized Running';
    Lgnd_cnt=Lgnd_cnt+1;
else
    scatter(Camera_time,Binarized_Running(1:length(Camera_time)),'filled','k');
    LgndTxt{Lgnd_cnt}='Binarized Running';
    Lgnd_cnt=Lgnd_cnt+1;
end
if exist('new_T_run','var')==1
    if ge(length(Camera_time),length(Binarized_Running))
        scatter(Camera_time(1:length(Binarized_Running)),Chunked_Running(1:length(Binarized_Running)),'filled','b');
        LgndTxt{Lgnd_cnt}='Analyzed Running';
        Lgnd_cnt=Lgnd_cnt+1;
    else
        scatter(Camera_time,Chunked_Running(1:length(Camera_time)),'filled','b');
        LgndTxt{Lgnd_cnt}='Analyzed Running';
        Lgnd_cnt=Lgnd_cnt+1;
    end
end
if any(RawData.Sol)==1
    if ge(length(Camera_time),length(Binarized_Running))
        if isnan(max(Left_Puff_Time))==0
            scatter(Camera_time(1:length(Binarized_Running)),Left_Puff_Time(1:length(Binarized_Running)), 'filled','g');
            LgndTxt{Lgnd_cnt}='Left Whisker Stim';
            Lgnd_cnt=Lgnd_cnt+1;
        end
        if isnan(max(Right_Puff_Time))==0
            scatter(Camera_time(1:length(Binarized_Running)),Right_Puff_Time(1:length(Binarized_Running)), 'filled','m');
            LgndTxt{Lgnd_cnt}='Right Whisker Stim';
            Lgnd_cnt=Lgnd_cnt+1;
        end
        if isnan(max(Noise_Puff_Time))==0
            scatter(Camera_time(1:length(Binarized_Running)),Noise_Puff_Time(1:length(Binarized_Running)),'filled','r');
            LgndTxt{Lgnd_cnt}='Noise Stim';
            Lgnd_cnt=Lgnd_cnt+1;
        end
    else
        if isnan(max(Left_Puff_Time))==0
            scatter(Camera_time,Left_Puff_Time(1:length(Camera_time)), 'filled','g');
            LgndTxt{Lgnd_cnt}='Left Whisker Stim';
            Lgnd_cnt=Lgnd_cnt+1;
        end
        if isnan(max(Right_Puff_Time))==0
            scatter(Camera_time,Right_Puff_Time(1:length(Camera_time)), 'filled','m');
            LgndTxt{Lgnd_cnt}='Right Whisker Stim';
            Lgnd_cnt=Lgnd_cnt+1;
        end
        if isnan(max(Noise_Puff_Time))==0
            scatter(Camera_time,Noise_Puff_Time(1:length(Camera_time)),'filled','r');
            LgndTxt{Lgnd_cnt}='Noise Stim';
            Lgnd_cnt=Lgnd_cnt+1;
        end
    end
end
if any(round(RawData.LED,0))==1
    if ge(length(Camera_time),length(Binarized_Running))
        scatter(Camera_time(1:length(Binarized_Running)),Laser_Time(1:length(Binarized_Running)),'filled','c');
        LgndTxt{Lgnd_cnt}='Opto Stim';
        Lgnd_cnt=Lgnd_cnt+1;
    else
        scatter(Camera_time,Laser_Time(1:length(Camera_time)),'filled','c');
        LgndTxt{Lgnd_cnt}='Opto Stim';
        Lgnd_cnt=Lgnd_cnt+1;
    end
end
ColorOrder={'b','r','g','c','m'};
for k=1:(size(ROI_name,1)-1)
    for q=1:length(coeffType)
        if q==1
            PlotColor=ColorMap((4+k),:);
        elseif q==2
            PlotColor=ColorMap((20-k),:);
        else
            NewMap=brewermap(11,'PRGn');
            PlotColor=NewMap((11-k),:);
        end
        plot(Camera_time,Normalizedrefl.(ROI_name{k}).(coeffType{q}),'Color',PlotColor);%ColorOrder{k});
        LgndTxt{Lgnd_cnt}=strrep([ROI_name{k} ' ' coeffType{q}],'_',' ');
        Lgnd_cnt=Lgnd_cnt+1;
    end
end
xlim([0 300]);
filnm=strrep(filename,'_', ' ');
title(['Whole Trial Normalized Reflectance ' filnm],'FontSize',14,'FontWeight','bold','FontName','Arial');
xlabel('Time (s)','FontSize',10,'FontWeight','bold','FontName','Arial');
ylabel('\Delta \muM HbT','FontSize',10,'FontWeight','bold','FontName','Arial');
legend(LgndTxt, 'Location','southeast',...
    'Orientation','vertical','FontSize',8,'FontWeight','bold','FontName','Arial');
h2=subplot(4,1,1);
plot(((1:length(BallAcc))/RawData.an_fs),abs(BallAcc),'k');
if max(abs(BallAcc))<1e-3
    ylim([0 5e-4]);
else
    ylim([0 (max(abs(BallAcc))+(max(abs(BallAcc))*0.1))]);
end
xlim([0 300]);
linkaxes([h1,h2],'x');
savefig([filename '_Normalized Reflectance Per Session W_locomotion']);
% lgnd_cnt=1;
% figure(33); hold on;
% for k=1:(size(ROI_name,1)-1)
% plot(Camera_time,(Normalizedrefl.(ROI_name{k}).HbO-Normalizedrefl.(ROI_name{k}).HbR),ColorOrder{k});
% lgndTxt{lgnd_cnt}=strrep([ROI_name{k} ' normalized residual HbO vs HbR'],'_',' ');
% lgnd_cnt=lgnd_cnt+1;
% end
% title(['Whole Trial Normalized difference between HbO and HbR extinction values ' filnm],'FontSize',14,'FontWeight','bold','FontName','Arial');
% xlabel('Time (s)','FontSize',10,'FontWeight','bold','FontName','Arial');
% ylabel('\Delta HbO-HbR normalized to HbO (%)','FontSize',10,'FontWeight','bold','FontName','Arial');
% legend(lgndTxt, 'Location','southeast',...
%     'Orientation','vertical','FontSize',8,'FontWeight','bold','FontName','Arial');
% savefig([filename '_Normalized difference between HbO and HbR extinction coefficients']);
close all;
end