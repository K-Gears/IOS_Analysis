function StimulusVideo_004(rawdatafile)
%Written by:
%Kyle Gheres
%Patrick Drew Lab
%Pennsylvania State University

%Description:
%This script takes binary video files collected from a Dalsa 1M30 camera
%and converts them in to an .avi file with indicators of stimulus periods

%Subscripts:
%ReadDalsaBinary_Matrix

%Usage:
%StimulusVideo_XXX('date_trial_dalsa.bin','animal_hem_date_trial_rawdata.mat');

animal=rawdatafile(1:9);%analogfile(1:6);
Hem=rawdatafile(11:12);%analogfile(8:9);
date=rawdatafile(14:19);%analogfile(11:16);
trial=rawdatafile(21:25);%analogfile(18:22);
misc=rawdatafile(27:30);
load(rawdatafile);
%% Find Optostim Artifacts
ROI_name=fieldnames(RawData.IOS);
ROI_num=find(strcmpi(ROI_name,'Optogenetics'));
FlashCatch=abs(diff(RawData.IOS.(ROI_name{ROI_num}).CBVrefl));
FlashFrames=find(abs(diff(RawData.IOS.(ROI_name{ROI_num}).CBVrefl))>=8)+1;
firstFlash=FlashFrames(diff(FlashFrames)==1);

%% Create Video Objects to write data
videofile=[date '_' trial '_' misc '_dalsa.bin'];
Window_Vid=VideoWriter([animal '_' date '_' trial '_CBVMovie.avi'],'Uncompressed AVI');
Window_Vid.FrameRate=RawData.dal_fr;
open(Window_Vid);

% Refl_Vid=VideoWriter([animal '_' date '_' trial '_ReflMovie.avi'],'Uncompressed AVI');
% Refl_Vid.FrameRate=RawData.dal_fr;
% open(Refl_Vid);
% 
% Ball_Vid=VideoWriter([animal '_' date '_' trial '_BallMovie.avi'],'Uncompressed AVI');
% Ball_Vid.FrameRate=RawData.dal_fr;
% open(Ball_Vid);

%% Open raw data files, extract velocity, reflectance, movie images
[Flash_Clear]=ReadDalsaBinary_Matrix(videofile,256,256);
Flash_Clear=double(Flash_Clear);
FrameTime=(1:size(Flash_Clear,3))/RawData.dal_fr;
%Avg_Frame_Diff=squeeze(mean(mean(diff(theimage,1,3),2),1));
%Flash_Inds=find(Avg_Frame_Diff>10)+1;
%Flash_Clear=theimage;
Flash_Clear(:,:,firstFlash)=NaN;
[Flash_Clear,~]=fillmissing(Flash_Clear,'spline',3);
CorrectedTime=FrameTime;
%CorrectedTime(Flash_Inds)=[];
Avg_Frame=mean(Flash_Clear,3);
Norm_Img=((Flash_Clear-Avg_Frame)./Avg_Frame)*100;
[z,p,k]=butter(3,0.5/15);
[sos,g]=zp2sos(z,p,k);
SmoothImg=permute(filtfilt(sos,g,permute(Norm_Img,[3 2 1])),[3 2 1]);
save([animal '_' date '_' trial '_Normalized_lowpass_movie.mat'],'SmoothImg','-v7.3');
TheMin=min(min(min(SmoothImg)));
ZeroedImg=SmoothImg+abs(TheMin);%sets darkest pixel value of all frames to 0
clear Flash_Clear Norm_Img; %Get rid of extra video variables to clear memory

%% Correct ROI Refl for LED Flashing
% ROI_name=fieldnames(RawData.IOS);
% FlashCatch=abs(diff(RawData.IOS.(ROI_name{ROI_num}).CBVrefl));
% FlashFrames=find(abs(diff(RawData.IOS.(ROI_name{ROI_num}).CBVrefl))>=10)+1;
% firstFlash=FlashFrames(diff(FlashFrames)==1);
if max(FlashCatch)>=10
    HoldRefl=RawData.IOS.(ROI_name{ROI_num}).CBVrefl-mean(RawData.IOS.(ROI_name{ROI_num}).CBVrefl);
%    Flash_Points=HoldRefl>=(3*std(RawData.IOS.(ROI_name{1}).CBVrefl));
    HoldRefl(firstFlash)=NaN;
    [Interp_Data,~]=fillmissing(HoldRefl,'spline');
    %% Correct Barrels ROI for Opto Stim
    Normalizedrefl.(ROI_name{ROI_num})=(((Interp_Data/mean(RawData.IOS.(ROI_name{ROI_num}).CBVrefl)))*100);
else
    Normalizedrefl.(ROI_name{ROI_num})=((((RawData.IOS.(ROI_name{ROI_num}).CBVrefl-mean(RawData.IOS.(ROI_name{ROI_num}).CBVrefl))/mean(RawData.IOS.(ROI_name{ROI_num}).CBVrefl)))*100);
end
%% Normalize and filter plotting data
[ball_b,ball_a]=butter(3,(RawData.dal_fr/(0.5*RawData.an_fs)),'low');
Ball_Velocity=downsample(filtfilt(ball_b,ball_a,RawData.vBall),round(RawData.an_fs/RawData.dal_fr,0));%lowpass filtered and downsampled ball velocity for plotting
[Refl_b,Refl_a]=butter(3,(1/(0.5*RawData.dal_fr)),'low');
CBV_Refl=filtfilt(Refl_b,Refl_a,Normalizedrefl.(ROI_name{ROI_num})); %lowpass filtered normalized reflectance for plotting
Plot_Time=(1:length(CBV_Refl))/RawData.dal_fr;
Plot_Vel(1:length(CBV_Refl))=NaN;
Plot_Vel(end)=0;
Plot_Refl(1:length(CBV_Refl))=NaN;
Plot_Refl(end)=0;
LEDTime=round(find(ceil(RawData.LED)==5)/RawData.an_fs,0);
stimpoint=1;
while stimpoint<length(LEDTime)
    StimWindow=find(LEDTime<=(LEDTime(stimpoint)+RawData.AcquistionParams.Laser_Duration));
    LEDTime(StimWindow((stimpoint+1):end))=[];
    stimpoint=stimpoint+1;
end
SolTime=round(find(ceil(RawData.Sol)>0)/RawData.an_fs,0);
stimpoint=1;
while stimpoint<length(SolTime)
    StimWindow=find(SolTime<=(SolTime(stimpoint)+(RawData.AcquistionParams.Solenoid_Duration*(RawData.AcquistionParams.Solenoid_Duty_Cycle*0.1))));
    SolTime(StimWindow((stimpoint+1):end))=[];
    stimpoint=stimpoint+1;
end
stimnum=1;
solnum=1;

%% Add colored circle to frames that occur during specified stimuli as indicator
% MatMov(1:size(SmoothImg,1),1:size(SmoothImg,1),(1:3),1:size(SmoothImg,3))=0;
for framenum=1:size(SmoothImg,3)
    frametime=CorrectedTime(framenum);
    if stimnum<=length(LEDTime)
        if LEDTime(stimnum)<=frametime
            if frametime<=(LEDTime(stimnum)+RawData.AcquistionParams.Laser_Duration)
                theframe=insertShape(repmat(SmoothImg(:,:,framenum),1,1,3),'FilledCircle',[20 20 15],'Color',[0 169 255]);
                theframe=theframe+abs(TheMin);%Zero after adding indicator to shift indicator intensity with rest of image.
                if frametime==(LEDTime(stimnum)+RawData.AcquistionParams.Laser_Duration)
                    stimnum=stimnum+1;
                end
            end
        else
            if solnum<=length(SolTime)
                if SolTime(solnum)<=frametime
                    if frametime<=(SolTime(solnum)+(RawData.AcquistionParams.Solenoid_Duration*(RawData.AcquistionParams.Solenoid_Duty_Cycle*0.01)))
                        theframe=insertShape(repmat(SmoothImg(:,:,framenum),1,1,3),'FilledCircle',[20 20 15],'Color',[255 0, 0]);
                        theframe=theframe+abs(TheMin);
                    end
                else
                    theframe=repmat(ZeroedImg(:,:,framenum),1,1,3);
                end
                if frametime>=(SolTime(solnum)+(RawData.AcquistionParams.Solenoid_Duration*(RawData.AcquistionParams.Solenoid_Duty_Cycle*0.01)))
                    solnum=solnum+1;
                end
            else
                theframe=repmat(ZeroedImg(:,:,framenum),1,1,3);
            end
            if frametime>(LEDTime(stimnum)+RawData.AcquistionParams.Laser_Duration)
                stimnum=stimnum+1;
            end
        end
    else
        if solnum<=length(SolTime)
            if SolTime(solnum)<=frametime
                if frametime<=(SolTime(solnum)+(RawData.AcquistionParams.Solenoid_Duration*(RawData.AcquistionParams.Solenoid_Duty_Cycle*0.01)))
                    theframe=insertShape(repmat(SmoothImg(:,:,framenum),1,1,3),'FilledCircle',[20 20 15],'Color',[255 0, 0]);
                    theframe=theframe+abs(TheMin);
                end
            else
                theframe=repmat(ZeroedImg(:,:,framenum),1,1,3);
            end
            if frametime>=(SolTime(solnum)+(RawData.AcquistionParams.Solenoid_Duration*(RawData.AcquistionParams.Solenoid_Duty_Cycle*0.01)))
                solnum=solnum+1;
            end
        else
            theframe=repmat(ZeroedImg(:,:,framenum),1,1,3);
        end
    end
%     MatMov(:,:,:,framenum)=theframe-abs(TheMin);
    Imgframe=theframe/255;%Video writer expects values [0 1] and assumes 8 bit dynamic range. Normalized values SHOULD fall in this range.
    Imgframe(Imgframe>1)=1;%correct any out of bounds pixel values
    Imgframe(Imgframe<0)=0;%correct any out of bounds pixel values
    writeVideo(Window_Vid,Imgframe);
%     Plot_Vel(1:framenum)=Ball_Velocity(1:framenum);
%     Plot_Refl(1:framenum)=CBV_Refl(1:framenum);
%     
%     %Plot the CBV reflectance
%     figure(99);plot(Plot_Time,Plot_Refl,'r','LineWidth',2);
%     xlim([0 300]);
%     ylim([-10 10]);
%     title('Window Reflectance');
%     xlabel('Time (sec)');
%     ylabel('Percent Change');
%     Reflframe=getframe(gcf);
%     writeVideo(Refl_Vid,Reflframe);
%     
%     %Plot the ball velocity
%     figure(100);plot(Plot_Time,Plot_Vel,'b','LineWidth',1);
%     xlim([0 300]);
%     ylim([-1 1]);
%     title('Locomotion velocity');
%     xlabel('Time (sec)');
%     ylabel('Percent Change');
%     Ballframe=getframe(gcf);
%     writeVideo(Ball_Vid,Ballframe);
end
% save([animal '_' date '_' trial '_CBVMovie.mat'],'MatMov','-v7.3');
close(Window_Vid);
% close(Ball_Vid);
% close(Refl_Vid);
end