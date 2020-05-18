function ProcessTrials_Intrinsic_Optogenetics_001(~)
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Mar 2013
%   Version 1
%   Revised for multiple ROIs by Kyle Gheres, Drew Lab, MCIBS,
%   Penn State University, Feb 2017
%   Revised to allow for processing of optogenetics stimuli Kyle Gheres Mar
%   2019
%
%   SUMMARY: Converts raw data from whisker trials into usable arrays and
%               saves them as a file in the current folder
%_______________________________________________________________
%   INPUTS:                     None
%_______________________________________________________________
%   OUTPUTS:                    Output is a file saved to the current
%                               folder
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%       - Declares dalsa frequency based on the status of "FlashTrial"
%       field from the .tdms file.
% -Includes changes to collect rotary encoder ball velocity (vBall)
%-Modification of RawData.Sol for use with Neonate Whisker Trial v5+ (all
%solenoids collected as single array with each channel having unique
%amplitude
%-Includes ability to create multiple ROI for independent analysis as well
%as pixel wise calculations within ROI
%_______________________________________________________________

%% Load computer specific variables
%[Q] = DetectMachine;

% Identify Files
% WorkingFolder=cd;
% Datafolders=dir(WorkingFolder);
% count=1;
% for foldernum=3:size(Datafolders,1)
%     cd([Datafolders(foldernum).folder '\' Datafolders(foldernum).name]);
%     findfile=dir('*dalsa.bin');
%     filenames{count}=findfile.name;
%     count=count+1;
% end
% cd(WorkingFolder);
filenames = uigetfile('*_dalsa.bin','MultiSelect', 'on');
RecordNeuro=1;%1 if yes, 0 if no
usefieldnames=1;%1 if yes, 0 if no
% Control for single file instead of a list
if iscell(filenames) == 0
    fileid = filenames;
    filenames = 1;
end

last_file = [];
for fil = 1:length(filenames)
    close all;
    Proc_Vars.fil=fil;
    % Adapt to list or single file
    if iscell(filenames) == 1
        indfile = filenames{fil};
    else
        indfile = fileid;
    end
    %% Import .tdms data (All channels);
    trialdata = ReadInTDMS_Trials_Intrinsic_Optogenetics([indfile(1:17) '.tdms'],usefieldnames);          %Index 1:17 is specific to this data type may require change in future
    animal = trialdata.Animal_Name;
    if strcmp(trialdata.Flash_Trial,'y')
        if ischar(trialdata.Dalsa_Frame_Rate)
        dal_fr=str2double(trialdata.Dalsa_Frame_Rate)/2;
        else
        dal_fr = trialdata.Dalsa_Frame_Rate/2;
        end
        display(['Dalsa Frame Rate has been divided by 2. Current Frame Rate is: ' num2str(dal_fr)]); %<- Temporary
    else
        if ischar(trialdata.Dalsa_Frame_Rate)
            dal_fr=str2double(trialdata.Dalsa_Frame_Rate);
        else
            dal_fr = trialdata.Dalsa_Frame_Rate;
        end
    end
    if ischar(trialdata.Analog_Sampling_Freq)
        an_fs=str2double(trialdata.Analog_Sampling_Freq);
    else
    an_fs = trialdata.Analog_Sampling_Freq;
    end
    hem = trialdata.Hemisphere_Recorded;
    curr_file = [animal trialdata.Session_Number];
    if exist([animal '_' hem '_' indfile(1:17) '_rawdata.mat'],'file') == 2
        disp('File already exists, continuing...'); continue;
    end
    vBall=trialdata.Data(1,:);%voltage signal from rotary encoder attached to ball 3-3-15 KG
    
    if strcmp(curr_file,last_file) == 1
        match = 1;
    else
        match = 0;
    end
    
    %% Import Dalsa images
    % Insert your Dalsa binary image to pixelwise and ROI averaged
    % conversion script here
    [CBVrefl,FlashData,Proc_Vars] = Bin2Intensity_MultiROI(indfile,animal,hem,trialdata.Flash_Trial,match,Proc_Vars);
   
    %% Create RawData structure for each trial for additional post processing
    ROI_names=fieldnames(CBVrefl); %Insert Structure containing ROI data here
    for k=1:size(ROI_names,1)
%         if strcmpi(ROI_names{k},'Pixelwise')==1
%             RawData.IOS.(ROI_names{k}).CBVrefl=CBVrefl.(ROI_names{k}).Refl;
%             RawData.IOS.(ROI_names{k}).PixelMap=CBVrefl.(ROI_name).Pixel_Map;
        %else
            RawData.IOS.(ROI_names{k}).CBVrefl=CBVrefl.(ROI_names{k}).Refl;  
            RawData.IOS.(ROI_names{k}).PixelMap=CBVrefl.(ROI_names{k}).Pixel_Map;
       % end
    end
    
    if  RecordNeuro==1
    
    RawData.FlashData = FlashData;
    RawData.Sol = trialdata.Data(3,:);
    RawData.LED=trialdata.Data(2,:);
    RawData.Neuro=trialdata.Data(4,:);
    RawData.LFP=trialdata.Data(5,:);
    RawData.MUA=trialdata.Data(6,:);
    RawData.dal_fr = dal_fr;
    RawData.an_fs = an_fs;
    RawData.Session = trialdata.Session_Number;
    RawData.FlashTrial = trialdata.Flash_Trial;
    RawData.WC_fr = trialdata.WebCam_Frame_Rate;
    RawData.vBall=vBall; %voltage signal from ball to be converted to velocity 3-3-15 KG
    RawData.AcquistionParams.Stim_Start=str2double(trialdata.Stim_offset); %Offset in seconds from start of trial for first stimulus trialdata.Stim_Start_Offse
    RawData.AcquistionParams.Time_Between_Stim=str2double(trialdata.Time_Between_Stim);%trialdata.Time_Between_Stim
    RawData.AcquistionParams.Laser_Duration=str2double(trialdata.LED_Duration);%trialdata.Laser_Duration
    RawData.AcquistionParams.Laser_Frequency=str2double(trialdata.LED_Frequency);
    RawData.AcquistionParams.Laser_Duty_Cycle=str2double(trialdata.LED_Duty_Cycle);
    RawData.AcquistionParams.Solenoid_Duration=str2double(trialdata.Solenoid_duration);
    RawData.AcquistionParams.Solenoid_Frequency=str2double(trialdata.Solenoid_Frequency);
    RawData.AcquistionParams.Solenoid_Duty_Cycle=str2double(trialdata.Solenoid_duty_cycle);
    else
    mergeSol=1;
    RawData.FlashData = FlashData;
    if mergeSol==1
        expectedlength=an_fs*(length(CBVrefl.(ROI_names{1}).Refl)/dal_fr);
        for solchan=2:size(trialdata.Data,1)
            channelnum=solchan-1;
            if length(trialdata.Data(solchan,:))<expectedlength
                trialdata.Data(solchan,(length(trialdata.Data(solchan,:)):expectedlength))=0;
            end
            Hold=round(trialdata.Data(solchan,(1:expectedlength)),0);
            Hold(Hold==5)=channelnum;
            Channeldata(channelnum,:)=Hold;
        end
    RawData.Sol=max(Channeldata,[],1);    
    else
    RawData.Sol = trialdata.Data(2,:);
    end
    RawData.dal_fr = dal_fr;
    RawData.an_fs = an_fs;
    RawData.Session = trialdata.Session_Number;
    RawData.FlashTrial = trialdata.Flash_Trial;
    RawData.WC_fr = trialdata.WebCam_Frame_Rate;
    RawData.vBall=vBall; %voltage signal from ball to be converted to velocity 3-3-15 KG
    RawData.AcquistionParams.Stim_Start=10; %Offset in seconds from start of trial for first stimulus
    RawData.AcquistionParams.Time_Between_Stim=30;
    RawData.AcquistionParams.Solenoid_Duration=1;
    RawData.AcquistionParams.Solenoid_Frequency=1;
    RawData.AcquistionParams.Solenoid_Duty_Cycle=0.5;
    end
    save([animal '_' hem '_' indfile(1:17) '_rawdata'],'RawData','-v7.3')
    last_file = curr_file;
    prev_flash = trialdata.Flash_Trial;
end