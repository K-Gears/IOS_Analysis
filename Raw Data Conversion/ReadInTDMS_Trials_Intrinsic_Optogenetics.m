% Adapted from Patrick Drew Code. 4/12/13 - Aaron Winder, Drew Lab, PSU
%Adapted for extracting LED stimulus times 4/28/19- Kyle Gheres, Drew Lab


function [TDMSFile]=ReadInTDMS_Trials_Intrinsic_Optogenetics(filename,usefieldnames)
% this function reads in .tdms files with convertTDMS and organizes the
% structure into 

%load in the struct
[TempStruct,ConvertVerion]=convertTDMS(0,filename);
%ADataStruct=struct('date_read_in', date);
%reorganize struct

TDMSFile.Data=zeros(length(TempStruct.Data.MeasuredData),length(TempStruct.Data.MeasuredData(1).Data));
for k=1:length(TempStruct.Data.MeasuredData)
    hold_data=TempStruct.Data.MeasuredData(1,k).Data;
TDMSFile.Data(k,:)=hold_data;
end

% Add trial properties to the structure <- FUTURE, DISPLAY TABLE OF TRIAL
% PROPERTIES
if usefieldnames==1
    Datafields=fieldnames(TempStruct.Data.Root);
    for fieldnum=1:size(Datafields,1)
        TDMSFile.(Datafields{fieldnum})=TempStruct.Data.Root.(Datafields{fieldnum});
    end
else
TDMSFile.Experimenter=TempStruct.Data.Root.Experimenter;
TDMSFile.Animal_Name=TempStruct.Data.Root.Animal_Name;
TDMSFile.Hemisphere_Recorded=TempStruct.Data.Root.Hemisphere_Recorded;
TDMSFile.Session_Number=TempStruct.Data.Root.Session_Number;
TDMSFile.Dalsa_Frame_Rate=str2num(TempStruct.Data.Root.Dalsa_Frame_Rate);
%TDMSFile.Basler_Frame_Rate=str2num(TempStruct.Data.Root.Basler_Frame_Rate);
TDMSFile.Analog_Sampling_Freq=str2num(TempStruct.Data.Root.Analog_Sampling_Freq);
TDMSFile.WebCam_Frame_Rate=str2num(TempStruct.Data.Root.WebCam_Frame_Rate);
TDMSFile.Isoflurane_Time=str2num(TempStruct.Data.Root.Isoflurane_Time);
TDMSFile.Flash_Trial=TempStruct.Data.Root.Flash_Trial;
TDMSFile.Stim_Start_Offset=TempStruct.Data.Root.Stim_offset;
TDMSFile.Time_Between_Stim=TempStruct.Data.Root.Time_Between_Stim;
TDMSFile.Laser_Duration=TempStruct.Data.Root.LED_Duration;
TDMSFile.Laser_Frequency=TempStruct.Data.Root.LED_Frequency;
TDMSFile.Laser_Duty_Cycle=TempStruct.Data.Root.LED_Duty_Cycle;
TDMSFile.Solenoid_Duration=TempStruct.Data.Root.Solenoid_duration;
TDMSFile.Solenoid_Frequency=TempStruct.Data.Root.Solenoid_Frequency;
TDMSFile.Solenoid_Duty_Cycle=TempStruct.Data.Root.Solenoid_duty_cycle;
end