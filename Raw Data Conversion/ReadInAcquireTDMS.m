function [ADataStruct]=ReadInAcquireTDMS(filename)
% this function reads in .tdms files with convertTDMS and organizes the
% structure into 

%load in the struct
[TempStruct,ConvertVerion]=convertTDMS(0,filename);
%ADataStruct=struct('date_read_in', date);
%reorganize struct

ADataStruct.Data=zeros(length(TempStruct.Data.MeasuredData),length(TempStruct.Data.MeasuredData(1).Data));
for k=1:length(TempStruct.Data.MeasuredData)
    hold_data=TempStruct.Data.MeasuredData(1,k).Data;
ADataStruct.Data(k,:)=hold_data;
end

%ADataStruct.Sampling_freq=str2num(TempStruct.Data.Root.Sampling_freq);
%ADataStruct.Experimenter=TempStruct.Data.Root.Experimenter;
%ADataStruct.Animal=TempStruct.Data.Root.Animal;
%ADataStruct.Ball_gain_V_per_Rad=str2num(TempStruct.Data.Root.Sampling_freq);
%ADataStruct.Ephys_gain=str2num(TempStruct.Data.Root.Ephys_gain);
