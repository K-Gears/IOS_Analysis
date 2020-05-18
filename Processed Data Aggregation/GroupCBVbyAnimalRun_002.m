function GroupCBVbyAnimalRun_002(~)
%Function to aggregate hemeodynamic data by animal and day.
%Written by Kyle Gheres, Nov 2014
% Identify Files
filename = uigetfile('*_procdata.mat','MultiSelect', 'on');
prompt='what age are these animals?';
Animal_Age=input(prompt,'s');
prompt='Were these injected?';
Injected=input(prompt,'s');
% Control for single file instead of a list
if iscell(filename) == 0
    filename = {filename};
end
for fil = 1:length(filename)
    close all;
    indfile = filename{fil};
    load(indfile);
    Fields=fieldnames(Proc_Data);
    m=0;
    for k=1:numel(Fields)
        if isempty(Proc_Data.(Fields{k}))==0
            m=m+1;
            Run_Type{m}=Fields{k};
        end
    end
    for k=1:numel(Run_Type)
        if isfield(Proc_Data.(Run_Type{k}),'Run_Events')==1
            Run_length=fieldnames(Proc_Data.(Run_Type{k}).Run_Events);
            for n=1:numel(Run_length)
                if isfield(Proc_Data.(Run_Type{k}).Run_Events,Run_length{n})==1
                    if fil==1
                        GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl(fil,:)=Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp;
                        GroupedData.(Run_Type{k}).(Run_length{n}).Single_AUC(fil)=(trapz(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2)/size(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2));
                    else
                        if isfield(GroupedData.(Run_Type{k}),Run_length{n})==1
                            count=size(GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl,1)+1;
                            if lt(size(GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl,2),size(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp,2))
                                GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl(count,:)=Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(1:(size(GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl,2)));
                                GroupedData.(Run_Type{k}).(Run_length{n}).Single_AUC(count)=(trapz(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2)/size(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2));
                            else
                                GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl=GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl(:,(1:size(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp,2)));
                                GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl(count,:)=Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp;
                                GroupedData.(Run_Type{k}).(Run_length{n}).Single_AUC(count)=(trapz(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2)/(size(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2)));
                            end
                        else
                            GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl(1,:)=Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp;
                            GroupedData.(Run_Type{k}).(Run_length{n}).Single_AUC(1)=(trapz(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2)/(size(Proc_Data.(Run_Type{k}).Run_Events.(Run_length{n}).Avg_Response_fp(60:end),2)));
                        end
                    end
                end
            end
        end
    end
end
Run_Type=fieldnames(GroupedData);
for k=1:numel(Run_Type)
    Run_length=fieldnames(GroupedData.(Run_Type{k}));
    for n=1:numel(Run_length)
        GroupedData.(Run_Type{k}).(Run_length{n}).Group_Avg_Refl=mean(GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl,1);
        GroupedData.(Run_Type{k}).(Run_length{n}).Group_Std_Refl=std(GroupedData.(Run_Type{k}).(Run_length{n}).Avg_Forepaw_Refl,0,1);
        GroupedData.(Run_Type{k}).(Run_length{n}).Group_Avg_AUC=(trapz(GroupedData.(Run_Type{k}).(Run_length{n}).Group_Avg_Refl(60:end),2)/size(GroupedData.(Run_Type{k}).(Run_length{n}).Group_Avg_Refl(60:end),2));
    end
end
save([Animal_Age '_' Injected '_GroupedRunDataAnimal'],'GroupedData');
end