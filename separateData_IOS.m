function [positiveHbT_data,negativeHbT_data]=separateData_IOS(AllData)

Run_Type={'Volitional','Evoked'};
for trigType=1:size(Run_Type,2)
    if ~isempty(AllData.(Run_Type{trigType}).All.IOS.forepaw.ChunkRefl)
        for eventNum=1:2
            if eventNum==1
                EventFlag=find(AllData.(Run_Type{trigType}).All.NegFlag==0);%Negative HbT response ==0
                EventLabel='negativeHbT';
            else
                EventFlag=find(AllData.(Run_Type{trigType}).All.NegFlag==1);
                EventLabel='positiveHbT';
            end  
            subfields=fieldnames(AllData.(Run_Type{trigType}).All);
            for fieldnum=1:size(subfields,1)
                if isstruct(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}))
                    finFields=fieldnames(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}));
                    for finnum=1:size(finFields,1)
                        if isstruct(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}))
                            termfields=fieldnames(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}));
                            for termNum=1:size(termfields,1)
                                if iscell(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}))
                                    for cellNum=1:length(EventFlag)
                                        Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}){cellNum}=AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}){EventFlag(cellNum)};
                                    end
                                else
                                    if size(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),3)>1
                                        GroupedLabel=['Avg' termfields{termNum}];
                                        Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})=AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})(:,:,EventFlag);
                                        Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(GroupedLabel)=mean(Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),3);
                                    else
                                        GroupedLabel=['Avg' termfields{termNum}];
                                        Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})=AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum})(EventFlag,:);
                                        Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(GroupedLabel)=mean(Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}).(termfields{termNum}),1);
                                    end
                                end
                            end
                        else
                            if iscell(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}))
                                for cellNum=1:length(EventFlag)
                                    Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}){cellNum}=AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}){EventFlag(cellNum)};
                                end
                            else
                                GroupedLabel=['Avg' finFields{finnum}];
                                if contains(GroupedLabel,'Spec')
%                                 if size(AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum}),3)>1
                                    Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum})=AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum})(:,:,EventFlag);
                                    Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(GroupedLabel)=median(Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}),3);
                                else
                                    Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum})=AllData.(Run_Type{trigType}).All.(subfields{fieldnum}).(finFields{finnum})(EventFlag,:);
                                    Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(GroupedLabel)=mean(Temp_Data.(EventLabel).(Run_Type{trigType}).(subfields{fieldnum}).(finFields{finnum}),1);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
positiveHbT_data=Temp_Data.positiveHbT;
negativeHbT_data=Temp_Data.negativeHbT;
end