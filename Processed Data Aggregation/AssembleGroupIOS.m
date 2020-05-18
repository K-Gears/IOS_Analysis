function [SortedData]=AssembleGroupIOS(Drivename)
filelist=dir(fullfile(Drivename,'**','*GroupedDataAnimal.mat'));
for filnum=1:size(filelist,1)
    load([filelist(filnum).folder '\' filelist(filnum).name]);
    StimTypes=fieldnames(GroupedData.PopulationData);
    SortedData.AnimalAge{filnum}=filelist(filnum).folder(4:end);
    for stimNum=1:size(StimTypes,1)
        ROIname=fieldnames(GroupedData.PopulationData.(StimTypes{stimNum}));
        for ROInum=1:size(ROIname,1)
            behaveName=fieldnames(GroupedData.PopulationData.(StimTypes{stimNum}).(ROIname{ROInum}));
            for behaveNum=1:size(behaveName,1)
                dataType=fieldnames(GroupedData.PopulationData.(StimTypes{stimNum}).(ROIname{ROInum}).(behaveName{behaveNum}));
                for dataNum=2:size(dataType)
                    dataFields=fieldnames(GroupedData.PopulationData.(StimTypes{stimNum}).(ROIname{ROInum}).(behaveName{behaveNum}).(dataType{dataNum}));
                    for fieldNums=1:size(dataFields,1)
                        SortedData.Population.(StimTypes{stimNum}).(ROIname{ROInum}).(behaveName{behaveNum}).(dataType{dataNum}).(dataFields{fieldNums})(:,filnum)=GroupedData.PopulationData.(StimTypes{stimNum}).(ROIname{ROInum}).(behaveName{behaveNum}).(dataType{dataNum}).(dataFields{fieldNums});
                        if ~strcmpi(dataType{dataNum},'IOS')
                            if isempty(strfind(dataFields{fieldNums},'std'))
                            SortedData.Animals.(StimTypes{stimNum}).(ROIname{ROInum}).(behaveName{behaveNum}).(dataType{dataNum}).(dataFields{fieldNums}){filnum}=GroupedData.AnimalData.(StimTypes{stimNum}).(ROIname{ROInum}).(behaveName{behaveNum}).(dataType{dataNum}).(dataFields{fieldNums});
                            end
                        end
                    end
                end
            end
        end
    end
end
save('DevelopmentNVC_IOS_Fig1','SortedData','-v7.3');
end
            
    