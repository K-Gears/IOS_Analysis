function AgeGroup_Volitional_IOS

driveNames={'E:\','F:\','G:\'};
for driveNum=1:size(driveNames,2)
    folderName='NeonateSleepIndividualAnimals';
    dirName=[driveNames{driveNum} folderName];
    cd(dirName);
    subfolders=dir(dirName);
    subfolders(~[subfolders.isdir])=[];
    tf=ismember({subfolders.name},{'.','..'});
    subfolders(tf)=[];
    for subNum=1:size(subfolders,1)
        if ~strcmpi(subfolders(subNum).name,'070919')
            subDir=[dirName '\' subfolders(subNum).name];
            cd(subDir);
            ageFolders=dir(subDir);
            ageFolders(~[ageFolders.isdir])=[];
            tf=ismember({ageFolders.name},{'.','..'});
            ageFolders(tf)=[];
            for ageNum=1:size(ageFolders.isdir,1)
                cd([ageFolders(ageNum).folder '\' ageFolders(ageNum).name]);
                animalFolders=dir([ageFolders(ageNum).folder '\' ageFolders(ageNum).name]);
                animalFolders(~[animalFolders.isdir])=[];
                tf=ismember({animalFolders.name},{'.','..'});
                animalFolders(tf)=[];
                for anNum=1:size(animalFolders,1)
                    workingDir=[animalFolders(anNum).folder '\' animalFolders(anNum).name];
                    cd(workingDir);
                    findfile=dir('*processedData*.mat');
                    for filnum=1:size(findfile,1)
                        dateVal(filnum)=findfile(filnum).datenum;
                    end
                    newestFil=max(dateVal);
                    NewInd=dateVal==newestFil;
                    filename=findfile(NewInd).name;
                    TheVars=whos('-file',filename);
                    for varNum=1:size(TheVars,1)
                        TheData=load(filename,TheVars(varNum).name);
                        evokeType=fieldnames(TheData);
                        if ~strcmpi(TheVars(varNum).name,'Puff_Data')
                            [SeparateData.positiveHbT_data,SeparateData.negativeHbT_data]=separateData_IOS(TheData.(TheVars(varNum).name));
                            outputNames=fieldnames(SeparateData);%{'positiveHbT_data','negativeHbT_data'};
                            for loopNum=1:size(outputNames,1)
                                behaviorStates=fieldnames(SeparateData.(outputNames{loopNum}));
                                for stateNum=1:size(behaviorStates,1)
                                    dataTypes=fieldnames(SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}));
                                    for dataNum=1:size(dataTypes,1)
                                        datafields=fieldnames(SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}).(dataTypes{dataNum}));
                                        avgFields=strfind(datafields,'Avg');
                                        for fieldnum=1:size(datafields,1)
                                            if isstruct(SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}).(dataTypes{dataNum}).(datafields{fieldnum}))
                                                finalFields=fieldnames(SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}).(dataTypes{dataNum}).(datafields{fieldnum}));
                                                avgFields=strfind(finalFields,'Avg');
                                                for q=1:size(finalFields,1)
                                                    if ~isempty(avgFields{q})
                                                        if anNum==1
                                                            count=1;
                                                        else
                                                            if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}),behaviorStates{stateNum})
                                                                if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}),outputNames{loopNum})
                                                                    if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}),dataTypes{dataNum})
                                                                        if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}),datafields{fieldnum})
                                                                            if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}).(datafields{fieldnum}),finalFields{q})
                                                                                count=size(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}).(datafields{fieldnum}).(finalFields{q}),1)+1;
                                                                            else
                                                                                count=1;
                                                                            end
                                                                        else
                                                                            count=1;
                                                                        end
                                                                    else
                                                                        count=1;
                                                                    end
                                                                else
                                                                    count=1;
                                                                end
                                                            else
                                                                count=1;
                                                            end
                                                        end
                                                        PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}).(datafields{fieldnum}).(finalFields{q})(count,:)=SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}).(dataTypes{dataNum}).(datafields{fieldnum}).(finalFields{q});
                                                    end
                                                end
                                            else
                                                if ~isempty(avgFields{fieldnum})
                                                    if anNum==1
                                                        count=1;
                                                    else
                                                        if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}),behaviorStates{stateNum})
                                                            if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}),outputNames{loopNum})
                                                                if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}),dataTypes{dataNum})
                                                                    if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}),datafields{fieldnum})
                                                                        count=size(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}).(datafields{fieldnum}),3)+1;
                                                                    else
                                                                        count=1;
                                                                    end
                                                                else
                                                                    count=1;
                                                                end
                                                            else
                                                                count=1;
                                                            end
                                                        else
                                                            count=1;
                                                        end
                                                    end
                                                    if ~isnan(SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}).(dataTypes{dataNum}).(datafields{fieldnum}))
                                                        PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(behaviorStates{stateNum}).(outputNames{loopNum}).(dataTypes{dataNum}).(datafields{fieldnum})(:,:,count)=SeparateData.(outputNames{loopNum}).(behaviorStates{stateNum}).(dataTypes{dataNum}).(datafields{fieldnum});
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        else
                            StimType=fieldnames(TheData.(evokeType{1}));
                            for stimNum=1:size(StimType,1)
                                if isstruct(TheData.(evokeType{1}).(StimType{stimNum}))
                                    stimFields=fieldnames(TheData.(evokeType{1}).(StimType{stimNum}));
                                    for fieldNum=1:size(stimFields,1)
                                        if isstruct(TheData.(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}))
                                            roiNames=fieldnames(TheData.(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}));
                                            for roiNum=1:size(roiNames,1)
                                                if ~strcmpi(roiNames{roiNum},'Pixelwise')
                                                    behaviorState=fieldnames(TheData.(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}));
                                                    for behavNum=1:size(behaviorState,1)
                                                        dataField=fieldnames(TheData.(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}).(behaviorState{behavNum}));
                                                        avgFields=strfind(dataField,'_');
                                                        for dataNum=1:size(dataField)
                                                            if ~isempty(avgFields{dataNum})
                                                                if anNum==1
                                                                    count=1;
                                                                else
                                                                    if isfield(PopulationData,ageFolders(ageNum).name)
                                                                        if isfield(PopulationData.(ageFolders(ageNum).name),evokeType{1})
                                                                            if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}),StimType{stimNum})
                                                                                if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(StimType{stimNum}),stimFields{fieldNum})
                                                                                    if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}),roiNames{roiNum})
                                                                                        if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}),behaviorState{behavNum})
                                                                                            if isfield(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}).(behaviorState{behavNum}),dataField{dataNum})
                                                                                                count=size(PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}).(behaviorState{behavNum}).(dataField{dataNum}),3)+1;
                                                                                            else
                                                                                                count=1;
                                                                                            end
                                                                                        else
                                                                                            count=1;
                                                                                        end
                                                                                    else
                                                                                        count=1;
                                                                                    end
                                                                                else
                                                                                    count=1;
                                                                                end
                                                                            else
                                                                                count=1;
                                                                            end
                                                                        else
                                                                            count=1;
                                                                        end
                                                                    else
                                                                        count=1;
                                                                    end
                                                                end
                                                                if ~isnan(TheData.(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}).(behaviorState{behavNum}).(dataField{dataNum}))
                                                                    PopulationData.(ageFolders(ageNum).name).(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}).(behaviorState{behavNum}).(dataField{dataNum})(:,:,count)...
                                                                        =TheData.(evokeType{1}).(StimType{stimNum}).(stimFields{fieldNum}).(roiNames{roiNum}).(behaviorState{behavNum}).(dataField{dataNum});
                                                                end
                                                            end
                                                        end
                                                    end
                                                end
                                            end
                                        end
                                    end
                                else
                                    if anNum==1
                                        PopulationData.(ageFolders(ageNum).name).(evokeType{1}).Params=TheData.(evokeType{1}).(StimType{stimNum});
                                    end
                                end
                            end
                        end
                        clear SeparateData
                    end
                end
            end
        end
    end
end
cd('F:\NeonateSleepPopulationAvg');
save('AgeGroupedStimulusTriggeredData','PopulationData','-v7.3');
end
    
