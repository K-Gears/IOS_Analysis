function BulkAnalyze_OptogeneticIOS(~)
Drivenames={'E:\','F:\','G:\'}; %
for driveNum=1:size(Drivenames,2)
    cd(Drivenames{driveNum});
    folders=dir;
    folders(~[folders.isdir])=[];
    tf=ismember({folders.name},{'.','..'});
    folders(tf)=[];
    tf2=ismember({folders.name},'NeonateSleepIndividualAnimals');
    folders(tf2==0)=[];
    %     for foldernum=1:length(folders)
    %         Test(foldernum)=strcmpi(folders(foldernum).name,'NeonateSleepIndividualAnimals');
    %     end
    %     Animalfolder=find(Test==1);
    cd([folders.folder '\' folders.name]);
    dates=dir;
    dates(~[dates.isdir])=[];
    tf3=ismember({dates.name},{'.','..'});
    dates(tf3)=[];
    %     if strcmpi(Drivenames{driveNum},'F:\')
    %         startfolder=6;
    %     else
    %         startfolder=3;
    %     end
    for dateNum=1:length(dates)
        if dates(dateNum).isdir==1
            if ~strcmpi(dates(dateNum).name,'070919')
                cd([dates(dateNum).folder '\' dates(dateNum).name]);
                Age=dir;
                Age(~[Age.isdir])=[];
                tf4=ismember({Age.name},{'.','..'});
                Age(tf4)=[];
                for ageNum=1:size(Age,1)
                    cd([Age(ageNum).folder '\' Age(ageNum).name]);
                    Animals=dir;
                    Animals(~[Animals.isdir])=[];
                    tf5=ismember({Animals.name},{'.','..'});
                    Animals(tf5)=[];
                    for animalNum=1:length(Animals)
                        cd([Animals(animalNum).folder '\' Animals(animalNum).name]);
                        thefiles=dir('*rawdata.mat');
                        for filnum=1:size(thefiles,1)
                            filenames{filnum}=thefiles(filnum).name;
                        end
                        animalAge=Age(ageNum).name;
                        ledType='M565L3';
                        bandfilterType='FB570-10';
                        cutfilterType='FEL0500';
                        Chunk_IOS_CBV_Optogenetics_Neuro_004(filenames,ledType,bandfilterType,cutfilterType,animalAge);
                        clear filenames
                    end
                end
            end
        end
    end
end
end