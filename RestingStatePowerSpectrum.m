function RestingStatePowerSpectrum(~)
Drivenames={'E:\','F:\','G:\'};
for driveNum=1:size(Drivenames,2)
    cd(Drivenames{driveNum});
    folders=dir;
    for foldernum=1:length(folders)
        Test(foldernum)=strcmpi(folders(foldernum).name,'NeonateSleepIndividualAnimals');
    end
    Animalfolder=find(Test==1);
    cd([folders(Animalfolder).folder '\' folders(Animalfolder).name]);
    dates=dir;
    for dateNum=3:length(dates)
        if dates(dateNum).isdir==1
        cd([dates(dateNum).folder '\' dates(dateNum).name]);
        Animals=dir;
        for animalNum=3:length(Animals)
            cd([Animals(animalNum).folder '\' Animals(animalNum).name]);
            thefiles=dir('*rawdata.mat');
            RestData=[];
            ledType='M565L3';
            bandfilterType='FB570-10';
            cutfilterType='FEL0500';
            [~,~,weightedcoeffHbT]=getHbcoeffs(ledType,bandfilterType,cutfilterType);
            params.Fs=30;
            params.tapers=[3 5];
            params.fpass=[0,15];
            params.trialave=1;
            params.err=[2 0.05];
            for filnum=1:size(thefiles,1)
                load(thefiles(filnum).name);
                [imp_bin]=velocity_binarize(RawData.vBall,RawData.an_fs,RawData.dal_fr,1e-4);
                LEDtrig=round(RawData.LED,0);
                Soltrig=round(RawData.Sol,0);
                runPts=find(imp_bin==1);
                runPerc=(length(runPts)/length(imp_bin))*100;
                if isfield(RawData,'IOS')
                if runPerc<10
                    if max(LEDtrig)==0
                        if max(Soltrig)==0
                            trialnum=size(RestData,2)+1;
                            RestData(:,trialnum)=detrend((weightedcoeffHbT*log(RawData.IOS.barrels.CBVrefl/mean(RawData.IOS.barrels.CBVrefl)))*1e6,'linear');
                            WhitenedData(:,trialnum)=diff(RestData(:,trialnum));
                        end
                    end
                end
                end
            end
            if ~isempty(RestData)
            [S,f,Serr]=mtspectrumc(RestData,params);
            [wS,wF,wSerr]=mtspectrumc(WhitenedData,params);
            figure;loglog(f,S);
            xlim([0 15]);
            procHold=date;
            procDate=strrep(procHold,'-','_');
            cd('D:\TempRestingData');
            save([Animals(animalNum).name '_RestingPowerSpectrum_' dates(dateNum).name '_' procDate],'S','f','Serr','trialnum','wS','wF','wSerr');
            else
                fprintf(['No resting data ' Animals(animalNum).name '\n'])
            end
        end
        end
    end
end
end

