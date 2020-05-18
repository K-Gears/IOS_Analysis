function AvgRestingIOSPowerSpectrum
cd('D:\Resting_IOS_PowerSpectrum');
Thefolders=dir;
for foldernum=3:length(Thefolders)
    if Thefolders(foldernum).isdir==1
    cd([Thefolders(foldernum).folder '\' Thefolders(foldernum).name]);
    thefiles=dir('*RestingPowerSpectrum*.mat');
    for filNum=1:length(thefiles)
        load(thefiles(filNum).name);
        
        PopulationPowerSpectrums.(Thefolders(foldernum).name).Animalname{filNum}=thefiles(filNum).name(1:9);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerSpectrum(:,filNum)=S;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrum(:,filNum)=wS;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).NormWhitenedPowerSpectrum(:,filNum)=wS./sum(wS);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).NormPowerSpectrum(:,filNum)=S./sum(S);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).trialnum(:,filNum)=trialnum;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).ErrorBars(:,:,filNum)=Serr;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).Frequency(:,filNum)=f;
        HighCut=find(f<=1,1,'last');
        out=PowerLawFit(f(1:HighCut),S(1:HighCut));%(1:HighCut)
        PopulationPowerSpectrums.(Thefolders(foldernum).name).NormWhitenedPowerSpectrumLow(:,filNum)=wS(1:HighCut)./sum(wS(1:HighCut));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrumLow(:,filNum)=wS(1:HighCut);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerPeak(:,filNum)=max(wS(1:HighCut));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerLoc(:,filNum)=f(find(wS(1:HighCut)==max(wS(1:HighCut))));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).NormPowerSpectrumLow(:,filNum)=S(1:HighCut)./sum(S(1:HighCut));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).NormPowerPeak(:,filNum)=PopulationPowerSpectrums.(Thefolders(foldernum).name).NormWhitenedPowerSpectrumLow(wS(1:HighCut)==max(wS(1:HighCut)),filNum);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).FrequencyLow(:,filNum)=f(1:HighCut);        
        PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.exponent(:,filNum)=out.LSERout.exponent;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.intercept(:,filNum)=out.LSERout.intercept;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.Rsqr(:,filNum)=1-(sum((out.LSERout.residual.^2))./sum(((log10(out.LSERout.power_resamp)-mean(log10(out.LSERout.power_resamp))).^2)));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.FitData(filNum)=out;
        
        PopulationPowerSpectrums.(Thefolders(foldernum).name).LowsumBounds=[0.005 0.05];
        StartInd=find(PopulationPowerSpectrums.(Thefolders(foldernum).name).FrequencyLow(:,filNum)<PopulationPowerSpectrums.(Thefolders(foldernum).name).LowsumBounds(1),1,'last')+1;
        EndInd=find(PopulationPowerSpectrums.(Thefolders(foldernum).name).FrequencyLow(:,filNum)>PopulationPowerSpectrums.(Thefolders(foldernum).name).LowsumBounds(2),1,'first')-1;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).Lowsum(filNum)=sum(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrumLow((StartInd:EndInd),filNum))...
            /length(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrumLow((StartInd:EndInd),filNum));
        
        PopulationPowerSpectrums.(Thefolders(foldernum).name).HighsumBounds=[0.1 1];
        StartInd=find(PopulationPowerSpectrums.(Thefolders(foldernum).name).FrequencyLow(:,filNum)<PopulationPowerSpectrums.(Thefolders(foldernum).name).HighsumBounds(1),1,'last')+1;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).Highsum(filNum)=sum(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrumLow((StartInd:end),filNum))...
            /length(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrumLow((StartInd:end),filNum));
        
        PopulationPowerSpectrums.(Thefolders(foldernum).name).HighLowRatio(filNum)=PopulationPowerSpectrums.(Thefolders(foldernum).name).Highsum(filNum)/PopulationPowerSpectrums.(Thefolders(foldernum).name).Lowsum(filNum);
    end
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PowerSpectrum=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerSpectrum,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormPowerSpectrum=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).NormPowerSpectrum,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.WhitenedPowerSpectrum=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrum,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.WhitenedPowerPeak=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerPeak,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.WhitenedPowerLoc=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerLoc,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormWhitenedPowerSpectrum=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).NormWhitenedPowerSpectrum,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormWhitenedPowerSpectrumLow=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).NormWhitenedPowerSpectrumLow,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormPowerPeak=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).NormPowerPeak,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.WhitenedPowerSpectrumLow=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).WhitenedPowerSpectrumLow,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormPowerSpectrumLow=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).NormPowerSpectrumLow,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PowerLawExponent=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.exponent,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PowerLawIntercept=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.intercept,2);
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PowerLawRsqr=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).PowerLaw.Rsqr,2);
        out=PowerLawFit(f(1:HighCut),PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PowerSpectrum(1:HighCut));%(1:HighCut)
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PopPowerLawExp=out.LSERout.exponent;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PopPowerLawInt=out.LSERout.intercept;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PopPowerLawRsqr=1-(sum((out.LSERout.residual.^2))./sum(((log10(out.LSERout.power_resamp)-mean(log10(out.LSERout.power_resamp))).^2)));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PopWhitenedPowerLoc=f(find(PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.WhitenedPowerSpectrumLow==max(PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.WhitenedPowerSpectrumLow)));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.Fit=out;
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.PopNormPowerPeak=PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormWhitenedPowerSpectrumLow(PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormWhitenedPowerSpectrumLow==max(PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.NormWhitenedPowerSpectrumLow));
        PopulationPowerSpectrums.(Thefolders(foldernum).name).AveragedData.HighLowRatio=mean(PopulationPowerSpectrums.(Thefolders(foldernum).name).HighLowRatio);
    end
end
cd('D:\Resting_IOS_PowerSpectrum');
save('Population_Resting_IOS_PowerSpectrum','PopulationPowerSpectrums','-v7.3');
%% Make Plot Color Maps
plotmaps(1:6,1:3)=NaN;
[cmap]=brewermap(11,'RdYlBu');
plotmaps(1,:)=cmap(11,:);
plotmaps(2,:)=cmap(10,:);
plotmaps(4,:)=cmap(3,:);
[cmap]=brewermap(6,'Dark2');
plotmaps(3,:)=cmap(6,:);
[cmap]=brewermap(1,'Set1');
plotmaps(5,:)=cmap(1,:);
[cmap]=brewermap(9,'YlOrRd');
plotmaps(6,:)=cmap(9,:);

%% Plot Data

Age=fieldnames(PopulationPowerSpectrums);
figure(87);hold on
for AgeNum=1:length(Age)
    loglog(PopulationPowerSpectrums.(Age{AgeNum}).Frequency(:,1),PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.PowerSpectrum);
    lgndTxt{AgeNum}=[Age{AgeNum} ' n=' num2str(length(PopulationPowerSpectrums.(Age{AgeNum}).Animalname))];
end
title('Population IOS resting state power spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (a.u)');
legend(lgndTxt);
set(gca,'YScale','Log');
set(gca,'XScale','Log');

figure(88);hold on
for AgeNum=1:length(Age)
    loglog(PopulationPowerSpectrums.(Age{AgeNum}).Frequency(:,1),PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.NormPowerSpectrum);
    lgndTxt{AgeNum}=[Age{AgeNum} ' n=' num2str(length(PopulationPowerSpectrums.(Age{AgeNum}).Animalname))];
end
title('Population IOS resting state normalized power spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (a.u)');
legend(lgndTxt);
set(gca,'YScale','Log');
set(gca,'XScale','Log');

figure(89);hold on
for AgeNum=1:length(Age)
    loglog(PopulationPowerSpectrums.(Age{AgeNum}).Frequency(:,1),PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.WhitenedPowerSpectrum);
    lgndTxt{AgeNum}=[Age{AgeNum} ' n=' num2str(length(PopulationPowerSpectrums.(Age{AgeNum}).Animalname))];
end
title('Population IOS resting state prewhitened power spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (a.u)');
legend(lgndTxt);
set(gca,'YScale','Log');
set(gca,'XScale','Log');
xlim([0 15]);

figure(90);hold on
for AgeNum=1:length(Age)
    loglog(PopulationPowerSpectrums.(Age{AgeNum}).Frequency(:,1),PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.NormWhitenedPowerSpectrum);
    lgndTxt{AgeNum}=[Age{AgeNum} ' n=' num2str(length(PopulationPowerSpectrums.(Age{AgeNum}).Animalname))];
end
title('Population IOS resting state prewhitened normalized power spectrum');
xlabel('Frequency (Hz)');
ylabel('Power (a.u)');
legend(lgndTxt);
set(gca,'YScale','Log');
set(gca,'XScale','Log');
xlim([0 15]);

figure(91);hold on
for AgeNum=1:length(Age)
    loglog(PopulationPowerSpectrums.(Age{AgeNum}).FrequencyLow(:,1),PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.NormWhitenedPowerSpectrumLow);
    lgndTxt{AgeNum}=[Age{AgeNum} ' n=' num2str(length(PopulationPowerSpectrums.(Age{AgeNum}).Animalname))];
end
title('Population IOS resting state prewhitened normalized power spectrum [DC:1]');
xlabel('Frequency (Hz)');
ylabel('Power (a.u)');
legend(lgndTxt);
set(gca,'YScale','Log');
set(gca,'XScale','Log');
xlim([0 1]);

figure(92);hold on
for AgeNum=1:length(Age)
    semilogx(PopulationPowerSpectrums.(Age{AgeNum}).FrequencyLow(:,1),10*log10(PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.WhitenedPowerSpectrumLow));
    lgndTxt{AgeNum}=[Age{AgeNum} ' n=' num2str(length(PopulationPowerSpectrums.(Age{AgeNum}).Animalname))];
end
title('Population IOS resting state prewhitened power spectrum [DC:1]');
xlabel('Frequency (Hz)');
ylabel('Power (dB)');
legend(lgndTxt);
%set(gca,'YScale','Log');
set(gca,'XScale','Log');
xlim([0 1]);

figure(93);hold on
for AgeNum=1:length(Age)
    xVals(1:length(PopulationPowerSpectrums.(Age{AgeNum}).PowerLaw.exponent))=AgeNum;
    scatter(xVals,PopulationPowerSpectrums.(Age{AgeNum}).PowerLaw.exponent,108,'filled','MarkerFaceColor',plotmaps(AgeNum,:),'MarkerEdgeColor','k');
    lgndTxt{AgeNum}=Age{AgeNum};
    clear xVals
    AvgExponent(AgeNum)=PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.PopPowerLawExp;
end
plot(AvgExponent,'Color','k','Marker','s','MarkerSize',10,'MarkerFaceColor','k');
title('Population IOS slope power law fit');
xlabel('Animal age');
ylabel('power law slope');
xticks(1:1:length(Age));
xticklabels(lgndTxt);

figure(94);hold on
for AgeNum=1:length(Age)
    xVals(1:length(PopulationPowerSpectrums.(Age{AgeNum}).WhitenedPowerLoc))=AgeNum;
    scatter(xVals,PopulationPowerSpectrums.(Age{AgeNum}).WhitenedPowerLoc,108,'filled','MarkerFaceColor',plotmaps(AgeNum,:),'MarkerEdgeColor','k');
    lgndTxt{AgeNum}=Age{AgeNum};
    clear xVals
    AvgPeak(AgeNum)=PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.PopWhitenedPowerLoc;
end
plot(AvgPeak,'Color','k','Marker','s','MarkerSize',10,'MarkerFaceColor','k');
set(gca,'YScale','Log');
title('Population IOS power spectrum peak');
xlabel('Animal age');
ylabel('Fequency (Hz)');
xticks(1:1:length(Age));
xticklabels(lgndTxt);
xlim([0 length(Age)+1]);
ylim([0 1]);

figure(95);hold on
for AgeNum=1:length(Age)
    xVals(1:length(PopulationPowerSpectrums.(Age{AgeNum}).NormPowerPeak))=AgeNum;
    scatter(xVals,(PopulationPowerSpectrums.(Age{AgeNum}).NormPowerPeak*100),108,'filled','MarkerFaceColor',plotmaps(AgeNum,:),'MarkerEdgeColor','k');
    lgndTxt{AgeNum}=Age{AgeNum};
    clear xVals
    AvgPeak(AgeNum)=PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.PopNormPowerPeak*100;
end
plot(AvgPeak,'Color','k','Marker','s','MarkerSize',10,'MarkerFaceColor','k');
set(gca,'YScale','Log');
title('Population IOS power spectrum peak Power');
xlabel('Animal age');
ylabel('Percentage total spectral power (%)');
xticks(1:1:length(Age));
xticklabels(lgndTxt);
xlim([0 length(Age)+1]);
ylim([0 1]);

figure(96);hold on
for AgeNum=1:length(Age)
    xVals(1:length(PopulationPowerSpectrums.(Age{AgeNum}).NormPowerPeak))=AgeNum;
    scatter(xVals,(PopulationPowerSpectrums.(Age{AgeNum}).HighLowRatio),108,'filled','MarkerFaceColor',plotmaps(AgeNum,:),'MarkerEdgeColor','k');
    lgndTxt{AgeNum}=Age{AgeNum};
    clear xVals
    AvgRatio(AgeNum)=PopulationPowerSpectrums.(Age{AgeNum}).AveragedData.HighLowRatio;
end
plot(AvgRatio,'Color','k','Marker','s','MarkerSize',10,'MarkerFaceColor','k');
title('Population IOS power spectrum peak Power');
xlabel('Animal age');
ylabel('Percentage total spectral power (%)');
xticks(1:1:length(Age));
xticklabels(lgndTxt);
xlim([0 length(Age)+1]);
ylim([0 1]);
