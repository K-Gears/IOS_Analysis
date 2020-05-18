plotTime=((1:length(lowRefl(ind,:)))/FrameRate)-leadtime;
peakScatter(1:length(lowRefl(ind,:)))=NaN;
peakScatter(locPosPeaks)=maxPosPeaks;
widthScatter(1:length(lowRefl(ind,:)))=NaN;
widthScatter([leftEdge,rightEdge])=maxPosPeaks*0.5;
figure(23);hold on;
plot(plotTime,lowRefl);
scatter(plotTime,peakScatter,'filled','r');
scatter(plotTime,widthScatter,'filled','k');
