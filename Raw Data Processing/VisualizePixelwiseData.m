function VisualizePixelwiseData(PixelwiseData,plotParams)
frameTimes=round(((Lead_Time-1):0.5:(Lead_Time+5.5))*ProcessedData.dal_fr,0);
labelTimes=(frameTimes/ProcessedData.dal_fr)-Lead_Time;
EmptyFrame((1:256),(1:256))=NaN;
PixMap=ProcessedData.pixelMaps.Pixelwise;
colnum=3;
rownum=ceil(numel(frameTimes)/colnum);
for n=1:size(Stim_Type,1)
    roiFrame=EmptyFrame;
    if strcmpi(Stim_Type{n},'Laser_Stim')
        roiFrame(round(SharVars.ROIs.Optogenetics.ycenter,0),round(SharVars.ROIs.Optogenetics.xcenter,0))=1;
        rotROI=rot90(roiFrame,3);
        [row,col]=find(rotROI==1);
        TheROI(:,1)=row;
        TheROI(:,2)=col;
    else
        for pic=1:length(SharVars.ROIs.barrels.xi)
            roiFrame(round(SharVars.ROIs.barrels.yi(pic),0),round(SharVars.ROIs.barrels.xi(pic),0))=1;
        end
        rotROI=rot90(roiFrame,3);
        [row,col]=find(rotROI==1);
        [TheROI(:,1),Inds]=sort(col);
         TheROI(:,2)=row(Inds);
         midRow=mean(TheROI(:,2));
         count=1;
         bump=1;
         for vert=1:length(newCol)
             if newRow(vert)<midRow
                 UpperInds(count,1)=newCol(vert);
                 UpperInds(count,2)=newRow(vert);
                 count=count+1;
             else
                 LowerInds(bump,1)=newCol(vert);
                 LowerInds(bump,2)=newRow(vert);
                 bump=bump+1;
             end
         end
         [newLower(:,1),Inds]=sort(LowerInds(:,1),'descend');
         newLower(:,2)=LowerInds(Inds,2);
         visROI=[UpperInds;newLower];
         visROI(length(visROI)+1,:)=visROI(1,:);
    end
    theStim=strrep(Stim_Type{n},'_',' ');
    for m=1:size(Run_State,2)
        if ~isempty(ProcessedData.(Stim_Type{n}).IOS.Pixelwise.(Run_State{m}).Refl)
        for framenum=1:length(frameTimes)
            theFrame(:,framenum)=mean(ProcessedData.(Stim_Type{n}).IOS.Pixelwise.(Run_State{m}).Avg_Refl(:,(frameTimes(framenum)-2):(frameTimes(framenum)+2)),2);
            MakePic(:,:,framenum)=EmptyFrame;
            for pixnum=1:size(theFrame,1)
                MakePic(PixMap(1,pixnum),PixMap(2,pixnum),framenum)=theFrame(pixnum,framenum);
            end
            plotLabel=['t= ' num2str(labelTimes(framenum))];
            if framenum==1
                theImg=EmptyFrame;
                for frame=1:5
                    for pixnum=1:size(RawData.IOS.Pixelwise.PixelMap,2)
                        theImg(PixMap(1,pixnum),PixMap(2,pixnum),frame)=RawData.IOS.Pixelwise.CBVrefl(pixnum,frame);
                    end
                    RefImg=mean(theImg,3);
                    [row,col]=find(~isnan(rot90(RefImg,3)));
                    row=sort(row); %find bounds for image display
                    col=sort(col); %find bounds for image display
   
                end
                figure(99);subplot(rownum,colnum,framenum);imagesc(rot90(RefImg,3));
                colormap gray;caxis([0 4095]); axis image; axis off; 
                xlim([col(1) col(end)]);ylim([row(1),row(end)]);
            end 
            figure(99);subplot(rownum,colnum,framenum+1);imagesc(rot90(MakePic(:,:,framenum),3));
            colormap parula; axis image; axis off; 
            xlim([col(1) col(end)]);ylim([row(1),row(end)]);
            title(plotLabel);
            if strcmpi(Stim_Type{n},'Laser_Stim')
               caxis([-2 2]);
               drawcircle('Center',[TheROI(2),TheROI(1)],'Radius',15,'Color',[1 1 1],'LineWidth',0.25,'FaceAlpha',0,'InteractionsAllowed','none');
            else
               caxis([-5 5]);
               drawpolygon('Position',visROI,'Color',[1 1 1],'FaceAlpha',0,'LineWidth',0.25,'InteractionsAllowed','none');
            end
            
        end
        runTitle={'still','running'};
        theTitle=['Average response to ' theStim ' while animal was ' runTitle{m}];
        sgtitle(theTitle);
        savefig([animal '_' date '_' Stim_Type{n} '_' Run_State{m} '_Pixelwise Average Stimulus Evoked Reflectance Change']);
        close;
        end
    end
    clear TheROI
end
close all;