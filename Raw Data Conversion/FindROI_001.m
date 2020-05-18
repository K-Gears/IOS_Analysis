function FindROI_001(filename,FPS,Limits)
[theimage]=ReadDalsaBinary_Matrix(filename,256,256);
Frequency_Cutoff=1/(30*0.5);
[z,p,k]=butter(6,Frequency_Cutoff);
[sos,g]=zp2sos(z,p,k);
  NewImage=double(theimage); %convert 12bit matrix to double to prevent round off errors
  NormImage=zeros(256,256,6000);
    clear theimage;
    AvgImage=mean(NewImage,3);
    for K=1:6000;
        NormImage(:,:,K)=((NewImage(:,:,K)-AvgImage)./AvgImage); %Subtract the mean from each frame prior to normalizing by mean
    end
    Filt_Image=permute(NormImage,[3,2,1]);
    clear NormImage;
    Filtered_Image=filtfilt(sos,g,Filt_Image);
    clear Filt_Image;
    NormImage=permute(Filtered_Image,[3,2,1]);
    NormImage=NormImage*100;
    figure(55);imagesc(NewImage(:,:,150));
    ylabel('Rostral','FontSize',10,'FontWeight','bold');
    xlabel('Lateral','FontSize',10,'FontWeight','bold');
    Frames_Per_Second=FPS;
    intensity_limits=Limits;
    handle=implay(NormImage,Frames_Per_Second);
    handle.Visual.ColorMap.UserRangeMin=intensity_limits(1);
    handle.Visual.ColorMap.UserRangeMax=intensity_limits(2);
    handle.Visual.ColorMap.UserRange=1;  
    ylabel('Rostral','FontSize',10,'FontWeight','bold');
    xlabel('Lateral','FontSize',10,'FontWeight','bold');
end