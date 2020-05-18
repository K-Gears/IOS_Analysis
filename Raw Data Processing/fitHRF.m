function fitHRF(~)




%% Identify and load files
if nargin==0
    filename = uigetfile('*.mat','MultiSelect', 'on');
    prompt='Calculate HRF then average HRFs or concatenate then solve for HRF?[Concatenate/Average]';
    Method=input(prompt,'s');
    prompt='What HRF Model Type do you want to fit with?[Numeric/pixelwise]';
    Model=input(prompt,'s');
end
% Control for single file instead of a list
if iscell(filename) == 0
    filename = {filename};
end
Pad_Time=5;%time in seconds to pad with zeros between concatenated trials
Downsample_Const=5; %Factor to downsample CBV data by before fitting.
Acc_Thresh=1e-5;
order=6;
Wn=1/15; %Bandpass filter properties [low freq, high freq]/nyquist freq
ftype='low';
[zeroa,poleb,gain]=butter(order,Wn,ftype);
[sos,g]=zp2sos(zeroa,poleb,gain);
for fil = 1:length(filename)
    close all;
    indfile = filename{fil};
    load(indfile)
    animal=indfile(1:5);
    date=indfile(10:15);
    hem=indfile(7:8);
    [imp_bin]=velocity_binarize(RawData.vBall,RawData.an_fs,RawData.dal_fr,Acc_Thresh);
    fprintf('Filtering CBV data\n')
    Pixelwise_Refl=filtfilt(sos,g,RawData.HRF.CBVrefl');
    Pixelwise_Refl=Pixelwise_Refl';
    if size(imp_bin,1)>size(Pixelwise_Refl,2)
        imp_bin=imp_bin(1:size(Pixelwise_Refl,2));
    else
        Pixelwise_Refl=Pixelwise_Refl(:,(1:size(imp_bin,1)));
    end
    Norm_Constant=mean(Pixelwise_Refl,2);
    Norm_Refl(size(Pixelwise_Refl,1),size(Pixelwise_Refl,2))=0;
    fprintf('Converting CBV to percent change\n')
    for w=1:size(Pixelwise_Refl,1)
        Norm_Refl(w,:)=((Pixelwise_Refl(w,:)/Norm_Constant(w))-1)*100;
    end
    if strcmpi(Method,'Concatenate')==1
        if fil==1
            Trial_Refl=downsample(Norm_Refl',3);
            Trial_Refl=Trial_Refl';
            Trial_Run=downsample(imp_bin',3);
        else
            fprintf('Concatenating CBV data\n')
            Pad_Refl=padarray(Norm_Refl,[0 (Pad_Time*RawData.dal_fr)],'pre');
            Down_Refl=downsample(Pad_Refl',Downsample_Const);
            Pad_Run=padarray(imp_bin',[0 (Pad_Time*RawData.dal_fr)],'pre');
            Down_Run=downsample(Pad_Run,Downsample_Const);
            Trial_Refl=horzcat(Trial_Refl,Down_Refl');
            Trial_Run=horzcat(Trial_Run,Down_Run);
        end
        
        clear Norm_Refl Pad_Refl Pixelwise_Refl
    else
        if strcmpi(Model,'Numeric')==1
            out = HRFNumerical(imp_bin',Norm_Refl,RawData.dal_fr);
            HRF_Fit.Average.Numeric.HRF{fil}=out;
        else
            Down_Refl=downsample(Norm_Refl',Downsample_Const);
            Down_Run=downsample(imp_bin,Downsample_Const);
            Time=(1:size(Down_Refl,2))/(RawData.dal_fr/Downsample_Const);
            [a, v, c, HRFa, HRFv, HRF, xFit_HRFa, xFit_HRFv, xFit_HRF, cc] = HRFModel(Down_Run,Down_Refl,Time);
            Red=getRGBcomponent(a);
            Blue=getRGBcomponent(v);
            HRF_Fit.Average.Pixelwise.ImgColor((1:256),(1:256),3)=0;
            for m=1:size(a,1)
                HRF_Fit.Average.Pixelwise.ImgColor(RawData.HRF.Pixel_Map(1,m),RawData.HRF.Pixel_Map(2,m),1)=Red(m);
                HRF_Fit.Average.Pixelwise.ImgColor(RawData.HRF.Pixel_Map(1,m),RawData.HRF.Pixel_Map(2,m),3)=Blue(m);
            end
            image(HRF_Fit.Average.Pixelwise.ImgColor);
            ax=gca;
            ax.XLabel.String='Lateral';
            ax.YLabel.String='Caudal';
            ax.Title.String='HRF Arterial/Venous pixel weights';
            
        end
    end
    
end
fprintf('Fitting HRF data\n')
if strcmpi(Method,'Concatenate')==1
    if strcmpi(Model,'Numeric')==1
        out = HRFNumerical(Trial_Run',Trial_Refl,(RawData.dal_fr/Downsample_Const));
        HRF_NumFit.Concatenate.HRF=out;
        clear out
        HRF_NumFit.Movie((1:256),(1:256),size(HRF_NumFit.Concatenate.HRF,2))=0;
        for u=1:size(HRF_NumFit.Concatenate.HRF,2)
            for w=1:size(RawData.HRF.Pixel_Map,2)
                HRF_Fit.Concatenate.Numeric.Movie(RawData.HRF.Pixel_Map(1,w),RawData.HRF.Pixel_Map(2,w),u)=HRF_NumFit.Concatenate.HRF(w,u);
            end
        end
    else
        Time=(1:size(Trial_Run,2))/(RawData.dal_fr/Downsample_Const);
        [a, v, c, HRFa, HRFv, HRF, xFit_HRFa, xFit_HRFv, xFit_HRF, cc] = HRFModel(Trial_Run,Trial_Refl,Time);
          Red=getRGBcomponent(a);
            Blue=getRGBcomponent(v);
            HRF_Fit.Concatenate.Pixelwise.ImgColor((1:256),(1:256),3)=0;
            for m=1:size(a,1)
                HRF_Fit.Concatenate.Pixelwise.ImgColor(RawData.HRF.Pixel_Map(1,m),RawData.HRF.Pixel_Map(2,m),1)=Red(m);
                HRF_Fit.Concatenate.Pixelwise.ImgColor(RawData.HRF.Pixel_Map(1,m),RawData.HRF.Pixel_Map(2,m),3)=Blue(m);
            end
            image(HRF_Fit.Concatenate.Pixelwise.ImgColor);
            ax=gca;
            ax.XLabel.String='Lateral';
            ax.YLabel.String='Caudal';
            ax.Title.String='HRF Arterial/Venous pixel weights';
            savefig([animal '_' date '_HRF Pixel Weights']);
    end
else
    for q=1:size(HRF_NumFit.Average.HRF,2)
        The_Data(:,:,q)=HRF_NumFit.Average.HRF{q}.HRF;
    end
    HRF_NumFit.Average.AVG_HRF=mean(The_Data,3);
    HRF_NumFit.Movie((1:256),(1:256),size(HRF_NumFit.Average.AVG_HRF,2))=0;
    for u=1:size(HRF_NumFit.Average.AVG_HRF,2)
        for w=1:size(RawData.HRF.Pixel_Map,2)
            HRF_NumFit.Movie(RawData.HRF.Pixel_Map(1,w),RawData.HRF.Pixel_Map(2,w),u)=HRF_NumFit.Average.AVG_HRF(w,u);
        end
    end
    
end
fprintf('Saving HRF data\n')
save([animal '_' date '_' hem '_HRF_Fit'],'HRF_Fit','-v7.3');
end