function [CBVrefl,flashdata,Proc_Vars] = Bin2Intensity_MultiROI(filename,animal,hem,flash,match,Proc_Vars)
%
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Mar 2013
%   Version 1
%   Revised by Kyle Gheres, Drew Lab, MCIBS, Penn State University, Feb
%   2017 ver 3
%
%   SUMMARY: Converts .bin file into a timeseries of average pixel
%   intensity values of a user defined ROI. Also converts ROI in to
%   pixelswise time series for HRF calculation
%_______________________________________________________________
%   INPUTS:                     filename - name of file with extension
%
%                               animal - animal name as a string
%
%                               hem - Recorded Hemisphere (LH/RH)
%
%                               flash - string (y/n) indicating whether
%                               LEDs were flashed during the trial
%
%                               machineVars - machine specific values
%                               which will store saved information about
%                               ROIs for the animal.
%
%                               match - 1 or 0 to indicate whether the
%                               current file is from the same session as
%                               the previous file. If so, the ROI does not
%                               need to be verified since the animal has
%                               not moved.
%_______________________________________________________________
%   OUTPUTS:                    CBVrefl - A timeseries of CBV reflectance
%                               values calculated from the user defined
%                               ROI.
%
%                               FlashData - A timeseries of average pixel
%                               intensity values from frames where LEDs
%                               were flashed.
%                               
%                               HRF- A structure containing [x,y] locations
%                               of pixels in CBVrefl.HRF defined by HRF ROI
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%       - Declares dalsa frequency based on the status of "FlashTrial"
%       field from the .tdms file.
%
%      - Capable of defining multiple ROIs with HRF ROI stored pixelwise to
%      allow for HRF calculation with behavior across entire window [pixel
%      x time]
%
%_____________________________________________________________

% Get Machine Specific Variables


% Convert Binary Files based on LED configuration for the trial
importstatus = 0;
while importstatus == 0
    if strcmp(flash,'y')
        [CBVframes,flashframes]=ReadDalsaBinary_Flash(filename, 256, 256);
        importstatus = 1;
    elseif strcmp(flash,'n')
        [CBVframes]=ReadDalsaBinary_WhiskTrials(filename,256,256);
        flashframes = {};
        importstatus = 1;
    end
end

nCBVframes = length(CBVframes);
nFlashframes = length(flashframes);
% Find 1st visible frame
framenum = 2;
display(['Animal: ' animal]);
display(['Hemisphere: ' hem]);
display(['Filename: ' filename]);
currimg = CBVframes{framenum};
% Define ROIs
if Proc_Vars.fil==1
Proc_Vars.numROI=input('How many ROIs do you wish to define?');
end
for d=1:Proc_Vars.numROI
    Proc_Vars.count=d;
    if match == 0
        if Proc_Vars.fil==1
            [mold,ROI_name] = GetROI_MultiROI(currimg,animal,hem,Proc_Vars);
            Proc_Vars.ROI_name{d}=ROI_name;
            if ~strcmpi(ROI_name,'Optogenetics')
                Proc_Vars.mold{d}=mold;
            else
                PixROI=strcmpi(Proc_Vars.ROI_name,'Pixelwise')==1;
                if max(PixROI)==1
                    Proc_Vars.mold{d}=Proc_Vars.mold{PixROI}.*mold;
                else
                    fprintf('Define ROI of whole window "Pixelwise"\n')
                    [pixMold,~] = GetROI_MultiROI(currimg,animal,hem,Proc_Vars);
                    Proc_Vars.ROI_name{d+1}='Pixelwise';
                    Proc_Vars.mold{d+1}=pixMold;
                    Proc_Vars.mold{d}=Proc_Vars.mold{d+1}.*mold;
                end
            end
        else
            ROI_name=Proc_Vars.ROI_name{d};
            mold=Proc_Vars.mold{d};
        end
    elseif match == 1
        SharVars = GetSharVars_2(animal,hem);
        disp('Previous ROI found for this animal and hemisphere, loading data...');
        ROI_name=Proc_Vars.ROI_name{d};
        xi = SharVars.ROIs.(ROI_name).xi;
        yi = SharVars.ROIs.(ROI_name).yi;
        mold = roipoly(currimg,xi,yi);
    end
    if strcmpi(ROI_name,'Pixelwise')==1
        [row,col]=find(mold==1);
        CBVrefl.(ROI_name).Pixel_Map(1,:)=row;
        CBVrefl.(ROI_name).Pixel_Map(2,:)=col;
        CBVrefl.(ROI_name).Refl((1:size(CBVrefl.(ROI_name).Pixel_Map,2)),(nCBVframes-framenum))=0;
        newMold=double(mold);
        newMold(newMold==0)=NaN; 
        for m=1:nCBVframes
            ImgMatrix=newMold.*double(CBVframes{m});
            ImgVector=reshape(ImgMatrix,[numel(ImgMatrix),1]);
            ImgVector(isnan(ImgVector))=[];
            CBVrefl.(ROI_name).Refl(:,m)=ImgVector;%pixels should be arranged based off of increasing column number
%         for m=1:nCBVframes
%             ThePixels=mold.*double(CBVframes{m});
%             for p=1:size(CBVrefl.(ROI_name).Pixel_Map,2)
%                 CBVrefl.(ROI_name).Refl(p,m)=ThePixels(CBVrefl.(ROI_name).Pixel_Map(1,p),CBVrefl.(ROI_name).Pixel_Map(2,p));
%             end
%         end
        end
    else
        [row,col]=find(mold==1);
        CBVrefl.(ROI_name).Pixel_Map(1,:)=row;
        CBVrefl.(ROI_name).Pixel_Map(2,:)=col;
        % Convert CBV frames into timeseries
        CBVrefl.(ROI_name).Refl = zeros(1,nCBVframes-framenum);
        for n = 1:nCBVframes
            mask = mold.*double(CBVframes{n});
            CBVrefl.(ROI_name).Refl(n) = mean(nonzeros(mask));
        end
    end
    
    % Convert Flash frames into timeseries, if exist
    if nFlashframes > 0
        flashdata = zeros(1,nCBVframes-framenum);
        for n = 1:nFlashframes
            mask = mold.*double(flashframes{n});
            flashdata(n) = mean(nonzeros(mask));
        end
    else
        flashdata = [];
    end
end
end




    