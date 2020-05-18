%%      function [mold] = GetROI(img, animal, hem)

%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Aug 2013
%   Version 1
%   Revised by Kyle Gheres, Drew Lab, MCIBS,Penn State University, Feb 2017
%   Ver 3

%   SUMMARY: Looks for existing ROI coordinates for a given animal and
%   recorded hemisphere. Loads the ROI if previous ROIs exist, calls
%   function to create a new ROI if none exist.
%_______________________________________________________________________
%   INPUTS:                     img - a sample image used for ROI
%                               reference

%                               animal - animal name as a string

%                               hem - hemisphere recorded as a string
%                                       (LH/RH)
%_______________________________________________________________________
%   OUTPUTS:                    mold - a logical matrix with size of the
%                               sample image.
%_______________________________________________________________________
%   REQUIRED SCRIPTS:           -[mold] = CreateROI(img)
%                               -[Q] = DetectMachine(curr_dir);
%_______________________________________________________________________
%   CALLED BY:                  -Bin2Refl_WhiskTrials.m
%_______________________________________________________________________
%   FUTURE VERSIONS:            -Detect images containing histology tracing
%                               of whiskers and automatically create ROI.
%_______________________________________________________________________
%   CHANGES FROM PREV VERS: Capable of switching between multiple pre-defined
%   ROIs.
%_______________________________________________________________________



function [mold,ROI_name] = GetROI_MultiROI(img, animal, hem,Proc_Vars)
isok= 'n';
%[Q] = DetectMachine_2;
SharVars = GetSharVars_2(animal,hem);
if isfield(SharVars,'ROIs')==0
    SharVars.ROIs = [];
end
if Proc_Vars.fil==1
    ROI_name=input('ROI name? (If pixelwise analysis desired define ROI "Pixelwise", if Opto stim both pixelwise and "Optogenetics" ROI needed) ','s');
else
    ROI_name=Proc_Vars.ROI_name{Proc_Vars.count};
end
if isfield(SharVars.ROIs,ROI_name)
    xi = SharVars.ROIs.(ROI_name).xi;
    yi = SharVars.ROIs.(ROI_name).yi;
    disp('Previous ROI found for this animal and hemisphere.')
    figure(99); imagesc(img); colormap(gray); axis image;
    xlabel('Caudal');
    if strcmp(hem,'LH')
        ylabel('Lateral');
    elseif strcmp(hem,'RH')
        ylabel('Medial')
    end
    impoly(gca,[xi,yi]);
    while strcmp(isok,'n') == 1
        if Proc_Vars.fil>1
            isok='y';
        else
            isok = input('Image and ROI okay? (y/n) ','s');
        end
        if strcmp(isok,'y')
            mold = roipoly(img,xi,yi);
            continue;
        end
        changeimage = input('Enter "1" change ROI, Enter "2" next trial.');
        switch changeimage
            case 1
                [mold] = CreateROI_MultiROI(img,animal,hem,ROI_name);
            case 2
                disp ('Moving to next trial. ')
                normrefl = [];
                err = 'y';
                return;
            otherwise
                disp('Invalid input:');
                ok2 = 'n';
        end
    end
else
    [mold]=CreateROI_MultiROI(img,animal,hem,ROI_name);
end
end
