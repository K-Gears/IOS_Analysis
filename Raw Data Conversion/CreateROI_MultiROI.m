%%      function [mold] = CreateROI(img)

%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Oct 2013
%   Version 1
%   Revised by Kyle Gheres, Drew Lab, MCIBS, Penn State University, Feb
%   2017. Ver. 3
%   SUMMARY: Prompts user to create an region of interest given an image.
%__________________________________________________________________________
%   INPUTS:                     img - image used as reference for ROI
%                               creation.
%__________________________________________________________________________
%   OUTPUTS:                    mold - a logical matrix with size of the 
%                               sample image.
%__________________________________________________________________________
%   REQUIRED SCRIPTS:           None
%__________________________________________________________________________
%   CALLED BY:                  [mold] = GetROI(img, animal, hem)
%__________________________________________________________________________
%   FUTURE VERSIONS:            - Create ROIs that can omit bubbles and junk
%__________________________________________________________________________
%   CHANGES FROM PREV VERS: Prompt to create multiple ROI to process multiple
%   regions simultaneously for CBV changes 2/2017 KG
%__________________________________________________________________________

function [mold] = CreateROI_MultiROI(img,animal,hem,ROI_name)
Q = DetectMachine_2;
[SharVars] = GetSharVars_2(animal,hem);
disp('Please select your region of interest.')
figure(99); imagesc(img); colormap(gray); axis image;
xlabel('Caudal');
if strcmp(hem,'LH')
    ylabel('Lateral');
elseif strcmp(hem,'RH')
    ylabel('Medial')
end
if ~strcmpi(ROI_name,'Optogenetics')
    [mold, xi, yi] = roipoly;
    SharVars.ROIs.(ROI_name).xi = xi;
    SharVars.ROIs.(ROI_name).yi = yi;
    impoly(gca,[xi,yi]);
    asksave = input('Save this ROI for future? (y/n) ','s');
    if asksave == 'y'
        save([Q.machineVars.SharedPath animal '_' hem '_SharVars.mat'], 'SharVars');
    end
else
    [mold,xi,yi,xcenter,ycenter]=DefineOptoROI(img);
    SharVars.ROIs.(ROI_name).xi = xi;
    SharVars.ROIs.(ROI_name).yi = yi;
    SharVars.ROIs.(ROI_name).xcenter = xcenter;
    SharVars.ROIs.(ROI_name).ycenter = ycenter;
    impoly(gca,[xi,yi]);
    asksave = input('Save this ROI for future? (y/n) ','s');
    if asksave == 'y'
        save([Q.machineVars.SharedPath animal '_' hem '_SharVars.mat'], 'SharVars');
    end 
end