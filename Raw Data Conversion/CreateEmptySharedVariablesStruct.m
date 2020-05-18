function CreateEmptySharedVariablesStruct(animal, hem)

%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Aug 2013
%   Version 1

%   SUMMARY: Creates an empty structure of shared variables so the
%   variables can be added.
%_______________________________________________________________
%   INPUTS:
%                           animal - animal name as a string

%                           hem - hemisphere recorded as a string
%                                   ('LH' 'RH')
%_______________________________________________________________
%   OUTPUTS:  None, the result is a structure saved to the shared variables
%   folder
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%           [Q] = DetectMachine(prevdir)
%_______________________________________________________________
%   CALLED BY:
%           DataNormalize.m
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

Q = DetectMachine_2;
SharVars.ROIs.Description = 'ROIs are coordinates of a spatial region of interest of some image';
SharVars.Baselines.Description = 'Baselines are reference values for normalization of some measurement';
SharVars.Thresholds.Description = 'Thresholds are a reference used to binarize some measurement';
SharVars.HRF.Description = 'HRFParams are parameters needed to reconstruct a hemodynamic response function (HRF)';
save([Q.machineVars.SharedPath animal '_' hem '_SharVars.mat'],'SharVars');