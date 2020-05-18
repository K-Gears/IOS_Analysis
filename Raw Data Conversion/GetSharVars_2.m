function SharVars = GetSharVars_2(animal,hem)
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Mar 2014
%   Version 1
%
%   SUMMARY: Loads animal and hemisphere-specific variables
%_______________________________________________________________
%   INPUTS:
%                           animal - [string] animal ID
%
%                           hem - [string] hemisphere recorded
%_______________________________________________________________
%   OUTPUTS:
%                           SharVars - [struct] contains variables specific
%                           to the animal and hemisphere
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%                           DetectMachine.m
%_______________________________________________________________
%   CALLED BY:
%                           GOF_Neur2CBV_WholeTrial
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

Q = DetectMachine_2;
if exist([Q.machineVars.SharedPath animal '_' hem '_SharVars.mat'],'file') == 2
    load([Q.machineVars.SharedPath animal '_' hem '_SharVars.mat']);
else
    CreateEmptySharedVariablesStruct(animal, hem);
    load([Q.machineVars.SharedPath animal '_' hem '_SharVars.mat']);
end

fieldnames = [{'Thresholds'},{'Baselines'},{'ROIs'},{'HRF'}];
for f = 1:length(fieldnames)
    if isfield(SharVars,fieldnames{f})==0
        SharVars.(fieldnames{f}) = [];
    end
end

