%% function DetectMachine - ATW - May 2013

% Summary: Function retrieves computer specific variables, such as paths to
%   files.

function [machineVars] = DetectMachine

varpath = userpath;
varpath = varpath(1:end-1);
if exist([varpath '\machineVars.mat'], 'file') == 2;
    machineVars = load([varpath '\machineVars.mat']);
else
    SetupMachineVars;
    machineVars = load([varpath '\machineVars.mat']);
end