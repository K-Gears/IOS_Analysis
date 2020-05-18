function SetupMachineVars

% Set up Shared Path
SharedPath = uigetdir('C:\','Select the location of your shared variables (ROIs, Thresholds etc.)');
NASpath = uigetdir('C:\', 'Select the location of your NAS files');
if or(isempty(SharedPath),isempty(NASpath))
    error('No folder selected...aborting.')
end
if strcmp(SharedPath(end),'\') 
    machineVars.SharedPath = SharedPath;
else
    machineVars.SharedPath = [SharedPath '\'];
end
if strcmp(NASpath(end),'\')
    machineVars.NASpath = NASpath;
else
    machineVars.NASpath = [NASpath '\'];
end
machineVars.nameind = 1:17;
varpath = userpath;
%varpath = varpath(1:end-1);
prevdir = cd(varpath);
save('machineVars.mat','machineVars');
cd(prevdir);