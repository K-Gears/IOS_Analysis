function [GlobVars] = GetGlobVars
%   Written by Aaron Winder, Drew Lab, ESM, Penn State University, Mar 2014
%   Version 1
%
%   SUMMARY: Loads Global Variables into Memory
%_______________________________________________________________
%   INPUTS:
%                           None
%_______________________________________________________________
%   OUTPUTS:
%                           GlobVars - [Struct] Variables that are applied
%                           across all analysis. 
%_______________________________________________________________
%   REQUIRED SCRIPTS:
%                           Q = DetectMachine;
%_______________________________________________________________
%   CALLED BY:
%_______________________________________________________________
%   FUTURE VERSIONS:
%_______________________________________________________________
%   CHANGES FROM PREV VERS:
%_______________________________________________________________

Q = DetectMachine;
load([Q.machineVars.SharedPath 'Global_Variables.mat'])