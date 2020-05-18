function ProcessIntriniscDirectories(Drivename,varargin)
if nargin==0
    Drivename=input('Select drive that contains data to be analyzed ("D:\")','s');
end
directoryContents=dir(Drivename);
folderNum=4:(size(directoryContents,1)-2);
theLED=input('Enter model of LED used for illumination e.g. M530L3. ','s');
bandFilter=input('Enter model of optical bandpass filter in light path. e.g FB570-10. ','s');
cutFilter=input('Enter model of optical cut off filter in light path. e.g FEL0500. ','s');
for k=1:length(folderNum)
    cd([directoryContents(folderNum(k)).folder  directoryContents(folderNum(k)).name]);
    tempDir=dir([directoryContents(folderNum(k)).folder  directoryContents(folderNum(k)).name]);
    for j=3:size(tempDir,1)
        cd([directoryContents(folderNum(k)).folder  directoryContents(folderNum(k)).name '\' tempDir(j).name]);
        thefiles=dir('*processeddata.mat');
        if isempty(thefiles)
        tempFiles=dir('*rawdata.mat');
        if ~isempty(tempFiles)
        for f=1:size(tempFiles,1)
            filenames{f}=tempFiles(f).name;
        end
        Chunk_IOS_CBV_Optogenetics_001(filenames,theLED,bandFilter,cutFilter);
        clear filenames
        else
            fprintf('No rawdata files in folder. Skipping\n');
        end
        else
            fprintf('Animal already analyzed. Skipping\n')
        end
    end
end
GroupCBVbyAnimal_WhiskerStim;
[SortedData]=AssembleGroupIOS(Drivename);
plotFig01(SortedData)
end


