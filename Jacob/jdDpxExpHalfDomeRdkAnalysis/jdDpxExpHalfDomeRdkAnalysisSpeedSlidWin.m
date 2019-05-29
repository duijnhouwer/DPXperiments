function jdDpxExpHalfDomeRdkAnalysisSpeedSlidWin(startMin,stepMin,winMin)
        
    % jdDpxExpHalfDomeRdkAnalysisSpeedSlidWin(startMin,stepMin,winMin)
    %
    % wrapper around jdDpxExpHalfDomeRdkAnalysisSpeed to perform a sliding
    % window analysis
    %
    % see also: jdDpxExpHalfDomeRdkAnalysisSpeed
    %
    % Jacob 20170314
    
    if nargin==0
        startMin=0;
        stepMin=1;
        winMin=15; % the width or duration of the time window in minutes
    end
    
    files=dpxUIgetFiles;
    disp([num2str(numel(files)) ' datafiles selected.']);
    if isempty(files)
        return;
    end
    
    for slid=0:1e30 % "infinite" steps
        timeWinSec=[startMin+slid*stepMin startMin+slid*stepMin+winMin]*60; 
        out=jdDpxExpHalfDomeRdkAnalysisSpeed(files,timeWinSec);
        if out.nFilesWithDataWithinWindow==0
            break; % the for loop
        end
    end
    
end