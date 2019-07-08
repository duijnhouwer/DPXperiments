function DPXD=datafiles_to_barebones_dpxd(files)
    
    % Analyze halfdome mouse on ball data
    % see also:
    %    jdDpxExpHalfDomeRdkAnalysisSpeedSlidWin
    %    jdDpxExpHalfDomeRdkAnalysisSpeedEarlyLate
    %

    if ~exist('files','var') || isempty(files)
        files=dpxUIgetFiles;
        disp([num2str(numel(files)) ' datafiles selected.']);
        if isempty(files)
            return;
        end
    end
    
    response_interval_sec=-1:1/60:3; % (60Hz)
    response_interval_sec([1 end])=[]; % remove first and last because often missing (would be nan value at 1 and 3 otherwise)
    
    E={};    
    for i=1:numel(files)
        fprintf('[%s] Loading %s...\n',mfilename,files{i});
        D=dpxdLoad(files{i});
        
        % check the subject id and fix if necessary
        id_dpxd=cell2mat(unique(D.exp_subjectId));
        [~,fname]=fileparts(files{i});
        id_filename = fname(find(fname=='-',1,'first')+1:find(fname=='-',1,'last')-1);
        if ~strcmpi(id_dpxd,id_filename)
            fprintf('WARNING! ID in DPXD is %s but in filename it''s %s!!!\n',id_dpxd,id_filename);
            if false
                fprintf('   Assuming it was entered wrong when recording and then corrected later in the fileNAME only\n');
                D.exp_subjectId=repmat({id_filename},1,D.N);
                fprintf('   Replaced the subject ID in the dpxd with the ID from the filename...\n');
            end
            continue
        end

        [D,suspect]=clarifyAndCheck(D,response_interval_sec);
        if ~suspect
            % Only include data files that are completely fine and have no suspicious
            % things happening, like poor correlations between the yaw measurements of
            % both computer mice on the ball (both measure yaw, should be about the
            % same)
            E{end+1}=D; %#ok<AGROW>
        else
            disp('clarifyAndCheck said file is suspect');
            fprintf(' ---> skipping file : %s\n', files{i} );
        end
        clear D;
    end
    if isempty(E)
        disp('None of the files passed the tests. Can''t continue.');
        return;
    else
        disp([num2str(numel(E)) ' out of ' num2str(numel(files)) ' data files passed the tests.']);
    end
    %
    % Merge all datafiles that we collected in cell array
    E=dpxdMerge(E); % E is now a DPXD
    E=barebonesify(E);
    DPXD=E; % Bare bones reverse phi data
    [fname,folder] = uiputfile('*.mat', 'Bare bones reverse phi data as...','barebonesreversephidata');
    if ~isnumeric(fname)
        save(fullfile(folder,fname),'DPXD');
    end
end

function [D,suspect]=clarifyAndCheck(D,response_interval_sec)
    % Make some changes to the DPXD that make the analysis easier to read;

    % Step 0 remove the first sample (outlier probably because mouse
    % cursor not moved to the center of the screen yet)
    for t=1:D.N
        D.resp_mouseBack_dxPx{t}(1)=[];
        D.resp_mouseSide_dxPx{t}(1)=[];
        D.resp_mouseBack_dyPx{t}(1)=[];
        D.resp_mouseSide_dyPx{t}(1)=[];
        D.resp_mouseBack_tSec{t}(1)=[];
        D.resp_mouseSide_tSec{t}(1)=[];
    end
    % Step 1, align time of traces to the start of the motion
    for t=1:D.N
        D.resp_mouseBack_tSec{t}=D.resp_mouseBack_tSec{t}-D.startSec(t)-D.rdk_motStartSec(t);
        D.resp_mouseSide_tSec{t}=D.resp_mouseSide_tSec{t}-D.startSec(t)-D.rdk_motStartSec(t);
    end
    % Step 2, remove offset from X value traces, because of monitor
    % settings in the Half Dome setup, the left-x is 0, and the control
    % computer starts at -1920. The Logitech mice are sampled on the
    % control monitor.
    for t=1:D.N
        D.resp_mouseBack_dxPx{t}=D.resp_mouseBack_dxPx{t}+1920;
        D.resp_mouseSide_dxPx{t}=D.resp_mouseSide_dxPx{t}+1920;
    end
    % Step 3, rename the mouse fields that code yaw. Also get the pitch
    % (forward speed) for each trial so we can filter on running speed. and
    % get the correlation between the yaw as measured with the sidemouse
    % and the mouse in the back. These correlation should be very high
    % because they measrue the same thing. After calculating the correlaton
    % apply interpolation to get the data on the specified reponse_interval
    for t=1:D.N
        D.resp_yaw_corr{t}=corr(D.resp_mouseBack_dyPx{t}(:),D.resp_mouseSide_dyPx{t}(:));
        D.resp_tSec{t}=response_interval_sec;
        tsec=mean([D.resp_mouseBack_tSec{t}(:) D.resp_mouseSide_tSec{t}(:)],2)'; 
        backyaw=interp1(tsec,D.resp_mouseBack_dyPx{t},response_interval_sec,'linear');
        sideyaw=interp1(tsec,D.resp_mouseSide_dyPx{t},response_interval_sec,'linear');
        D.resp_yaw{t}=mean([backyaw(:) sideyaw(:)],2)';
        D.resp_pitch{t}=interp1(tsec,D.resp_mouseBack_dxPx{t},response_interval_sec,'linear');
        D.resp_roll{t}=interp1(tsec,D.resp_mouseSide_dxPx{t},response_interval_sec,'linear');
    end
    % Step 4: Convert yaw pixels/frame to deg/s (added 20170710)
    scalar = jdDpxExpHalfDomeAuToDps;
    for i=1:numel(D.resp_yaw)
        D.resp_yaw{i}=D.resp_yaw{i}*scalar;
        D.resp_pitch{i}=D.resp_pitch{i}*scalar;
        D.resp_roll{i}=D.resp_roll{i}*scalar;
    end
    % See if the file is up to snuff
    suspect = false; % no checks performed as of yet
end

function B=barebonesify(E)
    % add a 'stimmode' field for the stimulus mode (phi, reverse phi,
    % unlimited lifetime)
    E_phi=dpxdSubset(E,(E.rdk_nSteps==1 & E.rdk_invertSteps==Inf));
    E_phi.stimmode=repmat('p',1,E_phi.N);
    E_revphi=dpxdSubset(E,(E.rdk_nSteps==1 & E.rdk_invertSteps==1));
    E_revphi.stimmode=repmat('r',1,E_revphi.N);
    E_unlim=dpxdSubset(E,(E.rdk_nSteps==Inf));
    E_unlim.stimmode=repmat('u',1,E_unlim.N);
    E=dpxdMerge({E_phi,E_revphi,E_unlim});
    clear('E_phi','E_revphi','E_unlim');
    %
    % Replace 'M003' style subject names into numbers
    E.exp_subjectId=cellfun(@(x)str2double(x(2:end)),E.exp_subjectId);
    %
    % Split the data by subject, stimmode, freezeframe,rdk_aziDps
    E=dpxdSplit(E,{'exp_subjectId','stimmode','rdk_freezeFlip','rdk_aziDps'});
    % 
    % Copy only the relevant data into B
    B=cell(size(E));
    for i=1:numel(E)
        B{i}.mouse=E{i}.exp_subjectId(1);
        B{i}.ff=E{i}.rdk_freezeFlip(1);
        B{i}.mode=E{i}.stimmode(1);
        B{i}.dps=E{i}.rdk_aziDps(1); % stimulus speed in deg per second
        B{i}.ms=E{i}.resp_tSec{1}(:)*1000;
        K=cellfun(@transpose,E{i}.resp_yaw,'UniformOutput',false);
        B{i}.yaw={[K{:}]}; % 1 cell with matrix with columns being trials, rows being time-samples of yaw
        K=cellfun(@transpose,E{i}.resp_pitch,'UniformOutput',false);
        B{i}.pitch={[K{:}]}; % 1 cell with matrix with columns being trials, rows being time-samples of pitch
        K=cellfun(@transpose,E{i}.resp_roll,'UniformOutput',false);
        B{i}.roll={[K{:}]}; % 1 cell with matrix with columns being trials, rows being time-samples of roll
        B{i}.framedropspersec={E{i}.nrMissedFlips./(E{i}.stopSec-E{i}.startSec)};
        B{i}.N=1;
    end 
    B=dpxdMerge(B);
end

