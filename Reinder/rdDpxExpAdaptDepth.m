function rdDpxExpAdaptDepth(varargin)
% Displays binocular rivalry stimulus, and after a set time brings up
% kinetic depth cylinders.
%
% Input Argument:
%   Type: 'diep' / 'bind'
%       respectively depth or binding experiment
        

%%%%%%%%%%%%%%%%%%
% STIMULUS INPUT %
%%%%%%%%%%%%%%%%%%
global IN

% adaptation
IN.adapSec     = 5;                 % time for adaptation
IN.warningSec  = 4;                 % stimulus warning, seconds before cylinder starts
IN.adapSize    = 8;

% cylinders
IN.cylRepeats  = 20;                % number of repeats of the cylinder stimulus
IN.cylOnSec    = 1;                 % Ton for cylinders
IN.cylOffSec   = 1.5;               % Toff for cylinder, in this case again adaptation
IN.disparities = [-.4 -.2 0 .2 .4]; % disparities used for the cylinders
IN.cylPosition = 1;                 % cylinder position, halve left, full right or vise versa
IN.rotSpeed    = [120 -120];        % speeds used by the cylinder in degree per second
IN.modes       = 'stereo';          % mode of depth in the stimulus, 'stereo' for disparity

%%%%%%%%%%%%%%%%%%%%%
%   START STUFF     %
%%%%%%%%%%%%%%%%%%%%%
if nargin~=1
  warning('wrong input.')
  disp('Running ''Diep''.\n')
end  

E=dpxCoreExperiment;
E.paradigm      = mfilename;
E.window.set('scrNr',0,'rectPx',[],'stereoMode','mirror'); % 'rectPx',[1440 0 1600+1440 1200]
E.window.set('distMm',1000,'interEyeMm',65,'widHeiMm',[394 295]);
E.window.set('gamma',0.49,'backRGBA',[.5 .5 .5 1],'skipSyncTests',1);

%prepare type specific stuff
switch lower(varargin{1})
    case 'diep'
        E.outputFolder  = 'C:\DPXDTemp\AdaptDiepte\';
        E.txtStart      = 'Diepte experiment\n\nStaar naar het kruisje. \nEerst ziet u enkele minuten een adaptatie stimulus, waarna een \ndraaiende cylinder verschijnt\n\nRapporteer diepte van linker cylinder\n\nPijltje omhoog  = Hol\nPijltje omlaag = bol ';
    case 'bind'
        E.outputFolder  = 'C:\DPXDTemp\AdaptDiepte\';
        E.txtStart      = 'Binding experiment\n\nStaar naar het kruisje. \nEerst ziet u enkele minuten een adaptatie stimulus, waarna een \ndraaiende cylinder verschijnt\n\nRapporteer bewegingsrichting van voorvlak van rechter cylinder\n\nPijltje omhoog = omhoog\nPijltje omlaag = omlaag.';
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%
%   FIRST ADAPTATION    %
%%%%%%%%%%%%%%%%%%%%%%%%%
adapC=dpxCoreCondition;
textC = dpxStimText;
set(textC,'str',sprintf('Cylinder stimulus starts in %d seconds',IN.warningSec),...
    'vAlign',1*E.window.deg2px,'onSec',IN.adapSec - IN.warningSec);
adapC.addStimulus(textC);
adapC = defineAdaptationStimulation(E.window,'adap',adapC);
adapC = defineCylinderStimulinder(false,adapC,0,0);

E.addCondition(adapC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   AFTER 1800 SEC CYLINDER STIMULUS  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for disparity = 1:numel(IN.disparities)
    for speeds = 1:numel(IN.rotSpeed)
    cylC = dpxCoreCondition;
    cylC = defineCylinderStimulinder(true,cylC,IN.disparities(disparity),IN.rotSpeed(speeds));
    cylC = defineAdaptationStimulation(E.window,'cyl',cylC);
    textC = dpxStimText; set(textC,'enabled',0'); cylC.addStimulus(textC); %placeholder
    E.addCondition(cylC);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%% 
% ALL STIMULI COLLECTED %
%%%%%%%%%%%%%%%%%%%%%%%%%
%        ! RUN !        %
%%%%%%%%%%%%%%%%%%%%%%%%%

cylSequence = repmat(2:numel(IN.disparities)+1,1,IN.cylRepeats);

E.conditionSequence=[1 cylSequence(randperm(numel(cylSequence)))];

E.run;

end

function C = defineAdaptationStimulation(W,state,C)
global IN

if strcmp(state,'adap'); 

C.durSec    = IN.adapSec;
S=dpxStimCross;
set(S,'wDeg',.2,'hDeg',.2,'lineWidDeg',.05,'onSec',0,'name','fix','visible',1);
C.addStimulus(S);

%left stimulus
ML = dpxStimMaskGaussian;
ML.name='MaskLeft';
ML.xDeg=0;
ML.hDeg = IN.adapSize*sqrt(2); %(50*sqrt(2))/W.deg2px;
ML.wDeg = IN.adapSize*sqrt(2); %(50*sqrt(2))/W.deg2px;
ML.sigmaDeg = ML.hDeg/8;
ML.durSec=IN.adapSec;
C.addStimulus(ML);

GL = dpxStimGrating;
GL.name = 'gratingLeft';
GL.xDeg=0;
GL.dirDeg=-45;
GL.contrastFrac=1;
GL.squareWave=false;
GL.cyclesPerSecond=0;
GL.cyclesPerDeg=2.5;
GL.wDeg= IN.adapSize; %(50)/W.deg2px;
GL.hDeg= IN.adapSize; %(50)/W.deg2px;
GL.durSec=IN.adapSec;
GL.buffer=0;
C.addStimulus(GL);

%right stimulus
MR = dpxStimMaskGaussian;
MR.name='MaskRite';
MR.xDeg=0;
MR.hDeg = IN.adapSize*sqrt(2); % (50*sqrt(2))/W.deg2px;
MR.wDeg = IN.adapSize*sqrt(2); % (50*sqrt(2))/W.deg2px;
MR.sigmaDeg = MR.hDeg/8;
MR.durSec=IN.adapSec;
C.addStimulus(MR);

GR = dpxStimGrating;
GR.name = 'gratingRight';
GR.xDeg=0;
GR.dirDeg=45;
GR.squareWave=false;
GR.cyclesPerSecond=0;
GR.cyclesPerDeg=2.5;
GR.wDeg= IN.adapSize; % (50)/W.deg2px;
GR.hDeg= IN.adapSize; % (50)/W.deg2px;
GR.durSec=IN.adapSec;
GR.buffer=1;
C.addStimulus(GR);

elseif strcmp(state,'cyl'); 
C.durSec    = IN.cylOnSec+IN.cylOffSec;
%left stimulus
ML = dpxStimMaskGaussian;
ML.name='MaskLeft';
ML.xDeg=0;
ML.hDeg = IN.adapSize*sqrt(2);
ML.wDeg = IN.adapSize*sqrt(2);
ML.sigmaDeg = ML.hDeg/8;
ML.onSec=IN.cylOnSec;
ML.durSec=IN.cylOffSec;
C.addStimulus(ML);

GL = dpxStimGrating;
GL.name = 'gratingLeft';
GL.xDeg=0;
GL.dirDeg=-45;
GL.contrastFrac=1;
GL.squareWave=false;
GL.cyclesPerSecond=0;
GL.cyclesPerDeg=2.5;
GL.wDeg=IN.adapSize;
GL.hDeg=IN.adapSize;
GL.onSec=IN.cylOnSec;
GL.durSec=IN.cylOffSec;
GL.buffer=0;
C.addStimulus(GL);

%right stimulus
MR = dpxStimMaskGaussian;
MR.name='MaskRite';
MR.xDeg=0;
MR.hDeg = IN.adapSize*sqrt(2);
MR.wDeg = IN.adapSize*sqrt(2);
MR.sigmaDeg = MR.hDeg/8;
MR.onSec=IN.cylOnSec;
MR.durSec=IN.cylOffSec;
C.addStimulus(MR);

GR = dpxStimGrating;
GR.name = 'gratingRight';
GR.xDeg=0;
GR.dirDeg=45;
GR.squareWave=false;
GR.cyclesPerSecond=0;
GR.cyclesPerDeg=2.5;
GR.wDeg= IN.adapSize;
GR.hDeg= IN.adapSize;
GR.onSec=IN.cylOnSec;
GR.durSec=IN.adapSec;
GR.buffer=1;
C.addStimulus(GR);
end
end

function C = defineCylinderStimulinder(state,C,disp,speed)
global IN

% The feedback stimulus for correct responses
S=dpxStimDot;
set(S,'wDeg',.25,'enabled',false,'durSec',.1,'RGBAfrac',[.75 .75 .75 .75],'name','fbCorrect');
C.addStimulus(S);

if state; visible = 1; C.durSec = IN.cylOnSec+IN.cylOffSec;
    % fixation cross must be on top
    S=dpxStimCross;
    set(S,'wDeg',.2,'hDeg',.2,'lineWidDeg',.05,'onSec',0,'name','fix','visible',1);
    C.addStimulus(S);
elseif ~state; visible = 0; C.durSec    = IN.adapSec;
end

% The full cylinder stimulus
S=dpxStimRotCylinder;
set(S,'dotsPerSqrDeg',12,'xDeg',IN.cylPosition*1.75,'wDeg',3,'hDeg',3,'dotDiamDeg',0.11 ...
    ,'rotSpeedDeg',speed,'disparityFrac',0,'sideToDraw','both' ...
    ,'onSec',0,'durSec',IN.cylOnSec,'stereoLumCorr',1,'fogFrac',0,'dotDiamScaleFrac',0 ...
    ,'name','fullTargetCyl','visible',visible);
C.addStimulus(S);

S=dpxStimRotCylinder;
set(S,'dotsPerSqrDeg',12,'xDeg',IN.cylPosition*-1.75,'wDeg',3,'hDeg',3,'dotDiamDeg',0.11 ...
    ,'rotSpeedDeg',speed,'disparityFrac',disp,'sideToDraw','front' ...
    ,'onSec',0,'durSec',IN.cylOnSec,'name','halfInducerCyl','visible',visible);
C.addStimulus(S);

% The response object
R=dpxRespKeyboard;
R.name='rightHand';
if state; R.allowAfterSec=0;
elseif ~state; R.allowAfterSec=IN.adapSec; end
R.kbNames='UpArrow,DownArrow';
R.correctStimName='fbCorrect';
R.correctKbNames='1';
R.correctEndsTrialAfterSec=inf;
R.wrongEndsTrialAfterSec=inf;
C.addResponse(R);
end


