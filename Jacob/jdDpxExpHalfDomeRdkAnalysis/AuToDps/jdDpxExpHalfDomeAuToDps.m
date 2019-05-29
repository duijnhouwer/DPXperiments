function scalar = jdDpxExpHalfDomeAuToDps
    
    % in the movie in this folder,
    % ...\DPXperiments\Jacob\jdDpxExpHalfDomeRdkAnalysis\AuToDps\IMG_8642.m4v
    % you can see the ball make almost 7 rolls, stopping short at
    % about 11 o'clock (if 12 o'clock is the top). This boils down to a
    % total revolution of 7*360-30 = 2490 degress
    
    totalDeg=2490;
    
    % the computer mouse output has been recorded too, load it now. It's u
    
    DPXD=dpxdLoad('jdDpxExpHalfDomeTransStaticNoise-WJ-20170704154103.mat');
    
    % this was the fourth trial of this session, and this axis of rotation
    % (roll) is recorded by DPXD.resp_mouseSide_dyPx. Load that now
    
    roll=DPXD.resp_mouseSide_dyPx{4};
    
    % I've looked at this data, and the ball start rolls between samples 385
    % to 916. Select that now
    
    ioi=385:916; % interval of interest
    roll=roll(ioi);
    
    % let's figure out how many seconds that was. looking from the movie,
    % it should be around 9 seconds. But let's make a more accurate
    % estimate using the sampling rate of the mice (same framerate as the
    % stimulus
    
    sampleHz=DPXD.window_measuredFrameRate(4);
    
    % that means the ball rolled for 
    
    rollSeconds = numel(ioi)/sampleHz; % 8.87
    
    % thus, the speed of the ball;
    
    ballDps = totalDeg/rollSeconds; % 281 dps
    
    % the speed in resp_mouseSide_dyPx is in Px per Frame. Now we have
    % enough info to convert px per frame (what we used to call arbitrary
    % units) to deg per second. This scalar is 
    
    scalar = ballDps/nanmean(roll); %1.77
    
    % i.e., multiply the pixel per frame value by scalar to get deg/s
    
    
    
    
    
    
    