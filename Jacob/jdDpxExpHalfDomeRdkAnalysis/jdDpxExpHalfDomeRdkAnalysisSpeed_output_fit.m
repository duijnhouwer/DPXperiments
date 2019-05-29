function [Speed,Log10Speed,Contrast,Log10Contrast,Yaw,Weights]=jdDpxExpHalfDomeRdkAnalysisSpeed_output_fit(DATA)
    
    % DATA = % output from jdDpxExpHalfDomeRdkAnalysisSpeed
    % 
    %   struct with fields:
    % 
    %      speeds: [1×140 double]
    %         yaw: [1×140 double]
    %      yawSem: [1×140 double]
    %     ctrlVar: [1×140 double]
    %       mouse: {1×140 cell}
    %           N: 140
    %
    % Fit a bi-variate, bi-cubic to Figures 2 and 3 of Laurens Mouse Psycho
    % half dome data
    
    
    DATA=dpxdSubset(DATA,strcmp(DATA.mouse,'MEAN'));
    
    Speed=DATA.speeds;
    Contrast=DATA.contrast;
    Yaw=DATA.yaw;
    Weights=1./DATA.yawSem;
    
    Log10Contrast=log10(Contrast);
    Log10Speed=log10(Speed);
    
    






   % fit 
   Bo=[-261 -76 323 -35 25 -90];
   Bhat=lsqcurvefit(@myBulge, Bo, {Log10Contrast Log10Speed}, Yaw, [], [], optimset('lsqcurvefit'));  
   
   
  
   cpsFindFig('yaw as function of speed and contrast'); 
   clf;
   plot3(Log10Contrast,Log10Speed,Yaw,'wo','MarkerFaceColor','k');
   hold on

   [XX,YY]=meshgrid(linspace(0,max(Log10Contrast),25),linspace(0,max(Log10Speed),25));
   ZZhat = myBulge(Bhat,{XX YY});
   surf(XX,YY,ZZhat);
   xlabel('Speed');
   ylabel('Contrast');
   zlabel('Yaw (deg/s)');
   
    
end
    

function Z=myBulge(B,XY)
    X=XY{1};
    Y=XY{2}; 
    Z= B(1) + B(2)*X + B(3)*Y + B(4)*X.^2 + B(5)*X.*Y + B(6)*Y.^2;
end
    
    

function [Speed, Contrast, Yaw, Weight]=addZeros(Speed, Contrast, Yaw, Weight)
    % add zero-yaw at zero-speed and at zero-contrast to further constrain
    % the bulge
    %
    % EXAMPLE:
    %     [Speed0, Contrast0, Yaw0, Weight0]=addZeros(Speed, Contrast, Yaw, Weight);
    
    nMice = numel(unique(DATA.mouse));
    nContrasts = numel(unique(DATA.contrast));
    nSpeeds = numel(unique(DATA.speeds));
    % add the zero-speeds
    Speed = [Speed zeros(1,nContrasts*nMice)];
    Contrast = [Contrast repmat(unique(Contrast),1,nMice)];
    Yaw = [Yaw zeros(1,nContrasts*nMice)];
    Weights = [Weights ones(1,nContrasts*nMice)*meanWeights];
    % add the zero-contrasts
    Speed = [Speed repmat(unique(DATA.speeds),1,nMice)];
    Contrast = [Contrast zeros(1,nSpeeds*nMice)];
    Yaw = [Yaw zeros(1,nSpeeds*nMice)];
    Weights = [Weights ones(1,nSpeeds*nMice)*meanWeights];  
end
    
