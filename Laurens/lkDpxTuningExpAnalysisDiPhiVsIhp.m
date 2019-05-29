function lkDpxTuningExpAnalysisDiPhiVsIhp(D,filterMode,normalize)
    
    % D should be the output from lkDpxTuningExpAnalysis
    [is,whynot] = dpxdIs(D);
    if ~is
        error(whynot)
    end
    
    if ~exist('filterMode') || isempty(filterMode)
        filterMode='includeall';
    end
     if ~exist('normalize') || isempty(normalize)
        normalize=true;
     end
        
    
    % Split D into a struct per cell recorded in each separate file (this
    % means that if the same cell was recorded in multiple files, they
    % will show up as separate units regardless, if this is a very common
    % tghing, we have to think up a way to merge them!
    C = dpxdSplit(D,'fileCellId');
    
    % calculate all Direction-Indices for phi (DIphi) and reverse phi
    % (DIihp) for each cell in each file
    DI.phi = [];
    DI.ihp = [];
    
    for i = 1:numel(C)
        [PHI,IHP] = dpxdSubset(C{i},strcmpi(C{i}.motType,'phi'));
        if normalize
            globalMean=nanmean([PHI.allDFoF{1}(:); IHP.allDFoF{1}(:)]);
            normPhi=PHI.allDFoF{1};
            %normPhi(:,1)=normPhi(:,1)./nanstd(
            
            
            normPhi=(PHI.allDFoF{1}-nanmean(all))./nanstd(PHI.allDFoF{1});
            normIhp=(IHP.allDFoF{1}-nanmean(all))./nanstd(IHP.allDFoF{1});            
            diPhi=diff(mean(normPhi,1));
            diIhp=diff(mean(normIhp,1));
        else
            diPhi=diff(PHI.meanDFoF{1});
            diIhp=diff(IHP.meanDFoF{1});
        end
        if strcmpi(filterMode,'phiSignificant')
            if ttest2(PHI.allDFoF{1}(:,1),PHI.allDFoF{1}(:,2),'alpha',0.1)
                DI.phi(end+1) = diPhi;
                DI.ihp(end+1) = diIhp;
            end
        elseif strcmpi(filterMode,'ihpsignificant')
            if ttest2(IHP.allDFoF{1}(:,1),IHP.allDFoF{1}(:,2),'alpha',0.1)
                DI.phi(end+1) = diPhi;
                DI.ihp(end+1) = diIhp;
            end
        elseif strcmpi(filterMode,'bothsignificant')
            if ttest2(PHI.allDFoF{1}(:,1),PHI.allDFoF{1}(:,2),'alpha',0.05,'vartype','unequal') && ttest2(IHP.allDFoF{1}(:,1),IHP.allDFoF{1}(:,2),'alpha',0.05,'vartype','unequal')
                DI.phi(end+1) = diPhi;
                DI.ihp(end+1) = diIhp;
            end
        elseif strcmpi(filterMode,'includeall')
            DI.phi(end+1) = diPhi;
            DI.ihp(end+1) = diIhp;
        else
            error('unknown filterMode');
        end
    end
    
    % plot the correlatiopn between DI.phi and DI.ihp. We hypothesize this
    % should be negative (because tuning to PHI and IHP should have
    % different sign)
    cpsFindFig(mfilename);
    plot(DI.phi,DI.ihp,'b.');
    axis equal
    a = axis;
    cpsRefLine('+','k--');
    h = refline;
    axis(a);
    set(h,'Color','r','LineWidth',2);
    [r,pVal]=corr(DI.phi(:),DI.ihp(:));
    h = cpsText({['Pearson r = ' num2str(r,'%.2f')],['N = ' num2str(numel(DI.phi))],['p = ' num2str(pVal,'%.3f')]});
    set(h,'FontSize',10)
    xlabel('phi[\theta] - phi[\theta+180] (\DeltaF/F)','FontSize',14);
    ylabel('revPhi[\theta] - revPhi[\theta+180] (\DeltaF/F)','FontSize',14);    
end
