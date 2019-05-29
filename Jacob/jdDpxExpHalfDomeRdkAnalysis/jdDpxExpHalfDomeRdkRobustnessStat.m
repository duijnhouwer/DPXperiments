function jdDpxExpHalfDomeRdkRobustnessStat(P)
    
    %     P has to be a :
    %
    %   struct with fields:
    %
    %       speeds: [1×192 double]
    %          yaw: [1×192 double]
    %       yawSem: [1×192 double]
    %      ctrlVar: [1×192 double]
    %        mouse: {1×192 cell}
    %     medSplit: {1×192 cell}
    %            N: 192
    %
    % Jacob 2017-10-31
    
    do_ttest(P,'between');
    do_ttest(P,'within');
    
end




function pval=do_ttest(P,split)

      
    cpsFindFig(split)
    
    if ~any(strcmp(split,{'between','within'}))
        error('split must be between or within');
    end
    late=dpxdSubset(P,strcmp(P.medSplit,[split 'Late']));
    early=dpxdSubset(P,strcmp(P.medSplit,[split 'Early']));
  
    
   subplot(2,1,1);
    contrasts=nan(1,late.N);
    yaw_late_array=contrasts;
    yaw_early_array=contrasts;
    for i=1:late.N
        yaw_late=late.yaw(i);
        K=dpxdSubset(early,strcmp(early.mouse,late.mouse{i}) & early.speeds==late.speeds(i));
        if K.N~=1
            error('not unique??');
        end
        yaw_early=K.yaw;
        contrasts(i)=(yaw_late-yaw_early)/(yaw_late+yaw_early);
        
        plot(K.speeds,yaw_late,'ro');
        hold on
        plot(K.speeds,yaw_early,'bo');
        
        yaw_late_array(i)=yaw_late;
        yaw_early_array(i)=yaw_early;
    end
    title(upper(split));
    
  % subplot(2,1,2);
  %  hist(contrasts,25);
  %  title([split ' this should look normal!'])
  %  xlabel('(OMRlate-OMRearly)/(OMRlate+OMRearly)');
  %  ylabel('N');
  %  [h,pval,ci,stats] = ttest(contrasts,0,'tail','both');
  %%  if h
   %     ttxt='were significantly different from zero';
   %% else
  %      ttxt='were not different from unity';
  %  end 
  %  sprintf('The early-vs-late contrasts in the %s session comparison %s (M=%.3f, SD=%.3f, t(%d)=%.3f, p=%.3f).',upper(split),ttxt,mean(contrasts),stats.sd,stats.df,stats.tstat,pval)
   
  subplot(2,1,2)
  plot(yaw_early_array,yaw_late_array,'+')
  cpsRefLine('/');
  [p,~,stats]=ranksum(yaw_early_array,yaw_late_array);
  text(.01,.9,sprintf('Two sided Wilcoxon rank sum test: zval=%.3f,ranksum=%d,p=%.3f',stats.zval,stats.ranksum,p),'Units','Normalized');
  xlabel('OMRearly (deg/s)')
  ylabel('OMRlate (deg/s)')
end
