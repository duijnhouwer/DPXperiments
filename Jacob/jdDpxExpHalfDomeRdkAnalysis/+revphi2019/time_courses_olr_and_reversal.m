function time_courses_olr_and_reversal
    
    options.datafile='barebonesreversephidata_resampled16ms.mat'; % created with datafiles_to_barebones_dpxd
    options.pitch_range=[1/(19.7*pi/360) Inf];  % X cm/second  = X/(19.7*pi/360) deg/second because circumference of ball is 19.7*pi cm
    options.roll_range=[-Inf Inf];
    options.max_framedrops_per_second=6;
    options.min_trials_per_condition=35;
    options.splineroughness=1e-9;
    options.splinestep_ms=100;
    options.min_spline_r2=0;
    options.detrend_per_trace_method='mean'; % 'mean','linear'
    options.detrend_per_mouse_method='mean_limited_lifetime_yaw'; % 'none','mean_unlimited_lifetime_yaw','mean_limited_lifetime_yaw','linear_prestimon'
    options.nrows=3;
    options.include_mice=1:9%1:9; % 0 and 1:9 mean all mice, 1 means mouse 1 etc, [1 2 4 5] means these mice
    options.freezeflips=[1 2 3 4 5 6 7]; % freezeflips to keep
    options.pool_freezeflips=false;
    options.reversal_measure='signed absolute product'; % tstat, linear, multiplicative, divisive, signrank, ttest
    options.ttest_interval=[1 2];
    options.ttest_alpha=0.05;
    
    ttt=tic;
    fprintf('[%s] loading %s ...',mfilename, options.datafile);
    D=dpxdLoad(options.datafile); % should be created with revphi2019.datafiles_to_barebones_dpxd.m
    toc(ttt);
    
    D=remove_bad_trials(D,options.max_framedrops_per_second);
    % remove the unlimited lifetime data
    % [D,UNLIM]=dpxdSubset(D,D.mode~='u');
    
    % remove the freezeflips that we don't want to analyse
    D=dpxdSubset(D,ismember(D.ff,options.freezeflips));
    if options.pool_freezeflips
        D.ff=repmat(666,1,D.N);
    end
    
    % keep only the mice we want to analyze (include_mice=0 means all)
    if options.include_mice
        D=dpxdSubset(D,ismember(D.mouse,options.include_mice));
    end
    
    D=subset_trials_by_pitch_and_roll(D,options);
    [D,options.include_mice]=remove_mice_conditions_with_too_few_n(D,options.min_trials_per_condition); % good to remove outliers now so no spline fitting required, save time
    
    D=convert_yaw_per_trial_to_yaw_per_condition(D,options);
    [D,options.include_mice]=remove_mice_conditions_with_too_few_n(D,options.min_trials_per_condition); % remove more outliers (if any)
    
    
    [D,UNLIM]=dpxdSubset(D,D.mode~='u');
    
    % Detrend the yaw data by removing the mean over all conditions for
    % each mouse.
    fprintf('[%s] Detrending yaw per mouse using method: %s',mfilename,options.detrend_per_mouse_method)
    if strcmpi(options.detrend_per_mouse_method,'mean_unlimited_lifetime_yaw')
        D=dpxdSplit(D,{'mouse'});
        for i=1:numel(D)
            this=dpxdSubset(UNLIM,UNLIM.mouse==D{i}.mouse(1) & UNLIM.dps==0);
            D{i}.yaw_mouse_trend=repmat(mean(this.yaw_mean,2),1,D{i}.N);
            D{i}.yaw_mean=D{i}.yaw_mean-D{i}.yaw_mouse_trend;
        end
        D=dpxdMerge(D);
    elseif strcmpi(options.detrend_per_mouse_method,'mean_limited_lifetime_yaw')
        D=dpxdSplit(D,{'mouse'});
        for i=1:numel(D)
            D{i}.yaw_mouse_trend=repmat(mean(D{i}.yaw_mean,2),1,D{i}.N);
            D{i}.yaw_mean=D{i}.yaw_mean-D{i}.yaw_mouse_trend;
        end
        D=dpxdMerge(D);
    elseif ~strcmpi(options.detrend_per_mouse_method,'none')
        error('unknown detrend_per_mouse_method');
    end
    
    % plot mean OLRs
    cpsFindFig([mfilename ' - ' options.datafile]);
    clf
    D=dpxdSplit(D,'ff');
    n_freeze_flips=numel(D);
    for i=1:n_freeze_flips
        [PHI,REV]=dpxdSubset(D{i},D{i}.mode=='p');
        h=subplot(options.nrows,n_freeze_flips,i);
        plot_olrs(h,PHI,options,'left');
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel('OLR (deg/s)');
        end
        h=subplot(options.nrows,n_freeze_flips,i+n_freeze_flips);
        plot_olrs(h,REV,options,'right');
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel('OLR (deg/s)');
        end
        h=subplot(options.nrows,n_freeze_flips,i+2*n_freeze_flips);
        [line_handles,olr_per_mouse(i)]=plot_right_minus_left(h,D{i},options);
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel('OLR (deg/s)');
            legend(line_handles,{'Phi','Reverse Phi'});
        end
    end
    cpsUnifyAxes
    
    % plot_diff_phi_contrast_inverted_revphi_contrast(D,options,4)
    suptitle(['Mouse' sprintf(' %d',options.include_mice)]);
    
    %
    plot_mean_olr_over_freezeframes(olr_per_mouse,n_freeze_flips,D)
end

function D=remove_bad_trials(D,maxframedroprate)
    fprintf('[%s] removing trials with completely still ball or too many framedrops...',mfilename);
    n_bad=0;
    n_total=0;
    for i=1:D.N
        ok=nanstd(D.yaw{i})>0 & D.framedropspersec{i}<maxframedroprate;
        D.yaw{i}=D.yaw{i}(:,ok);
        D.pitch{i}=D.pitch{i}(:,ok);
        D.roll{i}=D.roll{i}(:,ok);
        D.framedropspersec{i}=D.framedropspersec{i}(ok);
        %
        n_bad=n_bad+sum(~ok);
        n_total=n_total+numel(ok);
    end
    fprintf(' (%d/%d trials removed)\n',n_bad,n_total);
end

function D=subset_trials_by_pitch_and_roll(D,options)
    % pitch = forward running
    for i=1:D.N
        ok_pitch=nanmean(D.pitch{i},1)>min(options.pitch_range) & nanmean(D.pitch{i},1)<max(options.pitch_range);
        ok_roll=nanmean(D.roll{i},1)>min(options.roll_range) & nanmean(D.roll{i},1)<max(options.roll_range);
        ok=ok_pitch & ok_roll;
        D.yaw{i}=D.yaw{i}(:,ok);
        fprintf('[%s] removed %s from condition %d\n',mfilename,jdProp(~ok,'%'),i);
    end
end

function [D,good_mice]=remove_mice_conditions_with_too_few_n(D,min_trials_per_condition)
    D=dpxdSplit(D,'mouse');
    remove=false(size(D));
    good_mice=[];
    
    % 20191108, make a bargraph of the number of trails each mouse has
    % completed
    n_mice=numel(D);
    n_reverse_phi_trials_per_mouse=[];
    n_phi_trials_per_mouse=[];
    for i=1:n_mice
        PHI=dpxdSubset(D{i},D{i}.mode=='p');
        IHP=dpxdSubset(D{i},D{i}.mode=='r');
        n_phi_trials_per_mouse(i)=sum(cellfun(@(x)(size(x,2)),PHI.yaw));
        n_reverse_phi_trials_per_mouse(i)=sum(cellfun(@(x)(size(x,2)),IHP.yaw));
    end
    
    
    cpsFindFig(sprintf('[%s] remove_mice_conditions_with_too_few_n',mfilename));
    h=bar([n_phi_trials_per_mouse(:) n_reverse_phi_trials_per_mouse(:)]);
    legend(h,{'Phi','Reverse-Phi'},'FontSize',12);
    xlabel('Mouse ID','FontSize',16);
    ylabel('Nr of trials','FontSize',16);
    
    
    for i=1:numel(D)
        if any(cellfun(@(x)(size(x,2)),D{i}.yaw)<min_trials_per_condition)
            remove(i)=true;
            fprintf('[%s] discarding mouse %d because one or more conditions had fewer that %d remaining trials\n',mfilename,unique(D{i}.mouse),min_trials_per_condition);
        else
            good_mice(end+1)=unique(D{i}.mouse); %#ok<AGROW>
        end
    end
    D(remove)=[];
    D=dpxdMerge(D);
end


function D=convert_yaw_per_trial_to_yaw_per_condition(D,options)
    warning('off','SPLINES:CHCKXYWP:NaNs');
    ttt=tic;
    if strcmpi(options.detrend_per_trace_method,'mean')
        fprintf('[%s] removing pre-stim means from yaw traces... ',mfilename);
        for i=1:D.N
            mean_yaw_prestim=nanmean(D.yaw{i}(D.ms(:,i)>=-Inf & D.ms(:,i)<=0,:),1);
            D.yaw{i}=D.yaw{i}-mean_yaw_prestim;
        end
    elseif strcmpi(options.detrend_per_trace_method,'linear')
        fitObj=fittype('poly1');
        opts=fitoptions('poly1');
        opts.Robust='on';
        warning('off','curvefit:fit:iterationLimitReached');
        for i=1:D.N
            n_traces=size(D.yaw{i},2);
            fprintf('[%s] detrending %d yaw traces of condition %d based on slope and offset of prestim interval...\n',mfilename,n_traces,i);
            ioi=D.ms(:,i)>=-Inf & D.ms(:,i)<=0;
            xx=D.ms(ioi,i);
            for t=1:size(D.yaw{i},2)
                yy=D.yaw{i}(ioi,t);
                if any(isnan(yy))
                    opts.Weights=ones(size(yy));
                    opts.Weights(isnan(yy))=0;
                    yy(isnan(yy))=0;
                end
                [fitresult,gof] = fit(xx,yy,fitObj,opts);
                yHat=fitresult(D.ms(:,i)); % extrapolate to entire time range
                D.yaw{i}(:,t)=D.yaw{i}(:,t)-yHat(:);
            end
        end
        warning('on','curvefit:fit:iterationLimitReached');
    elseif ~strcmpi(options.detrend_per_trace_method,'none')
        error('unknown detrend_per_trace_method: %s',options.detrend_per_trace_method);
    end
    toc(ttt);
    %
    ttt=tic;
    fprintf('[%s] replacing yaw traces with splines...',mfilename);
    total_n_ok_r2=0;
    total_n=0;
    for i=1:D.N
        tt=D.ms(:,i);
        spline_tt=tt(1):options.splinestep_ms:tt(end);
        n_ok_r2=0;
        yaw_spline=nan(numel(spline_tt),size(D.yaw{i},2));
        for t=1:size(D.yaw{i},2)
            yy=D.yaw{i}(:,t);
            perfectspline=csaps(tt,yy,1,spline_tt);
            splineyaw=csaps(tt,yy,options.splineroughness,spline_tt);
            % calculate the R2 of this spline so we can filter out really
            % wild traces that can't be captured with the roughness
            % provided like they should (ball shouldn't be jumping around
            % but change course more or less smoothly)
            R2=1-sum((perfectspline-splineyaw).^2)/sum((perfectspline-mean(perfectspline)).^2);
            if R2>options.min_spline_r2
                n_ok_r2=n_ok_r2+1;
                yaw_spline(:,n_ok_r2)=splineyaw(:);
            end
            total_n=total_n+1;
        end
        D.yaw{i}=yaw_spline(:,1:n_ok_r2);
        total_n_ok_r2=total_n_ok_r2+n_ok_r2;
    end
    toc(ttt);
    %
    % replace the time axes with the (coarser) spline time axis), assuming
    % it's the same for all conditions ...
    D.ms=repmat(spline_tt(:),1,D.N);
    %
    fprintf('[%s] removed %d/%d traces with R2<%f at splineroughness %f...\n',mfilename,total_n-total_n_ok_r2,total_n,options.min_spline_r2,options.splineroughness);
    %
    ttt=tic;
    fprintf('[%s] calculate mean yaw with sems ... ',mfilename);
    [D.yaw_mean,D.yaw_sem]=deal(nan(size(D.ms)));
    for i=1:D.N
        D.yaw_mean(:,i)=nanmean(D.yaw{i},2);
        D.yaw_sem(:,i)=nanstd(D.yaw{i},[],2)/realsqrt(size(D.yaw{i},2));
        time_ok = D.ms(:,i)>min(options.ttest_interval*1000) & D.ms(:,i)<max(options.ttest_interval*1000);
        D.yaw_interval{i}=mean(D.yaw{i}(time_ok,:),1);
        if any(isnan(D.yaw_interval{i}))
            keyboard
        end
    end
    toc(ttt);
    warning('on','SPLINES:CHCKXYWP:NaNs')
end



function plot_olrs(h,D,options,tail)
    D=dpxdSubset(D,D.mode~='u'); % just in case the unlimited data is still in there
    if numel(unique(D.mode))>1
        error('more than 1 mode');
    end
    if unique(D.mode)=='p'
        line='-';
    elseif unique(D.mode)=='r'
        line='--';
    else
        error('unknown mode');
    end
    [L,R]=dpxdSubset(D,D.dps<0); % left , right
    t=L.ms(:,1)/1000;
    if numel(unique(L.mouse))>1
        err=std(L.yaw_mean,[],2)/sqrt(size(L.yaw_mean,2));
    else
        err=L.yaw_sem;
    end
    revphi2019.jdPlotBounded('axes',h','x',t,'y',mean(L.yaw_mean,2),'eu',err,'ed',err,'Color','r','LineStyle',line);
    hold on
    if numel(unique(R.mouse))>1
        err=std(R.yaw_mean,[],2)/sqrt(size(R.yaw_mean,2));
    else
        err=R.yaw_sem;
    end
    revphi2019.jdPlotBounded('axes',h','x',t,'y',mean(R.yaw_mean,2),'eu',err,'ed',err,'Color','b','LineStyle',line);
    set(h,'Xlim',[min(t) max(t)]);
    cpsRefLine('-','k--');
    
    if numel(unique(R.mouse))>1
        left_vals=mean(L.yaw_mean(t>min(options.ttest_interval)&t<max(options.ttest_interval),:));
        right_vals=mean(R.yaw_mean(t>min(options.ttest_interval)&t<max(options.ttest_interval),:));
        [~,p,~,stats]=ttest(left_vals(:),right_vals(:),'tail',tail);
        titstr=sprintf('paired t-test (one-sided):\nt(%d)=%.2f, p=%.3f',stats.df,stats.tstat,p);
    else
        left_vals=L.yaw_interval{1};
        right_vals=R.yaw_interval{1};
        [~,p,~,stats]=ttest2(left_vals(:),right_vals(:),'tail',tail);
        titstr=sprintf('two-sample t-test, (one-sided):\nt(%d)=%.2f, p=%.3f',stats.df,stats.tstat,p);
    end
    if p<0.05 && p>0.01
        text(mean(options.ttest_interval),0,'*','FontSize',15,'HorizontalAlignment','center');
    elseif p<0.01 && p>0.001
        text(mean(options.ttest_interval),0,'**','FontSize',15,'HorizontalAlignment','center');
    elseif p<0.001
        text(mean(options.ttest_interval),0,'***','FontSize',15,'HorizontalAlignment','center');
    end
    title(titstr);
end

function [line_h,olr_permouse_12sec]=plot_right_minus_left(h,D,options)
    
    % olr_permouse_12sec output added 20191116 to make frame duration vs OLR
    % plot
    
    D=dpxdSubset(D,D.mode~='u'); % just in case the unlimited data is still in there
    [PHI,IHP]=dpxdSubset(D,D.mode=='p');
    for m=1:2
        if m==1
            D=PHI;
            line='-';
            field='phi';
        else
            D=IHP;
            line='--';
            field='revphi';
        end
        D=dpxdSplit(D,'mouse');
        right_minus_left_mean=[];
        for i=1:numel(D)
            [L,R]=dpxdSubset(D{i},D{i}.dps<0);
            right_minus_left_mean(:,end+1)=R.yaw_mean(:)-L.yaw_mean(:); %#ok<AGROW>
            if numel(D)==1 % only one mouse selected
                % calculate the pooled SEM
                R_yaw_n=size(R.yaw{1},2);
                L_yaw_n=size(L.yaw{1},2);
                R_yaw_var=(R.yaw_sem*sqrt(R_yaw_n)).^2;
                L_yaw_var=(L.yaw_sem*sqrt(L_yaw_n)).^2;
                right_minus_left_sem = sqrt(R_yaw_var+L_yaw_var)./sqrt(R_yaw_n+L_yaw_n);
            else
                olr_permouse_12sec.(field)(i)=mean(right_minus_left_mean(R.ms>=1000&R.ms<2000,end));
                olr_permouse_12sec.(field)(i) =olr_permouse_12sec.(field)(i)/2; % divide by 2 to get mean OLR
            end
        end
        t=L.ms(:,1)/1000;
        y=mean(right_minus_left_mean,2)/2; % divide by 2, take mean of curves plotted as if motion is always to the right
        if numel(D)>1
            e=std(right_minus_left_mean,[],2)/sqrt(size(right_minus_left_mean,2));
        else % only one mouse selected
            e=right_minus_left_sem;
        end
        line_h(m)=revphi2019.jdPlotBounded('axes',h','x',t,'y',y,'eu',e,'ed',e,'Color','m','LineStyle',line);
        hold on
        set(h,'Xlim',[min(t) max(t)]);
        cpsRefLine('-','k--');
        % store mean YAW over interval of interest for t-test
        if m==1
            phi_vals=mean(right_minus_left_mean(t>min(options.ttest_interval)&t<max(options.ttest_interval),:));
        elseif m==2
            revphi_vals=mean(right_minus_left_mean(t>min(options.ttest_interval)&t<max(options.ttest_interval),:));
        else
            error('??');
        end
    end
    [~,p,~,stats]=ttest(phi_vals(:),revphi_vals(:),'tail','right');
    titstr=sprintf('two-sample t-test, (one-sided):\nt(%d)=%.2f, p=%.3f',stats.df,stats.tstat,p);
    if p<0.05 && p>0.01
        text(mean(options.ttest_interval),0,'*','FontSize',15,'HorizontalAlignment','center');
    elseif p<0.01 && p>0.001
        text(mean(options.ttest_interval),0,'**','FontSize',15,'HorizontalAlignment','center');
    elseif p<0.001
        text(mean(options.ttest_interval),0,'***','FontSize',15,'HorizontalAlignment','center');
    end
    title(titstr);
end


function plot_diff_phi_contrast_inverted_revphi_contrast(D,options,row)
    % plot reversal index, calculate them first
    D=dpxdMerge(D);
    D=dpxdSplit(D,{'mouse','ff'});
    for i=1:numel(D)
        [L,R]=dpxdSubset(D{i},D{i}.dps<0); % left , right;
        [Lphi,Lrevphi]=dpxdSubset(L,L.mode=='p');
        [Rphi,Rrevphi]=dpxdSubset(R,R.mode=='p');
        flipped=(Lphi.yaw_mean>Rphi.yaw_mean & Lrevphi.yaw_mean<Rrevphi.yaw_mean) | (Lphi.yaw_mean<Rphi.yaw_mean & Lrevphi.yaw_mean>Rrevphi.yaw_mean);
        amount=abs(Lphi.yaw_mean-Rphi.yaw_mean).*abs(Lrevphi.yaw_mean-Rrevphi.yaw_mean);
        reversal_idx=  sqrt(  amount  );
        reversal_idx(~flipped)=reversal_idx(~flipped)*-1; % positive means flipped
        ylabel_str='Reversal (deg/s)';
        D{i}.reversal=repmat(reversal_idx,1,D{i}.N);
    end
    D=dpxdMerge(D);
    if options.freezeflips
        D=dpxdSplit(D,'ff');
    else
        K{1}=D;
        D=K;
        clear K
    end
    
    n_freeze_flips=numel(D);
    for i=1:n_freeze_flips
        reverseaxes(i)=subplot(options.nrows,n_freeze_flips,i+(row-1)*n_freeze_flips); %#ok<AGROW>
        t=D{i}.ms(:,1);
        y=nanmean(D{i}.reversal,2);
        err=nanstd(D{i}.reversal,[],2)/sqrt(size(D{i}.reversal,2)/4);
        revphi2019.jdPlotBounded('axes',reverseaxes(i),'x',t,'y',y,'eu',err,'ed',err,'Color','k','LineStyle','-');
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel(ylabel_str);
        end
        set(reverseaxes(i),'Xlim',[min(t) max(t)]);
        cpsRefLine(reverseaxes(i),'-','k--');
    end
    cpsUnifyAxes(reverseaxes);
end


function plot_mean_olr_over_freezeframes(olr_per_mouse,n_freeze_flips,D)
    
    % olr_per_mouse is 2nd output of plot_left_vs_right
    fig=cpsFindFig('plot_mean_olr_over_freezeframes');
  
    ff=nan(1,n_freeze_flips);
    phi_mean=nan(1,n_freeze_flips);
    revphi_mean=nan(1,n_freeze_flips);
    phi_sem=nan(1,n_freeze_flips);
    revphi_sem=nan(1,n_freeze_flips);
    for i=1:n_freeze_flips
        ff(i)=unique(D{i}.ff);
        phi_mean(i)= mean(olr_per_mouse(i).phi);
        phi_sem(i)=std(olr_per_mouse(i).phi)/sqrt(numel(olr_per_mouse(i).phi));
        revphi_mean(i)=mean(olr_per_mouse(i).revphi);
        revphi_sem(i)=std(olr_per_mouse(i).revphi)/sqrt(numel(olr_per_mouse(i).revphi));
    end
    set(fig,'defaultLegendAutoUpdate','off');
    clf(fig);
    h(1)=errorbar(ff,phi_mean,phi_sem,'ok-','LineWidth',1,'MarkerFaceColor','k');
    hold on
    h(2)=errorbar(ff,revphi_mean,revphi_sem,'ok--','LineWidth',1,'MarkerFaceColor','w');
    xlabel('Frame duration (ms)','FontSize',14);
    ylabel('OLR (deg/s)','FontSize',14);
    legend(h,{'Phi motion','Reverse phi motion'},'FontSize',12);
    box off
    set(gca,'XLim',[0.5 n_freeze_flips+0.5]);
    set(gca,'XTickLabel',round((1:n_freeze_flips)*1000/60));
    set(gca,'YLim',[-20 20]);
    drawnow
    pause(1/10);
    cpsRefLine('-','k--')
end

