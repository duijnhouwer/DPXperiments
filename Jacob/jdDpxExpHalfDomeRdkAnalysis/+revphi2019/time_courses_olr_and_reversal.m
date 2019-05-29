function time_courses_olr_and_reversal
    
    options.pitch_range=[5 Inf];
    options.roll_range=[-Inf Inf];
    options.min_trials_per_condition=25;
    options.splineroughness=1/3;
    options.min_spline_r2=0;
    options.detrend_per_mouse=true;
    options.nrows=4;
    options.include_mice=1:9; % 0 and 1:9 mean all mice, 1 means mouse 1 etc, [1 2 4 5] means these mice
    options.freezeflips=[1 2 3 4 5 6 7]; % freezeflips to keep
    options.pool_freezeflips=false;
    options.reversal_measure='signed absolute product'; % tstat, linear, multiplicative, divisive, signrank, ttest
    
    load('barebonesreversephidata.mat','DPXD'); % created with revphi2019.datafiles_to_barebones_dpxd.m
    D=DPXD;
    clear DPXD;
    
    % remove the unlimited lifetime data
    D=dpxdSubset(D,D.mode~='u');
    
    % remove the freezeflips that we dont' want to analyse
    D=dpxdSubset(D,ismember(D.ff,options.freezeflips));
    if options.pool_freezeflips
        D.freezeflip=666;
    end
    
    D=subset_trials_by_pitch_and_roll(D,options);
    [D,options.include_mice]=remove_mice_conditions_with_too_few_n(D,options.min_trials_per_condition); % good to remove outliers now so no spline fitting required, save time
    D=convert_yaw_per_trial_to_yaw_per_condition(D,options);
    [D,options.include_mice]=remove_mice_conditions_with_too_few_n(D,options.min_trials_per_condition); % remove more outliers (if any)
        
    % Detrend the yaw data by removing the mean over all conditions for
    % each mouse
    if options.detrend_per_mouse
        D=dpxdSplit(D,{'mouse'});
        for i=1:numel(D)
            D{i}.yaw_mouse_trend=repmat(mean(D{i}.yaw_mean,2),1,D{i}.N);
            D{i}.yaw_mean=D{i}.yaw_mean-D{i}.yaw_mouse_trend;
        end
        D=dpxdMerge(D);
    end
    
    if options.include_mice
        D=dpxdSubset(D,ismember(D.mouse,options.include_mice)); % for debugging
    end
    
    
    % plot mean OLRs
    cpsFindFig(mfilename);
    clf
    D=dpxdSplit(D,'ff');
    n_freeze_flips=numel(D);
    for i=1:n_freeze_flips
        [PHI,REV]=dpxdSubset(D{i},D{i}.mode=='p');
        h=subplot(options.nrows,n_freeze_flips,i);
        plot_olrs(h,PHI);
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel('OLR (deg/s)');
        end
        h=subplot(options.nrows,n_freeze_flips,i+n_freeze_flips);
        plot_olrs(h,REV);
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel('OLR (deg/s)');
        end
        h=subplot(options.nrows,n_freeze_flips,i+2*n_freeze_flips);
        line_handles=plot_right_minus_left(h,D{i});
        if i==1
            xlabel('Time since motion onset (s)');
            ylabel([glyph('Delta') 'OLR (deg/s)']);
            legend(line_handles,{'Phi','Reverse Phi'});
        end
    end
    cpsUnifyAxes
    
    plot_diff_phi_contrast_inverted_revphi_contrast(D,options,4)
    suptitle(['Mouse' sprintf(' %d',options.include_mice)]);
    
    options
end

function D=subset_trials_by_pitch_and_roll(D,options)
    % pitch = forward running
    for i=1:D.N
        ok_pitch=median(D.pitch{i},1)>options.pitch_range(1) & median(D.pitch{i},1)<options.pitch_range(2);
        ok_roll=median(D.roll{i},1)>options.roll_range(1) & median(D.roll{i},1)<options.roll_range(2);
        ok=ok_pitch & ok_roll;
        D.yaw{i}=D.yaw{i}(:,ok);
        fprintf('[%s] removed %s from condition %d\n',mfilename,jdProp(~ok,'%'),i);
    end
end

function [D,good_mice]=remove_mice_conditions_with_too_few_n(D,min_trials_per_condition)
   D=dpxdSplit(D,'mouse');
   remove=false(size(D));
   good_mice=[];
   for i=1:numel(D)
       if any(cellfun(@(x)(size(x,2)),D{i}.yaw)<min_trials_per_condition)
           remove(i)=true;
           fprintf('[%s] Discarding mouse %d because one or more conditions had fewer that %d remaining trials\n',mfilename,unique(D{i}.mouse),min_trials_per_condition);
       else
           good_mice(end+1)=unique(D{i}.mouse); %#ok<AGROW>
       end
   end
   D(remove)=[];
   D=dpxdMerge(D);
end


function D=convert_yaw_per_trial_to_yaw_per_condition(D,options)
    ttt=tic;
    fprintf('[%s] removing pre-stim means... ',mfilename);
    for i=1:D.N
        mean_yaw_prestim=mean(D.yaw{i}(D.ms(:,i)>=-0.5 & D.ms(:,i)<=0,:),1);
        D.yaw{i}=D.yaw{i}-mean_yaw_prestim;
    end
    toc(ttt);
    %
    ttt=tic;
    fprintf('[%s] replacing yaw traces with splines...',mfilename);
    for i=1:D.N
        xx=D.ms(:,i);
        for t=1:size(D.yaw{i},2)
            thisyaw=D.yaw{i}(:,t);
           % if any(isnan(thisyaw)) %|| any(isnan(thispitch)) || any(isnan(thisroll))
           %     warning off
           % end
           dbstop if warning
            splineyaw=csaps(xx,thisyaw,options.splineroughness,xx);
            warning on
            % calculate the R2 of this spline so we can filter out really
            % wild traces that can't be captured with the roughness
            % provided like they should (ball shouldn't be jumping around
            % but change course more or less smoothly)
            R2=1-sum((thisyaw-splineyaw).^2)/sum((thisyaw-mean(thisyaw)).^2);
            if R2>options.min_spline_r2
                D.yaw{i}(:,t)=splineyaw;
            else
                D.yaw{i}(:,t)=Inf; % mark for removal
            end
        end
        
    end
    toc(ttt);
    %
    ttt=tic;
    fprintf('[%s] removing traces with R2<%f at splineroughness %f...',mfilename,options.min_spline_r2,options.splineroughness);
    n_removed=0;
    n_total=0;
    for i=1:D.N
        n_removed=n_removed+sum(mean(D.yaw{i})==Inf);
        n_total=n_total+size(D.yaw{i},2);
        D.yaw{i}(:,mean(D.yaw{i})==Inf)=[];
    end
    fprintf(' removed %d/%d traces. ',n_removed,n_total);
    toc(ttt);
    %
    ttt=tic;
    fprintf('[%s] calculate mean yaw with sems ... ',mfilename);
    [D.yaw_mean,D.yaw_sem]=deal(nan(size(D.ms)));
    for i=1:D.N
        D.yaw_mean(:,i)=nanmean(D.yaw{i},2);
        D.yaw_sem(:,i)=nanstd(D.yaw{i},[],2)./size(D.yaw{i},2);
    end
    toc(ttt);
end



function plot_olrs(h,D)
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
    t=L.ms(:,1);
    err=std(L.yaw_mean,[],2)/sqrt(size(L.yaw_mean,2));
    revphi2019.jdPlotBounded('axes',h','x',t,'y',mean(L.yaw_mean,2),'eu',err,'ed',err,'Color','r','LineStyle',line);
    hold on
    err=std(R.yaw_mean,[],2)/sqrt(size(R.yaw_mean,2));
    revphi2019.jdPlotBounded('axes',h','x',t,'y',mean(R.yaw_mean,2),'eu',err,'ed',err,'Color','b','LineStyle',line);
    set(h,'Xlim',[min(t) max(t)]);
    cpsRefLine('-','k--');
end

function line_h=plot_right_minus_left(h,D)
    
    [PHI,IHP]=dpxdSubset(D,D.mode=='p');
    for m=1:2
        if m==1
            D=PHI;
            line='-';
        else
            D=IHP;
            line='--';
        end
        D=dpxdSplit(D,'mouse');
        right_minus_left=[];
        for i=1:numel(D)
            [L,R]=dpxdSubset(D{i},D{i}.dps<0);
            right_minus_left(:,end+1)=R.yaw_mean(:)-L.yaw_mean(:); %#ok<AGROW>
        end
        t=L.ms(:,1);
        y=mean(right_minus_left,2);
        e=std(right_minus_left,[],2)./size(right_minus_left,2);
        line_h(m)=revphi2019.jdPlotBounded('axes',h','x',t,'y',y,'eu',e,'ed',e,'Color','m','LineStyle',line);
        hold on
        set(h,'Xlim',[min(t) max(t)]);
        cpsRefLine('-','k--');
    end
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

