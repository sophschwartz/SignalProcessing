function [lower_cluster, upper_cluster] = temporalwindow(data_folder,...
    subjectlists, time_start_ind, time_end_ind)

%%%%%%%%%
% Temporal component identification processing
% schwart2@bu.edu
%
% Inputs:
% subjectlists is a cell containing cells with subject string names
% example: 
%    subjlistA={'sub01', 'sub02', 'sub03'};
%    subjlistB={'sub04', 'sub05'};
%    subjlists={subjlistA, subjlistB};
% files inside data folder have:
% mERP_std and mERP_dev: 1x3 cell. Each cell stores a 3D array with 1 x
% time x trials. Cells are for Oddball, Integrated, and Segregated
% conditions, respectively.
% 
% time_start_ind, time_end_ind = Starting time point index, ending time point index for time points to
% consider window in
%%%%%%%%%

for list = 1:length(subjectlists)
    subjlist = subjectlists{list}
    nsub=length(sublist);
    for isub=1:nsub
        fname=subjlist{isub};
        if isempty(fname)
            error(['No file found for', subj]);
        else
            load([data_folder fname '_IntResults_Fz_03062018_final.mat']);
            for icond=1:3
                 std_trialm=mean(squeeze(mERP_std{icond}),2);
                 dev_trialm=mean(squeeze(mERP_dev{icond}),2);
                 fERP_std(:,isub,icond)=std_trialm;
                 fERP_dev(:,isub,icond)=dev_trialm;
                 fMMN(:,isub,icond)=std_trialm-dev_trialm;

            end
        end
    end

    % calculate t values and p values for each time point from 30 ms to 350 ms
    % (should get 320 points of interest) between the 150 or so deviant trials
    % and a random selection of 150 standard trials
    icond=1;
    oddball_dev=squeeze(fERP_dev(176:401,:,icond));
    oddball_mmn=squeeze(fMMN(176:401,:,icond));
    oddball_std=squeeze(fERP_std(176:401,:,icond));

    for i=1:225
        [h,p,ci,stats]=ttest(oddball_dev(i,:),oddball_std(i,:));
        tval(i)=stats.tstat;
        pval(i)=p;
    end

    % plot the t val distrubution (201 points) and determine a 5% threshold
    % point to consider a cluster as a 'cluster'
    % characterize that t value threshold
    
    plot(t(151:501), pval);hold on; plot([51 200],[0.05 0.05]);
    clear tclusterindex;
    clusterpieces=[];

    if ~isempty(find(pval)<0.05)
        tclusterindex=find(abs(tval)>2.99);

        firstindex=tclusterindex(1);
        lastindex=tclusterindex(end);
        jumpend=[]; jumpbeginning=[];
        for i=1:length(tclusterindex)-1
            priorindex=tclusterindex(i);
            nextindex=tclusterindex(i+1);
            if nextindex ~= priorindex+1
                jumpend=[jumpend priorindex];
                jumpbeginning=[jumpbeginning nextindex];
            end
        end
        jumpend=[jumpend tclusterindex(end)]

        % Sum each cluster of t values so you have multiple cluster sums
        clear tvalue;
        for i=1:length(jumpend)
            if i==1
                tvalue(i)=sum(tval(firstindex:jumpend(i)));
                clusterpieces{i}=firstindex:jumpend(i)
            else
                tvalue(i)=sum(tval(jumpbeginning(i-1):jumpend(i)));
                clusterpieces{i}=jumpbeginning(i-1):jumpend(i)
            end
        end
    else
        tvalue=0;
    end

    % Mock randomization for multiple comparisons
    tvaluemockfull=[]

    for i=1:1000
        allsamples=vertcat(oddball_std,oddball_dev);
        chosemocktrials=randperm(size(allsamples,1));
        mockdev=allsamples(chosemocktrials(1:(size(allsamples,1)/2)),:);
        mockzero=allsamples(chosemocktrials((size(allsamples,1)/2)+1:end),:);

        for k=1:225
            [h,p,ci,stats]=ttest(mockdev(k,:),mockzero(k,:));
            tvalmock(k)=stats.tstat;
            pvalmock(k)=p;
        end

        % Use the previously determined t value threshold to determine where
        % significant clusters are
        clear tclusterindex;
        if ~isempty(find(pvaluemock<0.01))
            tclusterindex=find(abs(tvalmock)>1);
            firstindex=tclusterindex(1);
            lastindex=tclusterindex(end);
            jumpend=[]; jumpbeginning=[];
            for d=1:length(tclusterindex)-1
                priorindex=tclusterindex(d);
                nextindex=tclusterindex(d+1);
                if nextindex ~= priorindex+1
                    jumpend=[jumpend priorindex];
                    jumpbeginning=[jumpbeginning nextindex];
                end
            end
            jumpend=[jumpend tclusterindex(end)];

        % Sum each cluster of t values so you have multiple cluster sums
            for r=1:length(jumpend)
                if r==1
                    tvaluemock(r)=sum(tvalmock(firstindex:jumpend(r)));
                else
                    tvaluemock(r)=sum(tvalmock(jumpbeginning(r-1):jumpend(r)));
                end
            end
        else
            tvaluemock=0;
        end
        tvaluemockfull=[tvaluemockfull tvaluemock];
    end

    % Identify clusters that are over 95% Confidence Interval of
    % distribution curve determined with mock sample
    
    lowerbound=mean(tvaluemockfull)-1.96*(std(tvaluemockfull));
    upperbound=mean(tvaluemockfull)+1.96*(std(tvaluemockfull));
    significanttvals_lower=find(tvalue<=lowerbound);
    significantclusterpieces_lower=clusterpieces(significanttvals_lower)
    significanttvals_upper=find(tvalue>=upperbound);
    significantclusterpieces_upper=clusterpieces(significanttvals_upper)
    
    %Output based on time sample
    tsub=t(176:351);
    lower_cluster = tsub(significantclusterpieces_lower{1});
    upper_cluster = tsub(significantclusterpieces_upper{1});
end    
end 