%%%%%%%%%
% Preprocessing of EEG data (128 channel EGI) from Oddball, Segregated, and
% Integrated Mismatch Negativity paradigms
% @schwart2@bu.edu
% 
% For each participant:
% Load in either .sb or .mat files
% Filtering: 0.3 - 35 Hz
% Downsampling
% Epoch data
% Artifact rejection
% Rereference to mastoids
% ICA for additional artifact rejection of eye blinks
% 100 ms pre-stimulus baselining rather than mean trial baselining
% Random permutation (without replacement) to select equal number of trials across conditions
%
%%%%%%%%%
%% List subjects by age group and load in.
% Deidentified: Subject IDs and folder names removed for sample code

clear all;
EGI_data_folder = '';
EGI_Results_folder = '';

subjlistA={''};

subjlistB={''};

subjlistC={''};

subjlistD={''};

subjlistE={''};

subjlists={subjlistA, subjlistB, subjlistC, subjlistD, subjlistE};


%% Parameter Settings
t = -100:451;
t_downsamp = downsample(t,4);
epoch_dur = 552;
artf_th = 100;
fs = 1000;
norder = fs;
cf1 = 0.3;
cf2 = 35;
Ny = fs/2;
Wn = [cf1 cf2]/(Ny);
h_bp = fir1(norder,Wn);
h_hp = fir1(norder*6,1/Ny,'high');
a = 1;

% cf1 = 1;
% cf2 = 59;
% Ny = fs/2;
% Wn = [cf1 cf2]/(Ny);
% h_bp_broad = fir1(norder,Wn);

% 10-20 system electrodes + left/right mastoid clusters (100/107, 56/57)
chan_list = [9,11,22,24,33,36,45,52,57,58,62,70,83,92,96,100,104,108,122,124];

intensity_trigger = {
    'DIN2', 'DIN4';
    'DIN6', 'DIN8';
    'DI10', 'DI12';
    'DI14', 'DI16';
    'DI18', 'DI20';
    'DI22', 'DI24';
    };
conds = {
    'Oddball Experiment'; % Condition 1. loud deviant, soft standard (select loud only later: icond*2 - 1)
    'Oddball Control'; % Condition 2. soft deviant, loud standard. (select loud only later: icond*2)
    'Integrated Experiment';
    'Integrated Control';
    'Segregated Experiment';
    'Segregated Control';
    };
stream_trigger = 'D150';
ncond = 6; nblk = 6;

load('ica_grp_input.mat');
load('file_proc_info_sample.mat');

for igroup=1:5
    %% Load in file
    clear subjlist
    subjlist=subjlists{igroup};
    
    for isub=1:length(subjlist)
        clearvars -except subjlist igroup isub EGI_data_folder conds intensity_trigger stream_trigger...
            chan_list h_hp h_bp a epoch_dur t t_downsamp artf_th subjlists ncond...
            EGI_Results_folder EGI_data_folder file_proc_info grp_proc_info_in h_bp_broad
        fname=subjlist{isub};
        isub
        if isempty(fname)
            error(['No file found for', subj]);
        else
            load([EGI_data_folder fname '.mat']);
            varlist=whos;
            
            nblk = 0;
            for ivar = 1:length(varlist)
                if ~isempty(strfind(varlist(ivar).name,fname))
                    nblk = nblk + 1;
                    order(nblk) = str2num(varlist(ivar).name(length(fname)+1:end));
                    eval(['blkdata{nblk} = ' varlist(ivar).name ';']);
                end
            end
            [y, I] = sort(order);
            blkdata = blkdata(I);
            
            
            %% Organization of all files
            for iblk = 1:nblk
                oddballstreamevents = cat(2,Oddball,Integrated,Segregated);
                index=find(cell2mat(oddballstreamevents(3,:))==iblk);
                alltriggers = oddballstreamevents(1,index);
                triggerlist = unique(alltriggers);
                len = [];
                for i = 1:length(triggerlist)
                    len(i) = length(strfind(cell2mat(alltriggers),triggerlist{i}));
                end
                [Y I] = sort(len);
                maintrigger{iblk} = cell2mat(triggerlist([I(end),I(end-1)]));
            end
            cond_trigger = sort(unique(maintrigger));
            if length(cond_trigger) ~= ncond
                error('number of conditions does not match the trigger list!');
            end
            [Y, I] = sort(mat2cell(cell2mat(intensity_trigger),[1 1 1 1 1 1],[8]));
            [Y, I] = sort(I);
            cond_trigger = cond_trigger(I);
            data = cell(6,1); 
            %data=cell(7,1); %switching from 6 to seven temporarily for extra block
            Oddball_bk = Oddball;
            Integrated_bk = Integrated;
            Segregated_bk = Segregated;
            for i = 1:length(cond_trigger)
                blklist{i} = find(~cellfun(@isempty,strfind(maintrigger,cond_trigger{i})));
                pre_blk_total_len_Odd = 0;
                pre_blk_total_len_Int = 0;
                pre_blk_total_len_Seg = 0;
                for j = 1:length(blklist{i})
                    data{i} = cat(2,data{i},blkdata{blklist{i}(j)});
                    [Oddball{3,(cell2mat(Oddball_bk(3,:))==blklist{i}(j))}] = deal(i);
                    [Integrated{3,(cell2mat(Integrated_bk(3,:))==blklist{i}(j))}] = deal(i);
                    [Segregated{3,(cell2mat(Segregated_bk(3,:))==blklist{i}(j))}] = deal(i);
                    if j > 1
                        ind_Oddball = cell2mat(Oddball_bk(3,:))==blklist{i}(j);
                        pre_blk_total_len_Odd = pre_blk_total_len_Odd + size(blkdata{blklist{i}(j-1)},2);
                        corrected_Oddballtiming = num2cell(cell2mat(Oddball_bk(2,(ind_Oddball)))+pre_blk_total_len_Odd);
                        [Oddball{2,ind_Oddball}] = deal(corrected_Oddballtiming{:});
                        
                        ind_Integrated = cell2mat(Integrated_bk(3,:))==blklist{i}(j);
                        pre_blk_total_len_Int = pre_blk_total_len_Int + size(blkdata{blklist{i}(j-1)},2);
                        corrected_Integratedtiming = num2cell(cell2mat(Integrated_bk(2,ind_Integrated))+pre_blk_total_len_Int);
                        [Integrated{2,ind_Integrated}] = deal(corrected_Integratedtiming{:});
                        
                        ind_Segregated = cell2mat(Segregated_bk(3,:))==blklist{i}(j);
                        pre_blk_total_len_Seg = pre_blk_total_len_Seg + size(blkdata{blklist{i}(j-1)},2);
                        corrected_Segregatedtiming = num2cell(cell2mat(Segregated_bk(2,ind_Segregated))+pre_blk_total_len_Seg);
                        [Segregated{2,ind_Segregated}] = deal(corrected_Segregatedtiming{:});
                    end
                end
            end
            [Oddball{4,:}] = deal(Oddball{2,:});
            [Integrated{4,:}] = deal(Integrated{2,:});
            [Segregated{4,:}] = deal(Segregated{2,:});
            
            %% Organization Complete. Now: Preprocessing: Artifact rejection, filtering, baselining.
            
            %% Set up epoching (not actually done here)
            for iblk = 1:ncond
                clearvars ind_cond evtlist evtime std_ind dev_ind std_time dev_time std_win dev_win...
                    TrialDataRelevant TrialDataChan_filt TrialData_filt std_alltr_allchan std_windowed_data...
                    dev_alltr_allchan dev_windowed_data std_bl dev_bl permuted_std_win data_std_resamp...
                    permuted_dev_win data_dev_resamp std_trials_clean dev_trials_clean...
                    std_trials_clean_m dev_trials_clean_m std_bl dev_bl
                switch iblk
                    case {1,2}
                        evtlist = Oddball;
                    case {3,4}
                        evtlist = Integrated;
                    case {5,6}
                        evtlist = Segregated;
                end
                
                ind_cond = cell2mat(evtlist(3,:));
                evtime = cell2mat(evtlist(2,:));
                
                evlabel = evtlist(1,:);
                std_ind = strcmp(intensity_trigger{iblk,1},evlabel);
                dev_ind = strcmp(intensity_trigger{iblk,2},evlabel);
                
                std_time = evtime(std_ind&(ind_cond==iblk));
                std_time(std_time>size(data{iblk},2)-epoch_dur) = [];
                dev_time = evtime(dev_ind&(ind_cond==iblk));
                dev_time(dev_time>size(data{iblk},2)-epoch_dur) = [];
                
                % Of note, the experiment was changed half way so some
                % participants had less trials for conditions that mattered
                % less (soft standards)
                std_win = cat(1,std_time,ones(551,size(std_time,2)));
                std_win = cumsum(std_win)-100;
                dev_win = cat(1,dev_time,ones(551,size(dev_time,2)));
                dev_win = cumsum(dev_win)-100;
                
                
                %% If testing out AEP window...
                stmtime = cell2mat(ECI_TCPIP_55513(2,strcmp(ECI_TCPIP_55513(1,:),'stm+')));
                ERP_win = cat(1,stmtime,ones(699,size(stmtime,2)));
                ERP_win = cumsum(ERP_win)-100;
        
                %% To optimize artifact rejection performance, a spatially
                % distributed subset of channels were selected (10-20 system + mastoid region).
                %TrialDataRelevant = data{iblk}(chan_list,:);
                TrialDataRelevant = data{iblk}(chan_list,:);
                
                % Filter data (0.1 - 35 Hz)
                % Do not downsample before filtering, unless performing
                % good anti-aliasing procedures.
                TrialData_filt=zeros(size(TrialDataRelevant));
                for i=1:size(TrialDataRelevant,1)
                    TrialData_filt(i,:)=filtfilt(h_bp,a,TrialDataRelevant(i,:));
                end
                
                % wICA and regular ICA - potentially do this on full data hz and then apply to 0.3-35 hz data?
                file_proc_info.beapp_filt_max_freq = 250;
                grp_proc_info_in.eeg = TrialDataRelevant;
                % grp_proc_info_in. =;
                
                file_proc_info.beapp_fname{1} = [fname '.mat'];
                grp_proc_info_in.beapp_fname_all{1} = [fname '.mat'];
                
                grp_proc_info_in_new = specialized_batch_beapp_ica(grp_proc_info_in, file_proc_in, TrialData_filt)
                
                % Segment data into epochs
                std_alltr_allchan=TrialData_filt(:,std_win);
                std_windowed_data=reshape(std_alltr_allchan,[size(TrialData_filt,1),size(std_win,1),size(std_win,2)]);
                
                dev_alltr_allchan=TrialData_filt(:,dev_win);
                dev_windowed_data=reshape(dev_alltr_allchan,[size(TrialData_filt,1),size(dev_win,1),size(dev_win,2)]);
                
                % Pre-stimulus baseline
                dev_bl=dev_windowed_data - repmat(mean(dev_windowed_data(:,1:100,:),2),[1,size(dev_win,1),1]); %size(dev_win,1) should be about 549
                std_bl=std_windowed_data - repmat(mean(std_windowed_data(:,1:100,:),2),[1,size(dev_win,1),1]);
                
                % Resample data to 250 Hz AFTER filtering to avoid aliasing higher frequencies
                srate_in = 1000; srate_out = 250;
                permuted_std_win = permute(std_bl, [2 3 1]);
                data_std_resamp(:,:,:)=permute((downsample(permuted_std_win(:,:,:),srate_in/srate_out,1)), [3 1 2]);
                
                permuted_dev_win = permute(dev_bl, [2 3 1]);
                data_dev_resamp(:,:,:)=permute((downsample(permuted_dev_win(:,:,:),srate_in/srate_out,1)), [3 1 2]);
                
                % Artifact rejection of 10-20 system with moving window. Trials with more than
                % 15% bad channels are rejected. Flag any participants with fewer than 80 trials per condition.
                % use artif_rej script if don't use ica.
                % Maybe use MARA ML Artifact rejection and ICA to say how eyeblinks were considered.
                perc_chan = 0.15; %85% of channels must be usable otherwise reject the trial
                allinclchans=1:size(dev_windowed_data,1);
                winms = 100;
                stepms = 50;
                std_trials_clean=artifactrej_movwin(data_std_resamp, allinclchans, artf_th, perc_chan, winms, stepms, srate_out);
                dev_trials_clean=artifactrej_movwin(data_dev_resamp, allinclchans, artf_th, perc_chan, winms, stepms, srate_out);
                
                % flag if trials less than we need, don't bother analyzing
                exit = 0;
                if (size(dev_trials_clean,3) < 80 && (iblk == 1 || iblk == 3 || iblk == 5))
                    p = [fname ' ' num2str(isub) ' ' num2str(iblk) ' NOT ENOUGH DEV TRIALS' num2str(size(dev_trials_clean,3))]
                    exit = 1;
                    break
                end
                if (size(std_trials_clean,3) < 80 && (iblk == 2 || iblk == 4 || iblk == 6))
                    p = [fname ' ' num2str(isub) ' ' num2str(iblk) ' NOT ENOUGH STD TRIALS' num2str(size(std_trials_clean,3))]
                    exit = 1;
                    break
                end
                
                % Joint Probability artifact rejection
%                 chan_labels = {EEG_tmp.chanlocs(beapp_indx{curr_epoch}).labels};
%                 EEG_tmp = pop_select(EEG_tmp,'channel', chan_labels);
%                 EEG_tmp = pop_jointprob(EEG_tmp,1,[1:length(chan_labels)],3,3,grp_proc_info_in.beapp_happe_seg_rej_plotting_on,0,...
%                     grp_proc_info_in.beapp_happe_seg_rej_plotting_on,[],0);
%                 EEG_tmp = eeg_rejsuperpose(EEG_tmp, 1, 0, 1, 1, 1, 1, 1, 1);
                
                % Rereference data to mastoids
                %ref_TrialData = ref_TrialData - repmat(mean(ref_TrialData([56,107],:)),[128,1]);  % 56,107 are the channels closest to L/R mastoids
                lmast = 10; rmast = 17;
                std_trials_clean_m = std_trials_clean - repmat(mean(std_trials_clean([lmast,rmast],:,:)),[size(std_trials_clean,1), 1, 1]);
                dev_trials_clean_m = dev_trials_clean - repmat(mean(dev_trials_clean([lmast,rmast],:,:)),[size(dev_trials_clean,1), 1, 1]);
                
                % Prestimulus baseline again
                clear dev_bl; clear std_bl;
                dev_bl=dev_trials_clean_m - repmat(mean(dev_trials_clean_m(:,1:25,:),2),[1,size(dev_trials_clean_m,2),1]);
                std_bl=std_trials_clean_m - repmat(mean(std_trials_clean_m(:,1:25,:),2),[1,size(std_trials_clean_m,2),1]);
                
                % save data for all blocks
                std_trials_tostore{iblk} = std_bl;
                dev_trials_tostore{iblk} = dev_bl;
                
            end %block loop
        end %loading if statement (not a loop)
        
        if exit == 1 %if the last loop exited, don't bother saving data from this subjectr
            continue
        end
        %Randomly select standards for subject to be equivalent to
        %number of deviants (randperm)
        trialsize = [];
        for icond = 1:3
            trialsize = [trialsize size(dev_trials_tostore{icond*2 - 1},3)];
        end
        dev_ntrials = min(trialsize);
        for icond = 1:3
            stdperm=randperm(size(std_trials_tostore{icond*2},3),dev_ntrials);
            std_ERP{icond}=std_trials_tostore{icond*2}(:,:,stdperm);
            
            devperm=randperm(size(dev_trials_tostore{icond*2-1},3),dev_ntrials);
            dev_ERP{icond}=dev_trials_tostore{icond*2-1}(:,:,devperm);
            
        end
        
        %     for icond=1:3
        %         %totalmin_std=min(size(std_ERP{icond*2},2));
        %         stdperm=randperm(size(std_trials_tostore{icond*2},3),size(dev_trials_tostore{icond},3));
        %         stdperm_tosave{icond}=stdperm;
        %         std_ERP{icond}=std_trials_tostore{icond*2}(:,:,stdperm);
        %         dev_ERP{icond}=dev_trials_tostore{icond*2-1};
        %     end
        
        % save full scalp data for each subject
        save([EGI_Results_folder fname '_Preprocessing_date.mat'],'fname','t_downsamp','std_ERP','dev_ERP');
    end %file loop
end %subject group loop
