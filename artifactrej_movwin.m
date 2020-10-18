function [trials_clean] = artifactrej(trials,allinclchans, artf_th, perc_chan, winms, stepms, fs)
%Artifact rejection of trial for all channels if trial bad in one channel,
%given 100 uV peak to peak
% New version of script 06/2020 that rejects trial if more than a certain
% percent of channels for that trial are above threshold.

triallist=1:size(trials,3);
%badchans=zeros(1,size(allinclchans,2));
%badtrials = zeros(1, size(trials,3));


winpnts  = floor(winms*fs/1000);
stepnts  = floor(stepms*fs/1000);
p1 = 1;
p2 = size(trials,2);

% c = squeeze(max(trials(:,:,:),[],2));
% d = squeeze(min(trials(:,:,:),[],2));
% peaktopeak = c - d;

for t = 1:size(trials,3)
    clear rej_win
    for j = p1:stepnts:p2-(winpnts-1)
        cur = trials(:, j:j+winpnts-1, t);
        c = squeeze(max(cur,[],2));
        d = squeeze(min(cur,[],2));
        peaktopeak_win = abs(c - d);
        rej_win = find(peaktopeak_win>artf_th);
        if length(rej_win) > floor(size(trials,1)*perc_chan)
            triallist(t)=0;
            break
        end
    end
end

trialsvalid=triallist>0;
trials_clean=trials(allinclchans,:,trialsvalid);
end

% function [trials_clean] = artifactrej(trials,allinclchans,artf_th, perc_chan)
% %Artifact rejection of trial for all channels if trial bad in one channel,
% %given 100 uV peak to peak
% 
% % trials=std_ERP_alltr_windowed;
% % allinclchans=Chans;
% % artf_th=100;
% 
% triallist=1:size(trials,3);
% badchans=zeros(1,size(allinclchans,2));
% for j=1:size(allinclchans,2)
%     chan=allinclchans(j);
%     peaktopeak=(squeeze(max(trials(chan,:,:),[],2))')-(squeeze(min(trials(chan,:,:),[],2))');
%     rej=find(peaktopeak>artf_th);
%     if ~isempty(rej)
%         badchans(j)=1;
%     end
%     triallist(rej)=0;
% end
% 
% trialsvalid=triallist>0;
% trials_clean=trials(allinclchans,:,trialsvalid);
% end