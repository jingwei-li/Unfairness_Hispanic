function [UPPS, UPPS_hdr, UPPS_colloquial] = ABCD_read_UPPS(subj_list, race, dohist, hist_dir, hist_fstem)

% [UPPS, UPPS_hdr, UPPS_colloquial] = ABCD_read_UPPS(subj_list, race, dohist, hist_dir, hist_fstem)
%
% Read and plot histograms of the mecessary measures from Modified UPPS-P for Children from PhenX.
%
% Inputs:
% - subj_list
%   List of subjects which passed fMRI prepreocessing quality control (full path). Default:
%   '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt'
%
% - race
%   A cell of strings. Each cell array corresponds to the ethnicity/race of one subject.
%   This cell can be obtained by `ABCD_read_race.m`.
% 
% - dohist
%   A 1/0 value determining whether the histograms are created or not. 
%   Default: 1, i.e. create plots.
%
% - hist_fname
%   Full path of output histogram filename.
%
% Outputs:
% - UPPS
%   A #subjects x 5 matrix. Each row corresponds to the behavioral scores of a subject.
%   The 5 columns correspond to Negative urgency, Positive urgency, Lack of planning, 
%   Lack of perseverance, Sensation seeking.
%
% - UPPS_hdr
%   A 1x5 cell. Headers of these 5 measures in ABCD csv file.
%  
% - UPPS_colloquial
%   A 1x5 cell. Colloquial names of these 5 measures.
%
% Example:
% [UPPS, UPPS_hdr, UPPS_colloquial] = ABCD_read_UPPS([], race, [], ...
%     '~/storage/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
% where "race" is obtained from
% race = ABCD_read_race([], [], ...
%     '~/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');
%
% Author: Jingwei Li

repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_dir, 'external')))
start_dir = pwd;

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

proj_dir = '/data/project/FairAI_Hispanic';
csv_dir = fullfile(proj_dir, 'data', 'inm7-superds', 'original', 'abcd', 'phenotype');
cd(csv_dir)
system('datalad get -n .')
system('git -C . config --local --add remote.datalad.annex-ignore true')
cd(fullfile(csv_dir, 'phenotype'))
system('datalad get -s inm7-storage abcd_mhy02.txt')

UPPS_csv = fullfile(csv_dir, 'phenotype', 'abcd_mhy02.txt');
UPPS_hdr = {'upps_y_ss_negative_urgency', 'upps_y_ss_positive_urgency', 'upps_y_ss_lack_of_planning', ...
    'upps_y_ss_lack_of_perseverance', 'upps_y_ss_sensation_seeking'};
UPPS_colloquial = {'Negative urgency', 'Positive urgency', 'Lack of planning', ...
    'Lack of perseverance', 'Sensation seeking'};
for c = 1:length(UPPS_colloquial)
    UPPS_col_plot{c} = regexprep(UPPS_colloquial{c}, ' +', '_');
end
subj_hdr = 'subjectkey';
event_hdr = 'eventname';

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = fullfile(proj_dir, 'scripts', 'lists', 'subjects_rs_censor.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_list);

% format in subj_list: e.g. sub-NDARINV007W6H7B
% format in csv: e.g. "NDAR_INV007W6H7B"
subjects_csv = cell(nsub, 1);
for s = 1:nsub
    subjects_csv{s} = [subjects{s}(5:8) '_' subjects{s}(9:end)];
end

d = readtable(UPPS_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');

% choose columns of selected UPPS measures
UPPS_read = [];
for c = 1:length(UPPS_hdr)
    curr_UPPS = d.(UPPS_hdr{c});
    UPPS_read = [UPPS_read curr_UPPS];
end

% select only the rows corresponding to required subjects
UPPS = cell(nsub, length(UPPS_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        UPPS(s,:) = UPPS_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, UPPS);
UPPS(empty_idx) = {'NaN'};
UPPS = cellfun(@str2num, UPPS);

if(dohist==1)
    binwidth = 1;
    nan_replace = 0;
    
    for c = 1:length(UPPS_hdr)
        UPPS_plot = UPPS(:,c);
        % assign NaN to 10 (an invalid number for this task), so that #subjects without scores can be plotted
        UPPS_plot(isnan(UPPS_plot)) = nan_replace; 
        
        %% histogram across all subjects
        h = histogram(UPPS_plot, 'binwidth', binwidth);
        box off
        set(gcf, 'Position', [0 0 1500 600])
        E = h.BinEdges;
        y = h.BinCounts; 
        start_idx = find(y~=0); start_idx = start_idx(2)-1;
        y(2:start_idx) = [];
        
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        xloc_txt = E(1:end-1); xloc_txt(2:start_idx) = [];
        text(xloc_txt, y+30, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%d', [nan round(xloc(2:end))]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [UPPS_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_UPPS = UPPS_plot(WA_filter);
        AA_UPPS = UPPS_plot(AA_filter);
        
        hc_WA = histcounts(WA_UPPS, E); hc_WA(2:start_idx) = [];
        hc_AA = histcounts(AA_UPPS, E); hc_AA(2:start_idx) = [];
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1500 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%d', [nan round(xloc(2:end))]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 (start_idx+1):end]);
        text(xloc-binwidth/2, hc_WA+30, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+30, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [UPPS_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop abcd_mhy02.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

