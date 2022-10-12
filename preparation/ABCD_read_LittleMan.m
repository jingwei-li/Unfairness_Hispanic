function [LMT, LMT_hdr, LMT_colloquial] = ABCD_read_LittleMan(subj_list, race, dohist, hist_dir, hist_fstem)

% [LMT, LMT_hdr, LMT_colloquial] = ABCD_read_LittleMan(subj_list, race, dohist, hist_dir, hist_fstem)
% 
% Read necessary measures from Little Man Task. Create histograms of behavioral distributions
% by race.
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
% - LMT
%   A #subjects x 3 matrix. Each row is the behavioral scores of a subject. The 3 columns 
%   correspond to Visuospatial accuracy, Visuospatial reaction time, Visuospatial efficiency.
%
% - LMT_hdr
%   A 1x3 cell. Headers of these 3 measures in ABCD csv file.
%
% - LMT_colloquial
%   A 1x3 cell. Colloquial names of these 3 measures.
% 
% Example:
% [LMT, LMT_hdr, LMT_colloquial] = ABCD_read_LittleMan([], race, [], ...
%     '~/storage/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
% where "race" is obtained from
% race = ABCD_read_race([], [], 
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
system('datalad get -s inm7-storage lmtp201.txt')

LMT_csv = fullfile(csv_dir, 'phenotype', 'lmtp201.txt');
LMT_hdr = {'lmt_scr_perc_correct', 'lmt_scr_rt_correct'}; % 'lmt_scr_efficiency' is removed because this column is all empty in the csv file from DCAN lab
LMT_colloquial = {'Visuospatial accuracy', 'Visuospatial reaction time'}; % therefore, 'Visuospatial efficiency' is removed
for c = 1:length(LMT_colloquial)
    LMT_col_plot{c} = regexprep(LMT_colloquial{c}, ' +', '_');
end
subj_hdr = 'subjectkey';
ses_hdr = 'eventname';
ses = 'baseline_year_1_arm_1';

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

d = readtable(LMT_csv);

% choose columns of selected LMT measures
LMT_read = [];
for c = 1:length(LMT_hdr)
    curr_LMT = d.(LMT_hdr{c});
    LMT_read = [LMT_read curr_LMT];
end

% select only the rows corresponding to required subjects
baseline_idx = strcmp(d.(ses_hdr), ses);
LMT = cell(nsub, length(LMT_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    tmp_idx = tmp_idx & baseline_idx;
    if(any(tmp_idx==1))
        LMT(s,:) = LMT_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, LMT);
LMT(empty_idx) = {'NaN'};
LMT = cellfun(@str2num, LMT);

if(dohist==1)
    binwidth = [0.1 200 4e-5];
    nan_replace = [-0.2 600 -5e-5];
    
    for c = 1:length(LMT_hdr)
        LMT_plot = LMT(:,c);
        % assign NaN to 10 (an invalid number for this task), so that #subjects without scores can be plotted
        LMT_plot(isnan(LMT_plot)) = nan_replace(c); 
        
        %% histogram across all subjects
        h = histogram(LMT_plot, 'binwidth', binwidth(c));
        box off
        set(gcf, 'Position', [0 0 1500 600])
        E = h.BinEdges;
        y = h.BinCounts; 
        start_idx = find(y~=0); start_idx = start_idx(2)-1;
        y(2:start_idx) = [];
        
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        xloc_txt = E(1:end-1); xloc_txt(2:start_idx) = [];
        text(xloc_txt, y+20, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%.1e', [nan xloc(2:end)]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [LMT_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_LMT = LMT_plot(WA_filter);
        AA_LMT = LMT_plot(AA_filter);
        
        hc_WA = histcounts(WA_LMT, E); hc_WA(2:start_idx) = [];
        hc_AA = histcounts(AA_LMT, E); hc_AA(2:start_idx) = [];
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1500 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%.1e', [nan xloc(2:end)]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 (start_idx+1):end]);
        text(xloc-binwidth(c)/2, hc_WA+10, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+10, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [LMT_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop lmtp201.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

