function [RAVLT, RAVLT_hdr, RAVLT_colloquial] = ABCD_read_RAVLT(subj_list, race, dohist, hist_dir, hist_fstem)

% [RAVLT, RAVLT_hdr, RAVLT_colloquial] = ABCD_read_RAVLT(subj_list, race, dohist, hist_dir, hist_fstem)
%
% Read necessary Rey Auditory Verbal Learning Test scores.
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
% - RAVLT
%   A #subjects x 2 matrix. Each row corresponds to the behavioral scores of a subject.
%   The 2 columns correspond to Short delay recall, Long delay recall.
%   
% - RAVLT_hdr
%   A 1x2 cell. Headers of these 2 measures in ABCD csv file.
%
% - RAVLT_colloquial
%   A 1x2 cell. Colloquial names of these 2 measures.
%
% Example:
% RAVLT = ABCD_read_RAVLT([], race, [], ...
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
% where "race" is obtained from
% race = ABCD_read_race([], [], ...
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');
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
system('datalad get -s inm7-storage abcd_ps01.txt')

RAVLT_csv = fullfile(csv_dir, 'phenotype', 'abcd_ps01.txt');
RAVLT_hdr = {'pea_ravlt_sd_trial_vi_tc', 'pea_ravlt_ld_trial_vii_tc'};
RAVLT_colloquial = {'Short delay recall', 'Long delay recall'};
for c = 1:length(RAVLT_colloquial)
    RAVLT_col_plot{c} = regexprep(RAVLT_colloquial{c}, ' +', '_');
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

d = readtable(RAVLT_csv);

% choose columns of selected RAVLT measures
RAVLT_read = [];
for c = 1:length(RAVLT_hdr)
    curr_RAVLT = d.(RAVLT_hdr{c});
    RAVLT_read = [RAVLT_read curr_RAVLT];
end

% select only the rows corresponding to required subjects
baseline_idx = strcmp(d.(ses_hdr), ses);
RAVLT = cell(nsub, length(RAVLT_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    tmp_idx = tmp_idx & baseline_idx;
    if(any(tmp_idx==1))
        RAVLT(s,:) = RAVLT_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, RAVLT);
RAVLT(empty_idx) = {'NaN'};
RAVLT = cellfun(@str2num, RAVLT);

if(dohist==1)
    for c = 1:length(RAVLT_hdr)
        RAVLT_plot = RAVLT(:,c);
        % assign NaN to -5 (an invalid number for this task), so that #subjects without scores can be plotted
        RAVLT_plot(isnan(RAVLT_plot)) = -5; 
        
        %% histogram across all subjects
        h = histogram(RAVLT_plot);
        box off
        set(gcf, 'Position', [0 0 800 600])
        E = h.BinEdges;
        y = h.BinCounts; y(2:5) = [];
        xloc = unique(RAVLT_plot);
        xloc_txt = E(1:end-1); xloc_txt(2:5) = [];
        text(xloc_txt, y+20, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%d', [nan; xloc(2:end)]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [RAVLT_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_RAVLT = RAVLT_plot(WA_filter);
        AA_RAVLT = RAVLT_plot(AA_filter);
        
        hc_WA = histcounts(WA_RAVLT, E); hc_WA(2:5) = [];
        hc_AA = histcounts(AA_RAVLT, E); hc_AA(2:5) = [];
        xloc = unique([WA_RAVLT; AA_RAVLT]);
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1300 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%d', [nan; xloc(2:end)]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 6:end]);
        text(xloc-0.4, hc_WA+10, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+10, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [RAVLT_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop abcd_ps01.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

