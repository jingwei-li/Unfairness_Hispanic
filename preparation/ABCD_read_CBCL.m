function [CBCL, CBCL_hdr, CBCL_colloquial] = ABCD_read_CBCL(subj_list, race, dohist, hist_dir, hist_fstem)

% [CBCL, CBCL_hdr, CBCL_colloquial] = ABCD_read_CBCL(subj_list, race, dohist, hist_dir, hist_fstem)
%
% Read and plot histograms of the necessary measures from Achenbach Child Behavior Check List.
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
% - CBCL
%   A #subjects x 8 matrix. Each row is the behavioral scores of a subject. The 8 columns 
%   correspond to Anxious,Depressed, Withdrawn,Depressed, Somatic complaints, Social problems, 
%   Thought problems, Attention problems, Rule-breaking behavior, Aggressive behavior
%
% - CBCL_hdr
%   A 1x8 cell. Headers of these 8 measures in ABCD csv file.
%
% - CBCL_colloquial
%   A 1x8 cell. Colloquial names of these 8 measures.
%
% Example:
% [CBCL, CBCL_hdr, CBCL_colloquial] = ABCD_read_CBCL([], race, [], ...
%     '~/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
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
system('datalad get -s inm7-storage abcd_cbcls01.txt')

CBCL_csv = fullfile(csv_dir, 'phenotype', 'abcd_cbcls01.txt');
CBCL_hdr = {'cbcl_scr_syn_anxdep_r', 'cbcl_scr_syn_withdep_r', 'cbcl_scr_syn_somatic_r', ...
    'cbcl_scr_syn_social_r', 'cbcl_scr_syn_thought_r', 'cbcl_scr_syn_attention_r', ...
    'cbcl_scr_syn_rulebreak_r', 'cbcl_scr_syn_aggressive_r'};
CBCL_colloquial = {'Anxious,Depressed', 'Withdrawn,Depressed', 'Somatic complaints', ...
    'Social problems', 'Thought problems', 'Attention problems', 'Rule-breaking behavior', ...
    'Aggressive behavior'};
for c = 1:length(CBCL_colloquial)
    CBCL_col_plot{c} = regexprep(CBCL_colloquial{c}, ' +', '_');
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

d = readtable(CBCL_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');

% choose columns of selected CBCL measures
CBCL_read = [];
for c = 1:length(CBCL_hdr)
    curr_CBCL = d.(CBCL_hdr{c});
    CBCL_read = [CBCL_read curr_CBCL];
end

% select only the rows corresponding to required subjects
CBCL = cell(nsub, length(CBCL_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s}) & base_event;
    if(any(tmp_idx==1))
        CBCL(s,:) = CBCL_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, CBCL);
CBCL(empty_idx) = {'NaN'};
CBCL = cellfun(@str2num, CBCL);

if(dohist==1)
    binwidth = 1;
    nan_replace = -2;
    
    for c = 1:length(CBCL_hdr)
        CBCL_plot = CBCL(:,c);
        % assign NaN to 10 (an invalid number for this task), so that #subjects without scores can be plotted
        CBCL_plot(isnan(CBCL_plot)) = nan_replace; 
        
        %% histogram across all subjects
        h = histogram(CBCL_plot, 'binwidth', binwidth);
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
        xticklabels(sprintfc('%d', [nan round(xloc(2:end))]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [CBCL_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_CBCL = CBCL_plot(WA_filter);
        AA_CBCL = CBCL_plot(AA_filter);
        
        hc_WA = histcounts(WA_CBCL, E); hc_WA(2:start_idx) = [];
        hc_AA = histcounts(AA_CBCL, E); hc_AA(2:start_idx) = [];
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1500 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%d', [nan round(xloc(2:end))]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 (start_idx+1):end]);
        text(xloc-binwidth/2, hc_WA+10, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+10, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [CBCL_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop abcd_cbcls01.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

