function [NIH, NIH_hdr, NIH_colloquial] = ABCD_read_NIH(subj_list, race, dohist, hist_dir, hist_fstem)

% [NIH, NIH_hdr, NIH_colloquial] = ABCD_read_NIH(subj_list, race, dohist, hist_dir, hist_fstem)
% 
% Read necessary NIH Toolbox scores. Create histograms of behavioral distributions by race.
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
% - NIH
%   A #subjects x 11 matrix. Each row corresponds to the behavioral scores of a subject.
%   The 11 columns correspond to Cognitive control/Attention (Flanker), Working memory (list sort),
%   Executive function (card sort), Reading (pronunciation), Processing speed, 
%   Visual episodic memory, Picture vocabulary, Fluid cognition, Crystallized cognition',
%   Overall cognition.
%
% - NIH_hdr
%   A 1x11 cell. Headers of these 11 measures in ABCD csv file.
%
% - NIH_colloquial
%   A 1x11 cell. Colloquial names of these 11 measures.
% 
% Example:
% [NIH, NIH_hdr, NIH_colloquial] = ABCD_read_NIH([], race, [], ...
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
system('datalad get -s inm7-storage abcd_tbss01.txt')

NIH_csv = fullfile(csv_dir, 'phenotype', 'abcd_tbss01.txt');
NIH_hdr = {'nihtbx_flanker_uncorrected', 'nihtbx_list_uncorrected', 'nihtbx_cardsort_uncorrected', ...
    'nihtbx_reading_uncorrected', 'nihtbx_pattern_uncorrected', 'nihtbx_picture_uncorrected', ...
    'nihtbx_picvocab_uncorrected', 'nihtbx_fluidcomp_uncorrected', 'nihtbx_cryst_uncorrected', ...
    'nihtbx_totalcomp_uncorrected'};
NIH_colloquial = {'Cognitive control, Attention (Flanker)', 'Working memory (list sort)', ...
    'Executive function (card sort)', 'Reading (pronunciation)', 'Processing speed', ...
    'Visual episodic memory', 'Picture vocabulary', 'Fluid cognition', 'Crystallized cognition', ...
    'Overall cognition'};
for c = 1:length(NIH_colloquial)
    NIH_col_plot{c} = regexprep(NIH_colloquial{c}, ' +', '_');
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

d = readtable(NIH_csv);

% choose columns of selected NIH measures
NIH_read = [];
for c = 1:length(NIH_hdr)
    curr_NIH = d.(NIH_hdr{c});
    NIH_read = [NIH_read curr_NIH];
end

% select only the rows corresponding to required subjects
baseline_idx = strcmp(d.(ses_hdr), ses);
NIH = cell(nsub, length(NIH_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    tmp_idx = tmp_idx & baseline_idx;
    if(any(tmp_idx==1))
        NIH(s,:) = NIH_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, NIH);
NIH(empty_idx) = {'NaN'};
NIH = cellfun(@str2num, NIH);

if(dohist==1)
    for c = 1:length(NIH_hdr)
        NIH_plot = NIH(:,c);
        % assign NaN to 10 (an invalid number for this task), so that #subjects without scores can be plotted
        NIH_plot(isnan(NIH_plot)) = 10; 
        
        %% histogram across all subjects
        h = histogram(NIH_plot, 'binwidth', 4);
        box off
        set(gcf, 'Position', [0 0 1000 600])
        E = h.BinEdges;
        y = h.BinCounts; 
        start_idx = find(y~=0); start_idx = start_idx(2)-1;
        y(2:start_idx) = [];
        
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        xloc_txt = E(1:end-1); xloc_txt(2:start_idx) = [];
        text(xloc_txt, y+20, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%d', [nan xloc(2:end)]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [NIH_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_NIH = NIH_plot(WA_filter);
        AA_NIH = NIH_plot(AA_filter);
        
        hc_WA = histcounts(WA_NIH, E); hc_WA(2:start_idx) = [];
        hc_AA = histcounts(AA_NIH, E); hc_AA(2:start_idx) = [];
        xloc = E(1:end-1) + diff(E)/2; xloc(2:start_idx) = [];
        bar(xloc, [hc_WA; hc_AA]')
        box off
        set(gcf, 'Position', [0 0 1300 600])
        
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
        xticklabels(sprintfc('%d', [nan xloc(2:end)]))
        
        legend({'WA', 'AA'}, 'FontSize', 13)
        legend boxoff
        
        xloc_txt = xloc([1 (start_idx+1):end]);
        text(xloc-2, hc_WA+10, string(hc_WA), 'FontSize', 13)
        text(xloc, hc_AA+10, string(hc_AA), 'FontSize', 13)
        
        fname2 = fullfile(hist_dir, [NIH_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop abcd_tbss01.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

