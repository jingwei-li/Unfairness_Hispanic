function [BISBAS, BISBAS_hdr, BISBAS_colloquial] = ABCD_read_BISBAS(subj_list, race, dohist, hist_dir, hist_fstem)

% [BISBAS, BISBAS_hdr, BISBAS_colloquial] = ABCD_read_BISBAS(subj_list, race, dohist, hist_dir, hist_fstem)
% 
% Read BIS/BAS scales from ABCD csv file. Create histograms of BIS/BAS scores by race.
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
% - BIABAS
%   A #subjects x 4 matrix. Each row is the behavioral scores of a subject. The 4 columns correspond
%   to Behavioral inhibition, BAS - Reward responsiveness, BAS - Drive, BAS - Fun seeking.
% 
% - BISBAS_hdr
%   A 1x4 cell. Headers of these 4 measures in ABCD csv file.
%
% - BISBAS_colloquial
%   A 1x4 cell. Colloquial names of these 4 measures.
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

BISBAS_csv = fullfile(csv_dir, 'phenotype', 'abcd_mhy02.txt');
BISBAS_hdr = {'bis_y_ss_bis_sum', 'bis_y_ss_bas_rr', 'bis_y_ss_bas_drive', 'bis_y_ss_bas_fs'};
BISBAS_colloquial = {'Behavioral inhibition ', 'BAS - Reward responsiveness', 'BAS - Drive', 'BAS - Fun seeking'};
for c = 1:length(BISBAS_colloquial)
    BISBAS_col_plot{c} = regexprep(BISBAS_colloquial{c}, ' +', '_');
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

d = readtable(BISBAS_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');

% choose columns of selected BISBAS measures
BISBAS_read = [];
for c = 1:length(BISBAS_hdr)
    curr_BISBAS = d.(BISBAS_hdr{c});
    BISBAS_read = [BISBAS_read curr_BISBAS];
end

% select only the rows corresponding to required subjects
BISBAS = cell(nsub, length(BISBAS_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s}) & base_event;
    if(any(tmp_idx==1))
        BISBAS(s,:) = BISBAS_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, BISBAS);
BISBAS(empty_idx) = {'NaN'};
BISBAS = cellfun(@str2num, BISBAS);

if(dohist==1)
    binwidth = 1;
    nan_replace = -2;
    
    for c = 1:length(BISBAS_hdr)
        BISBAS_plot = BISBAS(:,c);
        % assign NaN to 10 (an invalid number for this task), so that #subjects without scores can be plotted
        BISBAS_plot(isnan(BISBAS_plot)) = nan_replace; 
        
        %% histogram across all subjects
        h = histogram(BISBAS_plot, 'binwidth', binwidth);
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
        [imageData, alpha] = export_fig(fullfile(hist_dir, [BISBAS_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_BISBAS = BISBAS_plot(WA_filter);
        AA_BISBAS = BISBAS_plot(AA_filter);
        
        hc_WA = histcounts(WA_BISBAS, E); hc_WA(2:start_idx) = [];
        hc_AA = histcounts(AA_BISBAS, E); hc_AA(2:start_idx) = [];
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
        
        fname2 = fullfile(hist_dir, [BISBAS_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop abcd_mhy02.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

