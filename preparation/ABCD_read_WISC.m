function [WISC, WISC_hdr, WISC_colloquial] = ABCD_read_WISC(subj_list, race, dohist, hist_dir, hist_fstem)

% [WISC, WISC_hdr, WISC_colloquial] = ABCD_read_WISC(subj_list, race, dohist, hist_dir, hist_fstem)
%
% Read necessary measures from Wechsler Intelligence Scale for Children-V.
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
% - WISC
%   A #subjects x 1 vector. Each row corresponds to the Matrix reasoning score of a subject.
%
% - WISC_hdr
%   A 1x1 cell. Header of the behavioral measure in ABCD csv file.
%  
% - WISC_colloquial
%   A 1x1 cell. Colloquial name of the behavioral measure.
%
% Example:
% [WISC, WISC_hdr, WISC_colloquial] = ABCD_read_WISC([], race, [], ...
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist', '_pass_rs');
% where "race" is obtained from
% race = ABCD_read_race([], [], ...
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');

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

WISC_csv = fullfile(csv_dir, 'phenotype', 'abcd_ps01.txt');
WISC_hdr = {'pea_wiscv_trs'};
WISC_colloquial = {'Matrix reasoning'};
for c = 1:length(WISC_colloquial)
    WISC_col_plot{c} = regexprep(WISC_colloquial{c}, ' +', '_');
end
subj_hdr = 'subjectkey';

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

d = readtable(WISC_csv);

% choose columns of selected WISC measures
WISC_read = [];
for c = 1:length(WISC_hdr)
    curr_WISC = d.(WISC_hdr{c});
    WISC_read = [WISC_read curr_WISC];
end

% select only the rows corresponding to required subjects
WISC = cell(nsub, length(WISC_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        WISC(s,:) = WISC_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, WISC);
WISC(empty_idx) = {'NaN'};
WISC = cellfun(@str2num, WISC);

if(dohist==1)
    for c = 1:length(WISC_hdr)
        WISC_plot = WISC(:,c);
        % assign NaN to -5 (an invalid number for this task), so that #subjects without scores can be plotted
        WISC_plot(isnan(WISC_plot)) = -5; 
        
        %% histogram across all subjects
        h = histogram(WISC_plot);
        box off
        set(gcf, 'Position', [0 0 800 600])
        E = h.BinEdges;
        y = h.BinCounts; y(2:5) = [];
        xloc = unique(WISC_plot);
        xloc_txt = E(1:end-1); xloc_txt(2:5) = [];
        text(xloc_txt, y+20, string(y), 'FontSize', 13)
        set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
        xticklabels(sprintfc('%d', [nan; xloc(2:end)]))
        
        if(~exist(hist_dir, 'dir'))
            mkdir(hist_dir);
        end
        [imageData, alpha] = export_fig(fullfile(hist_dir, [WISC_col_plot{c} hist_fstem '.png']), '-png', '-nofontswap', '-a1');
        close(gcf)
        
        %% histogram for AA/WA separately
        WA_filter = strcmp(race, '1');
        AA_filter = strcmp(race, '2');
        WA_WISC = WISC_plot(WA_filter);
        AA_WISC = WISC_plot(AA_filter);
        
        hc_WA = histcounts(WA_WISC, E); hc_WA(2:5) = [];
        hc_AA = histcounts(AA_WISC, E); hc_AA(2:5) = [];
        xloc = unique([WA_WISC; AA_WISC]);
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
        
        fname2 = fullfile(hist_dir, [WISC_col_plot{c} hist_fstem '_WAvsAA.png']);
        [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
        close(gcf)
    end
end

system('datalad drop abcd_ps01.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

