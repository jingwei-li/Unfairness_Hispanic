function [race, race_hdr] = ABCD_read_race(subj_list, dohist, hist_fname)

% [race, race_hdr] = ABCD_read_race(subj_list, dohist, hist_fname)
%
% Read ethnicities/races of all subjects. Plot the distribution of ethnicity/race.
%
% Inputs:
% - subj_list
%   List of subjects which passed fMRI prepreocessing quality control (full path). Default:
%   '/mnt/eql/yeo13/data/ABCD/orig_scripts/release2.0/lists/subjects_pass_rs.txt'
%
% - dohist
%   A 1/0 value determining whether the histograms are created or not. 
%   Default: 1, i.e. create plots.
%
% - hist_fname
%   Full path of output histogram filename.
%
% Outputs:
% - race
%   A #subjects x 1 cell. Each entry contains the ethnicity/race of a subject.
%   
% - race_hdr
%   A string. Header of the ethnicity/race column in ABCD csv file.
%
% Example:
% race = ABCD_read_race([], [], ...
%     '/home/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');
%
% race = ABCD_read_race('/home/jingweil/storage/MyProject/fairAI/ABCD_race/scripts/lists/subjects_pass_rs_pass_pheno.txt',...
%    [], '/home/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs_pass_pheno.png');
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
system('datalad get -s inm7-storage acspsw03.txt')

race_csv = fullfile(csv_dir, 'phenotype', 'acspsw03.txt');
race_hdr = 'race_ethnicity';
subj_hdr = 'subjectkey';
ses_hdr = 'eventname';

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

race_labels = {'', '1', '2', '3', '4', '5'};
race_grps = {'Empty', 'White', 'Black', 'Hispanic', 'Asian', 'Other'};
d = readtable(race_csv);
ses = 'baseline_year_1_arm_1';
race = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s}) & strcmp(d.(ses_hdr), ses);
    if(any(tmp_idx==1))
        race(s) = d.(race_hdr)(tmp_idx);
    end
end

if(dohist)
    uniq_race = unique(race);
    [~,idx] = intersect(race_labels, uniq_race);
    xtl = race_grps(idx);
    race_plot = race;
    for l = 1:length(race_labels)
        idx = strcmp(race, race_labels{l});
        race_plot(idx) = {l-1};
    end
    race_plot = cell2mat(race_plot);
    
    h = histogram(race_plot, 'BinWidth', 0.9999);
    box off
    set(gcf, 'Position', [0 0 600 600])
    E = h.BinEdges;
    y = h.BinCounts;
    xloc = E(1:end-1)+diff(E)/2;
    xloc_txt = E(1:end-1)+diff(E)/3;
    text(xloc_txt, y+100, string(y), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
    xticklabels(xtl);
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close
end

system('datalad drop acspsw03.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

