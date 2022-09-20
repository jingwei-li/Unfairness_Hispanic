function [income, income_trans] = ABCD_read_income(subj_list, race, dohist, hist_fname)

% [income, income_trans] = ABCD_read_income(subj_list, race, dohist, hist_fname)
% 
% Read household income from ABCD csv file. Create histograms of income distributions
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
% - income
%   A #subjects x 1 cell. Each entry is the household income of a subject.
%
% Example:
% income = ABCD_read_income([], race, [], ...
%     '/home/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/income_pass_rs.png');
% where "race" is obtained from
% race = ABCD_read_race([], [], 
%     '/home/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/race_pass_rs.png');
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
system('datalad get -s inm7-storage pdem02.txt')

income_csv = fullfile(csv_dir, 'phenotype', 'pdem02.txt');
income_hdr = 'demo_comb_income_v2';
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

d = readtable(income_csv);
income = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        income(s) = d.(income_hdr)(tmp_idx);
    end
end
empty_idx = cellfun(@isempty, income);
if(any(empty_idx))
    income{empty_idx} = '999';
end

%% transform categorical data into median income of every group
income_trans = cellfun(@str2num, income);
nan_idx = income_trans==777 | income_trans==999;
income_trans(nan_idx) = nan;
% 1= Less than $5,000
idx = income_trans==1;
income_trans(idx) = 2500;
% 2=$5,000 through $11,999
idx = income_trans==2;
income_trans(idx) = 8500;
% 3=$12,000 through $15,999
idx = income_trans==3;
income_trans(idx) = 14000;
% 4=$16,000 through $24,999
idx = income_trans==4;
income_trans(idx) = 20500;
% 5=$25,000 through $34,999
idx = income_trans==5;
income_trans(idx) = 30000;
% 6=$35,000 through $49,999
idx = income_trans==6;
income_trans(idx) = 42500;
% 7=$50,000 through $74,999
idx = income_trans==7;
income_trans(idx) = 62500;
% 8= $75,000 through $99,999
idx = income_trans==8;
income_trans(idx) = 87500;
% 9=$100,000 through $199,999
idx = income_trans==9;
income_trans(idx) = 150000;
% 10=$200,000 and greater
idx = income_trans==10;
x = 1:9;
y = [2500 8500 14000 20500 30000 42500 62500 87500 150000];
f = @(b,x) b(1).*exp(b(2).*x)+b(3);                              % Objective Function
B = fminsearch(@(b) norm(y - f(b,x)), ones(3,1))
income_trans(idx) = round(f(B, 10)/100)*100;  % use exponentially fitted value

if(dohist==1)
    % 1= Less than $5,000; 2=$5,000 through $11,999; 3=$12,000 through $15,999; 4=$16,000 through $24,999; 
    % 5=$25,000 through $34,999; 6=$35,000 through $49,999; 7=$50,000 through $74,999; 
    % 8= $75,000 through $99,999; 9=$100,000 through $199,999; 10=$200,000 and greater. 
    % 999 = Don't know; 777 = Refuse to answer
    xtl = {'Refuse ans', 'Unknown', '<$5k', '$5k-12k', '$12k-16k', '$16k-25k', ...
        '$25k-35k', '$35k-50k', '$50k-75k', '$75k-100k', '$100k-200k', '>=$200k'};
    
    %% plot household income of all subjects
    income_plot = cellfun(@str2num, income);
    income_plot(income_plot==777) = -0.99;
    income_plot(income_plot==999) = 0;
    
    h = histogram(income_plot, 'BinWidth', 0.9999);
    box off
    set(gcf, 'Position', [0 0 1400 600])
    E = h.BinEdges;
    y = h.BinCounts;
    xloc = E(1:end-1)+diff(E)/2;
    xloc_txt = E(1:end-1) + diff(E)/3;
    text(xloc_txt, y+50, string(y), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 12, 'TickDir','out')
    xticklabels(xtl)
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
    
    %% plot household income for WA/AA separately
    WA_filter = strcmp(race, '1');
    AA_filter = strcmp(race, '2');
    WA_income = income_plot(WA_filter);
    AA_income = income_plot(AA_filter);
    
    hc_WA = histcounts(WA_income, E);
    hc_AA = histcounts(AA_income, E);
    bar(E(1:end-1), [hc_WA; hc_AA]')
    box off
    set(gcf, 'Position', [0 0 1300 600])
    xloc_txt = min(round(income_plot)):1:max(income_plot);
    set(gca, 'xtick', xloc_txt, 'linewidth', 2, 'fontsize', 12, 'TickDir', 'out')
    xticklabels(xtl)
    
    legend({'WA', 'AA'}, 'FontSize', 13)
    legend boxoff
    
    text(xloc_txt-0.4, hc_WA+30, string(hc_WA), 'FontSize', 13)
    text(xloc_txt, hc_AA+30, string(hc_AA), 'FontSize', 13)
    
    [~, outbase, outext] = fileparts(hist_fname);
    fname2 = fullfile(outdir, [outbase '_WAvsAA' outext]);
    [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
    close(gcf)
end

system('datalad drop pdem02.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

