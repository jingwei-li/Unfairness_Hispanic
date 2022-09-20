function [peduc, peduc_comb, peduc_avg, peduc_hdr, peduc_colloquial] = ABCD_read_Prt_Educ(subj_list, race, dohist, hist_fname)

% [peduc, peduc_comb, peduc_avg, peduc_hdr, peduc_colloquial] = ABCD_read_Prt_Educ(subj_list, race, dohist, hist_fname)
% 
% Read the highest degree of the parent and the parterner of the parent.
% 
% Based on ABCD documentation:
%  0 = Never attended/Kindergarten only; 1 = 1st grade; 2 = 2nd grade; 3 = 3rd grade; 4 = 4th grade; 
%  5 = 5th grade; 6 = 6th grade; 7 = 7th grade; 8 = 8th grade; 9 = 9th grade; 10 = 10th grade; 
%  11 = 11th grade; 12 = 12th grade; 13 = High school graduate; 14 = GED or equivalent Diploma; 
%  15 = Some college; 16 = Associate degree: Occupational; 17 = Associate degree: Academic Program; 
%  18 = Bachelor's degree (ex. BA); 19 = Master's degree (ex. MA); 20 = Professional School degree (ex. MD); 
%  21 = Doctoral degree (ex. PhD); 777 = Refused to answer Prefiero no responder ; 999 = Don't Know No
% 
% 777 & 999 are replaced with NaN.
% Then the average between the two parents are computed, returned as
% peduc_avg. Histogram of peduc_avg are plotted for the whole population,
% and for WA or AA separately.
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
% - peduc
%   A #subjects x 2 matrix. Each row is the education level of subject's parents.
% 
% - peduc_comb
%   A #subjects x 1 vector. Concatenatation of the two columns in "peduc".
%
% - peduc_avg
%   A #subjects x 1 vector. Average of the parents' education levels.
%
% - peduc_hdr
%   Header of the two education measures in ABCD csv files.
%
% - peduc_colloquial
%   Collquial names of the two education measures.
%
% Example:
% [peduc, peduc_comb, peduc_avg, peduc_hdr, peduc_colloquial] = ABCD_read_Prt_Educ([], race, [], ...
%     '/data/users/jingweil/storage/MyProject/fairAI/ABCD_race/figures/demo_hist/Prt_Educ_avg_pass_rs.png')
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

peduc_csv = fullfile(csv_dir, 'phenotype', 'pdem02.txt');
peduc_hdr = {'demo_prnt_ed_v2', 'demo_prtnr_ed_v2'};
peduc_colloquial = {'Parental degree', 'Parterner''s degree'};
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

d = readtable(peduc_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');

% choose columns of selected BISBAS measures
peduc_read = [];
for c = 1:length(peduc_hdr)
    curr_peduc = d.(peduc_hdr{c});
    peduc_read = [peduc_read curr_peduc];
end

% select only the rows corresponding to required subjects
peduc = cell(nsub, length(peduc_hdr));
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s});
    if(any(tmp_idx==1))
        tmp_idx = tmp_idx & base_event;
        peduc(s,:) = peduc_read(tmp_idx,:);
    end
end
empty_idx = cellfun(@isempty, peduc);
peduc(empty_idx) = {'NaN'};
peduc = cellfun(@str2num, peduc);

peduc(peduc>25) = nan;
peduc_comb = sprintfc('%03d', sort(peduc,2));
for l = 1:size(peduc_comb,1)
    peduc_comb{l,1} = strjoin(peduc_comb(l,:), '_');
end
peduc_comb(:,2) = [];
xtl = unique(peduc_comb);

peduc_avg = mean(peduc, 2);
nan_idx = find(isnan(peduc_avg));
for n = 1:length(nan_idx)
    if(~isnan(peduc(nan_idx(n),1)))
        peduc_avg(nan_idx(n)) = peduc(nan_idx(n),1);
    elseif(~isnan(peduc(nan_idx(n),2)))
        peduc_avg(nan_idx(n)) = peduc(nan_idx(n),2);
    end
end

if(dohist == 1)
    binwidth = 2;
    peduc_plot = peduc_avg;
    peduc_plot(isnan(peduc_plot)) = -2;
    
    E = (min(peduc_plot)-binwidth/2):binwidth:(max(peduc_plot)+binwidth/2);
    hc = histcounts(peduc_plot, E);
    bar(E(1:end-1) + diff(E)/2, hc');
    start_idx = find(hc~=0); start_idx = start_idx(2)-1;
    hc(2:start_idx) = [];
    
    box off
    set(gcf, 'Position', [0 0 800 600])
    xloc = E(1:end-1)+diff(E)/2;  xloc(2:start_idx) = [];
    xtl = sprintfc('%d', xloc);
    xloc_txt = E(1:end-1);   xloc_txt(2:start_idx) = [];
    text(xloc_txt, hc+50, string(hc), 'FontSize', 13)
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
    xticklabels(xtl)
    
    [outdir] = fileparts(hist_fname);
    if(~exist(outdir, 'dir'))
        mkdir(outdir);
    end
    [imageData, alpha] = export_fig(hist_fname, '-png', '-nofontswap', '-a1');
    close(gcf)
    
    %% plot AA/WA histograms across peduc_avg
    WA_filter = strcmp(race, '1');
    AA_filter = strcmp(race, '2');
    WA_peduc_avg = peduc_plot(WA_filter);
    AA_peduc_avg = peduc_plot(AA_filter);
    
    hc_WA = histcounts(WA_peduc_avg, E);
    hc_AA = histcounts(AA_peduc_avg, E);
    bar(E(1:end-1) + diff(E)/2, [hc_WA; hc_AA]')
    hc_WA(2:start_idx) = [];  hc_AA(2:start_idx) = [];
    box off
    set(gcf, 'Position', [0 0 1000 600])
    set(gca, 'xtick', xloc, 'linewidth', 2, 'fontsize', 13, 'TickDir', 'out')
    legend({'WA', 'AA'}, 'FontSize', 13, 'location', 'best')
    legend boxoff
    
    text(xloc-binwidth/2, hc_WA+15, string(hc_WA), 'FontSize', 13)
    text(xloc, hc_AA+15, string(hc_AA), 'FontSize', 13)
    xticklabels(xtl)
    
    [~, outbase, outext] = fileparts(hist_fname);
    fname2 = fullfile(outdir, [outbase '_WAvsAA' outext]);
    [imageData, alpha] = export_fig(fname2, '-png', '-nofontswap', '-a1');
    close(gcf)
end

system('datalad drop pdem02.txt')
cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end