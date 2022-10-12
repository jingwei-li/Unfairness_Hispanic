function ABCD_collect_RSFC(parc_dir, subj_censor, subj_ls, scale, out_fc)

% ABCD_collect_RSFC(parc_dir, subj_censor, out_rsfc)
%
% Average RSFC across runs within each subject. Concatenate RSFC across subjects.
%
% Input:
%   - parc_dir
%     Directory of parcellated timeseries and RSFC of each subject and each run.
%   - subj_censor
%     A .mat file containing the information of which runs of each subject have passed
%     the censoring criterion. It is the output file of 
%     https://github.com/jingwei-li/Parcellate_ABCD_DCANpreproc/blob/main/compute_RSFC_with_censor.m
%     It should contain at leaset the following variables:
%     `subjects`, a cell, list of subjects before censoring;
%     `subjects_pass`, a cell, list of subjects after censoring;
%     `pass_runs`, a cell with same length as `subjects`, which runs passed censoring criterion for each subject.
%   - subj_ls
%     A list of subjects of interest, which should be a subset of `subjects_pass`.
%   - scale
%     choose from 1 to 10. 
%     Scale 1: Schaefer parcellation with 100 areas + Tian parcellation at scale 1 (16 ROIs)
%     Scale 2: Schaefer parcellation with 200 areas + Tian parcellation at scale 2 (32 ROIs)
%     Scale 3: Schaefer parcellation with 300 areas + Tian parcellation at scale 3 (50 ROIs)
%     From 4 to 10: Schaefer parcellation with (scale * 100) areas + Tian parcellation at scale 4 (54 ROIs)
%
% Output:
%   - out_fc
%     Path of output RSFC mat file.

%% add external scripts folder to path
repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_dir, 'external')))

%% default parameters
if(isempty(parc_dir))
    parc_dir = '/data/project/parcellate_ABCD_preprocessed/data/parcellated_timeseries';
end

if(~isempty(subj_censor))
    subj_censor = '/data/project/parcellate_ABCD_preprocessed/scripts/lists/subjects_rs_censor.mat';
end

if(isempty(scale))
    scale = 2;
end

ses = 'ses-baselineYear1Arm1';
if(ischar(scale))
    scale = str2num(scale);
end
Schaefer_res = num2str(100*scale);
if(scale<4)
    Tian_res = num2str(scale);
else
    Tian_res = '4';
end

%% collect and average RSFC
load(subj_censor)
soi = CBIG_text2cell(subj_ls);
[~,~,idx1] = intersect(soi, subjects_pass, 'stable');
[~,~,idx2] = intersect(subjects_pass, subjects, 'stable');
corr_mat = [];
for i = 1:length(soi)
    s = soi{i};
    runs = pass_runs{idx2(idx1(i))};
    curr_fc = [];
    for j = 1:length(runs)
        fcname = fullfile(parc_dir, s, ses, 'func', [s '_' ses '_task-rest_' runs{j} ...
            '_RSFC_Schaefer' Schaefer_res '_Tian' Tian_res '.mat']);
        fc = load(fcname)
        curr_fc = cat(3, curr_fc, fc.corr_mat);
    end
    curr_fc = CBIG_StableAtanh(curr_fc);
    curr_fc = mean(curr_fc, 3);
    corr_mat = cat(3, corr_mat, curr_fc);
end

outdir = fileparts(out_fc)
if(~exist(outdir, 'dir'))
    mkdir(outdir)
end
save(out_fc, 'corr_mat', '-v7.3')

end