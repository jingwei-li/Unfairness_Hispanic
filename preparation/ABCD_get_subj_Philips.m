function [subj_philips, isPhilips] = ABCD_get_subj_Philips(subj_list)

% [subj_philips, isPhilips] = ABCD_get_subj_Philips(subj_list)
%
% Subjects collected on Philips scanners were problematic:
% https://github.com/ABCD-STUDY/fMRI-cleanup.
% These subjects need to be excluded from the study.
% This function checks if any subject in the given list was collected on
% Philips scanners.
%
% Inputs:
% - subj_list
%   List of subjects with all required phenotypes, and passed rs-fMRI quality 
%   control (full path).
%
% Outputs:
% - subj_philips
%   The subject IDs which were collected on Philips scanners, if any.
%
% - isPhilips
%   A #subjects x 1 boolean vector. An entry of True indicates that this 
%   subject was collected on Philips scanners.
%
% Author: Jingwei Li 

repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_dir, 'external')))
start_dir = pwd;

proj_dir = '/data/project/FairAI_Hispanic';
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

csv_dir = fullfile(proj_dir, 'data', 'inm7-superds', 'original', 'abcd', 'phenotype');
cd(csv_dir)
system('datalad get -n .')
system('git -C . config --local --add remote.datalad.annex-ignore true')
cd(fullfile(csv_dir, 'phenotype'))
system('datalad get -s inm7-storage abcd_mri01.txt')

mri_csv = fullfile(csv_dir, 'phenotype', 'abcd_mri01.txt');
subj_hdr = 'subjectkey';
scanner_hdr = 'mri_info_manufacturer';
event_hdr = 'eventname';

d = readtable(mri_csv);
base_event = strcmp(d.(event_hdr), 'baseline_year_1_arm_1');
scanner = cell(nsub,1);
for s = 1:nsub
    tmp_idx = strcmp(d.(subj_hdr), subjects_csv{s}) & base_event;
    if(any(tmp_idx==1))
        curr_scanner = d.(scanner_hdr)(tmp_idx);
        if(length(curr_scanner)>1)
            curr_scanner = curr_scanner(~cellfun('isempty', curr_scanner));
        end
        scanner(s) = curr_scanner;
    end
end

isPhilips = strcmp(scanner, 'Philips Medical Systems');
subj_philips = subjects(isPhilips);


end

