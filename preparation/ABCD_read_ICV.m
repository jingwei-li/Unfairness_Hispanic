function ICV = ABCD_read_ICV(subj_list, race, txtname, dohist, hist_fname, FSdir)

% ICV = ABCD_read_ICV(subj_list, race, dohist, hist_fname, FSdir)
%
% Read intracranial volume of all subjects from FreeSurfer recon-all results.
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
% - txtdir
%   The full-path filename of the output text file containing all ICV values.
% 
% - dohist
%   A 1/0 value determining whether the histograms are created or not. 
%   Default: 1, i.e. create plots.
%
% - hist_fname
%   Full path of output histogram filename.
%
% - FSdir
%   Full path of the directory storing all output files from FreeSurfer recon-all
%   when preprocessing structural MRI.
%
% Outputs:
% - ICV
%   A #subjects x 1 vector. Each entry is the intracranial volume of a subject.
%
% Author: Jingwei Li

repo_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_dir, 'external')))
start_dir = pwd;

if(~exist('dohist', 'var') || isempty(dohist))
    dohist = 1;
end

proj_dir = '/data/project/FairAI_Hispanic';

if(~exist('subj_list', 'var') || isempty(subj_list))
    subj_list = fullfile(proj_dir, 'scripts', 'lists', 'subjects_rs_censor.txt');
end
[subjects, nsub] = CBIG_text2cell(subj_list);

flag = 0;
if(exist(txtname, 'file'))
    ICV = dlmread(txtname);
    if(length(ICV) == nsub)
        flag = 1;
    end
end

if(flag == 0)
    if(~exist('FSdir', 'var') || isempty(FSdir))
        FSdir = fullfile(proj_dir, 'data', 'inm7-superds', 'original', 'abcd', 'derivatives', 'freesurfer-5.3.0-HCP');
    end
    cd(FSdir)
    system('datalad get -n .')
    system('git -C . config --local --add remote.datalad.annex-ignore true')

    ICV = nan(nsub,1);
    ses = 'ses-baselineYear1Arm1';
    for i = 1:nsub
        s = subjects{i}
        system(sprintf('datalad get -n %s', s))
        system(sprintf('git -C %s config --local --add remote.datalad.annex-ignore true', s))
        system(sprintf('datalad get -s inm7-storage %s/%s/stats/aseg.stats', s, ses))
        
        if(~exist(fullfile(FSdir, s, ses, 'stats', 'aseg.stats'), 'file'))
            warning('%s does not exist.\n', fullfile(FSdir, s, ses, 'stats', 'aseg.stats'))
            continue
        end

        FSstats = CBIG_text2cell(fullfile(FSdir, s, ses, 'stats', 'aseg.stats'));
        ICV_row = contains(FSstats, 'Intracranial Volume');
        if(exist('ICV_row', 'var') && any(ICV_row==1))
            ICV_row = FSstats{ICV_row==1};
            row_split = strsplit(ICV_row, ',');
            ICV(i) = str2double(row_split{4});
        end

        system(sprintf('datalad drop %s/%s/stats/aseg.stats', s, ses))
        system(sprintf('datalad uninstall %s', s))
    end

    txtdir = fileparts(txtname);
    if(~exist(txtdir, 'dir'))
        mkdir(txtdir)
    end
    dlmwrite(txtname, ICV)
end

if(dohist==1)
    binwidth = 5e4;
    E = (min(ICV)-binwidth/2):binwidth:(max(ICV)+binwidth/2);
    hc = histcounts(ICV, E);
    bar(E(1:end-1) + diff(E)/2, hc');
    
    box off
    set(gcf, 'Position', [0 0 800 600])
    xloc = E(1:end-1)+diff(E)/2; 
    xloc1 = xloc(1:4:length(xloc));
    xtl = sprintfc('%.2e', xloc1);
    xloc_txt = E(1:end-1);  
    text(xloc_txt, hc+50, string(hc), 'FontSize', 13)
    set(gca, 'xtick', xloc1, 'linewidth', 2, 'fontsize', 13, 'TickDir','out')
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
    WA_ICV = ICV(WA_filter);
    AA_ICV = ICV(AA_filter);
    
    hc_WA = histcounts(WA_ICV, E);
    hc_AA = histcounts(AA_ICV, E);
    bar(E(1:end-1) + diff(E)/2, [hc_WA; hc_AA]')
    box off
    set(gcf, 'Position', [0 0 1000 600])
    set(gca, 'xtick', xloc1, 'linewidth', 2, 'fontsize', 13, 'TickDir', 'out')
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

cd(start_dir)
rmpath(genpath(fullfile(repo_dir, 'external')))

end

