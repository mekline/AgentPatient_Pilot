

outputdir = '/mindhive/evlab/u/mekline/Documents/Projects/AgentPatient_Pilot/langloc_resp_ToM_20170317';

experiments(1)=struct(...
    'name','loc',...% language localizer 
    'pwd1','/mindhive/evlab/u/Shared/SUBJECTS',...  % path to the data directory
    'pwd2','firstlevel_ToMshort',... % path to the first-level analysis directory for the lang localizer
    'data',{{ ...
        '301_FED_20150708c_3T2',...
        '473_FED_20170210b_3T2',...
        '498_FED_20170210c_3T2'
        }}); % subject IDs
experiments(2)= struct(...
    'name','crit',...
    'pwd1','/mindhive/evlab/u/Shared/SUBJECTS',... 
    'pwd2','firstlevel_AgPat',...
    'data',{{ ...
        '301_FED_20161217b_3T2',...
        '473_FED_20170210b_3T2',...
        '498_FED_20170210c_3T2'
        }});% important that this is the same order as above

localizer_spmfiles={};
for nsub=1:length(experiments(1).data),
    localizer_spmfiles{nsub}=fullfile(experiments(1).pwd1,experiments(1).data{nsub},experiments(1).pwd2,'SPM.mat');
end

effectofinterest_spmfiles={};
for nsub=1:length(experiments(2).data),
    effectofinterest_spmfiles{nsub}=fullfile(experiments(2).pwd1,experiments(2).data{nsub},experiments(2).pwd2,'SPM.mat');
end

ss=struct(...
    'swd', outputdir,...   % output directory
    'EffectOfInterest_spm',{effectofinterest_spmfiles},...
    'Localizer_spm',{localizer_spmfiles},...
    'EffectOfInterest_contrasts',{{'agt','pat'}},... % contrasts of interest; here you just indicate the names of the four conditions
    'Localizer_contrasts',{{'bel-pho'}},...                     % localizer contrast 
    'Localizer_thr_type','percentile-ROI-level',...  % appendix F, paper by Alfonso & Ev, 2012, has all the options
    'Localizer_thr_p',0.1,... % proportion, not actually significance level (what proportion of ROI will be used, IF localizer_thr_type is percentile-ROI-level
    'type','mROI',...                                       % can be 'GcSS', 'mROI', or 'voxel'
    'ManualROIs','/users/evelina9/fMRI_PROJECTS/ROIS/ToMparcels_Mar2015.img',... % predefined, this is the one we usually use, ask beforehand
    'model',1,...    % can be 1 (one-sample t-test), 2 (two-sample t-test), or 3 (multiple regression), usually 1
    'estimation','OLS',... % basically always use, ordinary least squares, as opposed to ReML
    'overwrite',true,... % clears stuff from earlier toolbox analyses
    'ask','missing');    % can be 'none' (any missing information is assumed to take default values), 'missing' (any missing information will be asked to the user), 'all' (it will ask for confirmation on each parameter)

addpath('/mindhive/evlab/u/bpritche/Documents/fMRI_analyses/Toolbox/spm_ss_Apr4-2016'); %set spm_ss version here!

ss=spm_ss_design(ss);    % see help spm_ss_design for additional information
ss=spm_ss_estimate(ss);  % can be found in /software/spm_ss

% spm_ss_display and spm_ss_results allow you to interface with your
% results and visualize things
% spm_ss_watershed lets you visualize that probabilistic hilly space
