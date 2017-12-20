addpath /software/gablab/conn
addpath /users/evelina9/fMRI_PROJECTS/spm_ss_vE
addpath /software/spm12
addpath /software/spm_ss

addpath('/software/spm12');
addpath('/users/evelina9/fMRI_PROJECTS/spm_ss_vE/');
addpath('/mindhive/evlab/u/bpritche/Documents/fMRI_analyses/Toolbox/spm_ss_Apr4-2016'); %set spm_ss version here!

addpath /software/spm_ss/

% defines data sources
experiments(1)=struct(...
    'name','expt1',...% language localizer 
    'pwd1','/mindhive/evlab/u/Shared/SUBJECTS',...  % path to the data directory
    'pwd2','firstlevel_langlocSN',... % path to the first-level analysis directory for the lang localizer
    'data',{{
        '430_FED_20170823a_3T2_PL2017',...
        '632_FED_20170825a_3T2_PL2017',...
        '641_FED_20170824a_3T2_PL2017'
        }}); % subject IDs
experiments(2)=struct(...
    'name','expt2',...% non-lang expt
    'pwd1','/mindhive/evlab/u/Shared/SUBJECTS',...
    'pwd2','firstlevel_AgentPatient',...  % path to the first-level analysis directory for the critical task
    'data',{{
        '430_FED_20170823a_3T2_PL2017',...
        '632_FED_20170825a_3T2_PL2017',...
        '641_FED_20170824a_3T2_PL2017'
        }});

localizer_spmfiles={};
for nsub=1:length(experiments(1).data),
    localizer_spmfiles{nsub}=fullfile(experiments(1).pwd1,experiments(1).data{nsub},experiments(1).pwd2,'SPM.mat');
end

effectofinterest_spmfiles={};
for nsub=1:length(experiments(2).data),
    effectofinterest_spmfiles{nsub}=fullfile(experiments(2).pwd1,experiments(2).data{nsub},experiments(2).pwd2,'SPM.mat');
end

partition_names={'ODD_','EVEN_'};
condition_names={'agt', 'pat'};
mask_name=''; % note: leave empty (mask_name='';) if you do not want to use subject-specific masks within each ROI
mask_thr_type='';
mask_thr_value=1;

% rois to consider are obtained from a single ROI file possibly containing multiple ROIs (e.g. name: fROI.img)
roispath='/users/evelina9/fMRI_PROJECTS/ROIS';
roiindexes=[6,4,5,2,1,3];
clear rois roinames roinumbers;
for n1=1:length(roiindexes),
    rois{n1}=fullfile(roispath,['LangParcels_n220_LH.img']);
    roinames{n1}=['ROI#',num2str(roiindexes(n1))];
    roinumbers{n1}=roiindexes(n1);
end

% options
center=0;           % between-conditions centering (patterns are normalized to a mean of zero in each voxel across all conditions)
collapserois=0;     % set to 1 to collapse across voxels from all ROIs or 0 to analyze voxels within each ROI separately
minimumvoxels=10;   % spatial correlations across ROIs (after potential subject-specific masking) with less than minimumvoxels will not be computed (minimum value minimumvoxels=2)
filenameout=[mfilename,'.csv'];
measuretype=1;      % set to 1 to compute a correlation-based similarity measure between spatial patterns of responses (Haxby et al.)
                    % set to 2 to compute an unnormalized distance-based similarity measure between spatial patterns of responses (-log(distance)) 
                    % set to 3 to compute a normalized distance-based similarity measure between spatial patterns of responses (cosine similarity)    
                    
% locates appropriate subject-level files
disp('locating appropriate subject-level files');
spm_data=[];
maskfilename={};
if ~isempty(mask_name)
    if ischar(mask_name),
       for nsubject=1:numel(experiments(1).data),
		 current_spm=fullfile(experiments(1).pwd1,experiments(1).data{nsubject},experiments(1).pwd2,'SPM.mat');
            [spm_data,SPM]=spm_ss_importspm(spm_data,current_spm);
            maskfilename{nsubject}=char(spm_ss_createlocalizermask({SPM},mask_name,[],0,mask_thr_type,mask_thr_value)); % subjects
        end
    else
        maskfilename=reshape(mask_name,[],1);
    end
end
datafilename={};
if ~isempty(experiments(2))
   for nsubject=1:numel(experiments(2).data),
	  current_spm=fullfile(experiments(2).pwd1,experiments(2).data{nsubject},experiments(2).pwd2,'SPM.mat');
        [spm_data,SPM]=spm_ss_importspm(spm_data,current_spm);
        Cnames={SPM.xCon(:).name};
        ic=[];ok=1;
        for n1=1:numel(partition_names),
            for n2=1:numel(condition_names),
                temp=strmatch([partition_names{n1},condition_names{n2}],Cnames,'exact');if numel(temp)~=1,ok=0;break;else ic(n1,n2)=temp;end;
                datafilename{nsubject,n1,n2}=fullfile(fileparts(current_spm),['con_',num2str(ic(n1,n2),'%04d'),'.img']); % subjects x partitions x conditions
            end
        end
        if ~ok, error(['contrast name ',[partition_names{n1},condition_names{n2}],' not found at ',current_spm]); end
    end
end

% Computes within-condition and between-condition spatial correlations
% of effect sizes. Bivariate correlations are computed in both cases across
% different data partitions (e.g. ODD vs. EVEN)
disp('computing spatial correlations');

if collapserois,rois={char(rois)};end
%norig=[];n0=[];r_size=[];r_within=[];r_between=[];r_all=[];
R=[];
N=[];
for nroi=1:numel(rois),
    roi=rois{nroi};
    roinumber=roinumbers{nroi};
    %roi=''; % specify an ROI file for computing within-ROI spatial correlations, or leave this empty for computing whole-brain spatial correlations
    if isempty(roi),
        roi=fullfile(fileparts(which('spm')),'apriori','brainmask.nii');
    end

    if ~isempty(maskfilename)
        mask=rex(char(maskfilename),roi,'level','voxels','select_clusters',0,'selected_clusters',roinumber,'disregard_zeros',0); % subjects x voxels mask
    end
    npartitionpair=0;
    for part1=1:numel(partition_names),
        for part2=[1:part1-1,part1+1:numel(partition_names)],
            npartitionpair=npartitionpair+1;
            Data_part1=[];Data_part2=[];
            % loads data
            for ncondition=1:numel(condition_names),
                if ~isempty(datafilename)
                    data_part1=rex(char({datafilename{:,part1,ncondition}}),roi,'level','voxels','select_clusters',0,'selected_clusters',roinumber,'disregard_zeros',0); % subjects x voxels data
                    data_part2=rex(char({datafilename{:,part2,ncondition}}),roi,'level','voxels','select_clusters',0,'selected_clusters',roinumber,'disregard_zeros',0); % subjects x voxels data
                else
                    file1=fullfile(path_group,[partition_names{part1},condition_names{ncondition}],'SPM.mat');
                    file2=fullfile(path_group,[partition_names{part2},condition_names{ncondition}],'SPM.mat');
                    if ~spm_existfile(file1), error(['file ',file1,' not found']); end
                    if ~spm_existfile(file2), error(['file ',file2,' not found']); end
                    data_part1=rex(file1,roi,'level','voxels','select_clusters',0,'selected_clusters',roinumber,'disregard_zeros',0); % subjects x voxels data
                    data_part2=rex(file2,roi,'level','voxels','select_clusters',0,'selected_clusters',roinumber,'disregard_zeros',0); % subjects x voxels data
                end
                if size(data_part1,1)~=size(data_part2,1), error(['Different number of subjects in ',file1,' and ',file2]); end
                if size(data_part1,2)~=size(data_part2,2), error(['Different number of voxels in ',file1,' and ',file2]); end
                Data_part1=cat(3,Data_part1,data_part1); % subjects x voxels x conditions data
                Data_part2=cat(3,Data_part2,data_part2); % subjects x voxels x conditions data
            end
            % eliminates voxels with non-valid data (for any subject or condition)
            nsubjects=size(Data_part1,1);
            nvoxels=size(Data_part1,2);
            % centers (across conditions)
            if center, 
                Data_part1=Data_part1-repmat(mean(Data_part1,3),[1,1,size(Data_part1,3)]);
                Data_part2=Data_part2-repmat(mean(Data_part2,3),[1,1,size(Data_part2,3)]);
            end
            
            % computes correlations
            idxvoxels=find(~any(any(isnan(Data_part1),1),3)&~any(any(isnan(Data_part2),1),3));
            for nsubject=1:nsubjects,
                if ~isempty(maskfilename)
                    idxvoxels=find(mask(nsubject,:)>0&~any(isnan(Data_part1(nsubject,:,:)),3)&~any(isnan(Data_part2(nsubject,:,:)),3));
                end
                N(nsubject,npartitionpair,nroi)=numel(idxvoxels);
                for ncondition1=1:numel(condition_names),
                    for ncondition2=1:numel(condition_names),
                        if numel(idxvoxels)>=minimumvoxels,
                            if measuretype==1, r=corrcoef([Data_part1(nsubject,idxvoxels,ncondition1)',Data_part2(nsubject,idxvoxels,ncondition2)']);r=r(1,2);      % Spatial correlations: condition x condition x subjects x partitionpairs x rois
                            elseif measuretype==2, r=mean(abs(Data_part1(nsubject,idxvoxels,ncondition1)'-Data_part2(nsubject,idxvoxels,ncondition2)').^2,1);       % Spatial distance: condition x condition x subjects x partitionpairs x rois
                            elseif measuretype==3, r=(Data_part1(nsubject,idxvoxels,ncondition1)*Data_part2(nsubject,idxvoxels,ncondition2)')/sqrt(sum(abs(Data_part1(nsubject,idxvoxels,ncondition1)').^2,1))/sqrt(sum(abs(Data_part2(nsubject,idxvoxels,ncondition2)').^2,1)); end   % Spatial distance: condition x condition x subjects x partitionpairs x rois
                            R(ncondition1,ncondition2,nsubject,npartitionpair,nroi)=r;  
                        else 
                            R(ncondition1,ncondition2,nsubject,npartitionpair,nroi)=nan;
                        end
                    end
                end
            end
        end
    end
end
if measuretype==1, Z=atanh(R); elseif measuretype==2,  Z=-log10(R); elseif measuretype==3,  Z=atanh(R); end
Z=permute(Z,[1,2,5,3,4]); % Spatial correlations (fisher transformed corr coefficients): condition(partition 1) x condition(partition 2) x rois x subjects x numberofpartitionpairs
validsubjects=permute(all(N>1,2),[3,1,2]); % rois x subjects (invalid subjects are those where the subject-specific mask has one or zero voxels within the roi) 


% Computes correct classification (pairwise-comparisons: Table 1, Haxby et al.; and multiple-comparisons)
disp('computing identification rates');

acc_pairwise=zeros([size(Z,1),size(Z,2),size(Z,3),size(Z,4)]);% pairwise-conditions accuracy (within-condition correlation higher than ARBITRARY OTHER between-condition correlation, chance level 50%) for each subject (conditions x condition x rois x subjects)
acc_all=zeros([size(Z,1),size(Z,3),size(Z,4)]);% multiple-conditions accuracy (within-condition correlation higher than ALL OTHER between-condition correlations, chance level 100%/numberofconditions) for each subject (conditions x rois x subjects)
for ncondition1=1:size(Z,1),
    for ncondition2=[1:ncondition1-1,ncondition1+1:size(Z,1)],
        acc_pairwise(ncondition1,ncondition2,:,:)= mean( Z(ncondition1,ncondition1,:,:,:)>Z(ncondition1,ncondition2,:,:,:), 5);
    end
    acc_all(ncondition1,:,:)= shiftdim(mean( all(Z(ncondition1,ncondition1*ones(1,size(Z,2)-1),:,:,:)>Z(ncondition1,[1:ncondition1-1,ncondition1+1:size(Z,2)],:,:,:),2), 5),1);
end
acc_total=permute(sum(acc_pairwise,2)/(size(acc_pairwise,2)-1),[1,3,4,2]);% total accuracy for each subject (conditions x rois x subjects)
acc_total(:,find(~validsubjects))=nan;
acc_all(:,find(~validsubjects))=nan;

fh=fopen(filenameout,'wt');

acc_n=sum(~isnan(acc_pairwise),4);
acc_mean=nanmean(acc_pairwise,4);
acc_stderr=nanstd(acc_pairwise,0,4)./sqrt(acc_n);
chance_level=.5;
acc_T=(acc_mean-chance_level)./acc_stderr;
acc_p=1-tcdf(acc_T,acc_n-1);
fprintf(fh,'%s\n',['Accuracy of pairwise identification:    mean(standard error) p-value (chance level ',num2str(chance_level*100),'%)']);
fprintf(fh,'%s','ROI,# of voxels,# of valid subjects');
for n1=1:numel(condition_names),for n2=[1:n1-1,n1+1:numel(condition_names)],fprintf(fh,'%s',[',',condition_names{n1},' > ',condition_names{n2}]);end;end
fprintf(fh,'\n');
for nroi=1:size(acc_mean,3),
    fprintf(fh,'%s',[roinames{nroi},',']);
    fprintf(fh,'%s',[num2str(mean(mean(N(:,:,nroi),2),1),'%0.0f'),' (',num2str(std(mean(N(:,:,nroi),2),0,1)/sqrt(size(N,1)),'%0.0f'),'),']);
    fprintf(fh,'%s',[num2str(mean(mean(acc_n(:,:,nroi))),'%0.0f')]);
    for n1=1:size(acc_mean,1),for n2=[1:n1-1,n1+1:size(acc_mean,2)],fprintf(fh,'%s',[',',num2str(100*acc_mean(n1,n2,nroi),'%0.0f'),'(',num2str(100*acc_stderr(n1,n2,nroi),'%0.1f'),') p=',num2str(acc_p(n1,n2,nroi),'%0.3f')]);end;end
    fprintf(fh,'\n');
end
fprintf(fh,'\n');

acc_n=sum(~isnan(acc_total),3);
acc_mean=nanmean(acc_total,3);
acc_stderr=nanstd(acc_total,0,3)./sqrt(acc_n);
chance_level=.5;
acc_T=(acc_mean-chance_level)./acc_stderr;
acc_p=1-tcdf(acc_T,acc_n-1);
fprintf(fh,'%s\n',['Accuracy of identification (average):    mean(standard error) p-value (chance level ',num2str(chance_level*100),'%)']);
fprintf(fh,'%s','ROI,# of voxels,# of valid subjects');
for n1=1:numel(condition_names),fprintf(fh,'%s',[',',condition_names{n1}, ' > * (average)']);end
fprintf(fh,'\n');
for nroi=1:size(acc_mean,2),
    fprintf(fh,'%s',[roinames{nroi},',']);
    fprintf(fh,'%s',[num2str(mean(mean(N(:,:,nroi),2),1),'%0.0f'),' (',num2str(std(mean(N(:,:,nroi),2),0,1)/sqrt(size(N,1)),'%0.0f'),'),']);
    fprintf(fh,'%s',[num2str(mean(acc_n(:,nroi)),'%0.0f')]);
    for n1=1:size(acc_mean,1),fprintf(fh,'%s',[',',num2str(100*acc_mean(n1,nroi),'%0.0f'),'(',num2str(100*acc_stderr(n1,nroi),'%0.1f'),') p=',num2str(acc_p(n1,nroi),'%0.3f')]);end
    fprintf(fh,'\n');
end
fprintf(fh,'\n');

acc_n=sum(~isnan(acc_all),3);
acc_mean=nanmean(acc_all,3);
acc_stderr=nanstd(acc_all,0,3)./sqrt(acc_n);
chance_level=1/size(Z,2);
acc_T=(acc_mean-chance_level)./acc_stderr;
acc_p=1-tcdf(acc_T,acc_n-1);
fprintf(fh,'%s\n',['Accuracy of identification (across all conditions):    mean(standard error) p-value (chance level ',num2str(chance_level*100),'%)']);
fprintf(fh,'%s','ROI,# of voxels,# of valid subjects');
for n1=1:numel(condition_names),fprintf(fh,'%s',[',',condition_names{n1},' > * (all)']);end
fprintf(fh,'\n');
for nroi=1:size(acc_mean,2),
    fprintf(fh,'%s',[roinames{nroi},',']);
    fprintf(fh,'%s',[num2str(mean(mean(N(:,:,nroi),2),1),'%0.0f'),' (',num2str(std(mean(N(:,:,nroi),2),0,1)/sqrt(size(N,1)),'%0.0f'),'),']);
    fprintf(fh,'%s',[num2str(mean(acc_n(:,nroi)),'%0.0f')]);
    for n1=1:size(acc_mean,1),fprintf(fh,'%s',[',',num2str(100*acc_mean(n1,nroi),'%0.0f'),'(',num2str(100*acc_stderr(n1,nroi),'%0.1f'),') p=',num2str(acc_p(n1,nroi),'%0.3f')]);end
    fprintf(fh,'\n');
end
fprintf(fh,'\n');


Zmean=mean(Z,5); %condition(partition 1) x condition(partition 2) x rois x subjects
fprintf(fh,'%s\n',['Within- and Between- condition measures']);
fprintf(fh,'ROI,Measure');
for nsub=1:size(Zmean,4), fprintf(fh,',Subject#%d',nsub); end
fprintf(fh,',,Average,Standard error\n');
for nroi=1:size(Zmean,3),
    for ncondition=1:numel(condition_names),
        fprintf(fh,'%s,',roinames{nroi});
        fprintf(fh,'Within_%s',condition_names{ncondition});
        t0=Zmean(ncondition,ncondition,nroi,:);
        for nsub=1:size(Zmean,4),
            fprintf(fh,',%f',t0(nsub));
        end
        t1=nanmean(t0);
        t2=nanstd(t0);
        t3=sum(~isnan(t0));
        fprintf(fh,',,%f,%f\n',t1,t2./sqrt(t3));
    end
    for ncondition1=1:numel(condition_names),
        for ncondition2=[1:ncondition1-1,ncondition1+1:numel(condition_names)],
            fprintf(fh,'%s,',roinames{nroi});
            fprintf(fh,'Between_%s_%s',condition_names{ncondition1},condition_names{ncondition2});
            t0=Zmean(ncondition1,ncondition2,nroi,:);
            for nsub=1:size(Zmean,4),
                fprintf(fh,',%f',t0(nsub));
            end
            t1=nanmean(t0);
            t2=nanstd(t0);
            t3=sum(~isnan(t0));
            fprintf(fh,',,%f,%f\n',t1,t2./sqrt(t3));
        end
    end
end
fprintf(fh,'\n');

fprintf(fh,'%s\n',['Within- vs. Between- condition comparisons']);
fprintf(fh,'ROI,Measure');
for nsub=1:size(Zmean,4), fprintf(fh,',Subject#%d',nsub); end
fprintf(fh,',,Average,Standard error,T/F,dof,p\n');
for nroi=1:size(Zmean,3),
    for ncondition=1:numel(condition_names),
        fprintf(fh,'%s,',roinames{nroi});
        fprintf(fh,'Within_%s>0',condition_names{ncondition});
        t0=Zmean(ncondition,ncondition,nroi,:);
        for nsub=1:size(Zmean,4),
            fprintf(fh,',%f',t0(nsub));
        end
        t1=nanmean(t0);
        t2=nanstd(t0);
        t3=sum(~isnan(t0));
        t4=1-spm_Tcdf(t1./t2.*sqrt(t3),t3-1);
        fprintf(fh,',,%f,%f,%f,%d,%f\n',t1,t2./sqrt(t3),t1./t2.*sqrt(t3),t3-1,t4);
    end
    for ncondition1=1:numel(condition_names),
        for ncondition2=[1:ncondition1-1,ncondition1+1:numel(condition_names)],
            fprintf(fh,'%s,',roinames{nroi});
            fprintf(fh,'Within_%s>Between_%s_%s',condition_names{ncondition1},condition_names{ncondition1},condition_names{ncondition2});
            t0=Zmean(ncondition1,ncondition1,nroi,:)-Zmean(ncondition1,ncondition2,nroi,:);
            for nsub=1:size(Zmean,4),
                fprintf(fh,',%f',t0(nsub));
            end
            t1=nanmean(t0);
            t2=nanstd(t0);
            t3=sum(~isnan(t0));
            t4=1-spm_Tcdf(t1./t2.*sqrt(t3),t3-1);
            fprintf(fh,',,%f,%f,%f,%d,%f\n',t1,t2./sqrt(t3),t1./t2.*sqrt(t3),t3-1,t4);
        end
    end
    for ncondition1=1:numel(condition_names),
        for ncondition2=[ncondition1+1:numel(condition_names)],
            fprintf(fh,'%s,',roinames{nroi});
            fprintf(fh,'%s~=%s',condition_names{ncondition1},condition_names{ncondition2});
            t0=[squeeze(Zmean(ncondition1,ncondition1,nroi,:)-Zmean(ncondition1,ncondition2,nroi,:)),squeeze(Zmean(ncondition2,ncondition2,nroi,:)-Zmean(ncondition2,ncondition1,nroi,:))];
            for nsub=1:size(Zmean,4),
                fprintf(fh,',[%f;%f]',t0(nsub,1),t0(nsub,2));
            end
            t1=nanmean(t0);
            t2=nanstd(t0);
            t3=sum(~isnan(t0),1);
            t0(any(isnan(t0),2),:)=[];
            [h,f,p,dof]=glm(ones(size(t0,1),1),t0);
            fprintf(fh,',,[%f;%f],[%f;%f],%f,(%d;%d),%f\n',t1(1),t1(2),t2(1)./sqrt(t3(1)),t2(2)./sqrt(t3(2)),f,dof(1),dof(2),p);
        end
    end
end
fprintf(fh,'\n');

fclose(fh);
disp(['Results stored in ',filenameout,'.']);


% Creates plots (figure 4, Haxby et al.)
disp('displaying plots');

Z_n=sum(~isnan(mean(Z,5)),4);
Z_mean=nanmean(mean(Z,5),4);
Z_stderr=nanstd(mean(Z,5),0,4)./sqrt(Z_n); 
nrois=1:size(Z_mean,3);

figure('color','w');
for ncondition=1:size(Z_mean,1), % conditions
    subplot(size(Z_mean,1)+1,1,ncondition);hold on;
    z=shiftdim(Z_mean(ncondition,:,nrois),1);
    e=shiftdim(Z_stderr(ncondition,:,nrois),1);
    h1=bar(z);
    ztemp=z;ztemp(ncondition,:)=0;
    h2=bar(ztemp);
    xdata=[];
    for nroi=1:numel(h1), set(h1(nroi),'facecolor',[1,0,0]*nroi/numel(h1)); temp=get(get(h1(nroi),'children'),'xdata'); xdata(:,nroi)=mean(temp,1)'; end
    for nroi=1:numel(h2), set(h2(nroi),'facecolor',[0,0,1]*nroi/numel(h2)); end
    errorbar(xdata,z,e,'k.');
    hold off;
    axis tight;
    set(gca,'xlim',[0,size(Z_mean,2)+1],'xcolor','w');
    ylabel(condition_names{ncondition});
end
set(gca,'xcolor','k','xtick',1:numel(condition_names),'xticklabel',condition_names);
subplot(size(Z_mean,1)+1,3,size(Z_mean,1)*3+2);
c0=permute([1,0,0;0,0,1],[1,3,2]);c=[];
for nroi=1:size(Z_mean,3),c(:,nroi,:)=c0(:,1,:)*nroi/size(Z_mean,3); end
image(c);axis equal; axis tight;
ht=text(1:numel(roinames),size(c,1)+ones(1,numel(roinames)),roinames);set(ht,'rotation',-90);
if measuretype==1, set(gca,'xtick',[],'ytick',1:2,'yticklabel',{'Within-condition correlations','Between-condition correlations'},'yaxislocation','right');
elseif measuretype==2, set(gca,'xtick',[],'ytick',1:2,'yticklabel',{'Within-condition similarity','Between-condition similarity'},'yaxislocation','right'); 
elseif measuretype==3, set(gca,'xtick',[],'ytick',1:2,'yticklabel',{'Within-condition similarity','Between-condition similarity'},'yaxislocation','right'); end

saveas(gca, [mfilename '_corrs'], 'jpg');
