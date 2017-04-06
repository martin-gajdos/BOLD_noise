function BOLD_noise
% This is main function of the toolbox. The tool is designed to create
% simulated noise into BOLD data based on properties observed in real
% datasets. The simulated noise consists particularly from:
% 1. Noise caused by movements of subject's head
% 2. Noise caused by breathing and heartbeats
% 3. 1/f noise
%
% Output specification:
% 
% Output:
% 
% Toolbox requires SPM12
% Subpath source_data required
% -------------------------------------------------------------------------
% Version: 01
% Last update: 15.12.2016
% -------------------------------------------------------------------------
% Written by Martin Gajdos
%
% Supported by the Czech Science Foundation Grant Project No.14-33143S
% -------------------------------------------------------------------------

%% Specification of the output
no_of_scans=300;
TR=2;
breathing_const=1;
heartbeat_const=1;
movement_const=1;
vx_size=[3 3 3];



%% 

options.no_of_scans=no_of_scans;
options.TR=TR;
options.breathing_const=breathing_const;
options.heartbeat_const=heartbeat_const;
options.movement_const=movement_const;
options.vx_size=vx_size;
options.plot_movement_regressors=0;

main_fcn_path=mfilename('fullpath');
[toolbox_path,~]=fileparts(main_fcn_path);

output_path=toolbox_path;
output_name='test_sumu2.nii';

options.movement_maps_path=fullfile(toolbox_path,'source_data','movement_simulation');
options.movement_maps_mask_name='Mask.img';
options.movement_maps_mu_names={ % Order has to fit the order of variables in reg
    'movement_maps_mu_translationX.img';...
    'movement_maps_mu_translationY.img';...
    'movement_maps_mu_translationZ.img';...
    'movement_maps_mu_rotationX.img';...
    'movement_maps_mu_rotationY.img';...
    'movement_maps_mu_rotationZ.img';...
    'movement_maps_mu_translationX_squared.img';...
    'movement_maps_mu_translationY_squared.img';...
    'movement_maps_mu_translationZ_squared.img';...
    'movement_maps_mu_rotationX_squared.img';...
    'movement_maps_mu_rotationY_squared.img';...
    'movement_maps_mu_rotationZ_squared.img';...
    'movement_maps_mu_translationX_diference.img';...
    'movement_maps_mu_translationY_diference.img';...
    'movement_maps_mu_translationZ_diference.img';...
    'movement_maps_mu_rotationX_diference.img';...
    'movement_maps_mu_rotationY_diference.img';...
    'movement_maps_mu_rotationZ_diference.img';...
    'movement_maps_mu_translationX_diference_squared.img';...
    'movement_maps_mu_translationY_diference_squared.img';...
    'movement_maps_mu_translationZ_diference_squared.img';...
    'movement_maps_mu_rotationX_diference_squared.img';...
    'movement_maps_mu_rotationY_diference_squared.img';...
    'movement_maps_mu_rotationZ_diference_squared.img';...
    };
options.physio_path=fullfile(toolbox_path,'source_data','physio_simulation');

% options.inbrain_area_mask=fullfile(spm('dir'),'tpm','mask_ICV.nii'); % Intracranial mask
%% Movement noise simulation
[M_noise,movement_parameters]=get_movement_noise(options);

options.M_noise=M_noise;
% 
%% Breathing and heartbeats simulation
[P_noise,physio_parameters]=get_physio_noise(options);

%% 1/f noise
F_noise=get_1_f_noise(options);

BOLD_noise=(1/4).*M_noise+(2/4).*P_noise+(1/4).*F_noise;

%% Writing result

fname=fullfile(output_path,output_name);
% 
spm_unlink(fname)

vol=movement_parameters.mask_vol;
dt=spm_type('float32');
description='BOLD noise data';

write_4D_nifti_MG(BOLD_noise,vol,fname,dt,description);
% dim=size(BOLD_noise);
% dim=dim(1:3);
% for i=1:no_of_scans
%     BOLD_noise_struct(i) = deal(struct(...
% 		'fname',	[],...
% 		'dim',		dim,...
%         'dt',       [spm_type('int16') 0],...
% 		'mat',		movement_parameters.mask_vol.mat,...
% 		'pinfo',	[1 0 0]',...
% 		'descrip',	'BOLD noise data',...
%         'n',        [],...
%         'private',  []));
%     BOLD_noise_struct(i).fname   = fname;
%     BOLD_noise_struct(i) = spm_create_vol(BOLD_noise_struct(i));
%     spm_write_vol(BOLD_noise_struct(i),squeeze(BOLD_noise(:,:,:,i)));
% end

disp(' ')
disp('---------------------------')
disp(['BOLD noise file ',output_name,' successfully saved in ',output_path])
disp('---------------------------')
disp(' ')


function [M_noise,movement_parameters]=get_movement_noise(options)
% Function simulates movement of subject's head and reflection of
% this movement into the data. Simulation of the head movement is based on
% observed parameters. Reflection of movement in the MNI space is based on
% observed explained variability of the movement regressors in used source
% datasets.

%% Simulation of movement regressors
% Normal distribution of translations of head movement (differences)
mu_trX = -1.8058e-04;
mu_trY = 8.5625e-05;
mu_trZ = 0.0018;
sigma_trX = 0.0214;
sigma_trY = 0.1069;
sigma_trZ =0.0542;

% Normal distribution of rotations of head movement (differences)
mu_rotX = 1.4494e-05;
mu_rotY = 1.7250e-06;
mu_rotZ = -1.0947e-05;
sigma_rotX = 0.0010;
sigma_rotY = 4.0645e-04;
sigma_rotZ = 4.5307e-04;

% Limits for head positions (absolute values)
lim_min_trX=-1.5308;
lim_min_trY=-2.8109;
lim_min_trZ=-1.5869;
lim_max_trX=1.1740;
lim_max_trY=1.2805;
lim_max_trZ=7.4549;

lim_min_rotX=-0.0236;
lim_min_rotY=-0.0321;
lim_min_rotZ=-0.0272;
lim_max_rotX=0.0498;
lim_max_rotY=0.0490;
lim_max_rotZ=0.0477;

nscan=options.no_of_scans;
% Initial head position
translation_x=NaN(nscan,1);
translation_y=NaN(nscan,1);
translation_z=NaN(nscan,1);

rotation_x=NaN(nscan,1);
rotation_y=NaN(nscan,1);
rotation_z=NaN(nscan,1);

translation_x(1)=0;
translation_y(1)=0;
translation_z(1)=0;

rotation_x(1)=0;
rotation_y(1)=0;
rotation_z(1)=0;


for i=1:nscan-1
    
    diff_tr_x=random('Normal',mu_trX,sigma_trX);
    diff_tr_y=random('Normal',mu_trY,sigma_trY);
    diff_tr_z=random('Normal',mu_trZ,sigma_trZ);
    
    diff_rot_x=random('Normal',mu_rotX,sigma_rotX);
    diff_rot_y=random('Normal',mu_rotY,sigma_rotY);
    diff_rot_z=random('Normal',mu_rotZ,sigma_rotZ);
    
    translation_x(i+1) = translation_x(i)+diff_tr_x;
    translation_y(i+1) = translation_y(i)+diff_tr_y;
    translation_z(i+1) = translation_z(i)+diff_tr_z;
    
    rotation_x(i+1) = rotation_x(i)+diff_rot_x;
    rotation_y(i+1) = rotation_y(i)+diff_rot_y;
    rotation_z(i+1) = rotation_z(i)+diff_rot_z;
    
    if translation_x(i+1)<= lim_min_trX || translation_x(i+1)>= lim_max_trX
        translation_x(i+1) = translation_x(i)-diff_tr_x;
    end
    
    if translation_y(i+1)<= lim_min_trY || translation_y(i+1)>= lim_max_trY
        translation_y(i+1) = translation_y(i)-diff_tr_y;
    end
    
    if translation_z(i+1)<= lim_min_trZ || translation_z(i+1)>= lim_max_trZ
        translation_z(i+1) = translation_z(i)-diff_tr_z;
    end
    
    if rotation_x(i+1)<= lim_min_rotX || rotation_x(i+1)>= lim_max_rotX
        rotation_x(i+1) = rotation_x(i)-diff_rot_x;
    end
    
    if rotation_y(i+1)<= lim_min_rotY || rotation_y(i+1)>= lim_max_rotY
        rotation_y(i+1) = rotation_y(i)-diff_rot_y;
    end
    
    if rotation_z(i+1)<= lim_min_rotZ || rotation_z(i+1)>= lim_max_rotZ
        rotation_z(i+1) = rotation_z(i)-diff_rot_z;
    end
    
      
end

rp=[translation_x translation_y translation_z rotation_x rotation_y rotation_z];

rp2=rp.^2;

rp_diff=[zeros(1,size(rp,2));diff(rp)]; % First row are zeros, differences are following

rp_diff2=rp_diff.^2;

reg=[rp rp2 rp_diff rp_diff2];

if options.plot_movement_regressors==1
    figure
    plot(translation_x,'r')
    hold on
    plot(translation_y,'g')
    plot(translation_z,'b')
    legend('translace x','translace y','translace z')
    hold off
    
    figure
    plot(rotation_x,'r')
    hold on
    plot(rotation_y,'g')
    plot(rotation_z,'b')
    legend('rotace x', 'rotace y','rotace z')
    hold off
end

%% Simulation of 4D data
movement_maps_path=options.movement_maps_path;
movement_maps_mask_name=options.movement_maps_mask_name;
movement_maps_mu_names=options.movement_maps_mu_names;

Nmaps=size(movement_maps_mu_names,1);

mask_vol=spm_vol([fullfile(movement_maps_path,movement_maps_mask_name),',1']);
[mask_Y,mask_XYZ] = spm_read_vols(mask_vol);

val_coords_id=find(mask_Y);
Ncoords=numel(val_coords_id);
coords_mm=mask_XYZ(:,val_coords_id);
coords_mm=[coords_mm; ones(1,size(coords_mm,2))];
coords_vx=mask_vol.mat\coords_mm;
coords_vx=coords_vx(1:3,:);


M_noise=zeros(mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3),nscan);

rand_val=cell(Nmaps,1);
NaN(mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3),nscan);
for ii=1:Nmaps
    rand_val{ii,1}=zeros(mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3));
    
    tmp_mu_map_name=fullfile(movement_maps_path,movement_maps_mu_names{ii});
    tmp_mu_vol=spm_vol(tmp_mu_map_name);
    
    [tmp_mu,tmp_XYZ] = spm_read_vols(tmp_mu_vol);
    
    for jj=1:Ncoords
        tmp_vx=coords_vx(:,jj);
        tmp_mu_val=tmp_mu(tmp_vx(1),tmp_vx(2),tmp_vx(3)); % Doplnit linearni zmenseni v pomeru k TR

        rand_val{ii}(tmp_vx(1),tmp_vx(2),tmp_vx(3))=random('Exponential',tmp_mu_val);  % Statistic toolbox

    end
    
    tmp_reg=squeeze(reg(:,ii));
    tmp_reg=reshape(tmp_reg, [1, 1, 1, numel(tmp_reg)]);
    tmp_reg=repmat(tmp_reg, [mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3),1]);
    
    rand_val_rep=repmat(rand_val{ii},[1,1,1,nscan]);
    
    tmp_reg2=rand_val_rep.*tmp_reg;
    M_noise=M_noise+tmp_reg2;
  
end
M_noise=M_noise-mean(M_noise(:));
M_noise=M_noise./max(abs(M_noise(:)));

%%
% % for ii=1:Nmaps
% %     tmp_mu_map_name=fullfile(movement_maps_path,movement_maps_mu_names{ii});
% %     tmp_mu_vol=spm_vol(tmp_mu_map_name);
% %     
% %     [tmp_mu,tmp_XYZ] = spm_read_vols(tmp_mu_vol);
% %     
% %     for jj=1:Ncoords
% %         tmp_vx=coords_vx(:,jj);
% %         tmp_mu_val=tmp_mu(tmp_vx(1),tmp_vx(2),tmp_vx(3));
% % 
% %         rand_val=random('Exponential',tmp_mu_val);  % Statistic toolbox
% %         tmp_reg=rand_val*squeeze(reg(:,ii));
% %         tmp_reg=reshape(tmp_reg, [1, 1, 1, numel(tmp_reg)]);
% %         M_noise(tmp_vx(1),tmp_vx(2),tmp_vx(3),:)=M_noise(tmp_vx(1),tmp_vx(2),tmp_vx(3),:)+tmp_reg;
% %     end
% %     
% % end

movement_parameters.regressors=reg;
movement_parameters.map_names=movement_maps_mu_names;
movement_parameters.movement_maps_path=movement_maps_path;
movement_parameters.movement_maps_mask_name=movement_maps_mask_name;
movement_parameters.movement_maps_mu_names=movement_maps_mu_names;
movement_parameters.mask_vol=mask_vol;

disp(' ')
disp('Simulation of movement noise completed')
disp(' ')

function [P_noise,physio_parameters]=get_physio_noise(options)
% Function simulates breathing and heartbeats of subject and its reflection 
% into the data. It is based on observed impulse responses caused by breath
% and heart artifact, fitted by RETROICOR. The noise keeps also correlated
% structure in MNI space as observed in source datasets.

   
movement_maps_mask_name=fullfile(options.movement_maps_path,options.movement_maps_mask_name);
sim_mask_vol=spm_vol(movement_maps_mask_name);

no_of_scans=options.no_of_scans;    
TR=options.TR;
sig_length=no_of_scans*TR;

% Heartbeats
HB_T_mean = 0.8857;
HB_T_std = 0.1548;

% Respiration
Res_T_mean=3.7564;
Res_T_std=1.0256;

% Load basis function set, variable zset <1000 x 8>
load(fullfile(options.physio_path,'RETROICOR_zset.mat'))
zset_samples=size(zset,1);
% Load distribution of correlation structure
load(fullfile(options.physio_path,'breath_EKG_corr_distributions.mat'))

%% Heartbeats timecourse simulation
cond=1;
HB_T=[];

while cond==1
    HB_T(end+1,1)=HB_T_mean + HB_T_std*randn(1, 1);
    tmp=sum(HB_T);
    if sum(tmp>sig_length)
        cond=0;
    end
end

N_HB=numel(HB_T);
HB_t=NaN(numel(HB_T)+1,1); % Timecourse of onsets
HB_t(1)=abs(randn(1, 1).*HB_T_mean); % Random onset

for i=1:N_HB
    HB_t(i+1)=HB_t(i)+HB_T(i);
end

%% Betas for heartbeat impulse response
HB_beta_mean=[7928.48953290513 -12.8930714614321 -2.29230226396097 -2.54374938667226 -1.43773986459029 -1.75448400402841 -0.867089201423681 -2.94034030668576];
HB_beta_std=[3432.19301777188 163.854619079308 87.5674018078018 49.1680781869313 36.4988185002547 28.5509422602397 24.5989947188217 30.6070281543519];

HB_beta=NaN(8,1);

HB_beta(1)=abs(HB_beta_mean(1) + HB_beta_std(1)*randn(1,1)); % First beta is always positive
for i=2:8
    HB_beta(i)=HB_beta_mean(i) + HB_beta_std(i)*randn(1,1);
end
%% Heartbeat impulse response
HB_imp=zset*HB_beta;
HB_imp_t=HB_T_mean; % Length of impulse response
%% Heartbeat artifact timecourse
dt=HB_imp_t/zset_samples;
N_fine=round(sig_length/dt);
HB_onsets=zeros(N_fine,1);
for i=1:numel(HB_t)
    tmp_t=HB_t(i);
    tmp_N=round(tmp_t/dt);
    HB_onsets(tmp_N)=1;
end
N_fine=numel(HB_onsets);
dt2=round(N_fine/no_of_scans);

HB_art=conv(HB_imp,HB_onsets);
HB_art_down= downsample(HB_art,dt2);% Downsampled to TR, signal processing toolbox
HB_art_down=HB_art_down(1:no_of_scans);
HB_art_down=HB_art_down./max(abs(HB_art_down));
%% Heartbeat correlation structure
n_HB_pos=size(EKG_corr_mean,1);
HB_corr=ones(size(EKG_corr_mean));
for i=1:n_HB_pos
    for j=i:n_HB_pos
        if i~=j
            tmp=EKG_corr_mean(i,j) + EKG_corr_sdt(i,j)*randn(1,1);
            tmp=check_corr(tmp);
            HB_corr(i,j)=tmp;
            HB_corr(j,i)=tmp;
        end
    end
end

HB_sig=zeros(sim_mask_vol.dim(1),sim_mask_vol.dim(2),sim_mask_vol.dim(3),no_of_scans);

HB_pos_mm=EKG_corr_pos.mm;
HB_pos_mmx1=[HB_pos_mm; ones(1,size(HB_pos_mm,2))];
HB_vx=sim_mask_vol.mat\HB_pos_mmx1;
HB_vx=HB_vx(1:3,:);


%% Heartbeat 4D data
for i=1:n_HB_pos
    tmp_vx=round(HB_vx(:,i));
    corr_vect=HB_corr(i,:);
    corr_vect(i)=[];% Diagonal
    corr_vect(abs(corr_vect)<0.2)=[];
    tmp_corr=median(corr_vect);
    tmp_const=(1/tmp_corr)-1;
    tmp_rnd_sig=tmp_const*rand(no_of_scans,1);
    tmp_HB_sig=sign(tmp_corr).*HB_art_down+tmp_rnd_sig;
    tmp_HB_sig=tmp_HB_sig-mean(tmp_HB_sig);
    tmp_HB_sig=tmp_HB_sig./max(abs(tmp_HB_sig));
    HB_sig(tmp_vx(1),tmp_vx(2),tmp_vx(3),:)=tmp_HB_sig;
end
heart_noise=HB_sig;
heart_noise_const=0.5;

%% Generating respiration timecourse
cond=1;
Res_T=[];

while cond==1
    Res_T(end+1,1)=Res_T_mean + Res_T_std*randn(1,1);
    tmp=sum(Res_T);
    if sum(tmp>sig_length)
        cond=0;
    end
end

N_Res=numel(Res_T);
Res_t=NaN(numel(Res_T)+1,1); % Timecourse of onsets
Res_t(1)=abs(randn(1, 1).*Res_T_mean); % Random onset

for i=1:N_Res
    Res_t(i+1)=Res_t(i)+Res_T(i);
end

%% Betas for respiration impulse response
    
Res_beta_mean=[4.65104544833907 0.0692251295831512 -0.552516858219035 7.41935858988928 0.00135336583198837 -2.97490770527482 -0.313859699015594 -1.00166632428469];
Res_beta_std=[121.811152994733 27.6641874709205 22.2595038578895 187.817171503323 29.1896972144084 55.0845404754985 16.9292545023538 36.0814526698351];
Res_beta=NaN(8,1);

for i=1:8
    Res_beta(i)=Res_beta_mean(i) + Res_beta_std(i)*randn(1,1);
end


%% Respiration impulse response
Res_imp=zset*Res_beta;
Res_imp_t=Res_T_mean; % Length of impulse response
%% Respiration artifact timecourse
dt=Res_imp_t/zset_samples;
N_fine=round(sig_length/dt);
Res_onsets=zeros(N_fine,1);
for i=1:numel(Res_t)
    tmp_t=Res_t(i);
    tmp_N=round(tmp_t/dt);
    Res_onsets(tmp_N)=1;
end
N_fine=numel(Res_onsets);
dt2=round(N_fine/no_of_scans);

Res_art=conv(Res_imp,Res_onsets);
Res_art_down= downsample(Res_art,dt2);% Downsampled to TR, signal processing toolbox
Res_art_down=Res_art_down(1:no_of_scans);
Res_art_down=Res_art_down./max(abs(Res_art_down));
%% Respiration correlation structure
n_Res_pos=size(breath_corr_mean,1);
Res_corr=ones(size(breath_corr_mean));
for i=1:n_Res_pos
    for j=i:n_Res_pos
        if i~=j
            tmp=breath_corr_mean(i,j) + breath_corr_sdt(i,j)*randn(1,1);
            tmp=check_corr(tmp);
            Res_corr(i,j)=tmp;
            Res_corr(j,i)=tmp;
        end
    end
end

Res_sig=zeros(sim_mask_vol.dim(1),sim_mask_vol.dim(2),sim_mask_vol.dim(3),no_of_scans);

Res_pos_mm=breath_corr_pos.mm;
Res_pos_mmx1=[Res_pos_mm; ones(1,size(Res_pos_mm,2))];
Res_vx=sim_mask_vol.mat\Res_pos_mmx1;
Res_vx=Res_vx(1:3,:);


%% Respiration 4D data
for i=1:n_Res_pos
    tmp_vx=round(Res_vx(:,i));
    corr_vect=Res_corr(i,:);
    corr_vect(i)=[];% Diagonal
    corr_vect(abs(corr_vect)<0.2)=[];
    tmp_corr=median(corr_vect);
    tmp_const=(1/tmp_corr)-1;
    tmp_rnd_sig=tmp_const*rand(no_of_scans,1);
    tmp_Res_sig=sign(tmp_corr).*Res_art_down+tmp_rnd_sig;
    tmp_Res_sig=tmp_Res_sig-mean(tmp_Res_sig);
    tmp_Res_sig=tmp_Res_sig./max(abs(tmp_Res_sig));
    Res_sig(tmp_vx(1),tmp_vx(2),tmp_vx(3),:)=tmp_Res_sig;
end
breath_noise=Res_sig;

breath_noise_const=0.5;

P_noise=heart_noise_const*heart_noise+breath_noise_const*breath_noise;

physio_parameters.respiration_artifact=Res_art_down;
physio_parameters.heartbeat_artifact=HB_art_down;
physio_parameters.respiration_onsets=Res_onsets;
physio_parameters.heartbeat_onsets=HB_onsets;

disp(' ')
disp('Simulation of physiological noise completed')
disp(' ')



function corr=check_corr(corr)
% Function ensures, that generated correlations are in <-1 to 1>

if corr>1
    corr=0.5 + 0.01*randn(1,1);
%     corr=check_corr(corr);
end
        
if corr<-1
    corr=-0.5 + 0.01*randn(1,1);
%     corr=check_corr(corr);
end     
        

function F_noise=get_1_f_noise(options)
% Function simulates 1/f noise. It is based on noise observed in source
% datasets.

% Gamma distribution for hyperbola shape parameter
c_shape=1.7115;
c_scale=0.0008;


% Normal distribution for hyperbola shift in x axis
m_mean= -0.0032;
m_std=0.0074;


% Gamma distribution for hyperbola shift in y axis
n_shape=3.9462;
n_scale=0.0153;


%% Simulation of 4D data
nscan=options.no_of_scans;
mask_path=options.movement_maps_path;
mask_name=options.movement_maps_mask_name;
mask_vol=spm_vol([fullfile(mask_path,mask_name),',1']);
[mask_Y,mask_XYZ] = spm_read_vols(mask_vol);

val_coords_id=find(mask_Y);
Ncoords=numel(val_coords_id);
coords_mm=mask_XYZ(:,val_coords_id);
coords_mm=[coords_mm; ones(1,size(coords_mm,2))];
coords_vx=mask_vol.mat\coords_mm;
coords_vx=coords_vx(1:3,:);

TR=options.TR;

F_noise=zeros(mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3),nscan);

% xdata=1:nscan;
for jj=1:Ncoords
    tmp_vx=coords_vx(:,jj);
    %% hypebola preparation
    cond=0;
    while cond==0
        c = spm_gamrnd(c_shape, c_scale,1);
        
        m = m_mean + m_std*randn(1, 1);
        m = abs(m); % Subzero values are leading to unsuitable hyperbolas
        
        n = spm_gamrnd(n_shape, n_scale, 1);
        
%         y = c*(1./(xdata-m))+n;% Hyperbola
        % Random signal
        rand_sig=randn(1,nscan); % Normally distributed white noise
        % Spectra of random signal
        [spectrum_noise_tmp,spectrum_noise_phase,axis_label]=get_spectrum(rand_sig,TR);
        % Rescale amplitude spectrum to get 1 mean amplitude for every frequency
        spectrum_noise_tmp=spectrum_noise_tmp-mean(spectrum_noise_tmp)+1;
        
        y = c*(1./(axis_label+m))+n; % Only positive values of m leads to suitable hyperbolas - negative m leads to hyperbloas with negative values at (0,Inf)
        if sum(y>1.8)==0 && sum(y<0)==0 % Eliminate to high and negative hyperbolas
            cond=1;
        end
    end
    
    %% 1/f spectra
    % Spectrum is changed in first half of amplitude spectrum (0:Fs/s)
    fft_noise=spectrum_noise_tmp.*y;
    
    % Reconstruction of amplitude and phase spectra to full length
    N=nscan;
    if ~mod(N,2) % even length of signal - influences reconstruction of spectra
        spectrum_all=[fft_noise fft_noise(end-1:-1:2)];
        spectrum_ph_all=[spectrum_noise_phase -1*spectrum_noise_phase(end-1:-1:2)];
    else % uneven length of signal - influences reconstruction of spectra
        spectrum_all=[fft_noise fft_noise(end:-1:2)];
        spectrum_ph_all=[spectrum_noise_phase -1*spectrum_noise_phase(end:-1:2)];
    end
    % Conversion of amplitude and phase to imaginar and real numbers - necessary for Matlab ifft
    I=sin(spectrum_ph_all).*spectrum_all;
    R=cos(spectrum_ph_all).*spectrum_all;
    
    tmp_sig=N*real(ifft( R+I*1i ));
    
    tmp_sig=tmp_sig-mean(tmp_sig);
    tmp_sig=tmp_sig./max(abs(tmp_sig));
    
    F_noise(tmp_vx(1),tmp_vx(2),tmp_vx(3),:)=tmp_sig;
    
end

disp(' ')
disp('Simulation of 1/f noise completed')
disp(' ')


function [spectrum,spectrum_ph,axis_label]=get_spectrum(signal,TR)
% Function computes amplitude spectrum and creates label for this spectrum
% Variables:
% signal - time domain signal
% TR - reptition time [s]
% spectrum - amplitude spectrum 0:Fs/2; amplitudes are scalled to be 
% approximately equal to 'a' in equation y=a*sin(2*pi*f*t+phi); 
% in fft function
% axis_label- 0:Fs/2
%
% Written by MG
%
% Last update 8.1.2014 16:51
% 21.5.2014 Checked scaling and involved phase spectrum - to be able fully
% reconstruct signal, using complex fft representation in Matlab

Fs=1/TR; % Sampling frequency [Hz]

N=numel(signal);

Y = fft(signal)/N;
axis_label = Fs/2*linspace(0,1,length(signal)/2+1);
% % single-sided amplitude spectrum - source: Matlab fft help
spectrum=2*abs(Y(1:floor(length(signal)/2+1))); % Spektrum, jak je použito ve výpoètech
spectrum_ph=angle(Y(1:floor(length(signal)/2+1)));

% % %%%
% % 
% % 
% % %%
% % % function [spectrum,axis_label]=get_spectrum(signal,TR)
% % % % Function computes amplitude spectrum and creates label for this spectrum
% % % % Variables:
% % % % signal - time domain signal
% % % % TR - reptition time [s]
% % % % spectrum - amplitude spectrum 0:Fs/2; amplitudes are scalled to be 
% % % % approximately equal to 'a' in equation y=a*sin(2*pi*f*t+phi); length of
% % % % spectrum is not equal to 1/2 length of signal, because of using extension
% % % % of signal by NFFT (power of 2) in fft function
% % % % axis_label- 0:Fs/2
% % % %
% % % % Written by MG 17.12.2013
% % % %
% % % % Last update 8.1.2014 16:51
% % % 
% % % Fs=1/TR; % Sampling frequency [Hz]
% % % 
% % % N=numel(signal);
% % % 
% % % Y = fft(signal)/N;
% % % axis_label = Fs/2*linspace(0,1,length(signal)/2+1);
% % % % % single-sided amplitude spectrum - source: Matlab fft help
% % % spectrum=2*abs(Y(1:length(signal)/2+1));


function V4 = write_4D_nifti_MG(data,data_vol,fname,dt,description)
% Based on spm_file_merge
% Function writes 4D data to 4D NIFTI
% data      - data to write (4D)
% data_vol  - volume info for new NIFTI
% fname  - filename for output 4D volume [defaults: '4D.nii']
%          Unless explicit, output folder is the one containing first image
% dt     - datatype (see spm_type) [defaults: 0]
%          0 means same datatype than first input volume
%
% description - description into header
%
% V4     - spm_vol struct of the 4D volume

dim=size(data);
N_scans=dim(4);

for i=1:N_scans
    V(i)=data_vol;
end

% MG 10.1.2017
% ------------
% V4 = spm_file_merge(V,fname,dt)
% Concatenate 3D volumes into a single 4D volume
% FUNCTION V4 = spm_file_merge(V,fname,dt)
% V      - images to concatenate (char array or spm_vol struct)
% fname  - filename for output 4D volume [defaults: '4D.nii']
%          Unless explicit, output folder is the one containing first image
% dt     - datatype (see spm_type) [defaults: 0]
%          0 means same datatype than first input volume
%
% V4     - spm_vol struct of the 4D volume
%__________________________________________________________________________
%
% For integer datatypes, the file scale factor is chosen as to maximise
% the range of admissible values. This may lead to quantization error
% differences between the input and output images values.
%__________________________________________________________________________
% Copyright (C) 2009-2012 Wellcome Trust Centre for Neuroimaging

% John Ashburner
% $Id: spm_file_merge.m 5040 2012-11-07 12:36:23Z guillaume $

%-Input: V
%--------------------------------------------------------------------------
% % % if ~nargin
% % %     [V,sts] = spm_select([1 Inf],'image','Select images to concatenate');
% % %     if ~sts, return; end
% % % end
% % % if ischar(V)
% % %     V = spm_vol(V);
% % % elseif iscellstr(V)
% % %     V = spm_vol(char(V));
% % % end
ind  = cat(1,V.n);
N    = cat(1,V.private);

%-Set scalefactors and offsets
%==========================================================================
d = cat(1,V.dt); d = d(:,1);
s = cat(2,V.pinfo); s = s(1,:);
o = cat(2,V.pinfo); o = o(2,:);

%-Reuse parameters of input images if same scalefactor, offset and datatype
%--------------------------------------------------------------------------
if length(unique(s)) == 1 && length(unique(o)) == 1 ...
        && length(unique(d)) == 1 && d(1) == dt
    sf  = V(1).pinfo(1);
    off = V(1).pinfo(2);
else

    dmx  = spm_type(dt,'maxval');
    dmn  = spm_type(dt,'minval');

    %-Integer datatypes: scale to min/max of overall data
    %----------------------------------------------------------------------
    if isfinite(dmx)
        spm_progress_bar('Init',numel(V),'Computing scale factor','Volumes Complete');
        mx      = -Inf;
        mn      = Inf;
        for i=1:numel(V)
%             dat = V(i).private.dat(:,:,:,ind(i,1),ind(i,2));
            dat = data(:,:,:,i);
            dat = dat(isfinite(dat));
            mx  = max(mx,max(dat(:)));
            mn  = min(mn,min(dat(:)));
            spm_progress_bar('Set',i);
        end
        spm_progress_bar('Clear');
        if isempty(mx), mx = 0; end
        if isempty(mn), mn = 0; end

        if mx~=mn
            if dmn < 0
                sf = max(mx/dmx,-mn/dmn);
            else
                sf = mx/dmx;
            end
            off    = 0;
        else
            sf     = mx/dmx;
            off    = 0;
        end
    
    %-floating precison: no scaling
    %----------------------------------------------------------------------
    else
        sf         = 1;
        off        = 0;
    end
end

%-Create and write 4D volume image
%==========================================================================
spm_unlink(fname);

ni         = nifti;
ni.dat     = file_array(fname,...
                        [V(1).dim numel(V)],...
                        [dt spm_platform('bigend')],...
                        0,...
                        sf,...
                        off);
ni.mat     = N(1).mat;
ni.mat0    = N(1).mat;
ni.descrip = description;

create(ni);
spm_progress_bar('Init',size(ni.dat,4),'Saving 4D image','Volumes Complete');
for i=1:size(ni.dat,4)
    ni.dat(:,:,:,i) = data(:,:,:,i);
    spm_get_space([ni.dat.fname ',' num2str(i)], V(i).mat);
    spm_progress_bar('Set',i);
end
spm_progress_bar('Clear');

if nargout
    V4 = spm_vol(ni.dat.fname);
end
mat_fname=[fname(1:end-4),'.mat'];
spm_unlink(mat_fname)
