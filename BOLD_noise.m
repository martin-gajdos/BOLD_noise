function BOLD_noise
% This is main function of the toolbox. The tool is designed to create
% simulated noise into BOLD data based on properties observed in real
% datasets. The simulated noise consists particularly from:
% 1. Noise caused by movements of subject's head
% 2. Noise caused by breathing and heartbeats
% 3. 1/f noise
%
% Toolbox requires SPM12
%
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
output_name='test_sumu1.nii';

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

% options.inbrain_area_mask=fullfile(spm('dir'),'tpm','mask_ICV.nii'); % Intracranial mask

% %% Movement noise simulation
[M_noise,movement_parameters]=get_movement_noise(options);

options.M_noise=M_noise;

%% Breathing and heartbeats simulation
[P_noise,physio_parameters]=get_physio_noise(options);

%% 1/f noise
F_noise=get_1_f_noise(options);

BOLD_noise=M_noise+P_noise+F_noise;

%% Writing result



fname=fullfile(output_path,output_name);
% 
spm_unlink(fname)


dim=size(BOLD_noise);
dim=dim(1:3);
for i=1:no_of_scans
    BOLD_noise_struct(i) = deal(struct(...
		'fname',	[],...
		'dim',		dim,...
        'dt',       [spm_type('int16') 0],...
		'mat',		movement_parameters.mask_vol.mat,...
		'pinfo',	[1 0 0]',...
		'descrip',	'BOLD noise data',...
        'n',        [],...
        'private',  []));
    BOLD_noise_struct(i).fname   = fname;
    BOLD_noise_struct(i) = spm_create_vol(BOLD_noise_struct(i));
    spm_write_vol(BOLD_noise_struct(i),squeeze(BOLD_noise(:,:,:,i)));
end

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
    rand_val{ii,1}=NaN(mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3));
    
    tmp_mu_map_name=fullfile(movement_maps_path,movement_maps_mu_names{ii});
    tmp_mu_vol=spm_vol(tmp_mu_map_name);
    
    [tmp_mu,tmp_XYZ] = spm_read_vols(tmp_mu_vol);
    
    for jj=1:Ncoords
        tmp_vx=coords_vx(:,jj);
        tmp_mu_val=tmp_mu(tmp_vx(1),tmp_vx(2),tmp_vx(3)); 

        rand_val{ii}(tmp_vx(1),tmp_vx(2),tmp_vx(3))=random('Exponential',tmp_mu_val);  % Statistic toolbox

    end
    
    tmp_reg=squeeze(reg(:,ii));
    tmp_reg=reshape(tmp_reg, [1, 1, 1, numel(tmp_reg)]);
    tmp_reg=repmat(tmp_reg, [mask_vol.dim(1),mask_vol.dim(2),mask_vol.dim(3),1]);
    
    rand_val_rep=repmat(rand_val{ii},[1,1,1,nscan]);
    
    tmp_reg2=rand_val_rep.*tmp_reg;
    M_noise=M_noise+tmp_reg2;
  
end

movement_parameters.regressors=reg;
movement_parameters.map_names=movement_maps_mu_names;
movement_parameters.movement_maps_path=movement_maps_path;
movement_parameters.movement_maps_mask_name=movement_maps_mask_name;
movement_parameters.movement_maps_mu_names=movement_maps_mu_names;
movement_parameters.mask_vol=mask_vol;

function [P_noise,physio_parameters]=get_physio_noise(options)
% Function simulates breathing and heartbeats of subject and its reflection 
% into the data. It is based on observed impulse responses caused by breath
% and heart artifact, fitted by RETROICOR. The noise keeps also correlated
% structure in MNI space as observed in source datasets.

testing=1;
if testing==1
    P_noise=rand(size(options.M_noise));
    physio_parameters=[];
else
    %% Heartbeats timecourse simulation
    
    %% Generating heartbeat impulse response
    
    %% Preparing heartbeat 4D noise data
    heart_noise
    heart_noise_const 
    
    %% Breathing timecourse simulation
    
    %% Generating breathing impulse response
    
    %% Preparing breathing 4D noise data
    breath_noise
    breath_noise_const 
    
    P_noise=heart_noise_const*heart_noise+breath_noise_const*breath_noise;
end

function F_noise=get_1_f_noise(options)
% Function simulates 1/f noise. It is based on noise observed in source
% datasets.


testing=1;
if testing==1
    F_noise=rand(size(options.M_noise));
end





