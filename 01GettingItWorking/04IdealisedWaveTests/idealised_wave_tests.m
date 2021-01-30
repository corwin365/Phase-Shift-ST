clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the artificial waves from Hindley et al 2019 to test performance of 2D+1 ST
%
%Corwin Wright, c.wright@bath.ac.uk, 2020/07/22
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.NLevels = 100;%5:1:10;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the artificial wave field and tidy up the variable space
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iResolution = 1:1:numel(Settings.NLevels)
  
  %load the data
  load('specified_wave_field_workspace.mat');
  
  %reformat to look like an AIRS granule to the analysis routine
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %geolocation
  Field.l1_lat   = Workspace.Workspace.D1(:,:,1);
  Field.l1_lon   = Workspace.Workspace.D2(:,:,1);
  Field.ret_z    = squeeze(Workspace.Workspace.D3(1,1,:));
  
  %downsample field
  z = linspace(-100,100,Settings.NLevels(iResolution));
  W = reshape(Workspace.Workspace.Wavefield,201*201,201);
  W = reshape(interp1(Field.ret_z,W',z,'spline')',201,201,numel(z));
  Field.ret_z = z'; 

  %and create data fields
  Field.Tp       = W;
  Field.BG       = zeros(size(Field.Tp));
  Field.ret_temp = Field.Tp;

  clear Workspace z W
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% apply 3DST and both 2D+1 ST versions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  [ST3D,~] = gwanalyse_airs_3d(Field,                       ...
                               'NotAirsData',    true,       ...
                               'HeightScaling',  false,      ...
                               'ZRange',         [-100,100], ...
                               'Spacing',        [1,1,1],    ...
                               'TwoDPlusOne',    true);
                                

end