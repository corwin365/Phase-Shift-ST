clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%third draft of generalised 2D+1 ST analysis
%Corwin Wright, c.wright@bath.ac.uk, 20200624
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use the wavelength indices from the 3DST to pin down the waves, then use
%2DST phase differences to fine-tune their vertical wavelengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%minimum phase change to be meaningful
MindP = pi/10;

%data to analyse
Settings.DayScale         = datenum(2008,1,127);
Settings.Granules         = [56];

% Settings.DayScale  = datenum(2010,10,16);
% Settings.Granules  = [186];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.DayScale)

  Date = Settings.DayScale(iDay);
  
  for iGranule=1:1:numel(Settings.Granules)
    
    Gid = Settings.Granules(iGranule);
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% load data as a granule-pair
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %get granule
    Airs  = cat_struct(prep_airs_3d(Date,Gid,  'PreSmooth',[3,3,1]),   ...
                       prep_airs_3d(Date,Gid+1,'PreSmooth',[3,3,1]), ...
                       2,{'MetaData','Source','ret_z'});
                     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 3dst the granule-pair
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [ST3D,Airs] = gwanalyse_airs_3d(Airs,'TwoDPlusOne',false,'ZRange',[10,70]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2D+1 ST the granule pair
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %do a *****full***** 2D ST at each height
    textprogressbar('STing levels ')
    PhaseStore = NaN.*ST3D.C;
    for iLevel=1:1:numel(Airs.ret_z)
      
      %do the 2DST
      ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',[1,1].*.5);
      
      %retain the output field
      if iLevel == 1;
        STStore = ST.ST;
        freqs = ST.freqs;
      else  STStore = cat(5,STStore,ST.ST); 
      end
      
      textprogressbar(100.*iLevel./numel(Airs.ret_z))
    end
    clear iLevel ST F1 F2    
    textprogressbar('!')

    %%
    %convert to complex cospectra
    CC = STStore.*NaN;
    textprogressbar('CCing levels ')
    for iLevel=1:1:numel(Airs.ret_z)-2
      textprogressbar(iLevel./(numel(Airs.ret_z)-1).*100)
      CC(:,:,:,:,iLevel) = STStore(:,:,:,:,iLevel) .* conj(STStore(:,:,:,:,iLevel+2));
    end; clear iLevel  
    textprogressbar('!')

    %find the nearest wavelength combination to the 3DST fits
    F1_3D = ST3D.F1;  idx_F1_3D = NaN.*F1_3D;
    F2_3D = ST3D.F2;  idx_F2_3D = NaN.*F2_3D;
    for iPoint = 1:1:numel(F1_3D);
      idx_F1_3D(iPoint) = closest(freqs{1},F1_3D(iPoint));
      idx_F2_3D(iPoint) = closest(freqs{2},F2_3D(iPoint));
    end; clear iPoint F1_3D F2_3D
    
   
    %retain the phase difference of these fits
    sz = size(CC);
    CC        = reshape(CC,sz(1),sz(2),sz(3)*sz(4),sz(5));
    idx_F1_3D = reshape(idx_F1_3D,sz(3)*sz(4),sz(5));
    idx_F2_3D = reshape(idx_F2_3D,sz(3)*sz(4),sz(5));
    PhiStore  = NaN(sz(3)*sz(4),sz(5));
    for iLevel=1:1:sz(5)
      for iPoint=1:1:numel(sz(3)*sz(4));
        Phi = CC(idx_F1_3D(iPoint,iLevel),idx_F1_3D(iPoint,iLevel),iPoint,iLevel);
        PhiStore(iPoint,iLevel) = Phi;
      end
    end
    PhiStore = reshape(PhiStore,sz(3),sz(4),sz(5));
    clear sz CC idx_F1_3D idx_F2_3D iLevel iPoint Phi

    
    %convert each levels CC values to covarying phase
    AlldP = angle(PhiStore);
    clear PhiStore

    AlldP(abs(AlldP) < MindP) = NaN;
    
    %convert phase change to wavelength
    sz = size(AlldP);
    dZ = [0;diff(Airs.ret_z)+ circshift(diff(Airs.ret_z),1)]; %levels this makes wonky wll be removed later
    Lambda = permute(repmat(dZ,1,sz(1),sz(2)),[2,3,1])./AlldP.*2*pi;
    clear sz AlldP dZ

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% retain the stuff we want, for comparison
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        

    %and store
    ST3D.F3_2Dp1 = Lambda;

    
  end
end
