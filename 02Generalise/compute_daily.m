clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first draft of generalised 2D+1 ST analysis
%Corwin Wright, c.wright@bath.ac.uk, 20200624
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Settings.DayScale         = datenum(2008,1,127);
Settings.Granules         = [56];
Settings.BasisHeight      = 39; %km 
Settings.MinAmplitude     = 0.5; %cutoff amplitude to be used (arbitrary)
Settings.MinFracAmplitude = 0.1; %fraction of max amplitude below which wave dies
Settings.SampleIndices    = [128,38];

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
    
    ST3D = gwanalyse_airs_3d(Airs);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2D+1 ST the granule pair
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %do a *****full***** 2D ST at each height
    textprogressbar('STing levels ')
    for iLevel=1:1:numel(Airs.ret_z)
      ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',[1,1].*.5);
      
      %store
      if iLevel == 1; STStore = ST.ST; PointSpacing = ST.point_spacing;
      else            STStore = cat(5,STStore,ST.ST);
      end
      textprogressbar(iLevel./numel(Airs.ret_z).*100)
    end
    clear iLevel ST
    textprogressbar('!')
    
    %convert to complex cospectra
    CC = STStore;
    textprogressbar('CCing levels ')
    for iLevel=1:1:numel(Airs.ret_z)-2
      textprogressbar(iLevel./(numel(Airs.ret_z)-1).*100)
      CC(:,:,:,:,iLevel) = CC(:,:,:,:,iLevel) .* conj(CC(:,:,:,:,iLevel+2));
    end; clear iLevel
    textprogressbar('!')
    
    %remove extraneous level, and produce new height scale
    CC   = CC(:,:,:,:,1:end-1); 
    NewZ = Airs.ret_z(1:end-1) + diff(Airs.ret_z);
    dZ   = diff(Airs.ret_z)+ circshift(diff(Airs.ret_z),1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% find the strongest signal for each pixel and work out the vertical wavelength
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    %identify strongest spectral signal for each voxel
    sz = size(CC);
    CC = reshape(CC,sz(1)*sz(2),sz(3)*sz(4)*sz(5));
    [~,MaxIdx] = max(abs(CC),[],1);
    [idxX,idxY] = ind2sub(sz([1,2]),MaxIdx);
    CC = reshape(CC,sz);
    idxX = reshape(idxX,sz(3:end)); idxY = reshape(idxY,sz(3:end)); 
    
    %pull this signal out of the array
    for iX=1:1:sz(3);
      for iY=1:1:sz(4);
        for iZ=1:1:sz(5)
          CC(1,1,iX,iY,iZ) = CC(idxX(iX,iY,iZ),idxY(iX,iY,iZ),iX,iY,iZ);   
        end
      end
    end
    clear iX iY MaxIdx iZ
    CC = squeeze(CC(1,1,:,:,:));
    
    %convert each levels CC values to covarying amplitude and phase
    AllA  = sqrt(abs(CC));
    AlldP = angle(   CC);
    
    %convert phase change to wavelength
    sz = size(AlldP);
    Lambda = permute(repmat(dZ,1,sz(1),sz(2)),[2,3,1])./AlldP.*2*pi;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute the corresponding horizontal wavelengths
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    
    %produce along-track and cross-track cross-spectra at the same indices
    for iX=1:1:sz(3);
      for iY=1:1:sz(4);
        for iZ=1:1:sz(5)
          STStore(1,1,iX,iY,iZ) = STStore(idxX(iX,iY,iZ),idxY(iX,iY,iZ),iX,iY,iZ);   
        end
      end
    end
    STStore = squeeze(STStore(1,1,:,:,:));
    
    
    stop
    
    
  end
end
% % % % % %   else
% % % % %     
% % % % %   end
% % % % % 
% % % % % stop
% % % % % 
% % % % % % %smooth a bit
% % % % % % Airs.Tp = smoothn(Airs.Tp,[3,3,1]);
% % % % % 
% % % % % %do a *****full***** 2D ST at each height
% % % % % textprogressbar('STing levels ')
% % % % % for iLevel=1:1:numel(Airs.ret_z)
% % % % %   textprogressbar(iLevel./numel(Airs.ret_z).*100)
% % % % %   %do ST
% % % % %   ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',[1,1].*.5);
% % % % % 
% % % % %   %store
% % % % %   if iLevel == 1; STStore = ST.ST;
% % % % %   else            STStore = cat(5,STStore,ST.ST);
% % % % %   end
% % % % % end
% % % % % clear iLevel ST
% % % % % textprogressbar('!')
% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% compute complex cospectra relative to level above for each pair
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % CC = STStore;
% % % % % textprogressbar('CCing levels ')
% % % % % for iLevel=1:1:numel(Airs.ret_z)-1  
% % % % %   textprogressbar(iLevel./(numel(Airs.ret_z)-1).*100)
% % % % %   CC(:,:,:,:,iLevel) = CC(:,:,:,:,iLevel) .* conj(CC(:,:,:,:,iLevel+1));
% % % % % end; clear iLevel
% % % % % textprogressbar('!')
% % % % % 
% % % % % %remove extraneous level, and produce new height scale
% % % % % CC   = CC(:,:,:,:,1:end-1);
% % % % % NewZ = Airs.ret_z(1:end-1) + diff(Airs.ret_z)./2;
% % % % % dZ   = diff(Airs.ret_z);
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% at the basis level, find the peak frequency for each pixel
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %find the level
% % % % % zidx = closest(Airs.ret_z,Settings.BasisHeight);
% % % % % 
% % % % % %pull out this ST
% % % % % BasisST = STStore(:,:,:,:,zidx);
% % % % % 
% % % % % %reshape to facilitate max() operation
% % % % % sz = size(BasisST);
% % % % % BasisST = reshape(BasisST,sz(1)*sz(2),sz(3)*sz(4));
% % % % % 
% % % % % %take max(abs()) along spectral dimension
% % % % % [~,idx] = max(abs(BasisST),[],1);
% % % % % 
% % % % % %put BasisST back
% % % % % BasisST = reshape(BasisST,sz);
% % % % % 
% % % % % %convert to 2d indices
% % % % % [idx1,idx2] = ind2sub([sz(1),sz(2)],idx);
% % % % % idx1 = reshape(idx1,sz(3),sz(4));
% % % % % idx2 = reshape(idx2,sz(3),sz(4));
% % % % % 
% % % % % clear idx
% % % % % 
% % % % % 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% pull out amplitudes, phase at basis, and phase differences 
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % 
% % % % % %convert each levels CC values to covarying amplitude and phase
% % % % % AllA  = sqrt(abs(CC));
% % % % % AlldP = angle(CC);
% % % % % 
% % % % % %find basis level ST.A, for tuning
% % % % % A = NaN(sz(3),sz(4));
% % % % % for iX=1:1:sz(3);
% % % % %   for iY = 1:1:sz(4);
% % % % %     A(iX,iY) = abs(BasisST(idx1(iX,iY),idx2(iX,iY),iX,iY));
% % % % %   end
% % % % % end
% % % % % 
% % % % % %%
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %ok. let's go pixelwise for now for conceptual ease. we can speed this up later.
% % % % % sz = size(AllA);
% % % % % lz = NaN(sz(3:end));
% % % % % azs = lz;
% % % % % 
% % % % % for iX=2:1:sz(3)-1; %xt
% % % % %   for iY=2:1:sz(4)-1; %at
% % % % %     
% % % % %     %identify the peak frequency AT THE BASIS LEVEL for this pixel
% % % % %     jX = idx1(iX,iY);
% % % % %     jY = idx2(iX,iY);
% % % % %     
% % % % %     %if amplitude at the basis level is below cutoff, then skip
% % % % %     if AllA(jX,jY,iX,iY,zidx) < Settings.MinAmplitude; continue ; end
% % % % %     
% % % % %     %ok. find the covarying amplitude and phase with this frequency at each z level
% % % % %     
% % % % %     %take the mean of the nine SPATIAL pixels around it, rather than just the
% % % % %     %central value. SPECTRAL points remain distinct
% % % % %     pmiX = iX;%(-1:1:1) + iX;
% % % % %     pmiY = iY;%(-1:1:1) + iY;
% % % % %     
% % % % %     Az  = squeeze(nanmean(nanmean(AllA( jX,jY,pmiX,pmiY,:),3),4));
% % % % %     dPz = squeeze(nanmean(nanmean(AlldP(jX,jY,pmiX,pmiY,:),3),4));
% % % % %     azs(iX,iY,:) = Az;
% % % % %     
% % % % %     %hence, find the associated wavelengths
% % % % %     lz(iX,iY,:) = dZ ./( dPz ./ (2*pi));
% % % % % 
% % % % %     %and mask out any regions where the amplitude has dropped sharply
% % % % %     Bad = find(Az < 0.1 .* max(Az));
% % % % %     lz(iX,iY,Bad) = NaN;
% % % % % 
% % % % %     
% % % % %   end
% % % % % end
