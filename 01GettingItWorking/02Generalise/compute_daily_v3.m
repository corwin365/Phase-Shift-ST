clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%third draft of generalised 2D+1 ST analysis
%Corwin Wright, c.wright@bath.ac.uk, 20200624
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%let's go for something more final. try and match the wavelengths in a 
%sensible way between k/l and m, and then tidy up the code to do decent outputs


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
    
    ST3D = gwanalyse_airs_3d(Airs,'TwoDPlusOne',false);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% 2D+1 ST the granule pair
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    %do a *****full***** 2D ST at each height
    textprogressbar('STing levels ')
    for iLevel=1:1:numel(Airs.ret_z)
      
      %do the 2DST
      ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',[1,1].*.5);
      
      %retain the full complex ST, and *also* the **indices** of the fitted horizontal waves
      F1s = unique(ST.F1); F1 = ST.F1;
      for iF=1:1:numel(F1s); F1(F1 == F1s(iF)) = closest(ST.freqs{1},F1s(iF)); end
      F2s = unique(ST.F2); F2 = ST.F2;
      for iF=1:1:numel(F2s); F2(F2 == F2s(iF)) = closest(ST.freqs{2},F2s(iF)); end 
      clear F1s F2s iF
      
      
      if iLevel == 1;
        STStore = ST.ST;
        freqs = ST.freqs;
        F1Store = F1; F2Store = F2;
        kStore = ST.k; lStore = ST.l;
        AStore = ST.A; 
      else
        STStore = cat(5,STStore,ST.ST); AStore = cat(3,AStore,ST.A);
        F1Store = cat(3,F1Store,  F1); F2Store = cat(3,F2Store,F2  );
        kStore  = cat(3, kStore,ST.k); lStore  = cat(3, lStore,ST.l);
      end
      textprogressbar(iLevel./numel(Airs.ret_z).*100)
    end
    clear iLevel ST F1 F2
    textprogressbar('!')
    
    %convert to complex cospectra
    CC = STStore.*NaN;
    textprogressbar('CCing levels ')
    for iLevel=2:1:numel(Airs.ret_z)-1
      textprogressbar(iLevel./(numel(Airs.ret_z)-1).*100)
      CC(:,:,:,:,iLevel) = STStore(:,:,:,:,iLevel-1) .* conj(STStore(:,:,:,:,iLevel+1));
    end; clear iLevel
    textprogressbar('!')

    
    %pull the signal out of the array
    sz = size(CC);
    for iX=1:1:sz(3);
      for iY=1:1:sz(4);
        for iZ=1:1:sz(5)
          CC(1,1,iX,iY,iZ) = CC(F1Store(iX,iY,iZ),F2Store(iX,iY,iZ),iX,iY,iZ);   
        end
      end
    end
    clear iX iY iZ
    CC = squeeze(CC(1,1,:,:,:));
    
    %convert each levels CC values to covarying phase
    AlldP = angle(CC);
    clear CC

    AlldP(abs(AlldP) < MindP) = NaN;
    
    %convert phase change to wavelength
    sz = size(AlldP);
    dZ = [0;diff(Airs.ret_z)+ circshift(diff(Airs.ret_z),1)]; %levels this makes wonky wll be removed later
    Lambda = permute(repmat(dZ,1,sz(1),sz(2)),[2,3,1])./AlldP.*2*pi;
    clear sz AlldP dZ

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% retain the stuff we want, for comparison
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    %convert horizontal wavelength *indices* back to *wavelengths*
    fr = freqs{1}; for iF=1:1:numel(fr);F1Store(F1Store == iF) = fr(iF);end
    fr = freqs{2}; for iF=1:1:numel(fr);F2Store(F2Store == iF) = fr(iF);end
    clear fr iF freqs 
    
    %apply sign convention and store
    Positive = find(Lambda > 0);

    F1Store(Positive) = -1.*F1Store(Positive);
    F2Store(Positive) = -1.*F2Store(Positive);
    kStore( Positive) = -1.*kStore( Positive);
    lStore( Positive) = -1.*lStore( Positive);
    Lambda( Positive) = -1.*Lambda( Positive);
    clear Positive
    
    
    %trim to height range
    InZRange = find(Airs.ret_z >= 20 ...
                  & Airs.ret_z <= 60);
    
    %and store
    ST2D = struct();
    ST2D.F1 = F1Store(  :,:,InZRange); clear F1Store;
    ST2D.F2 = F2Store(  :,:,InZRange); clear F2Store;
    ST2D.k  = kStore(   :,:,InZRange);  clear kStore;
    ST2D.l  = lStore(   :,:,InZRange);  clear lStore;
    ST2D.F3 = 1./Lambda(:,:,InZRange);    
    ST2D.m  = 1./Lambda(:,:,InZRange); clear Lambda;
    ST2D.A  = AStore(   :,:,InZRange); clear AStore

    
  end
end
