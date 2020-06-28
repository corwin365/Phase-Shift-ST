clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% third attempt at using phase shift of 2D STs to get lambda_z
% collaborating with Laura Holt and Joan Alexander
%
%Corwin Wright, c.wright@bath.ac.uk, 2019/10/28
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.Date        = datenum(2008,1,127);
Settings.Granule     = 56;
Settings.BasisHeight = 39; %km 
Settings.MinAmplitude     = 1; %cutoff amplitude to be used (arbitrary)
Settings.MinFracAmplitude = 0.1; %fraction of max amplitude below which wave dies
Settings.SampleIndices = [174,34];
Settings.Mode = 2;

% Settings.Date        = datenum(2010,10,16);
% Settings.Granule     = 186;
% Settings.BasisHeight = 39; %km 
% Settings.MinAmplitude     = 0.5; %cutoff amplitude to be used (arbitrary)
% Settings.MinFracAmplitude = 0.1; %fraction of max amplitude below which wave dies
% Settings.Mode = 2;
% Settings.SampleIndices = [128,38];

% Settings.Date        = datenum(2007,1,13);
% Settings.Granule     = 122;
% Settings.BasisHeight = 39; %km 
% Settings.MinAmplitude     = 0.5; %cutoff amplitude to be used (arbitrary)
% Settings.MinFracAmplitude = 0.1; %fraction of max amplitude below which wave dies
% Settings.SampleIndices = [78,57];
% Settings.Mode = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data and ST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if Settings.Mode == 2;
  %get granule
  Airs1 = prep_airs_3d(Settings.Date,Settings.Granule);
  Airs2 = prep_airs_3d(Settings.Date,Settings.Granule+1);
  Airs = Airs1;
  Vars = {'l1_lat','l1_lon','l1_time','ret_temp','Tp','BG'};
  for iVar=1:1:numel(Vars)
    Airs.(Vars{iVar}) = cat(2,Airs1.(Vars{iVar}),Airs2.(Vars{iVar}));
  end
clear Airs1 Airs2 Vars iVar
else
  Airs= prep_airs_3d(Settings.Date,Settings.Granule);
end


% %smooth a bit
% Airs.Tp = smoothn(Airs.Tp,[3,3,1]);

%do a *****full***** 2D ST at each height
textprogressbar('STing levels ')
for iLevel=1:1:numel(Airs.ret_z)
  textprogressbar(iLevel./numel(Airs.ret_z).*100)
  %do ST
  ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',[1,1].*.5);

  %store
  if iLevel == 1; STStore = ST.ST;
  else            STStore = cat(5,STStore,ST.ST);
  end
end
clear iLevel ST
textprogressbar('!')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute complex cospectra relative to level above for each pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = STStore;
textprogressbar('CCing levels ')
for iLevel=1:1:numel(Airs.ret_z)-1  
  textprogressbar(iLevel./(numel(Airs.ret_z)-1).*100)
  CC(:,:,:,:,iLevel) = CC(:,:,:,:,iLevel) .* conj(CC(:,:,:,:,iLevel+1));
end; clear iLevel
textprogressbar('!')

%remove extraneous level, and produce new height scale
CC   = CC(:,:,:,:,1:end-1);
NewZ = Airs.ret_z(1:end-1) + diff(Airs.ret_z)./2;
dZ   = diff(Airs.ret_z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% at the basis level, find the peak frequency for each pixel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find the level
zidx = closest(Airs.ret_z,Settings.BasisHeight);

%pull out this ST
BasisST = STStore(:,:,:,:,zidx);

%reshape to facilitate max() operation
sz = size(BasisST);
BasisST = reshape(BasisST,sz(1)*sz(2),sz(3)*sz(4));

%take max(abs()) along spectral dimension
[~,idx] = max(abs(BasisST),[],1);

%put BasisST back
BasisST = reshape(BasisST,sz);

%convert to 2d indices
[idx1,idx2] = ind2sub([sz(1),sz(2)],idx);
idx1 = reshape(idx1,sz(3),sz(4));
idx2 = reshape(idx2,sz(3),sz(4));

clear idx


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull out amplitudes, phase at basis, and phase differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%convert each levels CC values to covarying amplitude and phase
AllA  = sqrt(abs(CC));
AlldP = angle(CC);

%find basis level ST.A, for tuning
A = NaN(sz(3),sz(4));
for iX=1:1:sz(3);
  for iY = 1:1:sz(4);
    A(iX,iY) = abs(BasisST(idx1(iX,iY),idx2(iX,iY),iX,iY));
  end
end

%%




%ok. let's go pixelwise for now for conceptual ease. we can speed this up later.
sz = size(AllA);
lz = NaN(sz(3:end));
azs = lz;

for iX=2:1:sz(3)-1; %xt
  for iY=2:1:sz(4)-1; %at
    
    %identify the peak frequency AT THE BASIS LEVEL for this pixel
    jX = idx1(iX,iY);
    jY = idx2(iX,iY);
    
    %if amplitude at the basis level is below cutoff, then skip
    if AllA(jX,jY,iX,iY,zidx) < Settings.MinAmplitude; continue ; end
    
    %ok. find the covarying amplitude and phase with this frequency at each z level
    
    %take the mean of the nine SPATIAL pixels around it, rather than just the
    %central value. SPECTRAL points remain distinct
    pmiX = iX;%(-1:1:1) + iX;
    pmiY = iY;%(-1:1:1) + iY;
    
    Az  = squeeze(nanmean(nanmean(AllA( jX,jY,pmiX,pmiY,:),3),4));
    dPz = squeeze(nanmean(nanmean(AlldP(jX,jY,pmiX,pmiY,:),3),4));
    azs(iX,iY,:) = Az;
    
    %hence, find the associated wavelengths
    lz(iX,iY,:) = dZ ./( dPz ./ (2*pi));

    %and mask out any regions where the amplitude has dropped sharply
    Bad = find(Az < 0.1 .* max(Az));
    lz(iX,iY,Bad) = NaN;

    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

set(gcf,'color','w')

subplot(2,3,1)
pc(Airs.Tp(:,:,zidx))
redyellowblue64
caxis([-1,1].*5)
title(['T''  at ',num2str(Settings.BasisHeight),'km'])

subplot(2,3,2)
pc(A)
redyellowblue64
caxis([0,1].*3)
title(['Cospectral amplitude  at ',num2str(Settings.BasisHeight+1.5),'km'])

subplot(2,3,3)
pc(lz(:,:,zidx))
caxis([-1,1].*30)
title(['\lambda_z  at ',num2str(Settings.BasisHeight+1.5),'km for A>',num2str(Settings.MinAmplitude),'K'])


subplot(2,3,4)
plot(squeeze(Airs.Tp(Settings.SampleIndices(2),Settings.SampleIndices(1),:)),Airs.ret_z,'k-o')
ylim([20 60])
title(['T'' at [',num2str(Settings.SampleIndices(1)),',',num2str(Settings.SampleIndices(2)),']'])

subplot(2,3,5)
plot(squeeze(azs(Settings.SampleIndices(2),Settings.SampleIndices(1),:)),NewZ)
ylim([20 60])
title(['A at [',num2str(Settings.SampleIndices(1)),',',num2str(Settings.SampleIndices(2)),']'])


subplot(2,3,6)
plot(smoothn(squeeze(lz(Settings.SampleIndices(2),Settings.SampleIndices(1),:)),[3,1]),NewZ)
ylim([20 60])
title(['\lambda_z at [',num2str(Settings.SampleIndices(1)),',',num2str(Settings.SampleIndices(2)),']'])
xlim([-1,1].*40)