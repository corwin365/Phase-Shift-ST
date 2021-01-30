clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first attempt at using phase shift of 2D STs to get lambda_z
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

% Settings.Date        = datenum(2010,10,16);
% Settings.Granule     = 186;
% Settings.BasisHeight = 39; %km 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data and ST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get granule
Airs1 = prep_airs_3d(Settings.Date,Settings.Granule);
Airs2 = prep_airs_3d(Settings.Date,Settings.Granule+1);
Airs = Airs1;
Vars = {'l1_lat','l1_lon','l1_time','ret_temp','Tp','BG'};
for iVar=1:1:numel(Vars)
  Airs.(Vars{iVar}) = cat(2,Airs1.(Vars{iVar}),Airs2.(Vars{iVar}));
end
clear Airs1 Airs2 Vars



% % %drop height extrema
% % InZRange = find(Airs.ret_z >= 20 & Airs.ret_z <= 60);
% % Airs.ret_z = Airs.ret_z(InZRange);
% % Airs.Tp  = Airs.Tp(:,:,InZRange);

%smooth a bit
Airs.Tp = smoothn(Airs.Tp,[3,3,3]);

%do a *****full***** 2D ST at each height
for iLevel=1:1:numel(Airs.ret_z)
  disp(iLevel)
  %do ST
  ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel),'FullST',true,'c',[1,1].*.5);

  %store
  if iLevel == 1; STStore = ST.ST;
  else            STStore = cat(5,STStore,ST.ST);
  end
end

%smooth the STStore space, separately in phase and amplitude
A = real(STStore);
B = imag(STStore);

A = smoothn(A,[5,5,1,1,1]);
B = smoothn(B,[5,5,1,1,1]);

STStore = complex(A,B);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ok. at the basis level, find the peak frequency for each pixel
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
%% compute complex cospectra relative to basis level
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CC = STStore;
for iLevel=1:1:numel(Airs.ret_z)  
  CC(:,:,:,:,iLevel) = CC(:,:,:,:,iLevel) .* conj(BasisST);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pull out amplitudes, phase at basis, and phase differences 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%basis level
A = NaN(sz(3),sz(4));
for iX=1:1:sz(3);
  for iY = 1:1:sz(4);
    A(iX,iY) = abs(BasisST(idx1(iX,iY),idx2(iX,iY),iX,iY));
  end
end

%other levels
dP = NaN(sz(3),sz(4),numel(Airs.ret_z));  % phase difference
CA = NaN(sz(3),sz(4),numel(Airs.ret_z));    %covarying amplitude
for iX=1:1:sz(3)
  for iY=1:1:sz(4)
    dP(iX,iY,:) = angle(CC(idx1(iX,iY),idx2(iX,iY),iX,iY,:));
    CA(iX,iY,:) = sqrt(abs( CC(idx1(iX,iY),idx2(iX,iY),iX,iY,:)));
  end
end

dP = dP ./ (2*pi); %scale into range -1 to 1

% %discard very small dP
% dP(abs(dP) < 0.05) = NaN;

% %smooth dP a bit
% dP = smoothn(dP,[5,5,1]);

%convert phase change to vertical wavelength
dz = diff(Airs.ret_z); dz(end+1) = dz(end); for iLev=2:1:numel(dz); dz(iLev) = dz(iLev - 1)+dz(iLev); end
dz = dz - dz(zidx);

lz = repmat(permute(dz,[2,3,1]),sz(3),sz(4),1) ./dP;

