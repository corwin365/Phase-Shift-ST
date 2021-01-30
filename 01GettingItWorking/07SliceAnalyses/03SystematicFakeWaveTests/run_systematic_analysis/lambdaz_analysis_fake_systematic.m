clearvars -except InfoFile
load(InfoFile)
delete(InfoFile)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure showing the 2D+1 analysis for four waves as 2D cuts
%
% Corwin Wright, c.wright@bath.ac.uk, 2020/08/24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%how many elements to store each side of wave centre
Settings.Length = 35;

%save file name
Settings.OutFile = SaveFile; clear SaveFile


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% primary loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


DataStore = struct(); %where our data will go
for iCombo = 1:1:size(ComboData,1);
  
  %prepare storage box
  Case = Cases.(['Case',num2str(ComboData(iCombo,6))]);
  Granules = Case{2};

  disp('=======================================================')
  disp(['Processing combination ',num2str(iCombo),' of ',num2str(size(ComboData,1))])
  disp('=======================================================')  
  tic

  for iGranule=1:1:numel(Granules)
    
    %load granule
    [Airs,Spacing] = prep_airs_3d(Case{1},Granules(iGranule),'PreSmooth',[3,3,1]);
    
    %keep only geolocation and zeros otherwise?
    if ComboData(iCombo,7) == 1; Airs.Tp(:) = 0; end
    
    %make the fake wave
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %create a sinusoid packet going in 1d over a few levels
    %parameters - wave
    Amplitude  = ComboData(iCombo,8); %K
    Lambda.x   = ComboData(iCombo,10); %km - along track
    Lambda.y   = ComboData(iCombo,11); %km - across track
    Lambda.z   = ComboData(iCombo,9); %km
    Rotation   =  [ComboData(iCombo,12),ComboData(iCombo,13),ComboData(iCombo,14)];%[30,30,30]; %degrees in x,y,z
    
    dx = Spacing(1); dy = Spacing(2); dz = 3;
    
    kx = (2*pi)./(Lambda.x./dx);
    ky = (2*pi)./(Lambda.y./dy);
    kz = (2*pi)./(Lambda.z./dz);
    
    %parameters - packets (values are fwhm of a gaussian, centred at granule centre)
    Width.x =  500;
    Width.y = 1000;
    Width.z =   40;
    
    %make wave
    [x,y,z] = ndgrid(1:1:500,1:1:500,1:1:200); %will trim to size later
    
    wave = Amplitude .* sin(kx.*x + ky.*y + kz.*z);
    
    %rotate wave
    wave = imrotate(wave,Rotation(1),'bilinear','crop'); %axis 1
    wave = permute(imrotate(permute(wave,[2,1,3]),Rotation(2),'bilinear','crop'),[2,1,3]); %axis 2
    wave = permute(imrotate(permute(wave,[1,3,2]),Rotation(2),'bilinear','crop'),[1,3,2]); %axis 3
    
    %trim to size, taking the middle bit only (this avoids bits rotated off the edges)
    wave = wave((1:1:size(Airs.l1_lat,1)) - floor(mean(size(Airs.l1_lat,1))./2) + 250,...
      (1:1:size(Airs.l1_lat,2)) - floor(mean(size(Airs.l1_lat,2))./2) + 250,...
      (1:1:size(Airs.ret_z, 1)) - floor(mean(size(Airs.ret_z, 1))./2) + 100);
    
    
    %packetise
    gauss = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;
    Gauss.x = gauss(1:1:size(Airs.l1_lat,1),size(Airs.l1_lat,1)./2,Width.x./dx./2.355,1,0);
    Gauss.y = gauss(1:1:size(Airs.l1_lat,2),size(Airs.l1_lat,2)./2,Width.y./dy./2.355,1,0);
    Gauss.z = gauss(1:1:size(Airs.ret_z, 1),size( Airs.ret_z,1)./2,Width.z./dz./2.355,1,0);
    [Gx,Gy,Gz] = ndgrid(Gauss.x,Gauss.y,Gauss.z); G = Gx.*Gy.*Gz;
    Wave = wave.*G;
    
    %add to the data]
    Airs.Tp = Airs.Tp + Wave;
    
    %tidy
    clear Amplitude Lambda dx kx ky kz Width x y z wave gauss Gx Gy Gz G Rotation
    
    
    
    %analyse the wave
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    TwoDSettings.c1          = ComboData(iCombo,1).*[1,1];
    TwoDSettings.c2          = ComboData(iCombo,2).*[1,1];
    TwoDSettings.NPeaks      = ComboData(iCombo,3);
    TwoDSettings.Steps       = ComboData(iCombo,4);
    TwoDSettings.Weight      = ComboData(iCombo,5);
    
    
    %3D S-Transform the granule
    ST = gwanalyse_airs_3d(Airs,'ZRange',[0 90],'TwoDPlusOne',true,'TwoDPlusOneSettings',TwoDSettings);
    
    %store granule
    Airs = rmfield(Airs,{'ret_temp','MetaData','Source'});
    Airs.Lz  = 1./ST.m;      Airs.A  = ST.A;
    Airs.Lz2 = 1./ST.m_2dp1; Airs.A2 = ST.A_2dp1;
    
    %smooth vertically
    %     Airs.Lz  = smoothn(Airs.Lz, [3,3,3]);
    %     Airs.Lz2 = smoothn(Airs.Lz2,[1,1,3]);
    
    if iGranule == 1; Data = Airs;
    else Data = cat_struct(Data,Airs,2,{'ret_z'});
    end
  end
  
  %ipull out the bit we want
  x = (-Settings.Length:1:Settings.Length) + Case{3};
  if min(x) <0; x = x - min(x) +1; end
  y = (-Settings.Length:1:Settings.Length) + Case{4};
  if min(y) <0; y = y - min(y) +1; end
  
  Data.l1_lat  = Data.l1_lat( y,x,:);
  Data.l1_lon  = Data.l1_lon( y,x,:);
  Data.l1_time = Data.l1_time(y,x,:);
  Data.Tp      = Data.Tp(     y,x,:);
  Data.BG      = Data.BG(     y,x,:);
  Data.A       = Data.A(      y,x,:);
  Data.Lz      = Data.Lz(     y,x,:);
  Data.A2      = Data.A2(     y,x,:);
  Data.Lz2     = Data.Lz2(    y,x,:);
  
  
  %store 4 eva
  DataStore.(['Combo',num2str(iCombo)]) = Data;
  
  
  %tidy up
  clear iGranule Granules Airs Case Data dy dz Gauss iCombo Spacing Wave x y TwoDSettings ST 
  


  %done! next loop
  toc
  
end; clear iCombo

%save
save(Settings.OutFile,'DataStore','Settings','ComboData','Cases')