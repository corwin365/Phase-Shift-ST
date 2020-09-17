clearvars

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

%cases to plot: {Date,Granules,XT row, AT row}
Cases.Case1 = {datenum(2008,1,127),1,67,45}; %we're using this one to get background, but the wave will be fake

%how many elements to plot each side of wave centre
Settings.Length = 35;

%flat-field or real noise?
Settings.Flat = 0;

%letters for labelling
Letters = 'abcdefghijklmnopqrstuvwxyz'; lettercount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% primary loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.018 0.005], [0.07,0.05], 0.05);


for iCase=1:1:numel(fieldnames(Cases))
  
  %% load data
  %%%%%%%%%%%%%%%%%%%%%%
  
  Case = Cases.(['Case',num2str(iCase)]);
  Granules = Case{2};
  Data = struct();
  
  for iGranule=1:1:numel(Granules)
    
    %load granule
    [Airs,Spacing] = prep_airs_3d(Case{1},Granules(iGranule),'PreSmooth',[3,3,1]);
    
    %keep only geolocation and zeros otherwise?
    if Settings.Flat == 1; Airs.Tp(:) = 0; end
    
    %make the fake wave
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %create a sinusoid packet going in 1d over a few levels
    %parameters - wave
    Amplitude  =           5; %K
    Lambda.x   =         600; %km - along track
    Lambda.y   =        1000; %km - across track
    Lambda.z   =          25; %km
    Rotation   =  [0,0,0];%[30,30,30]; %degrees in x,y,z
    
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
    
    %3D S-Transform the granule
    ST = gwanalyse_airs_3d(Airs,'ZRange',[15 65],'TwoDPlusOne',true);

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
  
  clear iGranule Granules Airs
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot xt and at slices of T', then AT slices of 3D ST-derived A and Lz
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for iPlot=1:1:5
    
    %prepare panel
    h = subplot(numel(fieldnames(Cases)),5,(iCase-1).*5 + iPlot);
    
    %get data and set settings
    switch iPlot
      case 1;
        ToPlot = squeeze(Data.Tp(:,Settings.Length+1,:))'; %at T'
        Colours = flipud(cbrewer('div','RdBu',33));
        Range = [-1,1].*8;
% %       case 2; 
% %         ToPlot = squeeze(Data.Tp(Settings.Length+1,:,:))'; %xt T'
% %         Colours = flipud(cbrewer('div','RdBu',33));
% %         Range = [-1,1].*8;
      case 2; 
        ToPlot = squeeze(Data.A( :,Settings.Length+1,:))'; %at A
        Colours = cbrewer('seq','Greens',33);
        Range = [0 10];
      case 3;
        ToPlot = squeeze(Data.Lz(:,Settings.Length+1,:))'; %at Lz 
        Colours = cbrewer('seq','Blues',33);
        Range = [15 35];
      case 4; 
        ToPlot = squeeze(Data.A2( :,Settings.Length+1,:))'; %at A
        Colours = cbrewer('seq','Greens',33);
        Range = [0 3];
      case 5;
        ToPlot = abs(squeeze(Data.Lz2(:,Settings.Length+1,:)))'; %at Lz
        Colours = cbrewer('seq','Blues',33);
        Range = [15 35];
    end
    

    %prepare x-axis
%     if iPlot == 2; xp = (-Settings.Length:1:Settings.Length).*Spacing(2);
%     else
      xp = (-Settings.Length:1:Settings.Length).*Spacing(1);
%     end
   
    %plot
%     contourf(xp,Data.ret_z,ToPlot,linspace(Range(1),Range(2),16),'edgecolor','none')
    imagesc(xp,Data.ret_z,ToPlot); hold on
    
    %tidy
    colormap(h,Colours)
    set(gca,'ydir','normal')
    caxis(Range)
    ylim([20 60])
  
    %labelling
    if iCase ==max(numel(fieldnames(Cases))); 
%       if iPlot ~= 2; xlabel('AT Distance [km]')
%     else
      xlabel('XT Distance [km]')
%       end
    end
    if iCase ==1;
      switch iPlot
        case 1; title('AT T''');
%         case 2; title('XT T''')
        case 2; title('3DST A')
        case 3; title('3DST \lambda z')
      end
    end
    if iPlot == 1; ylabel('Altitude [km]'); end
    if iPlot ~= 1; set(gca,'yticklabel',{});end
    if iCase ~= max(numel(fieldnames(Cases))); set(gca,'xticklabel',{});end;
    
    %lettering
    lettercount = lettercount+1;
    text(xp(1),24,['(',Letters(lettercount),')'],'fontweight','bold','fontsize',18)
    
    %done!
    colorbar
    drawnow
  end
  clear iPlot x xo y ToPlot
  

end; clear iCase