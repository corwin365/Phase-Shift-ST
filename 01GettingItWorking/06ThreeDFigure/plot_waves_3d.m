clear all
clf
set(gcf,'color','w')
set(gca,'fontsize',16)

for iWave = 2%1:1:3

  
  set(gcf,'position',[329.4,121.8,868,640]) %for size consistency
  cla
  hold on
  box on
  
  switch iWave
    case 1; 
      %Andes big 2008
      Airs = prep_airs_3d(datenum(2008,1,127),57);
      LatRange    = [-61,-41];
      LonRange    = [-80,-57];
      HeightRange = [20,60];
      ViewAngle   = [-36,28];
      IsoSurfaces = 4.5;%4:1:10;
      Alphas      = 0.9;%ones(size(IsoSurfaces)).*0.25;
      SmoothSize  = [5,5,1];   
    case 2;
      %Scandinavia big
      Airs = prep_airs_3d(datenum(2007,1,13),122);
      LatRange    = [45,70];
      LonRange    = [0,25];
      HeightRange = [20,52];
      ViewAngle   = [-60 16];
      IsoSurfaces = 1.5;%4:1:10;
      Alphas      = 0.9;%ones(size(IsoSurfaces)).*0.25;
      SmoothSize  = [5,5,1];               
    case 3; 
      %midwest convective
      Airs = prep_airs_3d(datenum(2005,7,8),85);
      Airs2 = prep_airs_3d(datenum(2005,7,8),86);
      Airs.l1_lat = cat(2,Airs.l1_lat,Airs2.l1_lat);
      Airs.l1_lon = cat(2,Airs.l1_lon,Airs2.l1_lon);      
      Airs.Tp     = cat(2,Airs.Tp,    Airs2.Tp);  
      clear Airs2
      LatRange    = [30,65];
      LonRange    = [-100,-80];
      HeightRange = [30,50];
      ViewAngle   = [155 50];
      IsoSurfaces = 1.5;%4:1:10;
      Alphas      = 0.9;%ones(size(IsoSurfaces)).*0.25;
      SmoothSize  = [3,3,1];             
    otherwise
      stop
  end
  
  
  stop
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%reformat for isosurface plotting
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Plot.Lon  = repmat(Airs.l1_lon,1,1,27);
  Plot.Lat  = repmat(Airs.l1_lat,1,1,27);
  Plot.Z    = repmat(permute(Airs.ret_z,[2,3,1]),size(Airs.l1_lon,1),size(Airs.l1_lon,2),1);
  Plot.Data = Airs.Tp;
  clear Airs
  
  
  %remove bad data at top
  Good = find(squeeze(nanmean(nanmean(Plot.Z))) >= HeightRange(1) ...
            & squeeze(nanmean(nanmean(Plot.Z))) <= HeightRange(2));
  Plot.Z    = Plot.Z(   :,:,Good);
  Plot.Lon  = Plot.Lon( :,:,Good);
  Plot.Lat  = Plot.Lat( :,:,Good);
  Plot.Data = Plot.Data(:,:,Good);
  clear Good HeightRange

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load and plot geolocation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %load topography
  [Topo,Map,~,~] = topo_etc(LonRange,LatRange,0,0,0,1);
  
  %smooth and vertically scale it a bit
  Topo.elev = smoothn(Topo.elev,[1,1].*3).*2;
  
  %axes
  axis([LonRange(1),LonRange(2),LatRange(1),LatRange(2),-.1,60])
  set(gca,'ytick',[-90:5:90],  'yticklabel',{'90S','85S','80S','75S','70S','65S','60S','55S','50S','45S','40S','35S','30S','25S','20S','15S','10S','5S','EQ','5N','10N','15N','20N','25N','30N','35N','40N','45N','50N','55N','60N','65N','70N','75N','80N','85N','90N'})
  set(gca,'xtick',[-180:5:180],'xticklabel',{'180W','175W','170W','165W','160W','155W','150W','145W','140W','135W','130W','125W','120W','115W','110W','105W','100W','95W','90W','85W','80W','75W','70W','65W','60W','55W','50W','45W','40W','35W','30W','25W','20W','15W','10W','5W','0W','5E','10E','15E','20E','25E','30E','35E','40E','45E','50E','55E','60E','65E','70E','75E','80E','85E','90E','95E','100E','105E','110E','115E','120E','125E','130E','135E','140E','145E','150E','155E','160E','165E','170E','175E','180E'})
  set(gca,'ztick',[0:20:60],   'zticklabel',{'0km','20km','40km','60km'})
  
  %prepare panel and plot terrain  
  hMap = surface(Topo.lons,Topo.lats,Topo.elev,Map.Map,...
    'FaceColor','texturemap','EdgeColor','none','CDataMapping','direct');
  set(hMap,'clipping','on')
  
  %set terrain to not reflect light specularly
  set(hMap,'DiffuseStrength',1,'SpecularStrength',0)
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %smooth
  Plot.Data = smoothn(Plot.Data,SmoothSize);
  
  %replace top and bottom layer with zeroes, to close isosurfaces
  Plot.Data(:,:,  1) = 0;    
  Plot.Data(:,:,end) = 0;  
  

  for iMode=1:1:2;
    
    if iMode == 1; %positive
      Colours = flipud(cbrewer('seq','Reds',numel(IsoSurfaces).*2));
    elseif iMode == 2; %negative
      Colours = flipud(cbrewer('seq','Blues',numel(IsoSurfaces).*2));
      IsoSurfaces = -1.*IsoSurfaces;
    end
    
    for iIso=1:1:numel(IsoSurfaces);
      %create unique ID for patch
      PatchId = ['Patch',num2str(iIso+iMode.*100)];
      
      %produce patch
      fv = isosurface(Plot.Lon,Plot.Lat,Plot.Z,Plot.Data,IsoSurfaces(iIso));
      eval([PatchId,' = patch(fv);'])
      
      % plot as transparent surfaces
      eval([PatchId,'.FaceColor = Colours(iIso,:);']); %
      eval([PatchId,'.EdgeColor = ''none'';']);
      eval([PatchId,'.FaceAlpha  = Alphas(iIso);']);
      
%       %fix reflection
%       eval(['set(',PatchId,',''DiffuseStrength'',1,''SpecularStrength'',0.5,''SpecularExponent'',15)'])
      
      
      drawnow
    end
    
  end
  


  camlight left; camlight  
  lighting gouraud
  view(ViewAngle)  
  axis square
  set(gca,'Projection','perspective')
  drawnow

  
  
  clearvars -except iWave
end