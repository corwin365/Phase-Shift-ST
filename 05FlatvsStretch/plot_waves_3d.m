%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot showing wave to working scale and real scale
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/14
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set(gcf,'position',[329.4,121.8,868,640],'color','w') %for size consistency
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = 198;
Airs = prep_airs_3d(datenum(2010,6,28),a-1);
Airs = cat_struct(Airs, ...
                  prep_airs_3d(datenum(2010,6,28),a), ...
                  2,{'ret_z','MetaData','Source'});
Airs = cat_struct(Airs, ...
                  prep_airs_3d(datenum(2010,6,28),a+1), ...
                  2,{'ret_z','MetaData','Source'});                

                
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
%% plotting settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LatRange    = [-75,-40];
LonRange    = [55,90];
HeightRange = [20,60];
ViewAngle   = [-82,40];
IsoSurfaces = 3.5;
Alphas      = 0.7;
SmoothSize  = [3,3,3];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% reformat for isosurface plotting
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

%smooth
Plot.Data = smoothn(Plot.Data,SmoothSize);

%replace top and bottom layer with zeroes, to close isosurfaces
Plot.Data(:,:,  1) = 0;
Plot.Data(:,:,end) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load geolocation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load topography
[Topo,Map,~,~] = topo_etc(LonRange,LatRange,0,0,0,1);

%smooth and vertically scale it a bit
Topo.elev = smoothn(Topo.elev,[1,1].*3).*2;



for iPlot=1:1:2;
  
  %prepare panel
  subplot(1,2,iPlot)

  if iPlot==2;
    %degree -> km conversion, but inverted for height
    %1 degree is about 112km
    Plot.Z    = Plot.Z    ./ 112;
    Topo.elev = Topo.elev ./112;
    
    %then scale to the panel
    Plot.Z    = Plot.Z    .*60./range(LonRange);
    Topo.elev = Topo.elev .*60./range(LonRange);
    
    Plot.Data = Plot.Data.*-1; %getting flipped somewhere as is clear from the output - this will "fix" it 
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot geolocation
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %axes
  axis([LonRange(1),LonRange(2),LatRange(1),LatRange(2),-.1,60])
  set(gca,'ytick',[-90:10:90],  'yticklabel',{'90S','85S','80S','75S','70S','65S','60S','55S','50S','45S','40S','35S','30S','25S','20S','15S','10S','5S','EQ','5N','10N','15N','20N','25N','30N','35N','40N','45N','50N','55N','60N','65N','70N','75N','80N','85N','90N'})
  set(gca,'xtick',[-80:10:180],'xticklabel',{'180W','175W','170W','165W','160W','155W','150W','145W','140W','135W','130W','125W','120W','115W','110W','105W','100W','95W','90W','85W','80W','75W','70W','65W','60W','55W','50W','45W','40W','35W','30W','25W','20W','15W','10W','5W','0W','5E','10E','15E','20E','25E','30E','35E','40E','45E','50E','55E','60E','65E','70E','75E','80E','85E','90E','95E','100E','105E','110E','115E','120E','125E','130E','135E','140E','145E','150E','155E','160E','165E','170E','175E','180E'})
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
      
      drawnow
    end
    
  end
  
  if iPlot==1;
    set(gca,'ztick',0:20:60)
  else
    set(gca,'ztick',[0,60]./112.*(60./range(LonRange)),'zticklabel',{' ','60km'})
    grid off
  end
  
  
  camlight left; camlight
  lighting gouraud
  view(ViewAngle)
  axis square
  set(gca,'Projection','perspective')
  drawnow
  
  
  
end
