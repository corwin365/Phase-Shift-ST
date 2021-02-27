clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do 2D+1 analysis of a real AIRS-observed waves
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%which granules?
Settings.Date = datenum(2007,1,13);
Settings.Gs   = [121:123];

%smooth fields?
Settings.SmoothSize = [3,3,1];

%overinterpolate for plotting? %all 1s to do nothing.
%smoothsize will be multipied by this as smoothing happens after interpolation
Settings.InterpFactor = [1,1,1];

%variables to plot
%these will be interpolated and smoothed.
%all others, except geolocation, will be removed.
Settings.PlotVars = {'IN','A','F3','A_2dp1','F3_2dp1',     ...
                     'k','l','k_2dp1','l_2dp1','BG','Lon', ...
                     'Lat','Z','T','BG',                   ...
                     'MF','MFx','MFy','MF_2dp1','MFx_2dp1','MFy_2dp1'};

%z range to analyse
Settings.ZRange = [20,55];

%number of colours
Settings.NColours = 15;

%map/plot settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%height level of the maps
Settings.Plotting.ZLevel = 33; %km

%latitude and longitude box
Settings.Plotting.LatRange = [45,70];
Settings.Plotting.LonRange = [-10,25];

%cross-track row to plot
Settings.Plotting.XTRow = 55 .*Settings.InterpFactor(1);

%along track elements to plot
Settings.Plotting.ATEls = (175:1:250) .*Settings.InterpFactor(1);

%amplitude cutoff settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%we want to remove areas which are not part of the wave
%to do this, we use a smoothed amplitude field
%so we need a (2D) smoothing size *and* a cutoff value
Settings.CutOff.A = 2.4;
Settings.CutOff.S = [1,1].*5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iGranule=1:1:numel(Settings.Gs)
  
  %get granule
  [Airs,Spacing] = prep_airs_3d(Settings.Date,Settings.Gs(iGranule));
  
  %store
  if iGranule == 1;
    Store = Airs;
  else
    Store = cat_struct(Store,Airs,2,{'ret_z','MetaData','Source'});
  end
  
end
Airs = Store; 
clear iGranule Store;
disp('Data loaded')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% S-Transform. This may take some time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ST,Airs] = gwanalyse_airs_3d(Airs,'ZRange',Settings.ZRange, ...
                                   'TwoDPlusOne',true,       ...
                                   'Spacing',Spacing         );
disp('S-Transform complete')

%copy vars we want in next step, including resizing to all have same form
ST.BG       = Airs.BG;
ST.T        = Airs.ret_temp;
ST.Lat      = repmat(Airs.l1_lat,1,1,numel(Airs.ret_z));
ST.Lon      = repmat(Airs.l1_lon,1,1,numel(Airs.ret_z));
ST.Z        = permute(repmat(Airs.ret_z,1,size(ST.Lat,1),size(ST.Lat,2)),[2,3,1]);
clear Airs

%work out momentum flux (easier than in the plotting loop below)
ST.MFx = 1000.*0.5 .* cjw_airdensity(h2p(ST.Z),ST.BG) .* (9.81/0.02).^2 .* (ST.A./ST.BG).^ 2 .* ST.k./ST.m;
ST.MFy = 1000.*0.5 .* cjw_airdensity(h2p(ST.Z),ST.BG) .* (9.81/0.02).^2 .* (ST.A./ST.BG).^ 2 .* ST.l./ST.m;
ST.MF = quadadd(ST.MFx,ST.MFy);

ST.MFx_2dp1 = 1000.*0.5 .* cjw_airdensity(h2p(ST.Z),ST.BG) .* (9.81/0.02).^2 .* (ST.A_2dp1./ST.BG).^ 2 .* ST.k_2dp1./ST.m_2dp1;
ST.MFy_2dp1 = 1000.*0.5 .* cjw_airdensity(h2p(ST.Z),ST.BG) .* (9.81/0.02).^2 .* (ST.A_2dp1./ST.BG).^ 2 .* ST.l_2dp1./ST.m_2dp1;
ST.MF_2dp1 = quadadd(ST.MFx_2dp1,ST.MFy_2dp1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% smooth and or overinterpolate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar=1:1:numel(Settings.PlotVars)
  
  %get var
  Var = ST.(Settings.PlotVars{iVar});
  x = 1:1:size(Var,1);
  y = 1:1:size(Var,2);
  z = 1:1:size(Var,3);
  
  
  if sum(Settings.InterpFactor) ~= 3;
    
    %overinterpolate var
    [xi,yi,zi] = meshgrid(x,y,z);
    
    x = 1:1/Settings.InterpFactor(1):size(Var,1);
    y = 1:1/Settings.InterpFactor(2):size(Var,2);
    z = 1:1/Settings.InterpFactor(3):size(Var,3);
    [xj,yj,zj] = meshgrid(x,y,z);
    
    Var = permute(interp3(xi,yi,zi,permute(Var,[2,1,3]),xj,yj,zj),[2,1,3]);
    
  end
  
  %smooth. Do not do this for geolocation
  if ~strcmp(Settings.PlotVars{iVar},'Lon') ...
  && ~strcmp(Settings.PlotVars{iVar},'Lat') ...
  && ~strcmp(Settings.PlotVars{iVar},'  Z');
   
     Smoother = Settings.SmoothSize.*Settings.InterpFactor;
     Var = smoothn(Var,make_odd(Smoother));
  end
  
  %store
  Data.(Settings.PlotVars{iVar}) = Var;
  clear x y z xi yi zi xj yj zj Var Smoother
  
end; clear iVar
clear ST

disp('Data overinterpolated and smoothed')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load background map and topography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[~,Map,~,Coasts] = topo_etc(Settings.Plotting.LonRange+[-5,5], ...
                            Settings.Plotting.LatRange+[-5,5], ...
                            0,0,1,1);

disp('Loaded topography and maps')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear the figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.03,0.015], 0.05, [0.05 0.01]);


%we want four rows:
%1 and 3: 3DST
%2 and 4: 2D+1 ST
%
%1 and 2: horizontal planes
%3 and 4: X-Z plane
%
%to share colour tables etc, we will plot these in a weird order!
% Letters = 'abcde fghijklmn opqr';
Letters = 'abcd  efg hijk  lmno';

for iP=1:1:20;

  %individual plot settings
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %plot settings - shared for 3D and 2D+1 ST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if iP ~=1 && iP ~= 11; Settings.NColours = Settings.NColours+3; end
  switch iP
    case {1,6};  %T' map      
      Type = 'm';  Range  =[-6,6];
      Colours  = flipud(cbrewer('div','RdBu',Settings.NColours));
    case {2,7}; %amplitude map
      Type = 'm';  Range  =[0,8];
      Colours  = cbrewer('seq','Greens',Settings.NColours);
    case {4,9}; %vertical wavelength map
      Type = 'm';  Range  =[10,30];
      Colours  = cbrewer('seq','Reds',Settings.NColours); 
    case {3,8} %horizontal wavelength map
      Type = 'm';  Range  =[100,1000];
      Colours  = cbrewer('seq','Oranges',Settings.NColours);       
    case {5,10} %momentum flux map
      Type = 'm';  Range  =[0,100];
      Colours  = cbrewer('seq','Purples',Settings.NColours);   
    case {11,16};  %T' slice
      Type = 's';  Range  =[-6,6];
      Colours  = flipud(cbrewer('div','RdBu',Settings.NColours));
    case {12,17}; %amplitude slice
      Type = 's';  Range  =[0,8];
      Colours  = cbrewer('seq','Greens',Settings.NColours);
    case {14,19}; %vertical wavelength slice
      Type = 's';  Range  =[10,30];
      Colours  = cbrewer('seq','Reds',Settings.NColours); 
    case {13,18} %horizontal wavelength slice
      Type = 's';  Range  =[100,1000];
      Colours  = cbrewer('seq','Oranges',Settings.NColours);       
    case {15,20} %momentum flux slice
      Type = 's';  Range  =[0,100];
      Colours  = cbrewer('seq','Purples',Settings.NColours);   
      
    otherwise; continue
  end
  
  if iP ~=1 && iP ~= 11; 
    Settings.NColours = Settings.NColours-3;
    Colours = Colours(4:end,:);
  end
    
  %load data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  switch iP
    case  {1,11}; PlotData = Data.IN; 
    case  {2,12}; PlotData = Data.A;
    case  {3,13}; PlotData = 1./quadadd(Data.k,Data.l);
    case  {4,14}; PlotData = 1./Data.F3;
    case  {5,15}; PlotData = Data.MF;
    case  {6,16}; continue; %identical to panel 1
    case  {7,17}; PlotData = Data.A_2dp1; 
    case  {8,18}; PlotData = 1./quadadd(Data.k_2dp1,Data.l_2dp1);
    case  {9,19}; PlotData = 1./Data.F3_2dp1;
    case {10,20}; PlotData = Data.MF_2dp1;
      
    otherwise; continue
  end  
  
  
  %impose data cutoff
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  if iP ~= 1 && iP ~= 11
    
    %compute and apply cutoff
    if iP <= 5 || (iP > 10 && iP < 15); Afield = smoothn(     Data.A,[Settings.CutOff.S,1]);
    else                                Afield = smoothn(Data.A_2dp1,[Settings.CutOff.S,1]);
    end
    Bad = find(Afield < Settings.CutOff.A);
    PlotData(Bad) = 0;
    clear Afield Bad
    
    %add white to the colour table to allow this effect
    Colours = cat(1,[1,1,1],Colours);
    
    %if we're plotting amplitude,change the primary colour table too
    if iP == 2 | iP == 12 | iP == 7 | iP == 17
      Levs = linspace(Range(1),Range(2),Settings.NColours);
      Colours(Levs < Settings.CutOff.A,:) = 1;
    end
    
  end

  %general prep
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  %prepare panel
  h = subplot(4,5,iP);  
  
  %plot either a map or a slice - different code for each
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if strcmp(Type,'m'); %map

    
    
    
    
    
    
    %%%%%%%%%%%%    
    %MAP
    %%%%%%%%%%%%
    
    %choose data level
    PlotData = squeeze(PlotData(:,:,closest(Data.Z(1,1,:),Settings.Plotting.ZLevel)));
    Lon      = squeeze(Data.Lon(:,:,1));
    Lat      = squeeze(Data.Lat(:,:,1));
    
    %plot a base map
    axis([Settings.Plotting.LonRange Settings.Plotting.LatRange])
    hold on
    imagesc(Map.LonScale,Map.LatScale,Map.Map)
    
    
    %plot data over it
    pcolor(Lon,Lat,PlotData); shading flat;
    
    %plot continents over top
    for iCoast=1:1:numel(Coasts)
      plot(Coasts(iCoast).Lon,Coasts(iCoast).Lat,'-','linewi',0.5,'color',[1,1,1].*0)
      hold on
    end
    
    %plot granule edges
    plot(Lon(  1,:),Lat(  1,:),'k:','linewi',2)
    plot(Lon(end,:),Lat(end,:),'k:','linewi',2)
    
    %add a grid
    for iLon = -180:10:180; plot([1,1].*iLon,[ -90, 90],':','color',[1,1,1].*0.4); end
    for iLat =  -90:10:90;  plot([-180,180],[1,1].*iLat,':','color',[1,1,1].*0.4); end
    
    %plot a line for the cuts below
    plot(Lon(Settings.Plotting.XTRow,Settings.Plotting.ATEls), ...
         Lat(Settings.Plotting.XTRow,Settings.Plotting.ATEls), ...
         'k:','linewi',1.5)
    %make clear where the line starts!
    plot(Lon(Settings.Plotting.XTRow(1),Settings.Plotting.ATEls(1)), ...
         Lat(Settings.Plotting.XTRow(1),Settings.Plotting.ATEls(1)), ...
         'ks','markerfacecolor','k')       
    
    %polish
    set(gca,'xtick',-180:10:180,'ytick',-90:10:90)     
        
    
    
    
  elseif strcmp(Type,'s') %slice through data
        
       
    %%%%%%%%%%%%
    %SLICE
    %%%%%%%%%%%%
    
    
    
    %choose data slice and range
    PlotData = squeeze(PlotData(Settings.Plotting.XTRow,Settings.Plotting.ATEls,:));
    
    %produce an along-track distance axis
    x = Settings.Plotting.ATEls - min(Settings.Plotting.ATEls);
    x = x.*Spacing(2);
    
    %prepare axes
    axis([min(x) max(x) min(Settings.ZRange) max(Settings.ZRange)])
    hold on
    
    %plot data
    imagesc(x,squeeze(Data.Z(1,1,:)),PlotData'); shading flat
    set(gca,'ydir','normal')
    
    %plot map height
    plot(minmax(x),[1,1].*Settings.Plotting.ZLevel,'k:','linewi',1.5)
    plot(min(   x),[  1].*Settings.Plotting.ZLevel,'ks','markerfacecolor','k')
    


  end
  
  %label
  xlims = get(gca,'xlim');  ylims = get(gca,'ylim');
  text(xlims(1)+0.01.*range(xlims), ...
       ylims(1)+0.07.*range(ylims), ...
       ['(',Letters(iP),')'])
    
  %finish off
  caxis(Range)
  colormap(h,Colours);  
  colorbar
  box on
  set(gca,'tickdir','out','fontsize',10,'ticklength',[1,1].*0.04)
  
  %different axes for each group of plots - handled here
  switch iP
    case {1,2,3,4,5,6,7,8,9,10}; set(gca,'xaxislocation','top')
    otherwise;                   set(gca,'xaxislocation','bottom')
  end
  switch iP
    case { 1, 2, 3, 4, 5}; xlabel('Longitude')
    case {16,17,18,19,20}; xlabel('Distance')
    otherwise
  end
  switch iP
    case { 1, 7}; ylabel('Latitude');
    case {11,17}; ylabel('Height');
    otherwise         ylabel('        '); 
  end
  switch iP
    case {1,2,3,4,5,17,18,19,20};
    case {7,8,9,10};             set(gca,'xticklabel','     ');
    case {12,13,14,15};          set(gca,'xticklabel','        ');
    case {17,18,19,20};          set(gca,'xticklabel','    ');
  end
  switch iP
    case {2,3,4,5,8,9,10,12,13,14,15,18,19,20}; set(gca,'yticklabel',' ');
    otherwise
  end
  

  drawnow
  
   %clear variable space
   clear CLevels Colours Lat Lon Panel PlotData Range Type x 
end