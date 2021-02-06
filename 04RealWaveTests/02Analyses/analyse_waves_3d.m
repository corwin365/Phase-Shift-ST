clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do 2D+1 analysis of some real AIRS-observed waves
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%smooth fields?
Settings.SmoothSize = [5,5,1];

%overinterpolate for plotting? %all 1s to do nothing.
Settings.InterpFactor = [3,3,3];

%variables to plot
Settings.PlotVars = {'IN','A','F3','A_2dp1','F3_2dp1'};

%z range to plot
Settings.ZRange = [20,60];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare figure, and go
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clf
set(gcf,'color','w')

for iWave = 1;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% load waves
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch iWave
    case 1;
      %Andes big 2008
      [Airs,Spacing] = prep_airs_3d(datenum(2008,1,127),57);
      CutLine        = 40;
      CutRange       = [1:1:90];
    case 2;
      %Scandinavia big
      [Airs,Spacing] = prep_airs_3d(datenum(2007,1,13),122);
      CutLine        = 80;
    case 3;
      %midwest convective
      [Airs,Spacing] = prep_airs_3d(datenum(2005,7,8),85);
      Airs2 = prep_airs_3d(datenum(2005,7,8),86);
      Airs.l1_lat = cat(2,Airs.l1_lat,Airs2.l1_lat);
      Airs.l1_lon = cat(2,Airs.l1_lon,Airs2.l1_lon);
      Airs.Tp     = cat(2,Airs.Tp,    Airs2.Tp);
      clear Airs2
      CutLine     = 110;
    otherwise
      stop
  end

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% analyse wave with 2D+1 and 3D ST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  [ST,Airs] = gwanalyse_airs_3d(Airs,'ZRange',Settings.ZRange, ...
    'TwoDPlusOne',true,       ...
    'Spacing',Spacing         );
  clear Spacing
  
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
      
      x = 1:1./Settings.InterpFactor(1):size(Var,1);
      y = 1:1./Settings.InterpFactor(2):size(Var,2);
      z = 1:1./Settings.InterpFactor(3):size(Var,3);
      [xj,yj,zj] = meshgrid(x,y,z);
      
      Var = permute(interp3(xi,yi,zi,permute(Var,[2,1,3]),xj,yj,zj),[2,1,3]);
      
    end
    
    %smooth
    Smoother = Settings.SmoothSize.*Settings.InterpFactor;
    Var = smoothn(Var,make_odd(Smoother));
    
    %store data
    Plot.(Settings.PlotVars{iVar}) =  Var;
    Plot.z = linspace(min(Airs.ret_z),max(Airs.ret_z),numel(z));
    Plot.x = x;
    Plot.y = y;
    
  end; clear xi yi zi Smoother Var xy jy zj x y z
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plotting - 2d cuts through the data
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %identify plotting location after over-interpolation
  x = CutLine.*Settings.InterpFactor(1);
  y = CutRange.*Settings.InterpFactor(1);
  
  
  %% input wave
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subplot(3,5,(3.*(iWave-1))+1)
  
  %get var
  Var = Plot.(Settings.PlotVars{1});
  
  %subset in space
  ToPlot = squeeze(Var(x,y,:));
    
  contourf(Plot.y(CutRange),Plot.z,ToPlot',-15:1.25:15,'edgecolor','none');
  colormap(flipud(cbrewer('div','RdBu',25)))
  caxis([-1,1].*15)
  hold on
  
  LineLevs = -20:4:20; LineLevs(LineLevs == 0) = [];
  [c,h] = contour(Plot.y(CutRange),Plot.z,ToPlot',LineLevs,'edgecolor','k');
  clabel(c,h)
  
  %% 3D wave amplitude
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subplot(3,5,(3.*(iWave-1))+2)
  
  %get var
  Var = Plot.(Settings.PlotVars{2});
  
  %subset in space
  ToPlot = squeeze(Var(x,y,:));
    
  contourf(Plot.y(CutRange),Plot.z,ToPlot',0:1.25:15,'edgecolor','none');
  colormap(flipud(cbrewer('div','RdBu',25)))
  caxis([0 15])
  hold on
  
  LineLevs = -20:4:20; LineLevs(LineLevs == 0) = [];
  [c,h] = contour(Plot.y(CutRange),Plot.z,ToPlot',LineLevs,'edgecolor','k');
  clabel(c,h)
  colorbar
  
   
  %% 2D+1 wave amplitude
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subplot(3,5,(3.*(iWave-1))+3)
  
  %get var
  Var = Plot.(Settings.PlotVars{4});
  
  %subset in space
  ToPlot = squeeze(Var(x,y,:));
    
  contourf(Plot.y(CutRange),Plot.z,ToPlot',0:1.25:15,'edgecolor','none');
  colormap(flipud(cbrewer('div','RdBu',25)))
  caxis(ColourRange)
  hold on
  
  LineLevs = -20:4:20; LineLevs(LineLevs == 0) = [];
  [c,h] = contour(Plot.y(CutRange),Plot.z,ToPlot',LineLevs,'edgecolor','k');
  clabel(c,h)
  colorbar 
  
  
  %% 3D wavelength
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  subplot(3,5,(3.*(iWave-1))+4)
  
  %get var
  Var = 1./Plot.(Settings.PlotVars{3});
  
  %subset in space
  ToPlot = squeeze(Var(x,y,:));
    
  contourf(Plot.y(CutRange),Plot.z,ToPlot',10:1:30,'edgecolor','none');
  colormap(flipud(cbrewer('div','RdBu',25)))
  caxis([10 30])
  hold on
  
  LineLevs = -20:4:20; LineLevs(LineLevs == 0) = [];
  [c,h] = contour(Plot.y(CutRange),Plot.z,ToPlot',LineLevs,'edgecolor','k');
  clabel(c,h)
  colorbar 

  
  
  
  stop
  
end