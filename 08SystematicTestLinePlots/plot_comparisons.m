clearvars

%ST tests carried out:
%for the basic test wave, varying as f(m):
%1. thinning makes essentially no differences
%2. height-weighting does make some difference - test when other settings pinned
%3. a step of 2 is better than a step of 3
%4. c1 = 1 and c2 = 0.25 works pretty well.
%5. thresholds above zero lose the short-vertical waves


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot systematic ideal wave variations with 3DST and 2D+1
%
% Corwin Wright, c.wright@bath.ac.uk, 2020/10/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%input file
% Settings.InFile = 'out_systematic.mat';
Settings.InFile = 'out_testing2.mat';
%average type we want to plot
Settings.AvType = 2;

%coding of lines
 %line 1 - 3D, flat-field
Settings.Line1.Colour = 'k';
Settings.Line1.Marker = '^';
 %line 2 - 3D, noise-field
Settings.Line2.Colour = 'k';
Settings.Line2.Marker = 'o';
 %line 3 - 2D+1, flat-field
Settings.Line3.Colour = 'r';
Settings.Line3.Marker = '^';
 %line 4 - 2D+1, noise-field
Settings.Line4.Colour = 'r'; 
Settings.Line4.Marker = 'o';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
Data = load(Settings.InFile);

%get list of the variables we saved
PlotVars = Data.Settings.OutVars;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over vars, and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Variables = fieldnames(Data.AllResults);

for iVar=1:1:numel(Variables)
  
  %create figure
  figure(iVar+20)
  clf
  set(gcf,'color','w')

  %information about the variable analysed
  VarName = Variables{iVar};
  XScale  = Data.Settings.Wave.(VarName);
  Output  = Data.AllResults.(Variables{iVar});  
  
  %sort data into order on x-axis, as first point is special in generator programme
  [~,idx] = sort(XScale,'asc');
  XScale = XScale(idx);
  Output = Output(:,:,:,:,idx);
   clear idx
  
  %loop over plotted variables
  for jVar=1:1:numel(PlotVars)
    
    %pull out data to plot
    ToPlot = squeeze(Output(:,:,jVar,Settings.AvType,:));
    
    %if the data is a wavenumber, flip it
    if any(strcmp(PlotVars{jVar},{'kh','m'})); ToPlot = 1./ToPlot; end
       
    %create panel
    subplot(2,ceil(numel(PlotVars)./2),jVar)
    
    %create axes
    axis([min(XScale) max(XScale) min(ToPlot(:)) max(ToPlot(:))])
    axis square
    hold on
    box on
    
    %label them
    xlabel(VarName)
    ylabel(PlotVars{jVar})
    
    %plot data
    plot(XScale,squeeze(ToPlot(1,1,:)), ...
         'color',Settings.Line1.Colour,'markerfacecolor',Settings.Line1.Colour, ...
         'marker',Settings.Line1.Marker)
    plot(XScale,squeeze(ToPlot(2,1,:)), ...
         'color',Settings.Line2.Colour,'markerfacecolor',Settings.Line2.Colour, ...
         'marker',Settings.Line2.Marker)    
    
    plot(XScale,squeeze(ToPlot(1,2,:)), ...
         'color',Settings.Line3.Colour,'markerfacecolor',Settings.Line3.Colour, ...
         'marker',Settings.Line3.Marker)
    plot(XScale,squeeze(ToPlot(2,2,:)), ...
         'color',Settings.Line4.Colour,'markerfacecolor',Settings.Line4.Colour, ...
         'marker',Settings.Line4.Marker)         
       
  end; clear jVar
  
  
  
  
  
  
end; clear iVar
