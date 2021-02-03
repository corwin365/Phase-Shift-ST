clearvars


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
% Settings.InFile = 'out_incfrac.mat';
Settings.InFile = 'out_testing2.mat';

%average type we want to plot
Settings.AvType = 2; %1 mean, 2 median, 3 mode


%coding of lines
 %line 1 - 3D, flat-field
Settings.Line1.Colour = [239, 86, 45]./255;
Settings.Line1.Marker = 'o';%'^';
 %line 2 - 3D, noise-field
Settings.Line2.Colour = [12, 76, 138]./255;
Settings.Line2.Marker = 'o';%'o';

 %line 3 - 2D+1, flat-field
Settings.Line3.Colour = [246, 210, 87]./255;
Settings.Line3.Marker = '^';%'s';
 %line 4 - 2D+1, noise-field
Settings.Line4.Colour = [136, 177, 75]./255; 
Settings.Line4.Marker = '^';%'d';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data, flip wavenumbers to wavelengths, rename variables to
%human-readable, and remove unwanted vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load data
Data = load(Settings.InFile);

% % %to reuse old code, overwrite the 'standard' variables with the
% % %'fractional' variables
% % Data.AllResults = Data.AllFracResults;

%get list of the input variables
Variables = fieldnames(Data.AllResults);

%get list of the output variables
PlotVars = Data.Settings.OutVars; 



%remove any input variables we don't want to plot
ToRemove = {'Rotationy'};%,'Lambday'};
List = [];
for iVar=1:1:numel(Variables)
  if ~any(strcmp(Variables{iVar},ToRemove));
    List(end+1) = iVar;
  end
end
Variables = Variables(List);
clear List

%variable flipping and renaming
for iVar=1:1:numel(PlotVars)
  if strcmp(PlotVars{iVar},'kh') |  strcmp(PlotVars{iVar},'F3')
    for jVar=1:1:numel(Variables)
      a = Data.AllResults.(Variables{jVar});
      a(:,:,iVar,:,:) = 1./a(:,:,iVar,:,:); 
      Data.AllResults.(Variables{jVar}) = a;
    end
  end
  if strcmp(PlotVars{iVar}, 'A'); PlotVars{iVar} = 'T'' [K]';       end
  if strcmp(PlotVars{iVar},'kh'); PlotVars{iVar} = '\lambda_h [km]'; end
  if strcmp(PlotVars{iVar},'F3'); PlotVars{iVar} = '\lambda_z [km]'; end
  if strcmp(PlotVars{iVar},'th'); PlotVars{iVar} = 'Propagation angle [deg]';    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% scale data to be fractional, then set plotting ranges
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% % % %find the range of each of the
% % % %output variables
% % % Ranges = NaN(2,numel(PlotVars)); Ranges(1,:) = 99e99; Ranges(2,:) = -99e99;
% % % for iV=1:1:numel(Variables)
% % %   a = Data.AllResults.(Variables{iV});
% % %   for iOut=1:1:numel(PlotVars);
% % %     mm = minmax(a(:,:,iOut,:,:));
% % %     if mm(1) < Ranges(1,iOut); Ranges(1,iOut) = mm(1);end
% % %     if mm(2) > Ranges(2,iOut); Ranges(2,iOut) = mm(2);end
% % %   end  
% % %   
% % % end; clear iV a iOut mm

%set manually for final version
Ranges(:,1) = [0,4.5];
Ranges(:,2) = [0 100];
Ranges(:,3) = [100 600];
Ranges(:,4) = [0 180];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over vars, and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare plotting positions
k = 0;
Positions = 1:1:numel(Variables).*numel(PlotVars);
Positions = reshape(Positions,numel(Variables),numel(PlotVars))';


%create figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.04 0.025], 0.06, [0.05,0.25]);
Letters = 'abcdefghijklmnopqrstuvwxyz';

%loop over INPUT variables
for iVar=1:1:numel(Variables)

  %information about the variable analysed
  VarName = Variables{iVar};
  XScale  = Data.Settings.Wave.(VarName);
  Output  = Data.AllResults.(Variables{iVar});  
  
  %sort data into order on x-axis, as first point is special in generator programme
  Basis = XScale(1);
  [~,idx] = sort(XScale,'asc');
  XScale = XScale(idx);
  Output = Output(:,:,:,:,idx);
  clear idx
  
  %find the name of the input variable
  switch VarName
    case 'Lambdax';   OutName = 'Along-track \lambda [km]';
    case 'Lambday';   OutName = 'Across-track \lambda [km]';
    case 'Lambdaz';   OutName = 'Vertical \lambda [km]';
    case 'Rotationx'; OutName = 'X-Rotation [deg]';      
    case 'Rotationy'; OutName = 'Y-Rotation [deg]';      
    case 'Rotationz'; OutName = 'Z-Rotation [deg]';      
    case 'Amplitude'; OutName = 'Amplitude [K]';   
    otherwise; OutName = VarName;
  end
  
  
  %loop over jOUTPUT variables
  for jVar=1:1:numel(PlotVars)
    
    %pull out data to plot
    ToPlot = squeeze(Output(:,:,jVar,Settings.AvType,:));
    
    %create panel
    k = k+1;
    subplot(numel(PlotVars),numel(Variables),Positions(k))
    
    %create axes
    axis([0 max(XScale) Ranges(:,jVar)'])
%     axis square
    hold on
    box on
    
    %label them
    if jVar == 1;                   xlabel(OutName); set(gca,'xaxislocation','top'); 
    elseif jVar == numel(PlotVars); xlabel(OutName);
    else;                           set(gca,'xticklabel',[]);
    end
    
    if     iVar == 1;                ylabel(PlotVars{jVar}); 
    elseif iVar == numel(Variables); ylabel(PlotVars{jVar}); set(gca,'yaxislocation','right'); 
    else          set(gca,'yticklabel',[]); 
    end
    
    
    %plot data
    plot(XScale,squeeze(ToPlot(1,1,:)), ...
         'color',Settings.Line1.Colour,'markerfacecolor',Settings.Line1.Colour, ...
         'marker',Settings.Line1.Marker,'linewi',1,'markersize',3)  
    plot(XScale,squeeze(ToPlot(2,1,:)), ...
         'color',Settings.Line2.Colour,'markerfacecolor',Settings.Line2.Colour, ...
         'marker',Settings.Line2.Marker,'linewi',1,'markersize',3)     
    
    plot(XScale,squeeze(ToPlot(1,2,:)), ...
         'color',Settings.Line3.Colour,'markerfacecolor',Settings.Line3.Colour, ...
         'marker',Settings.Line3.Marker,'linewi',1,'markersize',3)  
    plot(XScale,squeeze(ToPlot(2,2,:)), ...
         'color',Settings.Line4.Colour,'markerfacecolor',Settings.Line4.Colour, ...
         'marker',Settings.Line4.Marker,'linewi',1,'markersize',3)         
    
       
       
   %tidy up
   set(gca,'tickdir','out')
   set(gca,'color',[1,1,1].*0.9)
   set(gca,'TickLength',[1, 1].*0.05);
   text(0.05.*max(XScale),Ranges(1,jVar)+0.9.*range(Ranges(:,jVar)),['(',Letters(Positions(k)),')'])
       
  end; clear jVar
  
  

  
  drawnow
  
  
  
  
  
end; clear iVar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot key (relative to last panel)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xlims = get(gca,'xlim'); ylims = get(gca,'ylim');

plot(max(xlims)+[0.55,0.9].*range(xlims), min(ylims)+[1,1].*0.75.*range(ylims), ...
     'color',Settings.Line1.Colour,'markerfacecolor',Settings.Line1.Colour, ...
     'marker',Settings.Line1.Marker,'linewi',1,'clipping','off','markersize',3)  
text(max(xlims)+[1].*range(xlims),min(ylims)+0.75.*range(ylims),'3DST, no noise')

plot(max(xlims)+[0.55,0.9].*range(xlims), min(ylims)+[1,1].*0.60.*range(ylims), ...
     'color',Settings.Line2.Colour,'markerfacecolor',Settings.Line2.Colour, ...
     'marker',Settings.Line2.Marker,'linewi',1,'clipping','off','markersize',3)  
text(max(xlims)+[1].*range(xlims),min(ylims)+0.60.*range(ylims),'3DST, with noise')

plot(max(xlims)+[0.55,0.9].*range(xlims), min(ylims)+[1,1].*0.45.*range(ylims), ...
     'color',Settings.Line3.Colour,'markerfacecolor',Settings.Line3.Colour, ...
     'marker',Settings.Line3.Marker,'linewi',1,'clipping','off','markersize',3)  
text(max(xlims)+[1].*range(xlims),min(ylims)+0.45.*range(ylims),'2D+1 ST, no noise')

plot(max(xlims)+[0.55,0.9].*range(xlims), min(ylims)+[1,1].*0.30.*range(ylims), ...
     'color',Settings.Line4.Colour,'markerfacecolor',Settings.Line4.Colour, ...
     'marker',Settings.Line4.Marker,'linewi',1,'clipping','off','markersize',3)     
text(max(xlims)+[1].*range(xlims),min(ylims)+0.30.*range(ylims),'2D+1 ST, with noise')
   