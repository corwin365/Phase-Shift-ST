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
Cases.Case1 = {datenum(2008,1,127),[57],35,45};
Cases.Case2 = {datenum(2007,1,13),[122],80,55};
Cases.Case3 = {datenum(2005,7,8),[85,86],110,25};
Cases.Case4 = Cases.Case1; %placeholder

%how many elements to plot each side of wave centre
Settings.Length = 30;

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
    
    %3D S-Transform the granule
    ST = gwanalyse_airs_3d(Airs,'ZRange',[0 90]);

    %store granule
    Airs = rmfield(Airs,{'ret_temp','MetaData','Source'}); 
    Airs.Lz = 1./ST.m; Airs.A = ST.A;
    
    
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
        Range = [15 25];
      case 4; continue
      case 5; continue
    end
    

    %prepare x-axis
%     if iPlot == 2; xp = (-Settings.Length:1:Settings.Length).*Spacing(2);
%     else
      xp = (-Settings.Length:1:Settings.Length).*Spacing(1);
%     end
   
    %plot
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
    drawnow
  end
  clear iPlot x xo y ToPlot
  

end; clear iCase