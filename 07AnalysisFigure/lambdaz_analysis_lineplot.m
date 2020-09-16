clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% figure showing the 2D+1 analysis for four waves as 1D lines
%
% Corwin Wright, c.wright@bath.ac.uk, 2020/08/24
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cases to plot: {Date,Granules,XT row, AT row}
Cases.Case1 = {datenum(2008,1,127),[57],30,45};
Cases.Case2 = {datenum(2007,1,13),[122],80,55};
Cases.Case3 = {datenum(2005,7,8),[85,86],110,25};
% Cases.Case4 = Cases.Case1; %placeholder


%letters for labelling
Letters = 'abcdefghijklmnopqrstuvwxyz'; lettercount = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% primary loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%prepare figure
clf
set(gcf,'color','w')
subplot = @(m,n,p) subtightplot (m, n, p, [0.018 0.005], [0.07,0.05], 0.05);


for iCase=1:1:3%numel(fieldnames(Cases))
  
  %% load data
  %%%%%%%%%%%%%%%%%%%%%%
  
  Case = Cases.(['Case',num2str(iCase)]);
  Granules = Case{2};
  Data = struct();
  
  for iGranule=1:1:numel(Granules)
    
    %load granule
    [Airs,Spacing] = prep_airs_3d(Case{1},Granules(iGranule),'PreSmooth',[3,3,1]);
    
    %3D S-Transform the granule
    ST = gwanalyse_airs_3d(Airs,'ZRange',[0 90],'TwoDPlusOne',true);

    %store granule
    Airs = rmfield(Airs,{'ret_temp','MetaData','Source'}); 
    Airs.Lz  = 1./ST.m;      Airs.A  = ST.A;
    Airs.Lz2 = 1./ST.m_2dp1; Airs.A2 = ST.A_2dp1;     

% % %     %smooth vertically
% % %     Airs.Lz  = smoothn(Airs.Lz, [1,1,5]);
% %     Airs.Lz2 = smoothn(Airs.Lz2,[1,1,3]);
    
    if iGranule == 1; Data = Airs; 
    else Data = cat_struct(Data,Airs,2,{'ret_z'});
    end
  end
  
  %ipull out the bit we want
  x = Case{3};
  y = Case{4};
  
  Data.l1_lat  = Data.l1_lat( y,x,:);
  Data.l1_lon  = Data.l1_lon( y,x,:);
  Data.l1_time = Data.l1_time(y,x,:);
  Data.Tp      = Data.Tp(     y,x,:);
  Data.BG      = Data.BG(     y,x,:);
  Data.A       = Data.A(      y,x,:);
  Data.A2      = Data.A2(     y,x,:);
  Data.Lz      = Data.Lz(     y,x,:);
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
        ToPlot = squeeze(Data.Tp)'; %at T'
        Colours = flipud(cbrewer('div','RdBu',33));
        Range = [-1,1].*8;
      case 2; 
        ToPlot = squeeze(Data.A)'; %at A
        Range = [0 10];
      case 3;
        ToPlot = squeeze(Data.Lz)'; %at Lz 
        Range = [0 50];
      case 4; 
        ToPlot = squeeze(Data.A2)'; %at A
        Range = [0 3];
      case 5;
        ToPlot = abs(squeeze(Data.Lz2))'; %at Lz
        Range = [0 50];
    end

    %plot
    plot(ToPlot,Data.ret_z); hold on
    
    %tidy
    xlim(Range)
    ylim([20 60])
  
% % %     %labelling
% % %     if iCase ==max(numel(fieldnames(Cases))); 
% % % %       if iPlot ~= 2; xlabel('AT Distance [km]')
% % % %     else
% % %       xlabel('XT Distance [km]')
% % % %       end
% % %     end
% % %     if iCase ==1;
% % %       switch iPlot
% % %         case 1; title('AT T''');
% % % %         case 2; title('XT T''')
% % %         case 2; title('3DST A')
% % %         case 3; title('3DST \lambda z')
% % %       end
% % %     end
% % %     if iPlot == 1; ylabel('Altitude [km]'); end
% % %     if iPlot ~= 1; set(gca,'yticklabel',{});end
% % %     if iCase ~= max(numel(fieldnames(Cases))); set(gca,'xticklabel',{});end;
    
    %lettering
    lettercount = lettercount+1;
%     text(xp(1),24,['(',Letters(lettercount),')'],'fontweight','bold','fontsize',18)
    
    %done!
    drawnow
  end
  clear iPlot x xo y ToPlot
  

end; clear iCase