
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot data components of conceptual diagram
%to be stitched together in PPT
%
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%granules
Settings.Date     = datenum(2008,1,127);
Settings.Granules = 56:1:58;

%smoothing the data?
Settings.SmoothSize = [13,13,1];

%altitudes to plot
Settings.Z        = [28:7:60];

%exact range of data to plot (in along-track rows)
Settings.PlotRange = [70:1:250];

%colour and line levels for the plots
Settings.CLevs = -5:0.25:5;
Settings.LLevs = -99:2:100; 
Settings.LLevs(abs(Settings.LLevs) < 1) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load and stitch together data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iG=1:1:numel(Settings.Granules)
  
  Airs = prep_airs_3d(Settings.Date,Settings.Granules(iG));
  
  if iG == 1;
    Store.Tp    = Airs.Tp;
    Store.ret_z = Airs.ret_z;
  else
    Store.Tp = cat(2,Store.Tp,Airs.Tp);
  end
  
  
end; clear iG Airs
Airs = Store; clear Store

%smooth (after loading, to avoid inter-granule discontinuities)
Airs.Tp = smoothn(Airs.Tp,Settings.SmoothSize);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % % for iLevel=1:1:numel(Settings.Z)
% % % 
% % %   %prepare panel
% % %   figure(iLevel)
% % %   clf
% % %   set(gcf,'color','w')
% % %   
% % %  
% % %   %colour plot
% % %   ToPlot = Airs.Tp(:,Settings.PlotRange,closest(Airs.ret_z,Settings.Z(iLevel)));
% % %   ToPlot(ToPlot < min(Settings.CLevs)) = min(Settings.CLevs); 
% % %   contourf(ToPlot, ...
% % %            Settings.CLevs,'edgecolor','none')
% % %   caxis(minmax(Settings.CLevs))
% % %   colormap(flipud(cbrewer('div','RdBu',numel(Settings.CLevs))))  
% % %   hold on
% % %   
% % %   
% % %   %line plot
% % %   [c,h] = contour(ToPlot,Settings.LLevs,'k-');
% % %   clabel(c,h)
% % %   
% % %   %label
% % %   text(5,80,['z=',num2str(Settings.Z(iLevel)),'km'],'fontsize',15)
% % %   
% % %   
% % %   %tidy up and export
% % %   set(gca,'xtick',-30:30:1000,'ytick',-0:30:1000);%,'xticklabel',{},'yticklabel',{})
% % %   axis equal
% % %   grid on
% % %   drawnow
% % %   
% % %   
% % % %   export_fig(['Tp_',num2str(Settings.Z(iLevel)),'km'],'-png','-m2','-a4')
% % %   
% % %   
% % %   
% % % 
% % %   
% % % end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 2D STs of levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for iLevel=1:1:numel(Settings.Z)-1

  %prepare panel
  figure(iLevel)
  clf
  set(gcf,'color','w')
  
  %2D ST the level and the level below
  DataA = Airs.Tp(:,:,closest(Airs.ret_z,Settings.Z(iLevel)));
  STA = nph_ndst(DataA,{-20:1:20,-20:1:20},[1,1],[0.25,0.25],'full');
  
  DataB = Airs.Tp(:,:,closest(Airs.ret_z,Settings.Z(iLevel+1)));
  STB = nph_ndst(DataB,{-20:1:20,-20:1:20},[1,1],[0.25,0.25],'full');

  
  %find cospectrum
  CoSpec = STA.ST .* conj(STB.ST);
  
  %find phase difference and amplitude
  dPhi = atan2(imag(CoSpec),real(CoSpec));
  Amp = abs(CoSpec);
  
  dPhi = permute(dPhi,[3,4,1,2]);
  Amp  = permute( Amp,[3,4,1,2]);
  
  sz = size(dPhi);
  dPhi = reshape(dPhi,sz(1),sz(2),sz(3)*sz(4));
  Amp  = reshape( Amp,sz(1),sz(2),sz(3)*sz(4));
  
  
  %find phase difference at max
  PhiO = NaN(sz(1),sz(2));
  for iX=1:1:sz(1)
    for iY=1:1:sz(2);
      [~,idx] = max(Amp(iX,iY,:));
      PhiO(iX,iY) = dPhi(iX,iY,idx);
    end
  end
  
  
    
  %abs
  PhiO = abs(PhiO);  
  
  %unwrap
  PhiO = unwrap_phase(PhiO);
  
  %convert so that pi = 1
  PhiO = PhiO./pi;

  
  
  %colour contours
  contourf(PhiO(:,Settings.PlotRange),-5:.125:5,'edgecolor','none');
%   colormap(flipud(cbrewer('div','PRGn',60)))
  colormap(cbrewer('seq','Blues',60))
  hold on
  
  %line contours
  [c,h] = contour(PhiO(:,Settings.PlotRange),-2:.25:2,'edgecolor','k');  
  clabel(c,h)
  
  
  %tidy up and export
  set(gca,'xtick',0:30:1000,'ytick',0:30:1000,'xticklabel',{},'yticklabel',{})

  caxis([0,1]*1.5)
  axis equal
  grid on
  ylim([0 90])  
  drawnow

  export_fig(['dPhi_',num2str(Settings.Z(iLevel)),'km'],'-png','-m2','-a4')
  
  
  

  
end
