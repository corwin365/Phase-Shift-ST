clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% draw diagram showing how the artificial waves are characterised
%
% Corwin Wright, c.wright@bath.ac.uk, 2021/2/2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%background field to use: {Date,Granules,XT row, AT row}
%will be loaded regardless, but set to zero for flat-field analyses
%this should not have a nice wave in for these tests, just represent a
%reasonable case of background and noise level
Settings.Case = {datenum(2008,1,127),1}; 

%portion of the granule to average the result over
Settings.AvZone.X = (-25:1:25)+45;
Settings.AvZone.Y = (-25:1:25)+67;
Settings.AvZone.Z = ( -3:1: 3)+14;
Settings.Averages = {'nanmean','nanmedian','mode'};

%variables defining the wave we will create
Settings.Wave.Amplitude    =   25;
Settings.Wave.Lambdax      =  500;
Settings.Wave.Lambday      = 1000;
Settings.Wave.Lambdaz      =   25;
Settings.Wave.Rotationx    =   20;       
Settings.Wave.Rotationy    =   20;   
Settings.Wave.Rotationz    =   20;  
Settings.Wave.PacketWidthx = 1200;  %km width of packet in xt dir
Settings.Wave.PacketWidthy =  900; %km width of packet in at dir
Settings.Wave.PacketWidthz =   50; %km width of packet in z dir

%function to create gaussian packet envelope
gauss = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

%variables to plot, and the colour levels to show
Settings.PlotVars   = {'Wave'};%,'A_2dp1'};
Settings.ColourLevs = {[2],[5]}; %both negative and positve values will be shown

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the granule we'll be basing this on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Airs,Spacing] = prep_airs_3d(Settings.Case{1},Settings.Case{2},'PreSmooth',[3,3,1]);
Spacing(3) = 3; %3km grid in vertical in region of interest



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute the fake wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Params = Settings.Wave;

%create a sinusoid packet going in 1d over a few levels
kx = (2*pi)./(Params.Lambdax./Spacing(1));
ky = (2*pi)./(Params.Lambday./Spacing(2));
kz = (2*pi)./(Params.Lambdaz./Spacing(3));

%make wave
[x,y,z] = ndgrid(1:1:500,1:1:500,1:1:200); %will trim to size later
wave = Params.Amplitude .* sin(kx.*x + ky.*y + kz.*z + pi/2);
clear x y z

%rotate wave
wave = imrotate(wave,Params.Rotationx,'bilinear','crop'); %axis 1
wave = permute(imrotate(permute(wave,[2,1,3]),Params.Rotationy,'bilinear','crop'),[2,1,3]); %axis 2
wave = permute(imrotate(permute(wave,[1,3,2]),Params.Rotationz,'bilinear','crop'),[1,3,2]); %axis 3

%trim to size, taking the middle bit only (this avoids bits rotated off the edges)
wave = wave((1:1:size(Airs.l1_lat,1)) - floor(mean(size(Airs.l1_lat,1))./2) + 250,...
            (1:1:size(Airs.l1_lat,2)) - floor(mean(size(Airs.l1_lat,2))./2) + 250,...
            (1:1:size(Airs.ret_z, 1)) - floor(mean(size(Airs.ret_z, 1))./2) + 100);

%packetise
Gauss.x = gauss(1:1:size(Airs.l1_lat,1),size(Airs.l1_lat,1)./2,Params.PacketWidthx./Spacing(1)./2.355,1,0);
Gauss.y = gauss(1:1:size(Airs.l1_lat,2),size(Airs.l1_lat,2)./2,Params.PacketWidthy./Spacing(2)./2.355,1,0);
Gauss.z = gauss(1:1:size(Airs.ret_z, 1),size( Airs.ret_z,1)./2,Params.PacketWidthz./Spacing(3)./2.355,1,0);
[Gx,Gy,Gz] = ndgrid(Gauss.x,Gauss.y,Gauss.z); G = Gx.*Gy.*Gz;
G = G./0.7; G(G > 1) = 1; %make a flat region at the top

Wave = wave.*G;

clear x y z kx ky kz wave Gx Gy Gz

clear G gauss Gauss Params

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ST the wave, so we can plot output examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Airs.Tp = Wave;
ST = gwanalyse_airs_3d(Airs,'ZRange',[0 90],'TwoDPlusOne',true,'Spacing',Spacing);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make location matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sz = size(Wave);
[x,y,z] = meshgrid(1:sz(2),1:sz(1),1:sz(3));

% % x = repmat(1:1:90,sz(2),1,sz(3));
% % stop
% % y = permute(repmat(1:1:135,1,sz(1),sz(3)),[2,1,3]);
% % z = permute(repmat(Airs.ret_z,1,sz(1),sz(2)),[2,3,1]);
% % clear sz


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iPlotVar= 1%1:1:numel(Settings.PlotVars)
  
  
  %prepare figure
  clf
  set(gcf,'color','w')
  
  
  
  hold on
  
  %get variable
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch Settings.PlotVars{iPlotVar}
    case 'Wave';    PlotData = Wave;
    case 'A_2dp1';  PlotData = ST.A_2dp1;
  end
  
  
  %overinterpolate the wave
  %%%%%%%%%%%%%%%%%%%%%%%%%
  
  Oversample = 3;
  [xi,yi,zi] =  meshgrid(1:1./Oversample:sz(2),1:1./Oversample:sz(1),1:1./Oversample:sz(3));
  Wave2 = interp3(x,y,z,PlotData,xi,yi,zi);
  Wave2 = smoothn(Wave2,[9,9,9]);
  
  %scale to size of Airs granule
  xi = xi.*Spacing(1);
  yi = yi.*Spacing(2);
  zi = zi.*Spacing(3);
  
  %fill edges with zeros, to avoid truncated volumes
  Wave2([1,end],:,:) = NaN;
  Wave2(:,[1,end],:) = NaN;
  Wave2(:,:,[1,end]) = NaN;  
  
  %plot isosurfaces
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  IsoSurfaces= Settings.ColourLevs{iPlotVar};
  Alphas = ones(size(IsoSurfaces)).*0.5;
  
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
      fv = isosurface(xi,yi,zi,Wave2,IsoSurfaces(iIso));
      eval([PatchId,' = patch(fv);'])
      
      % plot as transparent surfaces
      eval([PatchId,'.FaceColor = Colours(iIso,:);']); %
      eval([PatchId,'.EdgeColor = ''none'';']);
      eval([PatchId,'.FaceAlpha  = Alphas(iIso);']);
      
      %fix reflection
      eval(['set(',PatchId,',''DiffuseStrength'',1,''SpecularStrength'',0.5,''SpecularExponent'',15)']);
      
      %done
      eval(['clear ',PatchId,';'])
      
      drawnow
    end
    
  end
  
  
  
  %plot averaging cage
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Cage = Settings.AvZone;
  Colour = [1,1,1]*0;
  LineStyle = '-';
  hold on
  
  Cage.X = Cage.X.*Spacing(1);
  Cage.Y = Cage.Y.*Spacing(2);
  Cage.Z = Cage.Z.*Spacing(3);
  
  plot3(Cage.Y([1,end,end,1,1]),Cage.X([  1,  1,end,end,  1]),Cage.Z([  1,  1,  1,  1,  1]),LineStyle,'color',Colour,'linewi',3)
  plot3(Cage.Y([1,end,end,1,1]),Cage.X([end,end,end,end,end]),Cage.Z([  1,  1,end,end,  1]),LineStyle,'color',Colour,'linewi',3)
  plot3(Cage.Y([1,end,end,1,1]),Cage.X([  1,  1,  1,  1,  1]),Cage.Z([  1,  1,end,end,  1]),LineStyle,'color',Colour,'linewi',3)
  plot3(Cage.Y([1,end,end,1,1]),Cage.X([  1,  1,end,end,  1]),Cage.Z([end,end,end,end,end]),LineStyle,'color',Colour,'linewi',3)
  
  
  %tidy up
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  xlabel('x')
  ylabel('y')
  zlabel('z')
  
  % axis([-5 135 -5 90 -5 31] + [20 -20 10 -10 3 -3])
  
%   camlight
  lighting gouraud
  view([40,20])
  axis square
  grid on
  % set(gca,'Projection','perspective')
  
  stpp
end