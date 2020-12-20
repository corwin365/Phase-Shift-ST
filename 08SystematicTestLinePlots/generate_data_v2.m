clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyse systematic ideal wave variations with 3DST and 2D+1
%
%modified to apply vertical weighting functions
%NOT USED FOR PAPER
%
%
% Corwin Wright, c.wright@bath.ac.uk, 2020/10/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output file
Settings.OutFile = 'out_testing_vert.mat';


%background field to use: {Date,Granules,XT row, AT row}
%will be loaded regardless, but set to zero for flat-field analyses
%this should not have a nice wave in for these tests, just represent a
%reasonable case of background and noise level
Settings.Case = {datenum(2008,1,127),1}; 

%portion of the granule to average the result over
Settings.AvZone.X = (-15:1:15)+45;
Settings.AvZone.Y = (-15:1:15)+67;
Settings.AvZone.Z = ( -3:1: 3)+14;
Settings.Averages = {'nanmean','nanmedian','mode'};

%variables defining the wave we will create
%only one will be varied at a time, and the first-listed value will be
%the one used while other variables are varied. Duplicates will be removed
%so don't worry about putting the default in twice.
Settings.Wave.Amplitude = 10;%[5,0.25:0.25:10];  %unique([5,1:1:5],'stable'); %K
Settings.Wave.Lambdax   = [100:100:1000];  %km - across track
Settings.Wave.Lambday   = 10000;%[1000,50:50:1500]; %km - along track
Settings.Wave.Lambdaz   = 60;%[10:10:80];%[25,2.5:2.5:100];  %km
Settings.Wave.Rotationx = 0;%[0:10:170];        %degrees in x to rotate wave
Settings.Wave.Rotationy = 0;%[0:10:170];        %degrees in y to rotate wave
Settings.Wave.Rotationz = 0;%[0:10:170];        %degrees in z to rotate wave

Settings.Wave.PacketWidthx =    500;  %km width of packet in xt dir
Settings.Wave.PacketWidthy =   1000; %km width of packet in at dir
Settings.Wave.PacketWidthz =     40; %km width of packet in z dir

%variable to retain for later analysis
Settings.OutVars = {'A','m','kh','th'};

%function to create gaussian packet envelope
gauss = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the granule we'll be working on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Airs,Spacing] = prep_airs_3d(Settings.Case{1},Settings.Case{2},'PreSmooth',[3,3,1]);
Spacing(3) = 3; %3km grid in vertical in region of interest


%instrument resolution as f(z)
%define the vertical resolution of the retrieval at each height
%taken from figure 5a of Hoffmann and Alexander 2009, cz=50km curve
%first column is altitude in km, second colum is resolution in km
%added a 30 to the top and bottom to represent the very large values off the ends
Res=[15,30;18,10.25;21,6.5;24,6.75;27,7.25;30,7.5;33,7.75;36,8.75;39,9;42,9;45,10;48,10;51,11;54,13.5;57,14.5;60,13.25;65,13.5;70,19.5;75,30];
Res(:,2) = Res(:,2) ./ 2.355; %fwhm -> stdev




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify which properties we want to vary over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find how many tests we need for each one
VarList = fieldnames(Settings.Wave);
NVaried = ones(size(VarList));
for iVar=1:1:numel(NVaried); 
  NVaried(iVar) = numel(Settings.Wave.(VarList{iVar}));
  if NVaried(iVar) > 1; Settings.Wave.(VarList{iVar}) = unique(Settings.Wave.(VarList{iVar}),'stable'); end
end

%sort so the ones with the fewest combinations go first
[~,idx] = sort(NVaried,'asc'); 
NVaried = NVaried(idx);
VarList = VarList(idx);

%and only process the ones with more than one setting
Vars = VarList(NVaried > 1); 

%tidy
clear iVar NVaried idx

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vary, compute and plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%loop over the variables we want to vary over
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar=1:1:numel(Vars)
  
  %pull out wave properties, both fixed and varing, from settings   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  Params = struct();
  for jVar=1:1:numel(VarList)
    if ~strcmp(Vars{iVar},VarList{jVar})
      v = Settings.Wave.(VarList{jVar});
      Params.(VarList{jVar}) = v(1);
    else; ToVary = Settings.Wave.(VarList{jVar});
    end
  end; clear jVar v
  
  %create results array   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  Results = NaN(2,2,                      ...  %first 2 is flat vs realistic, second 2 is 3D vs 2D+1
                numel(Settings.OutVars ), ...
                numel(Settings.Averages), ...
                numel(ToVary));
  
  %loop over this single variable, keeping all others fixed
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
  for iRun=1:1:numel(ToVary)
    
    disp(['==========================================='])    
    disp(['==========================================='])
    disp(['Varying ',Vars{iVar},', run ',num2str(iRun),' of ',num2str(numel(ToVary))])
    disp(['==========================================='])    
    disp(['==========================================='])    
    
    
    %copy over this run's param to the master list
    Params.(Vars{iVar}) = ToVary(iRun);
    
    %make the fake wave
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

    %apply 'weighting functions' to wave
    sz = size(Wave);

    Wave = reshape(Wave,sz(1)*sz(2),sz(3));
    Wave2 = Wave;
    for iLevel=1:1:size(Wave,2);
      z = Airs.ret_z(iLevel);
      r = interp1(Res(:,1),Res(:,2),[z],'linear','extrap');
      for iX=1:1:sz(1)*sz(2);
        w = filter(gausswin(8,r./3)./sum(gausswin(8,r./3)),1,squeeze(Wave(iX,:)),[],4);
        Wave2(iX,iLevel) = w(iLevel);
      end
    end
    Wave = reshape(Wave2,sz(1),sz(2),sz(3));

    clear Wave2 sz iLevel z r iX
    
    
    %do ST analysis, with and without flat-fielding
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    Tp = Airs.Tp; %temporarily  backup, so we can go back to it
    for Flat=[1,2];

      %generate field
      if Flat == 1; Airs.Tp = Airs.Tp + Wave; disp('===> Real noise');
      else          Airs.Tp = Wave;           disp('===> Flat-field');
      end

      
      Airs.Tp = smoothn(Airs.Tp,[3,3,1]);
        
      %do 3DST and 2D+1 ST
      ST = gwanalyse_airs_3d(Airs,'ZRange',[0 90],'TwoDPlusOne',true,'HeightScaling',false);

      %convert F1 and F2 to kx and th
      ST.kh      = quadadd(ST.F1,ST.F2); 
      ST.th      = atan2d(ST.F2,ST.F1);
      ST.kh_2dp1 = quadadd(ST.F1_2dp1,ST.F2_2dp1);
      ST.th_2dp1 = atan2d(ST.F2_2dp1,ST.F1_2dp1);

      %take average of region of interest
      for iAverage =1:1:numel(Settings.Averages)
        for iOutVar=1:1:numel(Settings.OutVars);
        fh = str2func(Settings.Averages{iAverage});
        VarA = ST.(Settings.OutVars{iOutVar});
        VarB = ST.([Settings.OutVars{iOutVar},'_2dp1']);

        Results(Flat,1,iOutVar,iAverage,iRun) = abs(fh(VarA(Settings.AvZone.X, ...
                                                            Settings.AvZone.Y, ...
                                                            Settings.AvZone.Z),'all'));
        Results(Flat,2,iOutVar,iAverage,iRun) = abs(fh(VarB(Settings.AvZone.X, ...
                                                            Settings.AvZone.Y, ...
                                                            Settings.AvZone.Z),'all'));  
        end; clear VarA VarB
      end;
      clear iAverage fh ST iOutVar
    end; 
    Airs.Tp = Tp; %put back in place, for future loops
    clear Tp Flat
        
  end; 
  
  %store results in master struct
  AllResults.(Vars{iVar}) = Results;
  save(Settings.OutFile,'Settings','AllResults')
  clear Results iRun Params
  
end; clear iVar 