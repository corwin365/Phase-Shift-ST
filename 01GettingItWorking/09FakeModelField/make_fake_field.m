clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make a big file of fake model data, which we will sample with AIRS
%weighting functions to make big waves
%
% Corwin Wright, c.wright@bath.ac.uk, 2020/12/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%output file
Settings.OutFile = 'fake_waves_2dp1.mat';

%background field to use: {Date,Granules,XT row, AT row}
%this is just to get an approximate area to simulate over
Settings.Case = {datenum(2008,1,127),21}; 

%'model' grid
Settings.Model.LonStep = 0.01; %deg
Settings.Model.LatStep = 0.01; %deg
Settings.Model.ZStep   = 0.5; %km 
Settings.Model.LonPad  = 5;
Settings.Model.LatPad  = 5;
Settings.Model.ZRange  = [0,100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% wave parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



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

%function to create gaussian packet envelope
gauss = @(x,mu,sig,amp,vo)amp*exp(-(((x-mu).^2)/(2*sig.^2)))+vo;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load the granule we'll be working on, and hence create the 'model' grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Airs,Spacing] = prep_airs_3d(Settings.Case{1},Settings.Case{2},'PreSmooth',[3,3,1]);
Spacing(3) = 3; %3km grid in vertical in region of interest

LonRange = [min(Airs.l1_lon(:)),max(Airs.l1_lon(:))] + [-1,1].*Settings.Model.LonPad;
LatRange = [min(Airs.l1_lat(:)),max(Airs.l1_lat(:))] + [-1,1].*Settings.Model.LatPad;

Model.Lon = min(LonRange):Settings.Model.LonStep:max(LonRange);
Model.Lat = min(LatRange):Settings.Model.LatStep:max(LatRange);
Model.Z   = min(Settings.Model.ZRange):Settings.Model.ZStep:max(Settings.Model.ZRange);


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
    [x,y,z] = ndgrid(1:1:numel(Model.Lon),1:1:numel(Model.Lat),1:1:numel(Model.Z));
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

    stop
    
    %store fields
    AllResults.(Vars{iVar}) = Results;
    save(Settings.OutFile,'Settings','AllResults')
    clear Results Params
  
  end; clear iRun
    
end; clear iVar 