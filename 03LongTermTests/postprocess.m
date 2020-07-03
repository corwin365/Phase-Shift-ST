clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%long-term test of 2D+1 ST
%
%Corwin Wright, c.wright@bath.ac.uk
%2020/07/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DataDir     = [LocalDataDir,'/corwin/airs_2dp1/'];
Settings.TimeScale   = datenum(2008,8,1):1:datenum(2008,8,30);
Settings.LatScale    = -90:5:90;
Settings.LonScale    = -180:5:180;
Settings.HeightScale = 20:3:60;
Settings.Vars        = {'A','k','l','m'};
Settings.Modes       ={'3D','2D'};
Settings.OutFile     = 'out_august2008.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%prepare grids
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Results.Zonal = NaN(numel(Settings.Vars),       ...
                    numel(Settings.TimeScale),  ...
                    numel(Settings.LatScale),   ...
                    numel(Settings.HeightScale), ...
                    2, ... %this 2 is for 3D 2+1
                    2);    %this 2 is for mean and median
                       
Results.Map40km = NaN(numel(Settings.Vars),       ...
                      numel(Settings.TimeScale), ...
                      numel(Settings.LatScale),  ...
                      numel(Settings.LonScale), ...
                      2, ... %this 2 is for 3D 2+1
                      2);    %this 2 is for mean and median

Results.Map30km = Results.Map40km; 
Results.Map50km = Results.Map40km;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load daily data, and grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.TimeScale)
  disp(datestr(Settings.TimeScale(iDay)))
  
  %identify and load file
  DayFile = [Settings.DataDir,'/gws_',num2str(Settings.TimeScale(iDay)),'.mat'];
  if ~exist(DayFile,'file'); clear DayFile; continue; end
  Data = load(DayFile); clear DayFile
  
  %pull out the geolocation
  Lon = double(Data.Results.Lon);
  Lat = double(Data.Results.Lat);
  Z   = double(Data.Results.Z);
  
  %and grid up 
  Lon = repmat(Lon,1,1,1,numel(Z));
  Lat = repmat(Lat,1,1,1,numel(Z));
  sz = size(Lon);
  Z   = permute(repmat(  Z,1,sz(1),sz(2),sz(3)),[2,3,4,1]);
  x = Lon(:); y = Lat(:); z = Z(:);
  clear sz Lon Lat Z
  
  %now loop over the vars, and grid them
  for iVar=1:1:numel(Settings.Vars);
    V = Data.Results.(Settings.Vars{iVar});
    for iMode=1:1:2;
      Vm = flatten(double(squeeze(V(:,:,:,:,iMode))));
     
      %zonal mean and median
      [xi,yi] = meshgrid(Settings.LatScale,Settings.HeightScale);
      Results.Zonal(iVar,iDay,:,:,iMode,1) = bin2mat(x,y,Vm,xi,yi,'@nanmean')';
      Results.Zonal(iVar,iDay,:,:,iMode,2) = bin2mat(x,y,Vm,xi,yi,'@nanmedian')';
      
      %30km map
      [xi,yi] = meshgrid(Settings.LonScale,Settings.LatScale);
      zval = 30; zidx= find(z == zval); clear zval
      Results.Map30km(iVar,iDay,:,:,iMode,1) =  bin2mat(x(zidx),y(zidx),Vm(zidx),xi,yi,'@nanmean');
      Results.Map30km(iVar,iDay,:,:,iMode,2) =  bin2mat(x(zidx),y(zidx),Vm(zidx),xi,yi,'@nanmedian');

      %40km map
      zval = 39; zidx= find(z == zval); clear zval
      Results.Map40km(iVar,iDay,:,:,iMode,1) =  bin2mat(x(zidx),y(zidx),Vm(zidx),xi,yi,'@nanmean');
      Results.Map40km(iVar,iDay,:,:,iMode,2) =  bin2mat(x(zidx),y(zidx),Vm(zidx),xi,yi,'@nanmedian');      

      %50km map
      zval = 51; zidx= find(z == zval); clear zval
      Results.Map50km(iVar,iDay,:,:,iMode,1) =  bin2mat(x(zidx),y(zidx),Vm(zidx),xi,yi,'@nanmean');
      Results.Map50km(iVar,iDay,:,:,iMode,2) =  bin2mat(x(zidx),y(zidx),Vm(zidx),xi,yi,'@nanmedian');         
      
      clear zidx xi yi
    end
    
  end
  
  
  
end

save(Settings.OutFile,'Settings','Results')