% clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%long-term test of 2D+1 ST
%
%Corwin Wright, c.wright@bath.ac.uk
%2020/07/02
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Settings.DayScale  = datenum(2008,1,127);%datenum(2010,10,8):1:datenum(2010,10,20);
Settings.OutDir = [LocalDataDir,'/corwin/airs_2dp1/'];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iDay=1:1:numel(Settings.DayScale)
  disp(datestr(Settings.DayScale(iDay)))
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% create results arrays
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  Results.A        = single(NaN(240,14,90,135,2));
  Results.k        = Results.A;
  Results.l        = Results.A;
  Results.m        = Results.A;
  Results.BG       = Results.A;

  Results.Lon      = single(NaN(240,90,135));
  Results.Lat      = Results.Lon;

  Middle = 68:68+135-1;
  
  
  OutFile = [Settings.OutDir,'/gws_',num2str(Settings.DayScale(iDay)),'.mat'];
  if exist(OutFile); disp('Already done'); continue; end
  
  for iGranule=1:1:240
    if mod(iGranule,60) == 0; disp(iGranule); end
    try
      
      %get granule metadata for this granule and the one before and after
      [Airs1,~,      Error,ErrorInfo] = prep_airs_3d(Settings.DayScale(iDay),iGranule,  'DayNightFlag',true,'Extrapolant','linear');
      if Error ~=0; continue; disp(ErrorInfo); end
      [Airs2,Spacing,Error,ErrorInfo] = prep_airs_3d(Settings.DayScale(iDay),iGranule+1,'DayNightFlag',true,'Extrapolant','linear');
      if Error ~=0; continue; disp(ErrorInfo); end
      clear Error
      
      
      %cat the granules together
      Vars = {'l1_lat','l1_lon','l1_time','ret_temp','Tp','BG','DayNightFlag'};
      Airs = Airs1;
      for iVar=1:1:numel(Vars)
        Airs.(Vars{iVar}) = cat(2,Airs.(Vars{iVar}),Airs2.(Vars{iVar}));
      end
      clear Vars Airs1 Airs2 Airs3 iVar
      
      %skip if we have no nighttime data
      if numel(find(Airs.DayNightFlag == 0)) == 0; continue; end
     
      %s-transform
      [ST,Airs,Error,ErrorInfo] = gwanalyse_airs_3d(Airs,'TwoDPlusOne',true);%,varargin)

      %discard daytime data for variables we will be storing
      Airs.l1_lon(Airs.DayNightFlag == 1) = NaN;
      Airs.l1_lat(Airs.DayNightFlag == 1) = NaN;  
      DNF = repmat(Airs.DayNightFlag,1,1,size(Airs.Tp,3));
      ST.A(DNF == 1) = NaN;
      ST.k(DNF == 1) = NaN;      
      ST.l(DNF == 1) = NaN;      
      ST.m(DNF == 1) = NaN;      
      ST.A(DNF == 1) = NaN;      
      Airs.BG(DNF == 1) = NaN;
      
      %store variables of interest
      Results.Lon(     iGranule,:,:)   = Airs.l1_lon(    :,Middle);
      Results.Lat(     iGranule,:,:)   = Airs.l1_lat(    :,Middle);
      Results.Z = Airs.ret_z;      
      
      Results.A(       iGranule,:,:,:,1) = permute(ST.A(   :,Middle,:),[3,1,2]);
      Results.k(       iGranule,:,:,:,1) = permute(ST.k(   :,Middle,:),[3,1,2]);
      Results.l(       iGranule,:,:,:,1) = permute(ST.l(   :,Middle,:),[3,1,2]);
      Results.m(       iGranule,:,:,:,1) = permute(ST.m(   :,Middle,:),[3,1,2]);
      Results.BG(      iGranule,:,:,:,1) = permute(Airs.BG(:,Middle,:),[3,1,2]);

      Results.A(       iGranule,:,:,:,2) = permute(ST.A_2dp1(   :,Middle,:),[3,1,2]);
      Results.k(       iGranule,:,:,:,2) = permute(ST.k_2dp1(   :,Middle,:),[3,1,2]);
      Results.l(       iGranule,:,:,:,2) = permute(ST.l_2dp1(   :,Middle,:),[3,1,2]);
      Results.m(       iGranule,:,:,:,2) = permute(ST.m_2dp1(   :,Middle,:),[3,1,2]);

    catch; end
  end
  
  
  save(OutFile,'Settings','Results');
end

