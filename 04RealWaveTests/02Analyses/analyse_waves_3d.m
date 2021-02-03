clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do 2D+1 analysis of some real AIRS-observed waves
%
%Corwin Wright, c.wright@bath.ac.uk, 2021/02/02
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%wave we want to analyse
Settings.WaveID = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load waves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


switch Settings.WaveID
  case 1;
    %Andes big 2008
    [Airs,Spacing] = prep_airs_3d(datenum(2008,1,127),57);
  case 2;
    %Scandinavia big
    [Airs,Spacing] = prep_airs_3d(datenum(2007,1,13),122);
  case 3;
    %midwest convective
    [Airs,Spacing] = prep_airs_3d(datenum(2005,7,8),85);
    Airs2 = prep_airs_3d(datenum(2005,7,8),86);
    Airs.l1_lat = cat(2,Airs.l1_lat,Airs2.l1_lat);
    Airs.l1_lon = cat(2,Airs.l1_lon,Airs2.l1_lon);
    Airs.Tp     = cat(2,Airs.Tp,    Airs2.Tp);
    clear Airs2
  otherwise
    stop
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% analyse wave with 2D+1 and 3D ST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ST = gwanalyse_airs_3d(Airs,'ZRange',[20 60],'TwoDPlusOne',true,'Spacing',Spacing);
clear Spacing
  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  



  
  
  
  