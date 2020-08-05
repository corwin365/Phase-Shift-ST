clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%new 2D+1 ST version, based on meeting of 24/07/2020
%Corwin Wright, c.wright@bath.ac.uk, 04/08/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%granule
%%%%%%%%%%%%%%%%%%%%%%

% % Settings.Day = datenum(2008,1,127);
% % Settings.GiD = [56];

% % Settings.Day = datenum(2010,10,16);
% % Settings.GiD = 186;

Settings.Day = datenum(2007,1,13);
Settings.GiD = [122];


%analysis
%%%%%%%%%%%%%%%%%%%%%%

Settings.BasisLevel = 40;            %km - nearest level taken
Settings.c1         = [1,1,1].*0.5;  %c for first pass
Settings.c2         = [1,1,1].*0.25; %c for second pass
Settings.NPeaks     = 3;             %number of spectral peaks to identify for each spatial point
Settings.Threshold  = 0.05;          %threshold in spectral space to be identified as a peak. this is a test and i am unsure of units yet, needs more work!


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Airs  = prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]);

% % Airs  = cat_struct(prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]),   ...
% %                    prep_airs_3d(Settings.Day,Settings.GiD+1,'PreSmooth',[3,3,1]), ...
% %                    2,{'MetaData','Source','ret_z'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ST basis level using coarse c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%identify basis level
zidx = closest(Airs.ret_z,Settings.BasisLevel);

%define scales
Scales{1} = 1:1:size(Airs.l1_lat,1)/2;
Scales{2} = 1:1:size(Airs.l1_lat,2)/2; Scales{2} = [reverse(-Scales{2}),Scales{2}];

Basis.ST = gwanalyse_airs_2d(Airs,Airs.ret_z(zidx), ...
                            'FullST',true,          ...
                            'c',Settings.c1(1:2),   ...
                            'Scales',Scales);

%find spectral peaks
PeakStore = NaN([size(Airs.l1_lat),Settings.NPeaks,3]);

%temporary warning as logic below doesn't work yet...
disp('NEED TO ADAPT PEAK FINDER TO SPECTRAL SPACE MAPPING')

for iX=1:1:size(Airs.l1_lat,1);
  for iY=1:1:size(Airs.l1_lat,2);
    
    %pull out the frame
    Frame = abs(Basis.ST.ST(:,:,iX,iY));
    
% %     %sometimes the dynamic range is very low, usually if there is a 
% %     %really big wave present. Scaling should increase the d-range.
% %     Frame = 10.^(Frame);
    
    %find maxima   
    cent = FastPeakFind(Frame);
    if numel(cent) == 0; continue; end    
    %convert from list to [x1,y1;x2,y2;...];
    Cent(:,1) = cent(1:2:end); %position on k_(along-track)  axis
    Cent(:,2) = cent(2:2:end); %position on k_(across-track) axis
    clear cent

    %find values of these indices
    for iZ=1:1:size(Cent,1); Cent(iZ,3) = Frame(Cent(iZ,2),Cent(iZ,1)); end
       
    %sort the data, select the number of peaks desired, and store
    [~,idx] = sort(Cent(:,3),'desc');
    if size(Cent,1) > Settings.NPeaks; Cent = Cent(idx(1:Settings.NPeaks),:);
    else;                              Cent = Cent(idx,:);
    end
    clear idx
    PeakStore(iX,iY,1:size(Cent,1),:) = Cent;
    clear Cent
    
  end
end
                          