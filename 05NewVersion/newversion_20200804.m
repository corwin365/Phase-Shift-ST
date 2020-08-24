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

Settings.Day = datenum(2008,1,127);
Settings.GiD = 56;

% % Settings.Day = datenum(2010,10,16);
% % Settings.GiD = 186;

% % Settings.Day = datenum(2007,1,13);
% % Settings.GiD = [122];

% % Settings.Day = datenum(2005,7,8);
% % Settings.GiD = 85;


%analysis
%%%%%%%%%%%%%%%%%%%%%%

Settings.BasisLevel = 40;                       %km - nearest level taken
Settings.c1         = [1,1].*0.5;               %c for  first pass (coarse in space, fine in freq)
Settings.c2         = [1,1].*0.25;              %c for second pass (coarse in freq, fine in space)
Settings.NPeaks     = 3;                        %number of spectral peaks to identify for each spatial point
Settings.Threshold  = 0;                        %threshold in spectral space to be identified as a peak. this is an empirical value from testing and I am unsure of units yet, needs more work!
Settings.Filt       = fspecial('gaussian',9,3); %characteristic size of point in spectral space. This is a little fatter than the default, to avoid very close peaks.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Airs  = prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]);

Airs  = cat_struct(prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]),   ...
                   prep_airs_3d(Settings.Day,Settings.GiD+1,'PreSmooth',[3,3,1]), ...
                   2,{'MetaData','Source','ret_z'});

%trim useless heights
InZRange = find(Airs.ret_z >= 20 & Airs.ret_z <= 60);
Airs.ret_temp = Airs.ret_temp(:,:,InZRange);
Airs.Tp       = Airs.Tp(      :,:,InZRange);
Airs.BG       = Airs.BG(      :,:,InZRange);
Airs.ret_z    = Airs.ret_z(       InZRange);
clear InZRange


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
clear Scales                         

%find spectral peaks
PeakStore = NaN([size(Airs.l1_lat),Settings.NPeaks,3]);

for iX=1:1:size(Airs.l1_lat,1);
  for iY=1:1:size(Airs.l1_lat,2);
    
    %pull out the frame
    Frame = abs(Basis.ST.ST(:,:,iX,iY));
    
    %find maxima   
    cent = FastPeakFind(Frame,0,Settings.Filt);
    if numel(cent) == 0; continue; end    
    %convert from list to [x1,y1;x2,y2;...];
    Cent(:,1) = cent(1:2:end); %position on k_(along-track)  axis
    Cent(:,2) = cent(2:2:end); %position on k_(across-track) axis
    clear cent
    
    %find values of these indices
    for iZ=1:1:size(Cent,1); Cent(iZ,3) = Frame(Cent(iZ,2),Cent(iZ,1)); end
    
    %remove small values
    Good = find(Cent(:,3) > Settings.Threshold);
    Cent = Cent(Good,:);
    clear Good
       
    %sort the data, select the number of peaks desired, and store
    [~,idx] = sort(Cent(:,3),'desc');
    if size(Cent,1) > Settings.NPeaks; Cent = Cent(idx(1:Settings.NPeaks),:);
    else;                              Cent = Cent(idx,:);
    end
    clear idx
    PeakStore(iX,iY,1:size(Cent,1),:) = Cent;
    clear Cent
    
    
% % % %     %convert to wavelengths (for sanity testing, not used in "production")
% % % %     f1 = Basis.ST.freqs{1};   f2 = Basis.ST.freqs{2}; 
% % % %     for iZ=1:1:Settings.NPeaks
% % % %       if ~isnan(PeakStore(iX,iY,iZ,3));
% % % %         PeakStore(iX,iY,iZ,1) = 1./f2(PeakStore(iX,iY,iZ,1));
% % % %         PeakStore(iX,iY,iZ,2) = 1./f1(PeakStore(iX,iY,iZ,2));
% % % %       end
% % % %     end

  end
end
clear iX iY iZ Frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% we now have the N strongest spectral signals at each geographic point,
% found using a coarse c to provide reasonable spectral fit.
%now we need to find how these extend to other levels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find all the unique scales of importance in the data
clear Scales
s1 = unique(PeakStore(:,:,:,1)); s1(isnan(s1)) = []; Scales{1} = s1; clear s1;
s2 = unique(PeakStore(:,:,:,2)); s2(isnan(s2)) = []; Scales{2} = s2; clear s2;

%find these scales at every level, using the tighter c values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%first, redo the basis level
Basis.ST2 = gwanalyse_airs_2d(Airs,Airs.ret_z(zidx), ...
                             'FullST',true,          ...
                             'c',Settings.c2(1:2),   ... %note this is c2, NOT c1
                             'Scales',Scales);

%now compute each other level and take the cospectrum, then store what we need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%results arrays
Store.A  = NaN([size(Airs.l1_lat),Settings.NPeaks,size(Airs.ret_z,1)]);
Store.L  = Store.A;
Store.F1 = Store.A;
Store.F2 = Store.A;

%frequencies from original ST, to fix horizontal wavelengths
f1 = Basis.ST.freqs{1};
f2 = Basis.ST.freqs{2};

for iLevel=1:1:size(Airs.ret_z,1);
  if iLevel == zidx; continue; end %the same by definition
  
  %do ST
  ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                         'FullST',true,          ...
                         'c',Settings.c2(1:2),   ...
                         'Scales',Scales);
 
  %take cospectrum
  CoSpectrum = ST.ST .* conj(Basis.ST2.ST);
  
  %pull out the right spectral points geographically, then take amplitude and phase difference
  for iX=1:1:size(Airs.l1_lat,1);
    for iY=1:1:size(Airs.l1_lat,2);
      for iZ=1:1:Settings.NPeaks;
        if ~isnan(PeakStore(iX,iY,iZ,3));
          
          %find new indices
          idx1 = closest(PeakStore(iX,iY,iZ,1),Scales{1});
          idx2 = closest(PeakStore(iX,iY,iZ,2),Scales{2});
          
          %cospectral value
          CS = CoSpectrum(idx1,idx2,iX,iY);
          
          %amplitude
          Store.A(iX,iY,iZ,iLevel) = sqrt(abs(CS));
          
          %vertical wavelength
          dZ = Airs.ret_z(iLevel) - Airs.ret_z(zidx);
          Store.L(iX,iY,iZ,iLevel) = 2.* pi .* abs(dZ ./angle(CS));
          
          %horizontal wavelengths
          Store.F1(iX,iY,iZ,iLevel) = f1(PeakStore(iX,iY,iZ,2));
          Store.F2(iX,iY,iZ,iLevel) = f2(PeakStore(iX,iY,iZ,1));
          
          
        end
      end
    end
  end
end
clear iLevel ST CoSpectrum iX iY iZ idx1 idx2 CS f1 f2 dZ zidx Scales
             

