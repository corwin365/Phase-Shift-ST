clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%new 2D+1 ST version, based on meeting of 07/08/2020
%Corwin Wright, c.wright@bath.ac.uk, 14/08/2020
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
% % Settings.GiD = 56;

% % Settings.Day = datenum(2010,10,16);
% % Settings.GiD = 186;

Settings.Day = datenum(2007,1,13);
Settings.GiD = [122];

% % Settings.Day = datenum(2005,7,8);
% % Settings.GiD = 85;


%analysis
%%%%%%%%%%%%%%%%%%%%%%

Settings.BasisLevel  = 40;                       %km - nearest level taken
Settings.c1          = [1,1].*0.5;               %c for  first pass (coarse in space, fine in freq)
Settings.c2          = [1,1].*0.25;              %c for second pass (coarse in freq, fine in space)
Settings.NPeaks      = 3;                        %number of spectral peaks to identify for each spatial point
Settings.Threshold   = 0;                        %threshold in spectral space to be identified as a peak. This is an amplitude for a lone 2DST (as opposed to a cospectrum)
Settings.Filt        = fspecial('gaussian',5,1); %characteristic size of point in spectral space. This is a little fatter than the default, to avoid very close peaks.
Settings.HeightRange = [20,60];                  %height range to analyse over
Settings.Thin        = 1;                        %thin out the number of scales (large runtime reduction, but changes the results)
Settings.Steps       = [0,1,2,3,4,5];           %number of steps to take phase difference over. '0' takes it from a basis level, defined above, while nonzero values use the phase shift with that many levels *above*

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Airs  = prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]);

Airs  = cat_struct(prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]),   ...
                   prep_airs_3d(Settings.Day,Settings.GiD+1,'PreSmooth',[3,3,1]), ...
                   2,{'MetaData','Source','ret_z'});

%trim useless heights
InZRange = find(Airs.ret_z >= min(Settings.HeightRange) & Airs.ret_z <= max(Settings.HeightRange));
Airs.ret_temp = Airs.ret_temp(:,:,InZRange);
Airs.Tp       = Airs.Tp(      :,:,InZRange);
Airs.BG       = Airs.BG(      :,:,InZRange);
Airs.ret_z    = Airs.ret_z(       InZRange);
clear InZRange

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2DST every level, and store as a sum and a weighted-sum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define scales
Scales{1} = 1:1:size(Airs.l1_lat,1)/2;
Scales{2} = 1:1:size(Airs.l1_lat,2)/2; Scales{2} = [reverse(-Scales{2}),Scales{2}];

%thin scales?
if Settings.Thin
  s1 = Scales{1}; Scales{1} = s1([1:1:5,6:2:end]);
  s2 = Scales{2}; Scales{2} = s2([1:1:5,6:2:end]);
  clear s1 s2
end

textprogressbar('Initial ST pass ')
for iLevel=1:1:numel(Airs.ret_z);
 
  %st level
  ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                         'FullST',true,           ...
                         'c',Settings.c1,         ...
                         'Scales',Scales);
  
  %compute weight
  CFac = exp(-(Airs.ret_z(iLevel)-42) ./ (2*scale_height(42)));  

  %store fields
  if iLevel == 1;
    %first level: just store
    Store.ST_Unweighted = ST.ST;
    Store.ST_Weighted   = ST.ST.*CFac;
    freqs = ST.freqs; %needed later
  else
    %subsequent level: integrate into existing average
    Store.ST_Unweighted = ((iLevel-1).*Store.ST_Unweighted + ST.ST        )./iLevel;
    Store.ST_Weighted   = ((iLevel-1).*Store.ST_Weighted   + ST.ST .* CFac)./iLevel;
  end


  %and loop
  textprogressbar(iLevel./numel(Airs.ret_z).*100)
end
textprogressbar('!')
clear iLevel CFac ST

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find spectral maxima in the two fields
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%find spectral peaks
PeakStore = NaN([size(Airs.l1_lat),Settings.NPeaks,3,2]);

textprogressbar('Finding peaks ')
for iX=1:1:size(Airs.l1_lat,1);
  for iY=1:1:size(Airs.l1_lat,2);
    
    %pull out the frame
    Working.U = abs(Store.ST_Unweighted(:,:,iX,iY));
    Working.W = abs(Store.ST_Weighted(  :,:,iX,iY));
    
    for iMode=1:1:2;
      switch iMode
        case 1; Frame = Working.U;
        case 2; Frame = Working.W;
      end
          
      %find maxima
      cent = FastPeakFind(Frame,0,Settings.Filt);
      if numel(cent) == 0; 
        %if we found no distinct peaks, just take the maximum as the largest "peak"
        [~,idx] = max(Frame(:));
        [x,y] = ind2sub(size(Frame),idx);
        cent(2) = x; cent(1) = y;
        clear idx x y
      end
      %convert from list to [x1,y1;x2,y2;...];
      Cent(:,1) = cent(1:2:end); %position on k_(along-track)  axis
      Cent(:,2) = cent(2:2:end); %position on k_(across-track) axis
    
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
      PeakStore(iX,iY,1:size(Cent,1),:,iMode) = Cent;
      clear Cent
      
    end
  end
  textprogressbar(iX./size(Airs.l1_lat,1).*100)
end
clear iX iY iZ Frame iMode Working cent
textprogressbar('!')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% we now have the N strongest spectral signals at each geographic point,
% found using a coarse c to provide reasonable spectral fit.
%now put back and do a finer fit to these frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find all the unique scales of importance in the data
clear Scales
s1 = unique(PeakStore(:,:,:,1,:)); s1(isnan(s1)) = []; Scales{1} = s1; clear s1;
s2 = unique(PeakStore(:,:,:,2,:)); s2(isnan(s2)) = []; Scales{2} = s2; clear s2;

%frequencies from original ST, to fix horizontal wavelengths
f1 = freqs{1};
f2 = freqs{2};

%create storage arrays for the fine STs
Store.STs  = complex(NaN([numel(Scales{1}),numel(Scales{2}), ...
                          size(Airs.l1_lat),size(Airs.ret_z,1)], ...
                          'single'));

textprogressbar('Computing fine STs ')
for iLevel=1:1:size(Airs.ret_z,1);
  
  %do ST
  ST = gwanalyse_airs_2d(Airs,Airs.ret_z(iLevel), ...
                         'FullST',true,          ...
                         'c',Settings.c2(1:2),   ...
                         'Scales',Scales);
  %store ST
  Store.STs(:,:,:,:,iLevel) = ST.ST;
  
  textprogressbar(iLevel ./ size(Airs.ret_z,1) .* 100);
  
end
clear ST iLevel
textprogressbar('!')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the requested cospectra, and store the resulting GW properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%results arrays
Store.A  = NaN([size(Airs.l1_lat),size(Airs.ret_z,1),2,numel(Settings.Steps),Settings.NPeaks]);
Store.L  = Store.A;
Store.F1 = Store.A;
Store.F2 = Store.A;

%find basis level (saves doing it every pass if used)
zidx = closest(Settings.BasisLevel,Airs.ret_z);


textprogressbar('Computing cospectra ')
for iLevel=1:1:size(Airs.ret_z,1);
  for iStep = 1:1:numel(Settings.Steps)
   
    %identify the level-pair to compute a cospectrum over
    ST1 = Store.STs(:,:,:,:,iLevel); %this level is one of them...
    if Settings.Steps(iStep) == 0;
      %The basis level
      if iLevel == zidx; continue; end
      ST2 = Store.STs(:,:,:,:,zidx);
      dZ =  Airs.ret_z(iLevel) - Airs.ret_z(zidx);
    else
      %choose the right level.
      if iLevel+Settings.Steps(iStep) <= size(Airs.ret_z,1);
        ST2 = Store.STs(:,:,:,:,iLevel+Settings.Steps(iStep));
        dZ = Airs.ret_z(iLevel) - Airs.ret_z(iLevel+Settings.Steps(iStep));
      else; continue; end
    end
    
    %ok, we have the levels. Take the cospectrum.
    CoSpectrum = ST1 .* conj(ST2);

    
    
    %then pull out the indices we need
    for iX=1:1:size(Airs.l1_lat,1);
      for iY=1:1:size(Airs.l1_lat,2);
        for iZ=1:1:Settings.NPeaks;
          for iMode=1:1:2;
            if ~isnan(PeakStore(iX,iY,iZ,3,iMode));
              
              %find new indices
              idx1 = closest(PeakStore(iX,iY,iZ,1,iMode),Scales{1});
              idx2 = closest(PeakStore(iX,iY,iZ,2,iMode),Scales{2});
              
              %cospectral value
              CS = CoSpectrum(idx1,idx2,iX,iY);
              
              %amplitude
              Store.A(iX,iY,iLevel,iMode,iStep,iZ) = sqrt(abs(CS));
              
              %vertical wavelength
              dPhi = angle(CS)./(2*pi);
              
              Store.L(iX,iY,iLevel,iMode,iStep,iZ) = dZ./dPhi;
              
              %horizontal wavelengths
              Store.F1(iX,iY,iLevel,iMode,iStep,iZ) = f1(PeakStore(iX,iY,iZ,2,iMode));
              Store.F2(iX,iY,iLevel,iMode,iStep,iZ) = f2(PeakStore(iX,iY,iZ,1,iMode));
              
            end
          end
        end
      end
    end
  end
  textprogressbar(iLevel ./ size(Airs.ret_z,1) .* 100);  
end; 
clear iLevel iX iY iZ iMode iStep CoSpectrum CS dPhi
clear idx1 idx2 dZ ST1 ST2 f1 f2 freqs PeakStore
Store = rmfield(Store,{'STs','ST_Unweighted','ST_Weighted'});
textprogressbar('!')
  

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %% phase unwrapping
% % % %only implemented for base-level relative, as:
% % % %(i) it's hard to do for the others
% % % %(ii) they're over short ranges so should be fine anyway
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % 
% % % for iStep=1:1:numel(Settings.Steps)  
% % %   if Settings.Steps(iStep) == 0;
% % %     
% % %     %pull out L
% % %     L = Store.L(:,:,:,:,iStep,:);
% % % 
% % %     %convert back to phases
% % %     sz = size(L);
% % %     dZ = permute(Airs.ret_z - Airs.ret_z(zidx),[2,3,1,4,5,6]);
% % %     dZ = repmat(dZ,sz(1),sz(2),1,sz(4),sz(5),sz(6));
% % %     L = L./dZ./(2*pi);
% % % 
% % %     %set the basis level to all-zeros
% % %     L(:,:,zidx,:,:,:) = 0;
% % %     
% % %     %unwrap
% % %     L = unwrap(L,[],3);
% % %     
% % %     %shift so it's zero at the basis level
% % %     Centre = L(:,:,zidx,:,:,:);
% % %     for iLevel=1:1:numel(Airs.ret_z);
% % %       L(:,:,iLevel,:,:,:) = L(:,:,iLevel,:,:,:) - Centre;
% % %     end
% % %     
% % %     %and convert back to wavelength
% % %     L = dZ .* L/(2*pi);
% % % % %     for iLevel=1:1:numel(Airs.ret_z);
% % % % %        dZ = Airs.ret_z(iLevel) - Settings.BasisLevel;
% % % % %        L(:,:,iLevel,:,:,:) = dZ./L(:,:,iLevel,:,:,:);
% % % % %     end
% % %     
% % %     %done! return to the pile
% % %     Store.L(:,:,:,:,iStep,:) = L;   
% % %   end
% % % end
% % % clear L iLevel Centre dZ iStep

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3dst?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ST3D = gwanalyse_airs_3d(Airs,'TwoDPlusOne',true);

%make dupe copy (w/uw) for comparison
Store.A2 =    repmat(permute(ST3D.A_2dp1,[1,2,3,6,4,5]),[1,1,1,2,1,1]);
Store.L2 = 1./repmat(permute(ST3D.m_2dp1,[1,2,3,6,4,5]),[1,1,1,2,1,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% done. save.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[yy,~,~] = datevec(Settings.Day);
dd = date2doy(Settings.Day);
OutFile = ['stout_',sprintf('%04d',yy),'d',sprintf('%03d',dd),'g',sprintf('%03d',Settings.GiD),'_thin',num2str(Settings.Thin),'.mat']


save(OutFile,'Settings','Airs','ST3D','Store','-v7.3')