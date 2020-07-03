clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%third draft of generalised 2D+1 ST analysis
%Corwin Wright, c.wright@bath.ac.uk, 20200624
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%let's go for something more final. try and match the wavelengths in a 
%sensible way between k/l and m, and then tidy up the code to do decent outputs


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%data to analyse
% Settings.Day = datenum(2008,1,127);
% Settings.GiD = [56];


Settings.Day = datenum(2010,10,16);
Settings.GiD = 186;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data as a granule-pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get granule
Airs  = cat_struct(prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]),   ...
                   prep_airs_3d(Settings.Day,Settings.GiD+1,'PreSmooth',[3,3,1]), ...
                   2,{'MetaData','Source','ret_z'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3dst the granule-pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ST3D,Airs] = gwanalyse_airs_3d(Airs,'TwoDPlusOne',true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot some maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 50;
zidx = closest(z,Airs.ret_z);
z = Airs.ret_z(zidx);


clf
set(gcf,'color','w')
redyellowblue32

subplot(2,2,1)
pcolor(Airs.l1_lon,Airs.l1_lat,squeeze(ST3D.A(:,:,zidx))); shading flat; colorbar
title(['Amplitude, 3DST, ',num2str(z),'km altitude'])


subplot(2,2,2)
pcolor(Airs.l1_lon,Airs.l1_lat,squeeze(ST3D.A_2dp1(:,:,zidx))); shading flat; colorbar
title(['Amplitude, 2D+1, ',num2str(z),'km altitude'])

subplot(2,2,3)
pcolor(Airs.l1_lon,Airs.l1_lat,1./squeeze(ST3D.m(:,:,zidx))); shading flat; colorbar
title(['\lambda_z, 3DST, ',num2str(z),'km altitude'])
caxis([10 40])


subplot(2,2,4)
pcolor(Airs.l1_lon,Airs.l1_lat,1./squeeze(ST3D.m_2dp1(:,:,zidx))); shading flat; colorbar
title(['\lambda_z, 2D+1, ',num2str(z),'km altitude'])
caxis([10 40])

sgtitle([datestr(Settings.Day),', g',num2str(Settings.GiD)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot some profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(1,3,1)
A.ST3D = ST3D.A; A.ST3D(A.ST3D < 1.5) = NaN;
plot(squeeze(nanmean(A.ST3D,[1,2])),Airs.ret_z,'color','r')
hold on
A.ST2D = ST3D.A_2dp1; A.ST2D(A.ST2D < 1.5) = NaN;
plot(squeeze(nanmean(A.ST2D,[1,2])),Airs.ret_z,'color','k')
title('Amplitude')

subplot(1,3,2)
A.ST3D = ST3D.A; m = ST3D.m; m(A.ST3D < 1.5) = NaN;
plot(1./squeeze(nanmean(m,[1,2])),Airs.ret_z,'color','r')
hold on
A.ST2D = ST3D.A_2dp1; m = ST3D.m_2dp1; m(A.ST2D < 1.5) = NaN;
plot(squeeze(1./nanmean(m,[1,2])),Airs.ret_z,'color','k')
title('Lz')
xlim([10 40])

subplot(1,3,3)
A.ST3D = ST3D.A; k = ST3D.k; k(A.ST3D < 1.5) = NaN;
plot(1./squeeze(nanmean(k,[1,2])),Airs.ret_z,'color','r')
hold on
A.ST2D = ST3D.A_2dp1; k = ST3D.k_2dp1; k(A.ST2D < 1.5) = NaN;
plot(squeeze(1./nanmean(k,[1,2])),Airs.ret_z,'color','k')
title('k^{-1}')
xlim([300 1000])
