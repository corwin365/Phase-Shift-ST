clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot height-longitude plots of old and new Lz output
%
%Corwin Wright, c.wright@bath.ac.uk, 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%get data
load('out_may2008.mat')
DayNumber = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% vertical wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1); set(gcf,'color','w')
redyellowblue16

%get amplitude data
Lambda = squeeze(Results.Zonal(4,:,:,:,:,:));

% %take time mean
% Lambda = 1./squeeze(nanmean(Lambda,1));
%take daily value
Lambda = 1./squeeze(Lambda(DayNumber,:,:,:,:,:));

%plot 1 - 3DST, mean
subplot(2,2,1)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,1,1))');
set(gca,'ydir','normal')
colorbar
title('Lz, 3DSt, mean')

%plot 2 - 3DST, median
subplot(2,2,2)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,1,2))');
set(gca,'ydir','normal')
colorbar
title('Lz, 3DSt, median')

%plot 3 - 2DST+1, mean
subplot(2,2,3)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,2,1))');
set(gca,'ydir','normal')
colorbar
title('Lz, 2D+1, mean')

%plot 2 - 2DST+1, median
subplot(2,2,4)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,2,2))');
set(gca,'ydir','normal')
colorbar
title('Lz, 2D+1, median')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% k-wavelength
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2); set(gcf,'color','w')
redyellowblue16

%get amplitude data
% Lambda = quadadd(squeeze(Results.Zonal(2,:,:,:,:,:)),squeeze(Results.Zonal(3,:,:,:,:,:)));
Lambda = squeeze(Results.Zonal(2,:,:,:,:,:));

% %take time mean
% Lambda = 1./squeeze(nanmean(Lambda,1));
%take daily value
Lambda = 1./squeeze(Lambda(DayNumber,:,:,:,:,:));

%plot 1 - 3DST, mean
subplot(2,2,1)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,1,1))');
set(gca,'ydir','normal')
colorbar; caxis([-1,1].*2500)
title('Lk, 3DSt, mean')


%plot 2 - 3DST, median
subplot(2,2,2)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,1,2))');
set(gca,'ydir','normal')
colorbar
title('Lk, 3DSt, median')


%plot 3 - 2DST+1, mean
subplot(2,2,3)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,2,1))');
set(gca,'ydir','normal')
colorbar
title('Lk, 2D+1, mean')

%plot 2 - 2DST+1, median
subplot(2,2,4)
imagesc(Settings.LatScale,Settings.HeightScale,squeeze(Lambda(:,:,2,2))');
set(gca,'ydir','normal')
colorbar
title('Lk, 2D+1, median')