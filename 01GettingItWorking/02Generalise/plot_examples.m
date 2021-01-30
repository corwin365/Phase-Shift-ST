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


Settings.Day = datenum(2008,1,127);
Settings.GiD = [56];
Settings.Fig = 9;

% % Settings.Day = datenum(2010,10,16);
% % Settings.GiD = 186;
% % Settings.Fig = 2;

% % Settings.Day = datenum(2007,1,13);
% % Settings.GiD = [122];
% % Settings.Fig = 3;

% % Settings.Day = datenum(2005,7,8);
% % Settings.GiD = [85];
% % Settings.Fig = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data as a granule-pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get granule
Airs  = cat_struct(prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]),   ...
                   prep_airs_3d(Settings.Day,Settings.GiD+1,'PreSmooth',[3,3,1]), ...
                   2,{'MetaData','Source','ret_z'});

% Airs = prep_airs_3d(Settings.Day,Settings.GiD,  'PreSmooth',[3,3,1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3dst the granule-pair
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ST3D,Airs] = gwanalyse_airs_3d(Airs,                  ...
                               'ZRange',      [10,70], ...
                               'TwoDPlusOne',    true, ...
                               'TwoDPlusOne_ind',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot some maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = 40;
zidx = closest(z,Airs.ret_z);
z = Airs.ret_z(zidx);

figure(Settings.Fig)
clf
set(gcf,'color','w')
redyellowblue32

subplot(2,2,1)
pcolor(Airs.l1_lon,Airs.l1_lat,squeeze(ST3D.A(:,:,zidx))); shading flat; colorbar
title(['Amplitude, 3DST, ',num2str(z),'km altitude'])


subplot(2,2,2)
pcolor(Airs.l1_lon,Airs.l1_lat,1./squeeze(ST3D.m_2dp1_ind(:,:,zidx))); shading flat; colorbar
title(['\lambda_z, 2D+1,ind fit, ',num2str(z),'km altitude'])
caxis([10 40])


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

%create figure
figure(10+Settings.Fig)
clf
set(gcf,'color','w')
redyellowblue32

%average wavelength over strong-signal area
NoSignal = find(ST3D.A < 2.5);

m.ST3D   = abs(ST3D.m);
m.ST2Dv1 = abs(ST3D.m_2dp1_ind);
m.ST2Dv2 = abs(ST3D.m_2dp1);

m.ST3D(  NoSignal) = NaN; m.ST3D   = squeeze(nanmean(m.ST3D,  [1,2]));
m.ST2Dv1(NoSignal) = NaN; m.ST2Dv1 = squeeze(nanmean(m.ST2Dv1,[1,2]));
m.ST2Dv2(NoSignal) = NaN; m.ST2Dv2 = squeeze(nanmean(m.ST2Dv2,[1,2]));




subplot(1,3,1)
axis([0,60,10,70])
hold on
plot(1./m.ST3D,  Airs.ret_z,'color','k')
plot(1./m.ST2Dv1,Airs.ret_z,'color','r')
plot(1./m.ST2Dv2,Airs.ret_z,'color','b')

subplot(1,3,2)
axis([0,10,10,70])
hold on
plot(squeeze(nanmean(ST3D.A,[1,2])),  Airs.ret_z,'color','k')
plot(squeeze(nanmean(ST3D.A_2dp1_ind,[1,2])),Airs.ret_z,'color','r')
plot(squeeze(nanmean(ST3D.A_2dp1,[1,2])),Airs.ret_z,'color','b')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot a slice
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identify which slice has the highest amplitude
[~,idx] = max(squeeze(nansum(ST3D.A,[2,3])));

%create figure
figure(20+Settings.Fig)
clf
set(gcf,'color','w')
redyellowblue16

%plot raw data 
subplot(2,2,1)
axis([0 size(ST3D.A,2) 10 70])
pcolor(1:1:size(ST3D.A,2),Airs.ret_z,squeeze(ST3D.IN(idx,:,:))'); 
shading flat; caxis([-1,1].*prctile(abs(flatten(squeeze(ST3D.IN(idx,:,Airs.ret_z > 25 & Airs.ret_z < 55)))),[97.5]))
title(['T'', track ',num2str(idx)])
colorbar

%plot 3DST output m field
subplot(2,2,2)
axis([0 size(ST3D.A,2) 10 70])
pcolor(1:1:size(ST3D.A,2),Airs.ret_z,squeeze(1./ST3D.m(idx,:,:))'); 
shading flat; caxis([0 40])
title(['\lambda_z, 3DST, track ',num2str(idx)])
colorbar

%plot ind output m field
subplot(2,2,3)
axis([0 size(ST3D.A,2) 10 70])
pcolor(1:1:size(ST3D.A,2),Airs.ret_z,squeeze(1./ST3D.m_2dp1_ind(idx,:,:))'); 
shading flat; caxis([0 40])
title(['\lambda_z, 2D+1, track ',num2str(idx)])
colorbar

%plot match output m field
subplot(2,2,4)
axis([0 size(ST3D.A,2) 10 70])
pcolor(1:1:size(ST3D.A,2),Airs.ret_z,squeeze(1./ST3D.m_2dp1(idx,:,:))'); 
shading flat; caxis([0 40])
title(['\lambda_z, 2D+1 ind fit, track ',num2str(idx)])
colorbar