clear all

[Airs,Spacing] = prep_airs_3d(datenum(2008,1,127),57);

% options = {'minwavelengths',[25 25 6],'maxwavelengths',[10000 10000 45]};
% ST2D = nph_2dst_plus1(Airs.Tp,100,[Spacing,3],ones(1,3).*0.25,options{:});

tic
ST = gwanalyse_airs_3d(Airs,'TwoDPlusOne',true,'Spacing',[Spacing,3]);
toc


figure
pc(ST.A(:,:,10))

figure
pc(ST.A_2dp1(:,:,10))