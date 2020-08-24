clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make plots for 2020/08/21 meetings base don program outputs
% Corwin Wright, c.wright@bath.ac.uk
% 2020/08/21
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % InFile = 'stout_2008d127g056_thin0.mat';
% % InFile = 'stout_2008d127g056_thin1.mat';
InFile = 'stout_2010d289g186_thin0.mat';
% % InFile = 'stout_2010d289g186_thin1.mat';
% % InFile = 'stout_2010d289g186_thin0.mat';
% % InFile = 'stout_2007d013g122_thin1.mat';


Peak = 1; %which peak in the data (1 = largest)

load(InFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% loop over vars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for iVar=1:1:2;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% get data and prep panel
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  switch iVar
    case 1; 
      Field1 = Airs.Tp; VarName1 = 'Tp'; CAxis1 = [-10,10];                  
      Field2 = Store.A; VarName2 = 'A';  XLim2  = [0,4];
      Field3 = ST3D.A;
    case 2;
      Field1 = abs(squeeze(Store.L(:,:,:,1,3,1)));
      Field1(abs(squeeze(Store.A(:,:,:,1,3,1))) < prctile(flatten(abs(squeeze(Store.A(:,:,:,1,3,1)))),80)) = NaN;
      VarName1 = 'Lz, Step 2, unweighted'; CAxis1 = [0,60];
      Field2 = abs(Store.L);  VarName2 = 'L_z'; XLim2  = [0 80];
      Field3 = 1./ST3D.m;
  end
  
 
  

  figure(iVar)
  clf
  set(gcf,'color','w')
  sgtitle(InFile,'interpreter','none')
  
  k = 0; %counter for panels

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot Tp field at various heights
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  for iZ=[30,40,50];
    
    zidx = closest(Airs.ret_z,iZ);
    k = k+1;subplot(3,3,k)
    
    pcolor(Field1(:,:,zidx));
    shading flat
    colorbar
    redyellowblue16
    caxis(CAxis1)
    title([VarName1,' @ ',num2str(Airs.ret_z(zidx)),'km'])
    drawnow
    
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% find most amplitudey column in 3DST
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% %   zidx = closest(Airs.ret_z,40);
% %   [~,Col] = max(flatten(smoothn(ST3D.A(:,:,zidx),[25,25])));
% %   [x,y] = ind2sub(size(Airs.l1_lat),Col);
  
  
  x = 42; y = 151;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot vertical profile of T'
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  k = k+1; subplot(3,3,[k,k+3]);
  plot(squeeze(Airs.Tp(x,y,:)),Airs.ret_z,'b-x')
  title(['Tp @ [',num2str(x),', ',num2str(y),']'])
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% plot fits
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  for iMode=1:1:2; %[unweighted,weighted]
    
    k = k+1; subplot(3,3,[k,k+3]);
    axis([XLim2 20 60])
    hold on
    
    %get data
    Data = squeeze(Field2(x,y,:,iMode,:,Peak));
    
    %loop over step sizes
    for iStep=0:1:numel(Settings.Steps) %0 omitted as very off-scale
      
      if iStep == 0; Colour = 'k'; Style = '-'; Width = 2; DisplayName = '3DST';
      else
        switch Settings.Steps(iStep)
          case 0; Colour = 'k'; Style = '-'; Width = 1;
          case 1; Colour = 'r'; Style = ':'; Width = 1;
          case 2; Colour = 'r'; Style = '-'; Width = 1;
          case 3; Colour = 'm'; Style = '-'; Width = 1;
          case 4; Colour = 'b'; Style = '-'; Width = 1;
          case 5; Colour = 'g'; Style = '-'; Width = 1;
        end
        if Settings.Steps(iStep) == 0; DisplayName = [num2str(Settings.BasisLevel),'km'];
        else DisplayName = [num2str(Settings.Steps(iStep)),'-step'];
        end
      end
      
      if iStep == 0; Line = squeeze(Field3(x,y,:)); else Line = squeeze(Data(:,iStep)); end
      
      plot(Line,Airs.ret_z,'color',Colour','linestyle',Style','linewidth',Width,'displayname',DisplayName)
      
    end
    
     legend
    drawnow
    if iMode == 1; title([VarName2,' @ [',num2str(x),', ',num2str(y),'], unweighted']);
    else           title([VarName2,' @ [',num2str(x),', ',num2str(y),'], weighted'])
    end
    
  end
 
% %   export_fig([InFile,'_',num2str(iVar)],'-png','-m1','-a3')
  
end
