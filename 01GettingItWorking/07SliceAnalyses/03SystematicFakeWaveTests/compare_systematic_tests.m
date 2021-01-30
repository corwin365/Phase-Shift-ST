clearvars


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% assess relative performance of method on systematically-varied parameters
%
%Corwin Wright, c.wright@bath.ac.uk, 17/SEP/2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%master data dir where tests are stored
Settings.DataDir = [LocalDataDir,'/corwin/2dp1/'];

%set of runs we want to assess here
Settings.RunID   = 'tuning1';

%INPUT variables we could have varied over. IN ORDER, do not change.
Settings.Variables = {'c1','c2','NPeaks','Steps','Weights','Case','RealNoise',  ...
                      'Amplitude','Lz','Lx','Ly','RotX','RotY','RotZ'};

%choose a baseline value for each variable.
%Comparisons will be made using the nearest value to this for each var
Settings.Baseline = [0.5,0.25,1,2,0,1,0,1,20,500,500,0,0,0];


%the OUTPUT variables we want to assess. 
%these will be plotted as individual panels varying against the INPUT variables
%we will average them in space, and plot how they vary with height
Settings.OutVars = {'Lz','Lz2','A','A2'};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate a complete list of attempted combinations in this set of runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%need to load files and glue together the ComboData variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%list of files. DO NOT CLEAR, as the order of the files listed will be used later
FileList = wildcardsearch([Settings.DataDir,'/',Settings.RunID,'/'],'*.mat');

%loop over and pull out the info
for iFile=1:1:numel(FileList)
  
  %load combodata list
  dat = load(FileList{iFile},'ComboData');
  
  %add file number to the end
  dat.ComboData = cat(2,dat.ComboData,ones(size(dat.ComboData,1),1).*iFile);
  
  %and the run number
  dat.ComboData = cat(2,dat.ComboData,(1:1:size(dat.ComboData,1))');
  
  %and store
  if ~exist('ComboData'); ComboData = dat.ComboData;
  else                    ComboData = cat(1,ComboData,dat.ComboData);
  end
  
end; clear iFile dat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find the closest value to baseline for each var, and replace the 'baseline' value with that
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iVar=1:1:numel(Settings.Baseline);
  Settings.Baseline(iVar) = ComboData(closest(ComboData(:,iVar),Settings.Baseline(iVar)),iVar);
end; clear iVar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify which variables we varied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Varied = NaN(numel(Settings.Variables),1);
for iVar = 1:1:numel(Settings.Variables); Varied(iVar) = numel(unique(ComboData(:,iVar))); end; clear iVar
Varied = find(Varied > 1); %these are the numbers of the possible variables we varied in this run


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, generate comparative plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%loop over variables
for iVar=1:1:numel(Varied)

  
  
  %extract all the combinations including the baseline for other variables
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Combos = 1:1:size(ComboData(:,1));
  for jVar=1:1:numel(Settings.Variables)
    
    %skip if this is the variable we're working on
    if jVar == Varied(iVar); continue; end
    
    %find all the cases where the variable is the baseline value
    IsBaseline = find(Settings.Baseline(jVar) == ComboData(:,jVar));
    
    %and intersect this with the remaining combolist
    Combos = intersect(Combos,IsBaseline);
    
  end; clear jVar IsBaseline
  
  %and reduce down the list of combinations to just this
  CB = ComboData(Combos,:);
  
  %identify the files containing these tests, and extract them
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Store = struct();
  OldFile = ''; %to speed up the programme by avoiding duplicate file loading
  for iJob=1:1:size(CB,1);
    
    %load the file
    File = FileList{CB(iJob,15)};
    if strcmp(File,OldFile) ~= 1; Data = load(File); OldFile = File; end
    
    %pull out the relevant run
    Store.(['Run',num2str(iJob)]) = Data.DataStore.(['Combo',num2str(CB(iJob,16))]);
    
  end; clear iJob OldFile Data File
  
  %extract the relevant output fields
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Out = NaN(numel(Store.Run1.ret_z),numel(Store),numel(Settings.OutVars));
  x = CB(:,Varied(iVar)); z = Store.Run1.ret_z;
  
  for iJob=1:1:numel(fieldnames(Store))
    for jVar=1:1:numel(Settings.OutVars)
      
      %get var
      Out(:,iJob,jVar) = squeeze(nanmean(Store.(['Run',num2str(iJob)]).(Settings.OutVars{jVar}),[1,2]));
      
% % %       %normalise to input value, then convert to a %age
% % %       switch Settings.OutVars{jVar}
% % %         case {'Lz','Lz2'}; InputValue = CB(iJob,9);
% % %         case { 'A', 'A2'}; InputValue = CB(iJob,8);
% % %         otherwise; disp('Outvar not specified, add as option'); stop
% % %       end
% % %       
% % %       Out(:,iJob,jVar) = Out(:,iJob,jVar)./InputValue.*100;
      
    end
  end
  
  %plot how they vary relative to the input value and height
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  figure(iVar+10)
  sgtitle(Settings.Variables{Varied(iVar)})
  
  
  for jVar=1:1:numel(Settings.OutVars)
    subplot(2,ceil(numel(Settings.OutVars)./2),jVar)
    
    imagesc(x,z,squeeze(Out(:,:,jVar)))
    set(gca,'ydir','normal')
    ylim([20 60])
    colorbar
% %     caxis([-1,1].*100)
    redyellowblue32
    title(Settings.OutVars{jVar})
    drawnow
    
  end
  
  
  
end; clear iVar