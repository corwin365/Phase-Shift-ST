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

%variables we could have varied over. IN ORDER, do not change.
Variables = {'c1','c2','NPeaks','Steps','Weights','Case','RealNoise',  ...
             'Amplitude','Lz','Lx','Ly','RotX','RotY','RotZ'};

%choose a baseline value for each variable.
%Comparisons will be made using the nearest value to this for each var
Baseline = [0.5,0.25,1,2,0,1,0,1,20,500,500,0,0,0];
           
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
  
  %and store
  if ~exist('ComboData'); ComboData = dat.ComboData;
  else                    ComboData = cat(1,ComboData,dat.ComboData);
  end
  
end; clear iFile dat 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% identify which variables we varied
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Varied = NaN(numel(Variables),1);
for iVar = 1:1:numel(Variables); Varied(iVar) = numel(unique(ComboData(:,iVar))); end; clear iVar
Varied = find(Varied > 1); %these are the numbers of the possible variables we varied in this run


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% now, generate comparative plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

