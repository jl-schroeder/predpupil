%% setup
addpath('/Users/scanlab/Documents/internship_luca/Toolboxes/fieldtrip-20240405')
ft_defaults

%% load data

% Define the relative directory
rel_dir = '/Users/scanlab/Documents/internship_luca/Pupil_data_Maastricht/pupilFieldtrip/blocks/predown/';


% Get a list of all .txt files in the directory
file_structs = dir(fullfile(rel_dir, '*.csv'));

% Extract the names of the files from the struct array
file_names = {file_structs.name};

% Sort the file names
sorted_file_names = sort(file_names);

% Display the sorted file names
% disp(sorted_file_names);

%% go through whole subject runs
for subj_run = 1:length(sorted_file_names)

% load data
continuous = 'yes'; % 'yes'/'no'
    
file2load = [rel_dir, sorted_file_names{subj_run}];
disp(sorted_file_names{subj_run})

T = readtable(file2load);

% set time, from time-accurate downsampling, as new time-column
T.time = T.time_new;

% move time to sec
T.time = T.time/1000;

%% put data into fieldtrip structure
data = struct;
data.label      = {'pupilData'};     % cell-array containing strings, Nchan*1

time_diff           = diff(T.time);
mean_time_diff      = mean(time_diff);
rounded_frequency   = round(1 / mean_time_diff,2);%5);*1000;

data.fsample    = rounded_frequency;  % sampling frequency in Hz, single number
data.trial      = {reshape(T.pupil, 1, [])};   % cell-array containing a data matrix for each
                                 % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
data.time       = {reshape(T.time, 1, [])};   % cell-array containing a time axis for each
                                 % trial (1*Ntrial), each time axis is a 1*Nsamples vector
data.trialinfo  = {reshape(T.Marker, 1, [])};
% data.sampleinfo = [T.time(1), T.time(end)];% optional array (Ntrial*2) containing the start and end
                % sample of each trial
data.dimord = 'chan_time';
             
% create header
hdr             = [];
hdr.Fs          = rounded_frequency;% sampling frequency
hdr.nChans      = 1;    % number of channels
hdr.nSamples    = length(data.trial);% number of samples per trial
hdr.nSamplesPre = 0;% number of pre-trigger samples in each trial
hdr.nTrials     = 144;  % number of trials
hdr.label       = data.label; % Nx1 cell-array with the label of each channel
hdr.chantype    = {'misc'};% Nx1 cell-array with the channel type, see FT_CHANTYPE
hdr.chanunit    = {'unknown'};% Nx1 cell-array with the physical units, see FT_CHANUNIT

data.hdr        = hdr; % save header to data

%% cut data into trials
if matches(continuous,'no')
    blk_start = find(strcmp(T.Marker, 'B'));
    blk_end = find(strcmp(T.Marker, 'E'));

    cfg.trl = zeros(length(blk_start), 3);
    cfg.trl(:, 1) = blk_start;
    cfg.trl(:, 2) = blk_end;
    cfg.trl(:, 3) = -1;

    trl = cfg.trl;

    trialdata = ft_redefinetrial(cfg, data);
else
    trialdata = data;
end

%% pre-processing
cfg = [];
cfg.continuous      = continuous;%'no';
cfg.channel         = 'pupilData';
cfg.lpfilter        = 'no';% 'no' or 'yes'
cfg.lpfreq          = 20;% lowpass  frequency in Hz
cfg.lpfilttype      = 'but';% digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')

preprodata = ft_preprocessing(cfg,trialdata);

%% artifact detection
% EOG
cfg            = [];
cfg.continuous = continuous;%'no';

% if threshhold file exists, load, else, manual selection
s_r_name = regexp(file2load, 'S(.*?)\.csv', 'match', 'once');

% set filename of threshoulding csv
fileName_tresh = [rel_dir,'tresh_storage/tresh_',s_r_name];

% check if file exists, otherwise set default value and manually check
if exist(fileName_tresh, 'file')
    cfg.artfctdef.zvalue.cutoff      =  readmatrix(fileName_tresh);
    cfg.artfctdef.zvalue.interactive = 'no';    % feedback
else
    tresh = 0.0015;
 
    cfg.artfctdef.zvalue.cutoff      = tresh;   % number, z-value threshold
    cfg.artfctdef.zvalue.interactive = 'yes';   % feedback
end

% channel selection, cutoff and padding
cfg.artfctdef.zvalue.channel     = 'pupilData';
cfg.artfctdef.zvalue.trlpadding  = 0;
cfg.artfctdef.zvalue.artpadding  = 0.1;
cfg.artfctdef.zvalue.fltpadding  = 0;

% algorithmic parameters
cfg.artfctdef.zvalue.bpfilter   = 'yes';
cfg.artfctdef.zvalue.bpfilttype = 'but';
cfg.artfctdef.zvalue.bpfreq     = [2 15];
cfg.artfctdef.zvalue.bpfiltord  = 4;
cfg.artfctdef.zvalue.hilbert    = 'yes';

cfg.artfctdef.zvalue.zscore = 'no';

[cfg, artifact_eog] = ft_artifact_zvalue(cfg, preprodata);

% save threshold value in file
tresh = cfg.artfctdef.zvalue.cutoff;
writematrix(tresh,fileName_tresh,'Delimiter',',');

%% reject artifacts
cfg                           = [];
cfg.artfctdef.reject          = 'nan'; % 'complete' rejects complete trials, use 'partial' if you want to do partial artifact rejection
cfg.artfctdef.eog.artifact    = artifact_eog;

data_no_artifacts = ft_rejectartifact(cfg,preprodata);

% replace values under certain size to nan, see e.g. S2-1
for dpt= 1:size(data_no_artifacts.trial{1},2)
    if data_no_artifacts.trial{1}(dpt) < 0.005
       data_no_artifacts.trial{1}(dpt) = NaN; 
    end
end

%% define int percentage based on missing data
% check % of missing data
blk_start = find(contains(data_no_artifacts.trialinfo{1}, 'B'));
blk_end = find(contains(data_no_artifacts.trialinfo{1}, 'E'));

% estimate and replace missing block end information
if size(blk_end,2) < 9
    blk_end2 = zeros(1,9);
    for blk = 1:8
        if blk_end(blk) > blk_start(blk+1)

            blk_end2(1:blk-1) = blk_end(1:blk-1);
            blk_end2(blk) = blk_start(2+1)-80;
            
            blk_end2(blk+1:end) = blk_end(blk:end);
        end
    end
    blk_end = sort(blk_end2);
end

% get amount of missing data in block
int_bool_blk = zeros(1,9);
for blk = 1:9
    idx_s = blk_start(blk);
    idx_e = blk_end(blk);
    
    nan_vals = sum(isnan(data_no_artifacts.trial{1}(idx_s:idx_e)));
    
    % set flag to 1 if percentage of NaN in block bigger than 35 percent
    if nan_vals/size(data_no_artifacts.trial{1}(idx_s:idx_e),2) > 0.35
        int_bool_blk(blk) = 1;
    end
    
    max_nan_len = 0;
    for p_val = idx_s:idx_e
        
        if isnan(data_no_artifacts.trial{1}(p_val))
            temp_max_nan_len = 1;
            
            if p_val < idx_e
                p_temp = p_val;
                while isnan(data_no_artifacts.trial{1}(p_temp+1))
                    temp_max_nan_len = temp_max_nan_len + 1;

                    p_temp = p_temp+1;

                    if p_temp == idx_e
                        break
                    end
                end
            end
            
            % take temporary block of nan len if longer than saved
            if temp_max_nan_len>max_nan_len
                max_nan_len = temp_max_nan_len;
            end
        end
    end
    
    % set flag to 1 if NaN - blocks longer than 1.5 seconds
    if max_nan_len > floor(1.5/mean(diff(data_no_artifacts.time{1})))
        int_bool_blk(blk) = 1;
    end
end

fileNameBlk = [rel_dir,'int_rateNaNs/int_',regexp(file2load, 'S(.*?)\.csv', 'match', 'once')];
writematrix(int_bool_blk,fileNameBlk,'Delimiter',',');

% end
%% interpolate missing data
cfg                         = [];
cfg.method                  = 'makima';%'spline';% string, interpolation method, see INTERP1 (default = 'linear') > spline, makima
cfg.prewindow               = 0.012;
cfg.postwindow              = 0.012;
cfg.feedback                = 'text';% string, 'no', 'text', 'textbar', 'gui' (default = 'text')

data_interpolated           = ft_interpolatenan(cfg, data_no_artifacts);

% downsample to 40 Hz
cfg             = [];
cfg.resamplefs  = 40;
cfg.detrend     = 'no';
data_40 = ft_resampledata(cfg, data_interpolated);

%% save data without downsampling
time_data = data_40.time{1};
trial_data = data_40.trial{1};

event_data_40 = strings(1,size(trial_data,2));
event_data = data_interpolated.trialinfo{1};

blk_start = find(contains(data_interpolated.trialinfo{1}, 'B'));
blk_end = find(contains(data_interpolated.trialinfo{1}, 'E'));

% estimate and replace missing block end information
if size(blk_end,2) < 9
    blk_end2 = zeros(1,9);
    for blk = 1:8
        if blk_end(blk) > blk_start(blk+1)

            blk_end2(1:blk-1) = blk_end(1:blk-1);
            blk_end2(blk) = blk_start(2+1)-80;

            blk_end2(blk+1:end) = blk_end(blk:end);
        end
    end
    blk_end = sort(blk_end2);
end

% get data from non-downsampled struct
for blk = 1:9
    t_s = data_interpolated.time{1}(blk_start(blk));   
    t_e = data_interpolated.time{1}(blk_end(blk));

    [minVal_s,idx_s] = min(abs(time_data-t_s));
    [minVal_e,idx_e] = min(abs(time_data-t_e));

    event_data_40(idx_s) = 'B';
    event_data_40(idx_e) = 'E';
end

% build output data
outputTable = table(reshape(time_data, [],1), reshape(trial_data, [],1), reshape(event_data_40, [],1), 'VariableNames',["time","pupil","event"] );

% save output data
match = regexp(file2load, 'S(.*?)\.csv', 'match', 'once');
fileName = [rel_dir,'fieldtrip/',regexp(file2load, 'S(.*?)\.csv', 'match', 'once')];
writetable(outputTable,fileName,'Delimiter',',');

% downsample data without interpolation
cfg             = [];
cfg.resamplefs  = 40;
cfg.detrend     = 'no';
data_noint_40 = ft_resampledata(cfg, preprodata);

% make and save plot to manually check correctness of interpolation
figure1 = figure;
axes1 = axes('Parent',figure1);
hold(axes1,'all');

plot(time_data, data_noint_40.trial{1})
plot(time_data, trial_data)
ylim([0 0.06])
legend('raw-data','inter-data')

int_plt_name = [rel_dir, 'pngs/ints/',sorted_file_names{subj_run}(1:end-4),'.png'];
saveas(figure1,int_plt_name)


hold off
close all
end