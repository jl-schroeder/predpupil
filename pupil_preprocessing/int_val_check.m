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


% subs = [1, 2, 71, 75, 106, 110, 124, 204, 270, 302, 304, 307, 328]; % 154,
subs = [1, 300, 303, 339, 365, 47, 153, 203, 237, 266, 261];
% 1-1, 4-5, 4-8, 6-11, 8-6, 11-8, 18-2, 20-3, 22-6, 24-3, 24-13
perc_min = 10;
perc_max = 60;

avg_int_data = zeros(size(subs,2),size(perc_min:perc_max,2));
avg_int_block_data = zeros(size(subs,2),100);

for s = 1:size(subs,2)%size(sorted_file_names,2)%
%     disp(subs(s)) 343

%% load data
% file2load = [rel_dir, sorted_file_names{subs(s)}];
file2load = [rel_dir, sorted_file_names{s}];

T = readtable(file2load);
% T = sortrows(T,'Var1','ascend');
T = sortrows(T,'timestamp','ascend');

% interpolate missing time data
% T.time = fillmissing(T.time,'linear');
T.time = T.time_new;

T.time = T.time/1000;

% plot runs
% plot(T.time, T.pupil)
% png_name = [rel_dir,'pngs/' , sorted_file_names{s}(1:end-4), '.png'];
% saveas(gcf,png_name)
% end
% blk_start   = find(contains(T.Marker, 'B'));
% blk_end     = find(contains(T.Marker, 'E'));
% 
% blk_pick    = 1;
% 
% orig_data   = T.pupil(blk_start(blk_pick):blk_end(blk_pick));
% time        = T.time(blk_start(blk_pick):blk_end(blk_pick));

% plot(time, orig_data)

%% make fieldtrip structure
data = struct;
data.label      = {'pupilData'};     % cell-array containing strings, Nchan*1

time_diff           = diff(time);
mean_time_diff      = mean(time_diff);
rounded_frequency   = round(1 / mean_time_diff,2);%5);*1000;

data.fsample    = rounded_frequency;  % sampling frequency in Hz, single number
data.trial      = {reshape(orig_data, 1, [])};   % cell-array containing a data matrix for each
                                 % trial (1*Ntrial), each data matrix is a Nchan*Nsamples matrix
data.time       = {reshape(time, 1, [])};   % cell-array containing a time axis for each
                                 % trial (1*Ntrial), each time axis is a 1*Nsamples vector
% data.trialinfo  = {reshape(T.Marker, 1, [])};
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
% 
data.hdr        = hdr; % save header to data


cfg = [];
cfg.continuous      = 'yes';%'no';
cfg.channel         = 'pupilData';

% if matches(continuous,'no')
cfg.lpfilter        = 'no';% 'no' or 'yes'
cfg.lpfreq          = 20;% lowpass  frequency in Hz
cfg.lpfilttype      = 'but';% digital filter type, 'but' or 'firws' or 'fir' or 'firls' (default = 'but')
% end
data = ft_preprocessing(cfg,data);


%% check random interpolation

int_vals = zeros(1,size(perc_min:perc_max,2));

for i = perc_min:perc_max
    
    indices_arr    = 1:size(orig_data,1);
    n              = length(indices_arr);
    nr             = round((n/100)*i);

    data2del       = indices_arr(randperm(n,nr));

    % cut out data
    data_cut = data;
    data_cut.trial{1}(data2del) = NaN;

    backup_pp = 1;
    while isnan(data_cut.trial{1}(1))
        backup_pp = backup_pp+1;
        data_cut.trial{1}(1) = data_cut.trial{1}(backup_pp);
    end

    backup_pp = size(data_cut.trial{1},2);
    while isnan(data_cut.trial{1}(end))
        backup_pp = backup_pp-1;
        data_cut.trial{1}(end) = data_cut.trial{1}(backup_pp);
    end

    % interpolate
    cfg                         = [];
    cfg.method                  = 'makima';
    cfg.prewindow               = 0.012;%0.005;
    cfg.postwindow              = 0.012;%0.005;
    cfg.feedback                = 'text';
    data_interpolated           = ft_interpolatenan(cfg, data_cut);

    % check corr
    % downsample before to 40 Hz
    cfg             = [];
    cfg.resamplefs  = 40;
    cfg.detrend     = 'no';
    data_40 = ft_resampledata(cfg, data);
    data_interpolated_40 = ft_resampledata(cfg, data_interpolated);

    corr = corrcoef(data_40.trial{1}, data_interpolated_40.trial{1});

    int_vals(i-(perc_min-1)) = corr(1,2);
end

% save in data over subjects
avg_int_data(s,:) = int_vals;

%% check longer missing data
% 1 time point is approx. 0.0044 ms
% check in 50 ms steps

int_block_vals = zeros(1,100);
for i = 1:100
    time2rem = 0.050*i;
    points2rem = round(time2rem/0.0044);
    
    % cut out data
    data_cut = data;
    start_point = 100;
    data_cut.trial{1}(start_point:(start_point+points2rem)) = NaN;
     
    % interpolate
    cfg                         = [];
    cfg.method                  = 'makima';
    cfg.prewindow               = 0.012;%0.005;
    cfg.postwindow              = 0.012;%0.005;
    cfg.feedback                = 'text';
    data_interpolated2           = ft_interpolatenan(cfg, data_cut);
    
    % downsampling
    cfg             = [];
    cfg.resamplefs  = 40;
    cfg.detrend     = 'no';

    data_40 = ft_resampledata(cfg, data);
    data_interpolated2_40 = ft_resampledata(cfg, data_interpolated2);
    % check corr
    corr = corrcoef(data_40.trial{1}, data_interpolated2_40.trial{1});
    
    int_block_vals(i) = corr(1,2);
end

avg_int_block_data(s,:) = int_block_vals;


end

figure(1)
plot(perc_min:perc_max,mean(avg_int_data,1))
title(['Average correlation (', num2str(size(subs,2)),' runs) for different random NaNs'])
xlabel('Random NaN (in %)') 
ylabel('correlation')

figure(2)
plot(.05:0.05:5,mean(avg_int_block_data,1))
title(['Average correlation (', num2str(size(subs,2)),' runs) for different lengths of NaN blocks'])
xlabel('Random NaN block len (in s)') 
ylabel('correlation')

% figure(3)
% plot(time, data.trial{1}, 'DisplayName','original')
% hold on
% plot(time, data_interpolated.trial{1}, 'DisplayName','35% rand')
% plot(time, data_interpolated2.trial{1} ,'DisplayName','1.5 sec NaN')
% colororder(["red";"blue";"green"])
% legend
% hold off