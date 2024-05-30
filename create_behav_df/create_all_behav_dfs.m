addpath('/Users/scanlab/Documents/internship_luca/Pupil_data_Maastricht')

cwd = '/Users/scanlab/Documents/internship_luca/Pupil_data_Maastricht/testing/eyetracking_data/data/';

for subj = 1:28
    % get filename
    file2load = [cwd,num2str(subj),'/',num2str(subj),'-mainpred.mat'];
    
    % load file
    load(file2load)
    
    % create behav data
    behav_table = creating_behav_df(segmentz, timingz);
    
    % subj name
    if subj <10
        sub_name = ['S0',num2str(subj)];
    else
        sub_name = ['S',num2str(subj)];
    end
    
    % apply drex values
    new_df = apply_drex_to_df(behav_table, sub_name);
end