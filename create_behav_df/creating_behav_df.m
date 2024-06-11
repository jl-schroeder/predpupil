function behav_table = creating_behav_df(segmentz, timingz, sub_name)
% Define column names
columnNames = {'frequencies', 'frequencies_oct', 'timing', 'timing_offset', 'run', 'block', 'segment','prob_cond','width_cond' , 'center_freq_a', 'center_freq_b', 'center_freq_a_oct', 'center_freq_b_oct', 'probability_a', 'probability_b', 'stimulus'}; 

% Create empty table
behav_table = table('Size',[size(timingz,2) numel(columnNames)], 'VariableTypes', repmat("double", 1, numel(columnNames)), 'VariableNames', columnNames);

behav_table.run              = timingz( 1, :)';
behav_table.block            = timingz( 2, :)';
behav_table.stimulus         = timingz( 3, :)';
behav_table.segment          = timingz( 4, :)';
% behav_table.prob_cond        = segmentz( 4, :)';
% behav_table.width_cond       = segmentz( 5, :)';
behav_table.timing           = timingz( 7, :)'; % temp 5, real 7
behav_table.timing_offset    = timingz( 8, :)'; % temp 6, real 8
behav_table.frequencies_oct  = timingz( 11, :)';
behav_table.frequencies      = timingz( 12, :)';


% empty for adding probability from data
segment_condition = zeros(size(timingz,2),1);
prob_condition = zeros(size(timingz,2),1);
width_condition = zeros(size(timingz,2),1);
probability_a = zeros(size(timingz,2),1);
probability_b = zeros(size(timingz,2),1);
% empty for adding center frequencies from data
center_freq_a = zeros(size(timingz,2),1);
center_freq_b = zeros(size(timingz,2),1);

for run = 1:max(behav_table.run)
    for block = 1:max(behav_table.block)
        for seg = 1:max(behav_table.segment)
            % probabilities
            prob_a = segmentz( 6,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            prob_b = segmentz( 7,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block ); 
            probability_a(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = prob_a(seg);
            probability_b(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = prob_b(seg);
            
            % conditions
            seg_cond   = segmentz( 3,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            prob_cond  = segmentz( 4,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            width_cond = segmentz( 5,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            segment_condition(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = seg_cond(seg);
            prob_condition(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = prob_cond(seg); 
            width_condition(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = width_cond(seg); 
            
            % center frequencies
            c_freq_a = segmentz( 9,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            c_freq_b = segmentz( 10,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            center_freq_a(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = c_freq_a(seg);
            center_freq_b(behav_table.block == block &  behav_table.run == run & behav_table.segment == seg) = c_freq_b(seg);
        end 
    end
end

behav_table.prob_cond        = prob_condition;
behav_table.width_cond       = width_condition;
behav_table.center_freq_a    = 2.^center_freq_a;
behav_table.center_freq_b    = 2.^center_freq_b;
behav_table.center_freq_a_oct= center_freq_a; 
behav_table.center_freq_b_oct= center_freq_b;
behav_table.probability_a    = probability_a; 
behav_table.probability_b    = probability_b;


general_cond_names = {'REG','RAND','d'};
cond_names = general_cond_names(segmentz( 3,  :));
columnNames_seg = {'cond_names', 'cond_vals' 'probability_a', 'probability_b','width','center_a','center_b','jitter_pre','jitter_post'}; 

% Create empty table
segment_table = table('Size',[size(segmentz,2) numel(columnNames_seg)], 'VariableTypes', repmat("double", 1, numel(columnNames_seg)), 'VariableNames', columnNames_seg);
segment_table.cond_names    = cond_names';
segment_table.cond_vals     = segmentz( 3,  :)';
segment_table.probability_a = segmentz( 6,  :)';
segment_table.probability_b = segmentz( 7,  :)';
segment_table.width         = segmentz( 5,  :)';
segment_table.center_a      = segmentz( 9,  :)';
segment_table.center_b      = segmentz( 10,  :)';
segment_table.jitter_pre    = segmentz( 11,  :)';
segment_table.jitter_post   = segmentz( 14,  :)';

writetable(segment_table,append('/Users/scanlab/Documents/internship_luca/Pupil_data_Maastricht/preproData/behav_data/segment_info/',sub_name,'_segInfo.csv'))

return