function behav_table = creating_behav_df(segmentz, timingz)
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
            prob_cond  = segmentz( 4,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
            width_cond = segmentz( 5,  segmentz( 1,  :)  == run & segmentz( 2,  :)  == block );
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

return