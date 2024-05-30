function new_df = apply_drex_to_df(old_df, subj)
    %% import modules
    addpath('/Users/scanlab/Documents/internship_luca/model-folder/drex-model/')
    
    % define parameters
    params = [];
    params.distribution = 'gmm';
    params.max_ncomp = 2;
    params.beta = 0.2; % scale down with higher value - probably
    params.D = 1;
    
    params.maxhyp = inf;
    params.memory = inf;
    
    % go through each block
    for run = 1:max(old_df.run)
    for blk = 1:max(old_df.block)
        idxs = find(old_df.block == blk & old_df.run == run);
        
        x = old_df.frequencies_oct(idxs(1):idxs(end));
        
        params.prior = estimate_suffstat(x,params);
        
        out = run_DREX_model(x,params);
        
        old_df.drex_surp(idxs(1):idxs(end)) = out.surprisal;
        
        % display figure of model
%         figure(blk); clf;
%         display_DREX_output(out,x)
    end
    end
    
    new_df = old_df;
    writetable(new_df,append('/Users/scanlab/Documents/internship_luca/Pupil_data_Maastricht/preproData/behav_data/',subj,'_beh.csv'))

    
return