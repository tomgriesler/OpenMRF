function [seq, TRID] = GE_add_TRID(seq, TRID, label_name, flag_GE)
    if flag_GE == 1

        if isempty(TRID)
            TRID.labels = {};
            TRID.all    = {};
        end

        % search if TRID label already exists
        idx = find(strcmp(TRID.labels, label_name), 1);

        % store new label to TRID list
        if isempty(idx)
            TRID.labels{end+1,1} = label_name;
            idx = numel(TRID.labels);
        end

        % store label to history
        TRID.all{end+1,1} = label_name;

        % add TRID segment label to seq object
        seq.addBlock(mr.makeLabel('SET', 'TRID', idx));
    
    end
end
