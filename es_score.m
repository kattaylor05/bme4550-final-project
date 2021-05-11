function [es_scores, max_es, index] = es_score(r_vals, num_hits, indexes, num_genes, p)
    n_r = 0;
    for i = 1:num_genes
        gene_index = indexes(i);
        if ismember(gene_index, num_hits)
            n_r = n_r + abs(r_vals(i))^p;
        end
    end
    hit_index = 0; miss_index = 0;
    es_scores = zeros(num_genes,1);
    for i = 1:num_genes
        gene_index = indexes(i);
        if ismember(gene_index, num_hits)
            hit_index = hit_index + abs(r_vals(i))^p;
        else
            miss_index = miss_index + 1;
        end
        p_hit = hit_index / n_r;
        p_miss = miss_index / (num_genes - 89);
        ES = p_hit - p_miss;
        es_scores(i) = ES;
    end
    [max_es, index] = max(abs(es_scores));
end