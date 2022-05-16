function v = sym2vec(M)

    n = size(M,1);

    v = nan(((n+1)*n)/2,1);

    ii = 1;

    for col_index = 1:n

        for row_index = col_index:n

            v(ii) = M(row_index,col_index);

            ii = ii + 1;

        end

    end

end