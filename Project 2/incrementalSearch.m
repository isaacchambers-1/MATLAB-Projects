function [x_l, x_u] = incrementalSearch(f, root_number)
    root_located = 0;
    x_i = -1000;
    max_iter = 1e6;
    count = 0;

    while (root_located ~= root_number && count < max_iter)
        first = f(x_i);
        second = f(x_i+1);

        if first == 0
            root_located = root_located + 1;
            x_l = x_i;
            x_u = x_i;
        elseif first * second < 0
            root_located = root_located + 1;
            x_l = x_i;
            x_u = x_i + 1;
        end

        x_i = x_i + 1;
        count = count + 1;
    end

    if count == max_iter
        error('not found');
    end
end
