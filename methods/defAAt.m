function [A, At] = defAAt(ind, res)
% def sampling operator
    function y = Adef(x, ind)
    y = x(ind);
    end
    function y = Atdef(x, ind, res)
    y = gpuArray.zeros(res);
    y(ind) = x;
    end

A = @(x) Adef(x, ind);
At = @(x) Atdef(x, ind, res);
end
