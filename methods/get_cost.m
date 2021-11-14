function cost = get_cost(p, ev, epsmin, A, x, b, lambda)
%cost computations (of previous iterate)
if p == 0
    shatten = 0.5 * sum(log(ev+epsmin));
else
    shatten = (1 / p) * sum((ev + epsmin).^(p / 2));
end
diff = A(x) - b;
if (lambda == 0) %equality constrained
    cost = shatten;
else
    cost = norm(diff(:)).^2 + lambda * shatten;
end

end