function [xf, y] = newt(f,xo);
% a newton method that will solve f=0, where f is a function that
% returns two parametes,f(x), f'(x) AS A CELL ARRAY!
% singleton functions only?

xf = xo;
for iter = 1:10
    g = f(xf);
    xf = xf - g{1}/g{2};
end
