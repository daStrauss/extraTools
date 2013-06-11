function f = factorizedLU(A)
  % f = factorizedLU(A)
% can I produce a factorized function just like I had in python?
% YES! wow.
  % returns a function that solves the system Ax = b. that is f(b) =
  % x;


[L U P Q R] = lu(A);

f = @(x) Q*(U\(L\(P*(R\x))));
