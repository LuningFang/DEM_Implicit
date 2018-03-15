function y = IsDiagonalDominant(A)

rowSum = sum(abs(A), 2);
diagnal = diag(abs(A));

if all(2*diagnal > rowSum)
    y = 1;
else
    y = 0;
end

end