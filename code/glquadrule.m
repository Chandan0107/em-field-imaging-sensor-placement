function [weights,quadpts] = glquadrule(n)
% Gauss-Legendre quadrature points and weights on [-1,1] for an n-point rule.
% Numerical (Golub-Welsch) version of the original symbolic implementation by
% Uday Khankhoje (23 Jan 2020); gives the same points and weights to machine
% precision without needing the Symbolic Math Toolbox.
%   weights : 1 x n row vector, quadpts : n x 1 column vector (ascending)
k = 1:n-1;
beta = k ./ sqrt(4*k.^2 - 1);
[V, D] = eig(diag(beta, 1) + diag(beta, -1));
[quadpts, order] = sort(diag(D));
weights = 2 * (V(1, order).^2);
end
