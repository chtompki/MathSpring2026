format long;

syms w x y z;
q = [x y z]';
A = [16 -28 0; -28 53 10; 0 10 29];
q'*A*q

r = [w x y z]';
B = [1 -2 3 -2; -2 20 -2 8; 3 -2 11 -5; -2 8 -5 9]
r'*B*r

K = cholesky(A)
L = cholesky(B)

% >> ps2prob7
%
% ans =
%
% conj(y)*(53*y - 28*x + 10*z) + conj(x)*(16*x - 28*y) + conj(z)*(10*y + 29*z)
%
%
% B =
%
%      1    -2     3    -2
%     -2    20    -2     8
%      3    -2    11    -5
%     -2     8    -5     9
%
%
% ans =
%
% conj(w)*(w - 2*x + 3*y - 2*z) - conj(x)*(2*w - 20*x + 2*y - 8*z) + conj(y)*(3*w - 2*x + 11*y - 5*z) - conj(z)*(2*w - 8*x + 5*y - 9*z)
%
%
% K =
%
%      4     0     0
%     -7     2     0
%      0     5     2
%
%
% L =
%
%      1     0     0     0
%     -2     4     0     0
%      3     1     1     0
%     -2     1     0     2