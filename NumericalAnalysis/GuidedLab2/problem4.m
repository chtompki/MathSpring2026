n=7;
A=magic(n);
x=[1:n]';
b=A*x;
[Q,R]=gs_factor(A);
x1 = qrSolve(Q,R,b);
norm(x - x1)

% problem4
%
% ans =
%
%      4.111978925869794e-14