function [x,r] = GmresR(A,b,M)
% The generalized minimum waste method solves linear system A x = b
% fahad Mostafa @ Texas tech, 21 OCT, 2020
%==========================================================
% A: non-symmetrical square matrix (n by n)
% k: natural number, total iterations 
% r: rule of (b - Ax), estimation of the order of error committed
% x: approximate solution
%  
n = length(A);
Q = zeros(n,M);  
H = zeros(M,M-1);
% Initial "solution" is zero.
r(1) = norm(b);
v = b/norm(b);
k=M;
[Q,H] = Arnoldis_iter(A,k,v);
  
  % Solve the GMRES problem with the current Q.
  z = (A*Q(:,1:M)) \ b;
  x = Q(:,1:M)*z;
  r = norm( A*x - b );
  
end