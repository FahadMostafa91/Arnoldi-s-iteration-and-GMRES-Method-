function [Q,H] = Arnoldis_iter(A,k_tot,q1)
% Arnoldi’s iteration for square matrix A, k-tot is tolerance
% fahad Mostafa @ Texas tech, 21 OCT, 2020
%==========================================================
% Matrix of dimensions n by n starting with a normalized q1 vector of dimension n
% For k < n returns a Q array of dimensions n by (k+1) with columns 
% orthonormal and an H matrix of dimensions (k+1) per k called 
% Hessenberg matrix which is such a top triangular matrix that
% A*Q(:,1:k) = Q(:,1:k)*H(1:k,1:k) + H(k+1,k)*Q(:,k+1)*E_k'
% Denote E_k is the k-th column of the k-sized identity matrix.
%
% A: square matrix of the linear system (n by n)
% q1: initial normalized vector, assumed not null
% k=k_tot: natural number, total iterations
% Q: Krylov Space Orthonormal Base (n per k+1)
% H: Hessenberg's top triangular matrix (k+1 per k)st A*Q(:,1:k)= Q*H 
% z: is a vector error estimate st AQ = QH + z*E_k

n =length(A);
Q = zeros(n,k_tot); 
q1 = q1/norm(q1);
Q(:,1) = q1;
H = zeros(min(k_tot+1,k_tot),n);
for k=1:k_tot
    z = A*Q(:,k);
    for i=1:k
        H(i,k) = Q(:,i)'*z;
        z = z - H(i,k)*Q(:,i);
    end
    if k < n
       H(k+1,k) = norm(z);
       if H(k+1,k) == 0 
           return, 
       end
       Q(:,k+1) = z/H(k+1,k);
   end
end
       
