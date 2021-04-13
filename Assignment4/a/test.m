% Demonstration of the Procrustes method
%
% CISC 472, Queen's University

% Let P = a set of 10 observations, each in a column

n = 10;
   
P = [];
for i = [1:n]
  P = [ P  [rand(); rand(); rand()] ];
end

% Set up some rotation and translation

R1 = makehgtform('axisrotate',[1,2,3],20);
R1 = R1(1:3,1:3);

t1 = [ 0.2; 0.1; -0.1 ];
t1 = repmat( t1, 1, n ); 

% Let Q = rotation and translation applied to P.
%
% Note that the rotation is not necessarily about the centroid of P.
R1
P
t1
Q = R1 * P + t1

% Use Procrustes method to get the rotation

meanP = repmat( mean(P,2), 1, n );
meanQ = repmat( mean(Q,2), 1, n );

Pcentred = P - meanP;
Qcentred = Q - meanQ;

[U, S, V] = svd( Qcentred * Pcentred' );

R2 = U * V';

% Calculate translation

t2 = meanQ - R2 * meanP;

% Apply rotation and translation

Q

P2 = R2 * P + t2

% Check the result, which should have a Frobenius norm of zero

difference = norm( P2 - Q, 'fro' )