format long; % change format to long for greater precision

%% Algorithm 7.1 Jacobi Iterative

% COMPARED TO OTHER METHODS:
% Jacobi Iterative is an iterative technique to solve the nxn linear system Ax=b that starts 
% with an initial approximation vector x(0) to the solution x and generates
% a sequence of vectors {x(k)} (from k = 0:infinity) that converges to x. A
% residual vector is defined as r = b - Ax', where x' is an
% approximation to the solution of the linear system defined by A=bx. This
% is what we use to check the TOL.

% To solve Ax=b given an initial approximation x(0):

% INPUT: the number of equations and unknowns n; the entries a(i,j), 1<=i,
% j<=n of the matrix A; the entries b(i), 1<=i<=n of bl the entries XO(i),
% 1<=i<=n of XO = x(0); tolerance TOL; max number of iterations N.

% OUTPUT: the approximate solution x(1), ..., x(n) or a message that the
% number of iterations was exceeded.


clear all; % clear any existing variables
data; % uncomment the data you would like to use in the data file
TOL = input ('Enter to the COMMAND WINDOW, TOL,  the tolerance \n'); % TOL allows us to specify when we should stop with the iterations
N = input ('Enter to the COMMAND WINDOW, N,  the maximum number of iterations \n'); % set maximum number of iterations, so if it doesnt converge and never gets to limit the program will stop

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
while k <= N % While k is less than or equal to our max iterations N:
    for i = 1:n % For all entries 1 to n: 
        for j = 1:n % For all entries 1 to n:
            if j ~= i % If j is not equal to i, since we exclude this term in our approximations:
                s = s + a(i,j).*XO(j); % compute s
            end % end if loop
        end % end for loop
        x(i) = (1./a(i,i)).*(-s + b(i)); % calculate the approximate ith iteration solution for vector x
        s = 0; % reset s back to zero
    end % end for loop
    if max( x' - XO ) < TOL % If our change between iterations is less than our tolerance:
        x % print x
        break; % break the while loop
    end % end if loop
    k = k+1; % Otherwise, our tolerance has not been reached, so increment iterative count by one to go through the loop again
    for i = 1:n % For all entries 1 to n:
        XO(i) = x(i)'; % set our x vector equal to XO for the next iteration
    end % end for loop
end % end while loop

fprintf('Maximum number of iterations exceeded \n') % output message to user
XO % print the final approximation we get either within our TOL or our max iterations

% Now we need to calculate the error
A = [a b]; % cocatenate the system of equations
R = rref(A);  % first we need the actual solution to the system, so we calculate the reduced row echelon form of A
for i = 1:n % For entries 1 to n:
    XA(i)=R(i,n+1); % compute the actual solution
end % end for loop

disp('Actual solution: ')
XA = XA' % transpose XA and print the actual solution
% and calculate relative error
ER = max(abs(XA-XO))./max(XA); % calculate error
fprintf('ER: %.6f (reached after %d iterations)\n', ER, k) % print error