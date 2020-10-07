format long;
%% Algorithm 7.2 Gauss-Seidel Iterative

% COMPARED TO OTHER METHODS:
%A possible improvement in Algorithm 7.1 can be seen by reconsidering our
% equations. The components of x(k-1) are used to compute the componenets
% of x(k). But for i > 1, the compoenent of x1(k),...,xi-1(k) of x(k) have
% already been computed. Using these instead of x(k-1), we can get a better
% approximation.

% To solve Ax=b given an initial approximation x(0):

% INPUT: the number of equations and unknowns n; the entries a(i,j), 1<=i,
% j<=n of the matrix A; the entries b(i), 1<=i<=n of bl the entries XO(i),
% 1<=i<=n of XO = x(0); tolerance TOL; max number of iterations N.

% OUTPUT: the approximate solution x(1), ..., x(n) or a message that the
% number of iterations was exceeded.

clear all; % clear any existing variables
data; % uncomment the data you would like to use in the data file
a = [1 0 -2; -.5 1 -.25; 1 -.5 1];
n = rank(a); %dimensions of our matrix
b = [.2; -1.425; 2];
XO = [0; 0; 0];
TOL = input ('Enter to the COMMAND WINDOW, TOL,  the tolerance \n'); % TOL allows us to specify when we should stop with the iterations
N = input ('Enter to the COMMAND WINDOW, N,  the maximum number of iterations \n'); % set maximum number of iterations, so if it doesnt converge and never gets to limit the program will stop

k = 1; % let k be 1 for the first iteration
s = 0; % let s be 0. We will use this variable for summation functions
p = 0; % let p be 0. We will use this variable for summation functions
while k <= N % While k is less than or equal to our max iterations N
    for i = 1:n % For entries 1 to n:
        for j = 1:i-1 % For entries 1 to i-1:
            s = s + a(i,j).*x(j); % compute the first summation
        end % end for loop
        for j = i+1:n % For entries i+1 to n:
            p = p + a(i,j).*XO(j); % compute second summation
        end % end for loop
        x(i) = (1./a(i,i)).*(-s - p + b(i)); % calculate the values in our approximated x vector
        s = 0; % reset s back to zero
        p = 0; % reset p back to zero
    end % end for loop
    if max(x' - XO) < TOL % If our change between iterations is less than our tolerance:
        x % print the approximation
        break; % break while loop
    end % end if loop
    k = k+1; % Otherwise, our tolerance has not been reached, so increment iterative count by one to go through the loop again
    for i = 1:n % For entries 1 to n:
        XO(i) = x(i)'; % set our x vector equal to XO for the next iteration
    end % end for loop
end % end while loop
fprintf('Maximum number of iterations exceeded \n') % output message to user
XO % print final approximation we get with our tolerance and max iterations

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