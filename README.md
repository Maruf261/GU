# GU
Gauss Elimination
%Created on 5/4/18 2:40pm
clc
clear all
syms x y z w
%% Converting Linear Simultaneous Equations to Matrix form
eqns=[-x+2*y+2*z-3*w==-1,
    x+2*z+3*w==1,
    6*x+2*y+2*z+4*w==1,
    y+z+4*w==2];
vars=[ x y z w];
% a is the matrix of co-efficient
% b is the matrix of constants
[a,b]=equationsToMatrix(eqns,vars)
% m is the row
% c is the column
[m,c]=size(a);
j=1; %j th column
k=1;%k th row
n=1;% n th pivoting row
z=2;% Row operation always starts from 2
%% Forward Elimination
% Shifting Loop for column operation
for i=1:c-1
    %% Partial Pivoting
    [x MErow]=max(abs(a(j:m,i)))
    % x is the hightest element of i th column
    % MErow is the iteration of getting Maximum element
    MErow=MErow+j-1 % Maximum element is in the MErow th row
    a([j MErow],:)=a([MErow j],:)
    b([j MErow],:)=b([MErow j],:)
    %% Making each desired column component zero
for r=z:m
    % Checking if any Lower triangle is Already Zero or not
    if a(r,j)==0
        % If any is zero left the entire row as it is and skip the step
        a(r,:)=a(r,:);b(r,:)=b(r,:);
    else
        b(r,:)=((a(r,j)/a(k,j))*b(n,:))-b(r,:)
        a(r,:)=((a(r,j)/a(k,j))*a(n,:))-a(r,:)
    end 
end
k=k+1;% Changing row after completion of inner loop
n=n+1;% Changing the pivoting row
z=z+1;% Setting a new condition for inner loop
j=j+1;% Changing column after completion of inner loop
end
%% Performing Back Substitution
y=linsolve(a,b)
% y is the matrix of the required result
