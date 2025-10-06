close all; clear all;

% program description
disp('This program computes the area under the curve of a function y = e^(dt)*sin(ct) using a chosen numerical method.');

% user's desired method
M = input("Choose an option: \n'1' - for Trapezoidal Rule\n'2' - for Simpson's Rule\n\nWhat is your choice?: ");
if M == 1
   disp('Trapezoidal Rule chosen');
elseif M == 2
   disp("Simpson's Rule chosen");
end

% user value inputs
c = input('Enter value of c for y = e^(dt)*sin(ct) function: ');
d = input('Enter value of d for y = e^(dt)*sin(ct) function: ');
a = input('Enter lower limit for the integration: ');
b = input('Enter upper limit for the integration: ');
n = input('Enter number of intervals: ');

dt = (b - a) / n; % compute delta t
t = a:dt:b;              % t-values
y = exp(d*t) .* sin(c*t); % function values with t-values

% zero crossing points
zc = find(diff(sign(y)));  % finds indices where y changes sign
t(zc); %t-values of the zero crossing points of y

% displays the zero crossing points
fprintf('Zero crossings at t = ');  
fprintf('%.4f ', t(zc));            
fprintf('\n');       





