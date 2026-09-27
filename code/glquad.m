%Uday Khankhoje, 23 Jan 2020
function [val,w_new,z_new] = glquad(varargin)
%computes integral of fun from a to b using a n-point GL rule
%a change of variables is needed from [a,b] to [-1,1], gives
%val=(b-a)/2(sum_1^n w_n f( (b-a)z_n/2 + (b+a)/2 )
%w_n,z_n are standard weights and nodes, obtained by calling glquadrule
%Two ways of calling this function, either ask it to calc [w,z],
%or supply it yourself

%these inputs are common
a = varargin{1};
b = varargin{2};
fun = varargin{3};

if nargin == 4
    n = varargin{4};
    [w,z] = glquadrule(n);
elseif nargin == 5
    w = varargin{4};
    z = varargin{5};
    n = length(w);
else
    error('Call the function as (a,b,fun,n) or (a,b,fun,w,z)');
end
%now compute
val = 0;
w_new = w .* (b-a)/2;
z_new = ((b-a).*z+(b+a))/2;
% for i=1:n
%     x = ((b-a)*z(i)+(b+a))/2;
%     val = val + w(i)*fun(x);
% end
% val = (b-a)/2 * val;
for i=1:n
    val = val + w_new(i)*fun(z_new(i));
end
end