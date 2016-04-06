function [ almost_eq ] = almost_equals( a,b )
%EQUALS Compares two floating-point numbers
treshold = 1;  % Good enough for IGRA coordinates
if (abs(a-b) < treshold)
    almost_eq = 1;
else
    almost_eq = 0;
end
end