clear all
close all
clc

a = [1 2 3 5 6]; 
b = [4 5 6 7];
c = [7 8 9];

nc = length(c);
d=a'*b;
e = cat(3,d*c(1),d*c(2),d*c(3));

%e= repmat (d,nc,1);
