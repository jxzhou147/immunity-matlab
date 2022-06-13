clear;clc

eq1 = 'DV=I-0.1*V';
eq2 = 'DI=V-M*I-0.1*I';
eq3 = 'DC=M+I-0.1*C';
eq4 = 'DM=I+C-0.1*M';

[V, I, C, M] = dsolve(eq1, eq2, eq3, eq4, 'V(0)=1,I(0)=0,M(0)=0,C(0)=0');