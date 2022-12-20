%creates the Lat matrix (matrix of connectivities) to use in fortran
%micromodel.f90
close all
clear all

nodes=7; %total number of nodes in one row of the lattice. 5 for Q1, 7 for Q2, 8 for Q3

I=eye(nodes^2); %creates main diagonal
D1=diag(ones(nodes^2-nodes,1),nodes); %creates right-most diagonal of ones
D2=diag(ones(nodes^2-nodes,1),-nodes); %creates left-most diagonal

dvect=ones(nodes^2-1,1);
for i=1:length(dvect)
    if(mod(i,nodes)==0)
        dvect(i)=0;
    end
end
D3=diag(dvect,1); %creates above-diagonal
D4=diag(dvect,-1); %creates below-diagonal

Lat=D1+D2+D3+D4+I;


%%%%Be sure to change the name of the saved file when you change "nodes"!!
save -ascii LatQ2.dat Lat