function [f,U]=find_force(N,D,NN,NL,L,pbc,r)
% N: number of atoms in the system
% D: space-dimension, which should be 3
% NN(N,1): NN(i) is the number of neighbors of atom i
% NL(N,:): NL(i,k) is the index of the k-th neighbor of atom i
% L(1,3): L(d) is the box length in the d-th direction
% pbc(1,3): pbc(d)=1(0) means periodic (free) in the d-th direction
% r(N,3): r(i,d) is the position of atom i in the d-th direction
% f(N,3): f(i,d) is the total force on atom i in the d-th direction
% U: total potential energy of the system
EPSILON=1.032e-2; % in units of eV (only for Argon)
SIGMA=3.405; % in units of Angstrom (only for Argon)
sigma_6=SIGMA^6;sigma_12=SIGMA^12;L_times_pbc=L.*pbc;
U=0; % initialize the total potential energy of the system
f=zeros(N,D); % initialize the total forces on each atom
for n1=1:N-1 % loop over the atoms
    for m=1:NN(n1) % loop over the neighbors of atom n1
        n2=NL(n1,m); % (n1, n2) is a pair of atoms
        if n2<n1;continue;end % Use Newton's 3rd law to speed up
        r12=r(n2,:)-r(n1,:);
        r12=r12-round(r12./L).*L_times_pbc; % minimum image convention
        d12_square=sum(r12.*r12);
        d_6=d12_square^3;d_8=d12_square*d_6;d_12=d_6*d_6;d_14=d_6*d_8;
        f12=(sigma_6/d_8-2.0*sigma_12/d_14)*24.0*EPSILON;
        f(n1,:)=f(n1,:)+f12*r12;
        f(n2,:)=f(n2,:)-f12*r12; % Newton's 3rd law used here
        U=U+4.0*EPSILON*(sigma_12/d_12-sigma_6/d_6); % accumulate energy
    end
end

