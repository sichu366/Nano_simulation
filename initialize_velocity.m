function [v]=initialize_velocity(K_B,N,D,T,m)
% K_B: Boltzmann's constant 8.6173303 x 10-5 eV K-1
% N : number of atoms in the system
% D : dimension, which is 3
% T : temperature in kelvin
% m(i) is the mass of atom i
% v(i,d) is the velocity of atom i in the d-th direction
v=rand(N,3)-0.5;
momentum_average=zeros(1,D);
for d=1:D
    momentum_average(d)=sum(v(:,d).*m)/N;
end
for n=1:N
    v(n,:)=v(n,:)-momentum_average/m(n);
end
v=v*sqrt(T*D*K_B*N/sum(m.*sum(v.^2,2))); % scale velocity
