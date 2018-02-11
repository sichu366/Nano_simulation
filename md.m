function [E,Position,Num_Neighbor,Neighbor_List]=md(Dimension,n0,nxyz,lattice_cons,pbc,Ne,Np,Ns,Cutoff_dis,dt,Temp_k,time2test)
% unit in my system: energy(eV); length(Angstrom); mass(atomic mass) 
% Dimension:  which should be 3
% n0: which is 4 ,it is number of atoms in a cubic unit cell
% nxyz(d) : number of unit cells in the d-th direction
% lattice_cons(1,3): lattice_cons(d) is the lattice constant in the d-th direction
% pbc(1,3): pbc(d)=1(0) means periodic (free) in the d-th direction
% Ne: number of time steps in the equilibration stage
% Np: number of time steps in the production stage
% Cutoff_dis: cutoff distance
% dt: time step of integration in units of fs
% Temp_k: temperature prescribed in units of K
% E(:,1): potenital eneryg per atom at regular time points
% E(:,2): kinetic eneryg per atom at regular time points
% E(:,3): total eneryg per atom at regular time points
K_B=8.625e-5; % Boltzmann's constant in my unit system  
TIME_UNIT_CONVERSION=10.18; % from fs to my unit system
N=n0*nxyz(1)*nxyz(2)*nxyz(3) % number of atoms
b =lattice_cons ;
lattice_cons = lattice_cons*[1,1,1];
BoxSize=lattice_cons.*nxyz; % box size (Angstrom)
dt=dt/TIME_UNIT_CONVERSION % time step in my unit system
m=ones(N,1)*40; % mass for Argon atom in my unit system
Position=initialize_position(N,Dimension,n0,nxyz,lattice_cons); % intial positions
h = figure;
plot3(Position(:,1),Position(:,2),Position(:,3),'o')
grid on

str = 'Init position,T=';
Tstr = num2str(Temp_k);
str = strcat(str,Tstr);
str = strcat(str,',Box side len=');
boxStr = num2str(b);
str = strcat(str,boxStr);
title(str,'fontsize',18);
basepath = pwd;
path_save = strcat(basepath,num2str(time2test));
path_save = strcat(path_save,str);
path_save = strcat(path_save,'.jpg');

saveas(h,path_save);

h = figure;
plot3(Position(1,1),Position(1,2),Position(1,3),'ro')
plot3(Position(2,1),Position(2,2),Position(2,3),'go')
plot3(Position(3,1),Position(3,2),Position(3,3),'bo')
plot3(Position(4,1),Position(4,2),Position(4,3),'co')
plot3(Position(5,1),Position(5,2),Position(5,3),'yo')
plot3(Position(6,1),Position(6,2),Position(6,3),'ko')
hold on
velocity=initialize_velocity(K_B,N,Dimension,Temp_k,m) ;% initial velocities
[Num_Neighbor,Neighbor_List]=find_neighbor(N,BoxSize,pbc,Cutoff_dis,Position) % fixed neighbor list

[f]=find_force(N,Dimension,Num_Neighbor,Neighbor_List,BoxSize,pbc,Position); % initial forces
E=zeros(Np/Ns,3); % energy data to be computed
for step=1:(Ne+Np) % time-evolution started
    for d=1:Dimension % step 1 of Velocity-Verlet
        velocity(:,d)=velocity(:,d)+(f(:,d)./m)*(dt*0.5); 
        Position(:,d)=Position(:,d)+velocity(:,d)*dt; 
    end
    [Num_Neighbor,Neighbor_List]=find_neighbor(N,BoxSize,pbc,Cutoff_dis,Position);% fixed neighbor list
    [f,U]=find_force(N,Dimension,Num_Neighbor,Neighbor_List,BoxSize,pbc,Position); % update forces
    for d=1:Dimension % step 2 of Velocity-Verlet
        velocity(:,d)=velocity(:,d)+(f(:,d)./m)*(dt*0.5);
    end 
    if step<=Ne; % control temperature in the equilibration stage
        velocity=velocity*sqrt(Temp_k*Dimension*K_B*N/sum(m.*sum(velocity.^2,2))); % scale velocity  should be the velocity under this tempeture.
    elseif mod(step,Ns)==0 % measure in the production stage
        E((step-Ne)/Ns,1)=U/N; % potential (per atom)
        E((step-Ne)/Ns,2)=0.5*sum(m.*sum(velocity.^2,2))/N; % kinetic energy
    end
    if mod(step,8)==0
        plot3(Position(1,1),Position(1,2),Position(1,3),'ro')
        plot3(Position(2,1),Position(2,2),Position(2,3),'go')
        plot3(Position(3,1),Position(3,2),Position(3,3),'bo')
        plot3(Position(4,1),Position(4,2),Position(4,3),'co')
        plot3(Position(5,1),Position(5,2),Position(5,3),'yo')
        plot3(Position(6,1),Position(6,2),Position(6,3),'ko')
    end

end % time-evolution completed
grid on
legend('atom 1','atom 2','atom 3','atom 4','atom 5','atom 6');
str = '   Phase Diagram,T=';
Tstr = num2str(Temp_k);
str = strcat(str,Tstr);
str = strcat(str,',Box side len=');
boxStr = num2str(b);
str = strcat(str,boxStr);
title(str,'fontsize',18);
basepath = pwd;
path_save = strcat(basepath,num2str(time2test));
path_save = strcat(path_save,str);
path_save = strcat(path_save,'.jpg');

saveas(h,path_save);

E(:,3)=E(:,1)+E(:,2); % total enegy (per atom)
