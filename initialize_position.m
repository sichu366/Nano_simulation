function [Position]=initialize_position(N,Dimension,n0,nxyz,lattice_cons) % FCC crystal

pos0=[0,0,0;0,0.5,0.5;0.5,0,0.5;0.5,0.5,0];
Position=zeros(N,Dimension);
index=0;
for nx=0:nxyz(1)-1
    for ny=0:nxyz(2)-1
        for nz=0:nxyz(3)-1
            for m=1:n0
                index=index+1;
                Position(index,:)=lattice_cons.*([nx,ny,nz]+pos0(m,:));
            end
        end
    end
end
