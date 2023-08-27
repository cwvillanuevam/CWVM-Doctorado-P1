function [PQori,PQmod]=DAmarket(Psys,lambda,Pbmin,Pbmax,PEMselect,T)
B=length(lambda);
PQori=zeros(T,1);
PQmod=zeros(T,1);
for t=1:T
    for b=1:B
        if Pbmin(b,1)<=Psys(t,1) && Pbmax(b,1)>Psys(t,1)
            PQori(t,1)=lambda(b,1);
        elseif Pbmax(b,1)<=Psys(t,1)
            PQori(t)=lambda(B,1);
        end
    end
    for b=1:B
        if Pbmin(b,1)<=Psys(t,1)+PEMselect(t,1)/1000 && Pbmax(b,1)>Psys(t,1)+PEMselect(t,1)/1000
            PQmod(t,1)=lambda(b,1);
        elseif Pbmax(b,1)<=Psys(t,1)+PEMselect(t,1)/1000
            PQmod(t,1)=lambda(B,1);
        end
    end
end
end