function evainvf(TT,SPEED,N,T,V,S,B,lambda,Pbmin,Pbmax,Psys,R,X,...
    SoCinit,CES,mES,Pvmax,nchg,ndsg,BCES,SoCmin,SoCmax,pis)
%% Ingreso de datos

%% SETTING SETS
ttset.name='tt';
ttset.uels={1:TT};

speedset.name='speed';
speedset.uels={1:SPEED};

nset.name='n';
nset.uels={1:N};

vset.name='v';
vset.uels={1:V};

tset.name='t';
tset.uels={1:T};

tiniset.name='tini';
tiniset.uels={1:1};
  
guel = @(s,v) strcat(s,strsplit(num2str(v)));
tinivalidval=[1,1];
tinivalidset.name='tinivalid';
tinivalidset.type='set';
tinivalidset.dim=2;
tinivalidset.form='sparse';
tinivalidset.uels={guel('',1:T),guel('',1:1)};
tinivalidset.val=tinivalidval;


t1=[];
t0=[];
tprevval=[];
for i=2:T
    t1(1,i-1)=i-1;
    t0(1,i-1)=i;
    tprevval(i-1,1)=i;
    tprevval(i-1,2)=i-1;
end
t1set.name='t1';
t1set.uels={t1};
t0set.name='t0';
t0set.uels={t0};
tprevset.name='tprev';
tprevset.type='set';
tprevset.dim=2;
tprevset.form='sparse';
tprevset.uels={guel('',1:T),guel('',1:T)};
tprevset.val=tprevval;

sset.name='s';
sset.uels={1:S};

bset.name='b';
bset.uels={1:B};
%% SETTING PARAMETERS
% deltaups.name='deltaup';
% deltaups.val=deltaupval;
% deltaups.form='full';
% deltaups.type='parameter';
% deltaups.uels=nbusset.uels;

lambdas.name='lambda';
lambdas.val=lambda';
lambdas.form='full';
lambdas.type='parameter';
lambdas.uels=bset.uels;

Pbmins.name='Pbmin';
Pbmins.val=Pbmin';
Pbmins.form='full';
Pbmins.type='parameter';
Pbmins.uels=bset.uels;

Pbmaxs.name='Pbmax';
Pbmaxs.val=Pbmax';
Pbmaxs.form='full';
Pbmaxs.type='parameter';
Pbmaxs.uels=bset.uels;

Psyss.name='Psys';
Psyss.val=Psys';
Psyss.form='full';
Psyss.type='parameter';
Psyss.uels=tset.uels;

Rs.name='R';
Rs.val=R;
Rs.form='full';
Rs.type='parameter';
Rs.uels={ttset.uels speedset.uels nset.uels tset.uels vset.uels sset.uels};

Xs.name='X';
Xs.val=X;
Xs.form='full';
Xs.type='parameter';
Xs.uels={ttset.uels speedset.uels nset.uels tset.uels vset.uels sset.uels};

SoCinits.name='SoCinit';
SoCinits.val=SoCinit;
SoCinits.form='full';
SoCinits.type='parameter';
SoCinits.uels={vset.uels sset.uels};

CESs.name='CES';
CESs.val=CES';
CESs.form='full';
CESs.type='parameter';
CESs.uels=vset.uels;

mESs.name='mES';
mESs.val=mES';
mESs.form='full';
mESs.type='parameter';
mESs.uels=vset.uels;

% Pvmaxs.name='Pvmax';
% Pvmaxs.val=Pvmax;
% Pvmaxs.form='full';
% Pvmaxs.type='scalar';
% 
% nchgs.name='nchg';
% nchgs.val=nchg;
% nchgs.form='full';
% nchgs.type='scalar';
% 
% ndsgs.name='ndsg';
% ndsgs.val=ndsg;
% ndsgs.form='full';
% ndsgs.type='scalar';
% 
% BCESs.name='BCES';
% BCESs.val=BCES;
% BCESs.form='full';
% BCESs.type='scalar';
% 
% SoCmins.name='SoCmin';
% SoCmins.val=SoCmin;
% SoCmins.form='full';
% SoCmins.type='scalar';
% 
% SoCmaxs.name='SoCmax';
% SoCmaxs.val=SoCmax;
% SoCmaxs.form='full';
% SoCmaxs.type='scalar';

piss.name='pis';
piss.val=pis;
piss.form='full';
piss.type='parameter';
piss.uels={ttset.uels speedset.uels nset.uels tset.uels sset.uels};

%% SAVING VALUES
% wgdx('in',nbusset,nbusPVset,nbusPQset,nlinesset,nrset,nlset,nopset,noppset,nbussumset,noleset,nrerset,...
%     ngenset,ngenbusset,vups,vlos,deltaups,deltalos,Pds,Qds,Qshs,Pgups,Pglos,Qgups,Qglos,as,Yms,ts,ags,bgs,cgs,...
%     basemvas,limSs,discls,ntlups,ntllos,cils,abs,Profas,tiops,Vmls,Pgls,ntlls);

wgdx('in.gdx',ttset,speedset,nset,vset,tset,tiniset,tinivalidset,t1set,...
    t0set,tprevset,sset,bset,lambdas,Pbmins,Pbmaxs,...
    Psyss,Rs,Xs,SoCinits,CESs,mESs,piss);%Pvmaxs,nchgs,ndsgs,...    BCESs,SoCmins,SoCmaxs,


end