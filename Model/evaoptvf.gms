SETs
tt
speed
n
v
t
tini(t)
tinivalid(t,t)
t1(t)
t0(t)
tprev(t,t)
s
b
;


$gdxin in
$load tt speed n v t tini tinivalid t1 t0 tprev s b
$gdxin

scalar Mbig /100000000000/;
scalar epsilon /0.0000001/;

parameters
lambda(b) Power cost PQP in b dolar per MW
Pbmin(b) Power lower quantity PQP in b MW
Pbmax(b) Power upper quantity PQP in b MW
Psys(t) System Power in t interval of time MW
R(tt,speed,n,t,v,s) Road energy consumption for s shortage MCS v ev and t time KWh
X(tt,speed,n,t,v,s) Road energy binary factor for s shortage MCS v ev and t time bin
SoCinit(v,s) Initial State of Charge for v ev and s shortage MCS KWh
CES(v) Cost of energy storage of v ev dolar per kWh
mES(v) energy storage degradation gradient of v ev pu
pis(tt,speed,n,t,s) shoratge probability for s MCS scenario
;

scalar
Pvmax /100/
*"Maximum charge rate kW"
nchg /0.9/
*"charging efficiency pu"
ndsg "discharging efficiency pu"   /0.9/
BCES /100/
*"battery capacity of energy storage kWh"
SoCmin /15/
*"minimum State of Charge kWh"
SoCmax /95/
*"maximum State of Charge kWh"
SScenarios /5/
;

$GDXIN in.gdx
$LOAD lambda Pbmin Pbmax Psys R X SoCinit CES mES pis
$GDXIN

Variables
C(tt,speed,n,t0) Cost in t time
CT(tt,speed,n) Total cost for EVA
CT2 Total cost for EVA adding slack cost
PEM(tt,speed,n,t0) Power traded in energy market
;

positive variables
PB2G(tt,speed,n,t0,v,s) Power traded from battery to grid
PB2R(tt,speed,n,t0,v,s) Power traded from battery to road
PG2B(tt,speed,n,t0,v,s) Power traded from grid to battery
PEMpos(tt,speed,n,t0) power traded with energy market positive for buying and negative for selling
PEMposefect(tt,speed,n,t0,b) PEM efective for b pqp energy market
PEMneg(tt,speed,n,t0) power traded with energy market positive for buying and negative for selling
PEMnegefect(tt,speed,n,t0,b) PEM efective for b pqp energy market
Ppos(tt,speed,n,t0,s) Power exedent sell in real time market
Pneg(tt,speed,n,t0,s) Power shortage buy from real time market to cover PEM
Pposefect(tt,speed,n,t0,s,b) Ppos efective for b pqp energy market
Pnegefect(tt,speed,n,t0,s,b) Pneg efective for b pqp energy market
soc(tt,speed,n,t,v,s) State of charge
SoCplot(tt,speed,n,t) State of charge
socdeg(tt,speed,n,t0,v,s) Degradation of state of charge
hpos(tt,speed,n,t0,v,s) slack variable
hneg(tt,speed,n,t0,v,s) slack variable
;

binary variables
pqp(tt,speed,n,t,b) binary activation of PEM in b pqp in t time
pqppos(tt,speed,n,t,s,b) binary activation of  Ppos in b pqp  for s MCS scenario t time
pqpneg(tt,speed,n,t,s,b) binary activation of P neg in b pqp for s MCS scenario t time
;

EQUATIONS
OBJ2,OBJ(tt,speed,n),timecost(tt,speed,n,t0),PEMposefectalg1(tt,speed,n,t0,b),PEMposefectalg2(tt,speed,n,t0,b),PEMposefectalg3(tt,speed,n,t0,b),PEMnegefectalg1(tt,speed,n,t0,b),PEMnegefectalg2(tt,speed,n,t0,b),PEMnegefectalg3(tt,speed,n,t0,b),Pposefectalg1(tt,speed,n,t0,s,b),Pposefectalg2(tt,speed,n,t0,s,b),Pposefectalg3(tt,speed,n,t0,s,b),Pnegefectalg1(tt,speed,n,t0,s,b),Pnegefectalg2(tt,speed,n,t0,s,b),Pnegefectalg3(tt,speed,n,t0,s,b),Pposuplim(tt,speed,n,t0,s),Pneguplim(tt,speed,n,t0,s),EMposbal(tt,speed,n,t0,s),EMnegbal(tt,speed,n,t0,s),SoCeqpos(tt,speed,n,t0,t1,v,s),SoCeqneg(tt,speed,n,t0,t1,v,s),SoCuplim(tt,speed,n,t0,v,s),SoClolim(tt,speed,n,t0,v,s),Road(tt,speed,n,t0,v,s),Pvmaxchg(tt,speed,n,t0,v,s),Pvmaxdsg(tt,speed,n,t0,v,s),InitSoC(tt,speed,n,t,tini,v,s),SoCdegeq(tt,speed,n,t0,t1,v,s),pqplocup(tt,speed,n,t0,b),pqploclo(tt,speed,n,t0,b),pqpposlocup(tt,speed,n,t0,s,b),pqpposloclo(tt,speed,n,t0,s,b),pqpneglocup(tt,speed,n,t0,s,b),pqpnegloclo(tt,speed,n,t0,s,b),PEMdeg(tt,speed,n,t0),SoCploteq(tt,speed,n,t);
OBJ(tt,speed,n).. CT(tt,speed,n)=E=sum(t0,C(tt,speed,n,t0));
OBJ2.. CT2=E=sum(tt,sum(speed,sum(n,CT(tt,speed,n))));
timecost(tt,speed,n,t0).. C(tt,speed,n,t0)=E=0.5*sum(b,lambda(b)*(PEMposefect(tt,speed,n,t0,b)-PEMnegefect(tt,speed,n,t0,b))/1000)-0.5*sum(s,pis(tt,speed,n,t0,s)*sum(b,lambda(b)*Pposefect(tt,speed,n,t0,s,b)/1000))+0.5*sum(s,pis(tt,speed,n,t0,s)*sum(b,lambda(b)*Pnegefect(tt,speed,n,t0,s,b)/1000))+sum(v,mES(v)/100*sum(s,socdeg(tt,speed,n,t0,v,s)/SScenarios*CES(v)))+sum(v,sum(s,Mbig*(hpos(tt,speed,n,t0,v,s)+hneg(tt,speed,n,t0,v,s))));
PEMposefectalg1(tt,speed,n,t0,b).. PEMposefect(tt,speed,n,t0,b)=L=PEMpos(tt,speed,n,t0);
PEMposefectalg2(tt,speed,n,t0,b).. PEMposefect(tt,speed,n,t0,b)=G=PEMpos(tt,speed,n,t0)-Mbig*(1-pqp(tt,speed,n,t0,b));
PEMposefectalg3(tt,speed,n,t0,b).. PEMposefect(tt,speed,n,t0,b)=L=pqp(tt,speed,n,t0,b);
PEMnegefectalg1(tt,speed,n,t0,b).. PEMnegefect(tt,speed,n,t0,b)=L=PEMneg(tt,speed,n,t0);
PEMnegefectalg2(tt,speed,n,t0,b).. PEMnegefect(tt,speed,n,t0,b)=G=PEMneg(tt,speed,n,t0)-Mbig*(1-pqp(tt,speed,n,t0,b));
PEMnegefectalg3(tt,speed,n,t0,b).. PEMnegefect(tt,speed,n,t0,b)=L=pqp(tt,speed,n,t0,b);
Pposefectalg1(tt,speed,n,t0,s,b).. Pposefect(tt,speed,n,t0,s,b)=L=Ppos(tt,speed,n,t0,s);
Pposefectalg2(tt,speed,n,t0,s,b).. Pposefect(tt,speed,n,t0,s,b)=G=Ppos(tt,speed,n,t0,s)-Mbig*(1-pqppos(tt,speed,n,t0,s,b));
Pposefectalg3(tt,speed,n,t0,s,b).. Pposefect(tt,speed,n,t0,s,b)=L=pqppos(tt,speed,n,t0,s,b);
Pnegefectalg1(tt,speed,n,t0,s,b).. Pnegefect(tt,speed,n,t0,s,b)=L=Pneg(tt,speed,n,t0,s);
Pnegefectalg2(tt,speed,n,t0,s,b).. Pnegefect(tt,speed,n,t0,s,b)=G=Pneg(tt,speed,n,t0,s)-Mbig*(1-pqpneg(tt,speed,n,t0,s,b));
Pnegefectalg3(tt,speed,n,t0,s,b).. Pnegefect(tt,speed,n,t0,s,b)=L=pqpneg(tt,speed,n,t0,s,b);
Pposuplim(tt,speed,n,t0,s).. Ppos(tt,speed,n,t0,s)=L=sum(v,PB2G(tt,speed,n,t0,v,s)*ndsg);
Pneguplim(tt,speed,n,t0,s).. Pneg(tt,speed,n,t0,s)=L=sum(v,PB2G(tt,speed,n,t0,v,s)*ndsg);
EMposbal(tt,speed,n,t0,s).. PEMpos(tt,speed,n,t0)-sum(v,PG2B(tt,speed,n,t0,v,s))=G=0;
EMnegbal(tt,speed,n,t0,s).. PEMneg(tt,speed,n,t0)=L=sum(v,PB2G(tt,speed,n,t0,v,s))*ndsg-Ppos(tt,speed,n,t0,s)+Pneg(tt,speed,n,t0,s);
SoCeqpos(tt,speed,n,t0,t1,v,s)$(tprev(t0,t1)).. soc(tt,speed,n,t0,v,s)-soc(tt,speed,n,t1,v,s)-0.5*(PG2B(tt,speed,n,t0,v,s)*nchg-PB2G(tt,speed,n,t0,v,s)-PB2R(tt,speed,n,t0,v,s))+hpos(tt,speed,n,t0,v,s)=L=0;
SoCeqneg(tt,speed,n,t0,t1,v,s)$(tprev(t0,t1)).. soc(tt,speed,n,t0,v,s)-soc(tt,speed,n,t1,v,s)-0.5*(PG2B(tt,speed,n,t0,v,s)*nchg-PB2G(tt,speed,n,t0,v,s)-PB2R(tt,speed,n,t0,v,s))-hneg(tt,speed,n,t0,v,s)=G=0;
SoCuplim(tt,speed,n,t0,v,s).. soc(tt,speed,n,t0,v,s)=L=SoCmax;
SoClolim(tt,speed,n,t0,v,s).. soc(tt,speed,n,t0,v,s)=G=SoCmin;
Road(tt,speed,n,t0,v,s).. 0.5*PB2R(tt,speed,n,t0,v,s)*ndsg=G=R(tt,speed,n,t0,v,s);
Pvmaxchg(tt,speed,n,t0,v,s).. PG2B(tt,speed,n,t0,v,s)+PB2G(tt,speed,n,t0,v,s)=L=Pvmax*(1-X(tt,speed,n,t0,v,s));
Pvmaxdsg(tt,speed,n,t0,v,s).. PB2R(tt,speed,n,t0,v,s)=L=Pvmax*X(tt,speed,n,t0,v,s);
InitSoC(tt,speed,n,t,tini,v,s)$(tinivalid(t,tini)).. soc(tt,speed,n,t,v,s)=G=SoCinit(v,s);
SoCdegeq(tt,speed,n,t0,t1,v,s)$(tprev(t0,t1)).. socdeg(tt,speed,n,t0,v,s)=G=soc(tt,speed,n,t1,v,s)-soc(tt,speed,n,t0,v,s);
pqplocup(tt,speed,n,t0,b).. Pbmax(b)*pqp(tt,speed,n,t0,b)-PEMpos(tt,speed,n,t0)/1000+PEMneg(tt,speed,n,t0)/1000-Psys(t0)=G=epsilon;
pqploclo(tt,speed,n,t0,b).. Pbmin(b)*pqp(tt,speed,n,t0,b)=L=PEMpos(tt,speed,n,t0)/1000-PEMneg(tt,speed,n,t0)/1000+Psys(t0);
pqpposlocup(tt,speed,n,t0,s,b).. Pbmax(b)*pqppos(tt,speed,n,t0,s,b)-PEMpos(tt,speed,n,t0)/1000+PEMneg(tt,speed,n,t0)/1000-Psys(t0)+Ppos(tt,speed,n,t0,s)/1000=G=epsilon;
pqpposloclo(tt,speed,n,t0,s,b).. Pbmin(b)*pqppos(tt,speed,n,t0,s,b)=L=PEMpos(tt,speed,n,t0)/1000-PEMneg(tt,speed,n,t0)/1000+Psys(t0)-Ppos(tt,speed,n,t0,s)/1000;
pqpneglocup(tt,speed,n,t0,s,b).. Pbmax(b)*pqpneg(tt,speed,n,t0,s,b)-PEMpos(tt,speed,n,t0)/1000+PEMneg(tt,speed,n,t0)/1000-Psys(t0)-Pneg(tt,speed,n,t0,s)/1000=G=epsilon;
pqpnegloclo(tt,speed,n,t0,s,b).. Pbmin(b)*pqpneg(tt,speed,n,t0,s,b)=L=PEMpos(tt,speed,n,t0)/1000-PEMneg(tt,speed,n,t0)/1000+Psys(t0)+Pneg(tt,speed,n,t0,s)/1000;
PEMdeg(tt,speed,n,t0).. PEM(tt,speed,n,t0)=E=PEMpos(tt,speed,n,t0)-PEMneg(tt,speed,n,t0);
SoCploteq(tt,speed,n,t).. SoCplot(tt,speed,n,t)=E=sum(v,sum(s,soc(tt,speed,n,t,v,s)));

*Initial values
*socdeg.l(t,v,s)=0;
soc.l(tt,speed,n,t,v,s)=22;
PEMpos.l(tt,speed,n,t0)=Psys(t0);
PB2G.l(tt,speed,n,t0,v,s)=Pvmax/2;
*PB2G.l(t0,v,s)=Psys(t0)/1000;
*PB2R.l(t0,v,s)=Psys(t0)/1000;
*PG2B.l(t0,v,s)=Psys(t0)/1000;

*Lower values
*soc.lo(t,v,s)=R(t,v,s)*(1-X(t,v,s))/48;

*Upper values
*PB2G.up(tt,speed,n,t,v,s)=Pvmax/2;
*PB2R.up(t,v,s)=40;
*PG2B.up(t,v,s)=40;
PEMpos.up(t0)=10000;
*PEMposefect.up(t0,b)=18000;
PEMneg.up(t0)=10000;
*PEMnegefect.up(t0,b)=18000;
*Ppos.up(t0,s)=18000;
*Pneg.up(t0,s)=18000;
*Pposefect.up(t0,s,b)=18000;
*Pnegefect.up(t0,s,b)=18000;

MODEL EVAOPT /ALL/;
option reslim=150;
option optca=0.05;
option optcr=0.05;
option MIP=Cplex;
SOLVE EVAOPT USING MIP MINIMIZING CT2;


EXECUTE_UNLOAD 'out.gdx',CT,C,PEM,PB2G,PG2B,PB2R,Ppos,Pneg,SoCplot,hpos,hneg;
