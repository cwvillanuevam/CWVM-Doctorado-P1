SETs
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
$load v t tini tinivalid t1 t0 tprev s b
$gdxin

scalar Mbig /10000/;
scalar epsilon /0.001/;

parameters
lambda(b) Power cost PQP in b dolar per MW
Pbmin(b) Power lower quantity PQP in b MW
Pbmax(b) Power upper quantity PQP in b MW
Psys(t) System Power in t interval of time MW
R(t,v,s) Road energy consumption for s shortage MCS v ev and t time KWh
X(t,v,s) Road energy binary factor for s shortage MCS v ev and t time bin
SoCinit(v,s) Initial State of Charge for v ev and s shortage MCS KWh
CES(v) Cost of energy storage of v ev dolar per kWh
mES(v) energy storage degradation gradient of v ev pu
pis(t,s) shoratge probability for s MCS scenario
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
C(t0) Cost in t time
CT Total cost for EVA
*CT2 Total cost for EVA adding slack cost
PEM(t0) Power traded in energy market
;

positive variables
PB2G(t0,v,s) Power traded from battery to grid
PB2R(t0,v,s) Power traded from battery to road
PG2B(t0,v,s) Power traded from grid to battery
PEMpos(t0) power traded with energy market positive for buying and negative for selling
PEMposefect(t0,b) PEM efective for b pqp energy market
PEMneg(t0) power traded with energy market positive for buying and negative for selling
PEMnegefect(t0,b) PEM efective for b pqp energy market
Ppos(t0,s) Power exedent sell in real time market
Pneg(t0,s) Power shortage buy from real time market to cover PEM
Pposefect(t0,s,b) Ppos efective for b pqp energy market
Pnegefect(t0,s,b) Pneg efective for b pqp energy market
soc(t,v,s) State of charge
SoCplot(t) State of charge
socdeg(t0,v,s) Degradation of state of charge
hpos(t0,v,s) slack variable
hneg(t0,v,s) slack variable
hsoclo(t0,v,s) slack variable
hlocup(t0,b) slack variable
hposlocup(t0,s,b) slack variable
hneglocup(t0,s,b) slack variable
;

binary variables
pqp(t,b) binary activation of PEM in b pqp in t time
pqppos(t,s,b) binary activation of  Ppos in b pqp  for s MCS scenario t time
pqpneg(t,s,b) binary activation of P neg in b pqp for s MCS scenario t time
;

EQUATIONS
OBJ,timecost(t0),PEMposefectalg1(t0,b),PEMposefectalg2(t0,b),PEMposefectalg3(t0,b),PEMnegefectalg1(t0,b),PEMnegefectalg2(t0,b),PEMnegefectalg3(t0,b),Pposefectalg1(t0,s,b),Pposefectalg2(t0,s,b),Pposefectalg3(t0,s,b),Pnegefectalg1(t0,s,b),Pnegefectalg2(t0,s,b),Pnegefectalg3(t0,s,b),Pposuplim(t0,s),Pneguplim(t0,s),EMposbal(t0,s),EMnegbal(t0,s),SoCeqpos(t0,t1,v,s),SoCeqneg(t0,t1,v,s),Road(t0,v,s),Pvmaxchg(t0,v,s),Pvmaxdsg(t0,v,s),SoCdegeq(t0,t1,v,s),pqplocup(t0,b),pqploclo(t0,b),pqpposlocup(t0,s,b),pqpposloclo(t0,s,b),pqpneglocup(t0,s,b),pqpnegloclo(t0,s,b),PEMdeg(t0),SoCploteq(t),InitSoC(t,tini,v,s),SoClolim(t0,v,s);
*,OBJ2,SoCuplim(t0,v,s)
OBJ.. CT=E=sum(t0,C(t0));
*OBJ2.. CT2=E=CT-sum(t,sum(b,pqp(t,b)+sum(s,pqppos(t,s,b)+pqpneg(t,s,b))));
timecost(t0).. C(t0)=E=0.5*sum(b,(lambda(b)-lambda(b-1))*(PEMposefect(t0,b)-PEMnegefect(t0,b))/1000)-0.5*sum(s,pis(t0,s)*sum(b,(lambda(b)-lambda(b-1))*Pposefect(t0,s,b)/1000))+0.5*sum(s,pis(t0,s)*sum(b,(lambda(b)-lambda(b-1))*Pnegefect(t0,s,b)/1000))+sum(v,mES(v)/100*sum(s,socdeg(t0,v,s)/SScenarios*CES(v)))+sum(v,sum(s,Mbig*(hpos(t0,v,s)+hneg(t0,v,s))));
*+Mbig*sum(b,hlocup(t0,b)+sum(s,hposlocup(t0,s,b)+hneglocup(t0,s,b)))
Pposuplim(t0,s).. Ppos(t0,s)=L=sum(v,PB2G(t0,v,s)*ndsg);
Pneguplim(t0,s).. Pneg(t0,s)=L=sum(v,PB2G(t0,v,s)*ndsg);
EMposbal(t0,s).. PEMpos(t0)-sum(v,PG2B(t0,v,s))=G=0;
EMnegbal(t0,s).. PEMneg(t0)=L=sum(v,PB2G(t0,v,s))*ndsg-Ppos(t0,s)+Pneg(t0,s);
SoCeqpos(t0,t1,v,s)$(tprev(t0,t1)).. soc(t0,v,s)-soc(t1,v,s)-0.5*(PG2B(t0,v,s)*nchg-PB2G(t0,v,s)-PB2R(t0,v,s))+hpos(t0,v,s)=L=0;
SoCeqneg(t0,t1,v,s)$(tprev(t0,t1)).. soc(t0,v,s)-soc(t1,v,s)-0.5*(PG2B(t0,v,s)*nchg-PB2G(t0,v,s)-PB2R(t0,v,s))-hneg(t0,v,s)=G=0;
*SoCuplim(t0,v,s).. soc(t0,v,s)=L=SoCmax;
SoClolim(t0,v,s).. soc(t0,v,s)=G=SoCmin;
*+hsoclo(t0,v,s)
Road(t0,v,s).. 0.5*PB2R(t0,v,s)*ndsg=G=R(t0,v,s);
Pvmaxchg(t0,v,s).. PG2B(t0,v,s)+PB2G(t0,v,s)=L=Pvmax*(1-X(t0,v,s));
Pvmaxdsg(t0,v,s).. PB2R(t0,v,s)=L=Pvmax*X(t0,v,s);
InitSoC(t,tini,v,s)$(tinivalid(t,tini)).. soc(t,v,s)=E=SoCinit(v,s);
SoCdegeq(t0,t1,v,s)$(tprev(t0,t1)).. socdeg(t0,v,s)=G=soc(t1,v,s)-soc(t0,v,s);
*pqplocup(t0,b).. Pbmax(b)*pqp(t0,b)-PEMpos(t0)/1000+PEMneg(t0)/1000-Psys(t0)+Mbig*(1-pqp(t0,b))=G=epsilon;
*pqploclo(t0,b).. Pbmin(b)*pqp(t0,b)=L=PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0);
*pqpposlocup(t0,s,b).. Pbmax(b)*pqppos(t0,s,b)-PEMpos(t0)/1000+PEMneg(t0)/1000-Psys(t0)+Mbig*(1-pqppos(t0,s,b))+Ppos(t0,s)/1000=G=epsilon;
*pqpposloclo(t0,s,b).. Pbmin(b)*pqppos(t0,s,b)=L=PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)-Ppos(t0,s)/1000;
*pqpneglocup(t0,s,b).. Pbmax(b)*pqpneg(t0,s,b)-PEMpos(t0)/1000+PEMneg(t0)/1000-Psys(t0)+Mbig*(1-pqpneg(t0,s,b))-Pneg(t0,s)/1000=G=epsilon;
*pqpnegloclo(t0,s,b).. Pbmin(b)*pqpneg(t0,s,b)=L=PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)+Pneg(t0,s)/1000;
pqplocup(t0,b).. PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)=L=Pbmin(b)+Mbig*pqp(t0,b);
*-hlocup(t0,b)
pqploclo(t0,b).. PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)=G=Pbmin(b)-Mbig*(1-pqp(t0,b));
pqpposlocup(t0,s,b).. PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)-Ppos(t0,s)/1000=L=Pbmin(b)+Mbig*pqppos(t0,s,b);
*-hposlocup(t0,s,b)
pqpposloclo(t0,s,b).. PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)-Ppos(t0,s)/1000=G=Pbmin(b)-Mbig*(1-pqppos(t0,s,b));
pqpneglocup(t0,s,b).. PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)+Pneg(t0,s)/1000=L=Pbmin(b)+Mbig*pqpneg(t0,s,b);
*-hneglocup(t0,s,b)
pqpnegloclo(t0,s,b).. PEMpos(t0)/1000-PEMneg(t0)/1000+Psys(t0)+Pneg(t0,s)/1000=G=Pbmin(b)-Mbig*(1-pqpneg(t0,s,b));
PEMposefectalg1(t0,b).. PEMposefect(t0,b)=L=PEMpos(t0);
PEMposefectalg2(t0,b).. PEMposefect(t0,b)=G=PEMpos(t0)-Mbig*(1-pqp(t0,b));
PEMposefectalg3(t0,b).. PEMposefect(t0,b)=L=Mbig*pqp(t0,b);
PEMnegefectalg1(t0,b).. PEMnegefect(t0,b)=L=PEMneg(t0);
PEMnegefectalg2(t0,b).. PEMnegefect(t0,b)=G=PEMneg(t0)-Mbig*(1-pqp(t0,b));
PEMnegefectalg3(t0,b).. PEMnegefect(t0,b)=L=Mbig*pqp(t0,b);
Pposefectalg1(t0,s,b).. Pposefect(t0,s,b)=L=Ppos(t0,s);
Pposefectalg2(t0,s,b).. Pposefect(t0,s,b)=G=Ppos(t0,s)-Mbig*(1-pqppos(t0,s,b));
Pposefectalg3(t0,s,b).. Pposefect(t0,s,b)=L=Mbig*pqppos(t0,s,b);
Pnegefectalg1(t0,s,b).. Pnegefect(t0,s,b)=L=Pneg(t0,s);
Pnegefectalg2(t0,s,b).. Pnegefect(t0,s,b)=G=Pneg(t0,s)-Mbig*(1-pqpneg(t0,s,b));
Pnegefectalg3(t0,s,b).. Pnegefect(t0,s,b)=L=Mbig*pqpneg(t0,s,b);
PEMdeg(t0).. PEM(t0)=E=PEMpos(t0)-PEMneg(t0);
SoCploteq(t).. SoCplot(t)=E=sum(v,sum(s,soc(t,v,s)));

*Initial values
*pqp.l(t0,b)=1;
*pqppos.l(t0,s,b)=1;
*pqpneg.l(t0,s,b)=1;
*socdeg.l(t,v,s)=0;
*soc.l(t,v,s)=22;
*PEMpos.l(t0)=Psys(t0);
*PB2G.l(t0,v,s)=Pvmax/2;
*PB2G.l(t0,v,s)=Psys(t0)/1000;
*PB2R.l(t0,v,s)=Psys(t0)/1000;
*PG2B.l(t0,v,s)=Psys(t0)/1000;

*Lower values
*soc.lo(t,v,s)=R(t,v,s)*(1-X(t,v,s))/48;
*soc.lo(t0,v,s)=SoCmin;

*Upper values
*PB2G.up(t,v,s)=40;
*PB2R.up(t,v,s)=40;
*PG2B.up(t,v,s)=40;
*PEMpos.up(t0)=10000;
*PEMposefect.up(t0,b)=18000;
*PEMneg.up(t0)=10000;
*PEMnegefect.up(t0,b)=18000;
*Ppos.up(t0,s)=18000;
*Pneg.up(t0,s)=18000;
*Pposefect.up(t0,s,b)=18000;
*Pnegefect.up(t0,s,b)=18000;
soc.up(t,v,s)=SoCmax;

MODEL EVAOPT /ALL/;
option reslim=500;
option optca=0.05;
option optcr=0.05;
option MIP=Cplex;
SOLVE EVAOPT USING MIP MINIMIZING CT;


EXECUTE_UNLOAD 'out.gdx',CT,C,PEM,PB2G,PG2B,PB2R,Ppos,Pneg,SoCplot,hpos,hneg,soc,hsoclo,hlocup,hposlocup,hneglocup;
