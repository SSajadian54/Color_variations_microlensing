#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include "VBBinaryLensingLibrary_n.h"
using namespace std;

///**************** constants *********************//
const int YZ=3615;
const int Num=2000;
const double MaxD=20.0;///kpc
const double RA=180.0/M_PI;
const double step=MaxD/(double)Num/1.0;///step in kpc
const int Ntem=int(24979);  
const int M=5;
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
const double pi= M_PI;
const double Pc=3.0857*pow(10.0,16);///[m]
const double G= 6.67408*pow(10.0,-11);/// [m^3/kg.s^2]
const double velocity =299792458.0;///[m/s]
const double Msun=1.989*pow(10.0,30.0);///[kg]
const double Rsun=6.9551*pow(10.0,8);///[m]
const double Au =1.495978707*pow(10.0,11);///[m]
const double binary_fraction=double(2.0/3.0);
const double vro_sun=226.0;
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
const double R_sun=8.0;
const double rho0[8]={4.0,7.9,6.2,4.0,5.8,4.9,6.6,3.96};///considering WD
const double d0[8]= {0.073117,0.0216524,0.0217405,0.0217901,0.0218061,0.0218118,0.0218121,0.0218121};
const double epci[8]={0.014,0.0268,0.0375,0.0551,0.0696,0.0785,0.0791,0.0791};
const double corr[8]={1.0,7.9/4.48419, 6.2/3.52112, 4.0/2.27237, 5.8/3.29525, 4.9/2.78402, 6.6/3.74991, 3.96/2.24994};
const double Rv[4]={3.1,2.5,3.1,3.1};
const int N1=19744,N2=35000, N3=4404, N4=2597;///CMD_BESANCON_new, thinD, bulge, thickD, halo
const double Avks=double(8.20922); 
const double FWHM[M]={0.50, 0.49 ,0.48 , 0.47 , 0.46};//Skottfeld et al. (2013)
const double sigma[M]={0.022, 0.022, 0.02, 0.023 ,0.025};//MOAاستفاده از مقاله کاردلی 
const double AlAv[M]={1.555, 1.332, 1.009, 0.841, 0.600};///From besancon model[UBVI]
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX DANISH Telescope XXXXXXXXXXXXXXXXXXXXXXXXXXXX
const double H_plank= 6.626068*pow(10.,-34.); //in [m^2Kg/s]
const double area=M_PI*pow(1.54/2.0,2.0);//dimeter of telescope is 1.5 meter
const double muskyI= 19.48;//Skotfelt 2015 
const double muskyV= 21.75;
const double mzeroI=fabs(2.5*log10(H_plank*0.790/0.2/(area*2.55*pow(10.0,-23.0))));//Bessell1979,VII
const double mzeroV=fabs(2.5*log10(H_plank*0.551/0.1/(area*3.64*pow(10.0,-23.0))));
const double OmegaV=M_PI*0.48*0.48/4.0; //arcs^2
const double OmegaI=M_PI*0.46*0.46/4.0; ///arcs^2
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


struct lens{
    double RE, tE, Vt, Dl, Ml;
    double t0, ksi, u0, pt1, pt2,tstar;
    double proj,dt;
    double vl, vs, xls;
    double rhomaxl;
    int numl, struc;
    double q, dis; 
    double tcross;
};
struct source{
    int nums, struc, cl; 
    double lon, lat; 
    double mass,Ds, logl,logg;
    double Teff, Rs, ro_star, type;
    double rho_disk[Num],rho_ThD[Num],rho_halo[Num],rho_stars[Num],rho_bulge[Num];
    double Rostar0[Num],Rostari[Num],Nstari[Num];
    double Nstart,Rostart,Romaxs;
    double FI, TET; 
    double Mab[M],Map[M], col;
    double Fluxb[M],magb[M]; 
    double blend[M],nsbl[M]; 
    double ext[M], Lstar[M];
    double slop; 
};
struct CMD{
  double Teff_d[N1],logl_d[N1],Mab_d[M][N1],mass_d[N1], type_d[N1], gra_d[N1],Rs_d[N1], slop_d[N1];int cl_d[N1];//thin_disk
  double Teff_b[N2],logl_b[N2],Mab_b[M][N2],mass_b[N2], type_b[N2], gra_b[N2],Rs_b[N2], slop_b[N2];int cl_b[N2];// bulge
  double Teff_t[N3],logl_t[N3],Mab_t[M][N3],mass_t[N3], type_t[N3], gra_t[N3],Rs_t[N3], slop_t[N3];int cl_t[N3];//thick_disk
  double Teff_h[N4],logl_h[N4],Mab_h[M][N4],mass_h[N4], type_h[N4], gra_h[N4],Rs_h[N4], slop_h[N4];int cl_h[N4];// halo
};
struct extinc{
   double dis[100];///distance 
   double Extks[100];///ks-band extinction
   double Ai[M],Av,Aks,AI;
   double exks;
};
struct Ellip{
    double xo, yo, zo, dx, dy,  tetc, xs, ys, zs, phi, e;
    double Rp, Req, i, c2, geff, geff0, bs, as, tete, tet;
    int flag1, flag2;
    double Vcrit2, omega;  
    double Iie, Iic; 
    double geff_r, geff_t; 
    double R; 
};
struct error{
   double errV, errI, errR, ered, evis;  
   double mvis, mred;
};
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_lens(source & s ,lens & l);
void func_source(source & s, CMD & cm, extinc & ex);
double RandN(double sigma, double N);
double RandR(double down, double up);
void read_cmd(CMD & cm);
void vrel(source & s,lens & l);
void Disk_model(source & s);
void Func_ellip(Ellip  &  p, source & s, lens & l);
double TETA(double x, double y,int flag);
int Extinction(extinc & ex,source & s);
double Interpol(double ds, extinc & ex);
void CauCross(lens & l, int icon); 
void ErrorCal(error & er, double mgR, double mgV, double mgI, double slope, double ratio); 
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    time_t _timeNow;
    unsigned int _randVal;
    unsigned int _dummyVal;
    FILE * _randStream;
///=========================================
int main(int, char **)
{


//================================================
   time(&_timeNow);
   _randStream = fopen("/dev/urandom", "r");
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
   _dummyVal = fread(&_randVal, sizeof(_randVal), 1, _randStream);
   srand(_randVal);
///=================================================
   printf("START_TIME >>>>>>>> %s",ctime(&_timeNow));
   lens l;
   source s;
   CMD cm;
   Ellip p;
   extinc ex;
   error er; 
   read_cmd(cm);



   VBBinaryLensing vbb;
   vbb.Tol=1.e-3;
   vbb.a1=0.0; 
   vbb.LoadESPLTable("./files/ESPL.tbl");
   FILE* magni;
   FILE* danish; 
   FILE* ellips; 
   
   FILE* param;
   param=fopen("./files/MONTEllipB/EllipChroB_all_a.txt","w");
   fclose(param);
   
   FILE * tems;
   tems= fopen("./files/FstarT_Bessel_t.dat","r"); 
   double Tem[Ntem], mm[M][Ntem];
   for(int i=0; i<Ntem; ++i){
   fscanf(tems,"%lf   %lf   %lf   %lf   %lf   %lf\n",&Tem[i],&mm[0][i],&mm[1][i],&mm[2][i],&mm[3][i],&mm[4][i]);}
   fclose(tems); 



    int nx, ny, i_op, imax, jmax, j1, ndat; 
    char filename1[40], filename3[40], filename2[40]; 
    double lonn, Tstar;
    double xlens, ylens, dVIe, dBRe, dVIc, dBRc;
    double Delta, rho0;
    double ux, uy, u;
    double RR, ellip, res[2], rr, circ, A, B,C, mini, t, rstar; 
    double yy, xx, zz, tet0,  Astar1, Astar2, Lstar[M]; 
    double Se0[M], Sc0[M];
    double chi0, chi1, morn, weather, tim, cade; 
    double ffe, col, col0, MV, MI, dchi;
    double ttr, deltaA; 
    int flad;
    double Astara, Astarb, ulens;
    double SqN_V, SqN_I, ratio, MRR, sigc; 
    double interval,  cVI, cBR, deltt; 
    int ftime, flagff;
    double Amax=0.0;   



   s.lat=-4.0;
   lonn=1.0;
   cout<<"latitude: "<<s.lat<<"\t longtitude: "<<lonn<<endl;
   if(lonn<=0.0)   s.lon=360.0+lonn;
   else            s.lon=lonn;
   s.TET=(360.0-s.lon)/RA;///radian
   s.FI=s.lat/RA;///radian
   Disk_model(s);
   if(Extinction(ex,s)==1){
   for(int icon=10; icon<200000; ++icon){ 
   do{
   func_source(s,cm, ex);
   Func_lens(s,l);
   ffe=double((double)rand()/(double)(RAND_MAX+1.));
   }while(l.tE<=0.5 or l.tE>300.0 or ffe>s.blend[4] or s.magb[4]>21.0 or s.magb[4]<12.0 or s.Teff<7000.0 or s.Teff>24999.0); 
   Func_ellip(p, s, l);
   vbb.PrintCau(l.dis,l.q,icon, 4);
   l.dt=double(2.0/60.0/24.0);
   CauCross(l,icon); 
   int Nlens=int( (l.pt2-l.pt1)*1.00/l.dt +3.00 ); 
   cout<<"Nlens:  "<<Nlens<<endl;
  
    
   double **mage=new double*[Nlens];
   double **magc=new double*[Nlens];
   double   **Se=new double*[Nlens];
   double   **Sc=new double*[Nlens];
   for(int i=0; i<Nlens; ++i){
   mage[i] =new double[5];
   magc[i] =new double[5];
   Se[i]   =new double[5];
   Sc[i]   =new double[5]; } 
   
  
   
   rho0= sqrt(p.dx*p.dy/M_PI)*s.Rs*l.proj;
   cout<<"rho0:  "<<rho0<<endl;
   
   for(int k=0;k<M;  ++k){
   for(int i=0; i<Nlens; ++i){
   mage[i][k]=magc[i][k]=0.0;
   Se[i][k]=Sc[i][k]=0.0; 
   Se0[k]=Sc0[k]=0.0;}}



 
    
    double aveL=0.0; 
    for(int i=0;  i<Ntem;  ++i){
    if(( (s.Teff-Tem[i])*(s.Teff-Tem[i-1])<0.0  and  i>0) or s.Teff==Tem[i]){
    for(int j=0; j<M; ++j)  s.Lstar[j]= pow(10.0,-0.4*mm[j][i]);
    break;}}
    if(s.Lstar[0]==0.0  or s.Lstar[1]==0.0  or s.Lstar[2]==0.0 or s.Lstar[3]==0 or s.Lstar[4]==0 or s.Teff<21.0 or s.Teff>24999){  
    cout<<"ERROR Tstar: "<<s.Teff<<"\t Lstar0: "<<s.Lstar[0]<<endl;
    cout<<"Lstar1:  "<<s.Lstar[1]<<"\t Lstar2:  "<<s.Lstar[2]<<endl;   int yye; cin>>yye;}
    aveL= double(s.Lstar[0]+s.Lstar[1]+s.Lstar[2]+s.Lstar[3]+s.Lstar[4])/5.0; 
    for(int i=0; i<M; ++i)    s.Lstar[i]=double(s.Lstar[i]/aveL); 
    
    flagff=0; 
    if(icon<250){
    flagff=1; 
    sprintf(filename1,"./files/MONTEllipB/%c%c%d.dat",'m','_',icon);
    magni=fopen(filename1,"w");
    sprintf(filename2,"./files/MONTEllipB/%c%c%d.dat",'d','_',icon);
    danish=fopen(filename2,"w");
    sprintf(filename3,"./files/MONTEllipB/%c%c%d.dat",'e','_',icon);
    ellips=fopen(filename3,"w"); }
    

  
    xlens = (l.tcross-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
    ylens = (l.tcross-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
    double ucross= sqrt(xlens*xlens+ylens*ylens );
    
    
    cout<<"The parameters: *********************************************"<<endl;
    cout<<"Ml:  "<<l.Ml<<"\t Dl:  "<<l.Dl<<"\t tE:  "<<l.tE<<"\t u0:  "<<l.u0<<endl;
    cout<<"numl:  "<<l.numl<<"\t strucl: "<<l.struc<<"\t Ds:  "<<s.Ds<<"\t l.Vt:  "<<l.Vt<<endl;
    cout<<"Tstar:  "<<s.Teff<<"\t logg:  "<<s.logg<<endl;
    cout<<"Req:  "<<p.Req<<"\t Rp:  "<<p.Rp<<endl;
    cout<<"omega:  "<<p.omega<<"\t incli:  "<<p.i<<"\t p.e:  "<<p.e<<endl;
    cout<<"Lstar[0]:  "<<s.Lstar[0]<<"\t Lstar[1]:  "<<s.Lstar[1]<<endl; 
    cout<<"Lstar[2]:  "<<s.Lstar[2]<<"\t Lstar[3]:  "<<s.Lstar[3]<<"\t Lstar[4]:  "<<s.Lstar[4]<<endl;
    cout<<"ro_star:  "<<s.ro_star<<"\t RE:  "<<l.RE/Au<<endl;
    cout<<"l.dt(min):  "<<l.dt*24.0*60.0 <<"\t cade(min):  "<<cade<<endl; 
    cout<<"dx:  "<<p.dx<<"\t dy:   "<<p.dy<<endl;
    cout<<"**************************************************************"<<endl;
    
    
///===========================================================================================================================
    nx=ny=0;
    for(p.xo=0.0;  p.xo<=p.Req;    p.xo+=p.dx){
    RR=sqrt(p.Req*p.Req-p.xo*p.xo); ++nx; ny=0;
    if(fabs(p.xo-p.Req)<p.dx)      RR=0.0;
    for(p.yo=-RR;  p.yo<=RR;       p.yo +=p.dy){
    p.flag1=p.flag2=0;     res[0]=res[1]=0.0; ++ny;
    ellip = fabs(p.xo*p.xo/p.Req/p.Req + p.yo*p.yo/p.c2/p.c2);
    rr=sqrt(p.xo*p.xo+p.yo*p.yo);
    circ=rr*rr/(p.R*p.R);
    p.tetc=TETA(p.xo,p.yo,0);
    if(circ==1.0 or circ<1.0)  p.flag2=1;
    
    if(ellip<1.0||ellip==1.0){
    A= p.Req*p.Req+p.Rp*p.Rp*cos(p.i)*cos(p.i)/(sin(p.i)*sin(p.i));
    B=-2.0*p.yo*p.Rp*cos(p.i)/(sin(p.i)*sin(p.i));
    C=p.xo*p.xo-p.Req*p.Req+p.yo*p.yo/(sin(p.i)*sin(p.i));
    Delta=B*B-4.0*A*C;
    if(Delta>0.0 || Delta==0.0 || fabs(Delta)<0.000000001){
    p.flag1=1; mini=1000000.0; i_op=-1;
    if(fabs(Delta)<0.000000001) Delta=0.0;
    res[0]=(-B+sqrt(Delta))/(2.0*A);
    res[1]=(-B-sqrt(Delta))/(2.0*A);
    if(Delta==0.0) imax=1;
    else imax=2;
    for(int i=0; i<imax; ++i){
    p.zo=p.Req*p.Req*(1.0-res[i]*res[i])+ p.Rp*p.Rp*res[i]*res[i]-p.xo*p.xo-p.yo*p.yo;
    if(p.zo>0.0 || p.zo==0.0 || (fabs(p.zo)<0.0000000001 && p.zo<0.0)){
    if(fabs(p.zo)<0.0000000001 && p.zo<0.0)  p.zo=0.0;
    if(p.zo==0.0) jmax=1;
    else jmax=2;
    for(int j=0; j<jmax;  ++j){
    if(j==0) p.zo= 1.0*sqrt(p.zo);
    if(j==1) p.zo=-1.0*sqrt(p.zo);
    p.xs= p.xo;
    p.ys= cos(p.i)*p.yo+sin(p.i)*p.zo;
    p.zs=-sin(p.i)*p.yo+cos(p.i)*p.zo;
    p.phi=TETA(p.zs,p.xs,0);
    if(p.phi==0.0 || p.phi==pi || (p.xs==0.0 && p.ys==0.0))                tet0=TETA(p.ys/p.Rp, p.zs/(p.Req*cos(p.phi)),0);
    else if(p.phi==pi/2.0 || p.phi==3.0*pi/2.0 ||(p.zs==0.0 && p.ys==0.0)) tet0=TETA(p.ys/p.Rp,p.xs/(p.Req*sin(p.phi)),0);
    else{
    yy=p.ys/p.Rp;   zz=p.zs/(p.Req*cos(p.phi));   tet0= TETA(yy,zz,0);}
    if(fabs(res[i]-cos(tet0))<1e-10){  p.tet=tet0;  mini=0.0;  j=jmax+2;  i=imax+2; }
    else if(fabs(res[i]-cos(tet0))<mini){mini=fabs(res[i]-cos(tet0));  p.tet=tet0; i_op=i;}}}}
    if(mini>0.01) {p.tet= acos(res[i_op]);}
    
    rstar=sqrt(p.Req*p.Req*sin(p.tet)*sin(p.tet)+p.Rp*p.Rp*cos(p.tet)*cos(p.tet))*s.Rs;
    p.geff_r= -s.mass/(rstar*rstar)+ p.omega*p.omega*p.Vcrit2*rstar*sin(p.tet)*sin(p.tet);
    p.geff_t= p.omega*p.omega*p.Vcrit2*rstar*sin(p.tet)*cos(p.tet);
    p.geff=sqrt(p.geff_r*p.geff_r+p.geff_t*p.geff_t)/p.geff0;
    
    p.bs=sqrt(p.yo*p.yo+p.xo*p.xo*p.c2*p.c2/(p.Req*p.Req));
    p.as=p.Req*p.bs/p.c2;
    p.tete=TETA(p.xo/p.as/p.as,p.yo/p.bs/p.bs,0);///normal to ellipce at the point (xo,yo)
    if(fabs(p.tetc-p.tete)*RA>20.0){
    cout<<"tetc: "<<p.tetc*RA<<"\t tete: "<<p.tete*RA<<"\t xo: "<<p.xo<<"\t yo: "<<p.yo<<endl;
    cout<<"big error (stopped)"<<endl;  int yyye; cin>>yyye;}}///end of if B^2-4AC
    else if((B*B-4.0*A*C)<0.0){
    cout<<"BIG ERROR in ellipce B^2-4AC is negative (Stopped): "<<B*B-4.0*A*C<<endl;
    cout<<"xo: "<<p.xo<<"\t yo: "<<p.yo<<"\t rr: "<<rr<<"\t ellip:"<<ellip<<endl; int eew; cin>>eew;}}///end of if ellipc
    
    if(p.flag1==1 or p.flag2==1){
    Tstar=s.Teff*pow(p.geff,0.25); 
    j1=-1; 
    if(Tstar<Tem[0] or Tstar==Tem[0])  j1=0;   
    else if(Tstar>Tem[Ntem-1])         j1=int(Ntem-1); 
    else{
    for(int k=1; k<Ntem; ++k){
    if( (Tstar-Tem[k])*(Tstar-Tem[k-1])<0.0 or Tstar==Tem[k]){ j1=k;   break;}}}
    
    for(int k=0; k<M; ++k){ Lstar[k]=pow(10.0,-0.4*mm[k][j1])/aveL;} 
    if(j1<0 or Tstar<21 or Tstar>25000.0){ cout<<"j1: "<<j1<<"\t Tstar:  "<<Tstar<<endl;  int uue;  cin>>uue;}
    
    
    for(int i=0; i<Nlens; ++i){
    t=double(l.pt1 + i*l.dt);
    xlens = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
    ylens = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi);
    ulens=sqrt(xlens*xlens + ylens*ylens);
    uy= ylens-p.yo*s.Rs*l.proj;
    
    
    if(double(fabs(ulens-ucross)/s.ro_star)<2.0)  flad=1;
    else flad=0;
    ux= xlens-p.xo*s.Rs*l.proj;
    if(flad>0)   Astar1=vbb.BinaryMag2(l.dis,l.q, ux, uy,rho0);
    else         Astar1=vbb.BinaryMag0(l.dis,l.q, ux, uy);    
    ux= xlens +p.xo*s.Rs*l.proj;
    if(flad>0)    Astar2=vbb.BinaryMag2(l.dis,l.q, ux, uy,rho0);
    else          Astar2=vbb.BinaryMag0(l.dis,l.q, ux, uy);    

    
    for(int k=0; k<M; ++k){
    if(i==0) Se0[k] +=2.0*  Lstar[k]*p.flag1;///ellipce
    if(i==0) Sc0[k] +=2.0*s.Lstar[k]*p.flag2;///cirlce
    Se[i][k] += (Astar1+Astar2)*  Lstar[k]*p.flag1;
    Sc[i][k] += (Astar1+Astar2)*s.Lstar[k]*p.flag2;}}
    
    
    if(p.flag1==1 and flagff>0){
    fprintf(ellips,"%-12.5lf  %-12.5lf   %-12.8lf   %-12.8lf  %-12.8lf  %.5lf\n",+1.0*p.xo,p.yo,p.geff,Lstar[2],Lstar[4],Tstar);
    fprintf(ellips,"%-12.5lf  %-12.5lf   %-12.8lf   %-12.8lf  %-12.8lf  %.5lf\n",-1.0*p.xo,p.yo,p.geff,Lstar[2],Lstar[4],Tstar);}}
    } cout<<"xo:  "<<p.xo<<"\t yo:  "<<p.yo<<endl;
    cout<<"*****************************************"<<endl;}//end of xo, yo loo
    for(int i=0; i<Nlens; i++){
    for(int k=0;k<M; ++k){
    mage[i][k]= s.blend[k] * double(Se[i][k]/Se0[k]) + 1.0 - s.blend[k];
    magc[i][k]= s.blend[k] * double(Sc[i][k]/Sc0[k]) + 1.0 - s.blend[k]; } }
///===========================================================================================================================


    tim=0.0;
    chi0=chi1=0.0;  
    ndat=0; 
    morn=fabs((double)rand()*24.0/(double)(RAND_MAX+1.)); 
    MV=MI=col=col0=0.0;
    
    interval=cVI=cBR=deltt=0.0; 
    ftime=0;   Amax=0.0;
    for(int k=0; k<Nlens; ++k){
    t=double(l.pt1 + k*l.dt);
    xlens = (t-l.t0)/l.tE * cos(l.ksi) - l.u0 * sin(l.ksi);
    ylens = (t-l.t0)/l.tE * sin(l.ksi) + l.u0 * cos(l.ksi); 
    cade= 40.0*sqrt(3.0/sqrt(5.0)/fabs(mage[k][4]))/60.0;///hour    
    if(cade<  double(2.0/60.0)) cade=double(2.0/60.0); ///hour
    if(cade>double(120.0/60.0)) cade=double(120.0/60.0);///hour 
    tim += double(l.dt*24.0);///hour 
    morn+= double(l.dt*24.0);///hour
    if(morn>24.0) morn -=24.0; 
    MV= s.magb[2]-2.5*log10( mage[k][2] );
    MI= s.magb[4]-2.5*log10( mage[k][4] ); 
    MRR=s.magb[3]-2.5*log10( mage[k][3] );  
    weather=fabs((double)rand()*100.0/(double)(RAND_MAX+1.));
    if(mage[k][4]>Amax)   Amax=fabs(mage[k][4]);
    if(tim>cade and weather<80.0 and MI>=11.0 and morn<6.5){    
    col0=double(MV-MI);
    ttr= RandN(1.0,3.0); 
    SqN_V=sqrt( fabs(pow(10.0,-0.4*MV) + pow(10.0,-0.4*muskyV)*OmegaV )*pow(10.0,0.4*mzeroV)*1.0 );  
    SqN_I=sqrt( fabs(pow(10.0,-0.4*MI) + pow(10.0,-0.4*muskyI)*OmegaI )*pow(10.0,0.4*mzeroI)*1.0 );
    ratio=fabs( double(1.0/SqN_V)/double(1.0/SqN_I) ); //(DF/F)_V / (DF/F)_I
    ErrorCal(er, MRR, MV, MI, s.slop, ratio); 
    sigc= sqrt(er.errI*er.errI + er.errV* er.errV);
    col= col0+ ttr* sigc;
    deltaA= fabs(pow(10.0,-0.4*er.errI)-1.0)*mage[k][4]; ///error in magnification
    chi0 += (col- col0 )*(col- col0)/(sigc*sigc);
    chi1 += (col- s.col)*(col-s.col)/(sigc*sigc); 
    if(flagff>0) fprintf(danish,"%.5lf %.5lf  %.5lf  %.5lf  %.5lf\n",
    (t-l.t0)/l.tstar,double(col-s.col)*1000.0,double(mage[k][4]+ttr*deltaA),sigc*1000.0,deltaA);
    ndat+=1;}
    if(tim>cade) tim= tim-cade;


    dVIe= -2.5*log10(double(mage[k][2])/double(mage[k][4]))*1000.0;
    dBRe= -2.5*log10(double(mage[k][1])/double(mage[k][3]))*1000.0;
    dVIc= -2.5*log10(double(magc[k][2])/double(magc[k][4]))*1000.0;
    dBRc= -2.5*log10(double(magc[k][1])/double(magc[k][3]))*1000.0;
    if(flagff>0) fprintf(magni,"%.4lf  %.9lf   %.9lf  %.5lf   %.5lf  %.5lf   %.5lf  %.5lf %.5lf %.5lf %.5lf %.5lf %.5lf %e  %e  %e  %e\n", (t-l.t0)/l.tstar,xlens,ylens, mage[k][0],mage[k][1],mage[k][2],mage[k][3],mage[k][4],
   magc[k][0],magc[k][1],magc[k][2],magc[k][3],magc[k][4],dVIe,dBRe, dVIc, dBRc);//17
   
   if(fabs(dVIe)>1.0 and ftime==0) {interval=t;    ftime=1; }
   if((fabs(dVIe)<1.0 or k==int(Nlens-1)) and ftime==1){deltt+=double(t-interval);  ftime=0;}
   if(cVI<fabs(dVIe))  cVI=fabs(dVIe); 
   if(cBR<fabs(dBRe))  cBR=fabs(dBRe); 
    
   
   
   if(k%500==0){
   cout<<"t:  "<<t<<"\t xlens:  "<<xlens<<"\t ylens:  "<<ylens<<endl;
   cout<<"delVI_cirlce:  "<<dVIc<<"\t delBR_circle:  "<<dBRc<<endl;
   cout<<"delVI_ellip :  "<<dVIe<<"\t delBR_ellip :  "<<dBRe<<endl;
   cout<<"k:  "<<k<<"\t Nlens:  "<<Nlens<<endl;
   cout<<"magt[0]: "<<mage[k][0]<<"\t magt[1]:  "<<mage[k][1]<<"\t magt[2]:  "<<mage[k][2]<<endl;
   cout<<"magt[3]: "<<mage[k][3]<<"\t magt[4]:  "<<mage[k][4]<<endl;
   cout<<"*******************************************************"<<endl;}}
   
   if(flagff>0) { fclose(magni);   fclose(danish);  fclose(ellips); }
   dchi=fabs(chi0-chi1); 
   param=fopen("./files/MONTEllipB/EllipChroB_all_a.txt","a+");
   fprintf(param,"%d  %d  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.7lf   %.5lf  %.5lf  %d  %d  %d  %.5lf  %.5lf  %.5lf %.5lf %.5lf  %.5lf  %.8lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf   %.5lf  %.6lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %d   %.5lf  %d %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf  %.5lf   %.5lf \n",
   l.numl,l.struc,l.Ml,l.Dl,l.xls,l.RE/Au,l.Vt,l.vs,l.vl,l.tE,l.u0/s.ro_star,l.ksi*RA,l.tstar,//13
   s.nums,s.struc,s.cl,s.mass,s.Ds,s.logg,s.logl,s.Teff,s.Rs,s.ro_star,s.magb[4],s.Map[4],//25
   p.i*RA,p.Req,p.Rp,p.omega,p.e,s.Lstar[0], s.Lstar[1], s.Lstar[2], s.Lstar[3], s.Lstar[4],icon, dchi, ndat, s.blend[4], l.q, l.dis
   , s.slop,deltt,cVI, cBR, Amax);//46
   fclose(param);
    cout<<"******************************************************"<<endl;
    cout<<" counter:    "<<icon<<endl;
    cout<<"LENS:   mass[Msun]: "<<l.Ml<<"\t lens_dis[kpc]:  "<<l.Dl<<"\t RE[AU]: "<<l.RE/Au<<endl;
    cout<<"xls:   "<<l.xls<<"\t u0:  "<<l.u0<<"\t ksi:  "<<l.ksi<<endl;
    cout<<"tstar:  "<<l.tstar<<"\t pt1:  "<<l.pt1<<"\t pt2:   "<<l.pt2<<endl;
    cout<<"tE[days]: "<<l.tE<<"\t Vt:  "<<l.Vt<<"\t Vl:  "<<l.vl<<endl;
    cout<<"source:  mass:  "<<s.mass<<"Ds[Kpc]: "<<s.Ds<<endl;
    cout<<"magb[4]:  "<<s.magb[4]<<"\t Map[4]: "<<s.Map[4]<<"\t blend[4]:  "<<s.blend[4]<<endl;
    cout<<"Rstar[Rsun]: "<<s.Rs<<"\t Tstar[K]:  "<<s.Teff<<endl;
    cout<<"ro_star: "<<s.ro_star<<"\t u0: "<<l.u0/s.ro_star<<endl;
    cout<<"logg:  "<<s.logg<<"\t logl:  "<<s.logl<<"\t Teff:  "<<s.Teff<<endl;
    cout<<"incli(deg):  "<<p.i*RA<<"\t Req:  "<<p.Req<<"\t Rp:  "<<p.Rp<<endl;
    cout<<"omega:  "<<p.omega<<"\t Acceleration_N:  "<<p.Vcrit2*p.Req*s.Rs<<endl;
    cout<<"Lstar[0]:  "<<s.Lstar[0]<<"\t Lstar[1]:  "<<s.Lstar[1]<<"\t Lstar[2]:  "<<s.Lstar[2]<<endl;
    cout<<"Lstar[3]: "<<s.Lstar[3]<<"\t Lstar[4]:  "<<s.Lstar[4]<<endl; 
    cout<<"dchi:  "<<dchi<<"\t ndat:  "<<ndat<<endl;
    cout<<"***************************************************** "<<endl;
    
    
    for(int i = 0; i <Nlens; ++i){
    delete [] mage[i], magc[i], Se[i], Sc[i]; }
    delete [] mage, magc, Se, Sc;  
    
    
    }}//end of loop


    time(&_timeNow);
    printf("END_TIME >>>>>>>>  %s ",ctime(&_timeNow));
    return(0);
}

///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
void CauCross(lens & l, int icon){

   FILE *caustic;
   char filnam1[40];
   sprintf(filnam1,"./files/MONTEllipB/%c%c%c%c%d.dat",'c','a','u','_',icon);
   caustic=fopen(filnam1,"r");
   if(!caustic) {cout<<"File is not exit!!!! "<<"\t icoc:  "<<icon<<endl;}
   double xcau[800], ycau[800];
   for (int i=0; i<800; ++i) {
   fscanf(caustic, "%lf    %lf\n", &xcau[i],&ycau[i]);}
   fclose(caustic); 
   cout<<"caustic is read!!!  "<<endl;
 
   int rad; 
   do{ 
   rad= int(fabs((double)rand()/(double)(RAND_MAX+1.)*798.0));
   if(rad<0 or rad>799) {cout<<"rad:   "<<rad<<endl;   int rrw;  cin>>rrw; } 
   l.u0= sqrt( fabs(xcau[rad]*xcau[rad] + ycau[rad]*ycau[rad])); 
   }while(l.u0>0.5);
   l.ksi= atan2(double(-xcau[rad]),  double(ycau[rad])); 
   l.pt1=-1.1;
   l.pt2=+1.1;         
   cout<<"pt1:  "<<l.pt1<<"\t pt2:  "<<l.pt2<<endl;
      
}
///&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
void ErrorCal(error & er, double mgR, double mgV, double mgI, double slope, double  ratio){

   double a[3]={0.9163,  0.0837,  1.6673};
   double b[3]={-0.0169, 1.0169,  1.7526}; 
   double dVI, sume, mind, tt; 
   int ie; 

   if(mgR<12.0)  er.errR= -3.5;
   if(mgR<18.5) {er.errR=(0.2561987)*mgR+ double(-6.2816372);}
   else if(mgR>=18.5){er.errR=(0.3334653)*mgR+ double(-7.6907161);}
   er.errR= pow(10.0,er.errR); 
   if((er.errR<0.0006 and mgR>=12.0) or (er.errR>0.5 and mgR<22.0)) {
   cout<<"Error errR:  "<<er.errR<<"\t magR:  "<<mgR<<endl;  int yye;  cin>>yye;}
   dVI= (slope*er.errR);

    
   mind=100000000000.0; ie=-1; 
   for (int i=0; i<300; i++){
   er.errV=double(dVI*i/300.0);
   er.errI=sqrt(dVI*dVI-er.errV*er.errV);
   tt=fabs(er.errV - fabs(2.5*log10(fabs(1.0+ ratio*fabs(pow(10.0,0.4*er.errI)-1.0)))));
   if(tt<mind){mind=tt;  ie=i;}}
   er.errV=double(dVI*ie/300.0); 
   er.errI=sqrt(dVI*dVI-er.errV*er.errV);
   if( (ratio>1.0 and er.errV<er.errI) or er.errV<0.0  or er.errI<0.0 or  (ratio<1.0 and er.errI<er.errV)){
   cout<<"Error (2):  alfa: "<<ratio<<"\t Sum:  "<<endl;
   cout<<"errV:  "<<er.errV<<"\t errI:    "<<er.errI<<endl;
   cout<<"dVI:   "<<dVI<<"\t slope:    "<<slope<<"\t   errR:   "<<er.errR<<endl;
   int uuw; cin>>uuw; }
   
  
   er.mvis= double(b[1]*mgV- a[1]*mgI-(a[2]*b[1]- a[1]*b[2]))/(a[0]*b[1]-a[1]*b[0]); 
   er.mred= double(b[0]*mgV- a[0]*mgI-(a[2]*b[0]- a[0]*b[2]))/(a[1]*b[0]-a[0]*b[1]);
   er.evis=sqrt(pow(b[1]*er.errV,2.0)+pow(a[1]*er.errI,2.0))/fabs(a[0]*b[1]-a[1]*b[0]);
   er.ered=sqrt(pow(b[0]*er.errV,2.0)+pow(a[0]*er.errI,2.0))/fabs(a[1]*b[0]-a[0]*b[1]);  

   
   if(fabs(dVI-sqrt(er.errV*er.errV +er.errI*er.errI))>0.5){
   cout<<"ERRRR!!!!!   errR:   "<<er.errR<<"\t errV:  "<<er.errV<<"\t errI:   "<<er.errI<<endl;
   cout<<"R: "<<er.errR<<"\t sqrt(V2+I2)/slope:   "<<sqrt(er.errV*er.errV + er.errI*er.errI)/slope<<endl; 
   cout<<"evis:  "<<er.evis<<"\t ered:  "<<er.ered<<"\t magR:  "<<mgR<<"\t slope:  "<<slope<<endl;
   int uue; cin>> uue; }
   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_ellip                              ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_ellip(Ellip & p, source & s, lens & l)
{
    p.omega= double(((double)rand()/(double)(RAND_MAX+1.0))*0.7+0.2);    
    p.i=(double( (double)rand()/(double)(RAND_MAX+1.0) )*80.0+5.0)/RA;    
    p.Rp=1.0  ;/// Rp=c normalized to the Rstar
    p.Req=3.0*p.Rp*cos((pi+acos(p.omega))/3.0)/p.omega;  //// Req=a
    p.dx=p.dy=double(p.Rp/35.54830547);
    p.c2= sqrt(cos(p.i)*cos(p.i)*p.Rp*p.Rp+p.Req*p.Req*sin(p.i)*sin(p.i));
    p.Vcrit2=pow(2.0/(3.0*p.Rp*s.Rs),3.0)*s.mass;
    p.e=sqrt(1.0-p.c2*p.c2/(p.Req*p.Req));
    cout<<"Req: "<<p.Req<<"\t Vcirt: "<<p.Vcrit2<<"\t Req_projected:  "<<p.c2<<endl;
    cout<<"inclination:  "<<p.i*RA<<"\t omega:  "<<p.omega<<endl;
    p.geff0=fabs(s.mass/(p.Rp*p.Rp*s.Rs*s.Rs));///geff at tet=0
    p.R=p.Rp;//sqrt(p.c2*p.Req);   
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_source                             ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void func_source(source & s, CMD & cm, extinc & ex)
{
     int num,struc,nums;
    double rho,rf,Nblend[M];
    double Akv,Avk,Alv;
    double Ds,Ai[M],Av;
    double Map[M];
   

    double maxnb=0.0;
    for(int i=0; i<M; ++i){
    s.Fluxb[i]=Nblend[i]=0.0;
    Nblend[i]=s.Nstart*pow(FWHM[i]*0.5,2)*M_PI/(3600.0*3600.0);
    //Nblend[i]=Nblend[i]+RandN(sqrt(Nblend[i]),1.0);
    if(Nblend[i]<0.0)  Nblend[i]=0.0; 
    if(Nblend[i]<=1.0) Nblend[i]=1.0;
    s.nsbl[i]=  double(Nblend[i]);
    if(Nblend[i]>maxnb) maxnb=Nblend[i];  }



    for(int k=1; k<=int(maxnb); ++k){
    do{
    num=int(fabs((double)rand()/(double)(RAND_MAX+1.)*Num*1.0));
    rho=fabs((double)rand()/((double)(RAND_MAX+1.))*s.Romaxs);
    }while(rho>s.Rostari[num] || num<5);///distance larger than 50.0 
    Ds=(double)(num*step);///in kpc
    nums=num;
    if(k==1){s.Ds=Ds;  s.nums=nums;}




    rf=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[nums];
         if (rf<= s.rho_disk[nums]) struc=0;///thin disk
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums])) struc=1;/// bulge structure
    else if (rf<=(s.rho_disk[nums]+s.rho_bulge[nums]+s.rho_ThD[nums])) struc=2;///thick disk
    else if (rf<=s.Rostar0[nums]) struc=3;///halo
    if(k==1)    s.struc=struc;
    //cout<<"Ds:  "<<Ds<<"\t struc:  "<<struc<<endl;


   double Mab[M]; 
    if(struc==0){///thin disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N1-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_d[i][num];}
    if(k==1){ 
    s.type=cm.type_d[num];
    s.mass=cm.mass_d[num];
    s.Teff=cm.Teff_d[num]; 
    s.logl=cm.logl_d[num];
    s.cl=    cm.cl_d[num]; 
    s.logg= cm.gra_d[num]; 
    s.Rs=    cm.Rs_d[num];
    s.slop=cm.slop_d[num];}} 


    if(struc==1){///bulge
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N2-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_b[i][num];}
    if(k==1){ 
    s.type=cm.type_b[num];
    s.mass=cm.mass_b[num];
    s.Teff=cm.Teff_b[num];
    s.logl=cm.logl_b[num];
    s.cl=    cm.cl_b[num];
    s.logg =cm.gra_b[num]; 
    s.Rs =   cm.Rs_b[num];
    s.slop=cm.slop_b[num];}}


    if(struc==2){///thick disk
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N3-2.0));
    for(int i=0; i<M; ++i){Mab[i]=cm.Mab_t[i][num]; }
    if(k==1){ 
    s.type=cm.type_t[num];
    s.mass=cm.mass_t[num];
    s.Teff=cm.Teff_t[num];
    s.logl=cm.logl_t[num];
    s.cl=    cm.cl_t[num];
    s.logg =cm.gra_t[num];
    s.Rs =   cm.Rs_t[num];
    s.slop=cm.slop_t[num];}} 


    if(struc==3){///stellar halo
    num=int(((double)rand()/(double)(RAND_MAX+1.))*(N4-2.0));
    for(int i=0; i<M; ++i) { Mab[i]=cm.Mab_h[i][num];}
    if(k==1){
    s.cl  =  cm.cl_h[num]; 
    s.type=cm.type_h[num];
    s.mass=cm.mass_h[num];
    s.Teff=cm.Teff_h[num];
    s.logl=cm.logl_h[num];
    s.logg =cm.gra_h[num]; 
    s.Rs =   cm.Rs_h[num];
    s.slop=cm.slop_h[num];}} 
    
    if(s.type>8.0 or s.type<2.0 or s.Rs<0.0 or s.Rs>1000.0 or s.mass<0.0 or s.Teff<0.0 or s.cl==6 or s.logg<0.0 or s.cl<0  
    or s.cl>9 or  s.slop<1.5 or s.slop>2.5){
    cout<<"ERROR:  type: "<<s.type<<"\t struc: "<<struc<<"\t num: "<<num<<endl; int rre;  cin>>rre;}
    
   
    ex.Aks=Interpol(Ds,ex);///extinction in Ks-band
    Av=ex.Aks*Avks;
    if(Av>20.0 or Av<0.0 or Ds>20.0  or Ds<0.0){cout<<"ERROR Ds:  "<<Ds<<" \t Av:  "<<Av<<endl; int yyw;  cin>>yyw; }
    if(Av<0.0)   Av=0.0;


    for(int i=0; i<M; ++i){    
    Ai[i]=fabs(Av*AlAv[i])+RandN(sigma[i],1.0);
    if(Ai[i]<0.0) Ai[i]=0.0;
    Map[i]=Mab[i]+5.0*log10(Ds*100.0)+Ai[i];
    if(Nblend[i]>=k)  s.Fluxb[i]+=fabs(pow(10.0,-0.4*Map[i]));
    if(k==1){s.ext[i]=Ai[i];  s.Map[i]=Map[i];    s.Mab[i]= Mab[i];}
    if(Ai[i]<0.0  or Ai[i]>100.0 or Map[i]<Mab[i] or Ds<0.0  or Ds>20.0  or Map[i]<0.0 or s.Fluxb[i]<0.0){
    cout<<"ERROR filter:  "<<i<<"\t extinction:  "<<Ai[i]<<"\t App_mag:  "<<Map[i]<<"\t Abso_mag:  "<<Mab[i]<<endl;
    int rre; cin>>rre;}}
    }///loop 


    for(int i=0; i<M; ++i){
    if(s.Fluxb[i]<=0.0){ cout<<"BIG ERROR Flux is negative: "<<s.Fluxb[i]<<endl; int yye; cin>>yye; }
    s.magb[i] = -2.5*log10( s.Fluxb[i] );
    s.col=s.magb[2]-s.magb[4];//V-I
    s.blend[i]=double(pow(10.0,-0.4*s.Map[i])/s.Fluxb[i]);
    if(int(s.nsbl[i])<0 || s.nsbl[i]==0.0 || (Nblend[i]==1.0 && s.blend[i]<1.0) or s.blend[i]>1.0  or s.blend[i]<0.0){ 
    cout<<"BIGG ERRROR nsbl: "<<s.nsbl[i]<<"\t Nlend: "<<Nblend[i]<<"\t s.blend  "<<s.blend[i]<<endl; int uue; cin>>uue;}}
}
////XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                        Func_lens                               ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Func_lens(source & s , lens & l)
{
    double rholens[Num],test, f;
    l.rhomaxl=0.0;
      
    for(int k=1;k<=s.nums;++k){
    rholens[k]=0.0;
    l.Dl=k*step;
    l.xls=l.Dl/s.Ds;
    if(l.Dl>s.Ds){cout<<"ERROR (Dl>Ds) Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;  int yye; cin>>yye;}
    rholens[k]= sqrt((s.Ds-l.Dl)*l.Dl/s.Ds)*s.Rostar0[k];
    if(rholens[k]>l.rhomaxl) l.rhomaxl=rholens[k];}



    do{
    l.numl = (int)((double)rand()*1000./((double)(RAND_MAX+1.)*1000.)*(s.nums-2.0)+1.0);
    test = ((double)rand()/(double)(RAND_MAX+1.)*l.rhomaxl);
    if(rholens[l.numl]>l.rhomaxl){cout<<"ERROR: rholens[numl]: "<<rholens[l.numl]<<""<<l.rhomaxl<<endl;
    int ue; cin>>ue;}
    }while(test>rholens[l.numl]);



   double  randflag=((double)rand()/(double)(RAND_MAX+1.))*s.Rostar0[l.numl];
        if (randflag<=s.rho_disk[l.numl]) l.struc=0;///thin disk
   else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl])) l.struc=1; // bulge structure
   else if (randflag<=(s.rho_disk[l.numl]+s.rho_bulge[l.numl]+s.rho_ThD[l.numl])) l.struc=2; //thick disk
   else if (randflag<= s.Rostar0[l.numl]) l.struc=3;//halo
   else {cout<<"randflag: "<<randflag<<"\t rho_star0: "<<s.Rostar0[l.numl]<<endl;  int ye; cin>>ye;}



  if(l.struc==0){///thin disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(4.5-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*57.0);
  if(l.Ml<=1.0) f=pow(l.Ml,-1.6);
  if(l.Ml>1.0)  f=pow(l.Ml,-3.0);
  }while(test>f);}


  if(l.struc==1){///Galactic bulge
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(1.4-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*378.0)+0.3;
  f=pow(l.Ml,-2.35);
  }while(test>f);}


  if(l.struc==2){///thick disk
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(1.4-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*3.0)+0.8;
  f=pow(l.Ml,-0.5);
  }while(test>f);}


  if(l.struc==3){///stellar halo
  do{
  l.Ml=fabs((double)rand()/(double)(RAND_MAX+1.)*(0.8-0.08))+0.08;
  test=fabs((double)rand()/(double)(RAND_MAX+1.)*1.0)+0.8;
  f=pow(l.Ml,-0.5);
  }while(test>f);}


    l.q=   RandR(0.1,1.0);
    l.dis= RandR(0.4,1.4);

    l.Dl=l.numl*step;///kpc
    l.xls=l.Dl/s.Ds;
    l.RE=sqrt(4.0*G*l.Ml*(1.0+l.q)*Msun*s.Ds*Pc*1000.0)/velocity;
    l.RE=l.RE*sqrt(l.xls*(1.0-l.xls));///meter
    vrel(s,l);
    l.tE=l.RE/fabs(l.Vt*1000.0*3600.0*24.0);///in day
    l.proj=double(l.xls*Rsun/l.RE);
    s.ro_star=fabs(s.Rs*l.proj);
    l.t0=0.0;
    l.ksi=RandR(0.0,359.0)*pi/180.0;///[radian]
    l.u0= RandR(0.0,0.1);//*s.ro_star;///*s.ro_star;
    l.tstar= fabs(l.tE*s.ro_star);///days
    //l.pt1=-1.3*l.tE;// l.tstar;///days
    //l.pt2=+1.3*l.tE;//l.tstar;///days

}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double RandN(double sigma, double N){
    double p, fp, frand;
    do{
    p=(double)rand()/((double)(RAND_MAX)+(double)(1))*N*sigma; ///[-N:N]
    fp=exp(-p*p/(2.*sigma*sigma));
    frand= (double)rand()/((double)(RAND_MAX)+(double)(1));
    }while(frand>fp);
    double sign= (double)rand()/((double)(RAND_MAX)+(double)(1));
    if(sign<0.5)     return(p);
    else             return(-p);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double TETA(double x, double y,int flag)
{
double alfa=0.0;

   if(y>=0.0 && x>0.0 )      alfa=atan(fabs(y)/fabs(x));
   else if(y>0.0  && x<=0.0) alfa=atan(fabs(x)/fabs(y))+M_PI/2.0;
   else if(y<=0.0 && x<0.0 ) alfa=atan(fabs(y)/fabs(x))+M_PI;
   else if(y<0.0 &&  x>=0.0) alfa=atan(fabs(x)/fabs(y))+3.0*M_PI/2.0;
   if(x==0.0 && y==0.0) alfa=0.0;
   if(x==0.0 && y>0.0 ) alfa=M_PI/2.0;
   if(x==0.0 && y<0.0 ) alfa=3.0*M_PI/2.0;
   if(y==0.0 && x>0.0 ) alfa=0.0;
   if(y==0.0 && x<0.0 ) alfa=M_PI;///in radian
   if(flag==1) {
   alfa=alfa*0.5*RA;///in degree;فقط باید در ۱/۲ ضرب شود وهیچ تغییر دیگری نیاز ندارد
   if(alfa>360.0) alfa=alfa-360.0;
   if(alfa>360.0) alfa=alfa-360.0;
   if(alfa>180.0) alfa=alfa-180.0;
   if(alfa>90.0 ) alfa=alfa-180.0;}
   return(alfa);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
double RandR(double down, double up){
    double p =(double)rand()/((double)(RAND_MAX)+(double)(1.0));
    return(p*(up-down)+down);
}
///==============================================================//
///                                                              //
///                  Linear interpolarion                        //
///                                                              //
///==============================================================//
double Interpol(double ds, extinc & ex){
  double F=-1.0;
  if(ds<ex.dis[0])        F=ex.Extks[0];
  else if(ds>=ex.dis[99]) F=ex.Extks[99];
  else{ 
  for(int i=0; i<99; ++i){
  if(ex.dis[i]>=ex.dis[i+1]){
  cout<<"ERROR dis[i]: "<<ex.dis[i]<<"\t disi+1: "<<ex.dis[i+1]<<endl;  int yye; cin>>yye; }
  if(ds>=ex.dis[i] && ds<ex.dis[i+1]){
  F = ex.Extks[i]+(double)(ds-ex.dis[i])*(ex.Extks[i+1]-ex.Extks[i])/(ex.dis[i+1]-ex.dis[i]);
  break;}}}
  if(F==-1.0||F<0.0){cout<<"ERROR big Extinction(ds): "<<F<<"\t ds: "<<ds<<endl; int uut; cin>>uut;}
  return(F);
}
///==============================================================//
///                                                              //
///                  EXtinction                                 //
///                                                              //
///==============================================================//
int Extinction(extinc & ex,source & s)
{
     double sig,Lon,Lat;
     int uue, flag=0;
     if(s.lon<0.0){sig=-1.0;cout<<"Strange!!!!longtitude is negative:  s.lon "<<s.lon<<endl; cin>>uue;}
     else sig=1.0;
     double delt=fabs(s.lon)-floor(fabs(s.lon));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR longtitude: delt: "<<delt<<"\t s.lon: "<<s.lon<<endl;  cin>>uue; }
     else if(delt<0.25) Lon=(floor(fabs(s.lon))+0.00)*sig;
     else if(delt<0.50) Lon=(floor(fabs(s.lon))+0.25)*sig;
     else if(delt<0.75) Lon=(floor(fabs(s.lon))+0.50)*sig;
     else               Lon=(floor(fabs(s.lon))+0.75)*sig;
     if(s.lon==0.0 or Lon<0.0 or Lon<0.1  or s.lon<0.1) Lon=360.00;

     if(s.lat<0.0) sig=-1.0;
     else sig=1.0;
     delt=fabs(s.lat)-floor(fabs(s.lat));
     if(delt>1.0 || delt<0.0) {cout<<"ERROR latitude: delt: "<<delt<<"\t s.lon: "<<s.lat<<endl;  cin>>uue;}
     else if(delt<0.25)  Lat=(floor(fabs(s.lat))+0.00)*sig;
     else if(delt<0.50)  Lat=(floor(fabs(s.lat))+0.25)*sig;
     else if(delt<0.75)  Lat=(floor(fabs(s.lat))+0.50)*sig;
     else                Lat=(floor(fabs(s.lat))+0.75)*sig;
     if(Lon>360.000 || Lon<0.25 || fabs(Lat)>10.0 || (Lon>100 && Lon<260)){
     cout<<"BIG error (stopped program) s.lon: "<<Lon<<"\t s.lat: "<<Lat<<endl;   cin>>uue;}


     char filename[40];
     FILE *fpd;
     sprintf(filename,"./files/Ext/%c%c%c%.2lf%c%.2lf.dat",'E','x','t',Lat,'_',Lon);
     fpd=fopen(filename,"r");
     cout<<"Lat:  "<<Lat<<"\t Lon:    "<<Lon<<endl;

     double lonti,latit;
     if(!fpd){
     cout<<"cannot open (extinction) file long : "<<Lon<<"\t latit: "<<Lat<<endl;
     FILE *SD;
     SD=fopen("./files/Ext/saved_direction.txt","r");
     for(int i=0; i<64881; ++i) {
     fscanf(SD,"%lf %lf \n",&latit,&lonti);
     if(fabs(Lat-latit)<0.1 && fabs(Lon-lonti)<0.1){
     cout<<"ERROR  long : "<<Lon<<"\t latit: "<<Lat<<endl;
     cout<<"Saved: Latii: "<<latit<<"\t lonti: "<<lonti<<endl; cin>>uue;}}
     flag=-1;}
     else{
     flag=1;
     for(int i=0; i<100; ++i){
     fscanf(fpd,"%lf  %lf\n",&ex.dis[i],&ex.Extks[i]);////Just extinctin in [Ks-band]
     //cout<<"distance: "<<ex.dis[i]<<"\t Extks: "<<ex.Extks[i]<<endl;
     if(ex.dis[i]<0.2  || ex.dis[i]>50.0 || ex.Extks[i]<0.0){
     cout<<"dis: "<<ex.dis[i]<<"\t extI: "<<ex.Extks[i]<<"\t i: "<<i<<endl;
     cout<<"filename: "<<filename<<endl;  ex.Extks[i]=0.0; }
    // if(ex.dis[i]==8.25) ex.exks8= ex.Extks[i]; 
     }}
     cout<<">>>>>>>>>>>>>>> END OF EXTINCTION FUNCTION <<<<<<<<<<<<<<<<<<"<<endl;
     fclose(fpd);
     return(flag);
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       relative velocity                        ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void vrel(source & s, lens & l)
{
if (l.Dl==0.0) l.Dl=0.00034735;
 double pi=M_PI;
 double Rlc=sqrt(l.Dl*l.Dl*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*l.Dl*cos(s.TET)*cos(s.FI));
 double Rsc=sqrt(s.Ds*s.Ds*cos(s.FI)*cos(s.FI)+R_sun*R_sun-2.*R_sun*s.Ds*cos(s.TET)*cos(s.FI));
 if(Rlc==0.0) Rlc=0.00034346123;
 if(Rsc==0.0) Rsc=0.0004762654134;
 ///Source and Lens velocity components in Galactocentric cylindrical coordinates
 double SVT, SVR, SVZ, LVT, LVR, LVZ,SVt,LVt;
 ///Source and Lens velocity components in heliocenteric galactic coordinates
 double SVb, SVl, LVb, LVl;
 double fv, testfv;
 double VSunl,VSunt,VSunb,vls_b,vls_l;
 double betal,betas,deltal,deltas,deltao;



 double NN=3.0;
 double VSunR =-10.3;
 double VSunT =vro_sun*(1.00762+0.00712)+6.3;
 double VSunZ = 5.9;
 double sigma_R_Disk=43.0, sigma_T_Disk=27.8, sigma_Z_Disk=17.5;
 double sigma_R_TDisk= 67.0, sigma_T_TDisk= 51.0, sigma_Z_TDisk= 42.0;
 double sigma_R_halo= 131.0, sigma_T_halo= 106.0, sigma_Z_halo= 85.0;
 double sigma_R_Bulge = 113.0,sigma_T_Bulge = 115.0,sigma_Z_Bulge = 100.0;
 double Rho[8]={00.0}; double maxr=0.0;
 for(int i=0; i<8; ++i){
 Rho[i]=rho0[i]*corr[i]/d0[i]; maxr=maxr+ Rho[i];}




  double test= ((double)rand()/(double)(RAND_MAX+1.))*maxr; ///total ages
     if(test<=Rho[0])     {sigma_R_Disk=16.7; sigma_T_Disk=10.8; sigma_Z_Disk=6.0;}
else if(test<=(Rho[0]+Rho[1])) {sigma_R_Disk=19.8; sigma_T_Disk=12.8; sigma_Z_Disk=8.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]))  {sigma_R_Disk=27.2; sigma_T_Disk=17.6; sigma_Z_Disk=10.0;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]))
                           {sigma_R_Disk=30.2; sigma_T_Disk=19.5; sigma_Z_Disk=13.2;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]))
                           {sigma_R_Disk=36.7; sigma_T_Disk=23.7; sigma_Z_Disk=15.8;}
else if(test<=(Rho[0]+Rho[1]+Rho[2]+Rho[3]+Rho[4]+Rho[5]))
                           {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.4;}
else if(test<=maxr)        {sigma_R_Disk=43.1; sigma_T_Disk=27.8; sigma_Z_Disk=17.5;}
else  {cout<<"BIG ERROR "<<test<<"\t maxr: "<<maxr<<endl;  int yye; cin>>yye;}


///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
/// Generate Source velocity components in Glactocenteric cylindrical coordinates(x',y')
    if(s.struc==0){///Galactic disk
    do{
    SVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Disk;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_Disk*sigma_R_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
    do{
    SVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Disk;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_Disk*sigma_T_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    SVZ = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Disk;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_Disk*sigma_Z_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    SVT =SVT +vro_sun*(1.00762 * pow(Rsc/R_sun,0.0394) + 0.00712);}

    else if(s.struc==1){///Galactic bulge
    do{
    SVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Bulge;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_Bulge*sigma_Z_Bulge));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Bulge;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_Bulge*sigma_R_Bulge));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Bulge;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_Bulge*sigma_T_Bulge));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}

    else if(s.struc==2) {///thick disk
    do{
    SVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_TDisk;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_TDisk*sigma_R_TDisk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    SVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_TDisk;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_TDisk*sigma_T_TDisk));
    testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
    do{
    SVZ = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_TDisk;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_TDisk*sigma_Z_TDisk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    SVT =SVT+ vro_sun *(1.00762*pow(Rsc/R_sun,0.0394) + 0.00712); }

    else if(s.struc==3){///stellar halo
    do{
    SVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_halo;
    fv = exp(-1./2.*SVZ*SVZ/(sigma_Z_halo*sigma_Z_halo));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_halo;
    fv = exp(-1./2.*SVR*SVR/(sigma_R_halo*sigma_R_halo));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
    do{
    SVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_halo;
    fv = exp(-1./2.*SVT*SVT/(sigma_T_halo*sigma_T_halo));
    testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}

    l.vs=sqrt(SVR*SVR+SVT*SVT+SVZ*SVZ);


///======================================================================================
/// Generate Lens velocity components in Glactocenteric cylindrical coordinates(x',y')
if(l.struc==0){///Galactic disk
    do{
    LVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Disk;
    fv = exp(-1./2.*LVR*LVR/(sigma_R_Disk*sigma_R_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    LVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Disk;
    fv = exp(-1./2.*LVT*LVT/(sigma_T_Disk*sigma_T_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    do{
    LVZ =(-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Disk;
    fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_Disk*sigma_Z_Disk));
    testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
    LVT =LVT+ vro_sun *(1.00762 * pow(Rlc/R_sun,0.0394) + 0.00712);}

   else if(l.struc==1){///Galactic bulge
   do{
   LVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_Bulge;
   fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_Bulge*sigma_Z_Bulge));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_Bulge;
   fv = exp(-1./2.*LVR*LVR/(sigma_R_Bulge*sigma_R_Bulge));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_Bulge;
   fv = exp(-1./2.*LVT*LVT/(sigma_T_Bulge*sigma_T_Bulge));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}

   else if(l.struc==2){///thick disk
   do{
   LVR = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_TDisk;
   fv = exp(-1./2.*LVR*LVR/(sigma_R_TDisk*sigma_R_TDisk));
   testfv = (double)rand()/(double)(RAND_MAX+1.);}while(testfv>fv);
   do{
   LVT = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_TDisk;
   fv = exp(-1./2.*LVT*LVT/(sigma_T_TDisk*sigma_T_TDisk));
   testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
   do{
   LVZ = (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_TDisk;
   fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_TDisk*sigma_Z_TDisk));
   testfv = (double)rand()/(double)(RAND_MAX+1.); }while(testfv>fv);
   LVT =LVT+ vro_sun *(1.00762*pow(Rlc/R_sun,0.0394) + 0.00712); }

   else if(l.struc==3){///stellar halo
   do{
   LVZ= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_Z_halo;
   fv = exp(-1./2.*LVZ*LVZ/(sigma_Z_halo*sigma_Z_halo));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVR= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_R_halo;
   fv = exp(-1./2.*LVR*LVR/(sigma_R_halo*sigma_R_halo));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);
   do{
   LVT= (-1. + 2.*(double)rand()/(double)(RAND_MAX+1.))*NN*sigma_T_halo;
   fv = exp(-1./2.*LVT*LVT/(sigma_T_halo*sigma_T_halo));
   testfv =((double)rand()/(double)(RAND_MAX+1.)); }while(testfv>fv);}
   l.vl=sqrt(LVT*LVT+LVZ*LVZ+LVR*LVR);
///HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

   if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc-1.0)<0.01) betal=pi/2.0;
   else if(fabs(l.Dl*cos(s.FI)*sin(s.TET)/Rlc+1.0)<0.01) betal=-pi/2.0;
   else  betal=asin(l.Dl*cos(s.FI)*sin(s.TET)/Rlc);///lens[-pi/2,pi/2]
   if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc-1.0)<0.01) betas=pi/2.0;
   else if(fabs(s.Ds*cos(s.FI)*sin(s.TET)/Rsc+1.0)<0.01) betas=-pi/2.0;
   else  betas=asin(s.Ds*cos(s.FI)*sin(s.TET)/Rsc);///lens[-pi/2,pi/2]
    if(fabs(l.Dl*cos(s.FI)*sin(s.TET))>Rlc || fabs(s.Ds*cos(s.FI)*sin(s.TET))>Rsc || Rlc==0.0 || Rsc==0.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl;
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
    cout<<"sin(l): "<<l.Dl*cos(s.FI)*sin(s.TET)/Rlc<<"\t sin(s): "<<s.Ds*cos(s.FI)*sin(s.TET)/Rsc<<endl;
     //int ew; cin>>ew;
      }

       //betao=0.0; ///observer
    if(fabs(l.Dl*cos(s.FI))>sqrt(pow(R_sun,2.0)+pow(l.Dl*cos(s.FI)*sin(s.TET),2.0)) ) betal= pi-betal;
    if(fabs(s.Ds*cos(s.FI))>sqrt(pow(R_sun,2.0)+pow(s.Ds*cos(s.FI)*sin(s.TET),2.0)) ) betas= pi-betas;



    if(fabs((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI))-1.0)<0.01)   deltal=0.0;
    else if (fabs((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI))+1.0)<0.01) deltal=pi;
    else    deltal=acos((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI)));
    if(fabs((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI))-1.0)<0.01)   deltas=0.0;
    else if (fabs((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI))+1.0)<0.01) deltas=pi;
    else    deltas=acos((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI)));
   if(fabs((Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI)))>1.002 || fabs((Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI)))>1.002  || l.Dl==0.0 || s.Ds==0.0 || fabs(s.FI)==pi/2.0){
    cout<<"ERROR Dl: "<<l.Dl<<"\t Ds: "<<s.Ds<<endl;
    cout<<"betal: "<<betal<<"\t betas: "<<betas<<endl;
    cout<<"FI: "<<s.FI<<"\t TET: "<<s.TET<<endl;
    cout<<"Rlc: "<<Rlc<<"\t Rsc: "<<Rsc<<endl;
    cout<<"cos(dl): "<<(Rlc-R_sun*cos(betal))/(l.Dl*cos(s.FI))<<"\t cos(ds): "<<(Rsc-R_sun*cos(betas))/(s.Ds*cos(s.FI))<<endl;
    // int ew; cin>>ew;
     }


    deltao=pi/2.0;
    SVl=-SVR*    sin(deltas)+ SVT* cos(deltas);
    LVl=-LVR*    sin(deltal)+ LVT* cos(deltal);
    VSunl=-VSunR*sin(deltao)+VSunT*cos(deltao);

    SVt=  1.0*SVR*cos(deltas)+  SVT*sin(deltas);
    LVt=  1.0*LVR*cos(deltal)+  LVT*sin(deltal);
    VSunt=1.0*VSunR*cos(deltao)+VSunT*sin(deltao);

    SVb=-sin(s.FI)*(SVt) + cos(s.FI)* (SVZ);
    LVb=-sin(s.FI)*(LVt) + cos(s.FI)* (LVZ);
    VSunb=-sin(s.FI)*(VSunt)+cos(s.FI)*(VSunZ);

    vls_l= LVl-l.xls*SVl -(1.0-l.xls)*VSunl;
    vls_b= LVb-l.xls*SVb -(1.0-l.xls)*VSunb;
    l.Vt=sqrt(fabs(vls_l*vls_l+ vls_b*vls_b));

if (l.Vt<0.0 || l.Vt>1.0e6 ){
    cout<<" Vt is very large: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;   int yee; cin>>yee;}
//cout<<"Vt: "<<l.Vt<<"\t vl: "<<l.vl<<"\t Vs: "<<l.vs<<endl;
}
///==================================================================
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
///                                                                ///
///                       Glactic model                            ///
///                                                                ///
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
void Disk_model(source & s)
{
   double x,xb,yb,zb,r4,r2,rdi,Rb;
   double nn=0.4/0.8;
   double mBarre;///stars/pc^3
   double Rx0, Ry0, Rz0,Rc,cp,cn,rhoE,rhoS;
   double alfa=12.89/RA;
   double xf,yf,zf,rho;
   s.Romaxs=s.Nstart=s.Rostart=0.0;
   double fd=1.0; ///see the program mass_averaged.cpp. we do not apply any limitation
   double fb=1.0;///0.657066;////just stars brighter than V=11.5, but we change to consider all stars
   double fh=1.0;///No limitation
   double Rdd=2.17;///2.53;///2.17;
   double Rhh=1.33;///1.32;//1.33;
/*
I assume that the up-limit of the mass is indicated by the simulation. Because it depends on the structre, .... all details can be found in mass_averaged.cpp  code. */

   /*
    char filename[40];
    FILE *fill;
    sprintf(filename,"./files/density/%c%.2lf%c%.2lf.dat",'D',s.lat,'_',s.lon);
    fill=fopen(filename,"w");
    if(!fill){cout<<"cannot open file longtitude : "<<s.lon<<"\t latitude: "<<s.lat<<endl;  exit(0);}
   */



for(int i=1;i<Num;++i){
   s.Rostar0[i]=s.Rostari[i]=s.Nstari[i]=0.0;
   s.rho_disk[i]=s.rho_bulge[i]=s.rho_halo[i]=s.rho_ThD[i]=0.0;

   x=i*step;
   zb = sin(s.FI)*x;
   yb = cos(s.FI)*sin(s.TET)*x;
   xb = R_sun-x*cos(s.FI)*cos(s.TET);
   Rb=sqrt(xb*xb+yb*yb);


///========== Galactic Thin Disk =====================
   for(int ii=0; ii<8; ++ii){
   rdi=Rb*Rb+zb*zb/(epci[ii]*epci[ii]);
   if(ii==0)     rho=exp(-rdi/25.0)-exp(-rdi/9.0);
   else if(ii>0) rho=exp(-sqrt(0.25+rdi/(Rdd*Rdd)))-exp(-sqrt(0.25+rdi/(Rhh*Rhh)));
   s.rho_disk[i]=s.rho_disk[i]+ rho0[ii]*corr[ii]*0.001*rho/d0[ii];}///M_sun/pc^3
///=================================================


///========== Galactic Thick Disk =====================

  double rho00=1.34*0.001+3.04*0.0001;
  if(fabs(zb)<0.4) s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*(1.0-zb*zb/(0.4*0.8*(2.0+nn)));
  else s.rho_ThD[i]=(rho00/0.999719)*exp(-(Rb-R_sun)/2.5)*exp(nn)*exp(-fabs(zb)/0.8)/(1.0+0.5*nn);///M_sun/pc^3
///=================================================


///========== Galactic Stellar Halo=================
   rdi=sqrt(Rb*Rb+ zb*zb/(0.76*0.76));
   if( rdi <=0.5)  s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(0.5/R_sun,-2.44);
   else            s.rho_halo[i]=1.0*(0.932*0.00001/867.067)*pow(rdi/R_sun,-2.44);///M_sun/pc^3
///=================================================



///========== Galactic bulge =====================
   xf = xb*cos(alfa) + yb*sin(alfa);
   yf =-xb*sin(alfa) + yb*cos(alfa);
   zf = zb;
   Rx0=1.46, Ry0=0.49, Rz0=0.39; Rc=3.43; cp=3.007;  cn=3.329;  mBarre=35.45/(3.84723);
   r4=pow(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(fabs(r4),1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4));
   else       rhoS= mBarre*1.0/(cosh(-r4)*cosh(-r4))*exp(-4.0*(r2-Rc)*(r2-Rc));

   Rx0=4.44, Ry0=1.31, Rz0=0.80; Rc=6.83; cp=2.786; cn=3.917; mBarre=2.27/87.0;//85.3789;
   r4=pow(fabs(pow(fabs(xf/Rx0),cn)+pow(fabs(yf/Ry0),cn)),cp/cn)+pow(fabs(zf/Rz0),cp);
   r4=pow(r4,1.0/cp);
   r2=sqrt(fabs(xf*xf+yf*yf));
   if(r2<=Rc) rhoE= mBarre*exp(-r4);
   else       rhoE= mBarre*exp(-r4)*exp(-4.0*(r2-Rc)*(r2-Rc));
  s.rho_bulge[i]= fabs(rhoS)+fabs(rhoE);///M_sun/pc^3
///=================================================
///در اینجا اینکه تعداد ستاره ه
///ا به این بستگی دارد که ما  نمودار قدر رنگ مطلق ستاره ها را چگونه درست کرده باشیم.اگر هیچ گونه محدودیتی برای
///درست کردن آن  در نظر نگرفته ایم،  پس تعداد کل ستاره ها را نظر  میگیریم.
/// ولی بهتر است که ما رابطه بین قدر مطلق و جرم را تعیین کنیم. در این صورت می توانیم  ورودی قدر رنگ
///وارد شده به کد را خودمان محدود به ستاره های روشن کنیم تا سرعت اجرای برنامه بالارود.
///averaged mass are the same as the previous work!!! because we did not change the besancon model


s.Rostar0[i]=fabs(s.rho_disk[i])+fabs(s.rho_ThD[i])+fabs(s.rho_bulge[i])+fabs(s.rho_halo[i]);///[Msun/pc^3]
s.Rostari[i]=s.Rostar0[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Msun/deg^2]
s.Nstari[i]=binary_fraction*(s.rho_disk[i]*fd/0.403445+s.rho_ThD[i]*fh/0.4542+s.rho_halo[i]*fh/0.4542+s.rho_bulge[i]*fb/0.308571);////[Nt/pc^3]

s.Nstari[i]= s.Nstari[i]*x*x*step*1.0e9*(M_PI/180.0)*(M_PI/180.0);///[Ni/deg^2]

s.Nstart  +=s.Nstari[i];///[Nt/deg^2]
s.Rostart += s.Rostari[i];///[Mt/deg^2]
if(s.Rostari[i]>s.Romaxs) s.Romaxs=s.Rostari[i];///source selection
//fprintf(fill,"%e   %e   %e   %e   %e  %e   %e\n",x,s.rho_disk[i],s.rho_bulge[i],s.rho_ThD[i],s.rho_halo[i],s.Rostar0[i],s.Nstari[i]);
   }
 // fclose(fill);
}
///>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
///==============================================================//
///                                                              //
///                  READ CMD FILE                               //
///                                                              //
///==============================================================//
void read_cmd(CMD & cm){

////mass log10(T) log10(age) log10(L) log10(g) metal U G R I Z Bj Vj Rj Ij CL TYP
    int yye, uui,h, k1, k2, g; double metal,age, mk;   
    char filename[40];
    FILE *fp2;

    int number[70]={0};   int count[70]={0};   double Metal[70]={0.0}; 
    FILE *meta; 
    meta=fopen("./files/CMD_BESANCON_new/metal.txt","r"); 
    for(int i=0; i<70; ++i){
    fscanf(meta,"%lf %d  %d\n",&Metal[i],&count[i],&number[i]);    
    if((Metal[i]<Metal[i-1] and i>0) or float(Metal[i])<-0.001 or number[i]==0 or count[i]>YZ or 
       (abs(count[i-1]+number[i-1]-count[i])>2 and i>0)){
    cout<<"ERROR Metal[i]: "<<Metal[i]<<"\t count[i]: "<<count[i]<<"\t number[i]: "<<number[i]<<"\t i: "<<i<<endl; cin>>uui;} }
    fclose(meta);
    double Age[YZ]={0.0}; double B1[YZ]={0.0};  double M1[YZ]={0.0};   double mm[YZ]={0.0}; 
    FILE *ji; 
    ji=fopen("./files/CMD_BESANCON_new/RVI.txt", "r"); 
    for(int i=0; i<YZ; ++i){
    fscanf(ji,"%lf   %lf   %lf  %lf\n",&Age[i],&mm[i],&B1[i],&M1[i]); 
    if(Age[i]<0.0 or mm[i]<0.0 or fabs(B1[i])>1.7 or M1[i]<0.5 or Age[i]>18.0){   
    cout<<"ERROR Age(JI): "<<Age[i]<<"\t metal: "<<mm[i]<<"\t B[i]"<<B1[i]<<"\t M[i]: "<<M1[i]<<"\t i: "<<i<<endl;
    cin>>uui;}}
    fclose(ji); 



////=================================== THIN DISK ==============================

    int j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c%c.dat",'C','M','D','T','i');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTi.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_d[j],&cm.Teff_d[j],&age,&cm.logl_d[j],&cm.gra_d[j],&metal,&cm.Rs_d[j],
    &cm.Mab_d[0][j],&cm.Mab_d[1][j],&cm.Mab_d[2][j],&cm.Mab_d[4][j],&mk,&cm.cl_d[j],&cm.type_d[j]);
///*******************************************************
    h=-1; 
    if(metal<Metal[0] or metal==Metal[0])         h=0; 
    else if(metal>Metal[69] or  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.05){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_d[3][j]= double(B1[g]+M1[g]*(cm.Mab_d[2][j]+cm.Mab_d[4][j])*0.5); ///R-band versus (V+I)
    cm.slop_d[j] = double(2.0/M1[g]);//slope of V+ I versus R
    if(fabs(cm.Mab_d[4][j]-cm.Mab_d[3][j])>2.5 or fabs(age-Age[g])>3.0 or cm.slop_d[j]<1.5  or cm.slop_d[j]>2.5){
    cout<<"slope_d:   "<<cm.slop_d[j]<<endl;
    cout<<"ERROR: Mab_d(z-band): "<<cm.Mab_d[4][j]<<"\t Mab[3]:  "<<cm.Mab_d[3][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************   
    if(cm.mass_d[j]<0.0 or cm.mass_d[j]==0.0 or cm.gra_d[j]>6.0 or cm.Teff_d[j]<0.0 or metal>0.12 or age>10.0 or age<0.0 or cm.cl_d[j]>5 or cm.type_d[j]>=9.0 or cm.type_d[j]<2.0 or (cm.cl_d[j]==5 and int(cm.type_d[j])>7) or (cm.cl_d[j]==6 and int(cm.type_d[j])!=9) or (cm.cl_d[j]<5 and int(cm.type_d[j])==9)){
    cout<<"ERROR(reading cmd file) structure thin disk: "<<"\t Nfil: "<<j+1<<endl; 
    cout<<"type_d: "<<cm.type_d[j]<<"\t CL(thin disk): "<<cm.cl_d[j]<<endl;
    cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N1){
    cout<<"BIG ERRROR j: "<<j<<"\t N1: "<<N1<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thin disk):  No. rows file: "<<j<<endl;







////=================================== Galactic Bulge ===========================================================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c.dat",'C','M','D','b');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDB.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_b[j],&cm.Teff_b[j],&age,&cm.logl_b[j],&cm.gra_b[j],&metal,&cm.Rs_b[j],
    &cm.Mab_b[0][j],&cm.Mab_b[1][j],&cm.Mab_b[2][j],&cm.Mab_b[4][j],&mk,&cm.cl_b[j],&cm.type_b[j]);
///*******************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.4){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_b[3][j]= double(B1[g]+M1[g]*(cm.Mab_b[2][j]+cm.Mab_b[4][j])*0.5); ///R-band versus (V+I)/2
    cm.slop_b[j] = double(2.0/M1[g]);//slope of V+ I versus R
    if(fabs(cm.Mab_b[4][j]-cm.Mab_b[3][j])>2.5 or fabs(age-Age[g])>3.0 or cm.slop_b[j]<1.5 or cm.slop_b[j]>2.5){
    cout<<"ERROR slope_b:  "<<cm.slop_b[j]<<endl;
    cout<<"ERROR:  Mab_b(z-band): "<<cm.Mab_b[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************   
    if(cm.mass_b[j]<0.0 or cm.mass_b[j]==0.0 or cm.Teff_b[j]<0.0 or age>10.0 or metal>0.2 or cm.cl_b[j]>5 or cm.type_b[j]>8.0 or (cm.cl_b[j]==5 and cm.type_b[j]>8.0) or (cm.cl_b[j]==6 and int(cm.type_b[j])!=9) or (cm.cl_b[j]<5 and int(cm.type_b[j])==9)){
    cout<<"ERROR(reading cmd file) structure bulge: "<<"\t Nfil: "<<j+1<<endl;
    cout<<"type_b: "<<cm.type_b[j]<<"\t CL(bulge): "<<cm.cl_b[j]<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N2){cout<<"BIG ERRROR j: "<<j<<"\t N2: "<<N2<<endl;  cin>>yye;}
    cout<<"End of CMD reading (bulge):  No. rows file: "<<j<<endl;










////==================================== THICK DISK ==========================================================
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c%c.dat",'C','M','D','T','k');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDTk.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_t[j],&cm.Teff_t[j],&age,&cm.logl_t[j],&cm.gra_t[j],&metal,&cm.Rs_t[j],
    &cm.Mab_t[0][j],&cm.Mab_t[1][j],&cm.Mab_t[2][j],&cm.Mab_t[4][j],&mk,&cm.cl_t[j],&cm.type_t[j]);
///*******************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.4){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_t[3][j]= double(B1[g]+M1[g]*(cm.Mab_t[2][j]+cm.Mab_t[4][j])*0.5); ///R-band versus (V+I)/2
    cm.slop_t[j] = double(2.0/M1[g]);//slope of V+ I versus R
    if(fabs(cm.Mab_t[4][j]-cm.Mab_t[3][j])>2.5 or fabs(age-Age[g])>3.0 or cm.slop_t[j]<1.5  or cm.slop_t[j]>2.5){
    cout<<"Error slope_t: "<<cm.slop_t[j]<<endl;
    cout<<"ERROR: Mab_t(z-band): "<<cm.Mab_t[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************   
    if(cm.mass_t[j]<0.0||  cm.mass_t[j]==0.0 or cm.Teff_t[j]<0.0 or metal>0.06 || cm.cl_t[j]>5 || cm.type_t[j]>8.0 or 
    (cm.cl_t[j]==5 and float(cm.type_t[j])>8.0) or (cm.cl_t[j]==6 and int(cm.type_t[j])!=9) or (cm.cl_t[j]<5 and int(cm.type_t[j])==9)){
    cout<<"type_thick: "<<cm.type_t[j]<<"\t CL(thick): "<<cm.cl_t[j]<<endl;
    cout<<"mass: "<<cm.mass_t[j]<<"\t TefF:  "<<cm.Teff_t[j]<<"\t metal:  "<<metal<<endl;
    cout<<"ERROR(reading cmd file) structure thick disk: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N3){cout<<"BIG ERRROR j: "<<j<<"\t N3: "<<N3<<endl;  cin>>yye;}
    cout<<"End of CMD reading (thick disk):  No. rows file: "<<j<<endl;






////=================================== STELLAR HALO ============================================================ 
    j=0; 
    sprintf(filename,"./files/CMD_BESANCON_new/%c%c%c%c.dat",'C','M','D','h');
    fp2=fopen(filename,"r");
    if(!fp2){cout<<"cannot read CMDH.dat long: "<<endl; cin>>yye;}
    while(!feof(fp2)){
    fscanf(fp2,"%lf %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf  %lf %lf %lf %d %lf\n",
    &cm.mass_h[j],&cm.Teff_h[j],&age,&cm.logl_h[j],&cm.gra_h[j],&metal,&cm.Rs_h[j],
    &cm.Mab_h[0][j],&cm.Mab_h[1][j],&cm.Mab_h[2][j],&cm.Mab_h[4][j],&mk,&cm.cl_h[j],&cm.type_h[j]);
///*******************************************************
    h=-1; 
    if(metal<Metal[0] || metal==Metal[0])         h=0; 
    else if(metal>Metal[69] ||  metal==Metal[69]) h=69;
    else{ 
    for(int i=1; i<70; ++i){
    if((metal-Metal[i-1])*(metal-Metal[i])<0.0 || metal==Metal[i-1]){h=int(i-1); break;}}}
    if(h==-1 or h<0 or h>69){cout<<"ERROR1 h: "<<h<<"\t metal: "<<metal<<endl; cin>>uui;} 
    g=-1; k1=int(count[h]);   k2=int(count[h]+number[h]);
    if(k1==k2 or k1>k2 or k1<0 or k2>YZ or mm[k1]!=mm[k2-1] or number[h]==0 or fabs(metal-Metal[h])>0.4){
    cout<<"ERROR4: k1: "<<k1<<"\t k2: "<<k2<<"\t h: "<<h<<"\t count[h]: "<<count[h]<<"\t number[h]: "<<number[h]<<endl; 
    cout<<"age: "<<age<<"\t Age[k1]: "<<Age[k1]<<"\t Age[k2-1]: "<<Age[k2-1]<<endl; 
    cout<<"metal[k1]: "<<mm[k1]<<"\t metal[k2-1]: "<<mm[k2-1]<<endl;  cin>>uui;}
    if(age<Age[k1]  or age==Age[k1])  g=k1; 
    else if(age>Age[k2-1] or age==Age[k2-1]) g=int(k2-1); 
    else{
    for(int k=k1+1;  k<k2;  ++k){
    if(Age[k-1]>Age[k] or mm[k-1]!=mm[k]){
    cout<<"Bad error: Age[k-1]: "<<Age[k-1]<<"\t Age[k]: "<<Age[k]<<"\t mm[k-1]:"<<mm[k-1]<<"\t mm[k]"<<mm[k]<<endl; cin>>uui;} 
    if((age-Age[k-1])*(age-Age[k])<0.0  ||  age==Age[k-1]) {g=int(k-1);   break;}}}
    if(g==-1 or g<0 or g>(YZ-1) or g<k1 or g>k2){cout<<"ERROR2 g:"<<g<<"\t age:"<<age<<"\t metal: "<<metal<<endl; cin>>uui;} 
    cm.Mab_h[3][j]= double(B1[g]+M1[g]*(cm.Mab_h[2][j]+cm.Mab_h[4][j])*0.5); ///R-band versus (V+I)/2.
    cm.slop_h[j] = double(2.0/M1[g]);//slope of V+ I versus R
    if(fabs(cm.Mab_h[4][j]-cm.Mab_h[3][j])>2.5 or fabs(age-Age[g])>3.0  or cm.slop_h[j]<1.5  or cm.slop_h[j]>2.5){
    cout<<"slope_h:  "<<cm.slop_h[j]<<endl;
    cout<<"ERROR: Mab_h(z-band): "<<cm.Mab_h[4][j]<<"\t metal: "<<metal<<endl; 
    cout<<"age: "<<age<<"\t Age[g]: "<<Age[g]<<"\t B[g]: "<<B1[g]<<"\t M[g]: "<<M1[g]<<"\t g: "<<g<<"\t Metal:"<<Metal[h]<<endl;
    cin>>uui;}
///*******************************************************       
    if(cm.mass_h[j]<0.0 || cm.mass_h[j]==0.0 || age<13.0 or age>15.0 or cm.cl_h[j]<0  or cm.cl_h[j]>5  or  cm.Teff_h[j]<0.0 ||
    metal>0.04 || cm.cl_h[j]>8 || cm.type_h[j]>9 or (cm.cl_h[j]==5 and int(cm.type_h[j])>7) or (cm.cl_h[j]==6 and int(cm.type_h[j])!=9) or   (cm.cl_h[j]<5 and int(cm.type_h[j])==9)){
    cout<<"type_halo: "<<cm.type_h[j]<<"\t CL(halo): "<<cm.cl_h[j]<<endl;
    cout<<"ERROR(reading cmd file) structure halo: "<<"\t Nfil: "<<j+1<<endl; cin>>yye;}
    j++;} fclose(fp2);
    if(j!=N4){cout<<"BIG ERRROR j: "<<j<<"\t N4: "<<N4<<endl;  cin>>yye;}
    cout<<"End of CMD reading (halo):  No. rows file: "<<j<<endl;
    cout<<">>>>>>>>>>>>>>>>> END OF CMD READING <<<<<<<<<<<<<<<<<<<<<<<<<"<<endl;      
}
///XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
