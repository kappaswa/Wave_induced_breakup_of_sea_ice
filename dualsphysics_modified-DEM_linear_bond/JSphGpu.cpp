//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpu.cpp \brief Implements the class \ref JSphGpu.

#include "JSphGpu.h"
#include "JException.h"
#include "JSphGpu_ker.h"
#include "JSphGpu_mdbc_iker.h"
#include "JSphGpuSimple_ker.h"
#include "JCellDivGpu.h"
#include "JDsGpuInfo.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include "JDsMotion.h"
#include "JDsFixedDt.h"
#include "JDsSaveDt.h"
#include "JDsOutputTime.h"
#include "JWaveGen.h"
#include "JMLPistons.h"
#include "JRelaxZones.h"
#include "JChronoObjects.h"
#include "JDsFtForcePoints.h"
#include "JDsDamping.h"
#include "JDsAccInput.h"
#include <algorithm>
#include "JXml.h"
#include "JDsGaugeSystem.h"
#include "JSphInOut.h"
#include "JSphShifting.h"
#include "JSphMk.h"
#include "JDataArrays.h"
#include "JSpVtkData.h"

#include <cfloat>
#include <climits>
#include <fstream>
#include <iomanip>
#include <vector>

using namespace std;

namespace{

struct StClusterLevels{
  bool ok;
  string reason;
  vector<double> levels;
  vector<unsigned> index;
  double spacing;
  StClusterLevels():ok(false),spacing(0){}
};

struct StBondAreaFit{
  bool ok;
  bool line2d;
  string reason;
  unsigned nu,nv;
  unsigned normalaxis,uaxis,vaxis;
  double lu,lv;
  double area_target,iuu_target,ivv_target,iuv_target;
  double area_disc,iuu_disc,ivv_disc,iuv_disc;
  vector<double> u;
  vector<double> v;
  vector<double> aeff;
  StBondAreaFit():ok(false),line2d(false),nu(0),nv(0),normalaxis(UINT_MAX),uaxis(UINT_MAX),vaxis(UINT_MAX),lu(0),lv(0),area_target(0),iuu_target(0),ivv_target(0),iuv_target(0)
    ,area_disc(0),iuu_disc(0),ivv_disc(0),iuv_disc(0){}
};

inline double Dot3(const tdouble3& a,const tdouble3& b){
  return(a.x*b.x + a.y*b.y + a.z*b.z);
}

inline double Norm3(const tdouble3& v){
  return(sqrt(Dot3(v,v)));
}

inline tdouble3 Normalize3(const tdouble3& v){
  const double m=Norm3(v);
  return(m>ALMOSTZERO? v/m: TDouble3(0));
}

inline tdouble3 Cross3(const tdouble3& a,const tdouble3& b){
  return(TDouble3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x));
}

inline tdouble3 PickLeastAlignedAxis(const tdouble3& n){
  const double ax=fabs(n.x), ay=fabs(n.y), az=fabs(n.z);
  if(ax<=ay && ax<=az)return(TDouble3(1,0,0));
  if(ay<=az)return(TDouble3(0,1,0));
  return(TDouble3(0,0,1));
}

inline double Coord3(const tdouble3& p,unsigned axis){
  return(axis==0? p.x: (axis==1? p.y: p.z));
}

inline const char* AxisName(unsigned axis){
  return(axis==0? "x": (axis==1? "y": "z"));
}

bool SolveLinear3(double a[3][3],double b[3],double x[3]){
  for(unsigned i=0;i<3;i++){
    unsigned pivot=i;
    double maxv=fabs(a[i][i]);
    for(unsigned r=i+1;r<3;r++){
      const double v=fabs(a[r][i]);
      if(v>maxv){ maxv=v; pivot=r; }
    }
    if(maxv<=ALMOSTZERO)return(false);
    if(pivot!=i){
      for(unsigned c=i;c<3;c++)swap(a[i][c],a[pivot][c]);
      swap(b[i],b[pivot]);
    }
    const double diag=a[i][i];
    for(unsigned c=i;c<3;c++)a[i][c]/=diag;
    b[i]/=diag;
    for(unsigned r=0;r<3;r++)if(r!=i){
      const double factor=a[r][i];
      if(fabs(factor)<=ALMOSTZERO)continue;
      for(unsigned c=i;c<3;c++)a[r][c]-=factor*a[i][c];
      b[r]-=factor*b[i];
    }
  }
  for(unsigned i=0;i<3;i++)x[i]=b[i];
  return(true);
}

StClusterLevels ClusterLevels(const vector<double>& coords,double tol){
  StClusterLevels out;
  const unsigned n=unsigned(coords.size());
  out.index.assign(n,UINT_MAX);
  if(!n){
    out.reason="empty coordinate list";
    return(out);
  }
  vector<pair<double,unsigned> > sorted;
  sorted.reserve(n);
  for(unsigned c=0;c<n;c++)sorted.push_back(make_pair(coords[c],c));
  sort(sorted.begin(),sorted.end());
  unsigned start=0;
  while(start<n){
    unsigned end=start+1;
    double sum=sorted[start].first;
    while(end<n && fabs(sorted[end].first-sorted[end-1].first)<=tol){
      sum+=sorted[end].first;
      end++;
    }
    const double level=sum/double(end-start);
    const unsigned idx=unsigned(out.levels.size());
    out.levels.push_back(level);
    for(unsigned c=start;c<end;c++)out.index[sorted[c].second]=idx;
    start=end;
  }
  if(out.levels.size()<=1){
    out.spacing=0;
    out.ok=true;
    return(out);
  }
  double sumsp=0;
  unsigned nsp=0;
  for(unsigned c=1;c<out.levels.size();c++){
    const double diff=out.levels[c]-out.levels[c-1];
    if(diff<=tol){
      out.reason="collapsed coordinate levels";
      return(out);
    }
    sumsp+=diff;
    nsp++;
  }
  out.spacing=(nsp? sumsp/double(nsp): 0.0);
  out.ok=true;
  return(out);
}

bool SolveMomentWeights(const vector<double>& coords,double target_sum,double target_second,vector<double>& weights,string& reason){
  const unsigned n=unsigned(coords.size());
  if(n<2){
    reason="need at least two levels per rectangular direction";
    return(false);
  }
  const double w0=target_sum/double(n);
  double ata[3][3]={{0,0,0},{0,0,0},{0,0,0}};
  double rhs[3]={0,0,0};
  for(unsigned c=0;c<n;c++){
    const double x0=1.0;
    const double x1=coords[c];
    const double x2=x1*x1;
    const double xs[3]={x0,x1,x2};
    for(unsigned i=0;i<3;i++){
      rhs[i]+=xs[i]*w0;
      for(unsigned j=0;j<3;j++)ata[i][j]+=xs[i]*xs[j];
    }
  }
  rhs[0]-=target_sum;
  rhs[1]-=0.0;
  rhs[2]-=target_second;
  double lambda[3]={0,0,0};
  if(!SolveLinear3(ata,rhs,lambda)){
    reason="singular least-deviation system";
    return(false);
  }
  weights.assign(n,0.0);
  double wmin=DBL_MAX;
  for(unsigned c=0;c<n;c++){
    const double corr=lambda[0] + lambda[1]*coords[c] + lambda[2]*coords[c]*coords[c];
    weights[c]=w0-corr;
    wmin=min(wmin,weights[c]);
  }
  const double negtol=max(target_sum*1e-9,1e-14);
  if(wmin<-negtol){
    reason="least-deviation weights turned negative";
    return(false);
  }
  for(unsigned c=0;c<n;c++)if(weights[c]<0)weights[c]=0;
  return(true);
}

StBondAreaFit FitRectangularBondAreas(const vector<tdouble3>& mids,const vector<tdouble3>& deltas,double aface,double dp
  ,double iuu_target_user,double ivv_target_user){
  StBondAreaFit fit;
  const unsigned n=unsigned(mids.size());
  fit.u.assign(n,0.0);
  fit.v.assign(n,0.0);
  fit.aeff.assign(n,0.0);
  if(!n){
    fit.reason="no bond midpoints";
    return(fit);
  }

  tdouble3 center=TDouble3(0);
  tdouble3 normal=TDouble3(0);
  for(unsigned c=0;c<n;c++){
    center=center+mids[c];
    normal=normal+deltas[c];
  }
  center=center/double(n);
  normal=Normalize3(normal);
  if(Norm3(normal)<=ALMOSTZERO){
    for(unsigned c=0;c<n && Norm3(normal)<=ALMOSTZERO;c++)normal=Normalize3(deltas[c]);
  }
  if(Norm3(normal)<=ALMOSTZERO){
    fit.reason="cannot estimate interface normal";
    return(fit);
  }

  tdouble3 t1=Normalize3(Cross3(PickLeastAlignedAxis(normal),normal));
  if(Norm3(t1)<=ALMOSTZERO)t1=Normalize3(Cross3(TDouble3(0,0,1),normal));
  if(Norm3(t1)<=ALMOSTZERO){
    fit.reason="cannot build local interface basis";
    return(fit);
  }
  const tdouble3 t2=Normalize3(Cross3(normal,t1));

  double cxx=0,cxy=0,cyy=0;
  vector<double> a(n,0.0),b(n,0.0);
  for(unsigned c=0;c<n;c++){
    const tdouble3 d=mids[c]-center;
    a[c]=Dot3(d,t1);
    b[c]=Dot3(d,t2);
    cxx+=a[c]*a[c];
    cxy+=a[c]*b[c];
    cyy+=b[c]*b[c];
  }
  const double theta=0.5*atan2(2.0*cxy,cxx-cyy);
  const double ct=cos(theta), st=sin(theta);
  vector<double> uraw(n,0.0), vraw(n,0.0);
  double umax=0,vmax=0;
  for(unsigned c=0;c<n;c++){
    uraw[c]= ct*a[c] + st*b[c];
    vraw[c]=-st*a[c] + ct*b[c];
    umax=max(umax,fabs(uraw[c]));
    vmax=max(vmax,fabs(vraw[c]));
  }
  const double tol=max(max(umax,vmax)*1e-6, max(double(dp)*1e-3,1e-12));
  const StClusterLevels ulevels_raw=ClusterLevels(uraw,tol);
  if(!ulevels_raw.ok){
    fit.reason="u clustering failed: "+ulevels_raw.reason;
    return(fit);
  }
  const StClusterLevels vlevels_raw=ClusterLevels(vraw,tol);
  if(!vlevels_raw.ok){
    fit.reason="v clustering failed: "+vlevels_raw.reason;
    return(fit);
  }
  vector<double> ulevels=vlevels_raw.levels;
  vector<double> vlevels=ulevels_raw.levels;
  vector<unsigned> uindex=vlevels_raw.index;
  vector<unsigned> vindex=ulevels_raw.index;
  fit.nu=unsigned(ulevels.size());
  fit.nv=unsigned(vlevels.size());
  fit.lu=double(fit.nu-1u)*dp;
  fit.lv=double(fit.nv-1u)*dp;
  if(fit.lu<fit.lv){
    ulevels=ulevels_raw.levels;
    vlevels=vlevels_raw.levels;
    uindex=ulevels_raw.index;
    vindex=vlevels_raw.index;
    fit.nu=unsigned(ulevels.size());
    fit.nv=unsigned(vlevels.size());
    fit.lu=double(fit.nu-1u)*dp;
    fit.lv=double(fit.nv-1u)*dp;
  }
  if(fit.nu<2 || fit.nv<2){
    fit.reason="rectangular fit needs at least 2x2 bond points";
    return(fit);
  }
  if(fit.nu*fit.nv!=n){
    fit.reason="bond points do not form a full rectangular grid";
    return(fit);
  }
  vector<int> occupancy(fit.nu*fit.nv,-1);
  for(unsigned c=0;c<n;c++){
    const unsigned iu=uindex[c];
    const unsigned iv=vindex[c];
    const unsigned idx=iv*fit.nu+iu;
    if(occupancy[idx]>=0){
      fit.reason="multiple bonds mapped to the same rectangular grid point";
      return(fit);
    }
    occupancy[idx]=int(c);
  }
  for(unsigned c=0;c<occupancy.size();c++){
    if(occupancy[c]<0){
      fit.reason="missing rectangular grid points on interface";
      return(fit);
    }
  }
  if(fit.lu<=ALMOSTZERO || fit.lv<=ALMOSTZERO){
    fit.reason="invalid rectangular side lengths";
    return(fit);
  }
  const double ucenter=0.5*(ulevels.front()+ulevels.back());
  const double vcenter=0.5*(vlevels.front()+vlevels.back());
  vector<double> ucoords(fit.nu,0.0), vcoords(fit.nv,0.0);
  for(unsigned c=0;c<fit.nu;c++)ucoords[c]=ulevels[c]-ucenter;
  for(unsigned c=0;c<fit.nv;c++)vcoords[c]=vlevels[c]-vcenter;

  // wu and wv are dimensionless one-dimensional area fractions. The tensor
  // product AFace_Bond*wu[i]*wv[j] gives the area represented by the bond at
  // rectangular grid location (i,j). The u weights target Ivv because Ivv is the
  // integral of u^2 dA; the v weights target Iuu because Iuu is the integral of
  // v^2 dA.
  vector<double> wu,wv;
  string whyu, whyv;
  fit.area_target=aface;
  if(iuu_target_user>0.0 && ivv_target_user>0.0){
    fit.iuu_target=iuu_target_user;
    fit.ivv_target=ivv_target_user;
    if(!SolveMomentWeights(ucoords,1.0,fit.ivv_target/aface,wu,whyu)){
      fit.reason="u weights failed for explicit Ivv target: "+whyu;
      return(fit);
    }
    if(!SolveMomentWeights(vcoords,1.0,fit.iuu_target/aface,wv,whyv)){
      fit.reason="v weights failed for explicit Iuu target: "+whyv;
      return(fit);
    }
  }
  else{
    fit.iuu_target=aface*fit.lv*fit.lv/12.0;
    fit.ivv_target=aface*fit.lu*fit.lu/12.0;
    if(!SolveMomentWeights(ucoords,fit.lu,fit.lu*fit.lu*fit.lu/12.0,wu,whyu)){
      fit.reason="u weights failed: "+whyu;
      return(fit);
    }
    if(!SolveMomentWeights(vcoords,fit.lv,fit.lv*fit.lv*fit.lv/12.0,wv,whyv)){
      fit.reason="v weights failed: "+whyv;
      return(fit);
    }
  }
  fit.iuv_target=0.0;
  for(unsigned c=0;c<n;c++){
    const unsigned iu=uindex[c];
    const unsigned iv=vindex[c];
    fit.u[c]=ucoords[iu];
    fit.v[c]=vcoords[iv];
    fit.aeff[c]=(iuu_target_user>0.0 && ivv_target_user>0.0? aface*wu[iu]*wv[iv]: wu[iu]*wv[iv]);
    fit.area_disc+=fit.aeff[c];
    fit.iuu_disc+=fit.aeff[c]*fit.v[c]*fit.v[c];
    fit.ivv_disc+=fit.aeff[c]*fit.u[c]*fit.u[c];
    fit.iuv_disc+=fit.aeff[c]*fit.u[c]*fit.v[c];
  }
  fit.ok=true;
  return(fit);
}

StBondAreaFit FitLineBondAreas2D(const vector<tdouble3>& mids,double aface,double dp,double ivv_target){
  StBondAreaFit fit;
  fit.line2d=true;
  const unsigned n=unsigned(mids.size());
  fit.u.assign(n,0.0);
  fit.v.assign(n,0.0);
  fit.aeff.assign(n,0.0);
  if(!n){
    fit.reason="no bond midpoints";
    return(fit);
  }
  if(aface<=ALMOSTZERO){
    fit.reason="AFace_Bond must be positive for 2D bonded line fitting";
    return(fit);
  }
  if(ivv_target<=ALMOSTZERO){
    fit.reason="Ivv_Bond must be positive for 2D bonded line fitting";
    return(fit);
  }
  if(dp<=ALMOSTZERO){
    fit.reason="dp must be positive for 2D bonded line fitting";
    return(fit);
  }

  double cmin[3]={mids[0].x,mids[0].y,mids[0].z};
  double cmax[3]={mids[0].x,mids[0].y,mids[0].z};
  for(unsigned c=1;c<n;c++){
    cmin[0]=min(cmin[0],mids[c].x); cmax[0]=max(cmax[0],mids[c].x);
    cmin[1]=min(cmin[1],mids[c].y); cmax[1]=max(cmax[1],mids[c].y);
    cmin[2]=min(cmin[2],mids[c].z); cmax[2]=max(cmax[2],mids[c].z);
  }
  const double span[3]={cmax[0]-cmin[0],cmax[1]-cmin[1],cmax[2]-cmin[2]};
  const unsigned outaxis=1;
  const unsigned aaxis=0;
  const unsigned baxis=2;
  const unsigned uaxis=(span[aaxis]>=span[baxis]? aaxis: baxis);
  const unsigned normalaxis=(span[aaxis]>=span[baxis]? baxis: aaxis);
  const unsigned vaxis=outaxis;
  const double maxline=max(span[uaxis],double(dp));
  const double tol=max(maxline*1e-6,max(double(dp)*1e-3,1e-12));
  if(span[outaxis]>tol){
    fit.reason=fun::PrintStr("2D bond midpoints left the x-z plane: y spread=%g exceeds tolerance=%g",span[outaxis],tol);
    return(fit);
  }
  if(span[normalaxis]>tol){
    fit.reason=fun::PrintStr("2D bond midpoints are not axis-aligned along %s: spread=%g exceeds tolerance=%g"
      ,AxisName(normalaxis),span[normalaxis],tol);
    return(fit);
  }
  if(span[uaxis]<=ALMOSTZERO){
    fit.reason="2D bonded line needs at least two distinct bond points along the interface";
    return(fit);
  }

  vector<double> uraw(n,0.0);
  for(unsigned c=0;c<n;c++)uraw[c]=Coord3(mids[c],uaxis);
  const StClusterLevels ulevels_raw=ClusterLevels(uraw,tol);
  if(!ulevels_raw.ok){
    fit.reason="u clustering failed: "+ulevels_raw.reason;
    return(fit);
  }
  fit.nu=unsigned(ulevels_raw.levels.size());
  fit.nv=1;
  fit.normalaxis=normalaxis;
  fit.uaxis=uaxis;
  fit.vaxis=vaxis;
  if(fit.nu<2){
    fit.reason="2D bonded line needs at least two u-levels";
    return(fit);
  }
  if(fit.nu!=n){
    fit.reason="2D bond points do not form a one-bond-per-level line";
    return(fit);
  }
  fit.lu=ulevels_raw.levels.back()-ulevels_raw.levels.front();
  fit.lv=0.0;
  if(fit.lu<=ALMOSTZERO){
    fit.reason="invalid 2D bonded line length";
    return(fit);
  }
  const double ucenter=0.5*(ulevels_raw.levels.front()+ulevels_raw.levels.back());
  vector<double> ucoords(fit.nu,0.0);
  for(unsigned c=0;c<fit.nu;c++)ucoords[c]=ulevels_raw.levels[c]-ucenter;

  vector<double> wu;
  string whyu;
  fit.area_target=aface;
  fit.iuu_target=0.0;
  fit.ivv_target=ivv_target;
  fit.iuv_target=0.0;
  if(!SolveMomentWeights(ucoords,aface,fit.ivv_target,wu,whyu)){
    fit.reason="u weights failed for explicit 2D Ivv target: "+whyu;
    return(fit);
  }
  for(unsigned c=0;c<n;c++){
    const unsigned iu=ulevels_raw.index[c];
    fit.u[c]=ucoords[iu];
    fit.v[c]=0.0;
    fit.aeff[c]=wu[iu];
    fit.area_disc+=fit.aeff[c];
    fit.iuu_disc+=0.0;
    fit.ivv_disc+=fit.aeff[c]*fit.u[c]*fit.u[c];
    fit.iuv_disc+=0.0;
  }
  fit.ok=true;
  return(fit);
}

}

//==============================================================================
/// Constructor.
//==============================================================================
JSphGpu::JSphGpu(bool withmpi):JSph(1,withmpi),DivAxis(MGDIV_None){
  ClassName="JSphGpu";
  CellDiv=NULL;
  Arrays_Cpu=NULL;
  Arrays_Gpu=NULL;
  GpuInfo=new JDsGpuInfo;
  Timersg=new JDsTimersGpu;
  InitVars();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphGpu::~JSphGpu(){
  DestructorActive=true;
  FreeCpuMemoryParticles();
  FreeGpuMemoryParticles();
  FreeGpuMemoryFixed();
  FreeCpuMemoryFixed();
  delete Arrays_Cpu; Arrays_Cpu=NULL;
  delete Arrays_Gpu; Arrays_Gpu=NULL;
  delete GpuInfo;    GpuInfo=NULL;
  delete Timersg;    Timersg=NULL;
  // cudaDeviceReset();
}

//==============================================================================
/// Throws exception related to a CUDA error from a static method.
//==============================================================================
void JSphGpu::RunExceptioonCudaStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,cudaError_t cuerr,std::string msg)
{
  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
  throw JException(srcfile,srcline,"JSphGpu",method,msg,"");
}

//==============================================================================
/// Checks CUDA error and throws exception from a static method.
//==============================================================================
void JSphGpu::CheckCudaErroorStatic(const std::string& srcfile,int srcline
  ,const std::string& method,std::string msg)
{
  const cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCudaStatic(srcfile,srcline,method,cuerr,msg);
}

//==============================================================================
/// Throws exception related to a CUDA error.
//==============================================================================
void JSphGpu::RunExceptioonCuda(const std::string& srcfile,int srcline
  ,const std::string& classname,const std::string& method
  ,cudaError_t cuerr,std::string msg)const
{
  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
  throw JException(srcfile,srcline,classname,method,msg,"");
}

//==============================================================================
/// Checks CUDA error and throws exception.
/// Comprueba error de CUDA y lanza excepcion si lo hubiera.
//==============================================================================
void JSphGpu::CheckCudaErroor(const std::string& srcfile,int srcline
  ,const std::string& classname,const std::string& method
  ,std::string msg)const
{
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCuda(srcfile,srcline,classname,method,cuerr,msg);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphGpu::InitVars(){
  RunMode="";
  memset(&BlockSizes,0,sizeof(StBlockSizes));
  BlockSizesStr="";

  DivData=DivDataGpuNull();

  Np=Npb=NpbOk=0;
  NpfPer=0;
  NpbPer=0;
  NpfPerM1=0;
  NpbPerM1=0;

  Idp_c=NULL;
  Code_c=NULL;
  Dcell_c=NULL;
  Posxy_c=NULL;
  Posz_c=NULL;
  Velrho_c=NULL;
  DemStress_c=NULL;
  DemStrain_c=NULL;
  DemShearStress_c=NULL;

  AuxPos_c=NULL;
  AuxVel_c=NULL;
  AuxRho_c=NULL;
  FreeCpuMemoryParticles();

  Idp_g=NULL;
  Code_g=NULL;
  Dcell_g=NULL;
  Posxy_g=NULL;
  Posz_g=NULL;
  PosCell_g=NULL;
  Velrho_g=NULL;
  DemStress_g=NULL;
  DemStrain_g=NULL;
  DemShearStress_g=NULL;
  DemBond_g=NULL;
  DemBondPairCount=0;
  DemBondRefCount=0;
  DemBondStartg=NULL;
  DemBondNgbg=NULL;
  DemBondL0g=NULL;
  DemBondPairIdg=NULL;
  DemBondOwnerg=NULL;
  DemBondStateg=NULL;
  DemBondUtg=NULL;
  DemBondLsegg=NULL;
  DemBondKng=NULL;
  DemBondAeffg=NULL;
  DemBondPairId1h.clear();
  DemBondPairId2h.clear();
  DemBondPairL0h.clear();
  DemBondPairIfcIdh.clear();
  DemBondPairUh.clear();
  DemBondPairVh.clear();
  DemBondPairNbondsh.clear();
  DemBondPairLsegh.clear();
  DemBondPairAfaceh.clear();
  DemBondPairAeffh.clear();
  DemBondPairKnh.clear();
  DemBondIfcBlock1h.clear();
  DemBondIfcBlock2h.clear();
  DemBondIfcNuh.clear();
  DemBondIfcNvh.clear();
  DemBondIfcLuh.clear();
  DemBondIfcLvh.clear();
  DemBondIfcAfaceh.clear();
  DemBondIfcIuuh.clear();
  DemBondIfcIvvh.clear();
  DemBondIfcIuvh.clear();

  PeriParent_g=NULL;

  BoundNor_g=NULL;     //-mDBC
  MotionVel_g=NULL;    //-mDBC2  //<vs_m2dbc>
  MotionAce_g=NULL;    //-mDBC2  //<vs_m2dbc>
  BoundMode_g=NULL;    //-mDBC2  //<vs_m2dbc>
  TangenVel_g=NULL;    //-mDBC2  //<vs_m2dbc>

  VelrhoM1_g=NULL;     //-Verlet
  PosxyPre_g=NULL;     //-Symplectic
  PoszPre_g=NULL;      //-Symplectic
  VelrhoPre_g=NULL;    //-Symplectic

  ViscDt_g=NULL; 
  Ace_g=NULL;
  Ar_g=NULL;
  Delta_g=NULL;
  NoPenShift_g = NULL; //<vs_m2dbcNP>

  SpsTauRho2_g=NULL;   //-Laminar+SPS.
  Sps2Strain_g=NULL;   //-Laminar+SPS.

  ShiftPosfs_g=NULL;   //-Shifting.

  ShiftVel_g=NULL;     //-ShiftingAdvanced //<vs_advshift>
  FSType_g=NULL;       //-ShiftingAdvanced //<vs_advshift>
  FSMinDist_g=NULL;    //-ShiftingAdvanced //<vs_advshift>
  FSNormal_g=NULL;     //-ShiftingAdvanced //<vs_advshift>

  //<vs_flexstruc_ini>
  FlexStrucDatag=NULL;
  FlexStrucRidpg=NULL;
  PosCell0g=NULL;
  NumPairsg=NULL;
  PairIdxBufferg=NULL;
  PairIdxg=NULL;
  KerCorrg=NULL;
  DefGradg=NULL;
  BoundNor0g=NULL;
  FlexStrucDtg=NULL;
  FlexStrucDtMax=0;
  //<vs_flexstruc_end>

  PsiClean_g=NULL;       //<vs_divclean>
  PsiCleanPre_g=NULL;    //<vs_divclean>
  PsiCleanRhs_g=NULL;    //<vs_divclean>
  CsPsiClean_g=NULL;     //<vs_divclean>
  CsPsiCleanMax=0;       //<vs_divclean>
  
  FreeGpuMemoryParticles();

  //-Variables for moving and floating bodies (GPU memory).
  RidpMotg=NULL;
  RidpBoundFloatg=NULL;
  FtoMasspg=NULL;
  FtoDatpg=NULL;
  FtoCenterg=NULL;
  FtoAceg=NULL;
  NStmFloatings=0;
  for(unsigned c=0;c<MaxNStmFloatings;c++)StmFloatings[c]=NULL;
  //-DEM.
  DemDatag=NULL;
  DemData2g=NULL;
  DemData3g=NULL;
  DemData4g=NULL;
  FreeGpuMemoryFixed();

  //-Variables for floating bodies (CPU memory).
  FtoCenterc=NULL;
  FreeCpuMemoryFixed();
}

//==============================================================================
/// Frees fixed memory on CPU for moving and floating bodies.
/// Libera memoria fija en CPU para moving y floating.
//==============================================================================
void JSphGpu::FreeCpuMemoryFixed(){
  MemCpuFixed=0;
  delete[] FtoCenterc; FtoCenterc=NULL;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocCpuMemoryFixed(){
  MemCpuFixed=0;
  try{
    //-Allocates memory for floating bodies.
    if(CaseNfloat){
      FtoCenterc=new tdouble3[FtCount];  MemCpuFixed+=(sizeof(tdouble3)*FtCount);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Frees fixed memory on the GPU for moving and floating bodies.
/// Libera memoria fija en GPU para moving y floating.
//==============================================================================
void JSphGpu::FreeGpuMemoryFixed(){
  MemGpuFixed=0;
  //-Memory for moving and floating particles.
  if(RidpMotg)  cudaFree(RidpMotg);    RidpMotg=NULL;
  if(RidpBoundFloatg) cudaFree(RidpBoundFloatg); RidpBoundFloatg=NULL;
  //-Memory for floating bodies.
  if(FtoMasspg) cudaFree(FtoMasspg);   FtoMasspg=NULL;
  if(FtoDatpg)  cudaFree(FtoDatpg);    FtoDatpg=NULL;
  if(FtoCenterg)cudaFree(FtoCenterg);  FtoCenterg=NULL;
  if(FtoAceg)   cudaFree(FtoAceg);     FtoAceg=NULL;
  //-Memory for DEM coefficients.
  if(DemDatag)  cudaFree(DemDatag);    DemDatag=NULL;
  if(DemData2g) cudaFree(DemData2g);   DemData2g=NULL;
  if(DemData3g) cudaFree(DemData3g);   DemData3g=NULL;
  if(DemData4g) cudaFree(DemData4g);   DemData4g=NULL;
  if(DemBondStartg)  cudaFree(DemBondStartg);  DemBondStartg=NULL;
  if(DemBondNgbg)    cudaFree(DemBondNgbg);    DemBondNgbg=NULL;
  if(DemBondL0g)     cudaFree(DemBondL0g);     DemBondL0g=NULL;
  if(DemBondPairIdg) cudaFree(DemBondPairIdg); DemBondPairIdg=NULL;
  if(DemBondOwnerg)  cudaFree(DemBondOwnerg);  DemBondOwnerg=NULL;
  if(DemBondStateg)  cudaFree(DemBondStateg);  DemBondStateg=NULL;
  if(DemBondUtg)     cudaFree(DemBondUtg);     DemBondUtg=NULL;
  if(DemBondLsegg)   cudaFree(DemBondLsegg);   DemBondLsegg=NULL;
  if(DemBondKng)     cudaFree(DemBondKng);     DemBondKng=NULL;
  if(DemBondAeffg)   cudaFree(DemBondAeffg);   DemBondAeffg=NULL;
  DemBondPairCount=0;
  DemBondRefCount=0;
  //<vs_flexstruc_ini>
  if(FlexStrucDatag)    cudaFree(FlexStrucDatag);     FlexStrucDatag=NULL;
  if(FlexStrucRidpg)    cudaFree(FlexStrucRidpg);     FlexStrucRidpg=NULL;
  if(PosCell0g)         cudaFree(PosCell0g);          PosCell0g=NULL;
  if(NumPairsg)         cudaFree(NumPairsg);          NumPairsg=NULL;
  if(PairIdxBufferg)    cudaFree(PairIdxBufferg);     PairIdxBufferg=NULL;
  if(PairIdxg)          cudaFree(PairIdxg);           PairIdxg=NULL;
  if(KerCorrg)          cudaFree(KerCorrg);           KerCorrg=NULL;
  if(DefGradg)          cudaFree(DefGradg);           DefGradg=NULL;
  if(BoundNor0g)        cudaFree(BoundNor0g);         BoundNor0g=NULL;
  if(FlexStrucDtg)      cudaFree(FlexStrucDtg);       FlexStrucDtg=NULL;
  //<vs_flexstruc_end>
  //-Frees streams for floating bodies.
  for(unsigned c=0;c<MaxNStmFloatings;c++){
    if(StmFloatings[c])cudaStreamDestroy(StmFloatings[c]);
    StmFloatings[c]=NULL;
  }
  NStmFloatings=0;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocGpuMemoryFixed(){
  MemGpuFixed=0;
  //-Allocates memory for moving and floating particles.
  if(CaseNmoving || CaseNfloat){
    MemGpuFixed+=fcuda::Malloc(&RidpMotg,CaseNmoving+CaseNfloat);
  }
  if(CaseNbound){
    MemGpuFixed+=fcuda::Malloc(&RidpBoundFloatg,CaseNbound);
  }
  //-Allocates memory for floating bodies.
  if(CaseNfloat){
    MemGpuFixed+=fcuda::Malloc(&FtoMasspg ,FtCount);
    MemGpuFixed+=fcuda::Malloc(&FtoDatpg  ,FtCount);
    MemGpuFixed+=fcuda::Malloc(&FtoCenterg,FtCount);
    MemGpuFixed+=fcuda::Malloc(&FtoAceg   ,FtCount*2);
  }
  //-Allocates streams for floating bodies.
  NStmFloatings=(FtCount>1? min(MaxNStmFloatings,FtCount): 0);
  for(unsigned c=0;c<NStmFloatings;c++){
    cudaStreamCreate(StmFloatings+c);
  }
  //-Allocates memory for DEM coefficients.
  if(UseDEM){
    MemGpuFixed+=fcuda::Malloc(&DemDatag,DemDataSize);
    MemGpuFixed+=fcuda::Malloc(&DemData2g,DemDataSize);
    MemGpuFixed+=fcuda::Malloc(&DemData3g,DemDataSize);
    MemGpuFixed+=fcuda::Malloc(&DemData4g,DemDataSize);
  }
  Check_CudaErroor("Failed GPU memory allocation.");
}

//==============================================================================
/// Frees CPU memory for the particles.
/// Libera memoria en CPU para particulas.
//==============================================================================
void JSphGpu::FreeCpuMemoryParticles(){
  //-Free array objects.
  delete Idp_c;      Idp_c=NULL;
  delete Code_c;     Code_c=NULL;
  delete Dcell_c;    Dcell_c=NULL;
  delete Posxy_c;    Posxy_c=NULL;
  delete Posz_c;     Posz_c=NULL;
  delete Velrho_c;   Velrho_c=NULL;
  delete DemStress_c;      DemStress_c=NULL;
  delete DemStrain_c;      DemStrain_c=NULL;
  delete DemShearStress_c; DemShearStress_c=NULL;

  delete AuxPos_c;   AuxPos_c=NULL;
  delete AuxVel_c;   AuxVel_c=NULL;
  delete AuxRho_c;   AuxRho_c=NULL;
  
  //-Free CPU memory for array objects.
  CpuParticlesSize=0;
  if(Arrays_Cpu)Arrays_Cpu->Reset();
}

//==============================================================================
/// Allocates memory for the main particle data.
/// Reserva memoria para datos principales de particulas.
//==============================================================================
void JSphGpu::AllocCpuMemoryParticles(unsigned np){
  FreeCpuMemoryParticles();
  //-Calculate number of partices to allocate memory.
  if(np>0)np=np+PARTICLES_OVERMEMORY_MIN;
  CpuParticlesSize=np;
  //-Set size of arrays.
  Arrays_Cpu->SetArraySize(CpuParticlesSize);
  //-Create arrays for basic particle data.
  Idp_c   =new acuint    ("Idpc"   ,Arrays_Cpu,true);
  Code_c  =new actypecode("Codec"  ,Arrays_Cpu,true);
  Dcell_c =new acuint    ("Dcellc" ,Arrays_Cpu,true);
  Posxy_c =new acdouble2 ("Posxyc" ,Arrays_Cpu,true);
  Posz_c  =new acdouble  ("Poszc"  ,Arrays_Cpu,true);
  Velrho_c=new acfloat4  ("Velrhoc",Arrays_Cpu,true);
  if(CaseNfloat)DemStress_c=new acfloat("DemStressc",Arrays_Cpu,true);
  if(CaseNfloat)DemStrain_c=new acfloat("DemStrainc",Arrays_Cpu,true);
  if(CaseNfloat)DemShearStress_c=new acfloat("DemShearStressc",Arrays_Cpu,true);

  AuxPos_c=new acdouble3 ("AuxPosc",Arrays_Cpu,true);
  AuxVel_c=new acfloat3  ("AuxVelc",Arrays_Cpu,true);
  AuxRho_c=new acfloat   ("AuxRhoc",Arrays_Cpu,true);
}

//==============================================================================
/// Frees GPU memory for the particles.
/// Libera memoria en GPU para particulas.
//==============================================================================
void JSphGpu::FreeGpuMemoryParticles(){
  //-Free array objects.
  delete Idp_g;         Idp_g=NULL;
  delete Code_g;        Code_g=NULL;
  delete Dcell_g;       Dcell_g=NULL;
  delete Posxy_g;       Posxy_g=NULL;
  delete Posz_g;        Posz_g=NULL;
  delete PosCell_g;     PosCell_g=NULL;
  delete Velrho_g;      Velrho_g=NULL;
  delete DemStress_g;      DemStress_g=NULL;
  delete DemStrain_g;      DemStrain_g=NULL;
  delete DemShearStress_g; DemShearStress_g=NULL;
  delete DemBond_g;        DemBond_g=NULL;

  delete PeriParent_g;  PeriParent_g=NULL;

  delete BoundNor_g;    BoundNor_g=NULL;    //-mDBC
  delete MotionVel_g;   MotionVel_g=NULL;   //-mDBC2  //<vs_m2dbc>
  delete MotionAce_g;   MotionAce_g=NULL;   //-mDBC2  //<vs_m2dbc>
  delete BoundMode_g;   BoundMode_g=NULL;   //-mDBC2  //<vs_m2dbc>
  delete TangenVel_g;   TangenVel_g=NULL;   //-mDBC2  //<vs_m2dbc>

  delete VelrhoM1_g;    VelrhoM1_g=NULL;    //-Verlet
  delete PosxyPre_g;    PosxyPre_g=NULL;    //-Symplectic
  delete PoszPre_g;     PoszPre_g=NULL;     //-Symplectic
  delete VelrhoPre_g;   VelrhoPre_g=NULL;   //-Symplectic

  delete ViscDt_g;      ViscDt_g=NULL; 
  delete Ace_g;         Ace_g=NULL;
  delete Ar_g;          Ar_g=NULL;
  delete Delta_g;       Delta_g=NULL;
  delete NoPenShift_g;  NoPenShift_g=NULL;  //-NoPenetration //<vs_m2dbcNP>

  delete SpsTauRho2_g;  SpsTauRho2_g=NULL;  //-Laminar+SPS.
  delete Sps2Strain_g;  Sps2Strain_g=NULL;  //-Laminar+SPS.

  delete ShiftPosfs_g;  ShiftPosfs_g=NULL;  //-Shifting.

  delete ShiftVel_g;    ShiftVel_g=NULL;    //-ShiftingAdvanced //<vs_advshift>
  delete FSType_g;      FSType_g=NULL;      //-ShiftingAdvanced //<vs_advshift>
  delete FSMinDist_g;   FSMinDist_g=NULL;   //-ShiftingAdvanced //<vs_advshift>
  delete FSNormal_g;    FSNormal_g=NULL;    //-ShiftingAdvanced //<vs_advshift>

  delete PsiClean_g;    PsiClean_g=NULL;    //<vs_divclean>
  delete PsiCleanPre_g; PsiCleanPre_g=NULL; //<vs_divclean>
  delete PsiCleanRhs_g; PsiCleanRhs_g=NULL; //<vs_divclean>
  delete CsPsiClean_g;  CsPsiClean_g=NULL;  //<vs_divclean>

  //-Free GPU memory for array objects.
  GpuParticlesSize=0;
  if(Arrays_Gpu)Arrays_Gpu->Reset();
}

//==============================================================================
/// Allocates GPU memory for the particles.
/// Reserva memoria en Gpu para las particulas. 
//==============================================================================
void JSphGpu::AllocGpuMemoryParticles(unsigned np){
  FreeGpuMemoryParticles();
  //-Calculate number of partices to allocate memory.
  GpuParticlesSize=np+PARTICLES_OVERMEMORY_MIN;
  //-Checks maximum number.
  if(GpuParticlesSize>=unsigned(INT_MAX))
    Run_Exceptioon("Particle number for allocation is too large.");
  //-Set size of arrays.
  Arrays_Gpu->SetArraySize(GpuParticlesSize);
  //-Create arrays for basic particle data.
  Idp_g    =new aguint    ("Idpg"    ,Arrays_Gpu,true);
  Code_g   =new agtypecode("Codeg"   ,Arrays_Gpu,true);
  Dcell_g  =new aguint    ("Dcellg"  ,Arrays_Gpu,true);
  Posxy_g  =new agdouble2 ("Posxyg"  ,Arrays_Gpu,true);
  Posz_g   =new agdouble  ("Poszg"   ,Arrays_Gpu,true);
  PosCell_g=new agfloat4  ("PosCellg",Arrays_Gpu,true);
  Velrho_g =new agfloat4  ("Velrhog" ,Arrays_Gpu,true);
  if(CaseNfloat)DemStress_g=new agfloat("DemStressg",Arrays_Gpu,true);
  if(CaseNfloat)DemStrain_g=new agfloat("DemStraing",Arrays_Gpu,true);
  if(CaseNfloat)DemShearStress_g=new agfloat("DemShearStressg",Arrays_Gpu,true);
  if(CaseNfloat)DemBond_g  =new agfloat("DemBondg"  ,Arrays_Gpu,true);
  //-Arrays for mDBC.
  if(UseNormals){
    BoundNor_g=new agfloat3("BoundNorg",Arrays_Gpu,true);
    if(SlipMode>=SLIP_NoSlip){ //<vs_m2dbc_ini>
      MotionVel_g=new agfloat3("MotionVelg",Arrays_Gpu,true);
      MotionAce_g=new agfloat3("MotionAceg",Arrays_Gpu,true);
      BoundMode_g=new agbyte  ("BoundModeg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
      TangenVel_g=new agfloat3("TangenVelg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
    } //<vs_m2dbc_end>
  }
  //-Arrays for Verlet.
  if(TStep==STEP_Verlet){
    VelrhoM1_g=new agfloat4  ("VelrhoM1g",Arrays_Gpu,true);
  }
  //-Arrays for Symplectic.
  if(TStep==STEP_Symplectic){
    PosxyPre_g =new agdouble2("PosxyPreg" ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
    PoszPre_g  =new agdouble ("PoszPreg"  ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
    VelrhoPre_g=new agfloat4 ("VelrhoPreg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  }
  //-Arrays for forces computation.
  ViscDt_g    =new agfloat ("ViscDtg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  Ace_g       =new agfloat3("Aceg"   ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
  Ar_g        =new agfloat ("Arg"    ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
  Delta_g     =new agfloat ("Deltag" ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
  //-Arrays for No Pentration.
  if(TMdbc2==MDBC2_NoPen){
    NoPenShift_g = new agfloat4("NoPenShiftg", Arrays_Gpu, false); //-NO INITIAL MEMORY.
  }
  //-Arrays for Laminar+SPS.
  if(TVisco==VISCO_LaminarSPS){
    SpsTauRho2_g=new agsymatrix3f("SpsTauRho2g",Arrays_Gpu,true);
    Sps2Strain_g=new agsymatrix3f("Sps2Straing",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  }
  //-Arrays for Shifting.
  ShiftPosfs_g=new agfloat4("ShiftPosfsg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  //-Arrays for Advanced shifting. //<vs_advshift_ini>
  if(ShiftingAdv!=NULL){
    ShiftVel_g    =new agfloat4 ("ShiftVel_g" ,Arrays_Gpu,true); 
    FSType_g      =new aguint   ("FSType_g"   ,Arrays_Gpu,true);
    FSNormal_g    =new agfloat3 ("FSNormal_g" ,Arrays_Gpu,false);  //-NO INITIAL MEMORY.
    FSMinDist_g   =new agfloat  ("FSMinDist_g",Arrays_Gpu,false);  //-NO INITIAL MEMORY.
    if(PeriActive && !PeriParent_g){
      PeriParent_g=new aguint(  "PeriParentg" ,Arrays_Gpu,true);
    }
  }//<vs_advshift_end>

  //<vs_divclean_ini>
  if(DivCleaning){
    PsiClean_g    =new agfloat("PhiClean_g"   ,Arrays_Gpu,true);
    PsiCleanPre_g =new agfloat("PhiCleanPre_g",Arrays_Gpu,false);
    PsiCleanRhs_g =new agfloat("PhiCleanRhs_g",Arrays_Gpu,false);
    CsPsiClean_g  =new agfloat("CsPhiClean_g" ,Arrays_Gpu,false); 
  }
  //<vs_divclean_end>

  //-Check CUDA errors.
  Check_CudaErroor("Failed GPU memory allocation.");
}

//==============================================================================
/// Resizes GPU memory for particles saving current data (ndatagpu).
//==============================================================================
void JSphGpu::ResizeGpuMemoryParticlesData(unsigned ndatagpu
  ,unsigned npnew,unsigned npmin)
{
  npnew=npnew+PARTICLES_OVERMEMORY_MIN;
  //-Shows GPU memory allocation.
  const llong memgpuparticles=Arrays_Gpu->GetAllocMemoryGpu();
  const double mbparticle=(double(memgpuparticles)/MEBIBYTE)/GpuParticlesSize; //-MiB per particle.
  const string txover=(npmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(npnew)/npmin): "");
  Log->Printf("**JSphGpu: Requesting GPU memory for %s particles%s: %.1f MiB (%u times)."
    ,KINT(npnew),txover.c_str(),mbparticle*npnew,Arrays_Gpu->GetCountResizeDataCalls()+1);
  //-Resizes GPU memory allocation.
  Arrays_Gpu->SetDataArraySize(npnew,ndatagpu);
  //-Updates values.
  GpuParticlesSize=npnew;
}

//==============================================================================
/// Returns the allocated memory on the CPU.
/// Devuelve la memoria reservada en CPU.
//==============================================================================
llong JSphGpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //-Allocated for particle arrays on CPU.
  if(Arrays_Cpu)s+=Arrays_Cpu->GetAllocMemoryCpu();
  //-Allocated in AllocCpuMemoryFixed().
  s+=MemCpuFixed;
  //-Allocated in other objects.
  if(MLPistons)s+=MLPistons->GetAllocMemoryCpu();
  return(s);
}

//==============================================================================
/// Returns the allocated memory in the GPU.
/// Devuelve la memoria reservada en GPU.
//==============================================================================
llong JSphGpu::GetAllocMemoryGpu()const{  
  llong s=0;
  //-Allocated for particle arrays on GPU.
  if(Arrays_Gpu)s+=Arrays_Gpu->GetAllocMemoryGpu();
  //-Allocated in AllocGpuMemoryFixed().
  s+=MemGpuFixed;
  //-Allocated in their objects.
  if(MLPistons)s+=MLPistons->GetAllocMemoryGpu();
  return(s);
}

//==============================================================================
/// Displays the allocated memory.
/// Visualiza la memoria reservada.
//==============================================================================
void JSphGpu::PrintAllocMemory(llong mcpu,llong mgpu)const{
  if(!Nstep)Log->Print("");
  Log->Printf("%s allocated memory in CPU: %s (%.2f MiB)  (%u particle arrays)"
    ,(!Nstep? "Initial": "**Updated"),KINT(mcpu),double(mcpu)/MEBIBYTE,Arrays_Cpu->GetArrayCount());
  Log->Printf("%s allocated memory in GPU: %s (%.2f MiB)  (%u particle arrays)"
    ,(!Nstep? "Initial": "**Updated"),KINT(mgpu),double(mgpu)/MEBIBYTE,Arrays_Gpu->GetArrayCount());
  Arrays_Cpu->GetArrayCountUpdated();
  Arrays_Gpu->GetArrayCountUpdated();
}

//==============================================================================
/// Loads constant data for GPU in ctes.
/// Carga constantes para GPU en ctes.
//==============================================================================
void JSphGpu::GetConstantData(StCteInteraction& ctes)const{
  memset(&ctes,0,sizeof(StCteInteraction));
  //-Set mass particle values.
  SetCtegMass(ctes,MassBound,MassFluid);
  //-Set distance values depending on h and dp.
  SetCtegKsize(ctes,KernelH,KernelSize,KernelSize2,PosCellSize
    ,Eta2,float(Dp),Scell,MovLimit);
  //-Wendland constants are always computed since this kernel is used in some parts where other kernels are not defined (e.g. mDBC, inlet/outlet...).
  SetCtegKerWendland(ctes,KWend);
  //-Copies constants for other kernels.
  if(TKernel==KERNEL_Cubic)SetCtegKerCubic(ctes,KCubic);
  //-Set values for density and pressure.
  SetCtegRho(ctes,RhopZero,1.f/RhopZero,Gamma,float(Cs0),CteB);
  //-Set values for DDT.
  SetCtegDdt(ctes,DDTkh,DDTgz);
  //-Set values on formulation options.
  SetCtegOpts(ctes,unsigned(TBoundary));
  //-Set values for open periodic boundaries.
  SetCtegPeriodic(ctes,PeriActive,PeriXinc,PeriYinc,PeriZinc);
  //-Set values for map definition.
  SetCtegMap(ctes,MapRealPosMin,MapRealSize);
  //-Set values depending on the assigned domain (can change).
  //-MGDIV_Z is necessary to avoid errors in KerGetInteraction_Cells().
  SetCtegDomain(ctes,MGDIV_Z,DomCellCode,DomPosMin);
}

//==============================================================================
/// Uploads constants to the GPU.
/// Copia constantes a la GPU.
//==============================================================================
void JSphGpu::ConstantDataUp(){
  StCteInteraction ctes;
  GetConstantData(ctes);
  //-Copy constants to GPU memory.
  cusph::CteInteractionUp(&ctes);
  Check_CudaErroor("Failed copying constants to GPU.");
}

//==============================================================================
/// Uploads particle data to the GPU.
/// Sube datos de particulas a la GPU.
//==============================================================================
void JSphGpu::ParticlesDataUp(unsigned n,const tfloat3* boundnor){
  Idp_g   ->CuCopyFromHost(Idp_c,n);
  Code_g  ->CuCopyFromHost(Code_c,n);
  Dcell_g ->CuCopyFromHost(Dcell_c,n);
  Posxy_g ->CuCopyFromHost(Posxy_c,n);
  Posz_g  ->CuCopyFromHost(Posz_c,n);
  Velrho_g->CuCopyFromHost(Velrho_c,n);
  if(UseNormals){
    BoundNor_g->CuCopyFromHost2(boundnor,n);
  }
  Check_CudaErroor("Failed copying data to GPU.");
}

//==============================================================================
/// Recovers particle data from the GPU and returns the particle number that
/// are less than n if the paeriodic particles are removed.
/// - code: Recovers data of Codeg.
/// - onlynormal: Only retains the normal particles, removes the periodic ones.
///
/// Recupera datos de particulas de la GPU y devuelve el numero de particulas que
/// sera menor que n si se eliminaron las periodicas.
/// - code: Recupera datos de Codeg.
/// - onlynormal: Solo se queda con las normales, elimina las particulas periodicas.
//==============================================================================
unsigned JSphGpu::ParticlesDataDown(unsigned n,unsigned pini,bool code
  ,bool onlynormal,const byte* filterg,unsigned& npfilterdel,float* demstress,float* demstrain,float* demshearstress,float* dembond)
{
  unsigned num=n;
  //-Copy data from GPU memory.
  Idp_g   ->CuCopyToHostOffset(pini,Idp_c  ,0,n);
  Posxy_g ->CuCopyToHostOffset(pini,Posxy_c,0,n);
  Posz_g  ->CuCopyToHostOffset(pini,Posz_c,0,n);
  Velrho_g->CuCopyToHostOffset(pini,Velrho_c,0,n);
  if(demstress && DemStress_g){
    DemStress_g->CuCopyToHostOffset(pini,DemStress_c,0,n);
  }
  else if(demstress){
    for(unsigned p=0;p<n;p++)demstress[p]=0;
  }
  if(demstrain && DemStrain_g){
    DemStrain_g->CuCopyToHostOffset(pini,DemStrain_c,0,n);
  }
  else if(demstrain){
    for(unsigned p=0;p<n;p++)demstrain[p]=0;
  }
  if(demshearstress && DemShearStress_g){
    DemShearStress_g->CuCopyToHostOffset(pini,DemShearStress_c,0,n);
  }
  else if(demshearstress){
    for(unsigned p=0;p<n;p++)demshearstress[p]=0;
  }
  if(dembond && DemBond_g)cudaMemcpy(dembond,DemBond_g->cptr()+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  else if(dembond){
    for(unsigned p=0;p<n;p++)dembond[p]=0;
  }
  if(code || onlynormal){
    Code_g->CuCopyToHostOffset(pini,Code_c,0,n);
  }
  //-Creates simple CPU pointers.
  unsigned* idpc   =Idp_c->ptr(); 
  tdouble2* posxyc =Posxy_c->ptr(); 
  double*   poszc  =Posz_c->ptr(); 
  tfloat4*  velrhoc=Velrho_c->ptr();
  typecode* codec  =Code_c->ptr();
  float* demstressc=((demstress && DemStress_c)? DemStress_c->ptr(): NULL);
  float* demstrainc=((demstrain && DemStrain_c)? DemStrain_c->ptr(): NULL);
  float* demshearc=((demshearstress && DemShearStress_c)? DemShearStress_c->ptr(): NULL);

  //-Obtain filter data on CPU memory. //<vs_outpaarts_ini>
  byte* filter=NULL; 
  if(filterg){
    filter=(byte*)AuxRho_c->ptr();
    cudaMemcpy(filter,filterg+pini,sizeof(byte)*n,cudaMemcpyDeviceToHost);
  }//<vs_outpaarts_end>
  Check_CudaErroor("Failed copying data from GPU.");
  
  //-Eliminate non-normal particles (periodic & others).
  const bool usefilter=(filter!=NULL);
  unsigned nfilter=0;
  if(onlynormal || usefilter){
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      const bool isnormal=(!onlynormal || CODE_IsNormal(codec[p]));
      const bool selected=(isnormal && (!usefilter || (filter[p]&1)));
      if(ndel && selected){
        const unsigned p2=p-ndel;
        idpc   [p2]=idpc   [p];
        posxyc [p2]=posxyc [p];
        poszc  [p2]=poszc  [p];
        velrhoc[p2]=velrhoc[p];
        if(demstressc)demstressc[p2]=demstressc[p];
        if(demstrainc)demstrainc[p2]=demstrainc[p];
        if(demshearc)demshearc[p2]=demshearc[p];
        if(dembond)dembond[p2]=dembond[p];
        codec  [p2]=codec  [p];
      }
      if(!selected){
        ndel++;
        if(isnormal)nfilter++;//-Normal particles removed by filters.
      }
    }
    num-=ndel;
  }
  npfilterdel=nfilter;
  filter=NULL;
  if(demstressc){
    for(unsigned p=0;p<num;p++)demstress[p]=demstressc[p];
  }
  if(demstrainc){
    for(unsigned p=0;p<num;p++)demstrain[p]=demstrainc[p];
  }
  if(demshearc){
    for(unsigned p=0;p<num;p++)demshearstress[p]=demshearc[p];
  }
  //-Converts data to a simple format in AuxPos_c ,AuxVel_c and AuxRhop_c.
  Pos21Vel4ToPos3Vel31(num,posxyc,poszc,velrhoc
    ,AuxPos_c->ptr(),AuxVel_c->ptr(),AuxRho_c->ptr());
  return(num);
}

//==============================================================================
/// Initialises CUDA device and returns GPU id.
/// Inicializa dispositivo CUDA y devuelve id de GPU.
//==============================================================================
int JSphGpu::SelecDevice(int gpuid){
  //-Shows information on available GPUs.
  JDsGpuInfo::ShowGpusInfo(Log);
  //-GPU selection.
  Log->Print("[GPU Hardware]");
  GpuInfo->SelectGpu(gpuid);
  GpuInfo->ShowSelectGpusInfo(Log);
  cudaSetDevice(GpuInfo->GetGpuId());
  return(GpuInfo->GetGpuId());
}

//==============================================================================
/// Configures BlockSize for main interaction CUDA kernels.
//==============================================================================
void JSphGpu::ConfigBlockSizes(bool usezone,bool useperi){
  Log->Print(" ");
  BlockSizesStr="";
  if(CellMode==CELLMODE_Full || CellMode==CELLMODE_Half){
    BlockSizes.forcesbound=BlockSizes.forcesfluid=BlockSizes.forcesdem=BSIZE_FORCES;
    //-Collects kernel information.
    StKerInfo kerinfo;
    memset(&kerinfo,0,sizeof(StKerInfo));
    StDivDataGpu divdatag;
    memset(&divdatag,0,sizeof(StDivDataGpu));
    #ifndef DISABLE_BSMODES
      const StInterParmsg parms=StrInterParmsg(Simulate2D
        ,TKernel,FtMode
        ,TVisco,TDensity,ShiftingMode,TMdbc2 //<vs_m2dbc>
        ,false,false,false,false //<vs_advshift>
        ,0,0,0,0,100,0,0
        ,0,0,divdatag,NULL
        ,NULL,NULL,NULL
        ,NULL,NULL,NULL
        ,NULL,NULL,NULL //<vs_m2dbc>
        ,NULL //<vs_m2dbc> //SHABA4
        ,NULL,NULL,NULL
        ,NULL,NULL,NULL,NULL
        ,NULL
        ,NULL
        ,NULL // SHABA
        ,NULL,NULL     //<vs_advshift>
        ,NULL,NULL,NULL,0,false //<vs_divclean>
        ,NULL,&kerinfo);
      cusph::Interaction_Forces(parms);
      if(UseDEM)cusph::Interaction_ForcesDem(BlockSizes.forcesdem,CaseNfloat,divdatag
        ,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,0
        ,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL
        ,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL
        ,&kerinfo);
    #endif
    //Log->Printf("====> bound -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesbound_rg,kerinfo.forcesbound_bs,kerinfo.forcesbound_bsmax);
    //Log->Printf("====> fluid -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesfluid_rg,kerinfo.forcesfluid_bs,kerinfo.forcesfluid_bsmax);
    //Log->Printf("====> dem   -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesdem_rg  ,kerinfo.forcesdem_bs  ,kerinfo.forcesdem_bsmax);
    Log->Printf("BlockSize calculation mode: Fixed.");
    const string txb=fun::PrintStr("BsForcesBound=%d ",BlockSizes.forcesbound)+(kerinfo.forcesbound_rg? fun::PrintStr("(%d regs)",kerinfo.forcesbound_rg): string("(? regs)"));
    const string txf=fun::PrintStr("BsForcesFluid=%d ",BlockSizes.forcesfluid)+(kerinfo.forcesfluid_rg? fun::PrintStr("(%d regs)",kerinfo.forcesfluid_rg): string("(? regs)"));
    const string txd=fun::PrintStr("BsForcesDem=%d "  ,BlockSizes.forcesdem)  +(kerinfo.forcesdem_rg  ? fun::PrintStr("(%d regs)",kerinfo.forcesdem_rg  ): string("(? regs)"));
    Log->Print(string("  ")+txb);
    Log->Print(string("  ")+txf);
    if(UseDEM)Log->Print(string("  ")+txd);
    if(!BlockSizesStr.empty())BlockSizesStr=BlockSizesStr+" - ";
    BlockSizesStr=BlockSizesStr+txb+" - "+txf;
    if(UseDEM)BlockSizesStr=BlockSizesStr+" - "+txd;
  }
  else Run_Exceptioon("CellMode unrecognised.");
  Log->Print(" ");
}

//==============================================================================
/// Configures execution mode in the GPU.
/// Configura modo de ejecucion en GPU.
//==============================================================================
void JSphGpu::ConfigRunMode(){
  Hardware=GpuInfo->GetHardware();
  //-Defines RunMode.
  RunMode="";
  if(Stable)RunMode=RunMode+(!RunMode.empty()? " - ": "")+"Stable";
  RunMode=RunMode+(!RunMode.empty()? " - ": "")+"Pos-Cell";
  RunMode=RunMode+(!RunMode.empty()? " - ": "")+"Single-GPU";
  //-Shows RunMode.
  Log->Print(" ");
  Log->Print(fun::VarStr("RunMode",RunMode));
  Log->Print(" ");
}

//==============================================================================
/// Adjusts variables of particles of floating bodies for GPU execution.
//==============================================================================
void JSphGpu::InitFloatingsGpu(float* ftomasspg,float4* ftodatag
  ,double3* ftocenterg,float4* demdatag,float4* demdata2g,float* demdata3g)const
{
  //-Copies massp values to GPU.
  {
    float* massp=new float[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++)massp[cf]=FtObjs[cf].massp;
    cudaMemcpy(ftomasspg,massp,sizeof(float)*FtCount,cudaMemcpyHostToDevice);
    delete[] massp; massp=NULL;
  }
  //-Copies floating values to GPU.
  {
    typedef struct{
      unsigned pini;
      unsigned np;
      float radius;
      float massp;
    }stdata;
    stdata* datp=new stdata[FtCount];
    tdouble3* centerc=new tdouble3[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++){
      const StFloatingData& fobj=FtObjs[cf];
      datp[cf].pini=fobj.begin-CaseNfixed;
      datp[cf].np=fobj.count;
      datp[cf].radius=fobj.radius;
      datp[cf].massp=fobj.massp;
      centerc[cf]=fobj.center;
    }
    cudaMemcpy(ftodatag,datp,sizeof(float4)*FtCount,cudaMemcpyHostToDevice);
    cudaMemcpy(ftocenterg,centerc,sizeof(double3)*FtCount,cudaMemcpyHostToDevice);
    delete[] datp;    datp=NULL;
    delete[] centerc; centerc=NULL;
  }
  //-Copies data object for DEM to GPU.
  if(UseDEM){
    float4* ddata=new float4[DemDataSize];
    float4* ddata2=new float4[DemDataSize];
    float* ddata3=new float[DemDataSize];
    float2* ddata4=new float2[DemDataSize];
    const float dp=float(Dp);
    unsigned active_dem_slots=0;
    for(unsigned c=0;c<DemDataSize;c++){
      ddata[c].x=DemData[c].mass;
      ddata[c].y=DemData[c].tau;
      ddata[c].z=DemData[c].kfric;
      ddata[c].w=DemData[c].restitu;
      ddata2[c].x=DemData[c].young;
      ddata2[c].y=DemData[c].poisson;
      ddata2[c].z=DemData[c].sigma_c;
      ddata2[c].w=DemData[c].ksbond;
      ddata3[c]=DemData[c].cnbond;
      ddata4[c].x=(DemData[c].bondfailmode==2? DemData[c].eps_tension_c: DemData[c].sigma_c);
      ddata4[c].y=float(DemData[c].bondfailmode);

      // Only report populated DEM material slots. Unused slots are left zeroed.
      if(ddata2[c].x>0.f || ddata4[c].y>0.f){
        active_dem_slots++;
      }
    }
    Log->Printf("DEM debug (GPU upload): active material slots=%u, bonded normal inputs loaded from material file, dp=%g"
      ,active_dem_slots,dp);
    for(unsigned c=0;c<DemDataSize;c++){
      if(ddata2[c].x>0.f || ddata4[c].y>0.f){
        const float lsegbond=(DemData[c].lsegbond>0.f? DemData[c].lsegbond: 1.f);
        const float afacebond=(DemData[c].afacebond>0.f? DemData[c].afacebond: 1.f);
        const float iuubond=(DemData[c].iuubond>0.f? DemData[c].iuubond: -1.f);
        const float ivvbond=(DemData[c].ivvbond>0.f? DemData[c].ivvbond: -1.f);
        const char* failmode=(DemData[c].bondfailmode==2? "strain": (DemData[c].bondfailmode==1? "stress": "none"));
        if(ddata2[c].w>=0.f){
          if(iuubond>0.f && ivvbond>0.f)Log->Printf("  DEM material[%u] young=%g poisson=%g failure=%s threshold=%g lseg_bond=%g aface_bond=%g iuu_bond=%g ivv_bond=%g ks_bond_explicit=%g cn_bond=%g"
            ,c,ddata2[c].x,ddata2[c].y,failmode,ddata4[c].x,lsegbond,afacebond,iuubond,ivvbond,ddata2[c].w,ddata3[c]);
          else Log->Printf("  DEM material[%u] young=%g poisson=%g failure=%s threshold=%g lseg_bond=%g aface_bond=%g ks_bond_explicit=%g cn_bond=%g"
            ,c,ddata2[c].x,ddata2[c].y,failmode,ddata4[c].x,lsegbond,afacebond,ddata2[c].w,ddata3[c]);
        }
        else{
          if(iuubond>0.f && ivvbond>0.f)Log->Printf("  DEM material[%u] young=%g poisson=%g failure=%s threshold=%g lseg_bond=%g aface_bond=%g iuu_bond=%g ivv_bond=%g ks_bond_default=from_kn_pair cn_bond=%g"
            ,c,ddata2[c].x,ddata2[c].y,failmode,ddata4[c].x,lsegbond,afacebond,iuubond,ivvbond,ddata3[c]);
          else Log->Printf("  DEM material[%u] young=%g poisson=%g failure=%s threshold=%g lseg_bond=%g aface_bond=%g ks_bond_default=from_kn_pair cn_bond=%g"
            ,c,ddata2[c].x,ddata2[c].y,failmode,ddata4[c].x,lsegbond,afacebond,ddata3[c]);
        }
      }
    }
    cudaMemcpy(demdatag,ddata,sizeof(float4)*DemDataSize,cudaMemcpyHostToDevice);
    cudaMemcpy(demdata2g,ddata2,sizeof(float4)*DemDataSize,cudaMemcpyHostToDevice);
    cudaMemcpy(demdata3g,ddata3,sizeof(float)*DemDataSize,cudaMemcpyHostToDevice);
    cudaMemcpy(DemData4g,ddata4,sizeof(float2)*DemDataSize,cudaMemcpyHostToDevice);
    delete[] ddata; ddata=NULL;
    delete[] ddata2; ddata2=NULL;
    delete[] ddata3; ddata3=NULL;
    delete[] ddata4; ddata4=NULL;
  }
}

//==============================================================================
/// Builds the persistent DEM bond graph from the initial floating-particle geometry.
//==============================================================================
void JSphGpu::InitDemBondsGpu(){
  // Reset the persistent bond graph. Bonds are created once from the initial
  // boundary/floating geometry and then only their state/history evolves.
  DemBondPairCount=0;
  DemBondRefCount=0;
  if(DemBondStartg){ cudaFree(DemBondStartg); DemBondStartg=NULL; }
  if(DemBondNgbg){ cudaFree(DemBondNgbg); DemBondNgbg=NULL; }
  if(DemBondL0g){ cudaFree(DemBondL0g); DemBondL0g=NULL; }
  if(DemBondPairIdg){ cudaFree(DemBondPairIdg); DemBondPairIdg=NULL; }
  if(DemBondOwnerg){ cudaFree(DemBondOwnerg); DemBondOwnerg=NULL; }
  if(DemBondStateg){ cudaFree(DemBondStateg); DemBondStateg=NULL; }
  if(DemBondUtg){ cudaFree(DemBondUtg); DemBondUtg=NULL; }
  if(DemBondLsegg){ cudaFree(DemBondLsegg); DemBondLsegg=NULL; }
  if(DemBondKng){ cudaFree(DemBondKng); DemBondKng=NULL; }
  if(DemBondAeffg){ cudaFree(DemBondAeffg); DemBondAeffg=NULL; }
  DemBondPairId1h.clear();
  DemBondPairId2h.clear();
  DemBondPairL0h.clear();
  DemBondPairIfcIdh.clear();
  DemBondPairUh.clear();
  DemBondPairVh.clear();
  DemBondPairNbondsh.clear();
  DemBondPairLsegh.clear();
  DemBondPairAfaceh.clear();
  DemBondPairAeffh.clear();
  DemBondPairKnh.clear();
  DemBondIfcBlock1h.clear();
  DemBondIfcBlock2h.clear();
  DemBondIfcNuh.clear();
  DemBondIfcNvh.clear();
  DemBondIfcLuh.clear();
  DemBondIfcLvh.clear();
  DemBondIfcAfaceh.clear();
  DemBondIfcIuuh.clear();
  DemBondIfcIvvh.clear();
  DemBondIfcIuvh.clear();
  if(!UseDEM || !CaseNfloat)return;

  // Floating particles occupy ids [CaseNpb, CaseNpb+CaseNfloat). We build bonds
  // only from those particles to boundary/floating particles of a different DEM
  // material slot (tav = CODE_GetTypeAndValue(code)).
  const unsigned idini=CaseNpb;
  const unsigned idfin=CaseNpb+CaseNfloat;
  std::vector<tdouble3> posb(CaseNbound);
  std::vector<typecode> tavb(CaseNbound,0);
  std::vector<byte> validb(CaseNbound,0);
  const unsigned* idp=Idp_c->cptr();
  const typecode* code=Code_c->cptr();
  const tdouble3* posc=AuxPos_c->cptr();
  for(unsigned p=0;p<Np;p++){
    const unsigned id=idp[p];
    if(id<CaseNbound){
      // Convert from current particle order to stable boundary/floating ids.
      // The bond graph is stored in stable-id space so it remains valid after
      // later particle reordering by the cell divide.
      posb[id]=posc[p];
      tavb[id]=CODE_GetTypeAndValue(code[p]);
      validb[id]=1;
    }
  }
  for(unsigned id=0;id<CaseNbound;id++){
    if(!validb[id])Run_Exceptioon("There are noboundary/floating particles found while initialising DEM bonds.");
  }

  // StBondPair stores the unique physical bond.
  // StBondRef stores one directed adjacency entry for the CSR traversal of a
  // floating particle. Floating-floating pairs will later receive a mirrored ref.
  struct StBondPair{ unsigned id1,id2; double l0; };
  struct StBondRef{ unsigned fp1,ngbid,pairid,owner; double l0; };
  std::vector<StBondPair> pairs;
  std::vector<StBondRef> refs;
  std::vector<unsigned> counts(CaseNfloat,0);
  // Initial bond search radius: neighbouring particles closer than ~dp are
  // considered bonded at startup.
  const double form_radius=double(Dp)*1.05;
  const double form_radius2=form_radius*form_radius;
  for(unsigned fp1=0;fp1<CaseNfloat;fp1++){
    const unsigned id1=idini+fp1;
    const tdouble3 p1=posb[id1];
    const typecode tv1=tavb[id1];
    for(unsigned id2=0;id2<CaseNbound;id2++){
      if(id1==id2)continue;
      if(tv1==tavb[id2])continue;
      const double dx=p1.x-posb[id2].x;
      const double dy=p1.y-posb[id2].y;
      const double dz=p1.z-posb[id2].z;
      const double rr2=dx*dx+dy*dy+dz*dz;
      if(rr2<=form_radius2){
        if(id2>=idini && id2<idfin && id1<id2){
          // Floating-floating bond: create one physical pair and one owner ref.
          // The mirrored non-owner ref is added later once the CSR offsets exist.
          const unsigned pid=unsigned(pairs.size());
          StBondPair pb={id1,id2,sqrt(rr2)};
          pairs.push_back(pb);
          refs.push_back((StBondRef){fp1,id2,pid,1u,pb.l0});
          counts[fp1]++;
          counts[id2-idini]++;
        }
        else if(id2<idini){
          // Floating-bound bond: only the floating side needs a traversal entry.
          const unsigned pid=unsigned(pairs.size());
          StBondPair pb={id1,id2,sqrt(rr2)};
          pairs.push_back(pb);
          refs.push_back((StBondRef){fp1,id2,pid,1u,pb.l0});
          counts[fp1]++;
        }
      }
    }
  }
  if(pairs.empty()){
    Log->Print("DEM bonds: no initial pairs were created.");
    return;
  }

  struct StBondInterface{
    unsigned block1,block2;
    unsigned tav1,tav2;
    unsigned nbonds;
    double lseg;
    double aface;
    bool rectfit;
    bool expliciti;
    string rectreason;
    bool line2d;
    unsigned nu,nv;
    unsigned normalaxis,uaxis,vaxis;
    double lu,lv;
    double area_disc,iuu_target,iuu_disc,ivv_target,ivv_disc,iuv_disc;
    vector<unsigned> pairids;
    StBondInterface():block1(UINT_MAX),block2(UINT_MAX),tav1(UINT_MAX),tav2(UINT_MAX),nbonds(0),lseg(1.0),aface(1.0)
      ,rectfit(false),expliciti(false),line2d(false),nu(0),nv(0),normalaxis(UINT_MAX),uaxis(UINT_MAX),vaxis(UINT_MAX),lu(0),lv(0),area_disc(0),iuu_target(0),iuu_disc(0),ivv_target(0),ivv_disc(0),iuv_disc(0){}
  };
  std::vector<StBondInterface> interfaces;
  const unsigned mkcount=MkInfo->Size();
  std::vector<unsigned> ifidx(mkcount*mkcount,UINT_MAX);
  std::vector<unsigned> ifaceid(pairs.size(),UINT_MAX);
  for(unsigned pid=0;pid<pairs.size();pid++){
    const unsigned id1=pairs[pid].id1;
    const unsigned id2=pairs[pid].id2;
    const unsigned blockid1=MkInfo->GetMkBlockById(id1);
    const unsigned blockid2=MkInfo->GetMkBlockById(id2);
    if(blockid1>=mkcount || blockid2>=mkcount || blockid1==blockid2)
      Run_Exceptioon("Invalid rigid-body interface while initialising DEM bonds.");
    const unsigned low =(blockid1<blockid2? blockid1: blockid2);
    const unsigned high=(blockid1<blockid2? blockid2: blockid1);
    unsigned& iface=ifidx[low*mkcount+high];
    if(iface==UINT_MAX){
      // First bond seen for this unordered rigid-body pair. Under the current
      // model assumptions this pair corresponds to exactly one bonded interface.
      iface=unsigned(interfaces.size());
      interfaces.push_back(StBondInterface());
      interfaces.back().block1=low;
      interfaces.back().block2=high;
      const unsigned tav1=CODE_GetTypeAndValue(MkInfo->GetCodeById(id1));
      const unsigned tav2=CODE_GetTypeAndValue(MkInfo->GetCodeById(id2));
      interfaces.back().tav1=tav1;
      interfaces.back().tav2=tav2;
      // LSeg_Bond and AFace_Bond are currently treated as global/interface
      // constants stored in the material file. We average the two endpoint
      // materials so mixed-material interfaces still have one effective value.
      const double lseg1=(DemData[tav1].lsegbond>0.f? DemData[tav1].lsegbond: 1.f);
      const double lseg2=(DemData[tav2].lsegbond>0.f? DemData[tav2].lsegbond: 1.f);
      const double aface1=(DemData[tav1].afacebond>0.f? DemData[tav1].afacebond: 1.f);
      const double aface2=(DemData[tav2].afacebond>0.f? DemData[tav2].afacebond: 1.f);
      interfaces.back().lseg=0.5*(lseg1+lseg2);
      interfaces.back().aface=0.5*(aface1+aface2);
      const double iuu1=(DemData[tav1].iuubond>0.f? DemData[tav1].iuubond: -1.0);
      const double iuu2=(DemData[tav2].iuubond>0.f? DemData[tav2].iuubond: -1.0);
      const double ivv1=(DemData[tav1].ivvbond>0.f? DemData[tav1].ivvbond: -1.0);
      const double ivv2=(DemData[tav2].ivvbond>0.f? DemData[tav2].ivvbond: -1.0);
      if(Simulate2D){
        const double iuutf_tol=1e-12;
        if((iuu1>=0.0 && fabs(iuu1)>iuutf_tol) || (iuu2>=0.0 && fabs(iuu2)>iuutf_tol)){
          Run_Exceptioon("DEM bonded 2D interfaces require Iuu_Bond to be zero or unset on both endpoint materials.");
        }
        if(ivv1<=0.0 || ivv2<=0.0){
          Run_Exceptioon("DEM bonded 2D interfaces require positive Ivv_Bond on both endpoint materials.");
        }
        interfaces.back().iuu_target=0.0;
        interfaces.back().ivv_target=0.5*(ivv1+ivv2);
      }
      else{
        interfaces.back().iuu_target=(iuu1>0.0 && iuu2>0.0? 0.5*(iuu1+iuu2): (iuu1>0.0? iuu1: iuu2));
        interfaces.back().ivv_target=(ivv1>0.0 && ivv2>0.0? 0.5*(ivv1+ivv2): (ivv1>0.0? ivv1: ivv2));
        interfaces.back().expliciti=(interfaces.back().iuu_target>0.0 && interfaces.back().ivv_target>0.0);
      }
    }
    // Record which interface each physical bond belongs to and count how many
    // discrete springs share that same interface.
    ifaceid[pid]=iface;
    interfaces[iface].nbonds++;
    interfaces[iface].pairids.push_back(pid);
  }

  std::vector<double> pairaeff(pairs.size(),0.0);
  std::vector<double> pairu(pairs.size(),0.0);
  std::vector<double> pairv(pairs.size(),0.0);
  for(unsigned ifc=0;ifc<interfaces.size();ifc++){
    StBondInterface& iface=interfaces[ifc];
    const unsigned nbonds=(iface.nbonds? iface.nbonds: 1u);
    const double aeff_uniform=iface.aface/double(nbonds);
    std::vector<tdouble3> mids;
    std::vector<tdouble3> deltas;
    mids.reserve(iface.pairids.size());
    deltas.reserve(iface.pairids.size());
    for(unsigned c=0;c<iface.pairids.size();c++){
      const unsigned pid=iface.pairids[c];
      const StBondPair& pb=pairs[pid];
      const tdouble3 p1=posb[pb.id1];
      const tdouble3 p2=posb[pb.id2];
      mids.push_back((p1+p2)*0.5);
      deltas.push_back(p2-p1);
    }
    const StBondAreaFit fit=(Simulate2D
      ? FitLineBondAreas2D(mids,iface.aface,double(Dp),iface.ivv_target)
      : FitRectangularBondAreas(mids,deltas,iface.aface,double(Dp)
          ,(iface.expliciti? iface.iuu_target: -1.0),(iface.expliciti? iface.ivv_target: -1.0)));
    iface.rectfit=fit.ok;
    iface.rectreason=fit.reason;
    if(fit.ok){
      iface.line2d=fit.line2d;
      iface.nu=fit.nu;
      iface.nv=fit.nv;
      iface.normalaxis=fit.normalaxis;
      iface.uaxis=fit.uaxis;
      iface.vaxis=fit.vaxis;
      iface.lu=fit.lu;
      iface.lv=fit.lv;
      iface.area_disc=fit.area_disc;
      iface.iuu_target=fit.iuu_target;
      iface.iuu_disc=fit.iuu_disc;
      iface.ivv_target=fit.ivv_target;
      iface.ivv_disc=fit.ivv_disc;
      iface.iuv_disc=fit.iuv_disc;
      for(unsigned c=0;c<iface.pairids.size();c++){
        const unsigned pid=iface.pairids[c];
        pairaeff[pid]=fit.aeff[c];
        pairu[pid]=fit.u[c];
        pairv[pid]=fit.v[c];
      }
    }
    else{
      for(unsigned c=0;c<iface.pairids.size();c++){
        const unsigned pid=iface.pairids[c];
        pairaeff[pid]=aeff_uniform;
      }
    }
  }

  std::vector<unsigned> start(CaseNfloat+1,0);
  for(unsigned fp=0;fp<CaseNfloat;fp++)start[fp+1]=start[fp]+counts[fp];
  DemBondPairCount=unsigned(pairs.size());
  DemBondRefCount=start[CaseNfloat];
  std::vector<unsigned> ngb(DemBondRefCount);
  std::vector<double> l0(DemBondRefCount);
  std::vector<unsigned> pairid(DemBondRefCount);
  std::vector<unsigned> owner(DemBondRefCount);
  std::vector<unsigned> cursor(start.begin(),start.end()-1);
  for(unsigned r=0;r<refs.size();r++){
    // Owner refs are the only refs that update tangential history.
    const StBondRef& rb=refs[r];
    unsigned cp=cursor[rb.fp1]++;
    ngb[cp]=rb.ngbid; l0[cp]=rb.l0; pairid[cp]=rb.pairid; owner[cp]=rb.owner;
  }
  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    const StBondPair& pb=pairs[pid];
    if(pb.id2>=idini && pb.id2<idfin){
      // Mirrored non-owner ref for the other floating endpoint. This lets both
      // particles apply equal and opposite forces while only one side evolves ut.
      const unsigned fp2=pb.id2-idini;
      unsigned cp=cursor[fp2]++;
      ngb[cp]=pb.id1; l0[cp]=pb.l0; pairid[cp]=pid; owner[cp]=0u;
    }
  }
  DemBondPairId1h.resize(DemBondPairCount);
  DemBondPairId2h.resize(DemBondPairCount);
  DemBondPairL0h.resize(DemBondPairCount);
  DemBondPairIfcIdh.resize(DemBondPairCount);
  DemBondPairUh.resize(DemBondPairCount);
  DemBondPairVh.resize(DemBondPairCount);
  DemBondPairNbondsh.resize(DemBondPairCount);
  DemBondPairLsegh.resize(DemBondPairCount);
  DemBondPairAfaceh.resize(DemBondPairCount);
  DemBondPairAeffh.resize(DemBondPairCount);
  DemBondPairKnh.resize(DemBondPairCount);
  DemBondIfcBlock1h.resize(interfaces.size());
  DemBondIfcBlock2h.resize(interfaces.size());
  DemBondIfcNuh.resize(interfaces.size());
  DemBondIfcNvh.resize(interfaces.size());
  DemBondIfcLuh.resize(interfaces.size());
  DemBondIfcLvh.resize(interfaces.size());
  DemBondIfcAfaceh.resize(interfaces.size());
  DemBondIfcIuuh.resize(interfaces.size());
  DemBondIfcIvvh.resize(interfaces.size());
  DemBondIfcIuvh.resize(interfaces.size());
  for(unsigned ifc=0;ifc<interfaces.size();ifc++){
    const StBondInterface& iface=interfaces[ifc];
    DemBondIfcBlock1h[ifc]=iface.block1;
    DemBondIfcBlock2h[ifc]=iface.block2;
    DemBondIfcNuh[ifc]=iface.nu;
    DemBondIfcNvh[ifc]=iface.nv;
    DemBondIfcLuh[ifc]=iface.lu;
    DemBondIfcLvh[ifc]=iface.lv;
    DemBondIfcAfaceh[ifc]=iface.aface;
    DemBondIfcIuuh[ifc]=iface.iuu_target;
    DemBondIfcIvvh[ifc]=iface.ivv_target;
    DemBondIfcIuvh[ifc]=0.0;
  }
  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    // Store the bond rest geometry exactly as formed at t=0. Runtime extension
    // is always measured against this l0 rather than against LSeg_Bond.
    DemBondPairId1h[pid]=pairs[pid].id1;
    DemBondPairId2h[pid]=pairs[pid].id2;
    DemBondPairL0h[pid]=pairs[pid].l0;
    const unsigned ifc=ifaceid[pid];
    const StBondInterface& iface=interfaces[ifc];
    const unsigned nbonds=(iface.nbonds? iface.nbonds: 1u);
    const unsigned tav1=CODE_GetTypeAndValue(MkInfo->GetCodeById(pairs[pid].id1));
    const unsigned tav2=CODE_GetTypeAndValue(MkInfo->GetCodeById(pairs[pid].id2));
    const byte failmode1=DemData[tav1].bondfailmode;
    const byte failmode2=DemData[tav2].bondfailmode;
    if(failmode1 && failmode2 && failmode1!=failmode2){
      Run_Exceptioon(fun::PrintStr("DEM bonded pair with ids %u and %u mixes Sigma_C and EpsTension_C failure modes."
        ,pairs[pid].id1,pairs[pid].id2));
    }
    const double young_pair=0.5*(double(DemData[tav1].young)+double(DemData[tav2].young));
    const double lseg=(iface.lseg>ALMOSTZERO? iface.lseg: 1.0);
    const double aface=(iface.aface>ALMOSTZERO? iface.aface: 1.0);
    // Each persistent bond stores its own representative interface area.
    const double aeff=pairaeff[pid];
    // Continuum-consistent bonded normal spring stiffness kn = E * Aeff / LSeg.
    const double kn=young_pair*aeff/lseg;
    DemBondPairIfcIdh[pid]=ifc;
    DemBondPairUh[pid]=pairu[pid];
    DemBondPairVh[pid]=pairv[pid];
    DemBondPairNbondsh[pid]=nbonds;
    DemBondPairLsegh[pid]=lseg;
    DemBondPairAfaceh[pid]=aface;
    DemBondPairAeffh[pid]=aeff;
    DemBondPairKnh[pid]=kn;
  }
  MemGpuFixed+=fcuda::Malloc(&DemBondStartg,CaseNfloat+1);
  MemGpuFixed+=fcuda::Malloc(&DemBondNgbg,DemBondRefCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondL0g,DemBondRefCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondPairIdg,DemBondRefCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondOwnerg,DemBondRefCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondStateg,DemBondPairCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondUtg,DemBondPairCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondLsegg,DemBondPairCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondKng,DemBondPairCount);
  MemGpuFixed+=fcuda::Malloc(&DemBondAeffg,DemBondPairCount);
  cudaMemcpy(DemBondStartg,&start[0],sizeof(unsigned)*(CaseNfloat+1),cudaMemcpyHostToDevice);
  cudaMemcpy(DemBondNgbg,&ngb[0],sizeof(unsigned)*DemBondRefCount,cudaMemcpyHostToDevice);
  cudaMemcpy(DemBondL0g,&l0[0],sizeof(double)*DemBondRefCount,cudaMemcpyHostToDevice);
  cudaMemcpy(DemBondPairIdg,&pairid[0],sizeof(unsigned)*DemBondRefCount,cudaMemcpyHostToDevice);
  cudaMemcpy(DemBondOwnerg,&owner[0],sizeof(unsigned)*DemBondRefCount,cudaMemcpyHostToDevice);
  std::vector<float> state(DemBondPairCount,1.f);
  cudaMemcpy(DemBondStateg,&state[0],sizeof(float)*DemBondPairCount,cudaMemcpyHostToDevice);
  // Tangential history starts at zero because the initial configuration is the
  // undeformed bonded reference state.
  cudaMemset(DemBondUtg,0,sizeof(double3)*DemBondPairCount);
  cudaMemcpy(DemBondLsegg,&DemBondPairLsegh[0],sizeof(double)*DemBondPairCount,cudaMemcpyHostToDevice);
  // Runtime force evaluation reads precomputed kn and Aeff directly from GPU.
  cudaMemcpy(DemBondKng,&DemBondPairKnh[0],sizeof(double)*DemBondPairCount,cudaMemcpyHostToDevice);
  cudaMemcpy(DemBondAeffg,&DemBondPairAeffh[0],sizeof(double)*DemBondPairCount,cudaMemcpyHostToDevice);
  Check_CudaErroor("Failed initialising DEM bonds.");
  ValidateDemBondGpu();
  Log->Printf("DEM bonds: %u persistent pairs created (form radius=%g).",DemBondPairCount,float(form_radius));
  for(unsigned c=0;c<interfaces.size();c++){
    const StBondInterface& iface=interfaces[c];
    const JSphMkBlock* mk1=MkInfo->Mkblock(iface.block1);
    const JSphMkBlock* mk2=MkInfo->Mkblock(iface.block2);
    const double young_pair=0.5*(double(DemData[iface.tav1].young)+double(DemData[iface.tav2].young));
    const double nu_pair=0.5*(double(DemData[iface.tav1].poisson)+double(DemData[iface.tav2].poisson));
    const double aeff_avg=iface.aface/double(iface.nbonds? iface.nbonds: 1u);
    const double kn_pair=young_pair*aeff_avg/(iface.lseg>ALMOSTZERO? iface.lseg: 1.0);
    double ks_pair=0.0;
    if(DemData[iface.tav1].ksbond>=0.f || DemData[iface.tav2].ksbond>=0.f){
      if(DemData[iface.tav1].ksbond>=0.f && DemData[iface.tav2].ksbond>=0.f)ks_pair=0.5*(double(DemData[iface.tav1].ksbond)+double(DemData[iface.tav2].ksbond));
      else ks_pair=(DemData[iface.tav1].ksbond>=0.f? DemData[iface.tav1].ksbond: DemData[iface.tav2].ksbond);
    }
    else{
      const double kappa=2.0*(1.0+nu_pair);
      ks_pair=(kappa>ALMOSTZERO? kn_pair/kappa: kn_pair);
    }
    double aeffmin=DBL_MAX, aeffmax=0;
    for(unsigned p=0;p<iface.pairids.size();p++){
      const double aeff=pairaeff[iface.pairids[p]];
      aeffmin=min(aeffmin,aeff);
      aeffmax=max(aeffmax,aeff);
    }
    if(Simulate2D){
      Log->Printf("  DEM bonded interface mkbound:%u-%u mk:%u-%u shape=line2d nbonds=%u grid=%ux%u axes(n,u,v)=%s,%s,%s lseg=%g aface=%g area=%g aeff=[%g,%g] kn_avg=%g ks_pair=%g Iuu=%g/%g Ivv=%g/%g"
        ,mk1->MkType,mk2->MkType,mk1->Mk,mk2->Mk
        ,iface.nbonds,iface.nu,iface.nv,AxisName(iface.normalaxis),AxisName(iface.uaxis),AxisName(iface.vaxis)
        ,iface.lseg,iface.aface,iface.area_disc,aeffmin,aeffmax,kn_pair,ks_pair,iface.iuu_disc,iface.iuu_target,iface.ivv_disc,iface.ivv_target);
    }
    else if(iface.rectfit){
      Log->Printf("  DEM bonded interface mkbound:%u-%u mk:%u-%u nbonds=%u grid=%ux%u (u=long, v=short) lseg=%g aface=%g aeff=[%g,%g] kn_avg=%g ks_pair=%g Iuu=%g/%g Ivv=%g/%g (%s)"
        ,mk1->MkType,mk2->MkType,mk1->Mk,mk2->Mk,iface.nbonds,iface.nu,iface.nv,iface.lseg,iface.aface,aeffmin,aeffmax,kn_pair,ks_pair
        ,iface.iuu_disc,iface.iuu_target,iface.ivv_disc,iface.ivv_target,(iface.expliciti? "explicit I targets": "geometry inferred"));
    }
    else{
      Log->Printf("  DEM bonded interface mkbound:%u-%u mk:%u-%u nbonds=%u lseg=%g aface=%g aeff=%g kn_pair=%g ks_pair=%g (uniform Aeff fallback: %s)"
        ,mk1->MkType,mk2->MkType,mk1->Mk,mk2->Mk,iface.nbonds,iface.lseg,iface.aface,aeffmin,kn_pair,ks_pair,iface.rectreason.c_str());
    }
  }
}

//==============================================================================
/// Copies the uploaded DEM bond graph back once and checks index consistency.
//==============================================================================
void JSphGpu::ValidateDemBondGpu()const{
  if(!UseDEM || !DemBondPairCount || !DemBondRefCount)return;

  std::vector<unsigned> start(CaseNfloat+1,0);
  std::vector<unsigned> ngb(DemBondRefCount,UINT_MAX);
  std::vector<unsigned> pairid(DemBondRefCount,UINT_MAX);
  std::vector<unsigned> owner(DemBondRefCount,UINT_MAX);
  std::vector<float> state(DemBondPairCount,0.f);
  std::vector<double3> ut(DemBondPairCount);
  std::vector<double> kn(DemBondPairCount,0.0);
  std::vector<double> aeff(DemBondPairCount,0.0);

  cudaMemcpy(&start[0],DemBondStartg,sizeof(unsigned)*(CaseNfloat+1),cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond start offsets.");
  cudaMemcpy(&ngb[0],DemBondNgbg,sizeof(unsigned)*DemBondRefCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond neighbour ids.");
  cudaMemcpy(&pairid[0],DemBondPairIdg,sizeof(unsigned)*DemBondRefCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond pair ids.");
  cudaMemcpy(&owner[0],DemBondOwnerg,sizeof(unsigned)*DemBondRefCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond owner flags.");
  cudaMemcpy(&state[0],DemBondStateg,sizeof(float)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond state buffer.");
  cudaMemcpy(&ut[0],DemBondUtg,sizeof(double3)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond tangential-history buffer.");
  cudaMemcpy(&kn[0],DemBondKng,sizeof(double)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond normal-stiffness buffer.");
  cudaMemcpy(&aeff[0],DemBondAeffg,sizeof(double)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed validating DEM bond effective-area buffer.");

  if(start[0]!=0)Run_Exceptioon("DEM bond validation failed: start[0] must be zero.");
  if(start[CaseNfloat]!=DemBondRefCount){
    Run_Exceptioon(fun::PrintStr("DEM bond validation failed: start[%u]=%u but DemBondRefCount=%u."
      ,CaseNfloat,start[CaseNfloat],DemBondRefCount));
  }
  for(unsigned fp=0;fp<CaseNfloat;fp++){
    if(start[fp]>start[fp+1] || start[fp+1]>DemBondRefCount){
      Run_Exceptioon(fun::PrintStr("DEM bond validation failed: invalid CSR range for floating particle %u (%u,%u) with refcount=%u."
        ,fp,start[fp],start[fp+1],DemBondRefCount));
    }
  }
  for(unsigned c=0;c<DemBondRefCount;c++){
    if(ngb[c]>=CaseNbound){
      Run_Exceptioon(fun::PrintStr("DEM bond validation failed: neighbour id %u out of range at ref %u (CaseNbound=%u)."
        ,ngb[c],c,CaseNbound));
    }
    if(pairid[c]>=DemBondPairCount){
      Run_Exceptioon(fun::PrintStr("DEM bond validation failed: pair id %u out of range at ref %u (paircount=%u)."
        ,pairid[c],c,DemBondPairCount));
    }
    if(owner[c]>1u){
      Run_Exceptioon(fun::PrintStr("DEM bond validation failed: owner flag %u is invalid at ref %u."
        ,owner[c],c));
    }
  }
  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    if(!(state[pid]==0.f || state[pid]==1.f)){
      Run_Exceptioon(fun::PrintStr("DEM bond validation failed: unexpected initial state %g at pair %u."
        ,double(state[pid]),pid));
    }
    if(!(aeff[pid]>0.0) || !(kn[pid]>0.0)){
      Run_Exceptioon(fun::PrintStr("DEM bond validation failed: non-positive uploaded coefficients at pair %u (aeff=%g kn=%g)."
        ,pid,aeff[pid],kn[pid]));
    }
  }
}

//==============================================================================
/// Builds a per-particle bonded interface shear-force field for output.
/// Each bonded interface is identified by the pair of MK blocks connected by the
/// intact bonds. The total interface shear-force vector is summed once per interface
/// and assigned uniformly to all particles on each side of that interface. Particles
/// participating in multiple interfaces receive the average of their incident
/// interface vectors.
//==============================================================================
void JSphGpu::ComputeDemShearForceOutput(unsigned np,const unsigned* idp,tfloat3* shearforce,float* shearmag)const{
  if(shearforce)memset(shearforce,0,sizeof(tfloat3)*np);
  if(shearmag)memset(shearmag,0,sizeof(float)*np);
  if(!UseDEM || !np || !idp || !shearforce || !DemBondPairCount || !DemBondStateg || !DemBondUtg)return;

  // Map stable particle ids to the current saved-particle indices. Output is
  // written after filtering/periodic removal, so we work in saved-particle order.
  std::vector<unsigned> ridp(CaseNbound,UINT_MAX);
  for(unsigned p=0;p<np;p++){
    const unsigned id=idp[p];
    if(id<CaseNbound)ridp[id]=p;
  }

  std::vector<float> bondstate(DemBondPairCount,0.f);
  std::vector<double3> bondut(DemBondPairCount);
  std::vector<double> bondkn(DemBondPairCount,0.0);
  std::vector<double> bondaeff(DemBondPairCount,0.0);
  const tdouble3* posc=AuxPos_c->cptr();
  const tfloat3* velc=AuxVel_c->cptr();
  cudaMemcpy(&bondstate[0],DemBondStateg,sizeof(float)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying DEM bonded state data for output.");
  cudaMemcpy(&bondut[0],DemBondUtg,sizeof(double3)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying DEM bonded tangential-history data for output.");
  cudaMemcpy(&bondkn[0],DemBondKng,sizeof(double)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying DEM bonded normal-stiffness data for output.");
  cudaMemcpy(&bondaeff[0],DemBondAeffg,sizeof(double)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying DEM bonded effective-area data for output.");

  // Each interface is identified by the pair of MK blocks it connects. For the
  // output field we want one total interface force per side, not one value per bond.
  struct StInterfaceData{
    unsigned block1,block2;
    tfloat3 force1;
    std::vector<unsigned> ids1;
    std::vector<unsigned> ids2;
    StInterfaceData():block1(UINT_MAX),block2(UINT_MAX),force1(TFloat3(0)){}
  };
  std::vector<StInterfaceData> interfaces;

  const unsigned mkcount=MkInfo->Size();
  std::vector<unsigned> ifidx(mkcount*mkcount,UINT_MAX);

  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    if(bondstate[pid]<=0.5f)continue;
    const unsigned id1=DemBondPairId1h[pid];
    const unsigned id2=DemBondPairId2h[pid];
    if(id1>=CaseNbound || id2>=CaseNbound)continue;

    // Interface grouping:
    // - GetMkBlockById(id) returns the MK block that owns each particle id.
    // - We define one "interface" as one unordered pair of MK blocks.
    // - Therefore every intact bond between the same two MK blocks contributes
    //   to the same interface-total shear force, even if many particle-particle
    //   bonds exist across that physical face.
    // - low/high provides a unique order so block pair (A,B) and (B,A) map to
    //   the same bucket in ifidx[].
    const unsigned blockid1=MkInfo->GetMkBlockById(id1);
    const unsigned blockid2=MkInfo->GetMkBlockById(id2);
    if(blockid1>=mkcount || blockid2>=mkcount || blockid1==blockid2)continue;
    const unsigned low =(blockid1<blockid2? blockid1: blockid2);
    const unsigned high=(blockid1<blockid2? blockid2: blockid1);
    unsigned& iface=ifidx[low*mkcount+high];
    if(iface==UINT_MAX){
      // First bond seen for this MK-block pair: create a new interface record.
      iface=unsigned(interfaces.size());
      interfaces.push_back(StInterfaceData());
      interfaces.back().block1=low;
      interfaces.back().block2=high;
    }

    // Recompute the same bonded shear force used by the solver:
    //   Ft = -ks_pair * ut - ct_pair * v_tan
    const unsigned tav1=CODE_GetTypeAndValue(MkInfo->GetCodeById(id1));
    const unsigned tav2=CODE_GetTypeAndValue(MkInfo->GetCodeById(id2));
    const unsigned p1=ridp[id1];
    const unsigned p2=ridp[id2];
    if(p1==UINT_MAX || p2==UINT_MAX)continue;
    double ks_pair=0.0;
    if(DemData[tav1].ksbond>=0.f || DemData[tav2].ksbond>=0.f){
      if(DemData[tav1].ksbond>=0.f && DemData[tav2].ksbond>=0.f)ks_pair=0.5*(double(DemData[tav1].ksbond)+double(DemData[tav2].ksbond));
      else ks_pair=(DemData[tav1].ksbond>=0.f? DemData[tav1].ksbond: DemData[tav2].ksbond);
    }
    else{
      const double kn_pair=DemBondPairKnh[pid];
      const double nu_pair=0.5*(double(DemData[tav1].poisson)+double(DemData[tav2].poisson));
      const double kappa=2.0*(1.0+nu_pair);
      ks_pair=(kappa>ALMOSTZERO? kn_pair/kappa: kn_pair);
    }

    const tdouble3 dr=TDouble3(posc[p1].x-posc[p2].x,posc[p1].y-posc[p2].y,posc[p1].z-posc[p2].z);
    const double rr2=dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;
    if(rr2<ALMOSTZERO)continue;
    const double rad=sqrt(rr2);
    const double nx=dr.x/rad, ny=dr.y/rad, nz=dr.z/rad;
    const double dvx=double(velc[p1].x)-double(velc[p2].x);
    const double dvy=double(velc[p1].y)-double(velc[p2].y);
    const double dvz=double(velc[p1].z)-double(velc[p2].z);
    const double vn=dvx*nx+dvy*ny+dvz*nz;
    const double dvxt=dvx-vn*nx;
    const double dvyt=dvy-vn*ny;
    const double dvzt=dvz-vn*nz;
    const double cn_pair=0.5*(double(DemData[tav1].cnbond)+double(DemData[tav2].cnbond));
    const double3 ut=bondut[pid];
    const tfloat3 ft=TFloat3(float(-ks_pair*ut.x-cn_pair*dvxt),float(-ks_pair*ut.y-cn_pair*dvyt),float(-ks_pair*ut.z-cn_pair*dvzt));
    StInterfaceData& idata=interfaces[iface];
    // Store the interface-total force using the low/high block ordering as the
    // internal sign convention. The opposite side will receive the opposite vector.
    const tfloat3 force_low=(blockid1==low? ft: TFloat3(-ft.x,-ft.y,-ft.z));
    idata.force1=idata.force1+force_low;
    idata.ids1.push_back(blockid1==low? id1: id2);
    idata.ids2.push_back(blockid1==low? id2: id1);
  }

  std::vector<unsigned> counts(np,0);
  for(unsigned c=0;c<interfaces.size();c++){
    StInterfaceData& idata=interfaces[c];
    // A particle can touch the same neighbour through multiple bonds; for the
    // interface field we assign the interface-total vector once per particle.
    std::sort(idata.ids1.begin(),idata.ids1.end());
    idata.ids1.erase(std::unique(idata.ids1.begin(),idata.ids1.end()),idata.ids1.end());
    std::sort(idata.ids2.begin(),idata.ids2.end());
    idata.ids2.erase(std::unique(idata.ids2.begin(),idata.ids2.end()),idata.ids2.end());
    const tfloat3 force2=TFloat3(-idata.force1.x,-idata.force1.y,-idata.force1.z);
    // Every particle on side 1 gets the same interface-total vector force1.
    // Every particle on side 2 gets the opposite vector force2.
    for(unsigned i=0;i<idata.ids1.size();i++){
      const unsigned p=ridp[idata.ids1[i]];
      if(p!=UINT_MAX){
        shearforce[p]=shearforce[p]+idata.force1;
        counts[p]++;
      }
    }
    for(unsigned i=0;i<idata.ids2.size();i++){
      const unsigned p=ridp[idata.ids2[i]];
      if(p!=UINT_MAX){
        shearforce[p]=shearforce[p]+force2;
        counts[p]++;
      }
    }
  }

  // Corner/edge particles may belong to more than one interface. Average the
  // interface vectors so one particle still carries a single output vector.
  for(unsigned p=0;p<np;p++){
    if(counts[p]>1){
      const float inv=1.f/float(counts[p]);
      shearforce[p].x*=inv;
      shearforce[p].y*=inv;
      shearforce[p].z*=inv;
    }
  }

  if(shearmag){
    for(unsigned p=0;p<np;p++){
      const tfloat3 f=shearforce[p];
      shearmag[p]=sqrt(f.x*f.x+f.y*f.y+f.z*f.z);
    }
  }
}

//==============================================================================
/// Builds a three-state DEM bond field for visualisation:
/// - 2: particle was never part of the initial bond graph
/// - 1: particle still has at least one intact bond
/// - 0: particle had bonds initially, but all of them are now broken
//==============================================================================
void JSphGpu::ComputeDemBondOutput(unsigned np,const unsigned* idp,float* dembond)const{
  if(!dembond)return;
  for(unsigned p=0;p<np;p++)dembond[p]=2.f;
  if(!UseDEM || !np || !idp || !DemBondPairCount || !DemBondStateg)return;

  // Map stable particle ids to the current saved-particle indices. Output is
  // written after filtering/periodic removal, so we work in saved-particle order.
  std::vector<unsigned> ridp(CaseNbound,UINT_MAX);
  for(unsigned p=0;p<np;p++){
    const unsigned id=idp[p];
    if(id<CaseNbound)ridp[id]=p;
  }

  std::vector<float> bondstate(DemBondPairCount,0.f);
  cudaMemcpy(&bondstate[0],DemBondStateg,sizeof(float)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying DEM bonded state data for bond-visualisation output.");

  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    const unsigned id1=DemBondPairId1h[pid];
    const unsigned id2=DemBondPairId2h[pid];
    const bool intact=(bondstate[pid]>0.5f);
    if(id1<CaseNbound){
      const unsigned p=ridp[id1];
      if(p!=UINT_MAX){
        if(intact)dembond[p]=1.f;
        else if(dembond[p]>1.f)dembond[p]=0.f;
      }
    }
    if(id2<CaseNbound){
      const unsigned p=ridp[id2];
      if(p!=UINT_MAX){
        if(intact)dembond[p]=1.f;
        else if(dembond[p]>1.f)dembond[p]=0.f;
      }
    }
  }
}

//==============================================================================
/// Dumps detailed DEM bond diagnostics for the first simulation steps.
//==============================================================================
void JSphGpu::DumpDemBondDebug(unsigned maxsteps){
  if(!UseDEM || !DemBondPairCount || !DemBondStateg || Nstep<1 || unsigned(Nstep)>maxsteps)return;

  unsigned npfilterdel=0;
  const unsigned nnormal=ParticlesDataDown(Np,0,true,true,NULL,npfilterdel,NULL,NULL,NULL,NULL);
  const unsigned* idpc=Idp_c->cptr();
  const tdouble3* posc=AuxPos_c->cptr();
  const typecode* codec=Code_c->cptr();
  std::vector<unsigned> ridp(CaseNbound,UINT_MAX);
  for(unsigned p=0;p<nnormal;p++){
    const unsigned id=idpc[p];
    if(id<CaseNbound)ridp[id]=p;
  }

  std::vector<float> bondstate(DemBondPairCount,0.f);
  cudaMemcpy(&bondstate[0],DemBondStateg,sizeof(float)*DemBondPairCount,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying DEM bond state for debug dump.");

  const string fname=DirOut+"DemBondDebug.csv";
  const bool create=(Nstep==1);
  ofstream pf(fname.c_str(),(create? ios::out: (ios::out|ios::app)));
  if(!pf.is_open())Run_Exceptioon(fun::PrintStr("Cannot open DEM bond debug file '%s'.",fname.c_str()));
  pf.setf(ios::scientific);
  pf<<setprecision(9);
  if(create){
    // Write both bond-level geometry/state and interface-level stiffness data so
    // the first 20 steps can be checked against the continuum interpretation.
    pf<<"step,time,pair_id,iface_id,state,id1,id2,mkbound1,mkbound2,mk1,mk2,tv1,tv2,x1,y1,z1,x2,y2,z2,rad,l0,ext,u0,v0,nbonds_interface,lseg,aface,aeff,kn_pair,ks_pair,fn_trial,sigma_signed,sigma_tension\n";
    Log->AddFileInfo(fname,"Detailed DEM bond diagnostics for the first 20 simulation steps.");
  }
  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    const unsigned id1=DemBondPairId1h[pid];
    const unsigned id2=DemBondPairId2h[pid];
    if(id1>=CaseNbound || id2>=CaseNbound)continue;
    const unsigned p1=ridp[id1];
    const unsigned p2=ridp[id2];
    if(p1==UINT_MAX || p2==UINT_MAX)continue;
    const tdouble3 ps1=posc[p1];
    const tdouble3 ps2=posc[p2];
    const double dx=ps1.x-ps2.x;
    const double dy=ps1.y-ps2.y;
    const double dz=ps1.z-ps2.z;
    const double rad=sqrt(dx*dx+dy*dy+dz*dz);
    const double l0=DemBondPairL0h[pid];
    const double ext=rad-l0;
    const unsigned tav1=CODE_GetTypeAndValue(codec[p1]);
    const unsigned tav2=CODE_GetTypeAndValue(codec[p2]);
    const unsigned tv1=CODE_GetTypeValue(codec[p1]);
    const unsigned tv2=CODE_GetTypeValue(codec[p2]);
    const unsigned cb1=MkInfo->GetMkBlockById(id1);
    const unsigned cb2=MkInfo->GetMkBlockById(id2);
    const word mkbound1=(cb1<MkInfo->Size()? MkInfo->Mkblock(cb1)->MkType: 0);
    const word mkbound2=(cb2<MkInfo->Size()? MkInfo->Mkblock(cb2)->MkType: 0);
    const word mk1=(cb1<MkInfo->Size()? MkInfo->Mkblock(cb1)->Mk: 0);
    const word mk2=(cb2<MkInfo->Size()? MkInfo->Mkblock(cb2)->Mk: 0);
    // Debug stress output uses the same Aeff and kn that the runtime bonded
    // kernel uses for force and failure.
    const double area=DemBondPairAeffh[pid];
    const double kn_pair=DemBondPairKnh[pid];
    double ks_pair=0.0;
    if(DemData[tav1].ksbond>=0.f || DemData[tav2].ksbond>=0.f){
      if(DemData[tav1].ksbond>=0.f && DemData[tav2].ksbond>=0.f)ks_pair=0.5*(double(DemData[tav1].ksbond)+double(DemData[tav2].ksbond));
      else ks_pair=(DemData[tav1].ksbond>=0.f? DemData[tav1].ksbond: DemData[tav2].ksbond);
    }
    else{
      const double nu_pair=0.5*(double(DemData[tav1].poisson)+double(DemData[tav2].poisson));
      const double kappa=2.0*(1.0+nu_pair);
      ks_pair=(kappa>ALMOSTZERO? kn_pair/kappa: kn_pair);
    }
    const double fn_trial=-kn_pair*ext;
    const double sigma_signed=-fn_trial/area;
    const double sigma_tension=(fn_trial<0.0? -fn_trial/area: 0.0);
    pf<<Nstep<<","<<TimeStep<<","<<pid
      <<","<<DemBondPairIfcIdh[pid]
      <<","<<bondstate[pid]
      <<","<<id1<<","<<id2
      <<","<<mkbound1<<","<<mkbound2
      <<","<<mk1<<","<<mk2
      <<","<<tv1<<","<<tv2
      <<","<<ps1.x<<","<<ps1.y<<","<<ps1.z
      <<","<<ps2.x<<","<<ps2.y<<","<<ps2.z
      <<","<<rad<<","<<l0<<","<<ext
      <<","<<DemBondPairUh[pid]
      <<","<<DemBondPairVh[pid]
      <<","<<DemBondPairNbondsh[pid]
      <<","<<DemBondPairLsegh[pid]
      <<","<<DemBondPairAfaceh[pid]
      <<","<<DemBondPairAeffh[pid]
      <<","<<kn_pair<<","<<ks_pair<<","<<fn_trial
      <<","<<sigma_signed<<","<<sigma_tension<<"\n";
  }

  const string fnameif=DirOut+"DemBondInterfaceDebug.csv";
  ofstream pif(fnameif.c_str(),(create? ios::out: (ios::out|ios::app)));
  if(!pif.is_open())Run_Exceptioon(fun::PrintStr("Cannot open DEM bond interface debug file '%s'.",fnameif.c_str()));
  pif.setf(ios::scientific);
  pif<<setprecision(9);
  if(create){
    pif<<"step,time,iface_id,mkbound1,mkbound2,mk1,mk2,nbonds_alive,aface_target,area_alive,iuu_target,iuu_alive,ivv_target,ivv_alive,iuv_target,iuv_alive\n";
    Log->AddFileInfo(fnameif,"Per-interface DEM bonded section diagnostics for the first 20 simulation steps.");
  }
  vector<unsigned> nbonds_alive(DemBondIfcAfaceh.size(),0);
  vector<double> area_alive(DemBondIfcAfaceh.size(),0.0);
  vector<double> iuu_alive(DemBondIfcAfaceh.size(),0.0);
  vector<double> ivv_alive(DemBondIfcAfaceh.size(),0.0);
  vector<double> iuv_alive(DemBondIfcAfaceh.size(),0.0);
  for(unsigned pid=0;pid<DemBondPairCount;pid++){
    if(bondstate[pid]<=0.5f)continue;
    const unsigned ifc=DemBondPairIfcIdh[pid];
    if(ifc>=DemBondIfcAfaceh.size())continue;
    const double area=DemBondPairAeffh[pid];
    const double u=DemBondPairUh[pid];
    const double v=DemBondPairVh[pid];
    nbonds_alive[ifc]++;
    area_alive[ifc]+=area;
    iuu_alive[ifc]+=area*v*v;
    ivv_alive[ifc]+=area*u*u;
    iuv_alive[ifc]+=area*u*v;
  }
  for(unsigned ifc=0;ifc<DemBondIfcAfaceh.size();ifc++){
    const unsigned block1=DemBondIfcBlock1h[ifc];
    const unsigned block2=DemBondIfcBlock2h[ifc];
    const JSphMkBlock* mkb1=(block1<MkInfo->Size()? MkInfo->Mkblock(block1): NULL);
    const JSphMkBlock* mkb2=(block2<MkInfo->Size()? MkInfo->Mkblock(block2): NULL);
    pif<<Nstep<<","<<TimeStep<<","<<ifc
      <<","<<(mkb1? mkb1->MkType: 0)<<","<<(mkb2? mkb2->MkType: 0)
      <<","<<(mkb1? mkb1->Mk: 0)<<","<<(mkb2? mkb2->Mk: 0)
      <<","<<nbonds_alive[ifc]
      <<","<<DemBondIfcAfaceh[ifc]
      <<","<<area_alive[ifc]
      <<","<<DemBondIfcIuuh[ifc]
      <<","<<iuu_alive[ifc]
      <<","<<DemBondIfcIvvh[ifc]
      <<","<<ivv_alive[ifc]
      <<","<<DemBondIfcIuvh[ifc]
      <<","<<iuv_alive[ifc]<<"\n";
  }
}

//==============================================================================
/// Initialises arrays and variables for the execution.
/// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphGpu::InitRunGpu(){
  unsigned nfilter=0;
  ParticlesDataDown(Np,0,true,false,NULL,nfilter,NULL,NULL,NULL,NULL);
  InitRun(Np,Idp_c->cptr(),AuxPos_c->cptr());
  if(DemStress_g)DemStress_g->CuMemset(0,Np);
  if(DemStrain_g)DemStrain_g->CuMemset(0,Np);
  if(DemShearStress_g)DemShearStress_g->CuMemset(0,Np);
  if(DemBond_g)DemBond_g->CuMemset(0,Np);
  if(CaseNfloat)InitFloatingsGpu(FtoMasspg,FtoDatpg,FtoCenterg,DemDatag,DemData2g,DemData3g);
  InitDemBondsGpu();
  if(TStep==STEP_Verlet)VelrhoM1_g->CuCopyFrom(Velrho_g,Np);
  if(TVisco==VISCO_LaminarSPS)SpsTauRho2_g->CuMemset(0,Np);
  if(MotionVel_g)MotionVel_g->CuMemset(0,Np); //<vs_m2dbc>
  if(MotionAce_g)MotionAce_g->CuMemset(0,Np); //<vs_m2dbc>
  if(ShiftVel_g)ShiftVel_g->CuMemset(0,Np);   //<vs_advshift>
  if(ShiftVel_g)FSType_g->CuMemset(3,Np);     //<vs_advshift>
  if(PsiClean_g)PsiClean_g->CuMemset(0,Np);   //<vs_divclean>
  Check_CudaErroor("Failed initializing variables for execution.");
}

//==============================================================================
/// Prepares variables for interaction (pre-loop and forces).
/// Prepara variables para interaccion (pre-loop y fuerzas).
//==============================================================================
void JSphGpu::PreInteraction_Forces(TpInterStep interstep){
  Timersg->TmStart(TMG_CfPreForces,false);
  //-Assign memory.
  ViscDt_g->Reserve();
  Ar_g->Reserve();
  Ace_g->Reserve();
  if(DDTArray)Delta_g->Reserve();
  if(Shifting)ShiftPosfs_g->Reserve();
  if(TVisco==VISCO_LaminarSPS)Sps2Strain_g->Reserve();
  if(TMdbc2==MDBC2_NoPen)NoPenShift_g->Reserve(); // SHABA
  if(ShiftingAdv){ //<vs_advshift_ini>
    FSMinDist_g->Reserve();
    FSNormal_g->Reserve();
  } //<vs_advshift_end>

  //<vs_divclean_ini>
  if(DivCleaning){
    PsiCleanRhs_g->Reserve(); 
    CsPsiClean_g->Reserve();
    PsiCleanRhs_g->CuMemset(0,Np);
    CsPsiClean_g->CuMemset(0,Np);
    CsPsiCleanMax=0;
  }  
  //vs_divclean_end>

  //-Initialise arrays.
  const unsigned npf=Np-Npb;
  ViscDt_g->CuMemset(0,Np);                                         //ViscDtg[]=0
  Ar_g->CuMemset(0,Np);                                             //Arg[]=0
  Ace_g->CuMemset(0,Np);                                            //Aceg[]=(0)
  if(DemStress_g)DemStress_g->CuMemset(0,Np);                       //DemStressg[]=0
  if(DemStrain_g)DemStrain_g->CuMemset(0,Np);                       //DemStraing[]=0
  if(DemShearStress_g)DemShearStress_g->CuMemset(0,Np);             //DemShearStressg[]=0
  if(DemBond_g)DemBond_g->CuMemset(0,Np);                           //DemBondg[]=0
  if(AG_CPTR(Delta_g))Delta_g->CuMemset(0,Np);                      //Deltag[]=0
  if(AG_CPTR(Sps2Strain_g))Sps2Strain_g->CuMemsetOffset(Npb,0,npf); //Sps2Straing[]=(0)
  if(AG_CPTR(NoPenShift_g))NoPenShift_g->CuMemset(0,Np);            //NoPenShiftg[]=(0) //<vs_m2dbcNP
  
  //-Select particles for shifting.
  if(AC_CPTR(ShiftPosfs_g))Shifting->InitGpu(npf,Npb,Posxy_g->cptr()
    ,Posz_g->cptr(),ShiftPosfs_g->ptr());

  //<vs_advshift_ini>
  if(AC_CPTR(ShiftVel_g) && interstep==INTERSTEP_SymPredictor)
    ShiftVel_g->CuMemset(0,Np);
  if(AC_CPTR(FSMinDist_g))FSMinDist_g->CuMemset(0,Np);
  if(AC_CPTR(FSNormal_g))FSNormal_g->CuMemset(0,Np);
  //<vs_advshift_end>

  //-Adds variable acceleration from input configuration.
  if(AccInput)AccInput->RunGpu(TimeStep,Gravity,npf,Npb,Code_g->cptr()
    ,Posxy_g->cptr(),Posz_g->cptr(),Velrho_g->cptr(),Ace_g->ptr());

  //-Computes VelMax: Includes the particles from floating bodies and does not affect the periodic conditions.
  //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
  const unsigned pini=(DtAllParticles? 0: Npb);
  cusph::ComputeVelMod(Np-pini,Velrho_g->cptr()+pini,ViscDt_g->ptr());
  const float velmax=cusph::ReduMaxFloat(Np-pini,0,ViscDt_g->ptr()
    ,CellDiv->GetAuxMem(cusph::ReduMaxFloatSize(Np-pini)));
  VelMax=sqrt(velmax);

  ViscDt_g->CuMemset(0,Np); //ViscDtg[]=0
  ViscDtMax=0;
  //<vs_flexstruc_ini>
  if(FlexStruc){
    cudaMemset(FlexStrucDtg,0,sizeof(float)*CaseNflexstruc);  //FlexStrucDtg[]=0
    FlexStrucDtMax=0;
  }
  //<vs_flexstruc_end>
  Timersg->TmStop(TMG_CfPreForces,true);
  Check_CudaErroor("Failed calculating VelMax.");
}

//==============================================================================
/// Frees memory allocated from Arrays_Gpu in PreInteraction_Forces().
/// Libera memoria asignada de Arrays_Gpu en PreInteraction_Forces().
//==============================================================================
void JSphGpu::PosInteraction_Forces(){
  ViscDt_g->Free();
  Ar_g->Free();
  Ace_g->Free();
  Delta_g->Free();
  ShiftPosfs_g->Free();
  if(Sps2Strain_g)Sps2Strain_g->Free();
  if(BoundMode_g)BoundMode_g->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
  if(TangenVel_g)TangenVel_g->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
  if(NoPenShift_g)NoPenShift_g->Free(); //<vs_m2dbcNP>
  if(FSMinDist_g)FSMinDist_g->Free(); //<vs_advshift>
  if(FSNormal_g)FSNormal_g->Free();   //<vs_advshift>
  if(PsiCleanRhs_g)PsiCleanRhs_g->Free();  //<vs_divclean>
  if(CsPsiClean_g)CsPsiClean_g->Free() ;   //<vs_divclean>
}

//==============================================================================
/// Updates particles according to forces and dt using Verlet.
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphGpu::ComputeVerlet(double dt){  //pdtedom
  Timersg->TmStart(TMG_SuComputeStep,false);
  const bool shift=(ShiftingMode!=SHIFT_None);
  const bool inout=(InOut!=NULL);
  const float3* indirvel=(inout? InOut->GetDirVelg(): NULL);
  const byte* boundmode=AG_CPTR(BoundMode_g); //<vs_m2dbc>
  VerletStep++;
  //-Allocates memory to compute the displacement.
  agdouble2 movxyg("movxyg",Arrays_Gpu,true);
  agdouble movzg("movzg",Arrays_Gpu,true);
  //-Computes displacement, velocity and density.
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    cusphs::ComputeStepVerlet(WithFloating,shift,inout,TMdbc2,Np,Npb
      ,Velrho_g->cptr(),VelrhoM1_g->cptr(),boundmode,Ar_g->cptr()
      ,Ace_g->cptr(),ShiftPosfs_g->cptr(),indirvel,AG_CPTR(NoPenShift_g),dt,twodt
      ,RhopZero,RhopOutMin,RhopOutMax,Gravity,Code_g->ptr()
      ,movxyg.ptr(),movzg.ptr(),VelrhoM1_g->ptr(),NULL);
  }
  else{
    cusphs::ComputeStepVerlet(WithFloating,shift,inout,TMdbc2,Np,Npb
      ,Velrho_g->cptr(),Velrho_g->cptr(),boundmode,Ar_g->cptr()
      ,Ace_g->cptr(),ShiftPosfs_g->cptr(),indirvel,AG_CPTR(NoPenShift_g),dt,dt
      ,RhopZero,RhopOutMin,RhopOutMax,Gravity,Code_g->ptr()
      ,movxyg.ptr(),movzg.ptr(),VelrhoM1_g->ptr(),NULL);
    VerletStep=0;
  }
  //-Applies displacement to non-periodic fluid particles.
  cusph::ComputeStepPos(PeriActive,WithFloating,Np,Npb,movxyg.cptr(),movzg.cptr()
    ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    #ifndef MDBC2_KEEPVEL
      if(MotionVel_g)cusphs::CopyMotionVelFlexStruc(CaseNflexstruc,Code_g->cptr(),FlexStrucRidpg,MotionVel_g->cptr(),Velrho_g->ptr());
    #endif
    cusphs::ComputeStepFlexStrucSemiImplicitEuler(CaseNflexstruc,Velrho_g->cptr(),Code_g->cptr(),FlexStrucRidpg,Ace_g->cptr(),dt,Gravity,movxyg.ptr(),movzg.ptr(),VelrhoM1_g->ptr(),NULL);
    cusph::ComputeStepPosFlexStruc(CaseNflexstruc,FlexStrucRidpg,Posxy_g->cptr(),Posz_g->cptr(),movxyg.cptr(),movzg.cptr(),Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
    if(!cusph::FlexStrucStepIsValid(Npb,Code_g->cptr()))
      Run_Exceptioon("Issue with FlexStruc particle update step (either naturally moved out of the domain or blew up).");
    BoundChanged=true;
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  //-The new values are calculated in VelrhoM1_g.
  Velrho_g->SwapPtr(VelrhoM1_g); //-Exchanges Velrho_g and VelrhoM1_g.

  Timersg->TmStop(TMG_SuComputeStep,true);
}

//==============================================================================
/// Updates particles according to forces and dt using Symplectic-Predictor.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphGpu::ComputeSymplecticPre(double dt){
  Timersg->TmStart(TMG_SuComputeStep,false);
  const bool shift=false; //(ShiftingMode!=SHIFT_None); //-We strongly recommend running the shifting correction only for the corrector. If you want to re-enable shifting in the predictor, change the value here to "true".
  const bool inout=(InOut!=NULL);
  //-Assign memory to PRE variables.
  PosxyPre_g->Reserve();
  PoszPre_g->Reserve();
  VelrhoPre_g->Reserve();
  if(DivCleaning)PsiCleanPre_g->Reserve();  //<vs_divclean>
  //-Move current data to PRE variables for calculating the new data.
  PosxyPre_g->SwapPtr(Posxy_g);   //- PosxyPre_g[]   <= Posxy_g[]
  PoszPre_g->SwapPtr(Posz_g);     //- PoszPre_g[]    <= Posz_g[]
  VelrhoPre_g->SwapPtr(Velrho_g); //- VelrhoPre_g[]  <= Velrho_g[]
  if(DivCleaning)PsiCleanPre_g->SwapPtr(PsiClean_g); //-PhiCleanPre[]<= PhiClean_g[] //<vs_divclean>
  //-Allocate memory to compute the diplacement.
  agdouble2 movxyg("movxyg",Arrays_Gpu,true);
  agdouble  movzg("movzg",Arrays_Gpu,true);
  //-Compute displacement, velocity and density.
  const double dt05=dt*.5;
  const float3* indirvel=(InOut? InOut->GetDirVelg(): NULL);
  const byte* boundmode=AG_CPTR(BoundMode_g); //<vs_m2dbc>
  cusphs::ComputeStepSymplecticPre(WithFloating,shift,inout,TMdbc2,Np,Npb
    ,VelrhoPre_g->cptr(),boundmode,Ar_g->cptr(),Ace_g->cptr(),ShiftPosfs_g->cptr()
    ,indirvel,dt05,RhopZero,RhopOutMin,RhopOutMax,Gravity
    ,Code_g->ptr(),movxyg.ptr(),movzg.ptr(),Velrho_g->ptr()
    ,AG_PTR(PsiClean_g),AG_CPTR(PsiCleanPre_g),AG_CPTR(PsiCleanRhs_g),DivCleaning,NULL);  //<vs_divclean>

  //-Applies displacement to non-periodic fluid particles.
  cusph::ComputeStepPos2(PeriActive,WithFloating,Np,Npb
    ,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr()
    ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
  //-Copy previous position of the boundary particles.
  Posxy_g->CuCopyFrom(PosxyPre_g,Npb);
  Posz_g->CuCopyFrom(PoszPre_g,Npb);

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    #ifndef MDBC2_KEEPVEL
      if(MotionVel_g)cusphs::CopyMotionVelFlexStruc(CaseNflexstruc,Code_g->cptr(),FlexStrucRidpg,MotionVel_g->cptr(),VelrhoPre_g->ptr());
    #endif
    cusphs::ComputeStepFlexStrucSymplecticPre(CaseNflexstruc,VelrhoPre_g->cptr(),Code_g->cptr(),FlexStrucRidpg,Ace_g->cptr(),dt05,Gravity,movxyg.ptr(),movzg.ptr(),Velrho_g->ptr(),NULL);
    cusph::ComputeStepPosFlexStruc(CaseNflexstruc,FlexStrucRidpg,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr(),Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
    if(!cusph::FlexStrucStepIsValid(Npb,Code_g->cptr()))
      Run_Exceptioon("Issue with FlexStruc particle update step (either naturally moved out of the domain or blew up).");
    BoundChanged=true;
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  Timersg->TmStop(TMG_SuComputeStep,false);
}

//==============================================================================
/// Updates particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphGpu::ComputeSymplecticCorr(double dt){
  Timersg->TmStart(TMG_SuComputeStep,false);
  const bool shift=(ShiftingMode!=SHIFT_None);
  const bool inout=(InOut!=NULL);
  const bool shiftadv=(ShiftingAdv!=NULL); //<vs_advshift>
  //-Allocate memory to compute the diplacement.
  agdouble2 movxyg("movxyg",Arrays_Gpu,true);
  agdouble  movzg("movzg",Arrays_Gpu,true);
  //-Computes displacement, velocity and density.
  const double dt05=dt*.5;
  const float3* indirvel=(InOut? InOut->GetDirVelg(): NULL);
  const byte*   boundmode=AG_CPTR(BoundMode_g); //<vs_m2dbc>
  const float4* shiftvel=AG_CPTR(ShiftVel_g); //<vs_advshift>
  cusphs::ComputeStepSymplecticCor(WithFloating,shift,shiftadv,inout,TMdbc2,Np,Npb
    ,VelrhoPre_g->cptr(),boundmode,Ar_g->cptr(),Ace_g->cptr(),ShiftPosfs_g->cptr()
    ,indirvel,AG_CPTR(NoPenShift_g),shiftvel,dt05,dt,RhopZero,RhopOutMin,RhopOutMax,Gravity // SHABA
    ,Code_g->ptr(),movxyg.ptr(),movzg.ptr(),Velrho_g->ptr()
    ,AG_PTR(PsiClean_g),AG_CPTR(PsiCleanPre_g),AG_CPTR(PsiCleanRhs_g),DivCleaning,NULL);  //<vs_divclean>

  //-Applies displacement to non-periodic fluid particles.
  cusph::ComputeStepPos2(PeriActive,WithFloating,Np,Npb
    ,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr()
    ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    cusphs::ComputeStepFlexStrucSymplecticCor(CaseNflexstruc,VelrhoPre_g->cptr(),Code_g->cptr(),FlexStrucRidpg,Ace_g->cptr(),dt05,dt,Gravity,movxyg.ptr(),movzg.ptr(),Velrho_g->ptr(),NULL);
    cusph::ComputeStepPosFlexStruc(CaseNflexstruc,FlexStrucRidpg,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr(),Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
    if(!cusph::FlexStrucStepIsValid(Npb,Code_g->cptr()))
      Run_Exceptioon("Issue with FlexStruc particle update step (either naturally moved out of the domain or blew up).");
    BoundChanged=true;
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  //-Free memory assigned to PRE variables in ComputeSymplecticPre().
  PosxyPre_g->Free();
  PoszPre_g->Free();
  VelrhoPre_g->Free();
  if(DivCleaning)PsiCleanPre_g->Free();  //<vs_divclean>
  Timersg->TmStop(TMG_SuComputeStep,true);
}

//==============================================================================
/// Computes the variable Dt.
/// Calcula un Dt variable.
//==============================================================================
double JSphGpu::DtVariable(bool final){
  //-dt1 depends on force per unit mass.
  const double acemax=sqrt(double(AceMax));
  const double dt1=(AceMax? (sqrt(double(KernelH)/AceMax)): DBL_MAX); 
  //-dt2 combines the Courant and the viscous time-step controls.
  const double dt2=double(KernelH)/(max(Cs0,VelMax*10.)+double(KernelH)*ViscDtMax);
  //-dt3 uses the maximum speed of sound across all structure particles.
  const double dt3=(FlexStrucDtMax? double(KernelH)/FlexStrucDtMax: DBL_MAX); //<vs_flexstruc>
  //-dt new value of time step.
  #ifdef AVAILABLE_DIVCLEAN
    //-dt4 uses the maximum local speed of sound across all fluid particles.
    const double dt4=(DivCleaning? double(KernelH)/CsPsiCleanMax: DBL_MAX); //<vs_divclean>
    double dt=CFLnumber*min(dt1,min(dt2,min(dt3,dt4)));
  #else
  double dt=CFLnumber*min(dt1,min(dt2,dt3));
  #endif
  //-Use dt value defined by FixedDt object.
  if(FixedDt)dt=FixedDt->GetDt(TimeStep,dt);
  //-Check that the dt is valid.
  if(fun::IsNAN(dt) || fun::IsInfinity(dt)){
    if(FlexStruc)Run_Exceptioon(fun::PrintStr(
    "The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f, FlexStrucDtMax=%f) is NaN or infinity at nstep=%u."
    ,dt,AceMax,VelMax,ViscDtMax,FlexStrucDtMax,Nstep));
    else         Run_Exceptioon(fun::PrintStr(
    "The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f) is NaN or infinity at nstep=%u."
    ,dt,AceMax,VelMax,ViscDtMax,Nstep));
  }
  //-Correct dt with minimum allowed value.
  if(dt<double(DtMin)){ 
    dt=double(DtMin); DtModif++;
    if(DtModif>=DtModifWrn){
      Log->PrintfWarning("%d DTs adjusted to DtMin (t:%g, nstep:%u)",DtModif,TimeStep,Nstep);
      DtModifWrn*=10;
    }
  }
  //-Saves information about dt.
  if(final){
    if(PartDtMin>dt)PartDtMin=dt;
    if(PartDtMax<dt)PartDtMax=dt;
    //-Saves detailed information about dt in SaveDt object.
    if(SaveDt)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber
      ,dt3*CFLnumber,AceMax,ViscDtMax,FlexStrucDtMax,VelMax);
  }
  return(dt);
}

//==============================================================================
/// Computes final shifting distance for the particle position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphGpu::RunShifting(double dt){
  Timersg->TmStart(TMG_SuShifting,false);
  Shifting->RunGpu(Np-Npb,Npb,dt,Velrho_g->cptr(),ShiftPosfs_g->ptr());
  Timersg->TmStop(TMG_SuShifting,true);
}

//==============================================================================
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
void JSphGpu::CalcMotion(double stepdt){
  Timersg->TmStart(TMG_SuMotion,false);
  JSph::CalcMotion(stepdt);
  Timersg->TmStop(TMG_SuMotion,false);
}

//==============================================================================
/// Processes boundary particle movement.
/// Procesa movimiento de boundary particles.
//==============================================================================
void JSphGpu::RunMotion(double stepdt){
  Timersg->TmStart(TMG_SuMotion,false);
  BoundChanged=false;
  //-Add motion from automatic wave generation.
  if(WaveGen)CalcMotionWaveGen(stepdt);
  //-Process particles motion.
  if(DsMotion->GetActiveMotion()){
    BoundChanged=true;
    const unsigned nref=DsMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      if(m.type==MOTT_Linear){//-Linear movement.
        //Log->PrintfDbg("%u] t:%g  dt:%g  mov:(%g,%g,%g)  vel:(%g,%g,%g)"
        //  ,Nstep,TimeStep,stepdt,m.linmov.x,m.linmov.y,m.linmov.z,m.linvel.x,m.linvel.y,m.linvel.z);
        cusph::MoveLinBound(PeriActive,m.count,m.idbegin-CaseNfixed,m.linmov
          ,ToTFloat3(m.linvel),RidpMotg,Posxy_g->ptr(),Posz_g->ptr()
          ,Dcell_g->ptr(),Velrho_g->ptr(),Code_g->ptr());
      }
      if(m.type==MOTT_Matrix){//-Matrix movement (for rotations).
        cusph::MoveMatBound(PeriActive,Simulate2D,m.count,m.idbegin-CaseNfixed
          ,m.matmov,stepdt,RidpMotg,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr()
          ,Velrho_g->ptr(),Code_g->ptr(),AG_PTR(BoundNor_g));
      }      
    }
  }
  //-Management of Multi-Layer Pistons.
  if(MLPistons){
    BoundChanged=true;
    if(MLPistons->GetPiston1dCount()){//-Process motion for pistons 1D.
      MLPistons->CalculateMotion1d(TimeStep+MLPistons->GetTimeMod()+stepdt);
      cusph::MovePiston1d(PeriActive!=0,CaseNmoving,0,Dp,MLPistons->GetPoszMin()
        ,MLPistons->GetPoszCount(),MLPistons->GetPistonIdGpu()
        ,MLPistons->GetMovxGpu(),MLPistons->GetVelxGpu(),RidpMotg
        ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Velrho_g->ptr()
        ,Code_g->ptr());
    }
    for(unsigned cp=0;cp<MLPistons->GetPiston2dCount();cp++){//-Process motion for pistons 2D.
      const JMLPistons::StMotionInfoPiston2D mot=MLPistons->CalculateMotion2d(cp
        ,TimeStep+MLPistons->GetTimeMod()+stepdt);
      cusph::MovePiston2d(PeriActive!=0,mot.np,mot.idbegin-CaseNfixed,Dp
        ,mot.posymin,mot.poszmin,mot.poszcount,mot.movyz,mot.velyz,RidpMotg
        ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Velrho_g->ptr()
        ,Code_g->ptr());
    }
  }
  //<vs_m2dbc_ini>
  //-Copy motion velocity and compute acceleration of moving particles.
  if(MotionVel_g){
    cusph::CopyMotionVelAce(CaseNmoving,stepdt,RidpMotg,Velrho_g->cptr()
      ,MotionVel_g->ptr(),MotionAce_g->ptr());
  } //<vs_m2dbc_end>
  Timersg->TmStop(TMG_SuMotion,true);
}

//==============================================================================
/// Applies RelaxZone to selected particles.
/// Aplica RelaxZone a las particulas indicadas.
//==============================================================================
void JSphGpu::RunRelaxZone(double dt){
  Timersg->TmStart(TMG_SuMotion,false);
  agbyte  rzid("rzid",Arrays_Gpu,false);
  agfloat rzfactor("rzfactor",Arrays_Gpu,false);
  RelaxZones->SetFluidVelGpu(TimeStep,dt,Np-Npb,Npb,(const tdouble2*)Posxy_g->cptr()
    ,Posz_g->cptr(),Idp_g->cptr(),(tfloat4*)Velrho_g->ptr(),rzid.ptr(),rzfactor.ptr());
  Timersg->TmStop(TMG_SuMotion,true);
}

//==============================================================================
/// Applies Damping to indicated particles.
/// Aplica Damping a las particulas indicadas.
//==============================================================================
void JSphGpu::RunDamping(double dt){
  const typecode* codeptr=(CaseNfloat || PeriActive? Code_g->cptr(): NULL);
  Damping->ComputeDampingGpu(TimeStep,dt,Np-Npb,Npb,Posxy_g->cptr(),Posz_g->cptr()
    ,codeptr,Velrho_g->ptr());
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::SaveVtkNormalsGpu(std::string filename,int numfile,unsigned np,unsigned npb
  ,const double2* posxyg,const double* poszg,const unsigned* idpg,const float3* boundnorg)
{
  //-Allocates memory (simple way for debug method).
  const unsigned n=(UseNormalsFt? np: npb);
  tdouble3* pos=fcuda::ToHostPosd3(0,n,posxyg,poszg);
  unsigned* idp=fcuda::ToHostUint(0,n,idpg);
  tfloat3*  nor=fcuda::ToHostFloat3(0,n,boundnorg);
  //-Generates VTK file.
  SaveVtkNormals(filename,numfile,np,npb,pos,idp,nor,1.f);
  //-Frees memory.
  delete[] pos;
  delete[] idp;
  delete[] nor;
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const double2* posxyg,const double* poszg
  ,const typecode* codeg,const unsigned* idpg,const float4* velrhog,const float3* fsnormal)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3* posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode* codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode* codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned* idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4*  velrhoh=fcuda::ToHostFloat4(pini,np,velrhog);
  tfloat3*  fsnormalh=fcuda::ToHostFloat3(pini,np,fsnormal);
  //-Creates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,np,posh,codeh,idph,velrhoh,fsnormalh);
  //-Deallocates memory.
  delete[] posh;  
  delete[] codeh;
  delete[] idph;
  delete[] velrhoh;
  delete[] fsnormalh;
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const double2* posxyg,const double* poszg
  ,const typecode* codeg,const unsigned* idpg,const float4* velrhog)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3* posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode* codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode* codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned* idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4*  velrhoh=fcuda::ToHostFloat4(pini,np,velrhog);
  //-Creates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,np,posh,codeh,idph,velrhoh);
  //-Deallocates memory.
  delete[] posh;  
  delete[] codeh;
  delete[] idph;
  delete[] velrhoh;
}

//==============================================================================
/// Saves VTK file with particle data (debug).
/// Graba fichero vtk con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,unsigned cellcode,const double2* posxyg
  ,const double* poszg,const unsigned* idpg,const unsigned* dcelg
  ,const typecode* codeg,const float4* velrhog,const float4* velrhom1g
  ,const float3* aceg)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  //-Loads position.
  tfloat3* pos=new tfloat3[n];
  {
    tdouble2* pxy=new tdouble2[n];
    double* pz=new double[n];
    cudaMemcpy(pxy,posxyg+pini,sizeof(double2)*n,cudaMemcpyDeviceToHost);
    cudaMemcpy(pz,poszg+pini,sizeof(double)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)pos[p]=TFloat3(float(pxy[p].x),float(pxy[p].y),float(pz[p]));
    delete[] pxy;
    delete[] pz;
  }
  //-Loads idp.
  unsigned* idp=NULL;
  if(idpg){
    idp=new unsigned[n];
    cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  }
  //-Loads dcel.
  tuint3* dcel=NULL;
  if(dcelg){
    dcel=new tuint3[n];
    unsigned* aux=new unsigned[n];
    cudaMemcpy(aux,dcelg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)dcel[p]=TUint3(unsigned(DCEL_Cellx(cellcode,aux[p])),unsigned(DCEL_Celly(cellcode,aux[p])),unsigned(DCEL_Cellz(cellcode,aux[p])));
    delete[] aux;
  }
  //-Loads vel and rho.
  tfloat3* vel=NULL;
  float* rho=NULL;
  if(velrhog){
    vel=new tfloat3[n];
    rho=new float[n];
    tfloat4* aux=new tfloat4[n];
    cudaMemcpy(aux,velrhog+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ vel[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rho[p]=aux[p].w; }
    delete[] aux;
  }
  //-Loads velm1 and rhom1.
  tfloat3* velm1=NULL;
  float* rhom1=NULL;
  if(velrhom1g){
    velm1=new tfloat3[n];
    rhom1=new float[n];
    tfloat4* aux=new tfloat4[n];
    cudaMemcpy(aux,velrhom1g+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ velm1[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rhom1[p]=aux[p].w; }
    delete[] aux;
  }
  //-Loads ace.
  tfloat3* ace=NULL;
  if(aceg){
    ace=new tfloat3[n];
    cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  }
  //-Loads type.
  byte* type=NULL;
  if(codeg){
    type=new byte[n];
    typecode* aux=new typecode[n];
    cudaMemcpy(aux,codeg+pini,sizeof(typecode)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ 
      const typecode cod=aux[p];
      byte tp=99;
      if(CODE_IsFixed(cod))tp=0;
      else if(CODE_IsMoving(cod))tp=1;
      else if(CODE_IsFloating(cod))tp=2;
      else if(CODE_IsFluid(cod))tp=3;
      if(CODE_IsNormal(cod))tp+=0;
      else if(CODE_IsPeriodic(cod))tp+=10;
      else if(CODE_IsOutIgnore(cod))tp+=20;
      else tp+=30;
      type[p]=tp;
    }
    delete[] aux;
  }
  //-Saves VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos",n,pos,true);
  if(idp)   arrays.AddArray("Idp"  ,n,idp  ,true);
  if(dcel)  arrays.AddArray("Dcel" ,n,dcel ,true);
  if(vel)   arrays.AddArray("Vel"  ,n,vel  ,true);
  if(rho)   arrays.AddArray("Rhop" ,n,rho  ,true);
  if(velm1) arrays.AddArray("Velm1",n,velm1,true);
  if(rhom1) arrays.AddArray("Rhom1",n,rhom1,true);
  if(ace)   arrays.AddArray("Ace"  ,n,ace  ,true);
  if(type)  arrays.AddArray("Typex",n,type ,true);
  JSpVtkData::Save(fun::FileNameSec(filename,numfile),arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,bool idp,bool vel,bool rho,bool code)
{
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  if(Log->GetMpiRank()>=0)filename=string("p")+fun::IntStr(Log->GetMpiRank())+"_"+filename;
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned nfilter=0;
  ParticlesDataDown(n,pini,code,false,NULL,nfilter,NULL,NULL,NULL,NULL);
  tfloat3* pos=new tfloat3[n];
  for(unsigned p=0;p<n;p++)pos[p]=ToTFloat3(AuxPos_c->cptr()[p]);
  byte* type=new byte[n];
  for(unsigned p=0;p<n;p++){
    const typecode cod=Code_c->cptr()[p];
    byte tp=99;
    if(CODE_IsFixed(cod))tp=0;
    else if(CODE_IsMoving(cod))tp=1;
    else if(CODE_IsFloating(cod))tp=2;
    else if(CODE_IsFluid(cod))tp=3;
    if(CODE_IsNormal(cod))tp+=0;
    else if(CODE_IsPeriodic(cod))tp+=10;
    else if(CODE_IsOutIgnore(cod))tp+=20;
    else tp+=30;
    type[p]=tp;
  }
  //-Saves VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos",n,pos,true);
  if(idp) arrays.AddArray("Idp"  ,n,Idp_c->cptr()   ,false);
  if(type)arrays.AddArray("Typex",n,type            ,true);
  if(vel) arrays.AddArray("Vel"  ,n,AuxVel_c->cptr(),false);
  if(rho) arrays.AddArray("Rho"  ,n,AuxRho_c->cptr(),false);
#ifdef CODE_SIZE4
  if(code)arrays.AddArray("Code" ,n,(unsigned*)Code_c->cptr(),false);
#else
  if(code)arrays.AddArray("Code" ,n,(word*)    Code_c->cptr(),false);
#endif
  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const float3* posg,const byte* checkg
  ,const unsigned* idpg,const float3* velg,const float* rhog)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  tfloat3*  pos=new tfloat3[n];
  byte*     check=NULL;
  unsigned* idp=NULL;
  tfloat3*  vel=NULL;
  float*    rho=NULL;
  if(checkg)check=new byte[n];
  if(idpg)  idp=new unsigned[n];
  if(velg)  vel=new tfloat3[n];
  if(rhog)  rho=new float[n];
  //-Copies data from GPU.
  cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(checkg)cudaMemcpy(check,checkg+pini,sizeof(byte)*n,cudaMemcpyDeviceToHost);
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhog)cudaMemcpy(rho,rhog+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  //-Generates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,n,pos,check,idp,vel,rho);
  //-Frees memory.
  delete[] pos;
  delete[] check;
  delete[] idp;
  delete[] vel;
  delete[] rho;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const float3* posg
  ,const unsigned* idpg,const float3* velg,const float* rhog,const float* arg
  ,const float3* aceg,const float3* vcorrg)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned* idp=NULL;   if(idpg)idp=new unsigned[n];
  tfloat3*  pos=NULL;   if(posg)pos=new tfloat3[n];
  tfloat3*  vel=NULL;   if(velg)vel=new tfloat3[n];
  float*    rho=NULL;   if(rhog)rho=new float[n];
  float*    ar=NULL;    if(arg) ar=new float[n];
  tfloat3*  ace=NULL;   if(aceg)ace=new tfloat3[n];
  tfloat3*  vcorr=NULL; if(vcorrg)vcorr=new tfloat3[n];
  //-Copies data from GPU.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhog)cudaMemcpy(rho,rhog+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(arg)cudaMemcpy(ar,arg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(aceg)cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(vcorrg)cudaMemcpy(vcorr,vcorrg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying data from GPU.");
  //-Generates CSV file.
  DgSaveCsvParticlesCpu(filename,numfile,0,n,head,pos,idp,vel,rho,ar,ace,vcorr);
  //-Frees memory.
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rho;
  delete[] ar;
  delete[] ace;
  delete[] vcorr;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu2(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const float3* posg
  ,const unsigned* idpg,const float3* velg,const float* rhog
  ,const float4* pospresg,const float4* velrhog)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned* idp=NULL;     if(idpg)idp=new unsigned[n];
  tfloat3*  pos=NULL;     if(posg)pos=new tfloat3[n];
  tfloat3*  vel=NULL;     if(velg)vel=new tfloat3[n];
  float*    rho=NULL;     if(rhog)rho=new float[n];
  tfloat4*  pospres=NULL; if(pospresg)pospres=new tfloat4[n];
  tfloat4*  velrho=NULL;  if(velrhog)velrho=new tfloat4[n];
  //-Copies data from GPU.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhog)cudaMemcpy(rho,rhog+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(pospresg)cudaMemcpy(pospres,pospresg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  if(velrhog)cudaMemcpy(velrho,velrhog+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying data from GPU.");
  //-Generates CSV file.
  DgSaveCsvParticles2(filename,numfile,0,n,head,pos,idp,vel,rho,pospres,velrho);
  //-Frees memory.
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rho;
  delete[] pospres;
  delete[] velrho;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticles2(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const tfloat3* pos
  ,const unsigned* idp,const tfloat3* vel,const float* rho
  ,const tfloat4* pospres,const tfloat4* velrho)
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Generates CSV file.
  ofstream pf;
  pf.open(filename.c_str());
  if(pf){
    if(!head.empty())pf << head << endl;
    pf << "Num";
    if(idp)pf << ";Idp";
    if(pos)pf << ";PosX;PosY;PosZ";
    if(vel)pf << ";VelX;VelY;VelZ";
    if(rho)pf << ";Rhop";
    if(pospres)pf << ";Px;Py;Pz;Pp";
    if(velrho) pf << ";Vx;Vy;Vz;Vr";
    pf << endl;
    const char fmt1[]="%f"; //="%24.16f";
    const char fmt3[]="%f;%f;%f"; //="%24.16f;%24.16f;%24.16f";
    for(unsigned p=pini;p<pfin;p++){
      pf << fun::UintStr(p-pini);
      if(idp)pf << ";" << fun::UintStr(idp[p]);
      if(pos)pf << ";" << fun::Float3Str(pos[p],fmt3);
      if(vel)pf << ";" << fun::Float3Str(vel[p],fmt3);
      if(rho)pf << ";" << fun::FloatStr(rho[p],fmt1);
      if(pospres)pf << ";" <<  fun::FloatStr(pospres[p].x,fmt1) << ";" << fun::FloatStr(pospres[p].y,fmt1) << ";" << fun::FloatStr(pospres[p].z,fmt1) << ";" << fun::FloatStr(pospres[p].w,fmt1);
      if(velrho)pf << ";" <<  fun::FloatStr(velrho[p].x,fmt1) << ";" << fun::FloatStr(velrho[p].y,fmt1) << ";" << fun::FloatStr(velrho[p].z,fmt1) << ";" << fun::FloatStr(velrho[p].w,fmt1);
      pf << endl;
    }
    if(pf.fail())Run_ExceptioonFile("Failed writing to file.",filename);
    pf.close();
  }
  else Run_ExceptioonFile("File could not be opened.",filename);
}

const unsigned JSphGpu::MaxNStmFloatings;


#ifdef _WITHMR //<vs_vrres_ini>
void JSphGpu::DgSaveVtkParticlesGpuMultiRes(std::string filename,int numfile,unsigned pini,unsigned pfin
 ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg
 ,const float4 *velrhopg,const double* rcond,const float4 *shiftposg)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3 *posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode *codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode *codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned *idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4  *velrhoph=fcuda::ToHostFloat4(pini,np,velrhopg);
  tfloat4  *shiftpos=NULL;
  if(shiftposg)shiftpos=fcuda::ToHostFloat4(pini,np,shiftposg);
  double   *rcondh=NULL;
  if(rcond)rcondh=fcuda::ToHostDouble(pini,np,rcond);
  //-Creates VTK file.
  DgSaveVtkParticlesCpuMR(filename+"Fluid.vtk",numfile,0,np,posh,codeh,idph,velrhoph,NULL,NULL,NULL,shiftpos);
  DgSaveVtkParticlesCpuMRBuffer(filename+"Buffer.vtk",numfile,0,np,posh,codeh,idph,velrhoph,NULL,rcondh,shiftpos);
  DgSaveVtkParticlesCpuMRBound(filename+"Bound.vtk",numfile,0,np,posh,codeh,idph,velrhoph,NULL);
  //-Deallocates memory.
  delete[] posh;
  delete[] codeh;
  delete[] idph;
  delete[] velrhoph;
  delete[] rcondh;
}
#endif  //<vs_vrres_end>



