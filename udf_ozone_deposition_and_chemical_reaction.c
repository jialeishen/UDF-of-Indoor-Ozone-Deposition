/**************************************************************************
				UDF of ozone deposition & chemical reaction
@author:Jialei Shen
@e-mail:shenjialei1992@163.com
@latest:2016.10.20
This is an UDF file to simulate indoor ozone deposition in CFD, included 
the chemical reaction. The UDF file includes the following terms: 
1  Ozone sink term;
2  B source&sink term;
3  P source term;
**************************************************************************/

#include "udf.h"
#include "sg.h"

#define r 2.0e-5           //reaction probability (-)
#define v 360.0            //Boltzmann velocity (m/s)
#define vs (r*v/4)         //surface uptake velocity (m/s)
#define Dm 1.82e-5         //diffusion coefficient of ozone in air (m2/s)
#define rho 1.205          //density of air (kg/m3)
#define kb_A 5111.111       //second order rate constant for ozone (kg/m3s) (convert from 0.0184ppb-1h-1) 98000.1978
#define kb_B 5111.111        //second order rate constant for ozone (kg/m3s) (convert from 0.0184ppb-1h-1) 34527.4604
#define B_source 5.25e-10  //source of B (kg/m2s) (convert from 1.89mg/m2h)

/*********ozone sink term (ozone deposition & chemical reaction)**********/

DEFINE_SOURCE(ozone_sink_udf,c,t,dS,eqn)
{
	face_t f;
	Thread *tf;
	int n;
	real NV_VEC(A);
	real xc[ND_ND], xf[ND_ND],y0[ND_ND];
	real source,depo_rate;
	real dy0;
	C_CENTROID(xc,c,t);
	source=kb_A*C_YI(c,t,1)*rho;
	c_face_loop (c,t,n)
	{
		f=C_FACE(c,t,n);
		tf=C_FACE_THREAD(c,t,n);
		F_CENTROID(xf,f,tf);
		if (THREAD_TYPE(tf)==THREAD_F_WALL)
		{
			NV_VV(y0,=,xc,-,xf);
			dy0=NV_MAG(y0);
			F_AREA(A,f,tf);
			depo_rate=vs/(1+vs*dy0/Dm)*NV_MAG(A)/C_VOLUME(c,t);
			source+=depo_rate;
		}
	}
	source=-rho*source;
	dS[eqn]=source;
	source*=C_YI(c,t,0);
	return source;
}

/***********B source from floor and sink of chemical reaction*************/

DEFINE_SOURCE(B_source_sink_udf,c,t,dS,eqn)
{
	face_t f;
	Thread *tf;
	int n;
	real NV_VEC(A);
	real source;
	real sourceB;
	real xf[ND_ND];
	C_CENTROID(xf,c,t);
	source=-kb_B*C_YI(c,t,0)*rho*C_YI(c,t,1)*rho;

	c_face_loop (c,t,n)
	{
		f=C_FACE(c,t,n);
		tf=C_FACE_THREAD(c,t,n);
		//F_CENTROID(xf,f,tf);
		if (THREAD_TYPE(tf)==THREAD_F_WALL && xf[1]<=0.00058)
		{
			F_AREA(A,f,tf);
			sourceB=B_source*NV_MAG(A)/C_VOLUME(c,t);
			source+=sourceB;
		}
	}
	dS[eqn]=-kb_B*C_YI(c,t,0)*rho*rho;
	return source;
}

/*********************P source of chemical reaction***********************/

DEFINE_SOURCE(P_source_udf,c,t,dS,eqn)
{
	face_t f;
	Thread *tf;
	real source;

	source=kb_B*C_YI(c,t,0)*rho*C_YI(c,t,1)*rho;
	dS[eqn]=0;
	return source;
}