/**************************************************************************
						UDF of ozone deposition
@author:Jialei Shen
@e-mail:shenjialei1992@163.com
@latest:2016.10.20
This is an UDF file to simulate indoor ozone deposition in CFD, excluded 
the chemical reaction. The UDF file includes the following terms: 
1  Ozone sink term near the surfaces;
**************************************************************************/

#include "udf.h"
#include "sg.h"

#define r 2.0e-5           //reaction probability (-)
#define v 360.0            //Boltzmann velocity (m/s)
#define vs (r*v/4)         //surface uptake velocity (m/s)
#define Dm 1.82e-5         //diffusion coefficient of ozone in air (m2/s)
#define rho 1.205          //density of air (kg/m3)

/********************ozone sink term (ozone deposition)********************/

DEFINE_SOURCE(ozone_deposition_udf,c,t,dS,eqn)
{
	face_t f;
	Thread *tf;
	int n;
	real NV_VEC(A);
	real xc[ND_ND], xf[ND_ND],y0[ND_ND];
	real source,depo_rate;
	real dy0;
	C_CENTROID(xc,c,t);
	source=0.0;
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