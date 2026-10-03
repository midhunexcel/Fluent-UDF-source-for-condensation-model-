# Fluent-UDF-source-for-condensation-model-
Development of condensation model for dehumdification encountered in wet coil heat exchangers. Source term implementation in the species, continuity and energy eq are done seperately. 
#include "udf.h"
#include "math.h"
#define A_cc 23.1964
#define B 3816.44
#define C -46.13
DEFINE_SOURCE(condensation_continuity, c,t,dS,eqn)
{
/* variable declaration*/
/*real Sp = 0.0, Sc = 0.0;// initialisation NOT REQ IN CONTINUITY EQ*/
real source=0.0;
Thread *tf; // wall (face) -thread pointer represent the cold-wall under const Temp
face_t f; // int data type for face in particular face thread
cell_t c0; // int data type for cells in a particular cell thread
Thread *t0; // fluid thread for first layer of cells




/* Step 1: Identify if this cell belongs to the cold-wall-adjacent region */

Domain *domain = Get_Domain(1); // domain_ID =1 by default
int wall_1_id = 6, wall_2_id =7; // top and bottom symmetric cold wall id 8 & 9

/* wall zone 8 */

tf=Lookup_Thread(domain, wall_1_id);//GET THE THREAD OF WALL ID 8 REPRESENTING THE CONDENSATION SURFACE
begin_f_loop(f, tf) /* face loop over the wall_zone_ID 8 to extract the FIRST LAYER of CELLS*/
{
	c0=F_C0(f,tf);// go through the first layer faces of cells near wall, then get those cells and its thread into c0 and t0
	t0=THREAD_T0(tf);// ADJACENT CELL THREAD ON ONE SIDE
	if (c==c0 && t==t0)
	{
	    /* Step 2 : INCLUDE THE CONDITION FOR CONDENSATION */
	    					/* ie  Twall < DPT */
	   
	  /*cell centre values*/ 
	    real rho=C_R(c0,t0);// density of cell 
	    real Yc=C_YI(c0,t0,0);// w_v at cell centre 
		real Tc=C_T(c0,t0);// cell centre temperature 
		real R_g_v=461.504;// R_g/M_v chara gas constant J/kg.K
		real pv_Pa= Yc*rho*R_g_v*Tc;// partial pressure of water vapor in the cell [pa]
		/*Clausius Clapeyron correlation for DPT*/ 
		real pv_MPa = pv_Pa * 1.0e-6;
		real Tsat_C = B / (A_cc - log(pv_MPa)) - C;// T sat in Degree celcius
		real Tsat_K=Tsat_C+273.15;// DPT in [K]
		
		
		
		
		/* face values */
		real Tf=F_T(f,tf);// face MACROS to get the temperature of particular face centre of the cell c0, tf sharing with the wall 8
		
		/*implementation of the condition for condensation */
		
		if (Tf<Tsat_K)
		{
		real A [ND_ND];
		F_AREA(A,f,tf);// vector of the face 
		real area_mag= NV_MAG(A);// FACE AREA :magnitude of vector A
		real vol=C_VOLUME(c0,t0);// volume of each cell
		real Diff=0.26e-04;// diffusion coefficient [m2/s]
		
		
		
		
	 /* delta-n calculation- applicable for structured  and orthogonal mesh */
		real cell_cent[ND_ND];
		real face_cent[ND_ND];
		
		C_CENTROID(cell_cent,c0,t0);
		F_CENTROID(face_cent,f,tf);// tf instead of t0
		
		real dist=NV_D(cell_cent, face_cent);// delta-n calculation based on structured mesh
		
		
		/*source term definition */
		
		real Y_s=Yc; //OMEGA_C * i.e. the term that changes with the iterations by updating the new cell centre value from Yc
		real k;// variable clubbing term 
		k=(-1)*(rho*Diff*area_mag)/(dist*vol);// k with -ive sign 
		
		real T_w=F_T(f,tf);// wall temperature [K] of the interface of each cell c0,t0 with face f in tf
		real t_w_C=T_w-273.15;// T_w in  degree C
		real p_s_w=610.78 * exp((17.2694 * t_w_C) / (t_w_C + 238.3)); //partial pressure [Pa] corresponding to the saturation temp t_w_C
		real Yw= p_s_w/(rho*R_g_v*T_w);// wall mass fraction for the interface calculated from the equation
		//real Sp=k*((1-Yw)/((1-Y_s)*(1-Y_s)));// (dS/dw_C)wc*= Sp
		real S_s=k*((Y_s-Yw)/(1-Y_s)); // ie S* IN THE EXPRESSION
		/*Sc=S_s-(Sp*Y_s); S*-()  only Species UDF requires Sp, Sc, but here we only need S_s.*/
		/* simplification of source here as here there is no requirememt of Sc+Sp phi format source = S_s after simplification*/
		source=S_s;// TOTAL MASS OF WATER CONDENSATION RATE PER UNIT VOLUME kg/s-m3
		
		/**dS=0.0;// LINEARISATION NOT REQUIRES AS source = Sc , slope =0*/
		}
		/*condition for no condensation*/
		else
		{
			source = 0.0;
			/**dS=0.0;*/
		}
		
		
		
	}
	
}
end_f_loop(f,tf)

/* wall zone 9*/

tf=Lookup_Thread(domain, wall_2_id);

/* START face loop over the wall thread 9*/

begin_f_loop(f, tf) /* face loop over the wall_zone_ID 8 to extract the FIRST LAYER of CELLS*/
{
	c0=F_C0(f,tf);// go through the first layer faces of cells near wall, then get those cells and its thread into c0 and t0
	t0=THREAD_T0(tf);// ADJACENT CELL THREAD ON ONE SIDE
	if (c==c0 && t==t0)
	{
	    /* Step 2 : INCLUDE THE CONDITION FOR CONDENSATION */
	    					/* ie  Twall < DPT */
	   
	  /*cell centre values*/ 
	    real rho=C_R(c0,t0);// density of cell 
	    real Yc=C_YI(c0,t0,0);// w_v at cell centre 
		real Tc=C_T(c0,t0);// cell centre temperature 
		real R_g_v=461.504;// R_g/M_v chara gas constant J/kg.K
		real pv_Pa= Yc*rho*R_g_v*Tc;// partial pressure of water vapor in the cell [pa]
		/*Clausius Clapeyron correlation for DPT*/ 
		real pv_MPa = pv_Pa * 1.0e-6;
		real Tsat_C = B / (A_cc - log(pv_MPa)) - C;// T sat in Degree celcius
		real Tsat_K=Tsat_C+273.15;// DPT in [K]
		
		
		
		
		/* face values */
		real Tf=F_T(f,tf);// face MACROS to get the temperature of particular face centre of the cell c0, tf sharing with the wall 8
		
		/*implementation of the condition for condensation */
		
		if (Tf<Tsat_K)
		{
		real A [ND_ND];
		F_AREA(A,f,tf);// vector of the face 
		real area_mag= NV_MAG(A);// FACE AREA :magnitude of vector A
		real vol=C_VOLUME(c0,t0);// volume of each cell
		real Diff=0.26e-04;// diffusion coefficient [m2/s]
		
		
		
		
	 /* delta-n calculation- applicable for structured  and orthogonal mesh */
		real cell_cent[ND_ND];
		real face_cent[ND_ND];
		
		C_CENTROID(cell_cent,c0,t0);
		F_CENTROID(face_cent,f,tf);// tf instead of t0
		
		real dist=NV_D(cell_cent, face_cent);// delta-n calculation based on structured mesh
		
		
		/*source term definition */
		
		real Y_s=Yc; //OMEGA_C * i.e. the term that changes with the iterations by updating the new cell centre value from Yc
		real k;// variable clubbing term 
		k=(-1)*(rho*Diff*area_mag)/(dist*vol);// k with -ive sign 
		
		real T_w=F_T(f,tf);// wall temperature [K] of the interface of each cell c0,t0 with face f in tf
		real t_w_C=T_w-273.15;// T_w in  degree C
		real p_s_w=610.78 * exp((17.2694 * t_w_C) / (t_w_C + 238.3)); //partial pressure [Pa] corresponding to the saturation temp t_w_C
		real Yw= p_s_w/(rho*R_g_v*T_w);// wall mass fraction for the interface calculated from the equation
		//real Sp=k*((1-Yw)/((1-Y_s)*(1-Y_s)));// (dS/dw_C)wc*= Sp
		real S_s=k*((Y_s-Yw)/(1-Y_s)); // ie S* IN THE EXPRESSION
		/*Sc=S_s-(Sp*Y_s);S*-()Species UDF requires Sp, Sc, but here we only need S_s.*/
		/* simplification of source here as here there is no requirememt of Sc+Sp phi format source = S_s after simplification*/
		source=S_s;// TOTAL MASS OF WATER CONDENSATION RATE PER UNIT VOLUME kg/s-m3
		
		/**dS=0.0;// LINEARISATION NOT REQUIRES AS source = Sc , slope =0*/
		}
		/*condition for no condensation*/
		else
		{
			source = 0.0;
			/**dS=0.0;*/
		}
		
		
		
	}
	
}
end_f_loop(f,tf)
/* END of face loop over wall thread 9*/
*dS=0.0; //slope =0, ie no linearisation implemented
return source;

}
