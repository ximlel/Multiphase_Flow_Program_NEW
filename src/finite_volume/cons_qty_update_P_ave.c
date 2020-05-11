#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "../include/var_struc.h"

int cons_qty_update_corr_ave_P(struct cell_var * cv, const struct mesh_var * mv,
							   const struct flu_var * FV, const double tau)
{
	const int dim = (int)config[0];
	const int num_cell = (int)config[3];
	const int order = (int)config[9];
	const double eps = config[4];
	int ** cp = mv->cell_pt;

	double U_u_a, U_v_a;
	int p_p, p_n;
	double length, Z_a;
//	for(int k = (int)config[13]; k < num_cell-(int)config[13]; ++k)
	for(int k = 0; k < num_cell; ++k)
		{
			U_u_a = cv->U_phi[k]*cv->U_u[k]/cv->U_rho[k];
			U_v_a = cv->U_phi[k]*cv->U_v[k]/cv->U_rho[k];			
			if (order == 1)					
				Z_a = FV->Z_a[k];
			else if (order == 2)
				Z_a = FV->Z_a[k]-0.5*tau*(cv->U_u[k]*cv->gradx_z_a[k]+cv->U_v[k]*cv->grady_z_a[k])/cv->U_rho[k];
			for(int j = 0; j < cp[k][0]; j++)
				{
					if(j == cp[k][0]-1) 
						{
							p_p=cp[k][1];
							p_n=cp[k][j+1];
						}				  
					else
						{
							p_p=cp[k][j+2];
							p_n=cp[k][j+1];
						}
					if (dim == 1)
						length = cv->n_x[k][j];
					else if (dim == 2)
						length = sqrt((mv->X[p_p] - mv->X[p_n])*(mv->X[p_p]-mv->X[p_n]) + (mv->Y[p_p] - mv->Y[p_n])*(mv->Y[p_p]-mv->Y[p_n]));				
					cv->U_rho[k] += - tau*cv->F_rho[k][j] * length / cv->vol[k];
					cv->U_e[k]   += - tau*cv->F_e[k][j]   * length / cv->vol[k];	
					cv->U_e_a[k] += - tau*(cv->F_e_a[k][j] + Z_a*cv->P_star[k][j]) * length / cv->vol[k];
					cv->U_u[k]   += - tau*cv->F_u[k][j]   * length / cv->vol[k];
					U_u_a += - tau*(cv->U_qt_add_c[k][j] + Z_a*cv->U_qt_star[k][j]) * length / cv->vol[k];
					if (dim > 1)
						{													
							cv->U_v[k] += - tau*cv->F_v[k][j] * length / cv->vol[k];
							U_v_a += - tau*(cv->V_qt_add_c[k][j] + Z_a*cv->V_qt_star[k][j]) * length / cv->vol[k];
						}
					if (dim > 2)
						cv->U_w[k] += - tau*cv->F_w[k][j] * length / cv->vol[k];
					if ((int)config[2] == 2)						
						cv->U_phi[k] += - tau*cv->F_phi[k][j] * length / cv->vol[k];
				}
			cv->U_e_a[k] += (cv->U_phi[k]*cv->U_u[k]/cv->U_rho[k]-U_u_a)*cv->U_u[k]/cv->U_rho[k];
			cv->U_e_a[k] += (cv->U_phi[k]*cv->U_v[k]/cv->U_rho[k]-U_v_a)*cv->U_v[k]/cv->U_rho[k];
		}
	return 1;
}
