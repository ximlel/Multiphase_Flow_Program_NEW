void cons_qty_init(const struct cell_var * cv, const struct flu_var * FV);
int cons2prim(struct i_f_var * ifv);
int cons_qty_update(const struct cell_var * cv, const struct mesh_var * mv,
					const struct flu_var *  FV, const double tau);
int cons_qty_update_corr_ave_P(struct cell_var * cv, const struct mesh_var * mv,
					const struct flu_var * FV, const double tau);


struct cell_var cell_mem_init(const struct mesh_var * mv, struct flu_var * FV);
void vol_comp(const struct cell_var * cv, const struct mesh_var * mv);
void cell_pt_clockwise(const struct mesh_var * mv);
void cell_rel(const struct cell_var * cv, const struct mesh_var * mv);
void cell_centroid(const struct cell_var * cv, const struct mesh_var * mv);


void slope_limiter(const struct cell_var * cv,const struct mesh_var * mv, const struct flu_var * FV);


void cons_qty_copy_cv2ifv(struct i_f_var * ifv, const struct cell_var * cv, const int c);
void cons_qty_copy_ifv2cv(const struct i_f_var * ifv, struct cell_var * cv, const int c);
void prim_var_copy_ifv2FV(const struct i_f_var * ifv, const struct flu_var * FV,const int c);
void flux_copy_ifv2cv(const struct i_f_var * ifv, const struct cell_var *cv, const int k, const int j);


int fluid_var_update(struct flu_var *FV, struct cell_var *cv);
int interface_var_init(const struct cell_var * cv, const struct mesh_var * mv,
					   struct i_f_var * ifv, struct i_f_var * ifv_R,
					   const int k, const int j, const int i);
double tau_calc(const struct cell_var * cv, const struct mesh_var * mv);


void Roe_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void HLL_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void Riemann_exact_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R);
void GRP_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau);
void GRP_2D_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau);
void GRP_IT_scheme(struct i_f_var * ifv, struct i_f_var * ifv_R, const double tau, double * phase_loc_rel, double * phase_loc);

void finite_volume_scheme(struct flu_var * FV, const struct mesh_var * mv, const char * scheme, const char * problem);
