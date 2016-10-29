#include "paper_3.h"


typedef boost::minstd_rand base_generator_type;

// continuous trait

void gibbssampler2(double *result, int * numRows, int * numCols, int * numCols2, int * fam, 
	double *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, 
	double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, double * arg_a, double * arg_b, int rvinfo)
{

	int ITER_NUM = 10000;
	int burnin = 500;

	base_generator_type generator(24);

	typedef std::vector<double> Vec;
	typedef std::vector<short> Vec_i;
	typedef std::vector<Vec> Mat;
	typedef std::vector<Vec_i> Mat_i;

	typedef boost::bernoulli_distribution<> distribution_type;
	typedef boost::normal_distribution<> distribution_type2;
	typedef boost::gamma_distribution<> distribution_type3;
	typedef boost::uniform_01<> distribution_type4;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;

	// read input files

	Vec pheno(0);
	int NUM_SUB = (*numRows);
	Mat geno;
	int NUM_RV = (*numCols);
	Mat_i geno_i;
	Mat m_error;
	Mat_i m_missing;
	Mat cov_m;
	int NUM_COV = (*numCols2);
	Vec maf_info;
	// Mat kin_m;
	int q_yes = 0;
	int cov_yes = 0;
	int maf_yes = 0;

	double thres_q = 20;
	double thres_i = 0.1;
	
	double p_a = 1.3;
	double p_b = 0.04;

	double tau_r_a = 0.001;
	double tau_r_b = 0.001;
	
		double * p = mat1;
			for(int i = 0; i < NUM_SUB; i++)
			{
				double temp = *p++;
				pheno.push_back(temp);
			}
			
		if((*arg_i) > 0)
		{
			thres_i = (*arg_i);
		}

		if((*arg_t) > 0)
		{
			thres_q = (*arg_t);
		}

		if((*arg_n) > 0)
		{
			ITER_NUM = (*arg_n);
		}

		if((*arg_bu) > 0)
		{
			burnin = (*arg_bu);
		}
		
		if((*arg_a) > 0)
		{
			p_a = (*arg_a);
		}

		if((*arg_b) > 0)
		{
			p_b = (*arg_b);
		}
		
		double * p2 = mat2;
		for(int i = 0; i < NUM_SUB; i++)
		{
			Vec row_ge;
			for(int j = 0; j < NUM_RV; j++)
			{
				double temp = *p2++;
				row_ge.push_back(temp);
			}
			geno.push_back(row_ge);
		}

		// quality information
		double * p3 = mat3;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_me(NUM_RV);
				Vec_i row_gi(NUM_RV);
				Vec_i row_mi(NUM_RV);

				for(int j = 0; j < NUM_RV; j++)
				{
					   row_me[j] = *p3++;
					   if(row_me[j]< 1)
					   {
							row_mi[j] = 1;
							if((row_me[j]<thres_i)||(row_me[j]>(1-thres_i)))
							{
								row_gi[j] = 0;
								geno[i][j] = row_me[j];
							}
							else
							{
								row_gi[j] = 1;
							}
					   }
					   else
					   {		
						   row_mi[j] = 0;
						   if(row_me[j]<thres_q)
							{
								row_gi[j] = 1;
							}
							else
							{
								row_gi[j] = 0;
							}
							row_me[j] = 1-pow(10,(row_me[j]/(-10)));
					   }
				}
				geno_i.push_back(row_gi);
				m_missing.push_back(row_mi);
				m_error.push_back(row_me);
			}
			
			q_yes = 1;
		
		
		// covariates		

		if(NUM_COV > 0)
		{
			
			double * p4 = mat4;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				for(int j = 0; j < NUM_COV; j++)
				{
					double temp = *p4++;
					row_cov.push_back(temp);
				}
				cov_m.push_back(row_cov);
			}
			cov_yes = 1;
			
		}

		// SVD kinship matrix
		Mat_i comp_ind;
		Mat comp_val;

		if((*fam)==1)
		{
			Mat kin_m;
			double * p_kin = mat_kin;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				for(int j = 0; j < NUM_SUB; j++)
				{
					double temp = *p_kin++;
					row_cov.push_back(temp);
				}
				kin_m.push_back(row_cov);
			}

			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec_i comp_ind_r;
				Vec comp_val_r;
				for(int j = 0; j < NUM_SUB; j++)
				{
					if(kin_m[j][i] != 0)
					{
						comp_ind_r.push_back(j);
						comp_val_r.push_back(kin_m[j][i]);
					}
				}
				comp_ind.push_back(comp_ind_r);
				comp_val.push_back(comp_val_r);
			}
		}

		//MAF

		if((*arg_m)>0)
		{
			
			for(int i = 0; i < NUM_RV; i++)
			{
				double temp;
				temp = arg_m[i];
				temp = (1-2*temp)/(2*temp);
				maf_info.push_back(temp);
			}
			maf_yes = 1;
		}
	

	// Gibbs sampling


	// estimate maf

	if(maf_yes == 0)
	{
		for(int i = 0; i < NUM_RV; i++)
		{
			int num_sub = 0;
			double num_min = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
				if(m_missing[j][i]==0)
				{
					num_sub++;
					num_min += geno[j][i];
				}
			}
			if(num_min > 0)
				maf_info.push_back((num_sub-num_min)/num_min);
			else
				maf_info.push_back(num_sub-1);
		}
	}

	
	// compress unobserved genotype matrix

    Mat_i ind_m;
    Mat_i geno_t;
	Mat error_m;
    Vec_i row_size_t(NUM_SUB);
	Vec_i sub_ind;
	int num_sub_err;
	int total_miss = 0;
	// Mat geno_0;

    for(int i = 0; i < NUM_SUB; i++)
	{
		Vec_i ind_r;
		Vec_i geno_r;
		Vec error_r;
		// Vec geno_e(NUM_RV,0);
		// geno_0.push_back(geno_e);

		for(int j = 0; j < NUM_RV; j++)
		{
			if(geno_i[i][j]!=0)
			{
				ind_r.push_back(j);
				geno_r.push_back(short(geno[i][j]));
				error_r.push_back(m_error[i][j]);
				total_miss++;
			}
		}
		row_size_t[i] = ind_r.size();
		
		if(row_size_t[i] > 0)
		{
			sub_ind.push_back(i);
		} 
		
		ind_m.push_back(ind_r);
		geno_t.push_back(geno_r);
		error_m.push_back(error_r);
	}
	num_sub_err = sub_ind.size();

	Vec y(NUM_SUB);

	for(int i = 0; i < NUM_SUB; i++)
	{
		y[i] = pheno[i];
	}

	Vec beta_0(ITER_NUM+1);
	Vec tau(ITER_NUM+1);
	Vec tau_ran(ITER_NUM+1);
	Vec_i delta(ITER_NUM+1);
	Vec delta_rb(ITER_NUM+1);
	// rvinfo == 1
	Mat tau_g;
	Mat gamma_0;
	Mat cov_z;
	
	Vec gamma_t(NUM_RV);
	Vec cov_t(NUM_COV);
	Vec tau_g_t(NUM_RV);
	Vec mu_t(NUM_SUB);
	Vec mu_cov_t(NUM_SUB);
	Vec ran_t(NUM_SUB);
	Vec ransum_t(NUM_SUB);

	// Initialization
	beta_0[0] = 0;
	tau[0] = 0.01;
	tau_ran[0] = 0.1;
	delta[0] = 0;
	delta_rb[0] = 0;

	Vec row_g(NUM_RV);

	for(int i = 0; i < NUM_RV; i++)
	{
		row_g[i] = 0.0;
		tau_g_t[i] = 1;
		gamma_t[i] = row_g[i];
	}

	for(int i = 0; i < NUM_COV; i++)
	{
		cov_t[i] = 0;
	}

	if(rvinfo==1)
	{
		gamma_0.push_back(gamma_t);
		tau_g.push_back(tau_g_t);
	}
	if(NUM_COV>0)
		cov_z.push_back(cov_t);

	double tau_0 = 0.01;
	// double p_a = 1.3;
	// double p_b = 0.04;

	double beta_0_t = beta_0[0];
	double tau_t = tau[0];
	double tau_ran_t = tau_ran[0];
	int delta_t = delta[0];

	for(int i = 0; i < NUM_SUB; i++)
	{
		mu_t[i] = 0;
		mu_cov_t[i] = 0;
		ran_t[i] = 0;
		ransum_t[i] = 0;
	}

	// Iteration

	for(int i = 1; i < ITER_NUM + 1; i++)
	{

		// beta_0
		double temp = 0;
		if(delta_t > 0)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp += y[j] - mu_t[j] - mu_cov_t[j] - ransum_t[j];
			}
		}
		else
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp += y[j] - mu_cov_t[j] - ransum_t[j];
			}
		}
		gen_type2 die_gen_b0(generator, distribution_type2((1/(NUM_SUB+tau_0/tau_t))*temp, sqrt(1/(NUM_SUB*tau_t+tau_0))));
		boost::generator_iterator<gen_type2> die_b0(&die_gen_b0);
		beta_0_t = *die_b0++;

		// cov

		for(int j = 0; j < NUM_COV; j++)
		{
			double new_cov = 0;
			double g_sq = 0;
			double sum_g_g = 0;

			for(int k = 0; k < NUM_SUB; k++)
			{
				if(delta_t == 0)
					sum_g_g += (y[k]-beta_0_t-mu_cov_t[k]-ransum_t[k]+cov_t[j]*cov_m[k][j])*cov_m[k][j];
				else
					sum_g_g += (y[k]-mu_t[k]-beta_0_t-mu_cov_t[k]-ransum_t[k]+cov_t[j]*cov_m[k][j])*cov_m[k][j];
				g_sq += cov_m[k][j]*cov_m[k][j];
			}

			gen_type2 die_gen_cov(generator, distribution_type2(sum_g_g/(tau_0/tau_t+g_sq),1/(sqrt(tau_0+tau_t*g_sq))));
			boost::generator_iterator<gen_type2> die_cov(&die_gen_cov);
			new_cov = *die_cov++;

			for(int k = 0; k < NUM_SUB; k++)
			{
				mu_cov_t[k] += (new_cov-cov_t[j])*cov_m[k][j];
			}
			
			cov_t[j] = new_cov;
		}
		
		// gamma
		for(int j = 0; j < NUM_RV; j++)
		{
			double new_gamma = 0;
			double g_sq = 0;
			double sum_g_g = 0;
			Vec_i flag(NUM_SUB);

			for(int k = 0; k < NUM_SUB; k++)
			{
				// double sum_g = 0;
				flag[k] = -1;
				for(int l = 0; l < row_size_t[k]; l++)
				{
					if(j == ind_m[k][l])
						flag[k] = l;
				}

				if(delta_t > 0)
				{
					if(flag[k] == -1)
					{
						sum_g_g += (y[k]-mu_t[k]+gamma_t[j]*geno[k][j]-beta_0_t-mu_cov_t[k]-ransum_t[k])*geno[k][j];
						g_sq += geno[k][j]*geno[k][j];
					}
					else
					{
						sum_g_g += (y[k]-mu_t[k] + gamma_t[j]*geno_t[k][flag[k]]-beta_0_t-mu_cov_t[k]-ransum_t[k])*geno_t[k][flag[k]];
						g_sq += geno_t[k][flag[k]]*geno_t[k][flag[k]];
					}
				}
			}
			if(delta_t==0)
			{
				gen_type2 die_gen8(generator, distribution_type2(0, sqrt(1/(tau_g_t[j]))));
				boost::generator_iterator<gen_type2> die8(&die_gen8);
				new_gamma = *die8++;
			}
			else
			{
				gen_type2 die_gen8(generator, distribution_type2(sum_g_g/(tau_g_t[j]/tau_t+g_sq), sqrt(1/(tau_g_t[j]+tau_t*g_sq))));
				boost::generator_iterator<gen_type2> die8(&die_gen8);
				new_gamma = *die8++;
			}
				

			for(int k = 0; k < NUM_SUB; k++)
			{
				if(flag[k] == -1)
				{
					mu_t[k] += (new_gamma-gamma_t[j])*geno[k][j];
				}
				else
				{
					mu_t[k] += (new_gamma-gamma_t[j])*geno_t[k][flag[k]];
				}
			}
			
			gamma_t[j] = new_gamma;
		}


		// ran_t 

		if((*fam)==1)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double new_ran = 0;
				double g_sq = 0;
				double sum_g_g = 0;

				int kin_size = comp_ind[j].size();
				for(int k = 0; k < kin_size; k++)
				{
					int real_ind = comp_ind[j][k];	
					if(delta_t == 0)
						sum_g_g += (y[real_ind]-beta_0_t-mu_cov_t[real_ind]-ransum_t[real_ind]+ran_t[j]*comp_val[j][k])*comp_val[j][k];
					else
						sum_g_g += (y[real_ind]-mu_t[real_ind]-beta_0_t-mu_cov_t[real_ind]-ransum_t[real_ind]+ran_t[j]*comp_val[j][k])*comp_val[j][k];
					g_sq += comp_val[j][k]*comp_val[j][k];
				}

				gen_type2 die_gen_ran(generator, distribution_type2(sum_g_g/(tau_ran_t/tau_t+g_sq),1/(sqrt(tau_ran_t+tau_t*g_sq))));
				boost::generator_iterator<gen_type2> die_ran(&die_gen_ran);
				new_ran = *die_ran++;

				for(int k = 0; k < kin_size; k++)
				{
					int real_ind = comp_ind[j][k];
					ransum_t[real_ind] += (new_ran-ran_t[j])*comp_val[j][k];
				}
			
				ran_t[j] = new_ran;
			}

			// tau_ran_t

			double temp_ran = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp_ran += ran_t[j]*ran_t[j];
			}
			gen_type3 die_gen_tauran(generator, distribution_type3(tau_r_a+0.5*NUM_SUB, 1/(tau_r_b+0.5*temp_ran) ));
			boost::generator_iterator<gen_type3> die_tauran(&die_gen_tauran);
			tau_ran_t = *die_tauran++;
		
		}


		// tau

		double temp4 = 0;
		if(delta_t>0)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double y_b0_b1_ax = y[j] - mu_t[j] - beta_0_t - mu_cov_t[j] - ransum_t[j];
				temp4 += y_b0_b1_ax*y_b0_b1_ax;
			}
		}
		else
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double y_b0_b1_ax = y[j] - beta_0_t - mu_cov_t[j] - ransum_t[j];
				temp4 += y_b0_b1_ax*y_b0_b1_ax;
			}
		}
	
		gen_type3 die_gen9(generator, distribution_type3(0.001+0.5*NUM_SUB, 1/(0.001+0.5*temp4)));
		boost::generator_iterator<gen_type3> die9(&die_gen9);
		tau_t = *die9++;

		beta_0[i] = beta_0_t;
		tau[i] = tau_t;
		tau_ran[i] = tau_ran_t;
		delta[i] = delta_t;
		if(rvinfo==1)
		{
			gamma_0.push_back(gamma_t);
			tau_g.push_back(tau_g_t);
		}
		if(NUM_COV>0)
			cov_z.push_back(cov_t);

		// tau_g

		for(int j = 0; j < NUM_RV; j++)
		{
			gen_type3 die_gen10(generator, distribution_type3(p_a+0.5, 1/(p_b+0.5*gamma_t[j]*gamma_t[j])));
			boost::generator_iterator<gen_type3> die10(&die_gen10);
			tau_g_t[j] = *die10++;
		}

		// delta

		double mu_sum = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			mu_sum += (2*(y[j]-beta_0_t-mu_cov_t[j]-ransum_t[j])-mu_t[j])*mu_t[j];
		}

		double check_of = tau_t*mu_sum/2;
		double r_delta = 0;
		double p_delta = 0;
		if(check_of<100)
		{
			r_delta = exp(check_of);
			p_delta = r_delta/(1+r_delta);
		}
		else
		{
			p_delta = 1;
		}

		gen_type die_gen_delta(generator, distribution_type(p_delta));
		boost::generator_iterator<gen_type> die_delta(&die_gen_delta);
		delta_t = *die_delta++;
		delta_rb[i] = p_delta;
	
		// geno_t
		
		for(int s = 0; s < num_sub_err; s++)
		{
                int j = sub_ind[s];
				
				for(int k = 0; k < row_size_t[j]; k++)
				{
					double r_1_num = 0;
					double inv_prior = 0;
					if(m_missing[j][ind_m[j][k]] == 0)
					{
						r_1_num = geno[j][ind_m[j][k]]+error_m[j][k]-2*geno[j][ind_m[j][k]]*error_m[j][k];
						
						inv_prior = maf_info[ind_m[j][k]];
					}
					else
					{
						r_1_num = 0.5;
						inv_prior = (1-error_m[j][k])/error_m[j][k];
					}

					double r_1 = (1-r_1_num)/r_1_num;

					double r_2 = 1;
					
					if(delta_t>0)
					{
						r_2 = exp(0.5*tau_t*( gamma_t[ind_m[j][k]]*(2*(y[j]-(mu_t[j] + beta_0_t + mu_cov_t[j] + ransum_t[j] - gamma_t[ind_m[j][k]]*geno_t[j][k])) - gamma_t[ind_m[j][k]]) ));
					}

					double ratio = r_1*r_2/inv_prior;
				
					
					gen_type die_gen11(generator, distribution_type(ratio/(1+ratio)));
					boost::generator_iterator<gen_type> die11(&die_gen11);
					int new_geno = *die11++;
					if(new_geno != geno_t[j][k])
					{
						mu_t[j] += (new_geno - geno_t[j][k])*gamma_t[ind_m[j][k]];
						geno_t[j][k] = new_geno;
					}
					// geno_0[j][k] += geno_t[j][k];
				}
        }

	}
	
	double beta_0_re = 0;
	double tau_re = 0;
	double tau_ran_re = 0;
	double delta_re = 0;
	double delta_rb_re = 0;
	Vec cov_re(NUM_COV);
	Vec gamma_re(NUM_RV);
	Vec gamma_sd(NUM_RV);
	Vec tau_g_re(NUM_RV);

	for(int i = 0; i < NUM_RV; i++)
	{
		gamma_re[i] = 0;
		gamma_sd[i] = 0;
		tau_g_re[i] = 0;
	}
	for(int i = 0; i < NUM_COV; i++)
	{
		cov_re[i] = 0;
	}

	// int l_t_z = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_re += beta_0[i];
		tau_re += tau[i];
		delta_re += delta[i];
		delta_rb_re += delta_rb[i];
		
		for(int j = 0; j < NUM_COV; j++)
		{
			cov_re[j] += cov_z[i][j];
		}
	}

	if(rvinfo==1)
	{
		for(int i = burnin; i < ITER_NUM; i++)
		{
			for(int j = 0; j < NUM_RV; j++)
			{
				gamma_re[j] += gamma_0[i][j];
				tau_g_re[j] += tau_g[i][j];
			}
		}
	}

	if((*fam)==1)
	{
		for(int i = burnin; i < ITER_NUM; i++)
		{
			tau_ran_re += tau_ran[i];
		}
	}

	double iter_count = ITER_NUM - burnin;
	beta_0_re /= iter_count;
	tau_re /= iter_count;
	tau_ran_re /= iter_count;
	if(rvinfo==1)
	{
		for(int j = 0; j < NUM_RV; j++)
		{
			gamma_re[j] /= iter_count;
			tau_g_re[j] /= iter_count;
		}

		double equ = 0.01;
		for(int i = ITER_NUM/10; i < ITER_NUM; i++)
		{
			for(int j = 0; j < NUM_RV; j++)
			{
				gamma_sd[j] += (gamma_0[i][j] - gamma_re[j])*(gamma_0[i][j] - gamma_re[j]);
			}
		}

		for(int j = 0; j < NUM_RV; j++)
		{
			gamma_sd[j] /= delta_re;
			gamma_sd[j] = sqrt(gamma_sd[j]);
		}
	}

	for(int j = 0; j < NUM_COV; j++)
	{
		cov_re[j] /= iter_count;
	}

	
	//beta_1_po /= (1-0.01128341555)/0.01128341555;
	delta_re /= iter_count;
	delta_rb_re /= iter_count;

	result[0] = delta_re;
	result[1] = delta_rb_re;
	result[2] = beta_0_re;
	result[3] = 1/tau_re;
	result[4] = total_miss;
	if((*fam)==1)
	{
		result[5] = 1/tau_ran_re;
	}
	else
	{
		result[5] = 0;
	}

	for(int j = 0; j < NUM_COV; j++)
	{
		result[6+j] = cov_re[j];
	}
	
	if(rvinfo==1)
	{
		for(int j = 0; j < NUM_RV; j++)
		{
			result[6+NUM_COV+j] = gamma_re[j];
		}
		for(int j = 0; j < NUM_RV; j++)
		{
			result[6+NUM_COV+NUM_RV+j] = gamma_sd[j];
		}
	}

}

// binary trait

void gibbssampler2_bin(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, double * arg_a, double * arg_b)
{

	int ITER_NUM = 10000;
	int burnin = 500;

	base_generator_type generator(24);

	typedef std::vector<double> Vec;
	typedef std::vector<short> Vec_i;
	typedef std::vector<Vec> Mat;
	typedef std::vector<Vec_i> Mat_i;

	typedef boost::bernoulli_distribution<> distribution_type;
	typedef boost::normal_distribution<> distribution_type2;
	typedef boost::gamma_distribution<> distribution_type3;
	typedef boost::uniform_01<> distribution_type4;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;

	// read input files

	Vec_i pheno(0);
	int NUM_SUB = (*numRows);
	Mat geno;
	int NUM_RV = (*numCols);
	Mat_i geno_i;
	Mat m_error;
	Mat_i m_missing;
	Mat cov_m;
	int NUM_COV = (*numCols2);
	Vec maf_info;
	// Mat kin_m;
	int q_yes = 0;
	int cov_yes = 0;
	int maf_yes = 0;

	double thres_q = 20;
	double thres_i = 0.1;
	
	double p_a = 1.3;
	double p_b = 0.04;
	
		int * p = mat1;
			for(int i = 0; i < NUM_SUB; i++)
			{
				int temp = *p++;
				pheno.push_back(temp);
			}
			
		if((*arg_i) > 0)
		{
			thres_i = (*arg_i);
		}

		if((*arg_t) > 0)
		{
			thres_q = (*arg_t);
		}

		if((*arg_n) > 0)
		{
			ITER_NUM = (*arg_n);
		}

		if((*arg_bu) > 0)
		{
			burnin = (*arg_bu);
		}
		
		if((*arg_a) > 0)
		{
			p_a = (*arg_a);
		}

		if((*arg_b) > 0)
		{
			p_b = (*arg_b);
		}
		
		double * p2 = mat2;
		for(int i = 0; i < NUM_SUB; i++)
		{
			Vec row_ge;
			for(int j = 0; j < NUM_RV; j++)
			{
				double temp = *p2++;
				row_ge.push_back(temp);
			}
			geno.push_back(row_ge);
		}

		// quality information
		double * p3 = mat3;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_me(NUM_RV);
				Vec_i row_gi(NUM_RV);
				Vec_i row_mi(NUM_RV);

				for(int j = 0; j < NUM_RV; j++)
				{
					   row_me[j] = *p3++;
					   if(row_me[j]< 1)
					   {
							row_mi[j] = 1;
							if((row_me[j]<thres_i)||(row_me[j]>(1-thres_i)))
							{
								row_gi[j] = 0;
								geno[i][j] = row_me[j];
							}
							else
							{
								row_gi[j] = 1;
							}
					   }
					   else
					   {		
						   row_mi[j] = 0;
						   if(row_me[j]<thres_q)
							{
								row_gi[j] = 1;
							}
							else
							{
								row_gi[j] = 0;
							}
							row_me[j] = 1-pow(10,(row_me[j]/(-10)));
					   }
				}
				geno_i.push_back(row_gi);
				m_missing.push_back(row_mi);
				m_error.push_back(row_me);
			}
			
			q_yes = 1;
		
		
		// covariates		

		if(NUM_COV > 0)
		{
			
			double * p4 = mat4;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				for(int j = 0; j < NUM_COV; j++)
				{
					double temp = *p4++;
					row_cov.push_back(temp);
				}
				cov_m.push_back(row_cov);
			}
			cov_yes = 1;
			
		}

		// SVD kinship matrix
		Mat_i comp_ind;
		Mat comp_val;

		if((*fam)==1)
		{
			Mat kin_m;
			double * p_kin = mat_kin;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				
				for(int j = 0; j < NUM_SUB; j++)
				{
					double temp = *p_kin++;
					row_cov.push_back(temp);				
				}
				kin_m.push_back(row_cov);
			}

			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec_i comp_ind_r;
				Vec comp_val_r;
				for(int j = 0; j < NUM_SUB; j++)
				{
					if(kin_m[j][i] != 0)
					{
						comp_ind_r.push_back(j);
						comp_val_r.push_back(kin_m[j][i]);
					}
				}
				comp_ind.push_back(comp_ind_r);
				comp_val.push_back(comp_val_r);
			}
		}

		//MAF

		if((*arg_m)>0)
		{
			
			for(int i = 0; i < NUM_RV; i++)
			{
				double temp;
				temp = arg_m[i];
				temp = (1-2*temp)/(2*temp);
				maf_info.push_back(temp);
			}
			maf_yes = 1;
		}
	

	// Gibbs sampling


	// estimate maf

	if(maf_yes == 0)
	{
		for(int i = 0; i < NUM_RV; i++)
		{
			int num_sub = 0;
			double num_min = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
				if(m_missing[j][i]==0)
				{
					num_sub++;
					num_min += geno[j][i];
				}
			}
			if(num_min > 0)
				maf_info.push_back((num_sub-num_min)/num_min);
			else
				maf_info.push_back(num_sub-1);
		}
	}

	
	// compress unobserved genotype matrix

    Mat_i ind_m;
    Mat_i geno_t;
	Mat error_m;
    Vec_i row_size_t(NUM_SUB);
	Vec_i sub_ind;
	int num_sub_err;
	int total_miss = 0;
	// Mat geno_0;

    for(int i = 0; i < NUM_SUB; i++)
	{
		Vec_i ind_r;
		Vec_i geno_r;
		Vec error_r;
		// Vec geno_e(NUM_RV,0);
		// geno_0.push_back(geno_e);

		for(int j = 0; j < NUM_RV; j++)
		{
			if(geno_i[i][j]!=0)
			{
				ind_r.push_back(j);
				geno_r.push_back(short(geno[i][j]));
				error_r.push_back(m_error[i][j]);
				total_miss++;
			}
		}
		row_size_t[i] = ind_r.size();
		
		if(row_size_t[i] > 0)
		{
			sub_ind.push_back(i);
		} 
		
		ind_m.push_back(ind_r);
		geno_t.push_back(geno_r);
		error_m.push_back(error_r);
	}
	num_sub_err = sub_ind.size();

	Vec_i y(NUM_SUB);

	for(int i = 0; i < NUM_SUB; i++)
	{
		y[i] = pheno[i];
	}

	Vec beta_0(ITER_NUM+1);
	// Vec tau(ITER_NUM+1);
	Vec tau_ran(ITER_NUM+1);
	Vec_i delta(ITER_NUM+1);
	Vec delta_rb(ITER_NUM+1);
	// Mat tau_g;
	// Mat gamma_0;
	Mat cov_z;
	
	Vec gamma_t(NUM_RV);
	Vec cov_t(NUM_COV);
	Vec tau_g_t(NUM_RV);
	Vec mu_t(NUM_SUB);
	Vec mu_cov_t(NUM_SUB);
	Vec ran_t(NUM_SUB);
	Vec ransum_t(NUM_SUB);

	// augment variable
	Vec w_t(NUM_SUB,0);

	// Initialization
	beta_0[0] = 0;
	// tau[0] = 0.01;
	tau_ran[0] = 1;
	delta[0] = 0;
	delta_rb[0] = 0;

	Vec row_g(NUM_RV);

	for(int i = 0; i < NUM_RV; i++)
	{
		row_g[i] = 0.0;
		tau_g_t[i] = 1;
		gamma_t[i] = row_g[i];
	}

	for(int i = 0; i < NUM_COV; i++)
	{
		cov_t[i] = 0;
	}

	// gamma_0.push_back(gamma_t);
	// tau_g.push_back(tau_g_t);
	if(NUM_COV>0)
		cov_z.push_back(cov_t);

	double tau_0 = 0.01;
	double tau_r_a = 1;
	double tau_r_b = 0.622;
	// double p_a = 1.3;
	// double p_b = 0.04;

	double beta_0_t = beta_0[0];
	// double tau_t = tau[0];
	double tau_ran_t = tau_ran[0];
	int delta_t = delta[0];

	for(int i = 0; i < NUM_SUB; i++)
	{
		mu_t[i] = 0;
		mu_cov_t[i] = 0;
		ran_t[i] = 0;
		ransum_t[i] = 0;
	}

	// Iteration

	for(int i = 1; i < ITER_NUM + 1; i++)
	{

		// beta_0
		double temp = 0;
		if(delta_t > 0)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp += w_t[j] - mu_t[j] - mu_cov_t[j] - ransum_t[j];
			}
		}
		else
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp += w_t[j] - mu_cov_t[j] - ransum_t[j];
			}
		}
		gen_type2 die_gen_b0(generator, distribution_type2((1/(NUM_SUB+tau_0))*temp, sqrt(1/(NUM_SUB+tau_0))));
		boost::generator_iterator<gen_type2> die_b0(&die_gen_b0);
		beta_0_t = *die_b0++;

		// cov

		for(int j = 0; j < NUM_COV; j++)
		{
			double new_cov = 0;
			double g_sq = 0;
			double sum_g_g = 0;

			for(int k = 0; k < NUM_SUB; k++)
			{
				if(delta_t == 0)
					sum_g_g += (w_t[k]-beta_0_t-mu_cov_t[k]-ransum_t[k]+cov_t[j]*cov_m[k][j])*cov_m[k][j];
				else
					sum_g_g += (w_t[k]-mu_t[k]-beta_0_t-mu_cov_t[k]-ransum_t[k]+cov_t[j]*cov_m[k][j])*cov_m[k][j];
				g_sq += cov_m[k][j]*cov_m[k][j];
			}

			gen_type2 die_gen_cov(generator, distribution_type2(sum_g_g/(tau_0+g_sq),1/(sqrt(tau_0+g_sq))));
			boost::generator_iterator<gen_type2> die_cov(&die_gen_cov);
			new_cov = *die_cov++;

			for(int k = 0; k < NUM_SUB; k++)
			{
				mu_cov_t[k] += (new_cov-cov_t[j])*cov_m[k][j];
			}
			
			cov_t[j] = new_cov;
		}
		
		// gamma
		for(int j = 0; j < NUM_RV; j++)
		{
			double new_gamma = 0;
			double g_sq = 0;
			double sum_g_g = 0;
			Vec_i flag(NUM_SUB);

			for(int k = 0; k < NUM_SUB; k++)
			{
				// double sum_g = 0;
				flag[k] = -1;
				for(int l = 0; l < row_size_t[k]; l++)
				{
					if(j == ind_m[k][l])
						flag[k] = l;
				}

				if(delta_t > 0)
				{
					if(flag[k] == -1)
					{
						sum_g_g += (w_t[k]-mu_t[k]+gamma_t[j]*geno[k][j]-beta_0_t-mu_cov_t[k]-ransum_t[k])*geno[k][j];
						g_sq += geno[k][j]*geno[k][j];
					}
					else
					{
						sum_g_g += (w_t[k]-mu_t[k] + gamma_t[j]*geno_t[k][flag[k]]-beta_0_t-mu_cov_t[k]-ransum_t[k])*geno_t[k][flag[k]];
						g_sq += geno_t[k][flag[k]]*geno_t[k][flag[k]];
					}
				}
			}
			if(delta_t==0)
			{
				gen_type2 die_gen8(generator, distribution_type2(0, sqrt(1/(tau_g_t[j]))));
				boost::generator_iterator<gen_type2> die8(&die_gen8);
				new_gamma = *die8++;
			}
			else
			{
				gen_type2 die_gen8(generator, distribution_type2(sum_g_g/(tau_g_t[j]+g_sq), sqrt(1/(tau_g_t[j]+g_sq))));
				boost::generator_iterator<gen_type2> die8(&die_gen8);
				new_gamma = *die8++;
			}
				

			for(int k = 0; k < NUM_SUB; k++)
			{
				if(flag[k] == -1)
				{
					mu_t[k] += (new_gamma-gamma_t[j])*geno[k][j];
				}
				else
				{
					mu_t[k] += (new_gamma-gamma_t[j])*geno_t[k][flag[k]];
				}
			}
			
			gamma_t[j] = new_gamma;
		}


		// ran_t 

		if((*fam)==1)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double new_ran = 0;
				double g_sq = 0;
				double sum_g_g = 0;

				int kin_size = comp_ind[j].size();
				for(int k = 0; k < kin_size; k++)
				{
						int real_ind = comp_ind[j][k];
						if(delta_t == 0)
							sum_g_g += (w_t[real_ind]-beta_0_t-mu_cov_t[real_ind]-ransum_t[real_ind]+ran_t[j]*comp_val[j][k])*comp_val[j][k];
						else
							sum_g_g += (w_t[real_ind]-mu_t[real_ind]-beta_0_t-mu_cov_t[real_ind]-ransum_t[real_ind]+ran_t[j]*comp_val[j][k])*comp_val[j][k];
						g_sq += comp_val[j][k]*comp_val[j][k];
				}

				gen_type2 die_gen_ran(generator, distribution_type2(sum_g_g/(tau_ran_t+g_sq),1/(sqrt(tau_ran_t+g_sq))));
				boost::generator_iterator<gen_type2> die_ran(&die_gen_ran);
				new_ran = *die_ran++;

				for(int k = 0; k < kin_size; k++)
				{
					int real_ind = comp_ind[j][k];	
					ransum_t[real_ind] += (new_ran-ran_t[j])*comp_val[j][k];
				}
			
				ran_t[j] = new_ran;
			}

			// tau_ran_t

			double temp_ran = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp_ran += ran_t[j]*ran_t[j];
			}
			gen_type3 die_gen_tauran(generator, distribution_type3(tau_r_a+0.5*NUM_SUB, 1/(tau_r_b+0.5*temp_ran) ));
			boost::generator_iterator<gen_type3> die_tauran(&die_gen_tauran);
			tau_ran_t = *die_tauran++;
		
		}


		// tau
		/*
		double temp4 = 0;
		if(delta_t>0)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double y_b0_b1_ax = y[j] - mu_t[j] - beta_0_t - mu_cov_t[j] - ransum_t[j];
				temp4 += y_b0_b1_ax*y_b0_b1_ax;
			}
		}
		else
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double y_b0_b1_ax = y[j] - beta_0_t - mu_cov_t[j] - ransum_t[j];
				temp4 += y_b0_b1_ax*y_b0_b1_ax;
			}
		}
	
		gen_type3 die_gen9(generator, distribution_type3(0.001+0.5*NUM_SUB, 1/(0.001+0.5*temp4)));
		boost::generator_iterator<gen_type3> die9(&die_gen9);
		tau_t = *die9++;
		*/

		// w
		for(int j = 0; j < NUM_SUB; j++)
		{
			double eta = delta_t*mu_t[j]+beta_0_t+mu_cov_t[j]+ransum_t[j];
			// double eta = beta_1_t*alpha_sum[j]*x[j]+beta_0_t+var_c_t[j]+ var_a_m_t[j] + var_a_d1_t[j]+ var_a_d2_t[j]+mu_cov_t[j];
			gen_type2 die_gen_w(generator, distribution_type2(eta, 1));
			boost::generator_iterator<gen_type2> die_w(&die_gen_w);
			if(y[j]==1)
			{
				if(eta>0)
				{
					w_t[j] = *die_w++;
					while(w_t[j] < 0)
					{
						w_t[j] = *die_w++;
					}
				}
				else
				{
					gen_type4 die_gen_u1(generator, distribution_type4());
					boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
					gen_type2 die_gen_z(generator, distribution_type2(0, 1));
					boost::generator_iterator<gen_type2> die_z(&die_gen_z);
					// double z1 = 1;
					// double z2 = 1;
					double u1 = 2;
					double z = 0;
					// while(u1>exp((-0.5)*(z2-z1)))
					while(u1>exp(z*eta))
					{	
						z = *die_z++;
						while(z < 0)
						{
							z = *die_z++;
						}
						// z1 = z*z;
						// z2 = (z-eta)*(z-eta)-eta*eta;
					
						u1 = *die_u1++;
					}
					w_t[j] = z;
				}
			}
			else
			{
				if(eta<0)
				{
					w_t[j] = *die_w++;
					while(w_t[j] >= 0)
					{
						w_t[j] = *die_w++;
					}
				}
				else
				{
					gen_type4 die_gen_u1(generator, distribution_type4());
					boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
					gen_type2 die_gen_z(generator, distribution_type2(0, 1));
					boost::generator_iterator<gen_type2> die_z(&die_gen_z);
					
					double u1 = 2;
					double z = 0;
					while(u1>exp((-1)*z*eta))
					{	
						z = *die_z++;
						while(z >= 0)
						{
							z = *die_z++;
						}
						
						u1 = *die_u1++;
					}
					w_t[j] = z;
				}
			}
		}
		

		beta_0[i] = beta_0_t;
		// tau[i] = tau_t;
		tau_ran[i] = tau_ran_t;
		delta[i] = delta_t;
		// gamma_0.push_back(gamma_t);
		// tau_g.push_back(tau_g_t);
		if(NUM_COV>0)
			cov_z.push_back(cov_t);

		// tau_g

		for(int j = 0; j < NUM_RV; j++)
		{
			gen_type3 die_gen10(generator, distribution_type3(p_a+0.5, 1/(p_b+0.5*gamma_t[j]*gamma_t[j])));
			boost::generator_iterator<gen_type3> die10(&die_gen10);
			tau_g_t[j] = *die10++;
		}

		// delta

		double mu_sum = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			mu_sum += (2*(w_t[j]-beta_0_t-mu_cov_t[j]-ransum_t[j])-mu_t[j])*mu_t[j];
		}

		double check_of = mu_sum/2;
		double r_delta = 0;
		double p_delta = 0;
		if(check_of<100)
		{
			r_delta = exp(check_of);
			p_delta = r_delta/(1+r_delta);
		}
		else
		{
			p_delta = 1;
		}

		gen_type die_gen_delta(generator, distribution_type(p_delta));
		boost::generator_iterator<gen_type> die_delta(&die_gen_delta);
		delta_t = *die_delta++;
		delta_rb[i] = p_delta;
	
		// geno_t
		
		for(int s = 0; s < num_sub_err; s++)
		{
                int j = sub_ind[s];
				
				for(int k = 0; k < row_size_t[j]; k++)
				{
					double r_1_num = 0;
					double inv_prior = 0;
					if(m_missing[j][ind_m[j][k]] == 0)
					{
						r_1_num = geno[j][ind_m[j][k]]+error_m[j][k]-2*geno[j][ind_m[j][k]]*error_m[j][k];
						
						inv_prior = maf_info[ind_m[j][k]];
					}
					else
					{
						r_1_num = 0.5;
						inv_prior = (1-error_m[j][k])/error_m[j][k];
					}

					double r_1 = (1-r_1_num)/r_1_num;

					double r_2 = 1;
					
					if(delta_t>0)
					{
						r_2 = exp(0.5*( gamma_t[ind_m[j][k]]*(2*(w_t[j]-(mu_t[j] + beta_0_t + mu_cov_t[j] + ransum_t[j] - gamma_t[ind_m[j][k]]*geno_t[j][k])) - gamma_t[ind_m[j][k]]) ));
					}

					double ratio = r_1*r_2/inv_prior;
				
					
					gen_type die_gen11(generator, distribution_type(ratio/(1+ratio)));
					boost::generator_iterator<gen_type> die11(&die_gen11);
					int new_geno = *die11++;
					if(new_geno != geno_t[j][k])
					{
						mu_t[j] += (new_geno - geno_t[j][k])*gamma_t[ind_m[j][k]];
						geno_t[j][k] = new_geno;
					}
					// geno_0[j][k] += geno_t[j][k];
				}
        }

	}
	
	double beta_0_re = 0;
	//double tau_re = 0;
	double tau_ran_re = 0;
	double delta_re = 0;
	double delta_rb_re = 0;
	Vec cov_re(NUM_COV);
	Vec gamma_re(NUM_RV);
	Vec gamma_sd(NUM_RV);
	Vec tau_g_re(NUM_RV);

	for(int i = 0; i < NUM_RV; i++)
	{
		gamma_re[i] = 0;
		gamma_sd[i] = 0;
		tau_g_re[i] = 0;
	}
	for(int i = 0; i < NUM_COV; i++)
	{
		cov_re[i] = 0;
	}

	// int l_t_z = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_re += beta_0[i];
		// tau_re += tau[i];
		delta_re += delta[i];
		delta_rb_re += delta_rb[i];
		/*
		for(int j = 0; j < NUM_RV; j++)
		{
			gamma_re[j] += gamma_0[i][j];
			tau_g_re[j] += tau_g[i][j];
		}
		*/
		for(int j = 0; j < NUM_COV; j++)
		{
			cov_re[j] += cov_z[i][j];
		}
	}

	if((*fam)==1)
	{
		for(int i = burnin; i < ITER_NUM; i++)
		{
			tau_ran_re += tau_ran[i];
		}
	}

	double iter_count = ITER_NUM - burnin;
	beta_0_re /= iter_count;
	// tau_re /= iter_count;
	tau_ran_re /= iter_count;
	/*
	for(int j = 0; j < NUM_RV; j++)
	{
		gamma_re[j] /= iter_count;
		tau_g_re[j] /= iter_count;
	}
	*/
	for(int j = 0; j < NUM_COV; j++)
	{
		cov_re[j] /= iter_count;
	}

	/*double equ = 0.01;
	for(int i = ITER_NUM/10; i < ITER_NUM; i++)
	{
		for(int j = 0; j < NUM_RV; j++)
		{
			gamma_sd[j] += (gamma_0[i][j] - gamma_re[j])*(gamma_0[i][j] - gamma_re[j]);
		}
	}

	for(int j = 0; j < NUM_RV; j++)
	{
		gamma_sd[j] /= delta_re;
		gamma_sd[j] = sqrt(gamma_sd[j]);
	}*/
	
	delta_re /= iter_count;
	delta_rb_re /= iter_count;

	result[0] = delta_re;
	result[1] = delta_rb_re;
	result[2] = beta_0_re;
	// result[3] = 1/tau_re;
	result[3] = total_miss;
	if((*fam)==1)
	{
		result[4] = 1/tau_ran_re;
	}
	else
	{
		result[4] = 0;
	}

	for(int j = 0; j < NUM_COV; j++)
	{
		result[5+j] = cov_re[j];
	}

}


// ordinal trait

void gibbssampler2_ord(double *result, int * numRows, int * numCols, int * numCols2, int * fam, int *mat1, double *mat2, double * mat3, double * mat4, double * mat_kin, double * arg_m, double * arg_i, int * arg_n, int * arg_bu, double * arg_t, double * arg_a, double * arg_b)
{

	int ITER_NUM = 10000;
	int burnin = 500;

	base_generator_type generator(24);

	typedef std::vector<double> Vec;
	typedef std::vector<short> Vec_i;
	typedef std::vector<Vec> Mat;
	typedef std::vector<Vec_i> Mat_i;

	typedef boost::bernoulli_distribution<> distribution_type;
	typedef boost::normal_distribution<> distribution_type2;
	typedef boost::gamma_distribution<> distribution_type3;
	typedef boost::uniform_01<> distribution_type4;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;

	// read input files

	Vec_i pheno(0);
	int NUM_SUB = (*numRows);
	Mat geno;
	int NUM_RV = (*numCols);
	Mat_i geno_i;
	Mat m_error;
	Mat_i m_missing;
	Mat cov_m;
	int NUM_COV = (*numCols2);
	Vec maf_info;
	// Mat kin_m;
	int q_yes = 0;
	int cov_yes = 0;
	int maf_yes = 0;

	double thres_q = 20;
	double thres_i = 0.1;
	
	double p_a = 1.3;
	double p_b = 0.04;

	double tau_r_a = 1;
	double tau_r_b = 0.622;
	
		int * p = mat1;
			for(int i = 0; i < NUM_SUB; i++)
			{
				int temp = *p++;
				pheno.push_back(temp);
			}
			
		if((*arg_i) > 0)
		{
			thres_i = (*arg_i);
		}

		if((*arg_t) > 0)
		{
			thres_q = (*arg_t);
		}

		if((*arg_n) > 0)
		{
			ITER_NUM = (*arg_n);
		}

		if((*arg_bu) > 0)
		{
			burnin = (*arg_bu);
		}
		
		if((*arg_a) > 0)
		{
			p_a = (*arg_a);
		}

		if((*arg_b) > 0)
		{
			p_b = (*arg_b);
		}
		
		double * p2 = mat2;
		for(int i = 0; i < NUM_SUB; i++)
		{
			Vec row_ge;
			for(int j = 0; j < NUM_RV; j++)
			{
				double temp = *p2++;
				row_ge.push_back(temp);
			}
			geno.push_back(row_ge);
		}

		// quality information
		double * p3 = mat3;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_me(NUM_RV);
				Vec_i row_gi(NUM_RV);
				Vec_i row_mi(NUM_RV);

				for(int j = 0; j < NUM_RV; j++)
				{
					   row_me[j] = *p3++;
					   if(row_me[j]< 1)
					   {
							row_mi[j] = 1;
							if((row_me[j]<thres_i)||(row_me[j]>(1-thres_i)))
							{
								row_gi[j] = 0;
								geno[i][j] = row_me[j];
							}
							else
							{
								row_gi[j] = 1;
							}
					   }
					   else
					   {		
						   row_mi[j] = 0;
						   if(row_me[j]<thres_q)
							{
								row_gi[j] = 1;
							}
							else
							{
								row_gi[j] = 0;
							}
							row_me[j] = 1-pow(10,(row_me[j]/(-10)));
					   }
				}
				geno_i.push_back(row_gi);
				m_missing.push_back(row_mi);
				m_error.push_back(row_me);
			}
			
			q_yes = 1;
		
		
		// covariates		

		if(NUM_COV > 0)
		{
			
			double * p4 = mat4;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				for(int j = 0; j < NUM_COV; j++)
				{
					double temp = *p4++;
					row_cov.push_back(temp);
				}
				cov_m.push_back(row_cov);
			}
			cov_yes = 1;
			
		}

		// SVD kinship matrix
		Mat_i comp_ind;
		Mat comp_val;

		if((*fam)==1)
		{
			Mat kin_m;
			double * p_kin = mat_kin;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				
				for(int j = 0; j < NUM_SUB; j++)
				{
					double temp = *p_kin++;
					row_cov.push_back(temp);				
				}
				kin_m.push_back(row_cov);
			}

			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec_i comp_ind_r;
				Vec comp_val_r;
				for(int j = 0; j < NUM_SUB; j++)
				{
					if(kin_m[j][i] != 0)
					{
						comp_ind_r.push_back(j);
						comp_val_r.push_back(kin_m[j][i]);
					}
				}
				comp_ind.push_back(comp_ind_r);
				comp_val.push_back(comp_val_r);
			}
		}

		//MAF

		if((*arg_m)>0)
		{
			
			for(int i = 0; i < NUM_RV; i++)
			{
				double temp;
				temp = arg_m[i];
				temp = (1-2*temp)/(2*temp);
				maf_info.push_back(temp);
			}
			maf_yes = 1;
		}
	

	// Gibbs sampling


	// estimate maf

	if(maf_yes == 0)
	{
		for(int i = 0; i < NUM_RV; i++)
		{
			int num_sub = 0;
			double num_min = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
				if(m_missing[j][i]==0)
				{
					num_sub++;
					num_min += geno[j][i];
				}
			}
			if(num_min > 0)
				maf_info.push_back((num_sub-num_min)/num_min);
			else
				maf_info.push_back(num_sub-1);
		}
	}

	
	// compress unobserved genotype matrix

    Mat_i ind_m;
    Mat_i geno_t;
	Mat error_m;
    Vec_i row_size_t(NUM_SUB);
	Vec_i sub_ind;
	int num_sub_err;
	int total_miss = 0;
	// Mat geno_0;

    for(int i = 0; i < NUM_SUB; i++)
	{
		Vec_i ind_r;
		Vec_i geno_r;
		Vec error_r;
		// Vec geno_e(NUM_RV,0);
		// geno_0.push_back(geno_e);

		for(int j = 0; j < NUM_RV; j++)
		{
			if(geno_i[i][j]!=0)
			{
				ind_r.push_back(j);
				geno_r.push_back(short(geno[i][j]));
				error_r.push_back(m_error[i][j]);
				total_miss++;
			}
		}
		row_size_t[i] = ind_r.size();
		
		if(row_size_t[i] > 0)
		{
			sub_ind.push_back(i);
		} 
		
		ind_m.push_back(ind_r);
		geno_t.push_back(geno_r);
		error_m.push_back(error_r);
	}
	num_sub_err = sub_ind.size();

	Vec_i y(NUM_SUB);

	// number of levels
	int level = 0;
	Vec_i levels(0);

	// store all levels in an ascending order
	for(int i = 0; i < NUM_SUB; i++)
	{
		// y[i] = pheno[i];
		bool exist = false;
		int i_ord = 0;
		for(; i_ord < level; i_ord++)
		{
			if(pheno[i]==levels[i_ord])
				exist = true;
		}
		if(exist==false)
		{
			int j = 0;
			for(; j < level; j++)
			{
				if(levels[j]>pheno[i])
				{
					levels.push_back(levels[level-1]);
					level++;
					for(int k = level-1; k > j; k--)
						levels[k] = levels[k-1];
					levels[j] = pheno[i];
					break;
				}
			}
			if(j==level)
			{
				levels.push_back(pheno[i]);
				level++;
			}
		}
	}

	Mat_i level_list(level);
	// index phenotypes according to the levels
	for(int i = 0; i < NUM_SUB; i++)
	{
		int j = 0;
		for(; j < level; j++)
			if(pheno[i]==levels[j])
				break;
		y[i] = j;
		level_list[j].push_back(i);
	}


	Vec beta_0(ITER_NUM+1);
	// Vec tau(ITER_NUM+1);
	Vec tau_ran(ITER_NUM+1);
	Vec_i delta(ITER_NUM+1);
	Vec delta_rb(ITER_NUM+1);
	// Mat tau_g;
	// Mat gamma_0;
	Mat cov_z;
	
	Vec gamma_t(NUM_RV);
	Vec cov_t(NUM_COV);
	Vec tau_g_t(NUM_RV);
	Vec mu_t(NUM_SUB);
	Vec mu_cov_t(NUM_SUB);
	Vec ran_t(NUM_SUB);
	Vec ransum_t(NUM_SUB);

	// augment variable
	Vec w_t(NUM_SUB,0);
	// ordinal threshold
	Vec t_t(level+1,0);
	t_t[0] = -100;
	t_t[1] = 0;
	t_t[level] = 100;
	for(int i = 2; i < level; i++)
	{
		t_t[i] = (i-1)/double(level-2);
	}

	// Initialization
	beta_0[0] = 0;
	// tau[0] = 0.01;
	tau_ran[0] = 1;
	delta[0] = 0;
	delta_rb[0] = 0;

	Vec row_g(NUM_RV);

	for(int i = 0; i < NUM_RV; i++)
	{
		row_g[i] = 0.0;
		tau_g_t[i] = 1;
		gamma_t[i] = row_g[i];
	}

	for(int i = 0; i < NUM_COV; i++)
	{
		cov_t[i] = 0;
	}

	// gamma_0.push_back(gamma_t);
	// tau_g.push_back(tau_g_t);
	if(NUM_COV>0)
		cov_z.push_back(cov_t);

	double tau_0 = 0.01;
	// double p_a = 1.3;
	// double p_b = 0.04;

	double beta_0_t = beta_0[0];
	// double tau_t = tau[0];
	double tau_ran_t = tau_ran[0];
	int delta_t = delta[0];

	for(int i = 0; i < NUM_SUB; i++)
	{
		mu_t[i] = 0;
		mu_cov_t[i] = 0;
		ran_t[i] = 0;
		ransum_t[i] = 0;
	}

	// Iteration

	for(int i = 1; i < ITER_NUM + 1; i++)
	{

		// beta_0
		double temp = 0;
		if(delta_t > 0)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp += w_t[j] - mu_t[j] - mu_cov_t[j] - ransum_t[j];
			}
		}
		else
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp += w_t[j] - mu_cov_t[j] - ransum_t[j];
			}
		}
		gen_type2 die_gen_b0(generator, distribution_type2((1/(NUM_SUB+tau_0))*temp, sqrt(1/(NUM_SUB+tau_0))));
		boost::generator_iterator<gen_type2> die_b0(&die_gen_b0);
		beta_0_t = *die_b0++;

		// cov

		for(int j = 0; j < NUM_COV; j++)
		{
			double new_cov = 0;
			double g_sq = 0;
			double sum_g_g = 0;

			for(int k = 0; k < NUM_SUB; k++)
			{
				if(delta_t == 0)
					sum_g_g += (w_t[k]-beta_0_t-mu_cov_t[k]-ransum_t[k]+cov_t[j]*cov_m[k][j])*cov_m[k][j];
				else
					sum_g_g += (w_t[k]-mu_t[k]-beta_0_t-mu_cov_t[k]-ransum_t[k]+cov_t[j]*cov_m[k][j])*cov_m[k][j];
				g_sq += cov_m[k][j]*cov_m[k][j];
			}

			gen_type2 die_gen_cov(generator, distribution_type2(sum_g_g/(tau_0+g_sq),1/(sqrt(tau_0+g_sq))));
			boost::generator_iterator<gen_type2> die_cov(&die_gen_cov);
			new_cov = *die_cov++;

			for(int k = 0; k < NUM_SUB; k++)
			{
				mu_cov_t[k] += (new_cov-cov_t[j])*cov_m[k][j];
			}
			
			cov_t[j] = new_cov;
		}
		
		// gamma
		for(int j = 0; j < NUM_RV; j++)
		{
			double new_gamma = 0;
			double g_sq = 0;
			double sum_g_g = 0;
			Vec_i flag(NUM_SUB);

			for(int k = 0; k < NUM_SUB; k++)
			{
				// double sum_g = 0;
				flag[k] = -1;
				for(int l = 0; l < row_size_t[k]; l++)
				{
					if(j == ind_m[k][l])
						flag[k] = l;
				}

				if(delta_t > 0)
				{
					if(flag[k] == -1)
					{
						sum_g_g += (w_t[k]-mu_t[k]+gamma_t[j]*geno[k][j]-beta_0_t-mu_cov_t[k]-ransum_t[k])*geno[k][j];
						g_sq += geno[k][j]*geno[k][j];
					}
					else
					{
						sum_g_g += (w_t[k]-mu_t[k] + gamma_t[j]*geno_t[k][flag[k]]-beta_0_t-mu_cov_t[k]-ransum_t[k])*geno_t[k][flag[k]];
						g_sq += geno_t[k][flag[k]]*geno_t[k][flag[k]];
					}
				}
			}
			if(delta_t==0)
			{
				gen_type2 die_gen8(generator, distribution_type2(0, sqrt(1/(tau_g_t[j]))));
				boost::generator_iterator<gen_type2> die8(&die_gen8);
				new_gamma = *die8++;
			}
			else
			{
				gen_type2 die_gen8(generator, distribution_type2(sum_g_g/(tau_g_t[j]+g_sq), sqrt(1/(tau_g_t[j]+g_sq))));
				boost::generator_iterator<gen_type2> die8(&die_gen8);
				new_gamma = *die8++;
			}
				

			for(int k = 0; k < NUM_SUB; k++)
			{
				if(flag[k] == -1)
				{
					mu_t[k] += (new_gamma-gamma_t[j])*geno[k][j];
				}
				else
				{
					mu_t[k] += (new_gamma-gamma_t[j])*geno_t[k][flag[k]];
				}
			}
			
			gamma_t[j] = new_gamma;
		}


		// ran_t 

		if((*fam)==1)
		{
			for(int j = 0; j < NUM_SUB; j++)
			{
				double new_ran = 0;
				double g_sq = 0;
				double sum_g_g = 0;

				int kin_size = comp_ind[j].size();
				for(int k = 0; k < kin_size; k++)
				{
						int real_ind = comp_ind[j][k];
						if(delta_t == 0)
							sum_g_g += (w_t[real_ind]-beta_0_t-mu_cov_t[real_ind]-ransum_t[real_ind]+ran_t[j]*comp_val[j][k])*comp_val[j][k];
						else
							sum_g_g += (w_t[real_ind]-mu_t[real_ind]-beta_0_t-mu_cov_t[real_ind]-ransum_t[real_ind]+ran_t[j]*comp_val[j][k])*comp_val[j][k];
						g_sq += comp_val[j][k]*comp_val[j][k];
				}

				gen_type2 die_gen_ran(generator, distribution_type2(sum_g_g/(tau_ran_t+g_sq),1/(sqrt(tau_ran_t+g_sq))));
				boost::generator_iterator<gen_type2> die_ran(&die_gen_ran);
				new_ran = *die_ran++;

				for(int k = 0; k < kin_size; k++)
				{
					int real_ind = comp_ind[j][k];	
					ransum_t[real_ind] += (new_ran-ran_t[j])*comp_val[j][k];
				}
			
				ran_t[j] = new_ran;
			}

			// tau_ran_t

			double temp_ran = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
				temp_ran += ran_t[j]*ran_t[j];
			}
			gen_type3 die_gen_tauran(generator, distribution_type3(tau_r_a+0.5*NUM_SUB, 1/(tau_r_b+0.5*temp_ran) ));
			boost::generator_iterator<gen_type3> die_tauran(&die_gen_tauran);
			tau_ran_t = *die_tauran++;
		
		}


		// w_t
		for(int j = 0; j < NUM_SUB; j++)
		{
			double eta = delta_t*mu_t[j]+beta_0_t+mu_cov_t[j]+ransum_t[j];
			double l_b = t_t[y[j]];
			double u_b = t_t[y[j]+1];
			double interv = u_b - l_b;
			double max_w = eta;

			gen_type4 die_gen_u1(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);

			double prop = 0;

			if(eta<l_b)
			{
				max_w = l_b;
				if(interv<1)
				{
					prop = *die_u1++;
					prop = l_b + prop*(u_b-l_b);
					double rej = *die_u1++;
					while(rej>exp((-0.5)*(prop-max_w)*(prop+max_w-2*eta)))
					{
						prop = *die_u1++;
						prop = l_b + prop*(u_b-l_b);
						rej = *die_u1++;
					}
				}
				else
				{
					gen_type2 die_gen_w(generator, distribution_type2(max_w,1));
					boost::generator_iterator<gen_type2> die_w(&die_gen_w);
					prop = *die_w++;
					while((prop<max_w)||(prop>u_b))
					{
						prop = *die_w++;
					}
					double rej = *die_u1++;
					while(rej>exp((max_w-eta)*(max_w-prop)))
					{
						prop = *die_w++;
						while((prop<max_w)||(prop>u_b))
						{
							prop = *die_w++;
						}
						rej = *die_u1++;
					}
				}
			}
			else
			{
				if(eta>u_b)
				{
					max_w = u_b;
					if(interv<1)
					{
						prop = *die_u1++;
						prop = l_b + prop*(u_b-l_b);
						double rej = *die_u1++;
						while(rej>exp((-0.5)*(prop-max_w)*(prop+max_w-2*eta)))
						{
							prop = *die_u1++;
							prop = l_b + prop*(u_b-l_b);
							rej = *die_u1++;
						}
					}
					else
					{
						gen_type2 die_gen_w(generator, distribution_type2(max_w,1));
						boost::generator_iterator<gen_type2> die_w(&die_gen_w);
						prop = *die_w++;
						while((prop>max_w)||(prop<l_b))
						{
							prop = *die_w++;
						}
						double rej = *die_u1++;
						while(rej>exp((max_w-eta)*(max_w-prop)))
						{
							prop = *die_w++;
							while((prop>max_w)||(prop<l_b))
							{
								prop = *die_w++;
							}
							rej = *die_u1++;
						}
					}
				}
				else
				{
					if(interv<1)
					{
						prop = *die_u1++;
						prop = l_b + prop*(u_b-l_b);
						double rej = *die_u1++;
						while(rej>exp((-0.5)*(prop-max_w)*(prop+max_w-2*eta)))
						{
							prop = *die_u1++;
							prop = l_b + prop*(u_b-l_b);
							rej = *die_u1++;
						}
					}
					else
					{
						gen_type2 die_gen_w(generator, distribution_type2(max_w,1));
						boost::generator_iterator<gen_type2> die_w(&die_gen_w);
						prop = *die_w++;
						while((prop<l_b)||(prop>u_b))
						{
							prop = *die_w++;
						}
					}
				}
			}
		
			w_t[j] = prop;
			
		}

		// t_t

		Vec w_min(level);
		Vec w_max(level);
		for(int j = 0; j < level; j++)
		{
			double temp_min = w_t[level_list[j][0]];
			double temp_max = temp_min;
			for(int k = 1; k < level_list[j].size(); k++)
			{
				double temp_w = w_t[level_list[j][k]];
				if(temp_w<temp_min)
					temp_min = temp_w;
				else
				{
					if(temp_w>temp_max)
						temp_max = temp_w;
				}
			}
			w_min[j] = temp_min;
			w_max[j] = temp_max;
		}

		gen_type4 die_gen_u_t(generator, distribution_type4());
		boost::generator_iterator<gen_type4> die_u_t(&die_gen_u_t);
		for(int j = 2; j < level; j++)
		{
			double new_t = *die_u_t++;
			t_t[j] = w_max[j-1]+(w_min[j]-w_max[j-1])*new_t;
		}
		

		beta_0[i] = beta_0_t;
		// tau[i] = tau_t;
		tau_ran[i] = tau_ran_t;
		delta[i] = delta_t;
		// gamma_0.push_back(gamma_t);
		// tau_g.push_back(tau_g_t);
		if(NUM_COV>0)
			cov_z.push_back(cov_t);

		// tau_g

		for(int j = 0; j < NUM_RV; j++)
		{
			gen_type3 die_gen10(generator, distribution_type3(p_a+0.5, 1/(p_b+0.5*gamma_t[j]*gamma_t[j])));
			boost::generator_iterator<gen_type3> die10(&die_gen10);
			tau_g_t[j] = *die10++;
		}

		// delta

		double mu_sum = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			mu_sum += (2*(w_t[j]-beta_0_t-mu_cov_t[j]-ransum_t[j])-mu_t[j])*mu_t[j];
		}

		double check_of = mu_sum/2;
		double r_delta = 0;
		double p_delta = 0;
		if(check_of<100)
		{
			r_delta = exp(check_of);
			p_delta = r_delta/(1+r_delta);
		}
		else
		{
			p_delta = 1;
		}

		gen_type die_gen_delta(generator, distribution_type(p_delta));
		boost::generator_iterator<gen_type> die_delta(&die_gen_delta);
		delta_t = *die_delta++;
		delta_rb[i] = p_delta;
	
		// geno_t
		
		for(int s = 0; s < num_sub_err; s++)
		{
                int j = sub_ind[s];
				
				for(int k = 0; k < row_size_t[j]; k++)
				{
					double r_1_num = 0;
					double inv_prior = 0;
					if(m_missing[j][ind_m[j][k]] == 0)
					{
						r_1_num = geno[j][ind_m[j][k]]+error_m[j][k]-2*geno[j][ind_m[j][k]]*error_m[j][k];
						
						inv_prior = maf_info[ind_m[j][k]];
					}
					else
					{
						r_1_num = 0.5;
						inv_prior = (1-error_m[j][k])/error_m[j][k];
					}

					double r_1 = (1-r_1_num)/r_1_num;

					double r_2 = 1;
					
					if(delta_t>0)
					{
						r_2 = exp(0.5*( gamma_t[ind_m[j][k]]*(2*(w_t[j]-(mu_t[j] + beta_0_t + mu_cov_t[j] + ransum_t[j] - gamma_t[ind_m[j][k]]*geno_t[j][k])) - gamma_t[ind_m[j][k]]) ));
					}

					double ratio = r_1*r_2/inv_prior;
				
					
					gen_type die_gen11(generator, distribution_type(ratio/(1+ratio)));
					boost::generator_iterator<gen_type> die11(&die_gen11);
					int new_geno = *die11++;
					if(new_geno != geno_t[j][k])
					{
						mu_t[j] += (new_geno - geno_t[j][k])*gamma_t[ind_m[j][k]];
						geno_t[j][k] = new_geno;
					}
					// geno_0[j][k] += geno_t[j][k];
				}
        }

	}
	
	double beta_0_re = 0;
	// double tau_re = 0;
	double tau_ran_re = 0;
	double delta_re = 0;
	double delta_rb_re = 0;
	Vec cov_re(NUM_COV);
	Vec gamma_re(NUM_RV);
	Vec gamma_sd(NUM_RV);
	Vec tau_g_re(NUM_RV);

	for(int i = 0; i < NUM_RV; i++)
	{
		gamma_re[i] = 0;
		gamma_sd[i] = 0;
		tau_g_re[i] = 0;
	}
	for(int i = 0; i < NUM_COV; i++)
	{
		cov_re[i] = 0;
	}

	// int l_t_z = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_re += beta_0[i];
		// tau_re += tau[i];
		delta_re += delta[i];
		delta_rb_re += delta_rb[i];
		/*
		for(int j = 0; j < NUM_RV; j++)
		{
			gamma_re[j] += gamma_0[i][j];
			tau_g_re[j] += tau_g[i][j];
		}
		*/
		for(int j = 0; j < NUM_COV; j++)
		{
			cov_re[j] += cov_z[i][j];
		}
	}

	if((*fam)==1)
	{
		for(int i = burnin; i < ITER_NUM; i++)
		{
			tau_ran_re += tau_ran[i];
		}
	}

	double iter_count = ITER_NUM - burnin;
	beta_0_re /= iter_count;
	// tau_re /= iter_count;
	tau_ran_re /= iter_count;
	/*
	for(int j = 0; j < NUM_RV; j++)
	{
		gamma_re[j] /= iter_count;
		tau_g_re[j] /= iter_count;
	}
	*/
	for(int j = 0; j < NUM_COV; j++)
	{
		cov_re[j] /= iter_count;
	}

	/*double equ = 0.01;
	for(int i = ITER_NUM/10; i < ITER_NUM; i++)
	{
		for(int j = 0; j < NUM_RV; j++)
		{
			gamma_sd[j] += (gamma_0[i][j] - gamma_re[j])*(gamma_0[i][j] - gamma_re[j]);
		}
	}

	for(int j = 0; j < NUM_RV; j++)
	{
		gamma_sd[j] /= delta_re;
		gamma_sd[j] = sqrt(gamma_sd[j]);
	}*/
	
	delta_re /= iter_count;
	delta_rb_re /= iter_count;

	result[0] = delta_re;
	result[1] = delta_rb_re;
	result[2] = beta_0_re;
	// result[3] = 1/tau_re;
	result[3] = total_miss;
	if((*fam)==1)
	{
		result[4] = 1/tau_ran_re;
	}
	else
	{
		result[4] = 0;
	}

	for(int j = 0; j < NUM_COV; j++)
	{
		result[5+j] = cov_re[j];
	}

}
