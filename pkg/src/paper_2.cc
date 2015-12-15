#include "paper_2.h"

typedef boost::minstd_rand base_generator_type;

void gibbssampler(double *result,
int * numRows, int * numCols, int * numCols2,
double *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g) 
{

	int ITER_NUM = 30000;
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
	typedef boost::exponential_distribution<> distribution_type_ex;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;
	typedef boost::variate_generator<base_generator_type&, distribution_type_ex> gen_type_ex;

	// read input files

	Vec pheno(0);
	int NUM_SUB = (*numRows);
	Mat_i geno;
	int NUM_RV = (*numCols);
	Mat cov_m;
	int NUM_COV = (*numCols2);
	Vec maf_info;
	int cov_yes = 0;
	int maf_yes = 0;
	Vec_i fam_no;
	Vec_i zyg;
	Mat_i fam;
	Mat_i mz;
	Mat_i dz;
	int NUM_FAM = 0;
	int NUM_MZ = 0;
	int NUM_DZ = 0;
	int IND_DZ = 0;
	int model = 1;

	double tau_0 = 0.04;
	double thres = 0.2;
	double initial_1 = (*arg_i_b);
	double initial_2 = (*arg_i_g);

		if((*arg_v) > 0)
		{
			tau_0 = 1/(*arg_v);
		}

		if((*arg_i) > 0)
		{
			thres = (*arg_i);
		}

		if((*arg_n) > 0)
		{
			ITER_NUM = (*arg_n);
		}

		if((*arg_b) > 0)
		{
			burnin = (*arg_b);
		}
		
		double * p = mat1;
			for(int i = 0; i < NUM_SUB; i++)
			{
				double temp = *p++;
				pheno.push_back(temp);
				temp = *p++;
				fam_no.push_back(short(temp));
				temp = *p++;
				zyg.push_back(short(temp));
			}
			
			for(int j = 0; j < NUM_SUB; j++)
			{
				int temp = fam_no[j];
				int i = 0;
				for(i = 0; i < fam.size(); i++)
				{
					if(temp == fam[i][0])
						break;
				}
				if( i == fam.size())
				{
					Vec_i new_fam;
					new_fam.push_back(temp);
					new_fam.push_back(j);
					fam.push_back(new_fam);
					if(zyg[j]==1)
					{
						mz.push_back(new_fam);
					}
					else
					{
						dz.push_back(new_fam);
						IND_DZ++;
					}
				}
				else
				{
					fam[i].push_back(j);
					int k = 0;
					if(zyg[j]==1)
					{
						for(k = 0; k < mz.size(); k++)
						{
							if(mz[k][0]==temp)
								break;
						}
						mz[k].push_back(j);
					}
					else
					{
						for(k = 0; k < dz.size(); k++)
						{
							if(dz[k][0]==temp)
								break;
						}
						dz[k].push_back(j);
						IND_DZ++;
					}
				}
			}
			NUM_FAM = fam.size();
			NUM_MZ = mz.size();
			NUM_DZ = dz.size();
		
		int * p2 = mat2;
		for(int i = 0; i < NUM_SUB; i++)
		{
			Vec_i row_ge;
			for(int j = 0; j < NUM_RV; j++)
			{
				int temp = *p2++;
				row_ge.push_back(temp);
			}
			geno.push_back(row_ge);
		}

		
		// covariates		

		if(NUM_COV > 0)
		{		
			p = mat3;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				for(int j = 0; j < NUM_COV; j++)
				{
					double temp = *p++;
					row_cov.push_back(temp);
				}
				cov_m.push_back(row_cov);
			}
			cov_yes = 1;
		}
		
		//MAF
		/*
		if((argv[arg_i][0]=='-')&&(argv[arg_i][1]=='m'))
		{
			std::fstream filestr5;
			filestr5.open(argv[arg_i+1], std::fstream::in);
			if(filestr5.fail())
			{
				cout<<"Please input the correct maf file."<<endl;
				return 0;
			}

			for(int i = 0; i < NUM_RV; i++)
			{
				double temp;
				filestr5 >> temp;
				temp = (1-2*temp)/(2*temp);
				maf_info.push_back(temp);
			}
			filestr5.close();
			maf_yes = 1;
		}
		*/

		if((*arg_t) > 0)
		{
			model = (*arg_t);
		}


	// Gibbs sampling

	// construct x, z and pos

	Vec_i x(NUM_SUB);
	Vec_i z(NUM_SUB);
	Mat_i pos;
	Mat_i alpha_ind;
	Vec_i z_pos_mz;
	Vec_i z_pos_dz;
	int z_num_mz = 0;
	int z_num_dz = 0;

	for(int i = 0; i < NUM_RV; i++)
	{
		Vec_i temp(0);
		alpha_ind.push_back(temp);
	}

	for(int i = 0; i < NUM_SUB; i++)
	{
		x[i] = 0;
		z[i] = 0;
		Vec_i pos_ind;

		int num_posi = 0;

		for(int j = 0; j < NUM_RV; j++)
		{
			if(geno[i][j] == 1)
			{
				num_posi++;
				pos_ind.push_back(j);
				alpha_ind[j].push_back(i);
			}
		}

		if(num_posi>0)
		{
			x[i] = 1;		
		}
		if(num_posi>1)
		{
			z[i] = 1;
			if(zyg[i]==1)
			{
				int k = 0;
				for(k = 0; k < NUM_MZ; k++)
				{
					if(mz[k][0]==fam_no[i])
						break;
				}
				int size_z_mz = z_pos_mz.size();
				int j = 0;
				for(; j < size_z_mz; j++)
				{
					if(z_pos_mz[j] == k)
						break;
				}
				if(j==size_z_mz)
				{
					z_pos_mz.push_back(k);
				}
			}
			else
			{
				z_pos_dz.push_back(i);
			}

		}
		pos.push_back(pos_ind);
	}

	z_num_mz = z_pos_mz.size();
	z_num_dz = z_pos_dz.size();

	// check mz genotype
	int check_mz = 0; 
	for(int i = 0; i < NUM_MZ; i++)
	{
		int mz_size = mz[i].size();
		int mz_x = x[mz[i][1]];
		int mz_z = z[mz[i][1]];
		
		for(int j = 1; j < mz_size; j++)
		{
			if(mz_x!=x[mz[i][j]])
				check_mz = 1;
			if(mz_z!=z[mz[i][j]])
				check_mz = 1;
		}
	}

	// estimate maf

	if(maf_yes == 0)
	{
		for(int i = 0; i < NUM_RV; i++)
		{
			int num_sub = 0;
			int num_min = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
					num_sub++;
					if(geno[j][i]==1)
						num_min++;
			}
			maf_info.push_back((num_sub-num_min)/num_min);
		}
	}

	
	Vec y(NUM_SUB);

	for(int i = 0; i < NUM_SUB; i++)
	{
		y[i] = pheno[i];
	}

	Vec beta_0(ITER_NUM+1);
	Vec beta_1(ITER_NUM+1);
	Vec beta_2(ITER_NUM+1);
	Vec tau(ITER_NUM+1,0);
	Vec tau_c(ITER_NUM+1,0);
	Vec tau_a(ITER_NUM+1,0);
	Vec_i sum_alpha(ITER_NUM+1);
	Vec_i sum_gamma(ITER_NUM+1);
	Mat_i alpha;
	Mat_i gamma;
	Mat cov_z;
	
	Vec_i alpha_t(NUM_RV);
	Mat p_t;
	Vec_i gamma_t(NUM_SUB);
	Mat q_t;
	Vec cov_t(NUM_COV);
	Vec mu_cov_t(NUM_SUB);
	// Vec mu_t(NUM_SUB);
	Vec var_c_t(NUM_SUB,0);
	Vec var_c_f_t(fam.size(),0);
	Vec var_a_m_t(NUM_SUB,0);
	Vec var_a_d1_t(NUM_SUB,0);
	Vec var_a_d2_t(NUM_SUB,0);
	Vec var_a_m_f_t(mz.size(),0);
	Vec var_a_d_f_t(dz.size(),0);

	// Initialization
	beta_0[0] = 0;
	beta_1[0] = initial_1;
	beta_2[0] = initial_2;
	tau[0] = 0.01;
	tau_c[0] = 0.01;
	tau_a[0] = 0.01;
	sum_alpha[0] = 0;
	sum_gamma[0] = 0;

	for(int i = 0; i < NUM_RV; i++)
	{
		alpha_t[i] = 0;
		Vec prob_p(3);
		prob_p[0] = 1.0/3;
		prob_p[1] = 1.0/3;
		prob_p[2] = 1.0/3;
		p_t.push_back(prob_p);
	}

	for(int i = 0; i < NUM_COV; i++)
	{
		cov_t[i] = 0;
	}
	
	alpha.push_back(alpha_t);
	if(NUM_COV>0)
		cov_z.push_back(cov_t);

	
	double p_a = 0.01;
	double p_b = 0.01;

	double beta_0_t = beta_0[0];
	double beta_1_t = beta_1[0];
	double beta_2_t = beta_2[0];
	double tau_t = tau[0];
	double tau_c_t = tau_c[0];
	double tau_a_t = tau_a[0];

	for(int i = 0; i < NUM_SUB; i++)
	{
		mu_cov_t[i] = 0;
		// mu_t[i] = 0;
		gamma_t[i] = 0;
		Vec prob_q(3);
		prob_q[0] = 1.0/3;
		prob_q[1] = 1.0/3;
		prob_q[2] = 1.0/3;
		q_t.push_back(prob_q);
	}

	gamma.push_back(gamma_t);

	// Iteration

	for(int i = 1; i < ITER_NUM + 1; i++)
	{

		// beta_0
		double temp = 0;
		int not_zero = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			double main = 0;
			double inter = 0;
			if(x[j] > 0)
			{
				int pos_sum = 0;
				int pos_size = pos[j].size();
				for(int k = 0; k < pos_size; k++)
				{
					pos_sum += alpha_t[pos[j][k]];
				}
				main = beta_1_t*pos_sum;
				if(pos_sum!=0)
				{
					not_zero = 1;
				}
			}
			if(z[j] > 0)
			{
				inter = beta_2_t*gamma_t[j];
			}

			temp += y[j] - main - inter - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j];
		}

		if(not_zero==0)
		{
			sum_alpha[i] = 0;
		}
		else
		{
			sum_alpha[i] = 1;
		}

		gen_type2 die_gen_b0(generator, distribution_type2((1/(NUM_SUB+tau_0/tau_t))*temp, sqrt(1/(NUM_SUB*tau_t+tau_0))));
		boost::generator_iterator<gen_type2> die_b0(&die_gen_b0);
		beta_0_t = *die_b0++;
	

		// beta_1

		double temp_1 = 0;
		double temp_2 = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			double inter = 0;
			
			if(z[j] > 0)
			{
				inter = beta_2_t*gamma_t[j];
			}

			if(x[j] > 0)
			{
				int pos_sum = 0;
				int pos_size = pos[j].size();
				for(int k = 0; k < pos_size; k++)
				{
					pos_sum += alpha_t[pos[j][k]];
				}
				temp_1 += (y[j] - beta_0_t - inter - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j])*pos_sum;
				temp_2 += pos_sum*pos_sum;
			}
		}

		double post_mu = (1/(temp_2+tau_0/tau_t))*temp_1;
		double post_sd = sqrt(1/(temp_2*tau_t+tau_0));

		if(post_mu > 0)
		{
			gen_type2 die_gen_b1(generator, distribution_type2(post_mu, post_sd));
			boost::generator_iterator<gen_type2> die_b1(&die_gen_b1);
			beta_1_t = *die_b1++;
			while(beta_1_t < 0)
			{
				beta_1_t = *die_b1++;
			}
		}
		else
		{/*
			double mu_minus =(-1)*(post_mu/post_sd);
			double alpha_st = (mu_minus + sqrt(mu_minus*mu_minus+4))/2;
			double u1 = 0;
			double z1 = 0;
			double rou_z1 = 1;
			gen_type_ex die_gen_z1(generator, distribution_type_ex(alpha_st));
			boost::generator_iterator<gen_type_ex> die_z1(&die_gen_z1);
			gen_type4 die_gen_u1(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
			while(u1 < rou_z1)
			{
				
				z1 = *die_z1++;
				z1 += mu_minus;
				rou_z1 = exp((-0.5)*(alpha_st-z1)*(alpha_st-z1));
								
				u1 = *die_u1++;
			}
		
			beta_1_t = z1*post_sd+post_mu;
			*/
			gen_type4 die_gen_u1(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
			gen_type2 die_gen_z(generator, distribution_type2(0, post_sd));
			boost::generator_iterator<gen_type2> die_z(&die_gen_z);
			
			double u1 = 2;
			double z = 0;
			while(u1>exp((-1)*z*post_mu*(temp_2*tau_t+tau_0)))
			{	
				z = *die_z++;
				while(z < 0)
				{
					z = *die_z++;
				}
					
				u1 = *die_u1++;
			}
			beta_1_t = z;
		}

		// beta_2

		temp_1 = 0;
		temp_2 = 0;
		int z_num = z_num_mz+z_num_dz;
		if(z_num > 0)
		{
			
			for(int j = 0; j < NUM_SUB; j++)
			{
				if(z[j] > 0)
				{
				
					double main = 0;
			
					int pos_sum = 0;
					int pos_size = pos[j].size();
					for(int k = 0; k < pos_size; k++)
					{
						pos_sum += alpha_t[pos[j][k]];
					}
					main = beta_1_t*pos_sum;
			
					temp_1 += (y[j] - beta_0_t - main - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j])*gamma_t[j];
					temp_2 += gamma_t[j]*gamma_t[j];
				}
			}

			post_mu = (1/(temp_2+tau_0/tau_t))*temp_1;
			post_sd = sqrt(1/(temp_2*tau_t+tau_0));

			if(post_mu > 0)
			{
				gen_type2 die_gen_b2(generator, distribution_type2(post_mu, post_sd));
				boost::generator_iterator<gen_type2> die_b2(&die_gen_b2);
				beta_2_t = *die_b2++;
				while(beta_2_t < 0)
				{
					beta_2_t = *die_b2++;
				}
			}
			else
			{/*
				double mu_minus =(-1)*(post_mu/post_sd);
				double alpha_st = (mu_minus + sqrt(mu_minus*mu_minus+4))/2;
				double u2 = 0;
				double z2 = 0;
				double rou_z2 = 1;
				gen_type_ex die_gen_z2(generator, distribution_type_ex(alpha_st));
				boost::generator_iterator<gen_type_ex> die_z2(&die_gen_z2);
				gen_type4 die_gen_u2(generator, distribution_type4());
				boost::generator_iterator<gen_type4> die_u2(&die_gen_u2);
				while(u2 < rou_z2)
				{			
					z2 = *die_z2++;
					z2 += mu_minus;
					rou_z2 = exp((-0.5)*(alpha_st-z2)*(alpha_st-z2));
				
					u2 = *die_u2++;
				}
		
				beta_2_t = z2*post_sd+post_mu;
				*/
				gen_type4 die_gen_u1(generator, distribution_type4());
				boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
				gen_type2 die_gen_z(generator, distribution_type2(0, post_sd));
				boost::generator_iterator<gen_type2> die_z(&die_gen_z);
			
				double u1 = 2;
				double z = 0;
				while(u1>exp((-1)*z*post_mu*(temp_2*tau_t+tau_0)))
				{	
					z = *die_z++;
					while(z < 0)
					{
						z = *die_z++;
					}
					
					u1 = *die_u1++;
				}
				beta_2_t = z;
			}
		}

		// P	

		gen_type3 die_gen_p(generator, distribution_type3(0.1, 1));
		boost::generator_iterator<gen_type3> die_p(&die_gen_p);
		gen_type3 die_gen_p2(generator, distribution_type3(1.1, 1));
		boost::generator_iterator<gen_type3> die_p2(&die_gen_p2);

		for(int j = 0; j < NUM_RV; j++)
		{
			double p_1 = 0;
			double p_2 = 0;
			double p_3 = 0;
			
			p_1 = *die_p++;
			p_2 = *die_p++;		
			p_3 = *die_p2++;

			double sum_p = p_1 + p_2 + p_3;

			if(alpha_t[j]==1)
			{
				p_t[j][0] = p_3/sum_p;
				p_t[j][1] = p_1/sum_p;
				p_t[j][2] = p_2/sum_p;
			}
			else
			{
				if(alpha_t[j]==0)
				{
					p_t[j][1] = p_3/sum_p;
					p_t[j][0] = p_1/sum_p;
					p_t[j][2] = p_2/sum_p;
				}
				else
				{
					p_t[j][2] = p_3/sum_p;
					p_t[j][0] = p_1/sum_p;
					p_t[j][1] = p_2/sum_p;
				}
			}
		}


		// alpha

		for(int j = 0; j < NUM_RV; j++)
		{
			double temp_a_1 = 0;
			double temp_a_2 = 0;
			double temp_a_3 = 0;
			int num_sub = alpha_ind[j].size();
			for(int k = 0; k < num_sub; k++)
			{
				double second_p = 0;
				double inter = 0;
				int ind_num = alpha_ind[j][k];
				if(z[ind_num]==1)
				{
					inter = beta_2_t*gamma_t[ind_num];
					int pos_sum = 0;
					int pos_size = pos[ind_num].size();
					for(int l = 0; l < pos_size; l++)
					{
						pos_sum += alpha_t[pos[ind_num][l]];
					}
					second_p = beta_1_t*(pos_sum-alpha_t[j]);
				}

				double diff_a_2 = y[ind_num] - beta_0_t - second_p - inter - mu_cov_t[ind_num] - var_c_t[ind_num] - var_a_m_t[ind_num] - var_a_d1_t[ind_num] - var_a_d2_t[ind_num];
				double diff_a_1 = diff_a_2 - beta_1_t;
				double diff_a_3 = diff_a_2 + beta_1_t;

				temp_a_1 += diff_a_1*diff_a_1;
				temp_a_2 += diff_a_2*diff_a_2;
				temp_a_3 += diff_a_3*diff_a_3;
			}

			temp_a_1 = p_t[j][0]*exp((-0.5)*tau_t*temp_a_1);
			temp_a_2 = p_t[j][1]*exp((-0.5)*tau_t*temp_a_2);
			temp_a_3 = p_t[j][2]*exp((-0.5)*tau_t*temp_a_3);
			double sum_temp_a = temp_a_1+temp_a_2+temp_a_3;

			gen_type4 die_gen_alpha(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_alpha(&die_gen_alpha);
			double unif_01 = *die_alpha++;
			if(unif_01 < (temp_a_1/sum_temp_a))
			{
				alpha_t[j] = 1;
			}
			else
			{
				if(unif_01 < ((temp_a_1+temp_a_2)/sum_temp_a))
					alpha_t[j] = 0;
				else
				{
					alpha_t[j] = -1;
				}
			}
		}

		// Q_dz

		gen_type3 die_gen_q(generator, distribution_type3(0.1, 1 ));
		boost::generator_iterator<gen_type3> die_q(&die_gen_q);
		gen_type3 die_gen_q2(generator, distribution_type3(1.1, 1 ));
		boost::generator_iterator<gen_type3> die_q2(&die_gen_q2);

		for(int j = 0; j < z_num_dz; j++)
		{
			double q_1 = 0;
			double q_2 = 0;
			double q_3 = 0;
			
			q_1 = *die_q++;
			q_2 = *die_q++;
			q_3 = *die_q2++;

			double sum_q = q_1 + q_2 + q_3;
			int sub_no = z_pos_dz[j];

			if(gamma_t[sub_no]==1)
			{
				q_t[sub_no][0] = q_3/sum_q;
				q_t[sub_no][1] = q_1/sum_q;
				q_t[sub_no][2] = q_2/sum_q;
			}
			else
			{
				if(gamma_t[sub_no]==0)
				{
					q_t[sub_no][1] = q_3/sum_q;
					q_t[sub_no][0] = q_1/sum_q;
					q_t[sub_no][2] = q_2/sum_q;
				}
				else
				{
					q_t[sub_no][2] = q_3/sum_q;
					q_t[sub_no][0] = q_1/sum_q;
					q_t[sub_no][1] = q_2/sum_q;
				}
			}
		}

		// Q_mz

		for(int j = 0; j < z_num_mz; j++)
		{
			double q_1 = 0;
			double q_2 = 0;
			double q_3 = 0;
			
			q_1 = *die_q++;
			q_2 = *die_q++;
			q_3 = *die_q2++;

			double sum_q = q_1 + q_2 + q_3;
			int mz_ind = z_pos_mz[j];

			if(gamma_t[mz[mz_ind][1]]==1)
			{
				q_t[mz[mz_ind][1]][0] = q_3/sum_q;
				q_t[mz[mz_ind][1]][1] = q_1/sum_q;
				q_t[mz[mz_ind][1]][2] = q_2/sum_q;
			}
			else
			{
				if(gamma_t[mz[mz_ind][1]]==0)
				{
					q_t[mz[mz_ind][1]][1] = q_3/sum_q;
					q_t[mz[mz_ind][1]][0] = q_1/sum_q;
					q_t[mz[mz_ind][1]][2] = q_2/sum_q;
				}
				else
				{
					q_t[mz[mz_ind][1]][2] = q_3/sum_q;
					q_t[mz[mz_ind][1]][0] = q_1/sum_q;
					q_t[mz[mz_ind][1]][1] = q_2/sum_q;
				}
			}
		}

		// gamma_dz

		int flag_gamma = 0;
		for(int j = 0; j < z_num_dz; j++)
		{
			double main = 0;

			int pos_sum = 0;
			int sub_no = z_pos_dz[j];
			int pos_size = pos[sub_no].size();
			for(int k = 0; k < pos_size; k++)
			{
				pos_sum += alpha_t[pos[sub_no][k]];
			}
			main = beta_1_t*pos_sum;

			
			double temp_g_2 = y[sub_no] - beta_0_t - main - mu_cov_t[sub_no] - var_c_t[sub_no] - var_a_d1_t[sub_no] - var_a_d2_t[sub_no];
			double temp_g_1 = temp_g_2 - beta_2_t;
			double temp_g_3 = temp_g_2 + beta_2_t;
			temp_g_1 *= temp_g_1;
			temp_g_2 *= temp_g_2;
			temp_g_3 *= temp_g_3;

			temp_g_1 = q_t[sub_no][0]*exp((-0.5)*tau_t*temp_g_1);
			temp_g_2 = q_t[sub_no][1]*exp((-0.5)*tau_t*temp_g_2);
			temp_g_3 = q_t[sub_no][2]*exp((-0.5)*tau_t*temp_g_3);
			double sum_temp_g = temp_g_1+temp_g_2+temp_g_3;

			gen_type4 die_gen_gamma(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_gamma(&die_gen_gamma);
			double unif_01 = *die_gamma++;
			if(unif_01 < (temp_g_1/sum_temp_g))
			{
				gamma_t[sub_no] = 1;
				flag_gamma = 1;
			}
			else
			{
				if(unif_01 < ((temp_g_1+temp_g_2)/sum_temp_g))
					gamma_t[sub_no] = 0;
				else
				{
					gamma_t[sub_no] = -1;
					flag_gamma = 1;
				}
			}
		}

		// gamma_mz

		for(int j = 0; j < z_num_mz; j++)
		{
			double main = 0;

			int pos_sum = 0;
			int mz_ind = z_pos_mz[j];
			int pos_size = pos[mz[mz_ind][1]].size();
			for(int l = 0; l < pos_size; l++)
			{
				pos_sum += alpha_t[pos[mz[mz_ind][1]][l]];
			}
			main = beta_1_t*pos_sum;

			double temp_gm_1 = 0;
			double temp_gm_2 = 0;
			double temp_gm_3 = 0;
			
			double sum_1 = 0;
			double sum_2 = 0;
			double sum_3 = 0;
			for(int k = 1; k < mz[mz_ind].size(); k++)
			{
				temp_gm_2 = y[mz[mz_ind][k]] - beta_0_t - main - mu_cov_t[mz[mz_ind][k]] - var_c_t[mz[mz_ind][k]] - var_a_m_t[mz[mz_ind][k]] ;
				temp_gm_1 = temp_gm_2 - beta_2_t;
				temp_gm_3 = temp_gm_2 + beta_2_t;
				temp_gm_1 *= temp_gm_1;
				temp_gm_2 *= temp_gm_2;
				temp_gm_3 *= temp_gm_3;
				sum_1 += temp_gm_1;
				sum_2 += temp_gm_2;
				sum_3 += temp_gm_3;
			}

			sum_1 = q_t[mz[mz_ind][1]][0]*exp((-0.5)*tau_t*sum_1);
			sum_2 = q_t[mz[mz_ind][1]][1]*exp((-0.5)*tau_t*sum_2);
			sum_3 = q_t[mz[mz_ind][1]][2]*exp((-0.5)*tau_t*sum_3);
			double sum_temp = sum_1+sum_2+sum_3;

			gen_type4 die_gen_gamma_mz(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_gamma_mz(&die_gen_gamma_mz);
			double unif_01 = *die_gamma_mz++;
			if(unif_01 < (sum_1/sum_temp))
			{
				for(int k = 1; k < mz[mz_ind].size(); k++)
					gamma_t[mz[mz_ind][k]] = 1;
				flag_gamma = 1;
			}
			else
			{
				if(unif_01 < ((sum_1+sum_2)/sum_temp))
				{
					for(int k = 1; k < mz[mz_ind].size(); k++)
						gamma_t[mz[mz_ind][k]] = 0;
				}
				else
				{
					for(int k = 1; k < mz[mz_ind].size(); k++)
						gamma_t[mz[mz_ind][k]] = -1;
					flag_gamma = 1;
				}
			}
		}
		if(flag_gamma == 1)
		{
			sum_gamma[i] = 1;
		}
		else
		{
			sum_gamma[i] = 0;
		}


		// cov

		for(int j = 0; j < NUM_COV; j++)
		{
			double new_cov = 0;
			double g_sq = 0;
			double sum_g_g = 0;

			for(int k = 0; k < NUM_SUB; k++)
			{
				int pos_sum = 0;
				int pos_size = pos[k].size();
				for(int l = 0; l < pos_size; l++)
				{
					pos_sum += alpha_t[pos[k][l]];
				}
				sum_g_g += (y[k]-beta_1_t*pos_sum*x[k]-beta_2_t*gamma_t[k]*z[k]-beta_0_t-mu_cov_t[k]+cov_t[j]*cov_m[k][j]-var_c_t[k]- var_a_m_t[k] - var_a_d1_t[k] - var_a_d2_t[k])*cov_m[k][j];
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
		
		
		// tau

		double temp4 = 0;
		
		for(int j = 0; j < NUM_SUB; j++)
		{
			int pos_sum = 0;
			int pos_size = pos[j].size();
			for(int k = 0; k < pos_size; k++)
			{
				pos_sum += alpha_t[pos[j][k]];
			}
			double y_b0_b1_ax = y[j] - beta_1_t*pos_sum*x[j]-beta_2_t*gamma_t[j]*z[j] - beta_0_t - mu_cov_t[j] - var_c_t[j]- var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j];
			temp4 += y_b0_b1_ax*y_b0_b1_ax;
		}
		
		gen_type3 die_gen9(generator, distribution_type3(0.001+0.5*NUM_SUB, 1/(0.001+0.5*temp4) ));
		boost::generator_iterator<gen_type3> die9(&die_gen9);
		tau_t = *die9++;

		if(model == 3)
		{
			// var_c

			for(int j = 0; j < NUM_FAM; j++)
			{
				double temp = 0;
				int size = fam[j].size();
				for(int k = 1; k < size; k++)
				{
					int sub = fam[j][k];
					double main = 0;
					double inter = 0;
					if(x[sub] > 0)
					{
						int pos_sum = 0;
						int pos_size = pos[sub].size();
						for(int l = 0; l < pos_size; l++)
						{
							pos_sum += alpha_t[pos[sub][l]];
						}
						main = beta_1_t*pos_sum;
					}
					if(z[sub] > 0)
					{
						inter = beta_2_t*gamma_t[sub];
					}

					temp += y[sub] - beta_0_t - main - inter - mu_cov_t[sub] - var_a_m_t[sub] - var_a_d1_t[sub] - var_a_d2_t[sub];
				}

				gen_type2 die_gen_var_c(generator, distribution_type2((1/(size-1+tau_c_t/tau_t))*temp, sqrt(1/((size-1)*tau_t+tau_c_t))));
				boost::generator_iterator<gen_type2> die_var_c(&die_gen_var_c);
				var_c_f_t[j] = *die_var_c++;
				for(int k = 1; k < size; k++)
				{
					int sub = fam[j][k];
					var_c_t[sub] = var_c_f_t[j];
				}
			}

			// tau_c

			double temp_c = 0;
		
			for(int j = 0; j < NUM_FAM; j++)
			{
				temp_c += var_c_f_t[j]*var_c_f_t[j];
			}
		
			gen_type3 die_gen_tau_c(generator, distribution_type3(0.001+0.5*NUM_FAM, 1/(0.001+0.5*temp_c) ));
			boost::generator_iterator<gen_type3> die_tau_c(&die_gen_tau_c);
			tau_c_t = *die_tau_c++;

		}

		if(model > 1)
		{
			double temp_a = 0;

			// var_a_m

			for(int j = 0; j < NUM_MZ; j++)
			{
				double temp = 0;
				int size = mz[j].size();
				for(int k = 1; k < size; k++)
				{
					int sub = mz[j][k];
					double main = 0;
					double inter = 0;
					if(x[sub] > 0)
					{
						int pos_sum = 0;
						int pos_size = pos[sub].size();
						for(int l = 0; l < pos_size; l++)
						{
							pos_sum += alpha_t[pos[sub][l]];
						}
						main = beta_1_t*pos_sum;
					}
					if(z[sub] > 0)
					{
						inter = beta_2_t*gamma_t[sub];
					}

					temp += y[sub] - beta_0_t - main - inter - mu_cov_t[sub] - var_c_t[sub];
				}

				gen_type2 die_gen_var_a_m(generator, distribution_type2((1/(size-1+tau_a_t/tau_t))*temp, sqrt(1/((size-1)*tau_t+tau_a_t))));
				boost::generator_iterator<gen_type2> die_var_a_m(&die_gen_var_a_m);
				var_a_m_f_t[j] = *die_var_a_m++;
				temp_a += var_a_m_f_t[j]*var_a_m_f_t[j];
				for(int k = 1; k < size; k++)
				{
					int sub = mz[j][k];
					var_a_m_t[sub] = var_a_m_f_t[j];
				}
			}

			// var_a_d

			for(int j = 0; j < NUM_DZ; j++)
			{
				double temp = 0;
				int size = dz[j].size();
				Vec mu_temp(size,0);
				for(int k = 1; k < size; k++)
				{
					int sub = dz[j][k];
					double main = 0;
					double inter = 0;
					if(x[sub] > 0)
					{
						int pos_sum = 0;
						int pos_size = pos[sub].size();
						for(int l = 0; l < pos_size; l++)
						{
							pos_sum += alpha_t[pos[sub][l]];
						}
						main = beta_1_t*pos_sum;
					}
					if(z[sub] > 0)
					{
						inter = beta_2_t*gamma_t[sub];
					}

					mu_temp[k] = beta_0_t + main + inter;

					temp += y[sub] - mu_temp[k] - mu_cov_t[sub] - var_c_t[sub] - var_a_d2_t[sub];
				}

				gen_type2 die_gen_var_a_d1(generator, distribution_type2((1/(size-1+2*tau_a_t/tau_t))*temp, sqrt(1/((size-1)*tau_t+2*tau_a_t))));
				boost::generator_iterator<gen_type2> die_var_a_d1(&die_gen_var_a_d1);
				var_a_d_f_t[j] = *die_var_a_d1++;
				temp_a += 2*var_a_d_f_t[j]*var_a_d_f_t[j];
				for(int k = 1; k < size; k++)
				{
					int sub = dz[j][k];
					var_a_d1_t[sub] = var_a_d_f_t[j];
				}

				for(int k = 1; k < size; k++)
				{
					int sub = dz[j][k];					
					temp = y[sub] - mu_temp[k] - mu_cov_t[sub] - var_c_t[sub] - var_a_d1_t[sub];

					gen_type2 die_gen_var_a_d2(generator, distribution_type2((1/(1+2*tau_a_t/tau_t))*temp, sqrt(1/(tau_t+2*tau_a_t))));
					boost::generator_iterator<gen_type2> die_var_a_d2(&die_gen_var_a_d2);
					var_a_d2_t[sub] = *die_var_a_d2++;
					temp_a += 2*var_a_d2_t[sub]*var_a_d2_t[sub];
				}
			}

			// tau_a

			gen_type3 die_gen_tau_a(generator, distribution_type3(0.001+0.5*(NUM_MZ+NUM_DZ+IND_DZ), 1/(0.001+0.5*temp_a) ));
			boost::generator_iterator<gen_type3> die_tau_a(&die_gen_tau_a);
			tau_a_t = *die_tau_a++;
		}

		beta_0[i] = beta_0_t;
		beta_1[i] = beta_1_t;
		beta_2[i] = beta_2_t;
		alpha.push_back(alpha_t);
		tau[i] = tau_t;
		if(model > 1)
		{
			tau_a[i] = tau_a_t;
			if(model > 2)
			{
				tau_c[i] = tau_c_t;
			}
		}
		if(NUM_COV>0)
			cov_z.push_back(cov_t);

		
	}
	
	double beta_0_re = 0;
	double beta_1_re = 0;
	double beta_2_re = 0;
	double tau_re = 0;
	double tau_c_re = 0;
	double tau_a_re = 0;
	int beta_1_posi = 0;
	int beta_2_posi = 0;
	int zero_alpha = 0;
	int zero_gamma = 0;
	
	Vec alpha_re(NUM_RV);
	Vec cov_re(NUM_COV);

	for(int i = 0; i < NUM_RV; i++)
	{
		alpha_re[i] = 0;
	}
	for(int i = 0; i < NUM_COV; i++)
	{
		cov_re[i] = 0;
	}

	int l_t_z = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_re += beta_0[i];
		beta_1_re += beta_1[i];
		beta_2_re += beta_2[i];
		tau_re += tau[i];
		tau_c_re += tau_c[i];
		tau_a_re += tau_a[i];

		zero_alpha += sum_alpha[i];
		zero_gamma += sum_gamma[i];
		
		if((sum_alpha[i] == 1)&&(beta_1[i]>thres))
			beta_1_posi++;
		if((sum_gamma[i-1] == 1)&&(beta_2[i]>thres))
			beta_2_posi++;

		for(int j = 0; j < NUM_RV; j++)
		{
			alpha_re[j] += alpha[i][j];
		}
		for(int j = 0; j < NUM_COV; j++)
		{
			cov_re[j] += cov_z[i][j];
		}
	}

	double iter_count = ITER_NUM - burnin;
	beta_0_re /= iter_count;
	beta_1_re /= iter_count;
	beta_2_re /= iter_count;
	tau_re /= iter_count;
	tau_c_re /= iter_count;
	tau_a_re /= iter_count;

	double zero_alpha_re = 1 - double(zero_alpha)/iter_count;
	double zero_gamma_re = 1 - double(zero_gamma)/iter_count;

	double beta_0_sd = 0;
	double beta_1_sd = 0;
	double beta_2_sd = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_sd += (beta_0[i] - beta_0_re)*(beta_0[i] - beta_0_re);
		beta_1_sd += (beta_1[i] - beta_1_re)*(beta_1[i] - beta_1_re);
		beta_2_sd += (beta_2[i] - beta_2_re)*(beta_2[i] - beta_2_re);
	}
	beta_0_sd /= iter_count;
	beta_1_sd /= iter_count;
	beta_2_sd /= iter_count;

	beta_0_sd = sqrt(beta_0_sd);
	beta_1_sd = sqrt(beta_1_sd);
	beta_2_sd = sqrt(beta_2_sd);
	
	//beta_1_po /= (1-0.01128341555)/0.01128341555;
	result[0] = beta_0_re;
	result[1] = beta_1_re;
	result[2] = beta_2_re;
	result[3] = 1/tau_re;
	result[4] = beta_1_posi;
	result[5] = beta_2_posi;
	result[6] = z_num_mz+z_num_dz;
	result[7] = check_mz;
	result[8] = 1/tau_a_re;
	result[9] = 1/tau_c_re;
	result[10] = beta_0_sd;
	result[11] = beta_1_sd;
	result[12] = beta_2_sd;

	/*
	cout<<beta_0_re<<endl;
	cout<<beta_1_re<<endl;
	cout<<beta_2_re<<endl;
	cout<<1/tau_re<<endl;
	cout<<beta_1_posi<<endl;
	cout<<beta_2_posi<<endl;
	cout<<z_num_mz+z_num_dz<<endl;
	cout<<check_mz<<endl;
	cout<<1/tau_a_re<<endl;
	cout<<1/tau_c_re<<endl;

    for(int i = 0; i < NUM_RV; i++)
    {
        //cout<<gamma_re[i]<<endl;
        //cout<<gamma_sd[i]<<endl;
		//cout<<gamma_posi[i]<<endl;
    }*/
	for(int j = 0; j < NUM_RV; j++)
	{
		alpha_re[j] /= iter_count;
		result[13+j] = alpha_re[j];
	}
	for(int j = 0; j < NUM_COV; j++)
	{
		cov_re[j] /= iter_count;
		result[13+NUM_RV+j] = cov_re[j];
	}
	/*
	for(int i = 0; i < NUM_COV; i++)
    {
        result[10+i] = cov_re[i];

	// cout<<cov_re[i]<<endl;
        //cout<<gamma_sd[i]<<endl;
		//cout<<gamma_posi[i]<<endl;
    }
	*/
}

void gibbssampler_bin(double *result,
int * numRows, int * numCols, int * numCols2,
int *mat1, int *mat2, double * mat3, double * arg_v, double * arg_i, int * arg_n, int * arg_b, int * arg_t, double * arg_i_b, double * arg_i_g) 
{

	int ITER_NUM = 10000;
	int burnin = 500;

	const double PI = 3.141592653589793;

	base_generator_type generator(24);

	typedef std::vector<double> Vec;
	typedef std::vector<short> Vec_i;
	typedef std::vector<Vec> Mat;
	typedef std::vector<Vec_i> Mat_i;

	typedef boost::bernoulli_distribution<> distribution_type;
	typedef boost::normal_distribution<> distribution_type2;
	typedef boost::gamma_distribution<> distribution_type3;
	typedef boost::uniform_01<> distribution_type4;
	typedef boost::exponential_distribution<> distribution_type_ex;

	typedef boost::variate_generator<base_generator_type&, distribution_type> gen_type;
	typedef boost::variate_generator<base_generator_type&, distribution_type2> gen_type2;
	typedef boost::variate_generator<base_generator_type&, distribution_type3> gen_type3;
	typedef boost::variate_generator<base_generator_type&, distribution_type4> gen_type4;
	typedef boost::variate_generator<base_generator_type&, distribution_type_ex> gen_type_ex;

	// read input files

	Vec_i pheno(0);
	int NUM_SUB = (*numRows);
	Mat_i geno;
	int NUM_RV = (*numCols);
	Mat cov_m;
	int NUM_COV = (*numCols2);
	Vec maf_info;
	int cov_yes = 0;
	int maf_yes = 0;
	Vec_i fam_no;
	Vec_i zyg;
	Mat_i fam;
	Mat_i mz;
	Mat_i dz;
	int NUM_FAM = 0;
	int NUM_MZ = 0;
	int NUM_DZ = 0;
	int IND_DZ = 0;
	int model = 1;

	// default prior var, lambda, initials
	double tau_0 = 1;
	double thres = 0.2;
	double initial_1 = (*arg_i_b);
	double initial_2 = (*arg_i_g);

	// beta_1 must be non-negative
	if(initial_1<0)
		initial_1 = 0.5;

	if((*arg_v) > 0)
		{
			tau_0 = 1/(*arg_v);
		}

		if((*arg_i) > 0)
		{
			thres = (*arg_i);
		}

		if((*arg_n) > 0)
		{
			ITER_NUM = (*arg_n);
		}

		if((*arg_b) > 0)
		{
			burnin = (*arg_b);
		}
		
		int * p = mat1;
			for(int i = 0; i < NUM_SUB; i++)
			{
				int temp = *p++;
				pheno.push_back(temp);
				temp = *p++;
				fam_no.push_back(short(temp));
				temp = *p++;
				zyg.push_back(short(temp));
			}
			
			for(int j = 0; j < NUM_SUB; j++)
			{
				int temp = fam_no[j];
				int i = 0;
				for(i = 0; i < fam.size(); i++)
				{
					if(temp == fam[i][0])
						break;
				}
				if( i == fam.size())
				{
					Vec_i new_fam;
					new_fam.push_back(temp);
					new_fam.push_back(j);
					fam.push_back(new_fam);
					if(zyg[j]==1)
					{
						mz.push_back(new_fam);
					}
					else
					{
						dz.push_back(new_fam);
						IND_DZ++;
					}
				}
				else
				{
					fam[i].push_back(j);
					int k = 0;
					if(zyg[j]==1)
					{
						for(k = 0; k < mz.size(); k++)
						{
							if(mz[k][0]==temp)
								break;
						}
						mz[k].push_back(j);
					}
					else
					{
						for(k = 0; k < dz.size(); k++)
						{
							if(dz[k][0]==temp)
								break;
						}
						dz[k].push_back(j);
						IND_DZ++;
					}
				}
			}
			NUM_FAM = fam.size();
			NUM_MZ = mz.size();
			NUM_DZ = dz.size();
		
		int * p2 = mat2;
		for(int i = 0; i < NUM_SUB; i++)
		{
			Vec_i row_ge;
			for(int j = 0; j < NUM_RV; j++)
			{
				int temp = *p2++;
				row_ge.push_back(temp);
			}
			geno.push_back(row_ge);
		}

		
		// covariates

		double * p3 = 0;

		if(NUM_COV > 0)
		{		
			p3 = mat3;
			for(int i = 0; i < NUM_SUB; i++)
			{
				Vec row_cov;
				for(int j = 0; j < NUM_COV; j++)
				{
					double temp = *p3++;
					row_cov.push_back(temp);
				}
				cov_m.push_back(row_cov);
			}
			cov_yes = 1;
		}
		
		if((*arg_t) > 0)
		{
			model = (*arg_t);
		}

	double gamma_a_c_a = 0.1;
	double gamma_a_c_b = 0.1;

	// Gibbs sampling

	// construct x, z and pos

	Vec_i x(NUM_SUB);
	Vec_i z(NUM_SUB);
	Mat_i pos;
	Mat_i alpha_ind;
	Vec_i z_pos_mz;
	Vec_i z_pos_dz;
	int z_num_mz = 0;
	int z_num_dz = 0;

	for(int i = 0; i < NUM_RV; i++)
	{
		Vec_i temp(0);
		alpha_ind.push_back(temp);
	}

	for(int i = 0; i < NUM_SUB; i++)
	{
		x[i] = 0;
		z[i] = 0;
		Vec_i pos_ind;

		int num_posi = 0;

		for(int j = 0; j < NUM_RV; j++)
		{
			if(geno[i][j] == 1)
			{
				num_posi++;
				pos_ind.push_back(j);
				alpha_ind[j].push_back(i);
			}
		}

		if(num_posi>0)
		{
			x[i] = 1;		
		}
		if(num_posi>1)
		{
			z[i] = 1;
			if(zyg[i]==1)
			{
				int k = 0;
				for(k = 0; k < NUM_MZ; k++)
				{
					if(mz[k][0]==fam_no[i])
						break;
				}
				int size_z_mz = z_pos_mz.size();
				int j = 0;
				for(; j < size_z_mz; j++)
				{
					if(z_pos_mz[j] == k)
						break;
				}
				if(j==size_z_mz)
				{
					z_pos_mz.push_back(k);
				}
			}
			else
			{
				z_pos_dz.push_back(i);
			}

		}
		pos.push_back(pos_ind);
	}

	z_num_mz = z_pos_mz.size();
	z_num_dz = z_pos_dz.size();

	// check mz genotype
	int check_mz = 0; 
	for(int i = 0; i < NUM_MZ; i++)
	{
		int mz_size = mz[i].size();
		int mz_x = x[mz[i][1]];
		int mz_z = z[mz[i][1]];
		
		for(int j = 1; j < mz_size; j++)
		{
			if(mz_x!=x[mz[i][j]])
				check_mz = 1;
			if(mz_z!=z[mz[i][j]])
				check_mz = 1;
		}
	}

	// estimate maf

	if(maf_yes == 0)
	{
		for(int i = 0; i < NUM_RV; i++)
		{
			int num_sub = 0;
			int num_min = 0;
			for(int j = 0; j < NUM_SUB; j++)
			{
					num_sub++;
					if(geno[j][i]==1)
						num_min++;
			}
			maf_info.push_back((num_sub-num_min)/num_min);
		}
	}

	
	Vec_i y(NUM_SUB);

	for(int i = 0; i < NUM_SUB; i++)
	{
		y[i] = pheno[i];
	}

	Vec beta_0(ITER_NUM+1);
	Vec beta_1(ITER_NUM+1);
	Vec beta_2(ITER_NUM+1);
	// Vec tau(ITER_NUM+1,0);
	Vec tau_c(ITER_NUM+1,0);
	Vec tau_a(ITER_NUM+1,0);
	Vec_i sum_alpha(ITER_NUM+1);
	Vec_i sum_gamma(ITER_NUM+1);
	Mat_i alpha;
	Mat_i gamma;
	Mat cov_z;
	/*
	double tau_d = 10000;
	double tau_d_est = 0;
	int calib = 0;
	double relative = abs(tau_d-tau_d_est)/tau_d;

	while(relative > 0.1)
	{
	tau_d_est = 0;
	alpha.clear();
	cov_z.clear();
	gamma.clear();
	*/
	Vec_i alpha_t(NUM_RV);
	Mat p_t;
	Vec_i gamma_t(NUM_SUB);
	Mat q_t;
	Vec cov_t(NUM_COV);
	Vec mu_cov_t(NUM_SUB);
	// Vec mu_t(NUM_SUB);
	Vec var_c_t(NUM_SUB,0);
	Vec var_c_f_t(fam.size(),0);
	Vec var_a_m_t(NUM_SUB,0);
	Vec var_a_d1_t(NUM_SUB,0);
	Vec var_a_d2_t(NUM_SUB,0);
	Vec var_a_m_f_t(mz.size(),0);
	Vec var_a_d_f_t(dz.size(),0);

	// augment variable
	Vec w_t(NUM_SUB,0);

	// Initialization
	beta_0[0] = 0;
	beta_1[0] = initial_1;
	beta_2[0] = initial_2;
	// tau[0] = 0.01;
	tau_c[0] = 1;
	tau_a[0] = 1;
	sum_alpha[0] = 0;
	sum_gamma[0] = 0;

	for(int i = 0; i < NUM_RV; i++)
	{
		alpha_t[i] = 1;
		Vec prob_p(3);
		prob_p[0] = 1.0/3;
		prob_p[1] = 1.0/3;
		prob_p[2] = 1.0/3;
		p_t.push_back(prob_p);
	}

	for(int i = 0; i < NUM_COV; i++)
	{
		cov_t[i] = 0;
	}
	
	alpha.push_back(alpha_t);
	if(NUM_COV>0)
		cov_z.push_back(cov_t);

	
	double p_a = 0.01;
	double p_b = 0.01;

	double beta_0_t = beta_0[0];
	double beta_1_t = beta_1[0];
	double beta_2_t = beta_2[0];
	// double tau_t = tau[0];
	double tau_c_t = tau_c[0];
	double tau_a_t = tau_a[0];
	double tau_dz_t = 1;

	for(int i = 0; i < NUM_SUB; i++)
	{
		mu_cov_t[i] = 0;
		// mu_t[i] = 0;
		// gamma_t[i] = 0;
		gamma_t[i] = 1;
		Vec prob_q(3);
		prob_q[0] = 1.0/3;
		prob_q[1] = 1.0/3;
		prob_q[2] = 1.0/3;
		q_t.push_back(prob_q);
	}

	gamma.push_back(gamma_t);

	int iter_trace = 0;

	// Iteration

	for(int i = 1; i < ITER_NUM + 1; i++)
	{
		iter_trace = i;

		// beta_0
		double temp = 0;
		int not_zero = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			double main = 0;
			double inter = 0;
			if(x[j] > 0)
			{
				int pos_sum = 0;
				int pos_size = pos[j].size();
				for(int k = 0; k < pos_size; k++)
				{
					pos_sum += alpha_t[pos[j][k]];
				}
				main = beta_1_t*pos_sum;
				if(pos_sum!=0)
				{
					not_zero = 1;
					sum_alpha[i] = 1;
				}
				else
				{
					sum_alpha[i] = 0;
				}
			}
			
			if(z[j] > 0)
			{
				inter = beta_2_t*gamma_t[j];
			}
			
			temp += w_t[j] - main - inter - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j];
			

			// temp += w_t[j] - main - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j];
		}
		
		if(not_zero==0)
		{
			sum_alpha[i] = 0;
		}
		else
		{
			sum_alpha[i] = 1;
		}
		
		gen_type2 die_gen_b0(generator, distribution_type2((1/(NUM_SUB+tau_0))*temp, sqrt(1/(NUM_SUB+tau_0))));
		boost::generator_iterator<gen_type2> die_b0(&die_gen_b0);
		beta_0_t = *die_b0++;
	
		
		// beta_1
		
		double temp_1 = 0;
		double temp_2 = 0;
		for(int j = 0; j < NUM_SUB; j++)
		{
			double inter = 0;
			
			if(z[j] > 0)
			{
				inter = beta_2_t*gamma_t[j];
			}
			
			if(x[j] > 0)
			{
				int pos_sum = 0;
				int pos_size = pos[j].size();
				for(int k = 0; k < pos_size; k++)
				{
					pos_sum += alpha_t[pos[j][k]];
				}
				temp_1 += (w_t[j] - beta_0_t - inter - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j])*pos_sum;
				// temp_1 += (w_t[j] - beta_0_t - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j])*pos_sum;
				temp_2 += pos_sum*pos_sum;
			}
		}

		double post_mu = (1/(temp_2+tau_0))*temp_1;
		double post_sd = sqrt(1/(temp_2+tau_0));

		if(post_mu > 0)
		{
			gen_type2 die_gen_b1(generator, distribution_type2(post_mu, post_sd));
			boost::generator_iterator<gen_type2> die_b1(&die_gen_b1);
			beta_1_t = *die_b1++;
			while(beta_1_t < 0)
			{
				beta_1_t = *die_b1++;
			}
		}
		else
		{
			/*
			double mu_minus =(-1)*(post_mu/post_sd);
			double alpha_st = (mu_minus + sqrt(mu_minus*mu_minus+4))/2;
			double u1 = 0;
			double z1 = 0;
			double rou_z1 = 1;
			gen_type_ex die_gen_z1(generator, distribution_type_ex(alpha_st));
			boost::generator_iterator<gen_type_ex> die_z1(&die_gen_z1);
			gen_type4 die_gen_u1(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
			while(u1 < rou_z1)
			{
				
				z1 = *die_z1++;
				z1 += mu_minus;
				rou_z1 = exp((-0.5)*(alpha_st-z1)*(alpha_st-z1));
								
				u1 = *die_u1++;
			}
		
			beta_1_t = z1*post_sd+post_mu;
			*/
			
			gen_type4 die_gen_u1(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
			gen_type2 die_gen_z(generator, distribution_type2(0, post_sd));
			boost::generator_iterator<gen_type2> die_z(&die_gen_z);
			
			double u1 = 2;
			double z = 0;
			while(u1>exp((-1)*z*post_mu*(temp_2+tau_0)))
			{	
				z = *die_z++;
				while(z < 0)
				{
					z = *die_z++;
				}
					
				u1 = *die_u1++;
			}
			beta_1_t = z;
			
		}
		
		
		// beta_2

		temp_1 = 0;
		temp_2 = 0;
		int z_num = z_num_mz+z_num_dz;
		if(z_num > 0)
		{
			
			for(int j = 0; j < NUM_SUB; j++)
			{
				if(z[j] > 0)
				{
				
					double main = 0;
			
					int pos_sum = 0;
					int pos_size = pos[j].size();
					for(int k = 0; k < pos_size; k++)
					{
						pos_sum += alpha_t[pos[j][k]];
					}
					main = beta_1_t*pos_sum;
			
					temp_1 += (w_t[j] - beta_0_t - main - mu_cov_t[j] - var_c_t[j] - var_a_m_t[j] - var_a_d1_t[j] - var_a_d2_t[j])*gamma_t[j];
					temp_2 += gamma_t[j]*gamma_t[j];
				}
			}

			post_mu = (1/(temp_2+tau_0))*temp_1;
			post_sd = sqrt(1/(temp_2+tau_0));

			//if(post_mu > 0)
			//{
				gen_type2 die_gen_b2(generator, distribution_type2(post_mu, post_sd));
				boost::generator_iterator<gen_type2> die_b2(&die_gen_b2);
				beta_2_t = *die_b2++;
			/*	while(beta_2_t < 0)
				{
					beta_2_t = *die_b2++;
				}
			}
			else
			{
				
				double mu_minus =(-1)*(post_mu/post_sd);
				double alpha_st = (mu_minus + sqrt(mu_minus*mu_minus+4))/2;
				double u2 = 0;
				double z2 = 0;
				double rou_z2 = 1;
				gen_type_ex die_gen_z2(generator, distribution_type_ex(alpha_st));
				boost::generator_iterator<gen_type_ex> die_z2(&die_gen_z2);
				gen_type4 die_gen_u2(generator, distribution_type4());
				boost::generator_iterator<gen_type4> die_u2(&die_gen_u2);
				while(u2 < rou_z2)
				{			
					z2 = *die_z2++;
					z2 += mu_minus;
					rou_z2 = exp((-0.5)*(alpha_st-z2)*(alpha_st-z2));
				
					u2 = *die_u2++;
				}
		
				beta_2_t = z2*post_sd+post_mu;
				
				gen_type4 die_gen_u1(generator, distribution_type4());
				boost::generator_iterator<gen_type4> die_u1(&die_gen_u1);
				gen_type2 die_gen_z(generator, distribution_type2(0, post_sd));
				boost::generator_iterator<gen_type2> die_z(&die_gen_z);
				double z1 = 1;
				double z2 = 1;
				double u1 = 2;
				double z = 0;
				while(u1>exp((-1)*z*post_mu*(temp_2+tau_0)))
				{	
					z = *die_z++;
					while(z < 0)
					{
						z = *die_z++;
					}
					
					u1 = *die_u1++;
				}
				beta_2_t = z;
			}*/
		}
		
		// P	
		
		gen_type3 die_gen_p(generator, distribution_type3(0.1, 1));
		boost::generator_iterator<gen_type3> die_p(&die_gen_p);
		gen_type3 die_gen_p2(generator, distribution_type3(1.1, 1));
		boost::generator_iterator<gen_type3> die_p2(&die_gen_p2);

		for(int j = 0; j < NUM_RV; j++)
		{
			double p_1 = 0;
			double p_2 = 0;
			double p_3 = 0;
			
			p_1 = *die_p++;
			p_2 = *die_p++;		
			p_3 = *die_p2++;

			double sum_p = p_1 + p_2 + p_3;

			if(alpha_t[j]==1)
			{
				p_t[j][0] = p_3/sum_p;
				p_t[j][1] = p_1/sum_p;
				p_t[j][2] = p_2/sum_p;
			}
			else
			{
				if(alpha_t[j]==0)
				{
					p_t[j][1] = p_3/sum_p;
					p_t[j][0] = p_1/sum_p;
					p_t[j][2] = p_2/sum_p;
				}
				else
				{
					p_t[j][2] = p_3/sum_p;
					p_t[j][0] = p_1/sum_p;
					p_t[j][1] = p_2/sum_p;
				}
			}
		}


		// alpha

		for(int j = 0; j < NUM_RV; j++)
		{
			double temp_a_1 = 0;
			double temp_a_2 = 0;
			double temp_a_3 = 0;
			int num_sub = alpha_ind[j].size();
			for(int k = 0; k < num_sub; k++)
			{
				double second_p = 0;
				double inter = 0;
				int ind_num = alpha_ind[j][k];
				if(z[ind_num]==1)
				{
					inter = beta_2_t*gamma_t[ind_num];
					int pos_sum = 0;
					int pos_size = pos[ind_num].size();
					for(int l = 0; l < pos_size; l++)
					{
						pos_sum += alpha_t[pos[ind_num][l]];
					}
					second_p = beta_1_t*(pos_sum-alpha_t[j]);
				}

				double diff_a_2 = w_t[ind_num] - beta_0_t - second_p - inter - mu_cov_t[ind_num] - var_c_t[ind_num] - var_a_m_t[ind_num] - var_a_d1_t[ind_num] - var_a_d2_t[ind_num];
				// double diff_a_2 = w_t[ind_num] - beta_0_t - second_p - mu_cov_t[ind_num] - var_c_t[ind_num] - var_a_m_t[ind_num] - var_a_d1_t[ind_num] - var_a_d2_t[ind_num];
				double diff_a_1 = diff_a_2 - beta_1_t;
				double diff_a_3 = diff_a_2 + beta_1_t;

				temp_a_1 += diff_a_1*diff_a_1;
				temp_a_2 += diff_a_2*diff_a_2;
				temp_a_3 += diff_a_3*diff_a_3;
			}

			temp_a_1 = p_t[j][0]*exp((-0.5)*temp_a_1);
			temp_a_2 = p_t[j][1]*exp((-0.5)*temp_a_2);
			temp_a_3 = p_t[j][2]*exp((-0.5)*temp_a_3);
			double sum_temp_a = temp_a_1+temp_a_2+temp_a_3;

			gen_type4 die_gen_alpha(generator, distribution_type4());
			boost::generator_iterator<gen_type4> die_alpha(&die_gen_alpha);
			double unif_01 = *die_alpha++;
			if(unif_01 < (temp_a_1/sum_temp_a))
			{
				alpha_t[j] = 1;
			}
			else
			{
				if(unif_01 < ((temp_a_1+temp_a_2)/sum_temp_a))
					alpha_t[j] = 0;
				else
				{
					alpha_t[j] = -1;
				}
			}
		}
		
		
		// cov

		// global alpha sum
		Vec_i alpha_sum(NUM_SUB, 0);

		for(int k = 0; k < NUM_SUB; k++)
		{
			int pos_sum = 0;
			int pos_size = pos[k].size();
			for(int l = 0; l < pos_size; l++)
			{
				pos_sum += alpha_t[pos[k][l]];
			}
			alpha_sum[k] = pos_sum;
		}

		for(int j = 0; j < NUM_COV; j++)
		{
			double new_cov = 0;
			double g_sq = 0;
			double sum_g_g = 0;

			for(int k = 0; k < NUM_SUB; k++)
			{
				/*
				int pos_sum = 0;
				int pos_size = pos[k].size();
				for(int l = 0; l < pos_size; l++)
				{
					pos_sum += alpha_t[pos[k][l]];
				}
				alpha_sum[k] = pos_sum; */
				sum_g_g += (w_t[k]-beta_1_t*alpha_sum[k]*x[k]-beta_2_t*gamma_t[k]*z[k]-beta_0_t-mu_cov_t[k]+cov_t[j]*cov_m[k][j]-var_c_t[k]- var_a_m_t[k] - var_a_d1_t[k] - var_a_d2_t[k])*cov_m[k][j];
				// sum_g_g += (w_t[k]-beta_1_t*alpha_sum[k]*x[k]-beta_0_t-mu_cov_t[k]+cov_t[j]*cov_m[k][j]-var_c_t[k]- var_a_m_t[k] - var_a_d1_t[k] - var_a_d2_t[k])*cov_m[k][j];
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
		

		// w
		for(int j = 0; j < NUM_SUB; j++)
		{
			double eta = beta_1_t*alpha_sum[j]*x[j]+beta_2_t*gamma_t[j]*z[j]+beta_0_t+var_c_t[j]+ var_a_m_t[j] + var_a_d1_t[j]+ var_a_d2_t[j]+mu_cov_t[j];
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
					double z1 = 1;
					double z2 = 1;
					double u1 = 2;
					double z = 0;
					while(u1>exp((-0.5)*(z2-z1)))
					{	
						z = *die_z++;
						while(z < 0)
						{
							z = *die_z++;
						}
						z1 = z*z;
						z2 = (z-eta)*(z-eta)-eta*eta;
					
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
		

		if(model == 3)
		{
			// var_c

			for(int j = 0; j < NUM_FAM; j++)
			{
				double temp = 0;
				int size = fam[j].size();
				for(int k = 1; k < size; k++)
				{
					int sub = fam[j][k];
					double main = 0;
					double inter = 0;
					if(x[sub] > 0)
					{
						main = beta_1_t*alpha_sum[sub];
					}
					
					if(z[sub] > 0)
					{
						inter = beta_2_t*gamma_t[sub];
					}
					
					temp += w_t[sub] - beta_0_t - main - inter - mu_cov_t[sub] - var_a_m_t[sub] - var_a_d1_t[sub] - var_a_d2_t[sub];
					
					// temp += w_t[sub] - beta_0_t - main - mu_cov_t[sub] - var_a_m_t[sub] - var_a_d1_t[sub] - var_a_d2_t[sub];
				}

				gen_type2 die_gen_var_c(generator, distribution_type2((1/(size-1+tau_c_t))*temp, sqrt(1/((size-1)+tau_c_t))));
				boost::generator_iterator<gen_type2> die_var_c(&die_gen_var_c);
				var_c_f_t[j] = *die_var_c++;
				for(int k = 1; k < size; k++)
				{
					int sub = fam[j][k];
					var_c_t[sub] = var_c_f_t[j];
				}
			}

			// tau_c

			double temp_c = 0;
		
			for(int j = 0; j < NUM_FAM; j++)
			{
				temp_c += var_c_f_t[j]*var_c_f_t[j];
			}
		
			gen_type3 die_gen_tau_c(generator, distribution_type3(gamma_a_c_a+0.5*NUM_FAM, 1/(gamma_a_c_b+0.5*temp_c) ));
			boost::generator_iterator<gen_type3> die_tau_c(&die_gen_tau_c);
			tau_c_t = *die_tau_c++;

		}

		if(model > 1)
		{
			double temp_a = 0;

			// var_a_m

			for(int j = 0; j < NUM_MZ; j++)
			{
				double temp = 0;
				int size = mz[j].size();
				for(int k = 1; k < size; k++)
				{
					int sub = mz[j][k];
					double main = 0;
					double inter = 0;
					if(x[sub] > 0)
					{
						main = beta_1_t*alpha_sum[sub];
					}
					
					if(z[sub] > 0)
					{
						inter = beta_2_t*gamma_t[sub];
					}

					temp += w_t[sub] - beta_0_t - main - inter - mu_cov_t[sub] - var_c_t[sub];
					
					// temp += w_t[sub] - beta_0_t - main - mu_cov_t[sub] - var_c_t[sub];
				}

				gen_type2 die_gen_var_a_m(generator, distribution_type2((1/(size-1+tau_a_t))*temp, sqrt(1/((size-1)+tau_a_t))));
				boost::generator_iterator<gen_type2> die_var_a_m(&die_gen_var_a_m);
				var_a_m_f_t[j] = *die_var_a_m++;
				temp_a += var_a_m_f_t[j]*var_a_m_f_t[j];
				for(int k = 1; k < size; k++)
				{
					int sub = mz[j][k];
					var_a_m_t[sub] = var_a_m_f_t[j];
				}
			}

			// var_a_d

			for(int j = 0; j < NUM_DZ; j++)
			{
				double temp = 0;
				int size = dz[j].size();
				Vec mu_temp(size,0);
				for(int k = 1; k < size; k++)
				{
					int sub = dz[j][k];
					double main = 0;
					double inter = 0;
					if(x[sub] > 0)
					{
						main = beta_1_t*alpha_sum[sub];
					}
					
					if(z[sub] > 0)
					{
						inter = beta_2_t*gamma_t[sub];
					}

					mu_temp[k] = beta_0_t + main + inter;
					
					// mu_temp[k] = beta_0_t + main;

					temp += w_t[sub] - mu_temp[k] - mu_cov_t[sub] - var_c_t[sub] - var_a_d2_t[sub];
				}

				gen_type2 die_gen_var_a_d1(generator, distribution_type2((1/(size-1+2*tau_a_t))*temp, sqrt(1/((size-1)+2*tau_a_t))));
				boost::generator_iterator<gen_type2> die_var_a_d1(&die_gen_var_a_d1);
				var_a_d_f_t[j] = *die_var_a_d1++;
				temp_a += 2*var_a_d_f_t[j]*var_a_d_f_t[j];
				for(int k = 1; k < size; k++)
				{
					int sub = dz[j][k];
					var_a_d1_t[sub] = var_a_d_f_t[j];
				}
				/*
				for(int k = 1; k < size; k++)
				{
					int sub = dz[j][k];					
					temp = w_t[sub] - mu_temp[k] - mu_cov_t[sub] - var_c_t[sub] - var_a_d1_t[sub];

					gen_type2 die_gen_var_a_d2(generator, distribution_type2((1/(1+2*tau_d))*temp, sqrt(1/(1+2*tau_d))));
					boost::generator_iterator<gen_type2> die_var_a_d2(&die_gen_var_a_d2);
					var_a_d2_t[sub] = *die_var_a_d2++;
					// temp_a += 2*var_a_d2_t[sub]*var_a_d2_t[sub];
				}
				*/
			}

			// tau_a

			// gen_type3 die_gen_tau_a(generator, distribution_type3(gamma_a_c_a+0.5*(NUM_MZ+NUM_DZ+IND_DZ), 1/(gamma_a_c_b+0.5*temp_a) ));
			gen_type3 die_gen_tau_a(generator, distribution_type3(gamma_a_c_a+0.5*(NUM_MZ+NUM_DZ), 1/(gamma_a_c_b+0.5*temp_a) ));
			boost::generator_iterator<gen_type3> die_tau_a(&die_gen_tau_a);
			tau_a_t = *die_tau_a++;
		}

		beta_0[i] = beta_0_t;
		beta_1[i] = beta_1_t;
		beta_2[i] = beta_2_t;
		alpha.push_back(alpha_t);
		//tau[i] = tau_t;
		if(model > 1)
		{
			tau_a[i] = tau_a_t;
			//tau_d_est += tau_a_t;
			if(model > 2)
			{
				tau_c[i] = tau_c_t;
			}
		}
		if(NUM_COV>0)
			cov_z.push_back(cov_t);

		//cout<<i<<endl;
	}

	// calibration
	/*
	tau_d_est /= ITER_NUM;
	tau_d_est *= 2;
	if(model>1)
	{
		relative = abs(tau_d-tau_d_est)/tau_d;
		calib++;
		tau_d = tau_d_est;
	}
	else
	{
		relative = 0;
	}
	}
	*/
	
	double beta_0_re = 0;
	double beta_1_re = 0;
	double beta_2_re = 0;
	//double tau_re = 0;
	double tau_c_re = 0;
	double tau_a_re = 0;
	int beta_1_posi = 0;
	int beta_2_posi = 0;
	int zero_alpha = 0;
	int zero_gamma = 0;
	
	Vec alpha_re(NUM_RV);
	Vec cov_re(NUM_COV);

	for(int i = 0; i < NUM_RV; i++)
	{
		alpha_re[i] = 0;
	}
	for(int i = 0; i < NUM_COV; i++)
	{
		cov_re[i] = 0;
	}

	int l_t_z = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_re += beta_0[i];
		beta_1_re += beta_1[i];
		beta_2_re += beta_2[i];
		//tau_re += tau[i];
		tau_c_re += tau_c[i];
		tau_a_re += tau_a[i];

		zero_alpha += sum_alpha[i];
		//zero_gamma += sum_gamma[i];
		
		if((sum_alpha[i] == 1)&&(beta_1[i]>thres))
			beta_1_posi++;
		//if((sum_gamma[i-1] == 1)&&(beta_2[i]>thres))
		//	beta_2_posi++;
		if(fabs(beta_2[i])>thres)
			beta_2_posi++;

		for(int j = 0; j < NUM_RV; j++)
		{
			alpha_re[j] += alpha[i][j];
		}
		for(int j = 0; j < NUM_COV; j++)
		{
			cov_re[j] += cov_z[i][j];
		}
	}

	double iter_count = ITER_NUM - burnin;
	beta_0_re /= iter_count;
	beta_1_re /= iter_count;
	beta_2_re /= iter_count;
	//tau_re /= iter_count;
	tau_c_re /= iter_count;
	tau_a_re /= iter_count;

	double zero_alpha_re = 1 - double(zero_alpha)/iter_count;
	// double zero_gamma_re = 1 - double(zero_gamma)/iter_count;
	
	double beta_0_sd = 0;
	double beta_1_sd = 0;
	double beta_2_sd = 0;
	for(int i = burnin; i < ITER_NUM; i++)
	{
		beta_0_sd += (beta_0[i] - beta_0_re)*(beta_0[i] - beta_0_re);
		beta_1_sd += (beta_1[i] - beta_1_re)*(beta_1[i] - beta_1_re);
		beta_2_sd += (beta_2[i] - beta_2_re)*(beta_2[i] - beta_2_re);
	}
	beta_0_sd /= iter_count;
	beta_1_sd /= iter_count;
	beta_2_sd /= iter_count;

	beta_0_sd = sqrt(beta_0_sd);
	beta_1_sd = sqrt(beta_1_sd);
	beta_2_sd = sqrt(beta_2_sd);

	result[0] = beta_0_re;
	result[1] = beta_1_re;
	result[2] = beta_2_re;
	result[3] = zero_alpha_re;
	result[4] = beta_1_posi;
	result[5] = beta_2_posi;
	result[6] = z_num_mz+z_num_dz;
	result[7] = check_mz;
	result[8] = 1/tau_a_re;
	result[9] = 1/tau_c_re;
	result[10] = beta_0_sd;
	result[11] = beta_1_sd;
	result[12] = beta_2_sd;
	
	for(int j = 0; j < NUM_RV; j++)
	{
		alpha_re[j] /= iter_count;
		result[13+j] = alpha_re[j];
	}
	for(int j = 0; j < NUM_COV; j++)
	{
		cov_re[j] /= iter_count;
		result[13+NUM_RV+j] = cov_re[j];
	}

	//beta_1_po /= (1-0.01128341555)/0.01128341555;
	/*
	cout<<beta_0_re<<endl;
	cout<<beta_1_re<<endl;
	cout<<beta_2_re<<endl;
	//cout<<1/tau_re<<endl;
	cout<<beta_1_posi<<endl;
	cout<<beta_2_posi<<endl;
	cout<<z_num_mz+z_num_dz<<endl;
	cout<<check_mz<<endl;
	cout<<1/tau_a_re<<endl;
	cout<<1/tau_c_re<<endl;
	//cout<<calib<<endl;
	//cout<<1/tau_d<<endl;
	
    for(int i = 0; i < NUM_RV; i++)
    {
        //cout<<gamma_re[i]<<endl;
        //cout<<gamma_sd[i]<<endl;
		//cout<<gamma_posi[i]<<endl;
    }*/
	/*
	for(int i = 0; i < NUM_COV; i++)
    {
        result[10+i] = cov_re[i];
		//cout<<cov_re[i]<<endl;
        //cout<<gamma_sd[i]<<endl;
		//cout<<gamma_posi[i]<<endl;
    }
	*/

}