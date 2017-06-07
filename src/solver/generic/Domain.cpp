#include "../generic/Domain.h"

#include <cmath>

// *******************************************************************
// **** DOMAIN CLASS ************************************************

using namespace espreso;

Domain::Domain(const ESPRESOSolver &configuration, Instance *instance_in, eslocal domain_index_in, eslocal USE_HTFETI_in):
		configuration(configuration),
		instance(instance_in),

		K(instance_in->K[domain_index_in]),

		Kplus_R(instance_in->N1[domain_index_in]),
		Kplus_R2(instance_in->N2[domain_index_in]),

		_RegMat(instance_in->RegMat[domain_index_in]),

		f(instance_in->f[domain_index_in]),
		vec_c(instance_in->B1c[domain_index_in]),
		vec_lb(instance_in->LB[domain_index_in])

{
		domain_prim_size = K.cols;
		domain_index = domain_index_in;

//		//Kernel setup
//		Kplus_Rb  = Kplus_R;
//		Kplus_Rb2 = Kplus_R2;
//
//		//Constraints and Dirichlet boundary condition
//		B1 = instance_in->B1[domain_index_in];
//		B1.type = 'G';
//		B1t = B1;
//		B1t.MatTransposeCOO();
//		B1t.ConvertToCSRwithSort(1);
//
//		B1_scale_vec = instance->B1duplicity[domain_index_in];
//		lambda_map_sub = instance->B1subdomainsMap[domain_index_in];
//
//
//		// HTFETI Section
		USE_HFETI = USE_HTFETI_in;
//		if (USE_HFETI == 1) {
//			B0 = instance_in->B0[domain_index_in];
//			B0.type = 'G';
//			B0.ConvertToCSRwithSort(1);
//		}



}

void Domain::SetDomain() {

	//Kernel setup
	Kplus_Rb  = Kplus_R;
	Kplus_Rb2 = Kplus_R2;

	//Constraints and Dirichlet boundary condition
	B1 = instance->B1[domain_index];
	B1.type = 'G';
	B1t = B1;
	B1t.MatTransposeCOO();
	B1t.ConvertToCSRwithSort(1);

	B1_scale_vec = instance->B1duplicity[domain_index];
	lambda_map_sub = instance->B1subdomainsMap[domain_index];


	// HTFETI Section
	//USE_HFETI = USE_HTFETI_in;
	if (USE_HFETI == 1) {
		B0 = instance->B0[domain_index];
		B0.type = 'G';
		B0.ConvertToCSRwithSort(1);
	}


}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out, eslocal x_in_vector_start_index, eslocal y_out_vector_start_index) {
	switch (configuration.Ksolver) {
	case ESPRESO_KSOLVER::DIRECT_DP:
		Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
		break;
	case ESPRESO_KSOLVER::DIRECT_SP:
		Kplus.Solve(x_in, y_out, x_in_vector_start_index, y_out_vector_start_index);
		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Invalid KSOLVER value.";
		exit(EXIT_FAILURE);
	}


	//		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
	//		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
	//		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
	//
	//		Kplus.Solve(x_in, x, x_in_vector_start_index, 0);
	//		if (enable_SP_refinement) {
	//			for (eslocal step = 0; step <= esconfiguration.Ksolver_SP_iter_steps; step++) {
	//				K.MatVec(x,r,'N');
	//				for (eslocal i = 0; i < r.size(); i++)
	//					r[i] = x_in[i + x_in_vector_start_index] - r[i];
	//				Kplus.Solve(r, z, 0, 0);
	//				for (eslocal i = 0; i < r.size(); i++)
	//					x[i] = x[i] + z[i];
	//
	//				double norm = 0.0;
	//				for (eslocal i = 0; i < r.size(); i++)
	//					norm += r[i]*r[i];
	//
	//				norm = sqrt(norm);
	//
	//				if (norm < esconfiguration.Ksolver_SP_iter_norm) {
	//					std::cout.precision(20);
	//					std::cout << "Refinement steps: " << step << " | norm: " << norm << std::endl;
	//					break;
	//				}
	//
	//			}
	//		}
	//
	//		for (eslocal i = 0; i < r.size(); i++)
	//			y_out[y_out_vector_start_index + i] = x[i];

//	case 1: {
//		Kplus.SolveCG(K, x_in_y_out);
//		break;
//	}

}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in, SEQ_VECTOR <double> & y_out) {
	switch (configuration.Ksolver) {
	case ESPRESO_KSOLVER::DIRECT_DP:
		Kplus.Solve(x_in, y_out, 0, 0);
		break;
	case ESPRESO_KSOLVER::ITERATIVE:
		Kplus.SolveCG(K, x_in, y_out);
		break;
	case ESPRESO_KSOLVER::DIRECT_SP:
		Kplus.Solve(x_in, y_out, 0, 0);
		break;
	case ESPRESO_KSOLVER::DIRECT_MP: {

		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);

		bool success = false;

		Kplus.Solve(x_in, x, 0, 0);
		if (enable_SP_refinement) {
			for (size_t step = 0; step <= configuration.Ksolver_iterations; step++) {
				K.MatVec(x,r,'N');
				for (size_t i = 0; i < r.size(); i++) {
					r[i] = x_in[i] - r[i];
				}
				Kplus.Solve(r, z, 0, 0);
				//Kplus.SolveCG(K, r, z);
				for (size_t i = 0; i < r.size(); i++) {
					x[i] = x[i] + z[i];
				}

				double norm = 0.0;
				for (size_t i = 0; i < r.size(); i++) {
					norm += r[i]*r[i];
				}

				norm = sqrt(norm);

				if (norm < configuration.Ksolver_epsilon) {
					ESINFO(PROGRESS3) << " " << step;
					success = true;
					break;
				}
			}
		}

		if (!success) {
			ESINFO(PROGRESS3) << "FAILED";
		}

		for (size_t i = 0; i < r.size(); i++) {
			y_out[i] = x[i];
		}
		break;
	}
//	case 4:
//		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//		Kplus.Solve(x_in, x, 0, 0);
//		Kplus.SolveCG(K, x_in, y_out, x);
//		break;
	default:
		ESINFO(GLOBAL_ERROR) << "Invalid KSOLVER value.";
		exit(EXIT_FAILURE);
	}
}

void Domain::multKplusLocal(SEQ_VECTOR <double> & x_in_y_out) {
	switch (configuration.Ksolver) {
	case ESPRESO_KSOLVER::DIRECT_DP:
		Kplus.Solve(x_in_y_out);
		break;
	case ESPRESO_KSOLVER::DIRECT_SP:
		Kplus.Solve(x_in_y_out);
		break;
	case ESPRESO_KSOLVER::DIRECT_MP: {
		bool success = false;

		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> r (Kplus.m_Kplus_size, 0.0);
		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);

		Kplus.Solve(x_in_y_out, x, 0, 0);

		if (enable_SP_refinement) {
			for (size_t step = 0; step <= configuration.Ksolver_iterations; step++) {
				K.MatVec(x,r,'N');
				for (size_t i = 0; i < r.size(); i++) {
					r[i] = x_in_y_out[i] - r[i];
				}
				Kplus.Solve(r, z, 0, 0);
				for (size_t i = 0; i < r.size(); i++) {
					x[i] = x[i] + z[i];
				}

				double norm = 0.0;
				for (size_t i = 0; i < r.size(); i++) {
					norm += r[i]*r[i];
				}

				norm = sqrt(norm);

				if (norm < configuration.Ksolver_epsilon) {
					ESINFO(PROGRESS3) << " " << step;
					break;
				}

			}
		}

		if (!success) {
			ESINFO(PROGRESS3) << "FAILED";
		}

		for (size_t i = 0; i < r.size(); i++) {
			x_in_y_out[i] = x[i];
		}

		break;
	}
//	case 4: { // DIRECT MIX - 2xSP
//
//		SEQ_VECTOR<double> x (Kplus.m_Kplus_size, 0.0);
//		SEQ_VECTOR<double> z (Kplus.m_Kplus_size, 0.0);
//
//		Kplus.Solve(x_in_y_out, x, 0, 0);
//		Kplus.SolveCG(K, x_in_y_out, z, x);
//
//		for (eslocal i = 0; i < z.size(); i++)
//			x_in_y_out[i] = z[i];
//
//		break;
//	}
	case ESPRESO_KSOLVER::ITERATIVE:
		Kplus.SolveCG(K, x_in_y_out);
		break;
	default:
		ESINFO(ERROR) << "Invalid KSOLVER value.";
		exit(EXIT_FAILURE);
	}
}

// **** END - DOMAIN CLASS *******************************************
// *******************************************************************
