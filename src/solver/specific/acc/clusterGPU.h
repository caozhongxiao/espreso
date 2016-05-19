/*
 * clusterGPU.h
 *
 *  Created on: Feb 24, 2016
 *      Author: lriha
 */

#ifndef SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_
#define SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_


#include "../cluster.h"
#include <cuda.h>
#include <cuda_runtime.h>
namespace espreso {

class ClusterGPU: public ClusterBase
{

public:
	// Constructor
	ClusterGPU(eslocal cluster_index): ClusterBase(cluster_index) {};
	ClusterGPU(): ClusterBase() {};

	void Create_SC_perDomain( bool USE_FLOAT );
	void GetSchurComplement( bool USE_FLOAT, eslocal i );

	void SetupKsolvers ( );

	void multKplusGlobal_GPU   ( SEQ_VECTOR<SEQ_VECTOR<double> > & x_in );


	void multKplus_HF      (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out);
	void multKplus_HF_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);

	void multKplus_HF_Loop1 (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_CP    ();

	void multKplus_HF_Loop2_SC   (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in, SEQ_VECTOR<SEQ_VECTOR<double> > & y_out);
	void multKplus_HF_Loop2_SPDS (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);
	void multKplus_HF_Loop2_MIX  (SEQ_VECTOR<SEQ_VECTOR<double> > & x_in);

};

}



#endif /* SOLVER_SPECIFIC_ACC_CLUSTERGPU_H_ */