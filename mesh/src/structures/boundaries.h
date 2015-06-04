#ifndef BOUNDARIES_H_
#define BOUNDARIES_H_

#include "../elements/elements.h"
#include "mesh.h"

namespace mesh {

class Boundaries
{

public:

	friend std::ostream& operator<<(std::ostream& os, const Boundaries &f);

	Boundaries(const Mesh &mesh);
	~Boundaries() { };

	template<typename T>
	void create_B1_l(	std::vector < SparseIJVMatrix >      & B1_local, 
						std::vector < SparseIJVMatrix >      & B0_local,
						std::vector < std::vector <T> >      & l2g_vec,
						std::vector < std::vector <eslocal> >    & lambda_map_sub_clst,
						std::vector < std::vector <eslocal> >    & lambda_map_sub_B1,
						std::vector < std::vector <eslocal> >    & lambda_map_sub_B0,
						std::vector < std::vector <double> > & B1_l_duplicity,
						std::map < T, double >             & dirichlet_x,
						std::map < T, double >             & dirichlet_y,
						std::map < T, double >             & dirichlet_z,
						const eslocal domains_num) ;


private:
	/** @brief Reference to a mesh. */
	const Mesh &_mesh;

	/** @brief Keeps mapping of nodes to mesh parts. */
	std::vector<std::set<eslocal> > _boundaries;
};

}

#include "boundaries.hpp"

#endif /* BOUNDARIES_H_ */
