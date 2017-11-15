
#include "solution.h"

#include "../old/mesh/settings/property.h"
#include "../old/mesh/structures/mesh.h"
#include "../old/mesh/structures/coordinates.h"
#include "../old/mesh/structures/elementtypes.h"
#include "../basis/logging/logging.h"

using namespace espreso;

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, size_t DOFs, std::vector<std::vector<double> > &data)
: name(name), eType(eType), DOFs(DOFs), data(data), _statistic(eType, mesh, data, DOFs)
{

}

Solution::Solution(const Mesh &mesh, const std::string &name, ElementType eType, size_t DOFs)
: name(name), eType(eType), DOFs(DOFs), data(_data), _statistic(eType, mesh, data, DOFs)
{
	// TODO: MESH
//	_data.resize(mesh.parts());
//	for (size_t p = 0; p < mesh.parts(); p++) {
//		switch (eType) {
//		case ElementType::ELEMENTS:
//			_data[p].resize(properties.size() * (mesh.getPartition()[p + 1] - mesh.getPartition()[p]));
//			break;
//		case ElementType::NODES:
//			_data[p].resize(properties.size() * mesh.coordinates().localSize(p));
//			break;
//		default:
//			ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: invalid solution element type.";
//		}
//	}
}

void Solution::fill(double value)
{
	for (size_t p = 0; p < _data.size(); p++) {
		std::fill(_data[p].begin(), _data[p].end(), value);
	}
}

void Solution::computeStatisticalData(const Step &step)
{
	_statistic.compute(step);
}

double Solution::getStatisticalData(size_t DOF, StatisticalData data, const Region *region) const
{
	if (DOF > DOFs) {
		return _statistic.getMagnitude(region, data);
	} else {
		return _statistic.get(region, DOF, data);
	}

}


