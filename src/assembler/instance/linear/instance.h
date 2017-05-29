
#ifndef SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_
#define SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_

#include "../instance.h"

#include "../../../configuration/output.h"
#include "../../constraints/constraints.h"
#include "../../../output/resultstore/vtkxmlascii.h"
#include "../../../solver/generic/LinearSolver.h"

namespace espreso {

class ESPRESOSolver;

template <class TPhysics, class TConfiguration>
struct LinearInstance: public OldInstance
{
	LinearInstance(const TConfiguration &configuration, const OutputConfiguration &output, Mesh &mesh): OldInstance(mesh),
	_output(output),
	_configuration(configuration.physics_solver.load_steps_settings.at(1)->espreso),
	_constrains(configuration.physics_solver.load_steps_settings.at(1)->espreso, mesh),
	_physics(mesh, _constrains, configuration),
	_linearSolver(configuration.physics_solver.load_steps_settings.at(1)->espreso, _physics, _constrains),
	_store(_output, &mesh, "results")
	{
		_timeStatistics.totalTime.startWithBarrier();
	};

	virtual void init();
	virtual void solve(std::vector<std::vector<double> > &solution);
	virtual void finalize();

	virtual ~LinearInstance() {};

	virtual const OldPhysics& physics() const { return _physics; }
	virtual const Constraints& constraints() const { return _constrains; }
	virtual OldPhysics& physics() { return _physics; }
	virtual Constraints& constraints() { return _constrains; }

protected:
	const OutputConfiguration &_output;
	const ESPRESOSolver &_configuration;
	Constraints _constrains;
	TPhysics _physics;
	LinearSolver _linearSolver;
	output::VTKXMLASCII _store;
};

}

#include "instance.hpp"

#endif /* SRC_ASSEMBLER_INSTANCE_LINEAR_INSTANCE_H_ */
