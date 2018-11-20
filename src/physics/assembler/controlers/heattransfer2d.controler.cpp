
#include "heattransfer2d.controler.h"

#include "../kernels/heattransfer2d.kernel.h"

#include "../../../basis/containers/point.h"
#include "../../../basis/containers/serializededata.h"
#include "../../../basis/evaluator/evaluator.h"
#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/physics/heattransfer.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/elementsregionstore.h"
#include "../../../mesh/store/nodestore.h"

using namespace espreso;


HeatTransfer2DControler::HeatTransfer2DControler(Mesh &mesh, const Step &step, const HeatTransferGlobalSettings &gSettings, const HeatTransferStepSettings &sSettings)
: HeatTransferControler(mesh, step, gSettings, sSettings)
{
	_kernel = new HeatTransfer2DKernel(_globalSettings);

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 273.15;
	double defaultThickness = 1;

	_temperature.data = new serializededata<eslocal, double>(1, _nDistribution);
	_coordinates.data = new serializededata<eslocal, double>(2, _nDistribution);
	_gradient.data    = new serializededata<eslocal, double>(2, _mesh.elements->distribution);

	_motion.isConts = tryElementConstness(sSettings.translation_motions, defaultMotion) && defaultMotion.x == defaultMotion.y;
	_motion.data = new serializededata<eslocal, double>(2, _nDistribution, 0);

	_heat.isConts = tryElementConstness(sSettings.heat_source, defaultHeat);
	_heat.data = new serializededata<eslocal, double>(1, _nDistribution, defaultHeat);

	_thickness.isConts = tryElementConstness(gSettings.thickness, defaultThickness);
	_thickness.data = new serializededata<eslocal, double>(1, _nDistribution, defaultThickness);
}

HeatTransfer2DControler::~HeatTransfer2DControler()
{
	delete _kernel;
}

void HeatTransfer2DControler::initData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _coordinates.data->begin(t);
		for (auto n = _mesh.elements->procNodes->datatarray().begin(t); n != _mesh.elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
			c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
		}

		auto evaluate1D = [&] (const std::map<std::string, ECFExpression> &settings, tarray<double> &data) {
			for (auto it = settings.begin(); it != settings.end(); ++it) {
				ElementsRegionStore *region = _mesh.eregion(it->first);
				it->second.evaluator->evaluate(
						region->elements->datatarray().size(t),
						region->elements->datatarray().begin(t),
						_mesh.elements->procNodes->boundarytarray().begin(t),
						2, _coordinates.data->datatarray().data(), NULL, 0, data.data()
				);
			}
		};
		auto evaluateND = [&] (const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data) {
			for (auto it = settings.begin(); it != settings.end(); ++it) {
				ElementsRegionStore *region = _mesh.eregion(it->first);
				it->second.x.evaluator->evaluate(
						region->elements->datatarray().size(t), 2,
						region->elements->datatarray().begin(t),
						_mesh.elements->procNodes->boundarytarray().begin(t),
						2, _coordinates.data->datatarray().data(), NULL, 0, data.data()
				);
				it->second.y.evaluator->evaluate(
						region->elements->datatarray().size(t), 2,
						region->elements->datatarray().begin(t),
						_mesh.elements->procNodes->boundarytarray().begin(t),
						2, _coordinates.data->datatarray().data(), NULL, 0, data.data() + 1
				);
			}
		};

		evaluate1D(_globalSettings.thickness, _thickness.data->datatarray());
		evaluate1D(_stepSettings.heat_source, _heat.data->datatarray());
		evaluateND(_stepSettings.translation_motions, _motion.data->datatarray());
	}
}

void HeatTransfer2DControler::processElements(Matrices matrices, InstanceFiller &filler)
{
	auto enodes = _mesh.elements->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::Iterator iterator;

	size_t noffset = enodes->begin() - _mesh.elements->procNodes->datatarray().begin();
	auto temperature = _temperature.data->begin() + noffset;
	auto coordinates = _coordinates.data->begin() + noffset;
	auto motion      = _motion.data->begin() + noffset;
	auto heat        = _heat.data->begin() + noffset;
	auto thickness   = _thickness.data->begin() + noffset;

	for (eslocal e = filler.begin; e < filler.end; ++e, ++enodes) {
		iterator.element     = _mesh.elements->epointers->datatarray()[e];
		iterator.material    = _mesh.materials[_mesh.elements->material->datatarray()[e]];
		iterator.temperature = temperature->begin(); temperature += enodes->size();
		iterator.coordinates = coordinates->begin(); coordinates += enodes->size();
		iterator.motion      = motion->begin();      motion      += enodes->size();
		iterator.heat        = heat->begin();        heat        += enodes->size();
		iterator.thickness   = thickness->begin();   thickness   += enodes->size();

		_kernel->processElement(matrices, iterator, _step, filler.Ke, filler.Me, filler.Re, filler.fe);
		filler.insert(enodes->size());
	}
}

void HeatTransfer2DControler::processBoundary(Matrices matrices, InstanceFiller &filler)
{

}

