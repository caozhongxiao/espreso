
#include "heattransfer2d.controler.h"

#include "../kernels/heattransfer2d.kernel.h"

#include "../../instance.h"

#include "../../../basis/containers/serializededata.h"
#include "../../../config/ecf/environment.h"
#include "../../../config/ecf/physics/heattransfer.h"

#include "../../../mesh/mesh.h"
#include "../../../mesh/store/elementstore.h"
#include "../../../mesh/store/nodestore.h"

using namespace espreso;


HeatTransfer2DControler::HeatTransfer2DControler(
		Mesh &mesh, const Step &step,
		const HeatTransferGlobalSettings &gSettings,
		const HeatTransferStepSettings &sSettings,
		const HeatTransferOutputSettings &oSettings)
: HeatTransferControler(mesh, step, gSettings, sSettings, oSettings)
{
	_kernel = new HeatTransfer2DKernel(_globalSettings, _outputSettings);

	Point defaultMotion(0, 0, 0);
	double defaultHeat = 273.15;
	double defaultThickness = 1;

	_ntemperature.data = new serializededata<eslocal, double>(1, _nDistribution);
	_ncoordinate.data = new serializededata<eslocal, double>(2, _nDistribution);

	_nmotion.isConts = tryElementConstness(sSettings.translation_motions, defaultMotion) && defaultMotion.x == defaultMotion.y;
	_nmotion.data = new serializededata<eslocal, double>(2, _nDistribution, defaultMotion.x);

	_nheat.isConts = tryElementConstness(sSettings.heat_source, defaultHeat);
	_nheat.data = new serializededata<eslocal, double>(1, _nDistribution, defaultHeat);

	_nthickness.isConts = tryElementConstness(gSettings.thickness, defaultThickness);
	_nthickness.data = new serializededata<eslocal, double>(1, _nDistribution, defaultThickness);

//	_temperature = _mesh.nodes->appendData(1, { "TEMPERATURE" });
	if (_mesh.hasPhaseChange()) {
		_phaseChange = _mesh.nodes->appendData(1, { "PHASE" });
		_latentHeat = _mesh.nodes->appendData(1, { "LATENT_HEAT" });
	}

	if (_outputSettings.translation_motions && _stepSettings.translation_motions.size()) {
		_motion = _mesh.elements->appendData(2, { "TRANSLATION_MOTION", "TRANSLATION_MOTION_X", "TRANSLATION_MOTION_Y" });
	}

	if (_outputSettings.gradient || _globalSettings.diffusion_split) {
		_gradient = _mesh.elements->appendData(2, { "GRADIENT", "GRADIENT_X", "GRADIENT_Y" });
	}

	if (_outputSettings.flux) {
		_flux = _mesh.elements->appendData(2, { "FLUX", "FLUX_X", "FLUX_Y" });
	}
}

HeatTransfer2DControler::~HeatTransfer2DControler()
{
	delete _kernel;
}

void HeatTransfer2DControler::evaluate(const std::map<std::string, ECFExpression> &settings, tarray<double> &data)
{
	Controler::evaluate(settings, data, 2, _ncoordinate.data->datatarray().data());
}

void HeatTransfer2DControler::evaluate(const std::map<std::string, ECFExpressionVector> &settings, tarray<double> &data)
{
	Controler::evaluate(settings, data, 2, _ncoordinate.data->datatarray().data());
}

void HeatTransfer2DControler::initData()
{
	size_t threads = environment->OMP_NUM_THREADS;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		auto c = _ncoordinate.data->begin(t);
		for (auto n = _mesh.elements->procNodes->datatarray().begin(t); n != _mesh.elements->procNodes->datatarray().end(t); ++n, ++c) {
			c->at(0) = _mesh.nodes->coordinates->datatarray()[*n].x;
			c->at(1) = _mesh.nodes->coordinates->datatarray()[*n].y;
		}
	}

	evaluate(_globalSettings.initial_temperature, _ntemperature.data->datatarray());
	evaluate(_globalSettings.thickness, _nthickness.data->datatarray());
	evaluate(_stepSettings.heat_source, _nheat.data->datatarray());
	evaluate(_stepSettings.translation_motions, _nmotion.data->datatarray());

	nodeValuesToElements(_nmotion.data->datatarray(), *_motion->data);
}

void HeatTransfer2DControler::updateData()
{

}

void HeatTransfer2DControler::processElements(Matrices matrices, InstanceFiller &filler)
{
	auto enodes = _mesh.elements->procNodes->cbegin() + filler.begin;
	HeatTransfer2DKernel::Iterator iterator;

	size_t noffset = enodes->begin() - _mesh.elements->procNodes->datatarray().begin();
	auto temperature = _ntemperature.data->begin() + noffset;
	auto coordinates = _ncoordinate.data->begin() + noffset;
	auto motion      = _nmotion.data->begin() + noffset;
	auto heat        = _nheat.data->begin() + noffset;
	auto thickness   = _nthickness.data->begin() + noffset;

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

