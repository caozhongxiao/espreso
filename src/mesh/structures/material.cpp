
#include "material.h"
#include "coordinates.h"
#include "../settings/evaluator.h"
#include "../../config/configuration.h"

using namespace espreso;

Material::Material(const Coordinates &coordinates): _coordinates(coordinates), _models((size_t)PHYSICS::SIZE, MATERIAL_MODEL::SIZE)
{
	for (size_t i = 0; i < (size_t)MATERIAL_PARAMETER::SIZE; i++) {
		_values.push_back(new ConstEvaluator(0));
	}
}

Material::Material(const Coordinates &coordinates, const Configuration &configuration)
: _coordinates(coordinates)
{
	for (size_t i = 0; i < (size_t)MATERIAL_PARAMETER::SIZE; i++) {
		_values.push_back(new ConstEvaluator(0));
	}
	if (configuration.parameters.find("MODEL") != configuration.parameters.end()) {
		_models.resize((size_t)PHYSICS::SIZE, (MATERIAL_MODEL)configuration.parameters.find("MODEL")->second->index());
	} else {
		_models.resize((size_t)PHYSICS::SIZE, MATERIAL_MODEL::SIZE);
	}

	for (size_t p = 0; p < configuration.orderedParameters.size(); p++) {
		if (!StringCompare::caseInsensitiveEq(configuration.orderedParameters[p]->name, "MODEL")) {
			const std::string &value = configuration.orderedParameters[p]->get();
			if (StringCompare::contains(value, "xyzt")) {
				_values[configuration.orderedParameters[p]->index()] = new CoordinatesEvaluator(value, _coordinates);
			} else {
				espreso::Expression expr(value, {});
				_values[configuration.orderedParameters[p]->index()] = new ConstEvaluator(expr.evaluate({}));
			}
		}
	}
}

Material::~Material()
{
	for (size_t i = 0; i < _values.size(); i++) {
		delete _values[i];
	}
}

void Material::set(MATERIAL_PARAMETER parameter, const std::string &value)
{
	delete _values[(size_t)parameter];

	if (StringCompare::contains(value, "xyzt")) {
		_values[(size_t)parameter] = new CoordinatesEvaluator(value, _coordinates);
	} else {
		espreso::Expression expr(value, {});
		_values[(size_t)parameter] = new ConstEvaluator(expr.evaluate({}));
	}
}

void Material::set(MATERIAL_PARAMETER parameter, Evaluator* value)
{
	delete _values[(size_t)parameter];
	_values[(size_t)parameter] = value;
}



