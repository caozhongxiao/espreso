
#include <fstream>

#include "element.h"

#include "../../../basis/evaluators/evaluator.h"
#include "../structures/region.h"
#include "../structures/coordinates.h"

using namespace espreso;

void OldElement::store(std::ofstream& os, const Coordinates &coordinates, size_t part)
{
	eslocal value = vtkCode(), pSize = params(), p;
	os.write(reinterpret_cast<const char *>(&value), sizeof(eslocal));
	for (size_t n = 0; n < nodes(); n++) {
		eslocal index = coordinates.localIndex(node(n), part);
		os.write(reinterpret_cast<const char *>(&index), sizeof(eslocal));
	}
	os.write(reinterpret_cast<const char *>(&pSize), sizeof(eslocal));
	if (pSize) {
		p = param(OldElement::MATERIAL);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		p = param(OldElement::CONSTANT);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		p = param(OldElement::COORDINATES);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
		p = param(OldElement::BODY);
		os.write(reinterpret_cast<const char *>(&p), sizeof(eslocal));
	}
}

bool OldElement::isFaceSwapped(const OldElement* face) const
{
	for (size_t i = 0; i < face->_regions.size(); i++) {
		return false;
	}
	size_t matches = 0;
	for (size_t f = 0; f < this->faces(); f++) {
		const std::vector<eslocal>& faceNodes = this->faceNodes(f);

		if (face->nodes() != faceNodes.size()) {
			continue;
		}

		matches = 0;
		for (; matches < faceNodes.size(); matches++) {
			if (indices()[faceNodes[matches]] != face->node(matches)) {
				break;
			}
		}
		if (matches == faceNodes.size()) {
			return false;
		}
	}

	return true;
}

void OldElement::rotateOutside(const OldElement* parent, const Coordinates &coordinates, Point &normal) const
{
	Point eMid(0, 0, 0), mid(0, 0, 0);
	for (size_t i = 0; i < parent->coarseNodes(); i++) {
		eMid += coordinates[parent->node(i)];
	}
	eMid /= parent->coarseNodes();

	for (size_t i = 0; i < coarseNodes(); i++) {
		mid += coordinates[node(i)];
	}
	mid /= coarseNodes();

	Point outside = mid - eMid;
	if (outside.x * normal.x + outside.y * normal.y + outside.z * normal.z > 0) {
		normal.flip();
	}
}

bool OldElement::hasProperty(Property property, size_t step) const
{
	for (size_t i = 0; i < _regions.size(); i++) {
		if (step < _regions[i]->settings.size() && _regions[i]->settings[step].count(property)) {
			return true;
		}
	}
	return false;
}

double OldElement::sumProperty(Property property, size_t step, const Point &p, double time, double temperature, double defaultValue) const
{
	double result = 0;
	bool set = false;
	for (size_t i = 0; i < _regions.size(); i++) {
		if (step < _regions[i]->settings.size()) {
			auto it = _regions[i]->settings[step].find(property);
			if (it != _regions[i]->settings[step].end()) {
				for (size_t j = 0; j < it->second.size(); j++) {
					result += it->second[j]->evaluate(p, time, temperature);
				}
				set = true;
			}
		}
	}

	return set ? result : defaultValue;
}

double OldElement::getProperty(Property property, size_t step, const Point &p, double time, double temperature, double defaultValue) const
{
	for (size_t i = 0; i < _regions.size(); i++) {
		if (step < _regions[i]->settings.size()) {
			auto it = _regions[i]->settings[step].find(property);
			if (it == _regions[i]->settings[step].end()) {
				continue;
			}
			return it->second.back()->evaluate(p, time, temperature);
		}
	}

	return defaultValue;
}



