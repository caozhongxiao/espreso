
#ifndef SRC_MESH_STRUCTURES_REGION_H_
#define SRC_MESH_STRUCTURES_REGION_H_

#include <vector>
#include <map>

namespace espreso {

class OldElement;
class OldEvaluator;
class Coordinates;
enum class Property;
enum class ElementType;

struct Region {
	std::string name;
	ElementType eType;
	std::vector<std::map<Property, std::vector<OldEvaluator*> > > settings;
	mutable double area;

	std::vector<OldElement*>& elements() { return *_elements; }
	const std::vector<OldElement*>& elements() const { return *_elements; }

	Region(ElementType eType);
	Region(ElementType eType, std::vector<OldElement*> &element);

	~Region();

	void computeArea(const Coordinates &coordinates) const;

protected:
	std::vector<OldElement*> *_elements;
	bool _destroy;
};

}



#endif /* SRC_MESH_STRUCTURES_REGION_H_ */
