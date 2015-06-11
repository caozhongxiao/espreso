
#include "pyramid5.h"

using namespace permoncube;

eslocal Pyramid5::subelements = Pyramid5Subelements;

eslocal Pyramid5::subnodes[3] = {
		Pyramid5Subnodes,
		Pyramid5Subnodes,
		Pyramid5Subnodes,
};

Pyramid5::Pyramid5(const permoncube::Settings &settings): _settings(settings)
{
	eslocal nodes[3];
	Utils<Pyramid5>::clusterNodesCount(_settings, nodes);

	_projection.reserve(clusterNodesCount(_settings));
	eslocal local = 0;
	for (eslocal z = 0; z < nodes[2]; z++) {
		for (eslocal y = 0; y < nodes[1]; y++) {
			for (eslocal x = 0; x < nodes[0]; x++) {
				_projection.push_back(local);
				if (addPoint(x, y, z)) {
					local++;
				}
			}
		}
	}
}

void Pyramid5::addElements(mesh::Mesh &mesh, const eslocal indices[])
{
	eslocal pyramid[8];
	pyramid[0] = _projection[indices[18]];
	pyramid[1] = _projection[indices[20]];
	pyramid[2] = _projection[indices[2]];
	pyramid[3] = _projection[indices[0]];
	pyramid[4] = _projection[indices[13]];
	pyramid[5] = _projection[indices[13]];
	pyramid[6] = _projection[indices[13]];
	pyramid[7] = _projection[indices[13]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[20]];
	pyramid[1] = _projection[indices[26]];
	pyramid[2] = _projection[indices[8]];
	pyramid[3] = _projection[indices[2]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[26]];
	pyramid[1] = _projection[indices[24]];
	pyramid[2] = _projection[indices[6]];
	pyramid[3] = _projection[indices[8]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[24]];
	pyramid[1] = _projection[indices[18]];
	pyramid[2] = _projection[indices[0]];
	pyramid[3] = _projection[indices[6]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[18]];
	pyramid[1] = _projection[indices[24]];
	pyramid[2] = _projection[indices[26]];
	pyramid[3] = _projection[indices[20]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));

	pyramid[0] = _projection[indices[2]];
	pyramid[1] = _projection[indices[8]];
	pyramid[2] = _projection[indices[6]];
	pyramid[3] = _projection[indices[0]];
	//mesh.pushElement(new mesh::Pyramid5(hexa));
}

eslocal Pyramid5::clusterNodesCount(const permoncube::Settings &settings)
{
	eslocal count = Hexahedron8::clusterNodesCount(settings);
	eslocal cElems = Utils<Hexahedron8>::clusterElementsCount(settings);

	return count + cElems;
}

esglobal Pyramid5::globalNodesCount(const permoncube::Settings &settings)
{
	eslocal count = Hexahedron8::globalNodesCount(settings);
	eslocal cElems = Utils<Hexahedron8>::clusterElementsCount(settings);
	cElems *= settings.clusters[0] * settings.clusters[1] * settings.clusters[2];

	return count + cElems;
}

