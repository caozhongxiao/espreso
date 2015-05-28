
#include "tetrahedron10.h"

using namespace permoncube;

size_t Tetrahedron10::subelements = Tetrahedron10Subelements;

size_t Tetrahedron10::subnodes[3] = {
		Tetrahedron10Subnodes,
		Tetrahedron10Subnodes,
		Tetrahedron10Subnodes
};

std::vector<idx_t> Tetrahedron10::_coordinateMapping;

void Tetrahedron10::addElements(mesh::Mesh &mesh, const idx_t indices[])
{
	idx_t tetra[20];
	tetra[0] = _coordinateMapping[indices[2]];
	tetra[1] = _coordinateMapping[indices[6]];
	tetra[2] = _coordinateMapping[indices[0]];
	tetra[4] = _coordinateMapping[indices[20]];

	tetra[8] = _coordinateMapping[indices[4]];
	tetra[9] = _coordinateMapping[indices[3]];
	tetra[11] = _coordinateMapping[indices[1]];
	tetra[16] = _coordinateMapping[indices[11]];
	tetra[17] = _coordinateMapping[indices[13]];
	tetra[18] = _coordinateMapping[indices[10]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[6]];
	tetra[1] = _coordinateMapping[indices[0]];
	tetra[2] = _coordinateMapping[indices[20]];
	tetra[4] = _coordinateMapping[indices[18]];

	tetra[8] = _coordinateMapping[indices[3]];
	tetra[9] = _coordinateMapping[indices[10]];
	tetra[11] = _coordinateMapping[indices[13]];
	tetra[16] = _coordinateMapping[indices[12]];
	tetra[17] = _coordinateMapping[indices[9]];
	tetra[18] = _coordinateMapping[indices[19]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[24]];
	tetra[1] = _coordinateMapping[indices[6]];
	tetra[2] = _coordinateMapping[indices[20]];
	tetra[4] = _coordinateMapping[indices[18]];

	tetra[8] = _coordinateMapping[indices[15]];
	tetra[9] = _coordinateMapping[indices[13]];
	tetra[11] = _coordinateMapping[indices[22]];
	tetra[16] = _coordinateMapping[indices[21]];
	tetra[17] = _coordinateMapping[indices[12]];
	tetra[18] = _coordinateMapping[indices[19]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[6]];
	tetra[1] = _coordinateMapping[indices[26]];
	tetra[2] = _coordinateMapping[indices[24]];
	tetra[4] = _coordinateMapping[indices[20]];

	tetra[8] = _coordinateMapping[indices[16]];
	tetra[9] = _coordinateMapping[indices[25]];
	tetra[11] = _coordinateMapping[indices[15]];
	tetra[16] = _coordinateMapping[indices[13]];
	tetra[17] = _coordinateMapping[indices[23]];
	tetra[18] = _coordinateMapping[indices[22]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[8]];
	tetra[1] = _coordinateMapping[indices[26]];
	tetra[2] = _coordinateMapping[indices[6]];
	tetra[4] = _coordinateMapping[indices[20]];

	tetra[8] = _coordinateMapping[indices[17]];
	tetra[9] = _coordinateMapping[indices[16]];
	tetra[11] = _coordinateMapping[indices[7]];
	tetra[16] = _coordinateMapping[indices[14]];
	tetra[17] = _coordinateMapping[indices[23]];
	tetra[18] = _coordinateMapping[indices[13]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));

	tetra[0] = _coordinateMapping[indices[2]];
	tetra[1] = _coordinateMapping[indices[20]];
	tetra[2] = _coordinateMapping[indices[8]];
	tetra[4] = _coordinateMapping[indices[6]];

	tetra[8] = _coordinateMapping[indices[11]];
	tetra[9] = _coordinateMapping[indices[14]];
	tetra[11] = _coordinateMapping[indices[5]];
	tetra[16] = _coordinateMapping[indices[4]];
	tetra[17] = _coordinateMapping[indices[13]];
	tetra[18] = _coordinateMapping[indices[7]];
	mesh.pushElement(new mesh::Tetrahedron10(tetra));
}

void Tetrahedron10::addCoordinates(mesh::Mesh &mesh, const permoncube::Settings &settings, const size_t cluster[])
{
	mesh::Coordinates &coordinates = mesh.coordinates();

	size_t nodes[3];

	ElementGenerator<Tetrahedron10>::clusterNodesCount(settings, nodes);
	mesh.coordinates().resize(nodes[0] * nodes[1] * nodes[2]);

	idx_t global = 0;
	idx_t local = 0;
	idx_t s[3], e[3];
	double step[3];
	for (int i = 0; i < 3; i++) {
		s[i] = (settings.subdomainsInCluster[i] * (settings.elementsInSubdomain[i] * (1 + subnodes[i]))) * cluster[i];
		e[i] = (settings.subdomainsInCluster[i] * (settings.elementsInSubdomain[i] * (1 + subnodes[i]))) * (cluster[i] + 1);
		std::cout << s[i] << " -- " << e[i] << "\n";
	}
	for (int i = 0; i < 3; i++) {
		step[i] = settings.clusterLength[i] / (nodes[i] - 1);
	}

	ElementGenerator<Tetrahedron10>::globalNodesCount(settings, nodes);
	_coordinateMapping.resize(nodes[0] * nodes[1] * nodes[2]);

	for (idx_t z = 0; z < nodes[2]; z++) {
		for (idx_t y = 0; y < nodes[1]; y++) {
			for (idx_t x = 0; x < nodes[0]; x++) {
				_coordinateMapping[global] = local;
				if (s[2] <= z && z <= e[2] && s[1] <= y && y <= e[1] && s[0] <= x && x <= e[0]) {
					coordinates.add(global, mesh::Point(x * step[0], y * step[1], z * step[2]));
					local++;
				}
				global++;
			}
		}
	}
}

void Tetrahedron10::fixZeroPlanes(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	size_t nodes[3];
	ElementGenerator<Tetrahedron10>::globalNodesCount(settings, nodes);
	idx_t index = 0;
	for (idx_t z = 0; z < nodes[2]; z++) {
		for (idx_t y = 0; y < nodes[1]; y++) {
			for (idx_t x = 0; x < nodes[0]; x++) {
				if (z == 0) {
					dirichlet_z[index] = 0.0;
				}
				if (y == 0) {
					dirichlet_y[index] = 0.0;
				}
				if (x == 0) {
					dirichlet_x[index] = 0.0;
				}
				index++;
			}
		}
	}
}

void Tetrahedron10::fixBottom(
			const permoncube::Settings &settings,
			std::map<int, double> &dirichlet_x,
			std::map<int, double> &dirichlet_y,
			std::map<int, double> &dirichlet_z)
{
	size_t nodes[3];
	ElementGenerator<Tetrahedron10>::globalNodesCount(settings, nodes);
	idx_t index = 0;
	for (idx_t y = 0; y < nodes[1]; y++) {
		for (idx_t x = 0; x < nodes[0]; x++) {
			dirichlet_z[index] = 0.0;
			dirichlet_y[index] = 0.0;
			dirichlet_x[index] = 0.0;
			index++;
		}
	}
}

void Tetrahedron10::clear()
{
	_coordinateMapping.clear();
}


