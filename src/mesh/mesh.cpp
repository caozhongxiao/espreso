
#include "mesh.h"

#include "esinfo/ecfinfo.h"
#include "esinfo/meshinfo.h"
#include "esinfo/envinfo.h"
#include "esinfo/mpiinfo.h"

#include "store/statisticsstore.h"
#include "store/elementstore.h"
#include "store/nodestore.h"
#include "store/elementsregionstore.h"
#include "store/boundaryregionstore.h"
#include "store/surfacestore.h"
#include "store/contactstore.h"

#include "output/result/resultstore.h"

#include "preprocessing/meshpreprocessing.h"

#include "elements/elements.h"

#include "basis/utilities/utils.h"
#include "basis/utilities/communication.h"
#include "basis/utilities/parser.h"

#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

#include "basis/containers/serializededata.h"
#include "basis/containers/tarray.h"

using namespace espreso;


Mesh::Mesh()
: elements(new ElementStore(_eclasses)), nodes(new NodeStore()),
  FETIData(NULL),
  halo(new ElementStore(_eclasses)),
  surface(NULL), domainsSurface(NULL),
  contacts(NULL),
  preprocessing(new MeshPreprocessing(this)),

  _eclasses(info::env::OMP_NUM_THREADS),
  store(NULL),
  _withGUI(false)
{
	size_t threads = info::env::OMP_NUM_THREADS;

	dimension = 0;
	switch (info::ecf->physics) {
	case PHYSICS::HEAT_TRANSFER_2D:
	case PHYSICS::STRUCTURAL_MECHANICS_2D:
	case PHYSICS::SHALLOW_WATER_2D:
		dimension = 2;
		break;
	case PHYSICS::HEAT_TRANSFER_3D:
	case PHYSICS::STRUCTURAL_MECHANICS_3D:
		dimension = 3;
		break;
	}

	preferedDomains = 1;
	uniformDecomposition = true;

	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		std::vector<Element> eclasses;
		_eclasses[t] = new Element[static_cast<int>(Element::CODE::SIZE)];

		eclasses.push_back(Point1::create(_eclasses[t]));

		eclasses.push_back(Line2::create(_eclasses[t]));
		eclasses.push_back(Line3::create(_eclasses[t]));

		eclasses.push_back(Triangle3::create(_eclasses[t]));
		eclasses.push_back(Triangle6::create(_eclasses[t]));
		eclasses.push_back(Square4::create(_eclasses[t]));
		eclasses.push_back(Square8::create(_eclasses[t]));

		eclasses.push_back(Tetrahedron4::create(_eclasses[t]));
		eclasses.push_back(Tetrahedron10::create(_eclasses[t]));
		eclasses.push_back(Pyramid5::create(_eclasses[t]));
		eclasses.push_back(Pyramid13::create(_eclasses[t]));
		eclasses.push_back(Prisma6::create(_eclasses[t]));
		eclasses.push_back(Prisma15::create(_eclasses[t]));
		eclasses.push_back(Hexahedron8::create(_eclasses[t]));
		eclasses.push_back(Hexahedron20::create(_eclasses[t]));

		std::sort(eclasses.begin(), eclasses.end(), [] (const Element &e1, const Element &e2) { return e1.code < e2.code; });

		memcpy(_eclasses[t], eclasses.data(), eclasses.size() * sizeof(Element));
	}
}

ElementsRegionStore* Mesh::eregion(const std::string &name)
{
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (StringCompare::caseSensitiveEq(elementsRegions[r]->name, name)) {
			return elementsRegions[r];
		}
	}
	ESINFO(ERROR) << "Unknown region of elements with '" << name << "'.";
	return NULL;
}

ElementsRegionsIntersectionStore* Mesh::ieregion(const std::string &name)
{
	for (size_t r = 0; r < elementsRegionsIntersections.size(); r++) {
		if (StringCompare::caseSensitiveEq(elementsRegionsIntersections[r]->name, name)) {
			return elementsRegionsIntersections[r];
		}
	}
	ESINFO(ERROR) << "ESPRESO internal error: request for unknown intersection of element regions '" << name << "'.";
	return NULL;
}

BoundaryRegionStore* Mesh::bregion(const std::string &name)
{
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (StringCompare::caseSensitiveEq(boundaryRegions[r]->name, name)) {
			return boundaryRegions[r];
		}
	}
	ESINFO(ERROR) << "Unknown boundary region '" << name << "'.";
	return NULL;
}

BoundaryRegionsIntersectionStore* Mesh::ibregion(const std::string &name)
{
	for (size_t r = 0; r < boundaryRegionsIntersections.size(); r++) {
		if (StringCompare::caseSensitiveEq(boundaryRegionsIntersections[r]->name, name)) {
			return boundaryRegionsIntersections[r];
		}
	}
	ESINFO(ERROR) << "ESPRESO internal error: request for unknown intersection of boundary regions '" << name << "'.";
	return NULL;
}

bool Mesh::onAllElements(const std::string &eregion) const
{
	return StringCompare::caseInsensitiveEq(eregion, "ALL_ELEMENTS");
}

bool Mesh::hasPhaseChange() const
{
	for (size_t m = 0; m < materials.size(); m++) {
		if (materials[m]->phase_change) {
			return true;
		}
	}
	return false;
}

void Mesh::storeMesh()
{
	store->updateMesh();
}

void Mesh::storeSolution()
{
	store->updateSolution();
}

void Mesh::update()
{
	int isEmpty = 0, quit;
	if (elements->size == 0) {
		isEmpty = 1;
	}

	MPI_Allreduce(&isEmpty, &quit, 1, MPI_INT, MPI_MAX, info::mpi::comm);
	if (quit) {
		ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::YELLOW << "ESPRESO quit computation. There is a process with no elements.";
		MPI_Barrier(info::mpi::comm);
		exit(EXIT_SUCCESS);
	}

	auto hasBEM = [] (const PhysicsConfiguration &physics) {
		for (auto it = physics.discretization.begin(); it != physics.discretization.end(); ++it) {
			if (it->second == DISCRETIZATION::BEM) {
				return true;
			}
		}
		return false;
	};

	auto getPhysics = [&] () -> const PhysicsConfiguration& {
		switch (info::ecf->physics) {
		case PHYSICS::HEAT_TRANSFER_2D:
			return info::ecf->heat_transfer_2d;
		case PHYSICS::HEAT_TRANSFER_3D:
			return info::ecf->heat_transfer_3d;
		case PHYSICS::STRUCTURAL_MECHANICS_2D:
			return info::ecf->structural_mechanics_2d;
		case PHYSICS::STRUCTURAL_MECHANICS_3D:
			return info::ecf->structural_mechanics_3d;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented physics.";
			exit(0);
		}
	};

	auto is3D = [&] () {
		switch (info::ecf->physics) {
		case PHYSICS::HEAT_TRANSFER_2D:
			return false;
		case PHYSICS::HEAT_TRANSFER_3D:
			return true;
		case PHYSICS::STRUCTURAL_MECHANICS_2D:
			return false;
		case PHYSICS::STRUCTURAL_MECHANICS_3D:
			return true;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented physics.";
			exit(0);
		}
	};

	auto forEachSteps = [&] (std::function<bool(const LoadStepConfiguration &)> fnc) {
		bool ret = false;
		switch (info::ecf->physics) {
		case PHYSICS::HEAT_TRANSFER_2D:
			for (auto step = info::ecf->heat_transfer_2d.load_steps_settings.begin(); step != info::ecf->heat_transfer_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PHYSICS::HEAT_TRANSFER_3D:
			for (auto step = info::ecf->heat_transfer_3d.load_steps_settings.begin(); step != info::ecf->heat_transfer_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PHYSICS::STRUCTURAL_MECHANICS_2D:
			for (auto step = info::ecf->structural_mechanics_2d.load_steps_settings.begin(); step != info::ecf->structural_mechanics_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PHYSICS::STRUCTURAL_MECHANICS_3D:
			for (auto step = info::ecf->structural_mechanics_3d.load_steps_settings.begin(); step != info::ecf->structural_mechanics_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented physics.";
			exit(0);
		}
		return ret;
	};

	auto ntob = [&] (const std::string &rname, int dimension) {
		BoundaryRegionStore *region = bregion(rname);
		if (region->dimension == 0) {
			preprocessing->computeBoundaryElementsFromNodes(region, dimension);
		}
	};

	preprocessing->startPreprocessing();

	ESINFO(OVERVIEW) << "Preprocess mesh data.";
	materials.clear();
	std::map<std::string, int> matindex;
	for (auto mat = info::ecf->getPhysics()->materials.begin(); mat != info::ecf->getPhysics()->materials.end(); ++mat) {
		materials.push_back(&mat->second);
		matindex[mat->first] = materials.size() - 1;
	}

	for (auto mat = info::ecf->getPhysics()->material_set.begin(); mat != info::ecf->getPhysics()->material_set.end(); ++mat) {
		ElementsRegionStore *region = eregion(mat->first);
		if (matindex.find(mat->second) == matindex.end()) {
			ESINFO(GLOBAL_ERROR) << "Unknown material '" << mat->second << "'.";
		}
		int material = matindex.find(mat->second)->second;
		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
			elements->material->datatarray()[*e] = material;
		}
	}

	if (
			info::ecf->decomposition.separate_materials ||
			info::ecf->decomposition.separate_regions ||
			info::ecf->decomposition.separate_etypes) {
		uniformDecomposition = false;
	}

	if (info::ecf->decomposition.balance_clusters) {
		preprocessing->reclusterize();
	}

	uniformDecomposition = false;
	if (uniformDecomposition) {
		// implement uniform decomposition
	} else {
		preprocessing->partitiate(preferedDomains);
	}

	if (info::ecf->physics == PHYSICS::STRUCTURAL_MECHANICS_2D || info::ecf->physics == PHYSICS::STRUCTURAL_MECHANICS_3D) {
		const StructuralMechanicsConfiguration *sm;
		int dimension;
		if (info::ecf->physics == PHYSICS::STRUCTURAL_MECHANICS_2D) {
			sm = &info::ecf->structural_mechanics_2d;
			dimension = 1;
		}
		if (info::ecf->physics == PHYSICS::STRUCTURAL_MECHANICS_3D) {
			sm = &info::ecf->structural_mechanics_3d;
			dimension = 2;
		}

		for (auto ls = sm->load_steps_settings.begin(); ls != sm->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.normal_pressure.begin(); bc != ls->second.normal_pressure.end(); ++bc) {
				ntob(bc->first, dimension);
			}
		}
	}

	if (info::ecf->physics == PHYSICS::HEAT_TRANSFER_2D || info::ecf->physics == PHYSICS::HEAT_TRANSFER_3D) {
		const HeatTransferConfiguration *ht;
		int dimension;
		if (info::ecf->physics == PHYSICS::HEAT_TRANSFER_2D) {
			ht = &info::ecf->heat_transfer_2d;
			dimension = 1;
		}
		if (info::ecf->physics == PHYSICS::HEAT_TRANSFER_3D) {
			ht = &info::ecf->heat_transfer_3d;
			dimension = 2;
		}

		for (auto ls = ht->load_steps_settings.begin(); ls != ht->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.heat_flow.begin(); bc != ls->second.heat_flow.end(); ++bc) {
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.heat_flux.begin(); bc != ls->second.heat_flux.end(); ++bc) {
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.convection.begin(); bc != ls->second.convection.end(); ++bc) {
				ntob(bc->first, dimension);
			}
			for (auto bc = ls->second.diffuse_radiation.begin(); bc != ls->second.diffuse_radiation.end(); ++bc) {
				ntob(bc->first, dimension);
			}
		}
	}

	preprocessing->arrangeRegions();

	if (forEachSteps([] (const LoadStepConfiguration &step) {
		return step.solver == LoadStepConfiguration::SOLVER::FETI;
	})) {

		preprocessing->computeLocalIndices();
	}

	if (is3D() && (hasBEM(getPhysics()))) {
		preprocessing->computeDomainsSurface();
		preprocessing->triangularizeDomainSurface();
	}

	if (is3D() && _withGUI) {
		preprocessing->computeRegionsSurface();
		for (size_t r = 0; r < elementsRegions.size(); r++) {
			preprocessing->triangularizeSurface(elementsRegions[r]->surface);
		}
		for (size_t r = 0; r < boundaryRegions.size(); r++) {
			preprocessing->triangularizeBoundary(boundaryRegions[r]);
		}
	}

	if (is3D() && info::ecf->output.format == OutputConfiguration::FORMAT::STL_SURFACE) {
		preprocessing->computeBodiesSurface();
		preprocessing->triangularizeSurface(surface);
	}

	if (forEachSteps([] (const LoadStepConfiguration &step) {
		return step.solver == LoadStepConfiguration::SOLVER::FETI && step.feti.method == FETI_METHOD::HYBRID_FETI && step.feti.B0_type == FETI_B0_TYPE::KERNELS;
	})) {

		preprocessing->computeSharedFaceNodes();
	}

	if (forEachSteps([] (const LoadStepConfiguration &step) {
		return step.solver == LoadStepConfiguration::SOLVER::FETI && step.feti.method == FETI_METHOD::HYBRID_FETI && step.feti.B0_type == FETI_B0_TYPE::CORNERS;
	})) {

		preprocessing->computeCornerNodes();
	}

	if (info::ecf->physics == PHYSICS::STRUCTURAL_MECHANICS_2D || info::ecf->physics == PHYSICS::STRUCTURAL_MECHANICS_3D) {
		if (forEachSteps([] (const LoadStepConfiguration &step) {
			return step.solver == LoadStepConfiguration::SOLVER::FETI && step.feti.regularization == FETI_REGULARIZATION::ANALYTIC;
		})) {

			preprocessing->computeFixPoints();
			if (hasBEM(getPhysics())) {
				preprocessing->computeFixPointsOnSurface();
			}
		}
	}

	if (info::ecf->mesh_morphing.type == MORPHING_TYPE::RBF) {
		for (auto it = info::ecf->mesh_morphing.rbf.begin(); it != info::ecf->mesh_morphing.rbf.end(); ++it) {
			switch (info::ecf->physics) {
			case PHYSICS::HEAT_TRANSFER_2D:
				preprocessing->morphRBF(it->first, it->second, 2);
				break;
			case PHYSICS::HEAT_TRANSFER_3D:
				preprocessing->morphRBF(it->first, it->second, 3);
				break;
			case PHYSICS::STRUCTURAL_MECHANICS_2D:
				preprocessing->morphRBF(it->first, it->second, 2);
				break;
			case PHYSICS::STRUCTURAL_MECHANICS_3D:
				preprocessing->morphRBF(it->first, it->second, 3);
				break;
			default:
				ESINFO(GLOBAL_ERROR) << "Not implemented physics.";
				exit(0);
			}
		}
	}

	if (getPhysics().contact_interfaces) {
		preprocessing->computeBodiesSurface();
		contacts = new ContactStore(surface);
		preprocessing->computeSurfaceLocations();
		preprocessing->searchContactInterfaces();
	}

	if (info::ecf->output.verbose_level > 0) {
		printMeshStatistics();
	}
	if (info::ecf->output.verbose_level > 1) {
		printDecompositionStatistics();
	}
}

void Mesh::printMeshStatistics()
{
	size_t namesize = 25;

	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (namesize < boundaryRegions[r]->name.size() + 5) {
			namesize = boundaryRegions[r]->name.size() + 5;
		}
	}

	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (namesize < elementsRegions[r]->name.size() + 5) {
			namesize = elementsRegions[r]->name.size() + 5;
		}
	}

	auto ename = [] (int code) -> std::string {
		switch (static_cast<Element::CODE>(code)) {

		case Element::CODE::POINT1: return "POINT1";

		case Element::CODE::LINE2: return "LINE2";

		case Element::CODE::TRIANGLE3: return "TRIANGLE3";
		case Element::CODE::SQUARE4: return "SQUARE4";

		case Element::CODE::TETRA4: return "TETRA4";
		case Element::CODE::PYRAMID5: return "PYRAMID5";
		case Element::CODE::PRISMA6: return "PRISMA6";
		case Element::CODE::HEXA8: return "HEXA8";

		case Element::CODE::LINE3: return "LINE3";

		case Element::CODE::TRIANGLE6: return "TRIANGLE6";
		case Element::CODE::SQUARE8: return "SQUARE8";

		case Element::CODE::TETRA10: return "TETRA10";
		case Element::CODE::PYRAMID13: return "PYRAMID13";
		case Element::CODE::PRISMA15: return "PRISMA15";
		case Element::CODE::HEXA20: return "HEXA20";

		default:
			ESINFO(ERROR) << "ESPRESO internal error: unknown element code.";
			return "";
		}
	};

	auto esize = [&] (int code) {;
		if (elements->ecounters[code]) {
			ESINFO(OVERVIEW) << std::string(namesize - ename(code).size(), ' ') << ename(code) << " : " << elements->ecounters[code];
		}
	};

	auto totalesize = [] (std::vector<esint> &ecounters) {
		esint size = 0;
		for (size_t etype = 0; etype < ecounters.size(); etype++) {
			size += ecounters[etype];
		}
		return size;
	};

	auto header = [&] (const std::string &name) {
		return name + std::string(namesize - name.size(), ' ') + " : ";
	};

	ESINFO(OVERVIEW) << "============= Mesh statistics ==============";

	ESINFO(OVERVIEW) << header(" Number of nodes") << nodes->uniqueTotalSize;

	ESINFO(OVERVIEW) << header(" Number of elements") << totalesize(elements->ecounters);
	for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
		esize(etype);
	}

	ESINFO(OVERVIEW);

	ESINFO(OVERVIEW) << header(" Element regions size");
	for (size_t r = 0; r < elementsRegions.size(); r++) {
		if (StringCompare::caseInsensitiveEq(elementsRegions[r]->name, "NAMELESS_ELEMENT_SET")) {
			ESINFO(OVERVIEW) << Info::TextColor::YELLOW
					<< std::string(namesize - elementsRegions[r]->name.size(), ' ') << elementsRegions[r]->name << " : " << totalesize(elementsRegions[r]->ecounters);
		} else if (StringCompare::caseInsensitiveEq(elementsRegions[r]->name, "ALL_ELEMENTS")) {
			ESINFO(OVERVIEW) << std::string(namesize - elementsRegions[r]->name.size(), ' ') << elementsRegions[r]->name << " : " << totalesize(elements->ecounters);
		} else {
			ESINFO(OVERVIEW) << std::string(namesize - elementsRegions[r]->name.size(), ' ') << elementsRegions[r]->name << " : " << totalesize(elementsRegions[r]->ecounters);
		}
	}
	ESINFO(OVERVIEW) << header(" Face regions size");
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->dimension == 2) {
			ESINFO(OVERVIEW) << std::string(namesize - boundaryRegions[r]->name.size(), ' ') << boundaryRegions[r]->name << " : " << totalesize(boundaryRegions[r]->ecounters);
		}
	}
	ESINFO(OVERVIEW) << header(" Edge regions size");
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->dimension == 1) {
			ESINFO(OVERVIEW) << std::string(namesize - boundaryRegions[r]->name.size(), ' ') << boundaryRegions[r]->name << " : " << totalesize(boundaryRegions[r]->ecounters);
		}
	}
	ESINFO(OVERVIEW) << header(" Node regions size");
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->dimension == 0) {
			ESINFO(OVERVIEW) << std::string(namesize - boundaryRegions[r]->name.size(), ' ') << boundaryRegions[r]->name << " : " << boundaryRegions[r]->uniqueTotalSize;
		}
	}

	ESINFO(OVERVIEW) << "============================================\n";
}

void Mesh::printDecompositionStatistics()
{
	size_t namesize = 25;

	auto header = [&] (const std::string &name) {
		return name + std::string(namesize - name.size(), ' ') + " : ";
	};

	ESINFO(DETAILS) << "========= Decomposition statistics =========";

	ESINFO(DETAILS) << header(" NUMBER OF PROCESSES") << info::mpi::size;

	ESINFO(DETAILS);

	esint totalNeighbors = 0, minNeighbors = 0, maxNeighbors = 0, neighbors = neighbours.size();
	MPI_Reduce(&neighbors, &minNeighbors, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&neighbors, &totalNeighbors, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&neighbors, &maxNeighbors, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << header(" NUMBER OF NEIGHBORS") << totalNeighbors;
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minNeighbors << ", " << maxNeighbors << " (" << totalNeighbors / (double)info::mpi::size << ")";
	ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxNeighbors / (double)minNeighbors;

	ESINFO(DETAILS);

	esint totalClusters = 0, minClusters = 0, maxClusters = 0, clusters = elements->nclusters;
	MPI_Reduce(&clusters, &minClusters, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&clusters, &totalClusters, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&clusters, &maxClusters, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << header(" NUMBER OF CLUSTERS") << totalClusters;
	ESINFO(DETAILS) << header(" clusters per MPI");
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minClusters << ", " << maxClusters << " (" << totalClusters / (double)info::mpi::size << ")";
	ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxClusters / (double)minClusters;

	ESINFO(DETAILS);

	esint totalDomains = 0, minDomains = 0, maxDomains = 0, domains = elements->ndomains;
	MPI_Reduce(&domains, &minDomains, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&domains, &totalDomains, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&domains, &maxDomains, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << header(" NUMBER OF DOMAINS") << totalDomains;
	ESINFO(DETAILS) << header(" domains per MPI");
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minDomains << ", " << maxDomains << " (" << totalDomains / (double)info::mpi::size << ")";
	if (maxDomains / (double)minDomains > 3) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxDomains / (double)minDomains;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxDomains / (double)minDomains;
	}

	ESINFO(DETAILS) << header(" domains per cluster");

	esint cdomains = elements->ndomains;
	for (esint c = 0; c < elements->nclusters; c++) {
		domains = 0;
		for (esint d = 0; d < elements->ndomains; d++) {
			if (elements->clusters[d] == c) {
				++domains;
			}
		}
		cdomains = std::min(domains, cdomains);
	}
	MPI_Reduce(&domains, &minDomains, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	for (esint c = 0; c < elements->nclusters; c++) {
		domains = 0;
		for (esint d = 0; d < elements->ndomains; d++) {
			if (elements->clusters[d] == c) {
				++domains;
			}
		}
		cdomains = std::max(domains, cdomains);
	}
	MPI_Reduce(&domains, &maxDomains, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);

	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minDomains << ", " << maxDomains << " (" << totalDomains / (double)totalClusters << ")";
	if (maxDomains / (double)minDomains > 3) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxDomains / (double)minDomains;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxDomains / (double)minDomains;
	}

	ESINFO(DETAILS);

	esint totalElements = 0, minElements = 0, maxElements = 0;
	esint minelements = elements->elementsDistribution[1], maxelements = 0, avgelements = elements->size;
	MPI_Reduce(&avgelements, &minElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&avgelements, &totalElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&avgelements, &maxElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << header(" NUMBER OF ELEMENTS") << totalElements;
	ESINFO(DETAILS) << header(" elements per MPI");
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minElements << ", " << maxElements << " (" << totalElements / (double)info::mpi::size << ")";
	if (maxElements / (double)minElements > 1.5) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxElements / (double)minElements;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxElements / (double)minElements;
	}


	ESINFO(DETAILS) << header(" elements per cluster");

	esint coffset = elements->gatherClustersDistribution()[info::mpi::rank];

	esint celements = elements->size;
	esint mincindex, gmincindex;
	for (esint c = 0; c < elements->nclusters; c++) {
		avgelements = 0;
		for (esint d = 0; d < elements->ndomains; d++) {
			if (elements->clusters[d] == c) {
				avgelements += elements->elementsDistribution[d + 1] - elements->elementsDistribution[d];
			}
		}
		celements = std::min(avgelements, celements);
		if (avgelements == celements) {
			mincindex = coffset + c;
		}
	}
	MPI_Allreduce(&celements, &minElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, info::mpi::comm);
	if (minElements != celements) {
		mincindex = 0;
	}
	MPI_Reduce(&mincindex, &gmincindex, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	for (esint c = 0; c < elements->nclusters; c++) {
		avgelements = 0;
		for (esint d = 0; d < elements->ndomains; d++) {
			if (elements->clusters[d] == c) {
				avgelements += elements->elementsDistribution[d + 1] - elements->elementsDistribution[d];
			}
		}
		celements = std::max(avgelements, celements);
	}
	MPI_Reduce(&celements, &maxelements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << std::string(namesize - 18, ' ') << "MIN[#n], MAX (AVG)" << " : "
			<< minElements << "[" << gmincindex << "], " << maxElements << " (" << totalElements / (double)totalClusters << ")";

	if (maxElements / (double)minElements > 2) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxElements / (double)minElements;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxElements / (double)minElements;
	}

	ESINFO(DETAILS) << header(" elements per domain");

	minelements = maxelements = elements->elementsDistribution[1];
	avgelements = 0;
	for (esint d = 0; d < elements->ndomains; d++) {
		minelements = std::min(minelements, elements->elementsDistribution[d + 1] - elements->elementsDistribution[d]);
		maxelements = std::max(maxelements, elements->elementsDistribution[d + 1] - elements->elementsDistribution[d]);
		avgelements += elements->elementsDistribution[d + 1] - elements->elementsDistribution[d];
	}
	MPI_Reduce(&minelements, &minElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&avgelements, &totalElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&maxelements, &maxElements, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minElements << ", " << maxElements << " (" << totalElements / (double)totalDomains << ")";

	if (maxElements / (double)minElements > 1.5) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxElements / (double)minElements;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxElements / (double)minElements;
	}

	esint totalNodes = 0, minNodes = 0, maxNodes = 0;
	esint minnodes = nodes->size, maxnodes = nodes->size, avgnodes = nodes->size;
	MPI_Reduce(&minnodes, &minNodes, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&avgnodes, &totalNodes, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&maxnodes, &maxNodes, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);

	ESINFO(DETAILS);

	ESINFO(DETAILS) << header(" NUMBER OF NODES") << nodes->uniqueTotalSize;
	ESINFO(DETAILS) << header(" nodes per MPI") << totalNodes;
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minNodes << ", " << maxNodes << " (" << totalNodes / (double)info::mpi::size << ")";
	ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxNodes / (double)minNodes;
	if (maxNodes / (double)minNodes > 1.5) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxNodes / (double)minNodes;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxNodes / (double)minNodes;
	}


	minnodes = maxnodes = avgnodes = nodes->dintervals[0].back().DOFOffset + nodes->dintervals[0].back().end - nodes->dintervals[0].back().begin;
	for (esint d = 1; d < elements->ndomains; d++) {
		minnodes = std::min(minnodes, nodes->dintervals[d].back().DOFOffset + nodes->dintervals[d].back().end - nodes->dintervals[d].back().begin);
		maxnodes = std::max(maxnodes, nodes->dintervals[d].back().DOFOffset + nodes->dintervals[d].back().end - nodes->dintervals[d].back().begin);
		avgnodes += nodes->dintervals[d].back().DOFOffset + nodes->dintervals[d].back().end - nodes->dintervals[d].back().begin;
	}
	MPI_Reduce(&minnodes, &minNodes, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&avgnodes, &totalNodes, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&maxnodes, &maxNodes, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);

	ESINFO(DETAILS) << header(" nodes per domain") << totalNodes;
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minNodes << ", " << maxNodes << " (" << totalNodes / (double)totalDomains << ")";
	if (maxNodes / (double)minNodes > 1.5) {
		ESINFO(DETAILS) << Info::TextColor::YELLOW << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxNodes / (double)minNodes;
	} else {
		ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxNodes / (double)minNodes;
	}

	esint minUnique = 0, maxUnique = 0;
	MPI_Reduce(&nodes->uniqueSize, &minUnique, sizeof(esint), MPI_BYTE, MPITools::esintOperations().min, 0, info::mpi::comm);
	MPI_Reduce(&nodes->uniqueSize, &maxUnique, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);

	ESINFO(DETAILS) << header(" unique nodes per MPI");
	ESINFO(DETAILS) << std::string(namesize - 14, ' ') << "MIN, MAX (AVG)" << " : "
			<< minUnique << ", " << maxUnique << " (" << maxUnique / (double)nodes->uniqueTotalSize << ")";
	ESINFO(DETAILS) << std::string(namesize - 17, ' ') << "ratio (MAX / MIN)" << " : " << maxUnique / (double)minUnique;


	ESINFO(DETAILS);

	esint maxduplicity = 0, maxDuplicity = 0;
	for (auto d = nodes->idomains->cbegin(); d != nodes->idomains->cend(); ++d) {
		maxduplicity = std::max(maxduplicity, (esint)d->size());
	}

	MPI_Reduce(&maxduplicity, &maxDuplicity, sizeof(esint), MPI_BYTE, MPITools::esintOperations().max, 0, info::mpi::comm);
	ESINFO(DETAILS) << header(" MAX NODE DUPLICITY") << maxDuplicity;

	ESINFO(DETAILS) << header(" COMMUNICATION VOLUME") << "BOUNDARY, INNER (RATIO)";
	esint inner = 0, boundary = 0, totalInner = 0, totalBoundary = 0;
	auto iranks = nodes->iranks->cbegin();
	for (size_t i = 0; i < nodes->pintervals.size(); ++i, ++iranks) {
		if (iranks->front() != info::mpi::rank || iranks->back() != info::mpi::rank) {
			if (nodes->pintervals[i].sourceProcess == info::mpi::rank)
			boundary += nodes->pintervals[i].end - nodes->pintervals[i].begin;
		} else {
			inner += nodes->pintervals[i].end - nodes->pintervals[i].begin;
		}
	}
	MPI_Reduce(&inner, &totalInner, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&boundary, &totalBoundary, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	ESINFO(DETAILS) << std::string(namesize - 7, ' ') << "per MPI : " << totalBoundary << ", " << totalInner << " (" << totalBoundary / (double)totalInner << ")";

	inner = 0, boundary = 0, totalInner = 0, totalBoundary = 0;
	for (esint d = 0; d < elements->ndomains; d++) {
		for (size_t i = 0; i < nodes->dintervals[d].size(); i++) {
			if ((nodes->idomains->cbegin() + nodes->dintervals[d][i].pindex)->size() > 1) {
				if ((nodes->idomains->cbegin() + nodes->dintervals[d][i].pindex)->front() == d + elements->firstDomain) {
					boundary += nodes->dintervals[d][i].end - nodes->dintervals[d][i].begin;
				}
			} else {
				inner += nodes->dintervals[d][i].end - nodes->dintervals[d][i].begin;
			}
		}
	}
	MPI_Reduce(&inner, &totalInner, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	MPI_Reduce(&boundary, &totalBoundary, sizeof(esint), MPI_BYTE, MPITools::esintOperations().sum, 0, info::mpi::comm);
	ESINFO(DETAILS) << std::string(namesize - 10, ' ') << "per domain : " << totalBoundary << ", " << totalInner << " (" << totalBoundary / (double)totalInner << ")";

	ESINFO(DETAILS) << "============================================\n";
}

