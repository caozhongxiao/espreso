
#include "mesh.h"

#include "store/statisticsstore.h"
#include "store/elementstore.h"
#include "store/nodestore.h"
#include "store/elementsregionstore.h"
#include "store/boundaryregionstore.h"

#include "preprocessing/meshpreprocessing.h"

#include "elements/elements.h"

#include "../basis/utilities/utils.h"
#include "../basis/utilities/communication.h"
#include "../basis/utilities/parser.h"

#include "../config/ecf/ecf.h"

#include "../assembler/step.h"

#include "../old/mesh/structures/mesh.h"
#include "../old/mesh/structures/coordinates.h"
#include "../old/mesh/structures/region.h"
#include "../old/mesh/structures/elementtypes.h"
#include "../old/mesh/elements/element.h"

#include <iostream>
#include <vector>
#include <numeric>

#include "../basis/containers/serializededata.h"
#include "../basis/containers/tarray.h"
#include "../output/result/visualization/vtklegacy.h"


using namespace espreso;


Mesh::Mesh(const ECFConfiguration &configuration)
: elements(new ElementStore(_eclasses)), nodes(new NodeStore()),
  FETIData(NULL),
  halo(new ElementStore(_eclasses)),
  domainsSurface(NULL),
  preprocessing(new MeshPreprocessing(this)),

  configuration(configuration),
  _eclasses(environment->OMP_NUM_THREADS),
  mesh(new OldMesh())
{
	size_t threads = environment->OMP_NUM_THREADS;

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
		if (StringCompare::caseInsensitiveEq(elementsRegions[r]->name, name)) {
			return elementsRegions[r];
		}
	}
	ESINFO(ERROR) << "ESPRESO internal error: request for unknown region of elements with '" << name << "'.";
	return NULL;
}

BoundaryRegionStore* Mesh::bregion(const std::string &name)
{
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (StringCompare::caseInsensitiveEq(boundaryRegions[r]->name, name)) {
			return boundaryRegions[r];
		}
	}
	ESINFO(ERROR) << "ESPRESO internal error: request for unknown boundary region '" << name << "'.";
	return NULL;
}

void Mesh::load()
{
	size_t threads = environment->OMP_NUM_THREADS;

	neighbours = mesh->neighbours();
	neighboursWithMe = neighbours;
	neighboursWithMe.push_back(environment->MPIrank);
	std::sort(neighboursWithMe.begin(), neighboursWithMe.end());

	std::vector<eslocal> shrink(mesh->nodes().size());
	// LOAD NODES
	{
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh->nodes().size());
		std::vector<std::vector<Point> > coordinates(threads);
		std::vector<std::vector<eslocal> > IDs(threads);
		std::vector<std::vector<eslocal> > ranksBoundaries(threads);
		std::vector<std::vector<int> > ranksData(threads);
		std::vector<eslocal> tempty(threads);

		ranksBoundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			eslocal empty = 0;
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				if (mesh->nodes()[n]->parentElements().size()) {
					coordinates[t].push_back(mesh->coordinates()[n]);
					IDs[t].push_back(mesh->coordinates().globalIndex(n));
					ranksBoundaries[t].push_back(offset = offset + mesh->nodes()[n]->clusters().size());
					for (size_t c = 0; c < mesh->nodes()[n]->clusters().size(); c++) {
						ranksData[t].push_back(mesh->nodes()[n]->clusters()[c]);
					}
					shrink[n] = empty;
				} else {
					++empty;
					shrink[n] = -1;
				}

			}
			tempty[t] = empty;
		}

		Esutils::sizesToOffsets(tempty);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			for (size_t n = distribution[t]; n < distribution[t + 1]; n++) {
				if (shrink[n] != -1 ) {
					shrink[n] += tempty[t];
				}
			}
		}

		Esutils::threadDistributionToFullDistribution(ranksBoundaries);

		nodes->IDs = new serializededata<eslocal, esglobal>(1, IDs);

		nodes->size = nodes->IDs->structures();
		nodes->distribution = nodes->IDs->datatarray().distribution();

		nodes->coordinates = new serializededata<eslocal, Point>(1, coordinates);
		nodes->ranks = new serializededata<eslocal, int>(ranksBoundaries, ranksData);
	}

	auto geteoffset = [] (eslocal vtkCode) {
		switch (vtkCode) {
		case  3: return static_cast<int>(Element::CODE::LINE2);
		case  4: return static_cast<int>(Element::CODE::LINE3);

		case  5: return static_cast<int>(Element::CODE::TRIANGLE3);
		case  9: return static_cast<int>(Element::CODE::SQUARE4);
		case 22: return static_cast<int>(Element::CODE::TRIANGLE6);
		case 23: return static_cast<int>(Element::CODE::SQUARE8);

		case 10: return static_cast<int>(Element::CODE::TETRA4);
		case 12: return static_cast<int>(Element::CODE::HEXA8);
		case 13: return static_cast<int>(Element::CODE::PRISMA6);
		case 14: return static_cast<int>(Element::CODE::PYRAMID5);

		case 24: return static_cast<int>(Element::CODE::TETRA10);
		case 25: return static_cast<int>(Element::CODE::HEXA20);
		case 26: return static_cast<int>(Element::CODE::PRISMA15);
		case 27: return static_cast<int>(Element::CODE::PYRAMID13);
		}
		return -1;
	};

	{
		size_t esize = mesh->elements().size();
		Communication::exscan(esize, MPITools::operations().sizeToOffsets);
		std::vector<size_t> distribution = tarray<eslocal>::distribute(threads, mesh->elements().size());
		std::vector<std::vector<eslocal> > boundaries(threads), indices(threads);
		std::vector<std::vector<Element*> > epointers(threads);
		std::vector<std::vector<eslocal> > eIDs(threads);
		std::vector<std::vector<int> > body(threads), material(threads);

		boundaries.front().push_back(0);
		#pragma omp parallel for
		for (size_t t = 0; t < threads; t++) {
			size_t offset = 0;
			for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {
				epointers[t].push_back(&_eclasses[t][geteoffset(mesh->elements()[e]->vtkCode())]);
				boundaries[t].push_back(offset = offset + mesh->elements()[e]->nodes());
				for (size_t n = 0; n < mesh->elements()[e]->nodes(); n++) {
					indices[t].push_back(mesh->elements()[e]->node(n) - shrink[mesh->elements()[e]->node(n)]);
				}

				eIDs[t].push_back(e + esize);
				body[t].push_back(mesh->elements()[e]->param(OldElement::Params::BODY));
				material[t].push_back(mesh->elements()[e]->param(OldElement::Params::MATERIAL));
			}
		}

		Esutils::threadDistributionToFullDistribution(boundaries);

		elements->size = mesh->elements().size();
		elements->distribution = distribution;

		elements->IDs = new serializededata<eslocal, eslocal>(1, eIDs);
		elements->nodes = new serializededata<eslocal, eslocal>(std::move(tarray<eslocal>(boundaries)), std::move(tarray<eslocal>(indices)));

		elements->body = new serializededata<eslocal, int>(1, body);
		elements->material = new serializededata<eslocal, int>(1, material);
		elements->epointers = new serializededata<eslocal, Element*>(1, std::move(tarray<Element*>(epointers)));
	}

	for (size_t r = 0; r < mesh->regions().size(); r++) {
		std::vector<size_t> tdistributions = tarray<size_t>::distribute(threads, mesh->regions()[r]->elements().size());
		std::vector<std::vector<eslocal> > rdistribution(threads), rdata(threads);
		std::vector<std::vector<Element*> > epointers(threads);

		switch (mesh->regions()[r]->eType) {
		case ElementType::ELEMENTS:
			elementsRegions.push_back(new ElementsRegionStore(mesh->regions()[r]->name));
			if (mesh->regions()[r]->elements().size() == mesh->elements().size()) {
				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					rdata[t].resize(tdistributions[t + 1] - tdistributions[t]);
					std::iota(rdata[t].begin(), rdata[t].end(), tdistributions[t]);
				}
			} else {
				#pragma omp parallel for
				for (size_t t = 0; t < threads; t++) {
					for (size_t e = tdistributions[t]; e < tdistributions[t + 1]; e++) {
						rdata[t].push_back(std::find(mesh->elements().begin(), mesh->elements().end(), mesh->regions()[r]->elements()[e]) - mesh->elements().begin());
					}
				}
			}
			elementsRegions.back()->elements = new serializededata<eslocal, eslocal>(1, rdata);
			std::sort(elementsRegions.back()->elements->datatarray().begin(), elementsRegions.back()->elements->datatarray().end());
			break;
		case ElementType::FACES:
		case ElementType::EDGES:
			rdistribution[0].push_back(0);
			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t e = tdistributions[t]; e < tdistributions[t + 1]; e++) {
					for (size_t n = 0; n < mesh->regions()[r]->elements()[e]->nodes(); n++) {
						rdata[t].push_back(mesh->regions()[r]->elements()[e]->node(n) - shrink[mesh->regions()[r]->elements()[e]->node(n)]);
					}
					epointers[t].push_back(&_eclasses[t][geteoffset(mesh->regions()[r]->elements()[e]->vtkCode())]);
					rdistribution[t].push_back(rdata[t].size());
				}
			}
			Esutils::threadDistributionToFullDistribution(rdistribution);

			boundaryRegions.push_back(new BoundaryRegionStore(mesh->regions()[r]->name, _eclasses));
			boundaryRegions.back()->distribution = tdistributions;
			boundaryRegions.back()->elements = new serializededata<eslocal, eslocal>(rdistribution, rdata);
			boundaryRegions.back()->epointers = new serializededata<eslocal, Element*>(1, epointers);
			if (mesh->regions()[r]->eType == ElementType::FACES) {
				boundaryRegions.back()->dimension = 2;
			} else {
				boundaryRegions.back()->dimension = 1;
			}
			break;
		case ElementType::NODES:

			#pragma omp parallel for
			for (size_t t = 0; t < threads; t++) {
				for (size_t e = tdistributions[t]; e < tdistributions[t + 1]; e++) {
					if (shrink[mesh->regions()[r]->elements()[e]->node(0)] != -1) {
						rdata[t].push_back(mesh->regions()[r]->elements()[e]->node(0) - shrink[mesh->regions()[r]->elements()[e]->node(0)]);
					}
				}
			}

			boundaryRegions.push_back(new BoundaryRegionStore(mesh->regions()[r]->name, _eclasses));
			boundaryRegions.back()->nodes = new serializededata<eslocal, eslocal>(1, rdata);
			std::sort(boundaryRegions.back()->nodes->datatarray().begin(), boundaryRegions.back()->nodes->datatarray().end());
			break;
		}
	}

	update();
}


void Mesh::update()
{
	ESINFO(OVERVIEW) << "Preprocess mesh data.";
	materials.clear();
	std::map<std::string, int> matindex;
	for (auto mat = configuration.getPhysics()->materials.begin(); mat != configuration.getPhysics()->materials.end(); ++mat) {
		materials.push_back(&mat->second);
		matindex[mat->first] = materials.size() - 1;
	}

	for (auto mat = configuration.getPhysics()->material_set.begin(); mat != configuration.getPhysics()->material_set.end(); ++mat) {
		ElementsRegionStore *region = eregion(mat->first);
		if (matindex.find(mat->second) == matindex.end()) {
			ESINFO(GLOBAL_ERROR) << "Unknown material '" << mat->second << "'.";
		}
		int material = matindex.find(mat->second)->second;
		for (auto e = region->elements->datatarray().cbegin(); e != region->elements->datatarray().cend(); ++e) {
			elements->material->datatarray()[*e] = material;
		}
	}

	if (configuration.decomposition.balance_clusters || configuration.input == INPUT_FORMAT::WORKBENCH) {
		preprocessing->reclusterize();
	}

	auto ntob = [&] (const std::string &rname, int dimension) {
		BoundaryRegionStore *region = bregion(rname);
		if (region->dimension == 0) {
			preprocessing->computeBoundaryElementsFromNodes(region, dimension);
		}
	};

	if (configuration.physics == PHYSICS::STRUCTURAL_MECHANICS_2D || configuration.physics == PHYSICS::STRUCTURAL_MECHANICS_3D) {
		const StructuralMechanicsConfiguration *sm;
		int dimension;
		if (configuration.physics == PHYSICS::STRUCTURAL_MECHANICS_2D) {
			sm = &configuration.structural_mechanics_2d;
			dimension = 1;
		}
		if (configuration.physics == PHYSICS::STRUCTURAL_MECHANICS_3D) {
			sm = &configuration.structural_mechanics_3d;
			dimension = 2;
		}

		for (auto ls = sm->load_steps_settings.begin(); ls != sm->load_steps_settings.end(); ++ls) {
			for (auto bc = ls->second.normal_pressure.begin(); bc != ls->second.normal_pressure.end(); ++bc) {
				ntob(bc->first, dimension);
			}
		}
	}

	if (configuration.physics == PHYSICS::HEAT_TRANSFER_2D || configuration.physics == PHYSICS::HEAT_TRANSFER_3D) {
		const HeatTransferConfiguration *ht;
		int dimension;
		if (configuration.physics == PHYSICS::HEAT_TRANSFER_2D) {
			ht = &configuration.heat_transfer_2d;
			dimension = 1;
		}
		if (configuration.physics == PHYSICS::HEAT_TRANSFER_3D) {
			ht = &configuration.heat_transfer_3d;
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

	bool uniformDecomposition = false;
	if (configuration.input == INPUT_FORMAT::GENERATOR) {
		switch (configuration.generator.shape) {
		case INPUT_GENERATOR_SHAPE::GRID:
			if (configuration.generator.grid.uniform_decomposition) {
				uniformDecomposition = true;
			}
			break;
		case INPUT_GENERATOR_SHAPE::GRID_TOWER:
			uniformDecomposition = true;
			for (auto it = configuration.generator.grid_tower.grids.begin(); it != configuration.generator.grid_tower.grids.end(); ++it) {
				if (!it->second.uniform_decomposition) {
					uniformDecomposition = false;
				}
			}
			break;
		case INPUT_GENERATOR_SHAPE::SPHERE:
			if (configuration.generator.sphere.uniform_decomposition) {
				uniformDecomposition = true;
			}
			break;
		}
	}
	if (
			configuration.decomposition.separate_materials ||
			configuration.decomposition.separate_regions ||
			configuration.decomposition.separate_etypes) {
		uniformDecomposition = false;
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
		switch (configuration.physics) {
		case PHYSICS::HEAT_TRANSFER_2D:
			return configuration.heat_transfer_2d;
		case PHYSICS::HEAT_TRANSFER_3D:
			return configuration.heat_transfer_3d;
		case PHYSICS::STRUCTURAL_MECHANICS_2D:
			return configuration.structural_mechanics_2d;
		case PHYSICS::STRUCTURAL_MECHANICS_3D:
			return configuration.structural_mechanics_3d;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented physics.";
			exit(0);
		}
	};

	auto forEachSteps = [&] (std::function<bool(const LoadStepConfiguration &)> fnc) {
		bool ret = false;
		switch (configuration.physics) {
		case PHYSICS::HEAT_TRANSFER_2D:
			for (auto step = configuration.heat_transfer_2d.load_steps_settings.begin(); step != configuration.heat_transfer_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PHYSICS::HEAT_TRANSFER_3D:
			for (auto step = configuration.heat_transfer_3d.load_steps_settings.begin(); step != configuration.heat_transfer_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PHYSICS::STRUCTURAL_MECHANICS_2D:
			for (auto step = configuration.structural_mechanics_2d.load_steps_settings.begin(); step != configuration.structural_mechanics_2d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		case PHYSICS::STRUCTURAL_MECHANICS_3D:
			for (auto step = configuration.structural_mechanics_3d.load_steps_settings.begin(); step != configuration.structural_mechanics_3d.load_steps_settings.end(); ++step) {
				ret |= fnc(step->second);
			}
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "Not implemented physics.";
			exit(0);
		}
		return ret;
	};

	uniformDecomposition = false;
	if (uniformDecomposition) {
		// implement uniform decomposition
	} else {
		preprocessing->partitiate(configuration.decomposition.domains ? configuration.decomposition.domains : mesh->parts(),
				configuration.decomposition.separate_materials || hasBEM(getPhysics()), // BEM domain has to have only one material
				configuration.decomposition.separate_regions || hasBEM(getPhysics()),
				configuration.decomposition.separate_etypes);
	}

	if (hasBEM(getPhysics())) {
		preprocessing->computeDomainsSurface();
		preprocessing->triangularizeDomainSurface();
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

	if (configuration.physics == PHYSICS::STRUCTURAL_MECHANICS_2D || configuration.physics == PHYSICS::STRUCTURAL_MECHANICS_3D) {
		if (forEachSteps([] (const LoadStepConfiguration &step) {
			return step.solver == LoadStepConfiguration::SOLVER::FETI && step.feti.regularization == FETI_REGULARIZATION::ANALYTIC;
		})) {

			preprocessing->computeFixPoints();
			if (hasBEM(getPhysics())) {
				preprocessing->computeFixPointsOnSurface();
			}
		}
	}

	if (configuration.mesh_morphing.type == MORPHING_TYPE::RBF) {
		for (auto it = configuration.mesh_morphing.rbf.begin(); it != configuration.mesh_morphing.rbf.end(); ++it) {
			switch (configuration.physics) {
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

	printStatistics();
}

void Mesh::initNodeData()
{
	for (auto datait = nodes->data.begin(); datait != nodes->data.end(); ++datait) {
		NodeData* data = *datait;
		if (data->names.size() && data->decomposedData != NULL) {
			data->gatheredData.resize(data->dimension * nodes->uniqueSize);
			data->sBuffer.resize(neighbours.size());
			data->rBuffer.resize(neighbours.size());

			for (size_t n = 0; n < neighbours.size(); ++n) {
				if (neighbours[n] < environment->MPIrank) {
					data->sBuffer[n].resize(data->dimension * nodes->scouters[n]);
				} else {
					data->rBuffer[n].resize(data->dimension * nodes->scouters[n + 1]);
				}
			}
		}
	}
}

void Mesh::gatherNodeData()
{
	// TODO: NUMA + load balancing

	auto n2i = [ & ] (int neighbour) {
		return std::lower_bound(neighbours.begin(), neighbours.end(), neighbour) - neighbours.begin();
	};

	auto doffset = [&] (eslocal d, eslocal i) {
		return std::lower_bound(nodes->dintervals[d].begin(), nodes->dintervals[d].end(), i, [] (const DomainInterval &interval, eslocal i) { return interval.pindex < i; })->DOFOffset;
	};

	for (auto datait = nodes->data.begin(); datait != nodes->data.end(); ++datait) {
		NodeData* data = *datait;
		if (data->names.size() && data->decomposedData != NULL) {
			for (size_t i = 0; i < data->sBuffer.size(); i++) {
				std::fill(data->sBuffer[i].begin(), data->sBuffer[i].end(), 0);
			}
			std::fill(data->gatheredData.begin(), data->gatheredData.end(), 0);

			#pragma omp parallel for
			for (size_t i = 0; i < nodes->pintervals.size(); ++i) {
				auto domains = nodes->idomains->cbegin() + i;
				eslocal offset, soffset, noffset;

				if (nodes->pintervals[i].sourceProcess < environment->MPIrank) {
					noffset = n2i(nodes->pintervals[i].sourceProcess);
					for (auto d = domains->begin(); d != domains->end(); ++d) {
						if (elements->firstDomain <= *d && *d < elements->firstDomain + elements->ndomains) {
							offset = data->dimension * doffset(*d - elements->firstDomain, i);
							soffset = data->dimension * nodes->soffsets[i];
							for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n) {
								for (int dof = 0; dof < data->dimension; ++dof, ++soffset, ++offset) {
									data->sBuffer[noffset][soffset] += (*data->decomposedData)[*d - elements->firstDomain][offset];
								}
							}
						}
					}
				}
			}

			if (!Communication::receiveUpperKnownSize(data->sBuffer, data->rBuffer, neighbours)) {
				ESINFO(ERROR) << "ESPRESO internal error: gather results\n";
			}

			#pragma omp parallel for
			for (size_t i = 0; i < nodes->pintervals.size(); ++i) {
				auto idomains = nodes->idomains->cbegin() + i;
				auto ineighbors = nodes->ineighborOffsets->cbegin() + i;

				eslocal offset, goffset, noffset;

				if (nodes->pintervals[i].sourceProcess == environment->MPIrank) {
					for (auto d = idomains->begin(); d != idomains->end() && *d < elements->firstDomain + elements->ndomains; ++d) {
						goffset = data->dimension * (nodes->pintervals[i].globalOffset - nodes->uniqueOffset);
						offset = data->dimension * doffset(*d - elements->firstDomain, i);
						for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n) {
							for (int dof = 0; dof < data->dimension; ++dof, ++goffset, ++offset) {
								data->gatheredData[goffset] += (*data->decomposedData)[*d - elements->firstDomain][offset];
							}
						}
					}
					for (auto neigh = ineighbors->begin(); neigh != ineighbors->end(); ++neigh) {
						goffset = data->dimension * (nodes->pintervals[i].globalOffset - nodes->uniqueOffset);
						offset = data->dimension * neigh->offset;
						noffset = n2i(neigh->process);
						for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n) {
							for (int dof = 0; dof < data->dimension; ++dof, ++goffset, ++offset) {
								data->gatheredData[goffset] += data->rBuffer[noffset][offset];
							}
						}
					}
					goffset = data->dimension * (nodes->pintervals[i].globalOffset - nodes->uniqueOffset);
					for (eslocal n = nodes->pintervals[i].begin; n < nodes->pintervals[i].end; ++n) {
						for (int dof = 0; dof < data->dimension; ++dof, ++goffset) {
							data->gatheredData[goffset] /= idomains->size();
						}
					}
				}
			}
		}
	}
}

double Mesh::sumSquares(const std::vector<std::vector<double> > &data, const BoundaryRegionStore* region)
{
	size_t threads = environment->OMP_NUM_THREADS;

	double csum = 0, gsum;

	int dimension = data.front().size() / (nodes->dintervals.front().back().DOFOffset + nodes->dintervals.front().back().end - nodes->dintervals.front().back().begin);
	#pragma omp parallel for
	for (size_t t = 0; t < threads; t++) {
		double tsum = 0, value;
		eslocal offset;

		for (eslocal d = elements->domainDistribution[t]; d < elements->domainDistribution[t + 1]; d++) {
			for (size_t di = 0; di < nodes->dintervals[d].size(); di++) {
				offset = dimension * nodes->pintervals[nodes->dintervals[d][di].pindex].begin;
				for (eslocal i = dimension * region->nintervals[nodes->dintervals[d][di].pindex].begin; i < dimension * region->nintervals[nodes->dintervals[d][di].pindex].end; ++i) {
					value = data[d][nodes->dintervals[d][di].DOFOffset + region->nodes->datatarray()[i] - offset];
					tsum += value * value;
				}
			}
		}

		#pragma omp atomic
		csum += tsum;
	}

	MPI_Allreduce(&csum, &gsum, 1, MPI_DOUBLE, MPI_SUM, environment->MPICommunicator);
	return gsum;
}

static void _computeGatheredNodeStatistics(const Mesh *mesh, const NodeData *data, const std::vector<ProcessInterval> &intervals, const tarray<eslocal> &nodes, eslocal offset, Statistics *statistics, MPI_Comm communicator)
{
	eslocal goffset, index;

	if (data->names.size() == 1) {
		statistics->reset();
		for (size_t i = 0; i < intervals.size(); i++) {
			if (intervals[i].sourceProcess == environment->MPIrank) {
				goffset = mesh->nodes->pintervals[i].globalOffset - offset;
				for (auto n = nodes.cbegin() + intervals[i].begin; n != nodes.cbegin() + intervals[i].end; ++n) {
					index = goffset + *n - mesh->nodes->pintervals[i].begin;
					statistics->min = std::min(statistics->min, data->gatheredData[index]);
					statistics->max = std::max(statistics->max, data->gatheredData[index]);
					statistics->avg += data->gatheredData[index];
					statistics->norm += data->gatheredData[index] * data->gatheredData[index];
				}
			}
		}
	} else {
		double value;
		for (int d = 0; d <= data->dimension; d++) {
			(statistics + d)->reset();
		}

		for (size_t i = 0; i < intervals.size(); i++) {
			if (intervals[i].sourceProcess == environment->MPIrank) {
				goffset = mesh->nodes->pintervals[i].globalOffset - offset;
				for (auto n = nodes.cbegin() + intervals[i].begin; n != nodes.cbegin() + intervals[i].end; ++n) {
					value = 0;
					index = goffset + *n - mesh->nodes->pintervals[i].begin;
					for (int d = 0; d < data->dimension; d++) {
						value +=  data->gatheredData[index * data->dimension + d] * data->gatheredData[index * data->dimension + d];
						(statistics + d + 1)->min = std::min((statistics + d + 1)->min, data->gatheredData[index * data->dimension + d]);
						(statistics + d + 1)->max = std::max((statistics + d + 1)->max, data->gatheredData[index * data->dimension + d]);
						(statistics + d + 1)->avg += data->gatheredData[index * data->dimension + d];
						(statistics + d + 1)->norm += data->gatheredData[index * data->dimension + d] * data->gatheredData[index * data->dimension + d];
					}
					value = std::sqrt(value);
					statistics->min = std::min(statistics->min, value);
					statistics->max = std::max(statistics->max, value);
					statistics->avg += value;
					statistics->norm += value * value;
				}
			}
		}
	}

	std::vector<Statistics> global(data->names.size());
	MPI_Reduce(statistics, global.data(), sizeof(Statistics) * data->names.size(), MPI_BYTE, MPITools::operations().mergeStatistics, 0, communicator);
	memcpy(statistics, global.data(), sizeof(Statistics) * data->names.size());
}

void Mesh::computeGatheredNodeStatistic(const NodeData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const
{
	_computeGatheredNodeStatistics(this, data, region->nintervals, region->nodes->datatarray(), nodes->uniqueOffset, statistics, communicator);
	for (size_t i = 0; i < data->names.size(); i++) {
		(statistics + i)->avg /= region->uniqueTotalSize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}

void Mesh::computeGatheredNodeStatistic(const NodeData *data, const BoundaryRegionStore* region, Statistics *statistics, MPI_Comm communicator) const
{
	_computeGatheredNodeStatistics(this, data, region->nintervals, region->nodes->datatarray(), nodes->uniqueOffset, statistics, communicator);
	for (size_t i = 0; i < data->names.size(); i++) {
		(statistics + i)->avg /= region->uniqueTotalSize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}

void Mesh::computeElementStatistic(const ElementData *data, const ElementsRegionStore* region, Statistics *statistics, MPI_Comm communicator) const
{
	if (data->names.size() == 1) {
		statistics->reset();
		for (size_t i = 0; i < region->eintervals.size(); i++) {
			for (auto e = region->elements->datatarray().cbegin() + region->eintervals[i].begin; e != region->elements->datatarray().cbegin() + region->eintervals[i].end; ++e) {
				statistics->min = std::min(statistics->min, (*data->data)[*e]);
				statistics->max = std::max(statistics->max, (*data->data)[*e]);
				statistics->avg += (*data->data)[*e];
				statistics->norm += (*data->data)[*e] * (*data->data)[*e];
			}
		}
	} else {
		double value;
		for (int d = 0; d <= data->dimension; d++) {
			(statistics + d)->reset();
		}

		for (size_t i = 0; i < region->eintervals.size(); i++) {
			for (auto e = region->elements->datatarray().cbegin() + region->eintervals[i].begin; e != region->elements->datatarray().cbegin() + region->eintervals[i].end; ++e) {
				value = 0;
				for (int d = 0; d < data->dimension; d++) {
					value += (*data->data)[*e * data->dimension + d] * (*data->data)[*e * data->dimension + d];
					(statistics + d + 1)->min = std::min((statistics + d + 1)->min, (*data->data)[*e * data->dimension + d]);
					(statistics + d + 1)->max = std::max((statistics + d + 1)->max, (*data->data)[*e * data->dimension + d]);
					(statistics + d + 1)->avg += (*data->data)[*e * data->dimension + d];
					(statistics + d + 1)->norm += (*data->data)[*e * data->dimension + d] * (*data->data)[*e * data->dimension + d];
				}
				value = std::sqrt(value);
				statistics->min = std::min(statistics->min, value);
				statistics->max = std::max(statistics->max, value);
				statistics->avg += value;
				statistics->norm += value * value;
			}
		}
	}

	std::vector<Statistics> global(data->names.size());
	MPI_Reduce(statistics, global.data(), sizeof(Statistics) * data->names.size(), MPI_BYTE, MPITools::operations().mergeStatistics, 0, communicator);
	memcpy(statistics, global.data(), sizeof(Statistics) * data->names.size());
	for (size_t i = 0; i < data->names.size(); i++) {
		(statistics + i)->avg /= region->uniqueTotalSize;
		(statistics + i)->norm = std::sqrt((statistics + i)->norm);
	}
}

void Mesh::printStatistics()
{
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
			ESINFO(OVERVIEW) << std::string(21 - ename(code).size(), ' ') << ename(code) << " : " << elements->ecounters[code];
		}
	};

	auto totalesize = [] (std::vector<eslocal> &ecounters) {
		eslocal size = 0;
		for (int etype = 0; etype < ecounters.size(); etype++) {
			size += ecounters[etype];
		}
		return size;
	};

	ESINFO(OVERVIEW) << "============= Mesh statistics =============";

	ESINFO(OVERVIEW) << " Number of nodes      : " << nodes->uniqueTotalSize;

	ESINFO(OVERVIEW) << " Number of elements   : " << totalesize(elements->ecounters);
	for (int etype = 0; etype < static_cast<int>(Element::CODE::SIZE); etype++) {
		esize(etype);
	}

	ESINFO(OVERVIEW);

	ESINFO(OVERVIEW) << " Element regions size :";
	ESINFO(OVERVIEW) << std::string(21 - elementsRegions[0]->name.size(), ' ') << elementsRegions[0]->name << " : " << totalesize(elements->ecounters);
	for (size_t r = 1; r < elementsRegions.size(); r++) {
		ESINFO(OVERVIEW) << std::string(21 - elementsRegions[r]->name.size(), ' ') << elementsRegions[r]->name << " : " << totalesize(elementsRegions[r]->ecounters);
	}
	ESINFO(OVERVIEW) << " Face regions size    :";
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->dimension == 2) {
			ESINFO(OVERVIEW) << std::string(21 - boundaryRegions[r]->name.size(), ' ') << boundaryRegions[r]->name << " : " << totalesize(boundaryRegions[r]->ecounters);
		}
	}
	ESINFO(OVERVIEW) << " Edge regions size    :";
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->dimension == 1) {
			ESINFO(OVERVIEW) << std::string(21 - boundaryRegions[r]->name.size(), ' ') << boundaryRegions[r]->name << " : " << totalesize(boundaryRegions[r]->ecounters);
		}
	}
	ESINFO(OVERVIEW) << " Node regions size    :";
	for (size_t r = 0; r < boundaryRegions.size(); r++) {
		if (boundaryRegions[r]->dimension == 0) {
			ESINFO(OVERVIEW) << std::string(21 - boundaryRegions[r]->name.size(), ' ') << boundaryRegions[r]->name << " : " << boundaryRegions[r]->uniqueTotalSize;
		}
	}

	ESINFO(OVERVIEW);

	int totalClusters = 0, clusters = elements->nclusters;
	MPI_Reduce(&clusters, &totalClusters, 1, MPI_INT, MPI_SUM, 0, environment->MPICommunicator);
	ESINFO(OVERVIEW) << " Number of clusters   : " << totalClusters;

	int totalDomains = 0, domains = elements->ndomains;
	MPI_Reduce(&domains, &totalDomains, 1, MPI_INT, MPI_SUM, 0, environment->MPICommunicator);
	ESINFO(OVERVIEW) << " Number of domains    : " << totalDomains;

	ESINFO(OVERVIEW) << "============================================";
}

