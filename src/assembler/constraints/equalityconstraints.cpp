
#include "equalityconstraints.h"

using namespace espreso;

void EqualityConstraints::insertDirichletToB1(const std::vector<Element*> &nodes, const std::vector<Property> &DOFs)
{
	size_t threads = Esutils::getEnv<size_t>("CILK_NWORKERS");
	std::vector<size_t> distribution = Esutils::getDistribution(threads, nodes.size());

	// part x thread x Dirichlet
	std::vector<std::vector<std::vector<esglobal> > > dirichlet(_mesh.parts(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<double> > > dirichletValues(_mesh.parts(), std::vector<std::vector<double> >(threads));

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (nodes[i]->settings().isSet(DOFs[dof])) {
					const std::vector<eslocal>& indices = nodes[i]->DOFsIndices();
					double value = nodes[i]->settings(DOFs[dof]).back()->evaluate(i);
					for(size_t d = 0; d < nodes[i]->domains().size(); d++) {
						if (indices[d * DOFs.size() + dof] != -1) {
							dirichlet[nodes[i]->domains()[d]][t].push_back(indices[d * DOFs.size() + dof] + IJVMatrixIndexing);
							dirichletValues[nodes[i]->domains()[d]][t].push_back(value);
						}
					}
				}
			}

		}
	}

	std::vector<size_t> dirichletSizes(_mesh.parts());
	#pragma cilk grainsize = 1
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		size_t size = 0;
		for (size_t t = 0; t < threads; t++) {
			size += dirichlet[p][t].size();
		}
		dirichletSizes[p] = size;
	}

	size_t clusterOffset = 0;
	std::vector<eslocal> subdomainsWithDirichlet;
	for (size_t p = 0; p < _mesh.parts(); p++) {
		clusterOffset += dirichletSizes[p];
		if (dirichletSizes[p]) {
			subdomainsWithDirichlet.push_back(p);
		}
	}

	size_t clusterDirichletSize = clusterOffset;
	size_t globalDirichletSize = synchronizeOffsets(clusterOffset);

	clusterOffset += B1[0].rows;
	#pragma cilk grainsize = 1
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		B1[p].rows += globalDirichletSize;
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		B1[s].nnz += dirichletSizes[s];
		B1[s].I_row_indices.reserve(B1[s].nnz);
		B1[s].J_col_indices.reserve(B1[s].nnz);
		B1[s].V_values.resize(B1[s].nnz, 1);
	}


	Esutils::sizesToOffsets(dirichletSizes);
	#pragma cilk grainsize = 1
	cilk_for (size_t i = 0; i < subdomainsWithDirichlet.size(); i++) {
		size_t s = subdomainsWithDirichlet[i];
		for (size_t i = 0; i < B1[s].nnz; i++) {
			B1[s].I_row_indices.push_back(clusterOffset + dirichletSizes[s] + i + IJVMatrixIndexing);
		}
		for (size_t t = 0; t < threads; t++) {
			B1[s].J_col_indices.insert(B1[s].J_col_indices.end(), dirichlet[s][t].begin(), dirichlet[s][t].end());
			B1c[s].insert(B1c[s].end(), dirichletValues[s][t].begin(), dirichletValues[s][t].end());
		}
		B1duplicity[s].resize(B1[s].I_row_indices.size(), 1);
		for (size_t r = B1subdomainsMap[s].size(); r < B1[s].nnz; r++) {
			B1subdomainsMap[s].push_back(B1[s].I_row_indices[r] - 1);
		}
	}

	B1clustersMap.reserve(B1clustersMap.size() + clusterDirichletSize);
	for (esglobal i = clusterOffset; i < clusterOffset + clusterDirichletSize; i++) {
		B1clustersMap.push_back({ i, config::env::MPIrank });
	}

	ESINFO(DETAILS) << "Lambdas with Dirichlet in B1: " << B1[0].rows;
	ESTEST(MANDATORY) << "ESPRESO requires some nodes with Dirichlet condition." << (globalDirichletSize == 0 ? TEST_FAILED : TEST_PASSED);
}

std::vector<esglobal> EqualityConstraints::computeLambdasID(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	std::vector<esglobal> lambdasID(elements.size() * DOFs.size(), -1);

	auto n2i = [ & ] (size_t neighbour) {
		return std::lower_bound(_mesh.neighbours().begin(), _mesh.neighbours().end(), neighbour) - _mesh.neighbours().begin();
	};

	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, elements.size());
	std::vector<size_t> offsets(threads);

	// neighbours x threads x data
	std::vector<std::vector<std::vector<esglobal> > > sLambdas(_mesh.neighbours().size(), std::vector<std::vector<esglobal> >(threads));
	std::vector<std::vector<std::vector<esglobal> > > rLambdas(_mesh.neighbours().size(), std::vector<std::vector<esglobal> >(threads));


	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		size_t lambdasSize = 0;
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				size_t n = elements[e]->numberOfGlobalDomainsWithDOF(dof);
				if (n > 1 && !elements[e]->settings().isSet(DOFs[dof])) { // Dirichlet and inner nodes are not glued
					if (elements[e]->clusters()[0] == config::env::MPIrank) { // set lambda ID
						lambdasID[e * DOFs.size() + dof] = n * (n - 1) / 2;
						lambdasSize += n * (n - 1) / 2;
						for (size_t c = 1; c < elements[e]->clusters().size(); c++) { // send to higher clusters
							sLambdas[n2i(elements[e]->clusters()[c])][t].push_back(e * DOFs.size() + dof);
						}
					} else { // pick ID from lower cluster
						rLambdas[n2i(elements[e]->clusters()[0])][t].push_back(e * DOFs.size() + dof);
					}
				}
			}

		}
		offsets[t] = lambdasSize;
	}

	size_t numberOfClusterLambdas = Esutils::sizesToOffsets(offsets);
	size_t clusterOffset = numberOfClusterLambdas;
	size_t totalNumberOfLambdas = synchronizeOffsets(clusterOffset);

	for (size_t i = 0; i < offsets.size(); i++) {
		offsets[i] += clusterOffset + B1[0].rows;
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		esglobal offset = offsets[t];
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t dof = 0; dof < DOFs.size(); dof++) {
				if (lambdasID[e * DOFs.size() + dof] > 0) {
					offset += lambdasID[e * DOFs.size() + dof];
					lambdasID[e * DOFs.size() + dof] = offset - lambdasID[e * DOFs.size() + dof];
				}
			}

		}

		for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
			for (size_t i = 0; i < sLambdas[n][t].size(); i++) {
				sLambdas[n][t][i] = lambdasID[sLambdas[n][t][i]];
			}
		}
	}

	std::vector<std::vector<esglobal> > rBuffer(_mesh.neighbours().size());

	#pragma cilk grainsize = 1
	cilk_for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		size_t size = rLambdas[n][0].size();
		for (size_t t = 1; t < threads; t++) {
			sLambdas[n][0].insert(sLambdas[n][0].end(), sLambdas[n][t].begin(), sLambdas[n][t].end());
			size += rLambdas[n][t].size();
		}
		rBuffer[n].resize(size);
	}


	std::vector<MPI_Request> req(_mesh.neighbours().size());
	for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
		if (_mesh.neighbours()[n] > config::env::MPIrank) {
			MPI_Isend(sLambdas[n][0].data(), sizeof(esglobal) * sLambdas[n][0].size(), MPI_BYTE, _mesh.neighbours()[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
		if (_mesh.neighbours()[n] < config::env::MPIrank) {
			MPI_Irecv(rBuffer[n].data(), sizeof(esglobal) * rBuffer[n].size(), MPI_BYTE, _mesh.neighbours()[n], 1, MPI_COMM_WORLD, req.data() + n);
		}
	}

	MPI_Waitall(_mesh.neighbours().size(), req.data(), MPI_STATUSES_IGNORE);

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t e = distribution[t]; e < distribution[t + 1]; e++) {

			for (size_t n = 0; n < _mesh.neighbours().size(); n++) {
				size_t offset = 0;
				for (size_t i = 0; i < t; i++) {
					offset += rLambdas[n][i].size();
				}
				for (size_t i = 0; i < rLambdas[n][t].size(); i++) {
					lambdasID[rLambdas[n][t][i]] = rBuffer[n][offset + i];
				}
			}

		}
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		B1[p].rows += totalNumberOfLambdas;
	}

	return lambdasID;
}

void EqualityConstraints::insertElementGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	std::vector<esglobal> lambdasID = computeLambdasID(elements, DOFs);

	std::vector<eslocal> permutation(lambdasID.size());
	std::iota(permutation.begin(), permutation.end(), 0);

	std::sort(permutation.begin(), permutation.end(), [&] (eslocal i, eslocal j) {
		return (lambdasID[i] < 0) ? false : (lambdasID[j] < 0) ? true : lambdasID[i] < lambdasID[j];
	});

	auto it = std::find_if(permutation.begin(), permutation.end(), [&] (eslocal i) { return lambdasID[i] == -1; });
	permutation.resize(it - permutation.begin());


	size_t threads = config::env::CILK_NWORKERS;
	std::vector<size_t> distribution = Esutils::getDistribution(threads, permutation.size());

	// threads x domains x data
	std::vector<std::vector<std::vector<esglobal> > > rows(threads, std::vector<std::vector<esglobal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > cols(threads, std::vector<std::vector<eslocal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<eslocal> > > vals(threads, std::vector<std::vector<eslocal> >(_mesh.parts()));
	std::vector<std::vector<std::vector<double> > > dup(threads, std::vector<std::vector<double> >(_mesh.parts()));

	std::vector<std::vector<std::vector<esglobal> > > cMap(threads);

	auto findDomain = [&] (const Element *e, size_t d, size_t dof) -> eslocal {
		auto &DOFIndices = e->DOFsIndices();
		size_t c = 0, DOFs = DOFIndices.size() / e->domains().size();
		for (size_t i = 0; i < e->domains().size(); i++) {
			if (DOFIndices[i * DOFs + dof] != -1) {
				if (d == c++) {
					return e->domains()[i];
				}
			}
		}
		return 0;
	};

	#pragma cilk grainsize = 1
	cilk_for (size_t t = 0; t < threads; t++) {
		for (size_t i = distribution[t]; i < distribution[t + 1]; i++) {

			const Element *e = elements[permutation[i] / DOFs.size()];
			size_t dof = permutation[i] % DOFs.size();
			esglobal offset = 0;
			double duplicity = 1. / e->numberOfGlobalDomainsWithDOF(dof);

			for (auto c1 = e->clusters().begin(); c1 != e->clusters().end(); ++c1) {
				for (size_t d1 = 0; d1 < e->DOFCounter(*c1, dof); d1++) {

					for (auto c2 = c1; c2 != e->clusters().end(); ++c2) {
						for (size_t d2 = (*c1 == *c2 ? d1 + 1 : 0); d2 < e->DOFCounter(*c2, dof); d2++) {

							if (*c1 == config::env::MPIrank) {
								eslocal d = findDomain(e, d1, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + IJVMatrixIndexing);
								cols[t][d].push_back(e->DOFIndex(d, dof) + IJVMatrixIndexing);
								vals[t][d].push_back(1);
								dup[t][d].push_back(duplicity);
							}

							if (*c2 == config::env::MPIrank) {
								eslocal d = findDomain(e, d2, dof);
								rows[t][d].push_back(lambdasID[permutation[i]] + offset + IJVMatrixIndexing);
								cols[t][d].push_back(e->DOFIndex(d, dof) + IJVMatrixIndexing);
								vals[t][d].push_back(-1);
								dup[t][d].push_back(duplicity);
							}

							if (*c1 == config::env::MPIrank || *c2 == config::env::MPIrank) {
								cMap[t].push_back({ lambdasID[permutation[i]] + offset });
								if (*c1 == *c2) {
									cMap[t].back().push_back(*c1);
								} else if (*c1 == config::env::MPIrank) {
									cMap[t].back().push_back(*c1);
									cMap[t].back().push_back(*c2);
								} else {
									cMap[t].back().push_back(*c2);
									cMap[t].back().push_back(*c1);
								}
							}

							offset++;

						}
					}

				}

			}

		}
	}

	#pragma cilk grainsize = 1
	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t t = 0; t < threads; t++) {
			B1[p].I_row_indices.insert(B1[p].I_row_indices.end(), rows[t][p].begin(), rows[t][p].end());
			B1[p].J_col_indices.insert(B1[p].J_col_indices.end(), cols[t][p].begin(), cols[t][p].end());
			B1[p].V_values.insert(B1[p].V_values.end(), vals[t][p].begin(), vals[t][p].end());
			B1duplicity[p].insert(B1duplicity[p].end(), dup[t][p].begin(), dup[t][p].end());
		}
		B1[p].nnz = B1[p].I_row_indices.size();
		B1c[p].resize(B1[p].nnz, 0);
		for (size_t r = B1subdomainsMap[p].size(); r < B1[p].nnz; r++) {
			B1subdomainsMap[p].push_back(B1[p].I_row_indices[r] - 1);
		}
	}

	for (size_t t = 0; t < threads; t++) {
		B1clustersMap.insert(B1clustersMap.end(), cMap[t].begin(), cMap[t].end());
	}

	ESINFO(DETAILS) << "Lambdas in B1: " << B1[0].rows;
}

#ifdef HAVE_MORTAR

#include "mortar.h"

void EqualityConstraints::insertMortarGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	std::vector<int> rows;
	std::vector<int> columns;
	std::vector<double> values;

	std::vector<std::vector<int> > masterElements;
	std::vector<Point_3D> masterCoordinates;
	std::vector<std::vector<int> > slaveElements;
	std::vector<Point_3D> slaveCoordinates;
	std::vector<int> nodes;

	for (size_t i = 0; i < elements.size(); i++) {
		if (elements[i]->settings().isSet(Property::NONMATCHING_ELEMENT)) {
			if (config::env::MPIrank) {
				masterElements.push_back(std::vector<int>());
				for (size_t n = 0; n < elements[i]->nodes(); n++) {
					masterElements.back().push_back(elements[i]->node(n));
					nodes.push_back(elements[i]->node(n));
				}
			} else {
				slaveElements.push_back(std::vector<int>());
				for (size_t n = 0; n < elements[i]->nodes(); n++) {
					slaveElements.back().push_back(elements[i]->node(n));
					nodes.push_back(elements[i]->node(n));
				}
			}

		}
	}

	if (!masterElements.size() && !slaveElements.size()) {
		// no MORTAR interface founded
		return;
	}

	std::sort(nodes.begin(), nodes.end());
	Esutils::removeDuplicity(nodes);

	for (size_t n = 0; n < nodes.size(); n++) {
		if (config::env::MPIrank) {
			masterCoordinates.push_back(Point_3D());
			masterCoordinates.back().x = _mesh.coordinates()[nodes[n]].x;
			masterCoordinates.back().y = _mesh.coordinates()[nodes[n]].y;
			masterCoordinates.back().z = _mesh.coordinates()[nodes[n]].z;
		} else {
			slaveCoordinates.push_back(Point_3D());
			slaveCoordinates.back().x = _mesh.coordinates()[nodes[n]].x;
			slaveCoordinates.back().y = _mesh.coordinates()[nodes[n]].y;
			slaveCoordinates.back().z = _mesh.coordinates()[nodes[n]].z;
		}
	}

	std::vector<int> buffer;

	for (size_t e = 0; e < masterElements.size(); e++) {
		buffer.push_back(masterElements[e].size());
		for (size_t n = 0; n < masterElements[e].size(); n++) {
			masterElements[e][n] = std::lower_bound(nodes.begin(), nodes.end(), masterElements[e][n]) - nodes.begin();
			buffer.push_back(masterElements[e][n]);
		}
	}

	for (size_t e = 0; e < slaveElements.size(); e++) {
		for (size_t n = 0; n < slaveElements[e].size(); n++) {
			slaveElements[e][n] = std::lower_bound(nodes.begin(), nodes.end(), slaveElements[e][n]) - nodes.begin();
		}
	}

	if (config::env::MPIrank) {
		MPI_Send(buffer.data(), buffer.size() * sizeof(int), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(masterCoordinates.data(), masterCoordinates.size() * sizeof(Point_3D), MPI_BYTE, 0, 1, MPI_COMM_WORLD);
	} else {
		MPI_Status status;
		int size;

		// ELEMENTS
		MPI_Probe(1, 0, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		buffer.resize(size / sizeof(int));
		MPI_Recv(buffer.data(), size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		// COORDINATES
		MPI_Probe(1, 1, MPI_COMM_WORLD, &status);
		MPI_Get_count(&status, MPI_BYTE, &size);
		masterCoordinates.resize(size / sizeof(Point_3D));
		MPI_Recv(masterCoordinates.data(), size, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);

		for (size_t i = 0; i < buffer.size(); i++) {
			masterElements.push_back(std::vector<int>());
			for (size_t n = 0; n < buffer[i - n]; n++, i++) {
				masterElements.back().push_back(buffer[i + 1]);
			}
		}
	}

	if (!config::env::MPIrank && (masterElements.size() || slaveElements.size())) {
		computeMortarEqualityConstraints(rows, columns, values, masterElements, masterCoordinates, slaveElements, slaveCoordinates);
	}
}

#else

void EqualityConstraints::insertMortarGluingToB1(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	// TODO: improve!
	size_t cc = 0;
	for (size_t i = 0; i < elements.size(); i++) {
		if (elements[i]->settings().isSet(Property::NONMATCHING_ELEMENT)) {
			cc++;
		}
	}
	if (cc) {
		ESINFO(GLOBAL_ERROR) << "Gluing of non-matching grids is not supported. Link MORTAR library!";
	}
}

#endif


void EqualityConstraints::insertDomainGluingToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs)
{
	if (!elements.size()) {
		return;
	}

	size_t lambdas = B0[0].rows;

	for (size_t e = 0; e < elements.size(); e++) {
		for (size_t dof = 0; dof < DOFs.size(); dof++) {
			if (elements[e]->numberOfLocalDomainsWithDOF(dof) > 1) { // inner nodes are not glued
				const std::vector<eslocal> &DOFIndices = elements[e]->DOFsIndices();

				for (size_t d1 = 0, d2 = 1; d2 < elements[e]->domains().size(); d1++, d2++) {

					B0[elements[e]->domains()[d1]].I_row_indices.push_back(lambdas + IJVMatrixIndexing);
					B0[elements[e]->domains()[d1]].J_col_indices.push_back(DOFIndices[d1 * DOFs.size() + dof] + IJVMatrixIndexing);
					B0[elements[e]->domains()[d1]].V_values.push_back(1);

					B0[elements[e]->domains()[d2]].I_row_indices.push_back(lambdas + IJVMatrixIndexing);
					B0[elements[e]->domains()[d2]].J_col_indices.push_back(DOFIndices[d2 * DOFs.size() + dof] + IJVMatrixIndexing);
					B0[elements[e]->domains()[d2]].V_values.push_back(-1);

					lambdas++;
				}

			}
		}
	}


	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		B0[p].rows = lambdas;
		B0[p].nnz = B0[p].I_row_indices.size();

		B0subdomainsMap[p].reserve(B0[p].nnz);
		for (size_t i = B0subdomainsMap[p].size(); i < B0[p].nnz; i++) {
			B0subdomainsMap[p].push_back(B0[p].I_row_indices[i] - IJVMatrixIndexing);
		}
	}

	ESINFO(DETAILS) << "Average number of lambdas in B0 is " << Info::averageValue(lambdas);
}

void EqualityConstraints::insertKernelsToB0(const std::vector<Element*> &elements, const std::vector<Property> &DOFs, const std::vector<SparseMatrix> &kernel)
{
	std::vector<Element*> el(elements);

	std::sort(el.begin(), el.end(), [] (Element* e1, Element* e2) {
		if (e1->domains().size() != e2->domains().size()) {
			return e1->domains().size() < e2->domains().size();
		}
		return e1->domains() < e2->domains();
	});

	std::vector<size_t> part;
	part.push_back(std::lower_bound(el.begin(), el.end(), 2, [] (Element *e, int size) { return e->domains().size() < size; }) - el.begin());
	ESTEST(MANDATORY) << "There are not elements on the sub-domains interface." << ((elements.size() - part[0]) ? TEST_PASSED : TEST_FAILED);
	for (size_t i = part[0] + 1; i < el.size(); i++) {
		if (i && el[i - 1]->domains() != el[i]->domains()) {
			part.push_back(i);
		}
	}
	part.push_back(el.size());

	cilk_for (size_t p = 0; p < _mesh.parts(); p++) {
		for (size_t i = 0; i < part.size() - 1; i++) {
			const std::vector<eslocal> &domains = el[part[i]]->domains();
			int sign = domains[0] == p ? 1 : domains[1] == p ? -1 : 0;
			if (sign == 0) {
				continue;
			}


			std::vector<Element*> nodes;
			for (size_t e = part[i]; e < part[i + 1]; e++) {
				for (size_t n = 0; n < el[e]->nodes(); n++) {
					nodes.push_back(_mesh.nodes()[el[e]->node(n)]);
				}
			}
			std::sort(nodes.begin(), nodes.end());
			Esutils::removeDuplicity(nodes);

			for (size_t col = 0; col < kernel[domains[0]].cols; col++) {
				for (size_t n = 0; n < nodes.size(); n++) {
					for (size_t dof = 0; dof < DOFs.size(); dof++) {
						B0[p].I_row_indices.push_back(i * kernel[0].cols + col + IJVMatrixIndexing);
						B0[p].J_col_indices.push_back(nodes[n]->DOFIndex(p, dof) + IJVMatrixIndexing);
						B0[p].V_values.push_back(sign * kernel[domains[0]].dense_values[kernel[domains[0]].rows * col + nodes[n]->DOFIndex(domains[0], dof)]);
					}
				}
			}
		}


		B0[p].rows = kernel[0].cols * (part.size() - 1);
		B0[p].nnz = B0[p].I_row_indices.size();
		B0subdomainsMap[p].reserve(B0[p].nnz);
		for (size_t i = B0subdomainsMap[p].size(); i < B0[p].nnz; i++) {
			B0subdomainsMap[p].push_back(B0[p].I_row_indices[i] - IJVMatrixIndexing);
		}
	}
}



