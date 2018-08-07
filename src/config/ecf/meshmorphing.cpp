
#include "../configuration.hpp"
#include "root.h"

espreso::RBFTargetTransformationConfiguration::RBFTargetTransformationConfiguration(ECFRoot *ECFRoot)
: _ECFRoot(ECFRoot),
  offset(ECFMetaData::getboundaryconditionvariables()),
  translation(DIMENSION::D3, ECFMetaData::getboundaryconditionvariables(), "0"),
  scaling(DIMENSION::D3, { "TIME" }, "100")
{
	transformation = MORPHING_TRANSFORMATION::TRANSLATION;
	REGISTER(transformation, ECFMetaData()
		.setdescription({ "Transformation variant." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("FIXED").setdescription("Fixed position."))
		.addoption(ECFOption().setname("OFFSET").setdescription("OFFSET."))
		.addoption(ECFOption().setname("SCALING").setdescription("Scaling."))
		.addoption(ECFOption().setname("TRANSLATION").setdescription("Translation."))
		.addoption(ECFOption().setname("ROTATION").setdescription("ROTATION."))
		);

	addSeparator();

	REGISTER(offset, ECFMetaData()
		.setdescription({ "Offset size." })
		.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(scaling, ECFMetaData()
		.setdescription({ "Scale vector." }));

	REGISTER(translation, ECFMetaData()
		.setdescription({ "Translation vector." }));

	addSeparator();

	REGISTER(coordinate_system, ECFMetaData()
		.setdescription({ "Configuration of coordinate system." }));

	REGISTER(override, ECFMetaData()
		.setdescription({ "Turn morphing target override on/off." })
		.setdatatype({ ECFDataType::BOOL }));

	ECFRoot->getParameter(&ECFRoot->physics)->addListener(ECFParameter::Event::VALUE_SET, [&] (const std::string &value) {
		switch (ECFRoot->physics) {
		case PHYSICS::HEAT_TRANSFER_2D:
		case PHYSICS::STRUCTURAL_MECHANICS_2D:
			translation.dimension = DIMENSION::D2;
			break;
		case PHYSICS::HEAT_TRANSFER_3D:
		case PHYSICS::STRUCTURAL_MECHANICS_3D:
			translation.dimension = DIMENSION::D3;
			break;
		default:
			ESINFO(GLOBAL_ERROR) << "ESPREOS internal error: unknown physics while set RBFTargetTransformation.";
		}
	});
}
std::ostream& espreso::operator<<(std::ostream& os, const espreso::RBFTargetTransformationConfiguration &t)
{
	os<<"\t\tTRANSFORMATION\t\t:";
	switch (t.transformation) {
			case MORPHING_TRANSFORMATION::FIXED: {
				os<<" FIXED\n";
			} break;

			case MORPHING_TRANSFORMATION::TRANSLATION: {
				os<<" TRANSLATION\n";
				os<<"\t\t\tX\t\t: "<<t.translation.x.value<<"\n";
				os<<"\t\t\tY\t\t: "<<t.translation.y.value<<"\n";
				os<<"\t\t\tZ\t\t: "<<t.translation.z.value<<"\n";
			} break;

			case MORPHING_TRANSFORMATION::OFFSET: {
				os<<" OFFSET\n";
			} break;

			case MORPHING_TRANSFORMATION::SCALING: {
				os<<" SCALING\n";
				os<<"\t\t\tSCALING X\t: "<<t.scaling.x.value<<"\n";
				os<<"\t\t\tSCALING Y\t: "<<t.scaling.y.value<<"\n";
				os<<"\t\t\tSCALING Z\t: "<<t.scaling.z.value<<"\n";
				os<<"\n";
				os<<"\t\t\tCENTER X\t: "<<t.coordinate_system.center.x.value<<"\n";
				os<<"\t\t\tCENTER Y\t: "<<t.coordinate_system.center.y.value<<"\n";
				os<<"\t\t\tCENTER Z\t: "<<t.coordinate_system.center.z.value<<"\n";

			} break;

			case MORPHING_TRANSFORMATION::ROTATION: {
				os<<" ROTATION\n";
				os<<"\t\t\tROTATION X\t: "<<t.coordinate_system.rotation.x.value<<"\n";
				os<<"\t\t\tROTATION Y\t: "<<t.coordinate_system.rotation.y.value<<"\n";
				os<<"\t\t\tROTATION Z\t: "<<t.coordinate_system.rotation.z.value<<"\n";
				os<<"\n";
				os<<"\t\t\tCENTER X\t: "<<t.coordinate_system.center.x.value<<"\n";
				os<<"\t\t\tCENTER Y\t: "<<t.coordinate_system.center.y.value<<"\n";
				os<<"\t\t\tCENTER Z\t: "<<t.coordinate_system.center.z.value<<"\n";

			}break;
		 	 default:
		 		 ESINFO(GLOBAL_ERROR) << "ESPRESO internal error: implement mesh morphing tranformation.";
		}

	return os;
}


espreso::RBFTargetConfiguration::RBFTargetConfiguration(ECFRoot *ECFRoot)
: function({ "R" }, "R"),
  external_ffd(ECFRoot)
{
	solver = MORPHING_RBF_SOLVER::DIRECT;
	REGISTER(solver, ECFMetaData()
		.setdescription({ "Mesh Morphing solver." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("ITERATIVE").setdescription("Iterative"))
		.addoption(ECFOption().setname("DIRECT").setdescription("Direct")));

	solver_precision = 1e-5;
	REGISTER(solver_precision, ECFMetaData()
		.setdescription({ "Solver requested precision." })
		.setdatatype({ ECFDataType::FLOAT }));

	solver_max_iter = 600;
		REGISTER(solver_max_iter, ECFMetaData()
			.setdescription({ "Solver requested maximum number of iterations." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	REGISTER(function, ECFMetaData()
		.setdescription({ "Radial basis function." })
		.setdatatype({ ECFDataType::EXPRESSION }));

	addSeparator();

	REGISTER(target, ECFMetaData()
		.setdescription({ "Set of morphed elements." })
		.setdatatype({ ECFDataType::BOUNDARY_REGION }));

	REGISTER(morphers , ECFMetaData()
		.setdescription({ "Morphed region name.", "Target configuration." })
		.setdatatype({ ECFDataType::BOUNDARY_REGION })
		.setpattern({ "REGION" }),
		ECFRoot);

	REGISTER(external_ffd, ECFMetaData()
			.setdescription({ "Configurations defined by an external FFD file." }));
}

std::ostream& espreso::operator<<(std::ostream& os, const espreso::RBFTargetConfiguration &t)
{
	os<<"\tTARGET REGION\t\t\t: "<<t.target<<"\n";
	os<<"\tSOLVER\t\t\t\t: ";
	if (t.solver==espreso::MORPHING_RBF_SOLVER::DIRECT) os<<"DIRECT\n";
	else {
		os<<"ITERATIVE (prec. "<<t.solver_precision<<", max. iter. "<<t.solver_max_iter<<")\n";
	}
	os<<"\tRADIAL BASIS FUNCTION\t\t: "<<t.function.value<<"\n";

	return os;
}

espreso::ExternalFFDConfiguration::ExternalFFDConfiguration(ECFRoot *ECFRoot)
{
	path = "";
	REGISTER(path, ECFMetaData()
			.setdescription({ "Set path to external configuration file with regions." })
			.setdatatype({ ECFDataType::STRING }));

	REGISTER(morphers , ECFMetaData()
			.setdescription({ "Morphed region name.", "Target configuration." })
			.setdatatype({ ECFDataType::BOUNDARY_REGION })
			.setpattern({ "REGION" }),
			ECFRoot);
}

espreso::MeshMorphing::MeshMorphing(ECFRoot *ECFRoot)
{
	type = MORPHING_TYPE::NONE;
	REGISTER(type, ECFMetaData()
		.setdescription({ "Mesh Morphing type." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("NONE").setdescription("No Mesh Morphing."))
		.addoption(ECFOption().setname("RBF").setdescription("RBF (Radial Base Function) Mesh Morphing.")));

	REGISTER(rbf, ECFMetaData()
		.setdescription({ "Morphing name.", "Named RBF configuration." })
		.setdatatype({ ECFDataType::STRING })
		.setpattern({ "MORPHING_NAME" }),
		ECFRoot);
}


