
#include "ecf.h"
#include "../configuration.hpp"

espreso::RBFTargetTransformationConfiguration::RBFTargetTransformationConfiguration(ECFConfiguration *ECFRoot)
: _ECFRoot(ECFRoot), translation(DIMENSION::D3, true), scaling(DIMENSION::D3, false)
{
	transformation = MORPHING_TRANSFORMATION::TRANSLATION;
	REGISTER(transformation, ECFMetaData()
		.setdescription({ "Transformation variant." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("FIXED").setdescription("Fixed position."))
		.addoption(ECFOption().setname("OFFSET").setdescription("OFFSET."))
		.addoption(ECFOption().setname("ROTATION").setdescription("ROTATION."))
		.addoption(ECFOption().setname("TRANSLATION").setdescription("Translation.")));

	REGISTER(offset, ECFMetaData()
		.setdescription({ "Offset size." })
		.setdatatype({ ECFDataType::EXPRESSION })
		.setboundaryconditionvariables());

	REGISTER(translation, ECFMetaData()
		.setdescription({ "Translation vector." }));

	REGISTER(scaling, ECFMetaData()
		.setdescription({ "Scale vector." }));

	REGISTER(coordinate_system, ECFMetaData()
		.setdescription({ "Configuration of coordinate system." }));

	REGISTER(overriding, ECFMetaData()
		.setdescription({ "Turn morphing target override on/off." })
		.setdatatype({ ECFDataType::BOOL }));

	ECFRoot->getParameter(&ECFRoot->physics)->addListener(ECFParameter::Event::VALUE_SET, [&] () {
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

espreso::RBFTargetConfiguration::RBFTargetConfiguration(ECFConfiguration *ECFRoot)
{
	solver = MORPHING_RBF_SOLVER::DIRECT;
	REGISTER(solver, ECFMetaData()
		.setdescription({ "Mesh Morphing solver." })
		.setdatatype({ ECFDataType::OPTION })
		.addoption(ECFOption().setname("ITERATIVE").setdescription("Iterative"))
		.addoption(ECFOption().setname("DENSE").setdescription("Dense")));

	solver_precision = 1e-5;
	REGISTER(solver_precision, ECFMetaData()
		.setdescription({ "Solver requested precision." })
		.setdatatype({ ECFDataType::FLOAT }));

	function.value = "R";
	REGISTER(function, ECFMetaData()
		.setdescription({ "Radial basis function." })
		.setdatatype({ ECFDataType::EXPRESSION })
		.setvariables({ "R" }));

	REGISTER(targets , ECFMetaData()
		.setdescription({ "Morphed region name.", "Target configuration." })
		.setdatatype({ ECFDataType::REGION })
		.setpattern({ "REGION" }),
		ECFRoot);
}

espreso::MeshMorphing::MeshMorphing(ECFConfiguration *ECFRoot)
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


