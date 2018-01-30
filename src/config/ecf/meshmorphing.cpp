
#include "ecf.h"
#include "../configuration.hpp"

espreso::RBFTargetTransformation::RBFTargetTransformation(ECFConfiguration *ECFRoot)
: _ECFRoot(ECFRoot), translation(DIMENSION::D3, true)
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
		.setdatatype({ ECFDataType::EXPRESSION }));

	REGISTER(translation, ECFMetaData()
		.setdescription({ "Translation vector." }));

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

espreso::RBFTarget::RBFTarget(ECFConfiguration *ECFRoot)
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
		.setdescription({ "x." })
		.setdatatype({ ECFDataType::STRING }));

	REGISTER(targets , ECFMetaData()
		.setdescription({ "a.", "b." })
		.setdatatype({ ECFDataType::STRING })
		.setpattern({ "MY_TARGET" }),
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
		.setdescription({ "Region.", "Target." })
		.setdatatype({ ECFDataType::STRING })
		.setpattern({ "RBF morphing configuration." }),
		ECFRoot);
}


