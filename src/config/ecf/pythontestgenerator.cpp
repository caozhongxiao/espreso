
#include "pythontestgenerator.h"
#include "../configuration.hpp"

using namespace espreso;

PythonTestGeneratorTable::PythonTestGeneratorTable()
{
	REGISTER(columns, ECFMetaData()
			.setdescription({ "Column values." })
			.setdatatype({ ECFDataType::STRING }));
	REGISTER(rows, ECFMetaData()
			.setdescription({ "Row values." })
			.setdatatype({ ECFDataType::STRING }));
	REGISTER(value, ECFMetaData()
			.setdescription({ "Mesured value." })
			.setdatatype({ ECFDataType::STRING }));
}

PythonTestGenerator::PythonTestGenerator()
{
	output = "generatedtests";
	REGISTER(output, ECFMetaData()
			.setdescription({ "Generated tests output directory." })
			.setdatatype({ ECFDataType::STRING }));

	measure_repetition = 2;
	REGISTER(measure_repetition, ECFMetaData()
			.setdescription({ "Number of test repetition. The result table average the values." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER }));

	gather_level = 1;
	REGISTER(gather_level, ECFMetaData()
			.setdescription({ "Higher levels will be executed together." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER }));

	REGISTER(levels, ECFMetaData()
			.setdescription({ "Python list generator. List has to contain tuple (value, name).", "Level that is available in expression via 'L1' where 1 is level index." })
			.setdatatype({ ECFDataType::POSITIVE_INTEGER, ECFDataType::STRING })
			.setpattern({ "1", "[1, 2, 3, 4]" }));
	REGISTER(variables, ECFMetaData()
			.setdescription({ "List of variables available in expressions.", "Variable available in expressions." })
			.setdatatype({ ECFDataType::STRING, ECFDataType::STRING })
			.setpattern({ "NPROCS", "L1 * L1 * L1" }));

	REGISTER(env, ECFMetaData()
			.setdescription({ "Environment settings." })
			.setdatatype({ ECFDataType::STRING }));
        REGISTER(warmup, ECFMetaData()
                        .setdescription({ "Run before testing." })
                        .setdatatype({ ECFDataType::STRING }));
	REGISTER(run, ECFMetaData()
			.setdescription({ "Run command." })
			.setdatatype({ ECFDataType::STRING }));
	REGISTER(exe, ECFMetaData()
			.setdescription({ "Execution command." })
			.setdatatype({ ECFDataType::STRING }));

	addSeparator();

	REGISTER(args, ECFMetaData()
			.setdescription({ "Expressions for default values.", "Index of ARG" })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER, ECFDataType::STRING })
			.setpattern({ "0", "L1" }));

	REGISTER(tables, ECFMetaData()
			.setdescription({ "Expressions for default values.", "Table description." })
			.setdatatype({ ECFDataType::NONNEGATIVE_INTEGER })
			.setpattern({ "RESULTS" }));
}


