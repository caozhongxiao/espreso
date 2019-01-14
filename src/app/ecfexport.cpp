#include "mpi.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "config/ecf/root.h"


using namespace espreso;

std::map<ECFDataType, std::string> ECFDatatypes;

std::string spaces(const int indent = 0)
{
    std::stringstream ss;
    for (int i = 0; i < indent; i++) ss << " ";
    return ss.str();
}

void printKeyValuePair(const std::string& key, const std::string& val, const int indent = 0)
{
    for (int i = 0; i < indent; i++) printf(" ");
    std::cout << key << ": \"" << val << "\"";
}

void printKeyArrayPair(const std::string& key, const std::vector<std::string>& arr, const int indent = 0, const char* itemEnclosed = "\"")
{
    for (int i = 0; i < indent; i++) printf(" ");
    std::cout << key << ": " << "[";
    
    std::vector<std::string> aux_arr;
    for (auto e = arr.begin(); e != arr.end(); ++e)
    {
        aux_arr.push_back(itemEnclosed);
        aux_arr.push_back(*e);
        aux_arr.push_back(itemEnclosed);
        aux_arr.push_back(", ");
    }
    if (aux_arr.size()) aux_arr.pop_back();
    for (auto e = aux_arr.begin(); e != aux_arr.end(); ++e)
    {
        std::cout << *e;
    }

    std::cout << "]";
}

void printKeyObjectPair(const std::string& key, 
    const std::function<void(int)>& printObjContent, 
    const int indent = 0)
{
    std::cout << spaces(indent);
    std::cout << key << ": " << "{" << std::endl;
    
    printObjContent(indent + 2);

    std::cout << std::endl << spaces(indent + 2) << "}" << std::endl;
}

void printMetaData(ECFMetaData *md, const int indent = 0)
{
    int ind = indent + 2;
    printKeyArrayPair("description", md->description, ind);
    std::cout << "," << std::endl;
    
    auto datatypesToStringVector = [&](const std::vector<ECFDataType>& dts) -> std::vector<std::string>{
        std::vector<std::string> ret;
        for (auto d = dts.begin(); d != dts.end(); ++d)
        {
            ret.push_back(ECFDatatypes[*d]);
        }
        return ret;
    };
    printKeyArrayPair("datatype", datatypesToStringVector(md->datatype), ind);
    
    std::cout << "," << std::endl;
    
    printKeyArrayPair("pattern", md->pattern, ind);

    std::cout << "," << std::endl;

    auto printOptions = [&] () {
        std::vector<std::string> ret;
        for (auto option = md->options.begin(); option != md->options.end(); ++option)
        {
            std::stringstream ss;
            ss 
                << "{ name: \""
                << option->name
                << "\", description: "
                << option->description
                << "\" }";
            ret.push_back(ss.str());
        }

        return ret;
    };
    printKeyArrayPair("options", printOptions(), ind, "");
    
    std::cout << "," << std::endl;
    
    printKeyArrayPair("variables", md->variables, ind);

    std::cout << "," << std::endl;
    
    printKeyValuePair("unit", md->unit.unit(), ind);
}


void printValue(ECFParameter* val, const int indent = 0)
{
    printKeyValuePair("value", val->getValue(), indent + 2);
    std::cout << "," << std::endl;
    auto printMetadata = [&](const int p_indent = 0) {
        printMetaData(&val->metadata, p_indent);
    };
    printKeyObjectPair("metadata", printMetadata, indent + 2);
}

void printParameter(ECFParameter* param, const int indent = 0)
{
    std::cout 
            << std::endl
            << spaces(indent) 
            << param->name
            << ": {" << std::endl;

    if (param->isValue())
    {
        printValue(param, indent);
    }
    if (param->isObject())
    {
        ECFObject* obj = static_cast<ECFObject*>(param);
        // ECFParameter* pattern = obj->getPattern();
        // if (pattern) 
        // {
        //     pattern->name = "pattern";
        //     printParameter(pattern, indent + 2);
        // }
        
        auto p = obj->parameters.begin();
        if (p != obj->parameters.end()) {
            printParameter(*p, indent + 2);
            p++;
        }
        for (; p != obj->parameters.end(); ++p)
        {
            std::cout << "," << std::endl;
            printParameter(*p, indent + 2);
        }
    }

    std::cout << spaces(indent) << "}";
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    ECFDatatypes[ECFDataType::BOOL] = "bool";
	ECFDatatypes[ECFDataType::STRING] = "string";
	ECFDatatypes[ECFDataType::INTEGER] = "integer";
	ECFDatatypes[ECFDataType::POSITIVE_INTEGER] = "positive_integer";
	ECFDatatypes[ECFDataType::NONNEGATIVE_INTEGER] = "nonnegative_integer";
	ECFDatatypes[ECFDataType::FLOAT] = "float";
	ECFDatatypes[ECFDataType::ENUM_FLAGS] = "enum_flags";
	ECFDatatypes[ECFDataType::OPTION] = "option";
	ECFDatatypes[ECFDataType::REGION] = "region";
	ECFDatatypes[ECFDataType::BOUNDARY_REGION] = "boundary_region";
	ECFDatatypes[ECFDataType::ELEMENTS_REGION] = "elements_region";
	ECFDatatypes[ECFDataType::MATERIAL] = "material";
	ECFDatatypes[ECFDataType::LOAD_STEP] = "load_step";
	ECFDatatypes[ECFDataType::EXPRESSION] = "expression";
	ECFDatatypes[ECFDataType::TENSOR] = "tensor";
	ECFDatatypes[ECFDataType::INTERVAL] = "interval";
	ECFDatatypes[ECFDataType::SPACE] = "space";
	ECFDatatypes[ECFDataType::SEPARATOR] = "separator";

    ECFRoot ecf("benchmarks/diffusion2D/steadystate/simple/espreso.ecf");
    if (rank == 0)
        printParameter(&ecf);
        // printParameter(&ecf.heat_transfer_3d.materials["material"]);

    ecf.forEachParameters([&] (const ECFParameter *parameter) {
		if (parameter->metadata.condition->isset()) {
			std::cout << parameter->name << ": ";
			ecf.forEachParameters([&] (const ECFParameter *conditionalParameter) {
				if (parameter->metadata.condition->match(conditionalParameter->data())) {
					std::cout << parameter->metadata.condition->compose(conditionalParameter);
				}
			}, false);
			std::cout << "\n";
		}
	}, false);

    MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}



