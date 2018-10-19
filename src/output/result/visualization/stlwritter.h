
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_STLWRITER_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_STLWRITER_H_

#include <sstream>
#include <iomanip>

namespace espreso {

struct STLWriter {
	virtual void storeHeader(std::stringstream &os, const std::string &name) const =0;
	virtual void storeFooter(std::stringstream &os, const std::string &name) const =0;

	virtual void storeSize(std::stringstream &os, int size) const=0;

	virtual void beginFace(std::stringstream &os, float x, float y, float z) const=0;
	virtual void addVertex(std::stringstream &os, float x, float y, float z) const=0;
	virtual void endFace(std::stringstream &os) const=0;

	virtual ~STLWriter() {};
};

struct STLBinaryWriter: public STLWriter {

	inline void storeHeader(std::stringstream &os, const std::string &description) const
	{
		os << description << std::string(80 - description.size(), '\0');
	}

	inline void storeFooter(std::stringstream &os, const std::string &name) const
	{

	}

	inline void storeSize(std::stringstream &os, int size) const
	{
		os.write(reinterpret_cast<char*>(&size), sizeof(int));
	}

	inline void beginFace(std::stringstream &os, float x, float y, float z) const
	{
		os.write(reinterpret_cast<char*>(&x), sizeof(float));
		os.write(reinterpret_cast<char*>(&y), sizeof(float));
		os.write(reinterpret_cast<char*>(&z), sizeof(float));
	}

	inline void addVertex(std::stringstream &os, float x, float y, float z) const
	{
		os.write(reinterpret_cast<char*>(&x), sizeof(float));
		os.write(reinterpret_cast<char*>(&y), sizeof(float));
		os.write(reinterpret_cast<char*>(&z), sizeof(float));
	}

	inline void endFace(std::stringstream &os) const
	{
		int attribute = 0;
		os.write(reinterpret_cast<char*>(&attribute), 2);
	}
};

struct STLASCIIWriter: public STLWriter {

	inline void storeHeader(std::stringstream &os, const std::string &name) const
	{
		os << "solid " << name << "\n";
	}

	inline void storeFooter(std::stringstream &os, const std::string &name) const
	{
		os << "endsolid " << name << "\n";
	}

	inline void storeSize(std::stringstream &os, int size) const
	{

	}

	inline void beginFace(std::stringstream &os, float x, float y, float z) const
	{
		os << " facet normal " << x << " " << y << " " << z << "\n";
		os << "  outer loop\n";
	}

	inline void addVertex(std::stringstream &os, float x, float y, float z) const
	{
		os << "   vertex " << x << " " << y << " " << z << "\n";
	}

	inline void endFace(std::stringstream &os) const
	{
		os << "  endloop\n";
		os << " endfacet\n";
	}
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_STLWRITER_H_ */
