
#ifndef SRC_OUTPUT_RESULT_VISUALIZATION_ENSIGHTWRITER_H_
#define SRC_OUTPUT_RESULT_VISUALIZATION_ENSIGHTWRITER_H_

#include <sstream>
#include <iomanip>

namespace espreso {

struct EnsightWriter {
	virtual void storeFormat(std::stringstream &os) const =0;
	virtual void storeDescriptionLine(std::stringstream &os, const std::string &description) const =0;

	virtual void storeInt(std::stringstream &os, int value) const =0;
	virtual void storeFloat(std::stringstream &os, float value) const =0;

	virtual void eIndex(std::stringstream &os, int value) const =0;
	virtual void eEnd(std::stringstream &os) const =0;


	virtual ~EnsightWriter() {};
};

struct EnsightBinaryWriter: public EnsightWriter {
	inline void storeFormat(std::stringstream &os) const
	{
		storeDescriptionLine(os, "C Binary");
	}

	inline void storeDescriptionLine(std::stringstream &os, const std::string &description) const
	{
		os << description << std::string(80 - description.size(), '\0');
	}

	inline void storeInt(std::stringstream &os, int value) const
	{
		os.write(reinterpret_cast<char*>(&value), sizeof(int));
	}

	inline void storeFloat(std::stringstream &os, float value) const
	{
		os.write(reinterpret_cast<char*>(&value), sizeof(float));
	}

	inline void eIndex(std::stringstream &os, int value) const
	{
		os.write(reinterpret_cast<char*>(&value), sizeof(int));
	}

	inline void eEnd(std::stringstream &os) const { }
};

struct EnsightASCIIWriter: public EnsightWriter {
	inline void storeFormat(std::stringstream &os) const { }

	inline void storeDescriptionLine(std::stringstream &os, const std::string &description) const
	{
		os << description << "\n";
	}

	inline void storeInt(std::stringstream &os, int value) const
	{
		os << std::setw(10) << value << "\n";
	}

	inline void storeFloat(std::stringstream &os, float value) const
	{
		os << value << "\n";
	}

	inline void eIndex(std::stringstream &os, int value) const
	{
		os << std::setw(10) << value;
	}

	inline void eEnd(std::stringstream &os) const
	{
		os << "\n";
	}
};

}



#endif /* SRC_OUTPUT_RESULT_VISUALIZATION_ENSIGHTWRITER_H_ */
