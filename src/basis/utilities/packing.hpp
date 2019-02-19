
#include "packing.h"

namespace espreso {
namespace utils {

template<>
inline size_t packedSize(const std::string &data)
{
	return sizeof(size_t) + data.size();
}

template<typename Ttype>
inline size_t packedSize(const Ttype &data)
{
	return sizeof(Ttype);
}

template<>
inline size_t packedSize(const std::vector<std::string> &data)
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		size += packedSize(data[i]);
	}
	return size;
}

template<typename Ttype>
inline size_t packedSize(const std::vector<Ttype> &data)
{
	size_t size = sizeof(size_t);
	for (size_t i = 0; i < data.size(); i++) {
		size += packedSize(data[i]);
	}
	return size;
}

template<>
inline void pack(const std::string &data, char* &p)
{
	size_t size = data.size();
	std::memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	std::memcpy(p, data.data(), data.size());
	p += data.size();
}

template<typename Ttype>
inline void pack(const Ttype &data, char* &p)
{
	std::memcpy(p, &data, packedSize(data));
	p += packedSize(data);
}

template<>
inline void pack(const std::vector<std::string> &data, char* &p)
{
	size_t size = data.size();
	std::memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	for (size_t i = 0; i < data.size(); i++) {
		pack(data[i], p);
	}
}

template<typename Ttype>
inline void pack(const std::vector<Ttype> &data, char* &p)
{
	size_t size = data.size();

	std::memcpy(p, &size, packedSize(size));
	p += packedSize(size);

	for (size_t i = 0; i < data.size(); i++) {
		pack(data[i], p);
	}
}

template<>
inline void unpack(std::string &data, const char* &p)
{
	size_t size;
	std::memcpy(&size, p, packedSize(size));
	p += packedSize(size);
	data = std::string(p, size);
	p += size;
}

template<typename Ttype>
inline void unpack(Ttype &data, const char* &p)
{
	std::memcpy(&data, p, packedSize(data));
	p += packedSize(data);
}

template<>
inline void unpack(std::vector<std::string> &data, const char* &p)
{
	size_t size;
	std::memcpy(&size, p, packedSize(size));
	p += packedSize(size);

	if (size) {
		data.resize(size);
		for (size_t i = 0; i < data.size(); i++) {
			unpack(data[i], p);
		}
	}
}

template<typename Ttype>
inline void unpack(std::vector<Ttype> &data, const char* &p)
{
	size_t size;
	std::memcpy(&size, p, packedSize(size));
	p += packedSize(size);

	data.resize(size);
	for (size_t i = 0; i < data.size(); i++) {
		unpack(data[i], p);
	}
}

}
}
