
#ifndef SRC_MESH_TRANSFORMATION_TFLAGS_H_
#define SRC_MESH_TRANSFORMATION_TFLAGS_H_

#define ALLOWS_OPERATORS(enum_name) \
inline constexpr enum_name  operator& (enum_name  d1, enum_name  d2) { return enum_name( static_cast<int>(d1) & static_cast<int>(d2)); } \
inline constexpr enum_name  operator| (enum_name  d1, enum_name  d2) { return enum_name( static_cast<int>(d1) | static_cast<int>(d2)); } \
inline constexpr enum_name  operator^ (enum_name  d1, enum_name  d2) { return enum_name( static_cast<int>(d1) ^ static_cast<int>(d2)); } \
inline constexpr enum_name  operator~ (enum_name  d1)                { return enum_name(~static_cast<int>(d1)                       ); } \
\
inline const     enum_name& operator&=(enum_name &d1, enum_name &d2) { return d1 = d1 & d2; } \
inline const     enum_name& operator|=(enum_name &d1, enum_name &d2) { return d1 = d1 | d2; } \
inline const     enum_name& operator^=(enum_name &d1, enum_name &d2) { return d1 = d1 ^ d2; }

namespace espreso {

struct TFlags {

	enum SEPARATE: int {
		MATERIALS,
		DEGREES_OF_FREEDOM
	};

	enum class ELEVEL: int {
		ELEMENT,
		FACE,
		EDGE,
		NODE
	};
};

ALLOWS_OPERATORS(TFlags::SEPARATE)

}



#endif /* SRC_MESH_TRANSFORMATION_TFLAGS_H_ */
