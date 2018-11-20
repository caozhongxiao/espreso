
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_MATRIXINDEX_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_MATRIXINDEX_H_

namespace espreso {

struct IJ {
	eslocal row, column;
};

bool operator==(const IJ &left, const IJ &right)
{
	return left.row == right.row && left.column == right.column;
}

bool operator!=(const IJ &left, const IJ &right)
{
	return !(left == right);
}

bool operator<(const IJ &left, const IJ &right)
{
	return left.row == right.row ? left.column < right.column : left.row < right.row;
}

}



#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_MATRIXINDEX_H_ */
