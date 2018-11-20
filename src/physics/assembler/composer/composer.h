
#ifndef SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_
#define SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_

namespace espreso {

class Mesh;
class Step;
class Instance;
class Controler;
enum Matrices: int;

class Composer {

public:
	Composer(Mesh &mesh, Step &step, Instance &instance, Controler &controler);

	virtual void initDOFs() = 0;
	virtual void buildPatterns() = 0;

	void initData();
	virtual void assemble(Matrices matrices) = 0;

	virtual void setDirichlet() = 0;
	virtual void synchronize() = 0;

	virtual ~Composer() {}

protected:
	Mesh &_mesh;
	Step &_step;
	Instance &_instance;
	Controler &_controler;
};

}


#endif /* SRC_PHYSICS_ASSEMBLER_COMPOSER_COMPOSER_H_ */
