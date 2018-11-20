
#include "composer.h"

#include "../controlers/controler.h"

using namespace espreso;

Composer::Composer(Mesh &mesh, Step &step, Instance &instance, Controler &controler)
: _mesh(mesh), _step(step), _instance(instance), _controler(controler)
{

}

void Composer::initData()
{
	_controler.initData();
}

