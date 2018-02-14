
#include "../insituwrapper.h"

#include "../../../basis/logging/logging.h"

using namespace espreso;

InSituWrapper::InSituWrapper()
{
	ESINFO(ALWAYS_ON_ROOT) << Info::TextColor::YELLOW << "ESPRESO was built without any in situ wrapper. Set in situ wrapper and re-configure ESPRESO.";
}

void InSituWrapper::update()
{

}

InSituWrapper::~InSituWrapper()
{

}



