
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include <vector>
#include <string>
#include <map>

namespace espreso {

class Mesh;
struct DataHolder;
struct Physics;
class LinearSolver;
class ProviderOOLLDD;
class Composer;
class Controler;
class TimeStepSolver;
class LoadStepSolver;

struct ECFRoot;
struct LoadStepConfiguration;

class FactoryLoader {

	friend class APITestESPRESODataProvider;
public:
	virtual ~FactoryLoader();

	void preprocessMesh();
	virtual size_t loadSteps() const =0;

	virtual LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh) =0;

	template<class TLoadStepSettings>
	TLoadStepSettings& getLoadStepsSettings(size_t step, std::map<size_t, TLoadStepSettings> &setting) const
	{
		if (setting.find(step + 1) == setting.end()) {
			printError("Missing setting for LOAD STEP " + std::to_string(step + 1));
		}
		return setting.find(step + 1)->second;
	}

	LinearSolver* getLinearSolver(LoadStepConfiguration &settings, DataHolder *instance) const;

protected:
	void printError(const std::string &error) const;

	std::vector<DataHolder*> _instances;
	std::vector<LinearSolver*> _linearSolvers;
	std::vector<ProviderOOLLDD*> _provider;
	std::vector<Composer*> _composer;
	std::vector<TimeStepSolver*> _timeStepSolvers;
	std::vector<LoadStepSolver*> _loadStepSolvers;
};


class Factory {

	friend class ESPRESO;
	friend class APITestESPRESODataProvider;

protected:
	Factory(ECFRoot &configuration, Mesh &mesh);
	~Factory();

	void solve();

	FactoryLoader* createFactoryLoader(ECFRoot &configuration);

	Mesh *_mesh;

	FactoryLoader *_loader;
	std::vector<LoadStepSolver*> _loadSteps;
};

class ESPRESO {

public:
	static void run(int *argc, char ***argv);
};


}



#endif /* APP_FACTORY_FACTORY_H_ */
