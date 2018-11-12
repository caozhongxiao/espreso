
#ifndef APP_FACTORY_FACTORY_H_
#define APP_FACTORY_FACTORY_H_

#include <vector>
#include <string>
#include <map>

namespace espreso {

class Mesh;
struct Instance;
struct Step;
struct Physics;
class LinearSolver;
class Composer;
class TimeStepSolver;
class LoadStepSolver;
class ResultStore;

struct ECFRoot;
struct LoadStepConfiguration;

class FactoryLoader {

	friend class APITestESPRESODataProvider;
public:
	virtual ~FactoryLoader();

	void preprocessMesh();
	virtual size_t loadSteps() const =0;

	virtual LoadStepSolver* getLoadStepSolver(size_t step, Mesh *mesh, ResultStore *store) =0;

	template<class TLoadStepSettings>
	const TLoadStepSettings& getLoadStepsSettings(size_t step, const std::map<size_t, TLoadStepSettings> &setting) const
	{
		if (setting.find(step + 1) == setting.end()) {
			printError("Missing setting for LOAD STEP " + std::to_string(step + 1));
		}
		return setting.find(step + 1)->second;
	}

	LinearSolver* getLinearSolver(const LoadStepConfiguration &settings, Instance *instance) const;

protected:
	void printError(const std::string &error) const;

	std::vector<Instance*> _instances;
	std::vector<Physics*> _physics;
	std::vector<LinearSolver*> _linearSolvers;
	std::vector<Composer*> _composers;
	std::vector<TimeStepSolver*> _timeStepSolvers;
	std::vector<LoadStepSolver*> _loadStepSolvers;
};


class Factory {

	friend class ESPRESO;
	friend class APITestESPRESODataProvider;

protected:
	Factory(const ECFRoot &configuration, Mesh &mesh, ResultStore &store);
	~Factory();

	void solve();

	FactoryLoader* createFactoryLoader(const ECFRoot &configuration);

	Mesh *_mesh;
	ResultStore *_store;
	Step *_step;

	FactoryLoader *_loader;
	std::vector<LoadStepSolver*> _loadSteps;
};

class ESPRESO {

public:
	static void run(int *argc, char ***argv);
};


}



#endif /* APP_FACTORY_FACTORY_H_ */
