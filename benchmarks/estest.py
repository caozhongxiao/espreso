
import shutil, os, subprocess, copy

try:
    import requests, git
    snailwatch = True
except ImportError:
    snailwatch = False

class ESPRESOTest:

    root = os.path.dirname(os.path.dirname(__file__))
    espreso = os.path.join(root, "build", "espreso")
    mpirun = [ "mpirun", "--map-by", "socket:pe=2", "--bind-to", "core", "-n" ]

    env = dict(os.environ)
    env["MKL_NUM_THREADS"] = "1"
    env["CILK_NWORKERS"] = "1"
    env["OMP_NUM_THREADS"] = "4"
    env["SOLVER_NUM_THREADS"] = "4"
    env["PAR_NUM_THREADS"] = "4"
    env["OMP_PROC_BIND"] = "TRUE"

    path = ""
    ecf = "espreso.ecf"
    processes = 4
    args = []
    store_results = False

    _program = []

    @staticmethod
    def has_snailwatch():
        return snailwatch and "SNAILWATCH_URL" in os.environ and "SNAILWATCH_TOKEN"in os.environ

    @staticmethod
    def raise_error(error):
        raise Exception("\n {3} \n\nPath: {2}\nProgram: {1}\n{0}\n\n {3} \n".format(
                error, " ".join(ESPRESOTest._program), ESPRESOTest.path, "#" * 80))

    @staticmethod
    def clean():
        shutil.rmtree(os.path.join(ESPRESOTest.path, "results"), ignore_errors=True)

    @staticmethod
    def run_program(program):
        ESPRESOTest._program = copy.deepcopy(program)
        if not ESPRESOTest.store_results:
            program.append("--OUTPUT::RESULTS_STORE_FREQUENCY=NEVER")
            if ESPRESOTest.has_snailwatch():
                program.append("-m")

        return subprocess.Popen(program,
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                cwd=ESPRESOTest.path,
                                env=ESPRESOTest.env).communicate()

    @staticmethod
    def run():
        program = copy.deepcopy(ESPRESOTest.mpirun)
        program.append(str(ESPRESOTest.processes))
        program.append(ESPRESOTest.espreso)
        program.extend([ "-c", os.path.join(ESPRESOTest.path, ESPRESOTest.ecf) ])
        program.extend(map(str, ESPRESOTest.args))

        output, error = ESPRESOTest.run_program(program)
        if error != "":
            ESPRESOTest.raise_error(error)

    @staticmethod
    def compare(preset):
        def compare(v1, v2):
            if v1 != v2 and abs(float(v1) - float(v2)) > 1e-4 and abs((float(v1) - float(v2)) / float(v1)) > 1e-3:
                error = "various monitored results:\n"
                error += "  region={0}".format(table1[0][column])
                error += ", property={0}".format(table1[1][column])
                error += ", stat={0}".format(table1[2][column])
                error += ", step={0}".format(table1[row][0])
                error += ", substep={0}\n".format(table1[row][1])
                error += "  preset={0} != current={1}\n".format(v1, v2)
                ESPRESOTest.raise_error(error)

        def create_table(emr):
            return [ [ value.strip() for value in line.split(";") ] for line in open(emr, "r").readlines() if len(line.strip()) ]

        table1 = create_table(os.path.join(ESPRESOTest.path, preset))
        table2 = create_table(os.path.join(ESPRESOTest.path, "results", "last", ESPRESOTest.ecf.replace(".ecf", ".emr")))

        if len(table1) != len(table2):
            ESPRESOTest.raise_error("various time steps")
        if len(table1[0]) != len(table2[1]):
            ESPRESOTest.raise_error("various monitored properties")

        for row, (row1, row2) in enumerate(zip(table1, table2)):
            for column, (value1, value2) in enumerate(zip(row1, row2)):
                compare(value1, value2)

    @staticmethod
    def report(timereport):
        if not ESPRESOTest.has_snailwatch():
            return

        header = {
            "Content-Type": "application/json",
            "Authorization": os.environ['SNAILWATCH_TOKEN']
        }

        path = os.path.relpath(ESPRESOTest.path, os.path.dirname(__file__))

        benchmark = "-".join(path.split("/") + map(str, ESPRESOTest.args))
        commit = git.Repo(search_parent_directories=True).head.object.hexsha
        processes = str(ESPRESOTest.processes)
        threads = ESPRESOTest.env["OMP_NUM_THREADS"]

        log = os.path.join(ESPRESOTest.path, "results", "last", ESPRESOTest.ecf.replace(".ecf", ".log"))

        results = { }
        regions = ["Mesh preprocessing timing- Total", "Physics solver timing- Total"]
        with open(log, 'r') as file:
            for line in file:
                for region in regions:
                    if line.startswith(region):
                        results[region.replace(" ", "_")] = { "type": "time", "value": line.split("avg.:")[1].split()[0] }

        response = requests.post(
            os.environ['SNAILWATCH_URL'],
            headers=header,
            json={
                "benchmark": benchmark,
                "environment": { "commit": commit, "processes": processes, "threads": threads },
                "result": results
            })

        if response.status_code != 201:
            ESPRESOTest.raise_error("Cannot push to snailwatch")

