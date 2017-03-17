
import subprocess
import os
import sys
import re
import shutil

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ESPRESO_ROOT = os.path.dirname(ESPRESO_TESTS)
EXAMPLES = os.path.join(ESPRESO_TESTS, "examples")

ENV = {
    "MKL_NUM_THREADS": "1",
    "CILK_NWORKERS": "1",

    "OMP_NUM_THREADS": "4",
    "SOLVER_NUM_THREADS": "4",
    "PAR_NUM_THREADS": "4",

    "PARDISOLICMESSAGE": "1",

    "OMP_PROC_BIND": "TRUE", # there is problem with linking CUDA
}


class TestCaseCreator:

    @staticmethod
    def create_test(testClass, function, name, *args, **kwargs):
        def test_method(self):
            function(self, *args, **kwargs)

        setattr(testClass, 'test_' + name, test_method)
        test_method.__name__ = 'test_' + name

    @staticmethod
    def iterate(function, *args):
        next = [ True for arg in args]
        iterators = [ Iterator(arg) for arg in args ]

        while reduce(lambda x, y: x or y, next):
            function(*[ it.get() for it in iterators ])
            for i in range(0, len(next)):
                next[i] = iterators[i].next()
                if next[i]:
                    break

    @staticmethod
    def select(*directories):
        argv = [ sys.argv[0] ]
        origin = []
        selection = []

        for directory in directories:
            origin.append([ directory ])
            selection.append([])
            for arg in sys.argv[1:]:
                if arg.startswith("-") or arg.startswith("--"):
                    argv.append(arg)
                    continue
                if os.path.commonprefix([ directory, os.path.join(ESPRESO_ROOT, arg) ]) == directory:
                    selection[-1].append(os.path.join(ESPRESO_ROOT, arg))

        if len(sys.argv) == len(argv):
            if len(origin) == 1:
                return origin[0]
            return tuple(origin)
        sys.argv = argv
        if len(selection) == 1:
            return selection[0]
        return tuple(selection)

    @staticmethod
    def gather(folder, ext, omit = "$^"):
        def skip(files):
            for file in files:
                if file.endswith(".skip"):
                    print "skip folder " + root
                    for line in open(os.path.join(root, file)):
                        print "--> " + line
                    return True
            return False

        # recurse tree a gather all examples
        omit = re.compile(omit)
        examples = []
        for root, subFolders, files in os.walk(folder):
            if not skip(files):
                for file in files:
                    if file.endswith(ext) and not omit.match(file):
                        examples.append(( root, file.replace(ext, '')))

        examples.sort()

        result = []
        for example in examples:
            name = os.path.relpath(example[0], folder).replace('/', '_')
            if example[1] != "espreso":
                name += "_" + example[1]
            result.append((name, example[0], example[1] + ext))
        return result

class Iterator:

    def __init__(self, items):
        self.items = items
        self.keys = items.keys()
        self.pointers = [ 0 for i in items ]

    def next(self):
        self.pointers[-1] += 1
        for i in xrange(1, len(self.pointers)):
            if self.pointers[-i] == len(self.items[self.keys[-i]]):
                self.pointers[-i] = 0
                self.pointers[-i - 1] += 1

        if self.pointers[0] == len(self.items[self.keys[0]]):
            self.pointers[0] = 0
            return False

        return True

    def __getitem__(self, i):
        return self.get()[i]

    def values(self):
        return self.get().values()

    def get(self):
        result = {}
        for i in xrange(0, len(self.pointers)):
            result[self.keys[i]] = self.items[self.keys[i]][self.pointers[i]]
        return result

class RunInfo:

    def __init__(self, output):
        self._iterations = 0
        self._precision = 1
        self._oscilated = False
        self._max_oscilation = 0

        norm = 2
        min = 2
        for line in output.split("\n"):
            if line.find("CONVERGENCE:") != -1:
                try:
                    self._iterations = int(line.split()[1])
                except ValueError:
                    continue
                if norm < float(line.split()[2]):
                    self._oscilated = True
                norm = float(line.split()[2])
                if min > norm:
                    min = norm
                if min < norm:
                    self._max_oscilation = norm / min

        self._precision = norm

    def iterations(self, min, max=None):
        if self._iterations < min:
            raise EspresoError("Unexpected number of iterations: {0}".format(self._iterations))
        if max is not None and self._iterations > max:
            raise EspresoError("Number of iterations is too high: {0}".format(self._iterations))

    def precision(self, precision):
        if self._precision > precision:
            raise EspresoError("Example not convergated to a requested precision ({0} > {1}".format(self._precision, precision))

    def oscilation(self, allowed, max=None):
        if not allowed and self._oscilated:
            raise EspresoError("Not allowed oscilation")
        if max is not None and self._max_oscilation > max:
            raise EspresoError("Oscilation {0} is higher than {1}".format(self._max_oscilation, max))

class EspresoError:

    program = []

    def __init__(self, error):
        raise Exception("{0}\nProgram: {1}".format(error, self.program))

class Espreso:

    def __init__(self, path=ESPRESO_ROOT, config={}):
        self.path = path
        self.mpirun = [ "mpirun" ]

        bashrc = os.path.join(os.path.expanduser("~"), ".bashrc")
        if os.path.isfile(bashrc):
            for line in open(bashrc, 'r'):
                if line.find("alias") != -1 and line.find("mpirun") != -1:
                    self.mpirun = map(lambda x: x.strip(" '\""), line.split("=")[1].split())

        if path != ESPRESO_ROOT:
            self.path = os.path.join(ESPRESO_TESTS, path)

        self.env = os.environ
        self.append_env("LD_LIBRARY_PATH", os.path.join(self.path, "libs/"))
        self.append_env("LD_LIBRARY_PATH", "/usr/local/cuda-7.5/lib64")
        for key, value in ENV.items():
            self.set_env(key, value)

        if path != ESPRESO_ROOT:
            self.clone()

    def append_env(self, key, value):
        self.env[key] = value + ":" + self.env[key]

    def set_env(self, key, value):
        self.env[key] = value


    def clone(self):
        self.remove()

        subprocess.call(["git", "clone", "git@code.it4i.cz:mec059/espreso.git", self.path, "-q"], cwd=ESPRESO_TESTS)
        subprocess.call(["mkdir", "libs"], cwd=self.path)

        origin = os.path.join(ESPRESO_ROOT, "libs/")
        target = os.path.join(self.path, "libs/")
        if os.path.isfile(origin + "libpardiso500-INTEL120-X86-64.so"):
            subprocess.call(["cp", origin + "libpardiso500-INTEL120-X86-64.so", target])
        if os.path.isfile(origin + "libifcore.a"):
            subprocess.call(["cp", origin + "libifcore.a", target])

    def remove(self):
        shutil.rmtree(self.path, ignore_errors=True)


    def waf(self, args=[], config={}):
        return self.run_program([os.path.join(self.path, "waf")], self.path, config, args)

    def install(self, config={}):
        def check(result, error, method):
            success = False
            for line in result.splitlines():
                if line.find("'" + method + "' finished successfully") != -1:
                    success = True

            if error != "":
                raise EspresoError(error)
            if success == False:
                raise EspresoError(result)

        result, error = self.waf(["configure"], config)
        check(result, error, "configure")
        result, error = self.waf(["install"])
        check(result, error, "install")

    def run_program(self, program, cwd="", config={}, args=[]):
        program += [ str(x) for x in args ]
        for key, value in config.items():
            program += [ "--{0}={1}".format(key, value) ]

        if cwd:
            cwd = os.path.join(ESPRESO_ROOT, cwd)
        EspresoError.program = program
        result = subprocess.Popen(program,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                  cwd=cwd or self.path,
                                  env=self.env)

        return result.communicate()


    def run(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes),"--map-by", "core", os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            raise EspresoError(error)
        if output != "":
            raise EspresoError(output)

    def valgrind(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes), "valgrind", "-q", "--leak-check=full", "--suppressions={0}/espreso.supp".format(self.path), os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            skip = False
            warningless = ""
            for line in error.split("\n"):
                tokens = line.split(" ")
                if len(tokens) == 1:
                    continue
                if tokens[1] == "Warning:":
                    skip = True
                    continue
                elif tokens[1] == "" and skip:
                    continue
                else:
                    skip = False
                warningless += line + "\n"
            if warningless:
                raise EspresoError("\n" + warningless)
        if output != "":
            raise EspresoError(output)

    def decompose(self, *args, **kwargs):
        program = [ os.path.join(self.path, "decomposer") ]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            raise EspresoError(error)
        return output

    def output(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes), os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error != "":
            raise EspresoError(error)

        return output

    def fail(self, processes, *args, **kwargs):
        program = self.mpirun + [ "-n", str(processes), os.path.join(self.path, "espreso")]

        output, error = self.run_program(program, *args, **kwargs)
        if error == "":
            raise EspresoError("Expected fail, but run was correct.")






