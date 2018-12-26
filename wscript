
import commands
import sys
import os

reload(sys)
sys.setdefaultencoding('utf8')

from waflib.TaskGen import after_method, feature

@feature('cxx')
@after_method('apply_incpaths')
def insert_srcdir(self):
    path = os.path.join(self.bld.srcnode.abspath(), "src")
    self.env.prepend_value('INCPATHS', path)

def set_compiler(ctx):
    def trycompiler():
        try:
            ctx.find_program(ctx.options.mpicxx)
        except ctx.errors.ConfigurationError:
            return False
        return True

    if not trycompiler():
        if ctx.options.mpicxx == "mpiicpc":
            ctx.options.mpicxx = "mpic++"
        elif ctx.options.mpicxx == "mpic++":
            ctx.options.mpicxx = "mpiicpc"
        if not trycompiler():
            ctx.fatal("Cannot found MPI compiler. Set a correct one by 'mpicxx=' parameter.")

    if ctx.options.mpicxx == "mpiicpc":
        ctx.options.cxx = "icpc"
    if ctx.options.mpicxx == "mpic++":
        ctx.options.cxx = "g++"
    ctx.load(ctx.options.cxx)
    ctx.env.CXX = ctx.env.LINK_CXX =ctx.options.mpicxx

def set_openmp(ctx):
    if ctx.options.cxx == "icpc":
        ctx.env.append_unique("CXXFLAGS", [ "-qopenmp" ])
        ctx.env.append_unique("LINKFLAGS", [ "-qopenmp" ])
    if ctx.options.cxx == "g++":
        ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
        ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])

def set_metis(ctx):
    ctx.check_cxx(
        header_name="metis.h parmetis.h", lib=["metis", "parmetis"], uselib_store="METIS",
        define_name="",
        defines=["HAVE_METIS", "IDXTYPEWIDTH=" + str(ctx.options.intwidth), "REALTYPEWIDTH=64"],
        includes=ctx.options.parmetisroot + "/include",
        libpath=ctx.options.parmetisroot + "/lib",
        msg="Checking for library ParMETIS",
        errmsg="set 'parmetisroot'.")

def set_mkl(ctx):
    defines = [ "HAVE_MKL" ]
    libs = []
    if ctx.options.intwidth == 32:
        defines.append("MKL_INT=int")
        libs.append("mkl_intel_lp64")
    if ctx.options.intwidth == 64:
        defines.append("MKL_INT=long")
        libs.append("mkl_intel_ilp64")
    libs.append("mkl_core")
    if ctx.options.cxx == "icpc":
        libs.append("mkl_intel_thread")
    if ctx.options.cxx == "g++":
        libs.append("mkl_gnu_thread")

    ctx.check_cxx(
        header_name="mkl.h", lib=libs, uselib_store="MKL",
        define_name="", defines=defines,
        includes=ctx.options.mklroot + "/include",
        libpath=ctx.options.mklroot + "/lib",
        msg="Checking for library MKL",
        errmsg="set 'mklroot'.")

def try_hypre(ctx):
    ctx.check_cxx(
        header_name="HYPRE.h", lib="HYPRE", uselib_store="HYPRE",
        define_name="", defines="HAVE_HYPRE",
        includes=ctx.options.hypreroot + "/include",
        libpath=ctx.options.hypreroot + "/lib",
        mandatory=False,
        msg="Checking for library HYPRE",
        errmsg="set 'hypreroot' to use HYPRE solver.")

def try_catalyst(ctx):
    pass

def set_variables(ctx):
    if ctx.options.intwidth == 32:
        ctx.env.append_unique("DEFINES", [ "esint=int", "esint_mpi=MPI_INT" ])
    if ctx.options.intwidth == 64:
        ctx.env.append_unique("DEFINES", [ "esint=long", "esint_mpi=MPI_LONG" ])

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++11", "-Wall"])
    ctx.env.mode = ctx.options.mode
    if ctx.options.mode == "release":
        ctx.env.append_unique("CXXFLAGS", [ "-O3" ])
    if ctx.options.mode == "development":
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-g" ])
    if ctx.options.mode == "debug":
        ctx.env.append_unique("CXXFLAGS", [ "-O0", "-g" ])

def configure(ctx):
    set_compiler(ctx)
    set_openmp(ctx)
    set_mkl(ctx)
    set_metis(ctx)
    try_hypre(ctx)
    set_variables(ctx)

    ctx.msg("Setting compiler to", ctx.options.mpicxx)
    ctx.msg("Setting int width to", ctx.options.intwidth)
    ctx.msg("Setting build mode to", ctx.options.mode)
    ctx.msg("Setting solver to", ctx.options.solver)

fetisources= (
   "src/solver/generic/Domain.cpp",
   "src/solver/generic/SparseMatrix.cpp",
   "src/solver/generic/utils.cpp",
   "src/solver/generic/FETISolver.cpp",
   "src/solver/specific/cluster.cpp",
   "src/solver/specific/itersolver.cpp",
   "src/solver/specific/cpu/SparseSolverMKL.cpp",
   "src/solver/specific/cpu/clustercpu.cpp",
   "src/solver/specific/cpu/itersolvercpu.cpp",
   "src/solver/specific/cpu/DenseSolverMKL.cpp"
)

def build(ctx):

    ctx.objects(source=ctx.path.ant_glob('src/basis/**/*.cpp'),target="basis")
    ctx.objects(source=ctx.path.ant_glob('src/globals/**/*.cpp'),target="globals",defines="MODE="+ctx.env.mode.upper())
    ctx.objects(source=ctx.path.ant_glob('src/config/**/*.cpp'),target="config")
    ctx.objects(source=ctx.path.ant_glob('src/mesh/**/*.cpp'),target="mesh")
    ctx.objects(source=ctx.path.ant_glob('src/input/**/*.cpp'),target="input")
    ctx.objects(source=ctx.path.ant_glob('src/output/**/*.cpp'),target="output",includes="tools/ASYNC")
    ctx.objects(source=ctx.path.ant_glob('src/physics/**/*.cpp'),target="physics",defines=["SOLVER_MKL"])
    ctx.objects(source=ctx.path.ant_glob('src/linearsolver/**/*.cpp'),target="linearsolver")
    ctx.objects(source=ctx.path.ant_glob('src/wrappers/metis/**/*.cpp'),target="metis",use="METIS")
    ctx.objects(source=ctx.path.ant_glob('src/wrappers/math/**/*.cpp'),target="math",use="MKL")
    ctx.objects(source=ctx.path.ant_glob('src/wrappers/hypre/**/*.cpp'),target="hypre",use="HYPRE")
    ctx.objects(source=ctx.path.ant_glob('src/wrappers/catalyst/**/*.cpp'),target="catalyst",use="CATALYST")
    ctx.objects(source=fetisources,target="feti",defines=["SOLVER_MKL"])

    ctx.program(
        source="src/app/espreso.cpp",target="espreso",
        use="config basis globals mesh input output physics feti linearsolver math hypre catalyst metis",
    )

def options(opt):
    espreso = opt.add_option_group("ESPRESO library options")

    espreso.add_option("--mpicxx",
        action="store",
        type="string",
        default="mpiicpc",
        help="MPI compiler used for building of the library [default: %default]")

    espreso.add_option("--cxx",
        action="store",
        choices=["icpc", "g++"],
        metavar="icpc,g++",
        default="icpc",
        help="C++ compiler (set it in the case of non-standard 'mpicxx' settings) [default: %default]")

    espreso.add_option("--intwidth",
        action="store",
        default=32,
        choices=[32, 64],
        metavar="32,64",
        help="ESPRESO integer datatype width [default: %default]")

    modes=["release", "development", "debug"]
    espreso.add_option("-m", "--mode",
        action="store",
        default="release",
        choices=modes,
        help="ESPRESO build mode: " + ", ".join(modes) + " [default: %default]")

    solvers=["mkl"]
    espreso.add_option("--solver",
        action="store",
        default="mkl",
        choices=solvers,
        help="ESPRESO solver " + ", ".join(solvers) + " [default: %default]")

    espreso.add_option("--parmetisroot",
        action="store",
        type="string",
        default="",
        help="Path to ParMETIS.")

    espreso.add_option("--mklroot",
        action="store",
        type="string",
        default="",
        help="Path to MKL.")

    espreso.add_option("--hypreroot",
        action="store",
        type="string",
        default="",
        help="Path to HYPRE solver.")

