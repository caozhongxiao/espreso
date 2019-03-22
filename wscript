
import sys, os, logging, subprocess

reload(sys)
sys.setdefaultencoding('utf8')

from waflib import Logs
from waflib.Build import BuildContext
class ShowConfiguration(BuildContext):
    cmd = "show"
    fun = "show"

def set_compiler(ctx):
    def trycompiler():
        try:
            ctx.find_program(ctx.options.mpicxx, var="MPI_CXX")
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
    ctx.env.CXX = ctx.env.LINK_CXX = ctx.env.MPI_CXX

def set_openmp(ctx):
    ctx.env.append_unique("CXXFLAGS", [ "-fopenmp" ])
    ctx.env.append_unique("LINKFLAGS", [ "-fopenmp" ])

def set_metis(ctx):
    includes = []
    libpath = []
    if ctx.options.parmetis:
        includes = [ ctx.options.parmetis + "/includes" ]
        libpath = [ ctx.options.parmetis + "/lib" ]

    ctx.check_cxx(
        header_name="metis.h parmetis.h", lib=["parmetis", "metis"], uselib_store="METIS",
        define_name="",
        defines=["HAVE_METIS", "IDXTYPEWIDTH=" + ctx.options.intwidth, "REALTYPEWIDTH=32"],
        includes=includes, libpath=libpath,
        msg="Checking for library ParMETIS",
        errmsg="set 'parmetis'.")

    ctx.check_cxx(
        fragment='''
            #include <metis.h>
            #include "mpi.h"
            int main(int argc, char** argv) {{ MPI_Init(&argc, &argv); return IDXTYPEWIDTH != {0}; MPI_Finalize(); }}
        '''.format(ctx.options.intwidth),
        execute=True,
        includes=includes, libpath=libpath,
        msg="Checking for METIS IDXTYPEWIDTH="+ctx.options.intwidth)

def set_mkl(ctx):
    includes = []
    libpath = []
    if ctx.options.mkl:
        includes = [ ctx.options.mkl + "/include" ]
        libpath = [ ctx.options.mkl + "/lib" ]

    defines = [ "HAVE_MKL" ]
    libs = []
    if ctx.options.intwidth == "32":
        defines.append("MKL_INT=int")
        libs.append("mkl_intel_lp64")
    if ctx.options.intwidth == "64":
        defines.append("MKL_INT=long")
        libs.append("mkl_intel_ilp64")
    libs.append("mkl_core")
    if ctx.options.cxx == "icpc":
        libs.append("mkl_intel_thread")
        libs.append("mkl_blacs_intelmpi_lp64")
    if ctx.options.cxx == "g++":
        libs.append("mkl_gnu_thread")
        libs.append("mkl_blacs_openmpi_lp64")

    ctx.check_cxx(
        header_name="mkl.h", lib=libs, uselib_store="MKL",
        define_name="", defines=defines,
        includes=includes, libpath=libpath,
        msg="Checking for library MKL",
        errmsg="set 'mkl'.")

def set_info(ctx):
    ctx.env.append_unique("DEFINES_INFO", [ "__ESMODE__="+ctx.options.mode.upper() ])

def set_async(ctx):
    if ctx.options.without_async:
        ctx.env.append_unique("DEFINES_ASYNC", [ "WITHOUT_ASYNC" ])
    else:
        ctx.env.append_unique("INCLUDES_ASYNC", [ "tools/ASYNC" ])

def try_hypre(ctx):
    includes = []
    libpath = []
    if ctx.options.hypre:
        includes = [ ctx.options.hypre + "/include" ]
        libpath = [ ctx.options.hypre + "/lib" ]

    ctx.check_cxx(
        header_name="HYPRE.h", lib="HYPRE", uselib_store="HYPRE",
        define_name="", defines="HAVE_HYPRE",
        includes=includes, libpath=libpath,
        mandatory=False,
        msg="Checking for library HYPRE",
        errmsg="set 'hypre' to use HYPRE solver.")

def try_bem(ctx):
    includes = []
    libpath = []
    if ctx.options.bem:
        includes = [ ctx.options.bem + "/include" ]
        libpath = [ ctx.options.bem + "/lib" ]

    ctx.check_cxx(
        header_name="heatdtn.h", lib="heatdtn_int"+ctx.options.intwidth, uselib_store="BEM",
        define_name="", defines="HAVE_BEM",
        includes=includes, libpath=libpath,
        mandatory=False,
        msg="Checking for library HeatDTN",
        errmsg="set 'bem' to use BEM assembler.")

def try_catalyst(ctx):
    pass

def set_variables(ctx):
    if ctx.options.intwidth == "32":
        ctx.env.append_unique("DEFINES", [ "esint=int", "esint_mpi=MPI_INT" ])
    if ctx.options.intwidth == "64":
        ctx.env.append_unique("DEFINES", [ "esint=long", "esint_mpi=MPI_LONG" ])

    ctx.env.append_unique("CXXFLAGS", [ "-std=c++11", "-Wall" ])
    ctx.env.append_unique("CXXFLAGS", [ "-Wno-deprecated-declarations" ]) # MKL
    ctx.env.append_unique("CXXFLAGS", [ "-Wno-format-security" ]) # loggers
    ctx.env.append_unique("CXXFLAGS", ctx.options.cxxflags.split())
    if ctx.options.mode == "release":
        ctx.env.append_unique("CXXFLAGS", [ "-O3", "-g" ])
    if ctx.options.mode == "measurement":
        ctx.env.append_unique("CXXFLAGS", [ "-O3", "-g" ])
    if ctx.options.mode == "devel":
        ctx.env.append_unique("CXXFLAGS", [ "-O2", "-g" ])
    if ctx.options.mode == "debug":
        ctx.env.append_unique("CXXFLAGS", [ "-O0", "-g" ])

    ctx.env.append_unique("INCLUDES", "src")

def configure(ctx):
    set_compiler(ctx)
    set_variables(ctx)
    set_openmp(ctx)
    set_mkl(ctx)
    set_metis(ctx)
    set_async(ctx)
    set_info(ctx)
    try_hypre(ctx)
    try_bem(ctx)

    ctx.msg("Setting compiler to", ctx.options.mpicxx)
    ctx.msg("Setting int width to", ctx.options.intwidth)
    ctx.msg("Setting build mode to", ctx.options.mode)
    ctx.msg("Setting solver to", ctx.options.solver)

def show(ctx):
    ctx.logger = logging.getLogger('show')
    ctx.logger.handlers = Logs.log_handler()

    ctx.msg("CXX", " ".join(ctx.env.CXX))
    ctx.msg("INFO", " ".join(ctx.env.DEFINES_INFO))
    ctx.msg("DEFINES", " ".join(ctx.env.DEFINES))
    ctx.msg("CXXFLAGS", " ".join(ctx.env.CXXFLAGS))
    ctx.msg("MKL", "HAVE_MKL" in ctx.env.DEFINES_MKL)
    ctx.msg("METIS", "HAVE_METIS" in ctx.env.DEFINES_METIS)
    ctx.msg("HYPRE", "HAVE_HYPRE" in ctx.env.DEFINES_HYPRE)
    ctx.msg("BEM4I", "HAVE_BEM" in ctx.env.DEFINES_BEM)

fetisources= (
   "src/solver/generic/Domain.cpp",
   "src/solver/generic/SparseMatrix.cpp",
   "src/solver/generic/utils.cpp",
   "src/solver/generic/timeeval.cpp",
   "src/solver/generic/FETISolver.cpp",
   "src/solver/specific/cluster.cpp",
   "src/solver/specific/itersolver.cpp",
   "src/solver/specific/cpu/SparseSolverMKL.cpp",
   "src/solver/specific/cpu/clustercpu.cpp",
   "src/solver/specific/cpu/itersolvercpu.cpp",
   "src/solver/specific/cpu/DenseSolverMKL.cpp"
)

def build(ctx):
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"]).rstrip()
    ctx.env.append_unique("DEFINES_INFO", [ '__ESCOMMIT__=\"{0}\"'.format(commit) ])

    ctx.shlib(source=ctx.path.ant_glob('src/basis/**/*.cpp'),target="basis")
    ctx.shlib(source=ctx.path.ant_glob('src/esinfo/**/*.cpp'),target="esinfo",use="INFO")
    ctx.shlib(source=ctx.path.ant_glob('src/config/**/*.cpp'),target="config")
    ctx.shlib(source=ctx.path.ant_glob('src/mesh/**/*.cpp'),target="mesh")
    ctx.shlib(source=ctx.path.ant_glob('src/input/**/*.cpp'),target="input")
    ctx.shlib(source=ctx.path.ant_glob('src/output/**/*.cpp'),target="output",use="ASYNC")
    ctx.shlib(source=ctx.path.ant_glob('src/physics/**/*.cpp'),target="physics",defines=["SOLVER_MKL"])
    ctx.shlib(source=ctx.path.ant_glob('src/linearsolver/**/*.cpp'),target="linearsolver")
    ctx.shlib(source=ctx.path.ant_glob('src/wrappers/metis/**/*.cpp'),target="wmetis",use="METIS")
    ctx.shlib(source=ctx.path.ant_glob('src/wrappers/math/**/*.cpp'),target="wmath",use="MKL")
    ctx.shlib(source=ctx.path.ant_glob('src/wrappers/mklpdss/**/*.cpp'),target="wmklpdss",use="MKL")
    ctx.shlib(source=ctx.path.ant_glob('src/wrappers/hypre/**/*.cpp'),target="whypre",use="HYPRE")
    ctx.shlib(source=ctx.path.ant_glob('src/wrappers/bem/**/*.cpp'),target="wbem",use="BEM")
    ctx.shlib(source=ctx.path.ant_glob('src/wrappers/catalyst/**/*.cpp'),target="wcatalyst",use="CATALYST")
    ctx.shlib(source=fetisources,target="feti",defines=["SOLVER_MKL"])

    ctx.program(
        source="src/app/espreso.cpp",target="espreso",
        use="config basis esinfo mesh input output physics feti linearsolver wmath wmklpdss whypre wcatalyst wmetis wbem",
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

    espreso.add_option("--cxxflags",
        action="store",
        type="string",
        default="",
        help="C++ compiler flags (space separated list)")

    espreso.add_option("--intwidth",
        action="store",
        default="32",
        choices=["32", "64"],
        metavar="32,64",
        help="ESPRESO integer datatype width [default: %default]")

    modes=["release", "measurement", "devel", "debug"]
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

    espreso.add_option("--parmetis",
        action="store",
        type="string",
        metavar="ParMETIS_ROOT",
        default="",
        help="Path to ParMETIS.")

    espreso.add_option("--mkl",
        action="store",
        type="string",
        metavar="MKL_ROOT",
        default="",
        help="Path to MKL.")

    espreso.add_option("--hypre",
        action="store",
        type="string",
        metavar="HYPRE_ROOT",
        default="",
        help="Path to HYPRE solver.")

    espreso.add_option("--bem",
        action="store",
        type="string",
        metavar="BEM4I_ROOT",
        default="",
        help="Path to BEM4I assembler.")

    espreso.add_option("--without-async",
        action="store_true",
        default=False,
        help="Build without asynchronous output.")
