
import os
import sys

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(ESPRESO_TESTS, "utils"))

from testing import *

class ESPRESOTestEvaluator:

    def __init__(self, root):
        self.ROOT = root
        self.repetition = 0
        self.levels = []
        self.table = {}
        self.load()

    def get_values(self, table, position):
        if len(position):
            return self.get_values(table[position[0]], position[1:])
        else:
            return table

    def load(self):
        def fill_values(table, levels):
            if len(levels):
                for value in levels[0]:
                    table[value] = {}
                    fill_values(table[value], levels[1:])

        run = os.path.join(self.ROOT, "run")
        log = os.path.join(self.ROOT, "log")

        self.levels = []
        for path, dirs, files in os.walk(run):
            if any([file.endswith(".ecf") for file in files]):
                self.repetitions = len([file for file in files if file.endswith(".ecf")])
                break
            else:
                self.levels.append(dirs)

        fill_values(self.table, self.levels)

        files = len(os.listdir(log))
        for file in os.listdir(log):
            position = file.split(".")[0:len(self.levels)]
            data = {}
            for line in open(os.path.join(log, file), "r"):
                if "avg.:" in line:
                    key = line.split("avg.:")[0]
                    if key in data:
                        data[key].append(float(" ".join(line.split()).split("avg.: ")[1].split()[0]))
                    else:
                        data[key] = [ float(" ".join(line.split()).split("avg.: ")[1].split()[0]) ]

            values = self.get_values(self.table, position)
            for key in data:
                value = sum(data[key]) / len(data[key])
                if key in values:
                    values[key].append(value)
                else:
                    values[key] = [ value ]

        def average(table, levels):
            if len(levels):
                for value in levels[0]:
                    average(table[value], levels[1:])
            else:
                for value in table:
                    table[value] = sum(table[value]) / len(table[value])
        average(self.table, self.levels)

    def evaluate(self, tablename, rows, columns, **kwargs):
        ranges = []
        for level in range(1, len(self.levels) + 1):
            ranges.append({"L{0}".format(level): kwargs["L{0}".format(level)]})
        ranges.append({"VALUES": kwargs["VALUES"]})

        stats = os.path.join(self.ROOT, "stats")

        trows = []
        tcolumns = []
        freevars = []
        isrowvar = []
        iscolumnvar = []
        isargfree = []
        for level in ranges:
            for k, v in level.iteritems():
                if k in rows:
                    trows.append({k: v})
                    isrowvar.append(True)
                    iscolumnvar.append(False)
                    isargfree.append(False)
                if k in columns:
                    tcolumns.append({k: v})
                    isrowvar.append(False)
                    iscolumnvar.append(True)
                    isargfree.append(False)
                if k not in rows and k not in columns:
                    freevars.append({k: v})
                    isrowvar.append(False)
                    iscolumnvar.append(False)
                    isargfree.append(True)

        tables = {}
        def create_table(*args):
            args = [ arg.values()[0].replace(" ", "_") for arg in args ]
            name = "{0}.{1}".format(tablename, ".".join(args))
            tables[name] = { "file": open(os.path.join(stats, name + ".csv"), "w") }

        if len(freevars):
            TestCaseCreator.iterate(create_table, *freevars)
        else:
            tables[tablename] = { "file": open(os.path.join(stats, tablename + ".csv"), "w") }

        tablesrows = []
        tablescolumns = []
        def add_row(*args):
            tablesrows.append(".".join([ arg.values()[0].replace(" ", "_") for arg in args ]))
        def add_column(*args):
            tablescolumns.append(".".join([ arg.values()[0].replace(" ", "_") for arg in args ]))
        TestCaseCreator.iterate(add_row, *trows)
        TestCaseCreator.iterate(add_column, *tcolumns)

        for table in tables:
            tables[table]["data"] = {}
            for row in tablesrows:
                tables[table]["data"][row] = {}
                for column in tablescolumns:
                    tables[table]["data"][row][column] = 0;

        def write_value(*args):
            args = [ arg.values()[0] for arg in args ]
            values = self.get_values(self.table, args[:-1])
            value = None
            for key in values:
                if args[-1] in key:
                    value = values[key]
                    break

            table = ".".join([tablename] + [arg.replace(" ", "_") for i, arg in enumerate(args) if isargfree[i]])
            row = ".".join([arg.replace(" ", "_") for i, arg in enumerate(args) if isrowvar[i]])
            column = ".".join([arg.replace(" ", "_") for i, arg in enumerate(args) if iscolumnvar[i]])
            tables[table]["data"][row][column] = value

        TestCaseCreator.iterate(write_value, *ranges)


        cwidth = len(max(tablescolumns, key=lambda x: len(x)))
        rwidth = len(max(tablesrows, key=lambda x: len(x)))
        for table in tables:
            tables[table]["file"].write("{0:{width}} ; ".format("", width=rwidth))
            for column in sorted(tablescolumns):
                tables[table]["file"].write("{0:^{width}} ; ".format(column, width=cwidth))
            tables[table]["file"].write("\n")
            for row in sorted(tablesrows):
                tables[table]["file"].write("{0:<{width}} ; ".format(row, width=rwidth))
                for column in sorted(tablescolumns):
                    if tables[table]["data"][row][column] is not None:
                        tables[table]["file"].write("{0:>{width}.3f} ; ".format(tables[table]["data"][row][column], width=cwidth))
                    else:
                        tables[table]["file"].write("{0:>{width}} ; ".format("", width=cwidth))
                tables[table]["file"].write("\n")

