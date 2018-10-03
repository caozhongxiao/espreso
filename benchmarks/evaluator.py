
import os
import sys

ESPRESO_TESTS = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.abspath(__file__))

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

class ESPRESOTestEvaluator:

    def __init__(self, root):
        self.ROOT = root
        self.levels = []
        self.table = {}
        self.load()

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
                break
            else:
                self.levels.append(dirs)

        fill_values(self.table, self.levels)

        repetitions = max([int(file.split(".")[-2]) for file in os.listdir(log)]) + 1

        files = len(os.listdir(log))
        for file in os.listdir(log):
            position = file.split(".")[0:len(self.levels)]
            repetition = int(file.split(".")[-2])
            data = {}
            for line in open(os.path.join(log, file), "r"):
                if "avg.:" in line:
                    key = line.split("avg.:")[0]
                    if key not in data:
                        data[key] = []
                    data[key].append(float(" ".join(line.split()).split("avg.: ")[1].split()[0]))

            values = self.get_values(self.table, position)
            for key in data:
                value = sum(data[key]) / len(data[key])
                if key not in values:
                    values[key] = [None] * repetitions
                values[key][repetition] = value

        def average(table, levels):
            if len(levels):
                for value in levels[0]:
                    average(table[value], levels[1:])
            else:
                for value in table:
                    first = table[value][0]
                    rest = [v for v in table[value][1:] if v is not None]
                    table[value] = {
                        "first": first,
                        "rest": (sum(rest) / len(rest), None)[len(rest) == 0],
                        "count": len(rest),
                        "min": min(rest),
                        "max": max(rest),
                        "ratio": max(rest) / min(rest)}
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
        def add_files(tablename, stats):
            tables[tablename] = {}
            tables[tablename]["file"] = {}
            tables[tablename]["file"]["rest"] = open(os.path.join(stats, tablename + ".sum.csv"), "w")
            tables[tablename]["file"]["first"] = open(os.path.join(stats, tablename + ".first.csv"), "w")
            tables[tablename]["file"]["min"] = open(os.path.join(stats, tablename + ".min.csv"), "w")
            tables[tablename]["file"]["max"] = open(os.path.join(stats, tablename + ".max.csv"), "w")
            tables[tablename]["file"]["ratio"] = open(os.path.join(stats, tablename + ".ratio.csv"), "w")
            tables[tablename]["file"]["count"] = open(os.path.join(stats, tablename + ".count.csv"), "w")

        def create_table(*args):
            args = [ arg.values()[0].replace(" ", "_") for arg in args ]
            name = "{0}.{1}".format(tablename, ".".join(args))
            add_files(name, stats)

        if len(freevars):
            ESPRESOTestEvaluator.iterate(create_table, *freevars)
        else:
            add_files(tablename, stats)

        tablesrows = []
        tablescolumns = []
        def add_row(*args):
            tablesrows.append(".".join([ arg.values()[0].replace(" ", "_") for arg in args ]))
        def add_column(*args):
            tablescolumns.append(".".join([ arg.values()[0].replace(" ", "_") for arg in args ]))
        ESPRESOTestEvaluator.iterate(add_row, *trows)
        ESPRESOTestEvaluator.iterate(add_column, *tcolumns)

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

        ESPRESOTestEvaluator.iterate(write_value, *ranges)

        cwidth = len(max(tablescolumns, key=lambda x: len(x)))
        rwidth = len(max(tablesrows, key=lambda x: len(x)))
        for stat in [ "rest", "first", "min", "max", "count", "ratio" ]:
            for table in tables:
                tables[table]["file"][stat].write("{0:{width}} ; ".format("", width=rwidth))
                for column in sorted(tablescolumns):
                    tables[table]["file"][stat].write("{0:^{width}} ; ".format(column, width=cwidth))
                tables[table]["file"][stat].write("\n")
                for row in sorted(tablesrows):
                    tables[table]["file"][stat].write("{0:<{width}} ; ".format(row, width=rwidth))
                    for column in sorted(tablescolumns):
                        if tables[table]["data"][row][column] is not None and tables[table]["data"][row][column][stat] is not None:
                            tables[table]["file"][stat].write("{0:>{width}.3f} ; ".format(tables[table]["data"][row][column][stat], width=cwidth))
                        else:
                            tables[table]["file"][stat].write("{0:>{width}} ; ".format("", width=cwidth))
                    tables[table]["file"][stat].write("\n")

