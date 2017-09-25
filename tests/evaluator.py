
import os
import sys

class ESPRESOTestEvaluator:

    def evaluate(self, root, tablename, rows, columns, value):
        rows = set(rows)
        columns = set(columns)
        table = {}

        for row in rows:
            table[row] = {}
            for column in columns:
                table[row][column] = []

        run   = os.path.join(root, "run")
        log   = os.path.join(root, "log")
        stats = os.path.join(root, "stats")

        levels = []
        for path, dirs, files in os.walk(run):
            if "results" in dirs:
                repetitions = len([file for file in files if file.endswith(".ecf")])
                break
            else:
                levels.append(dirs)

        for file in os.listdir(log):
            names = file.split(".")[:-1] # .log
            if repetitions > 1:
                names = names[:-1] # .repetition

            row = column = None
            for name in names:
                if name in rows and name in columns:
                    print "Invalid table description. '{0}' is in both rows and columns.".format(name)
                    exit()
                if name in rows:
                    row = name
                if name in columns:
                    column = name

            if len(value):
                if row is None and len(rows) == 1:
                    row = list(rows)[0]
                if column is None and len(columns) == 1:
                    column = list(columns)[0]
                if row is None or column is None:
                    print "Unrecognized table format. A row or column not found."
                    exit()

                time = None
                for line in open(os.path.join(log, file), "r"):
                    if value in line:
                        time = float(" ".join(line.split()).split("avg.: ")[1].split()[0])
                        break
                if time is None:
                    print value, "NOT FOUND"
                    exit()

                table[row][column].append(time)

            else:
                if row is not None and column is not None:
                    print "Missing parameter value in table '{0}'".format(tablename)
                    exit()

                if row is None:
                    times = [None] * len(rows)
                if column is None:
                    times = [None] * len(columns)
                for line in open(os.path.join(log, file), "r"):
                    if row is None:
                        for i, rvalue in enumerate(sorted(rows)):
                            if rvalue in line:
                                times[i] = float(" ".join(line.split()).split("avg.: ")[1].split()[0])
                    if column is None:
                        for i, cvalue in enumerate(sorted(column)):
                            if cvalue in line:
                                times[i] = float(" ".join(line.split()).split("avg.: ")[1].split()[0])
                if row is None:
                    for i, rvalue in enumerate(sorted(rows)):
                        if times[i] is None:
                            print "Missing parameter '{0}' in table '{1}'".format(rvalue, tablename)
                        else:
                            table[rvalue][column].append(times[i])
                if column is None:
                    for i, cvalue in enumerate(sorted(column)):
                        if times[i] is None:
                            print "Missing parameter '{0}' in table '{1}'".format(cvalue, tablename)
                        else:
                            table[row][cvalue].append(times[i])


        cwidth = len(max(columns, key=lambda x: len(x)))
        rwidth = len(max(rows, key=lambda x: len(x)))
        statistics = open(os.path.join(stats, "{0}.csv".format(tablename)), "w")
        statistics.write("{0:{width}} ; ".format("", width=rwidth))
        for column in sorted(columns):
            statistics.write("{0:^{width}} ; ".format(column, width=cwidth))
        statistics.write("\n")
        for row in sorted(rows):
            statistics.write("{0:<{width}} ; ".format(row, width=rwidth))
            for column in sorted(columns):
                if table[row][column][0]:
                    statistics.write("{:>{width}.3f} ; ".format(sum(table[row][column]) / len(table[row][column]), width=cwidth))
                else:
                    statistics.write("{:>{width}.3f} ; ".format(0, width=cwidth))
            statistics.write("\n")
