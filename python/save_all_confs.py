import os
import importlib
import build_cut_table as tt
pt = importlib.import_module("plot-tet")


def main():
    args = lambda: None

    folder = os.path.join(os.getcwd(), "confs")

    if not os.path.exists(folder):
        os.makedirs(folder)


    old_file = os.path.join(os.getcwd(), "temp-plot.html")
    table = tt.CutTable()
    for i in range(len(table.table)):
        args.all = -i
        pt.main(args)
        if os.path.exists(old_file):
            os.rename(old_file, os.path.join(folder, "conf_{}.html".format(i)))


if __name__ == '__main__':
    main()
