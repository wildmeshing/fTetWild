import build_cut_table as tt
import argparse

import os


def parse_args():
    parser = argparse.ArgumentParser(description='Export table as a c file')
    parser.add_argument('--filename', type=str, help='name of the file')

    return parser.parse_args()


def main():
    class_name = "CutTable"

    args = parse_args()
    save_file = None
    if args.filename:
        fname = os.path.basename(args.filename)
        save_file = os.path.join(os.getcwd(), args.filename)
    else:
        fname = "test"

    table = tt.CutTable()



    table_ret = "std::vector<std::vector<Vector4i>>"
    edges_ret = "std::vector<std::vector<Vector2i>>"
    track_ret = "std::vector<std::vector<std::array<bool, 4>>>"
    origi_ret = "std::vector<std::vector<Vector4i>>"

    table_fun_name = "get_tet_confs(const int idx)"
    track_fun_name = "get_surface_conf(const int idx)"
    origi_fun_name = "get_face_id_conf(const int idx)"
    edges_fun_name = "get_diag_confs(const int idx)"

    table_fun_name1 = "get_tet_conf(const int idx, const int cfg)"
    track_fun_name1 = "get_surface_conf(const int idx, const int cfg)"
    origi_fun_name1 = "get_face_id_conf(const int idx, const int cfg)"


    hpp = "#pragma once\n\n"
    hpp += "#include <floattetwild/Types.hpp>\n\n"
    hpp += "#include <array>\n"
    hpp += "#include <vector>\n\n"

    hpp += "namespace floatTetWild {{\n\tclass {} {{\npublic:\n".format(class_name)

    hpp += "\t\tstatic const {}& {};\n".format(table_ret, table_fun_name)
    hpp += "\t\tstatic const {}& {};\n".format(edges_ret, edges_fun_name)
    hpp += "\t\tstatic const {}& {};\n".format(track_ret, track_fun_name)
    hpp += "\t\tstatic const {}& {};\n".format(origi_ret, origi_fun_name)

    hpp += "\t\tstatic inline const {}& {}{{ return get_tet_confs(idx)[cfg]; }}\n".format("std::vector<Vector4i>", table_fun_name1)
    hpp += "\t\tstatic inline const {}& {}{{ return get_surface_conf(idx)[cfg]; }}\n".format("std::vector<std::array<bool, 4>>", track_fun_name1)
    hpp += "\t\tstatic inline const {}& {}{{ return get_face_id_conf(idx)[cfg]; }}\n".format("std::vector<Vector4i>", origi_fun_name1)


    hpp += "\t};\n}\n"


    cpp = '#include "{}.hpp"\n\n'.format(fname)
    cpp += '#include <cassert>\n\n'.format(fname)

    cpp += "namespace floatTetWild {\n"




    ###############################################################################################
    cpp += "\tconst {}& {}::{} {{\n".format(table_ret, class_name, table_fun_name)

    cpp += "\t\tstatic const std::array<{}, {}> table= {{{{\n".format(table_ret, len(table.table))
    for i in range(len(table.table)):
        cpp += "\n"
        confs = table.table[i]
        if confs is None:
            cpp += "\t\t\t{},"
            continue

        cpp += "\t\t\t{\n"
        for c in confs:
            cpp += "\n"
            cpp += "\t\t\t\t{\n\t\t\t\t\t"

            max_v = -1
            has_8 = False
            for tet in c:
                max_v = max(max_v, max(tet))
                if tet[0] == 8 or tet[1] == 8 or tet[2] == 8 or tet[3] == 8:
                    has_8 = True

            decrease_9 = max_v == 9 and not has_8
            for tet in c:
                i0 = 8 if decrease_9 and tet[0] == 9 else tet[0]
                i1 = 8 if decrease_9 and tet[1] == 9 else tet[1]
                i2 = 8 if decrease_9 and tet[2] == 9 else tet[2]
                i3 = 8 if decrease_9 and tet[3] == 9 else tet[3]
                cpp += "Vector4i({}, {}, {}, {}),".format(i0, i1, i2, i3)
            cpp = cpp[:-1]
            cpp += "\n\t\t\t\t},"
        cpp = cpp[:-1]
        cpp += "\n\t\t\t},"

    cpp = cpp[:-1]
    cpp += "\n\t\t}};\n\n"

    cpp += "\t\tassert(!table[idx].empty());\n"
    cpp += "\t\treturn table[idx];\n"
    cpp += "\t}\n\n\n"
    ###############################################################################################

    ###############################################################################################
    cpp += "\tconst {}& {}::{} {{\n".format(edges_ret, class_name, edges_fun_name)

    cpp += "\t\tstatic const std::array<{}, {}> table= {{{{\n".format(edges_ret, len(table.edges_table))
    for i in range(len(table.edges_table)):
        cpp += "\n"
        confs = table.table[i]
        if confs is None:
            cpp += "\t\t\t{},"
            continue

        edges = table.edges_table[i]
        cpp += "\t\t\t{\n"
        for c in edges:
            cpp += "\n"
            cpp += "\t\t\t\t{\n\t\t\t\t\t"
            for e in c:
                cpp += "Vector2i({}, {}),".format(e[0], e[1])
            cpp = cpp[:-1]
            cpp += "\n\t\t\t\t},"
        cpp = cpp[:-1]
        cpp += "\n\t\t\t},"

    cpp = cpp[:-1]
    cpp += "\n\t\t}};\n\n"

    cpp += "\t\t//assert(!table[idx].empty());\n"
    cpp += "\t\treturn table[idx];\n"
    cpp += "\t}\n\n\n"
    ###############################################################################################

    ###############################################################################################
    cpp += "\tconst {}& {}::{} {{\n".format(track_ret, class_name, track_fun_name)

    cpp += "\t\tstatic const std::array<{}, {}> table= {{{{\n".format(track_ret, len(table.track_faces))
    for i in range(len(table.track_faces)):
        cpp += "\n"
        confs = table.table[i]
        if confs is None:
            cpp += "\t\t\t{},"
            continue

        track = table.track_faces[i]
        cpp += "\t\t\t{\n"
        for c in track:
            cpp += "\n"
            cpp += "\t\t\t\t{\n\t\t\t\t\t"
            for t in c:
                cpp += "{{{{{}, {}, {}, {}}}}},".format(str(t[0]).lower(), str(t[1]).lower(), str(t[2]).lower(), str(t[3]).lower())
            cpp = cpp[:-1]
            cpp += "\n\t\t\t\t},"
        cpp = cpp[:-1]
        cpp += "\n\t\t\t},"

    cpp = cpp[:-1]
    cpp += "\n\t\t}};\n\n"

    cpp += "\t\tassert(!table[idx].empty());\n"
    cpp += "\t\treturn table[idx];\n"
    cpp += "\t}\n\n\n"
    ###############################################################################################

    ###############################################################################################
    cpp += "\tconst {}& {}::{} {{\n".format(origi_ret, class_name, origi_fun_name)

    cpp += "\t\tstatic const std::array<{}, {}> table= {{{{\n".format(table_ret, len(table.original_faces))
    for i in range(len(table.original_faces)):
        cpp += "\n"
        confs = table.table[i]
        if confs is None:
            cpp += "\t\t\t{},"
            continue

        cpp += "\t\t\t{\n"
        for c in table.original_faces[i]:
            cpp += "\n"
            cpp += "\t\t\t\t{\n\t\t\t\t\t"
            for face in c:
                cpp += "Vector4i({}, {}, {}, {}),".format(face[0], face[1], face[2], face[3])
            cpp = cpp[:-1]
            cpp += "\n\t\t\t\t},"
        cpp = cpp[:-1]
        cpp += "\n\t\t\t},"

    cpp = cpp[:-1]
    cpp += "\n\t\t}};\n\n"

    cpp += "\t\tassert(!table[idx].empty());\n"
    cpp += "\t\treturn table[idx];\n"
    cpp += "\t}\n\n\n"
    ###############################################################################################

    cpp += "}\n"

    if save_file:
        with open(save_file + ".hpp", 'w') as f:
            f.write(hpp)

        with open(save_file + ".cpp", 'w') as f:
            f.write(cpp)
    else:
        print(cpp)


if __name__ == '__main__':
    main()
