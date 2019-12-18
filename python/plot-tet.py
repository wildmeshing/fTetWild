import argparse
import numpy as np

import build_cut_table as tt

import plotly.graph_objs as go
import plotly.offline as plotly


cols = [
    '#55efc4',
    '#00b894',
    '#ffeaa7',
    '#fdcb6e',
    '#81ecec',
    '#00cec9',
    '#fab1a0',
    '#e17055',
    '#74b9ff',
    '#0984e3',
    '#ff7675',
    '#d63031',
    '#a29bfe',
    '#6c5ce7',
    '#fd79a8',
    '#e84393',
    '#dfe6e9',
    '#b2bec3',
    '#636e72',
    '#2d3436'
]


def get_layout(a, size=None):
    layout = go.Layout(
        autosize=True,
        width=(None if size is None else size[0]),
        height=(None if size is None else size[1]),
        showlegend=False,
        scene=go.layout.Scene(
            aspectmode='data',
            annotations=a,
            xaxis=dict(
                autorange=True,
                # showgrid=False,
                # zeroline=False,
                # showline=False,
                # ticks='',
                # showticklabels=False
            ),
            yaxis=dict(
                autorange=True,
                showgrid=False,
                zeroline=False,
                # showline=False,
                # ticks='',
                # showticklabels=False
            ),
            zaxis=dict(
                autorange=True,
                # showgrid=False,
                # zeroline=False,
                # showline=False,
                # ticks='',
                # showticklabels=False
            )
        )
    )

    return layout


def parse_args():
    parser = argparse.ArgumentParser(description='Plot tet')
    parser.add_argument('--bin', type=str, help='binary representation')
    parser.add_argument('--edges', type=str, help='edges')
    parser.add_argument('--int', type=int, help='integer representation')
    parser.add_argument('--conf', default=0, type=int, help='configuration')
    parser.add_argument('--all', type=int, help='plot_all')
    parser.add_argument('--tets', type=str, help='load tets')
    parser.add_argument('--vrts', type=str, help='load vertices')

    return parser.parse_args()



def get_vertices(int_flag, table):
    vertices = table.vertices

    flags = table.int_to_flags(int_flag)
    for e in range(6):
        if flags[5-e] == '1':
            edge = table.edges[e, :]
            p = (vertices[edge[0], :]+vertices[edge[1], :]) / 2
            vertices = np.append(vertices, [p], axis=0)

    tmp = table.new_points[int_flag]
    if len(tmp) > 0:
        vertices = np.append(vertices, tmp, axis=0)
    return vertices


def tet_trace(table, v, tet, edges, track, original_faces, col, shrink):
    bary = v[tet, :].mean(0)
    vertices = v - bary
    vertices *= shrink
    vertices += bary

    # print(original_faces)

    connectivity = np.array([[tet[0], tet[1], tet[2]], [tet[0], tet[1], tet[3]], [tet[0], tet[2], tet[3]], [tet[1], tet[2], tet[3]]])

    mesh = go.Mesh3d(
        x=vertices[:, 0], y=vertices[:, 1], z=vertices[:, 2],
        i=connectivity[:, 0], j=connectivity[:, 1], k=connectivity[:, 2],
        hoverinfo='none',
        alphahull=-1, opacity=0.4, color=col)

    trace = [mesh]


    annotations = []

    for i in range(len(original_faces)):
        f = original_faces[i]
        annotations.append(
            dict(
                showarrow=False,
                x=vertices[tet[i], 0],
                y=vertices[tet[i], 1],
                z=vertices[tet[i], 2],
                text="{}".format(f),
                xanchor="left",
                xshift=10,
                opacity=0.9
            )
        )

    # ['solid', 'dot', 'dash', 'longdash', 'dashdot', 'longdashdot']

    for e in table.edges:
        v0 = tet[e[0]]
        v1 = tet[e[1]]

        dash = 'solid'

        for ee in edges:
            if ee[0] == v0 and ee[1] == v1:
                dash = 'longdash'
                break
            if ee[1] == v0 and ee[0] == v1:
                dash = 'longdash'
                break

        tmp = go.Scatter3d(
            x=[vertices[v0, 0], vertices[v1, 0]],
            y=[vertices[v0, 1], vertices[v1, 1]],
            z=[vertices[v0, 2], vertices[v1, 2]],
            marker=dict(color=col),
            line=dict(color=col, width=10,dash=dash),
            # hoverinfo='none'
        )
        trace.append(tmp)

    for i in range(len(track)):
        if track[i]:
            trace.append(
                go.Scatter3d(
                    x=[vertices[tet[i], 0]],
                    y=[vertices[tet[i], 1]],
                    z=[vertices[tet[i], 2]],
                    marker=dict(color='black'),
                    # hoverinfo='none'
                )
            )

    annotations = []
    return trace, annotations


def tets_trace(int_flag, conf_index, table):
    trace = []
    annotations = []
    vertices = get_vertices(int_flag, table)
    tets = table.table[int_flag][conf_index]
    edges  = table.edges_table[int_flag][conf_index]
    track  = table.track_faces[int_flag][conf_index]
    original_faces  = table.original_faces[int_flag][conf_index]

    print("---------------")
    print(str(int_flag))
    print("\tconfiguration {}".format(conf_index))
    print("\t"+ str(tets))
    print("\t"+ str(edges))
    print("\t"+ str(track))
    print("\t"+ str(original_faces))
    print("---------------\n")


    for i in range(len(tets)):
        tet = tets[i]

        tmp, a = tet_trace(table, vertices, tet, edges, track[i], original_faces[i], cols[i], .7)
        trace += tmp
        annotations += a

    return trace, annotations


def html_div_for_index(index, conf_index, table, include_js):
    trace, annotations = tets_trace(index, conf_index, table)

    fig = go.Figure(data=trace, layout=get_layout(annotations, [450, 450]))
    html  = "<div style='height: 480px;width: 450px; float: left'>"
    html += plotly.plot(fig, output_type='div', include_plotlyjs=include_js, show_link=False)
    html += "<div style='text-align: center; font-weight: bold; width: 450px'>{} - {} ({})</div>".format(index, table.int_to_flags(index), conf_index)
    html += "</div>"

    return html


def main(args):
    table = tt.CutTable()

    if args.all is not None:
        html = "<html><div style='display: inline;'>"
        number = 0
        if args.all < 0:
            index = -args.all
            tet = table.table[index]

            if tet is None:
                return
            for j in range(len(tet)):
                html += html_div_for_index(index, j, table, number == 0)
                number += 1

            # for i in range(len(table.table)):
            #     tet = table.table[i]

            #     if tet:
            #         html += html_div_for_index(i, args.conf, table, number == 0)
            #         number += 1
        else:
            indices = []
            if args.all == 0:
                indices.append(table.index_for_edges([]))
            elif args.all == 1:
                indices.append(table.index_for_edges([0, 2, 3]))
                # indices.append(table.index_for_edges([0, 1, 4]))
                # indices.append(table.index_for_edges([1, 2, 5]))
                # indices.append(table.index_for_edges([3, 4, 5]))
            elif args.all == 2:
                indices.append(table.index_for_edges([0, 1, 3, 5]))
                indices.append(table.index_for_edges([1, 2, 3, 4]))
                indices.append(table.index_for_edges([0, 2, 4, 5]))
            elif args.all == 3:
                indices.append(table.index_for_edges([0, 1]))
                indices.append(table.index_for_edges([0, 2]))
                indices.append(table.index_for_edges([0, 3]))
                indices.append(table.index_for_edges([0, 4]))

                indices.append(table.index_for_edges([1, 2]))
                indices.append(table.index_for_edges([1, 4]))
                indices.append(table.index_for_edges([1, 5]))

                indices.append(table.index_for_edges([2, 3]))
                indices.append(table.index_for_edges([2, 5]))

                indices.append(table.index_for_edges([3, 4]))
                indices.append(table.index_for_edges([3, 5]))

                indices.append(table.index_for_edges([4, 5]))

            for index in indices:
                tet = table.table[index]

                if tet:
                    for j in range(len(tet)):
                        html += html_div_for_index(index, j, table, number == 0)
                        number += 1

        html += "</div></html>"

        with open("temp-plot.html", "w") as f:
            f.write(html)
    elif args.tets is not None:
        assert(args.vrts is not None)

        tets = []
        vertices = []
        with open(args.tets, 'r') as f:
            line = f.readline()
            while line:
                tmp = line.split(' ')

                assert(len(tmp) == 4)
                tet = [int(tmp[0]), int(tmp[1]), int(tmp[2]), int(tmp[3])]
                tets.append(tet)
                line = f.readline()

        with open(args.vrts, 'r') as f:
            line = f.readline()
            while line:
                tmp = line.split(' ')

                assert(len(tmp) == 3)
                vertex = [float(tmp[0]), float(tmp[1]), float(tmp[2])]
                vertices.append(vertex)
                line = f.readline()

        vertices = np.array(vertices)

        trace = []
        annotations = []
        for i in range(len(tets)):
            tet = tets[i]
            tmp, a = tet_trace(table, vertices, tet, [], [], [], cols[i % len(cols)], .7)
            trace += tmp
            annotations += a

        fig = go.Figure(data=trace, layout=get_layout(annotations))
        plotly.plot(fig)

    else:
        int_flag = -1

        if args.bin:
            int_flag = 0
            mask = 1

            n_args = len(args.bin)
            assert(n_args == table.n_flags)

            for i in range(table.n_flags):
                c = args.bin[table.n_flags-1-i]
                if c == '1':
                    int_flag = int_flag | mask

                mask = mask << 1
        elif args.edges:
            etmp = []
            for c in args.edges:
                etmp.append(int(c))

            int_flag = table.index_for_edges(etmp)
        elif args.int:
            int_flag = args.int


        if int_flag < 0:
            vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
            tet = [0, 1, 2, 3]
            trace, annotations = tet_trace(table, vertices, tet, [], [], [], cols[0], 1)
        else:
            trace, annotations = tets_trace(int_flag, args.conf, table)


        fig = go.Figure(data=trace, layout=get_layout(annotations))

        plotly.plot(fig)


if __name__ == '__main__':
    args = parse_args()
    main(args)
