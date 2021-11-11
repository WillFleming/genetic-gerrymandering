import faulthandler
import argparse
import functools
import itertools
from random import shuffle

import geopandas
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import nxmetis
import pandas as pd
import pysal
from deap import base, creator, tools
from deap.algorithms import eaMuPlusLambda
from nxmetis.types import MetisOptions
from shapely.geometry.multipolygon import MultiPolygon
from shapely.geometry.polygon import Polygon

plt.switch_backend('agg')
faulthandler.enable()


def read_shp(fn):
    df = geopandas.read_file(fn)
    usecols = [
        'STATEFP10',
        'COUNTYFP10',
        'VTDST10',
        'GEOID10',
        'VTDI10',
        'NAME10',
        'NAMELSAD10',
        'VAP',
        'USCDV2010', 'USCRV2010',
        'geometry'
    ]
    df = df[usecols]
    df['GEOID10'] = df.GEOID10.apply(lambda x: x.encode("utf-8"))
    return(df)


def shp_to_queen(fn):
    A, idx = pysal.queen_from_shapefile(fn, idVariable='GEOID10').full()
    return A, idx


def adj_to_graph(A, df):
    G = nx.from_numpy_matrix(A)
    for i, v in list(enumerate(df.VAP.values)):
        G.node[i]['VAP'] = v
    return(G)


def random_balanced_partition(G, n):
    w = np.random.randint(0, 100, len(G.edges.keys()))
    rand_edge_weights = {k: v for k, v in zip(G.edges.keys(), w)}
    nx.set_edge_attributes(G, rand_edge_weights, 'rand_edge_weight')
    recursive = False  # np.random.random() < .5
    options = MetisOptions(contig=True)
    objval, parts = nxmetis.partition(
        G,
        nparts=n,
        node_weight='VAP',
        node_size='size',
        edge_weight='rand_edge_weight',  # 'weight',
        tpwgts=None,
        ubvec=None,
        options=options,
        recursive=recursive
    )
    part_nums = []
    for i in range(len(parts)):
        for j in parts[i]:
            part_nums.append((i, j))
    part_nums.sort(key=lambda x: x[1])
    part_nums = [i[0] for i in part_nums]
    # partition = {geo_id: part_num for geo_id, part_num in zip(geo_ids,part_nums)} dict version
    return(part_nums)


def save_map(partition, data, highlight, name='test'):
    if highlight:
        color_parts = [p if i not in highlight else -
                       i-1 for i, p in enumerate(partition)]
    else:
        color_parts = partition
    districts = data.dissolve(color_parts)
    f, ax = plt.subplots(1, figsize=(16, 16), dpi=200)
    districts['part'] = districts.index.values.astype(str)
    png = districts.plot(column='part', axes=ax,
                         edgecolor='black', linewidth=.1)
    plt.savefig('map_'+str(name)+'.png')


def mu_swap_geo(ind, data, A, n):
    per = data.VAP.sum()/len(set(ind))
    current = data.groupby(ind)['VAP'].sum()
    thres = (current - per).abs().mean()
    arr = np.array(ind)
    mat = (arr[:, None] == arr[None, :]).astype(int)
    mat = 1-mat
    check = mat * A
    a, b = (check > 0).nonzero()
    neighbors = list(zip(a, b))
    np.random.shuffle(neighbors)
    found = 0
    for pair in neighbors:
        if found > 0:
            break
        ind_ = list(ind)
        ind_[pair[0]] = ind[pair[1]]
        new = data.groupby(ind_)['VAP'].sum()
        imb = (new - per).abs().mean()
        if imb < thres:
            ind[:] = ind_
            found += 1
            continue
        else:
            ind_ = list(ind)
            ind_[pair[1]] = ind[pair[0]]
            new = data.groupby(ind_)['VAP'].sum()
            imb = (new - per).abs().mean()
        if imb < thres:
            ind[:] = ind_
            found += 1
        else:
            continue
    assert len(set(list(ind))) == n
    return(ind,)


def fast_check_intersect(df1, df2):
    arr1 = df1.geometry.values.copy()
    arr2 = df2.geometry.values.copy()
    np.random.shuffle(arr1)
    np.random.shuffle(arr2)
    pairs = itertools.product(arr1, arr2)
    for p1, p2 in pairs:
        if p1 is None or p2 is None:
            continue
        elif p1.intersects(p2):
            return(True)
    return(False)


def cx_dissolve_difference(ind1, ind2, data, G, n):
    MAX_FOUND = 15
    df = data[['GEOID10', 'geometry']].copy(deep=True)
    inds = (list(ind1), list(ind2))

    df['part'] = list(ind1)
    df1 = df.copy()

    df['part'] = list(ind2)
    df2 = df.copy()

    dfs = (df1, df2)
    disjoint1 = []
    disjoint2 = []
    combined = []
    d1 = [(0, i) for i in range(max(ind1))]
    d2 = [(1, i) for i in range(max(ind2))]
    np.random.shuffle(d1)
    np.random.shuffle(d2)
    choices = list(itertools.chain.from_iterable(zip(d1, d2)))
    n_found = 0
    for district in choices:
        ind = inds[district[0]]
        df = dfs[district[0]]
        df = df[df.part == district[1]]
        if len(combined) > 0:
            intersection = geopandas.sjoin(
                df, combined, op='intersects', how='inner')
        else:
            intersection = []
        found_disjoint = True if len(intersection) == 0 else False
        if found_disjoint:
            if district[0] == 0:
                disjoint1.append(district[1])
            if district[0] == 1:
                disjoint2.append(district[1])
            if len(combined) == 0:
                combined = geopandas.GeoDataFrame()
            combined = pd.concat([combined, df])
            combined.reset_index(drop=True)
            n_found += 1
        else:
            pass
        if n_found >= MAX_FOUND:
            break

    data['ind1'] = list(ind1)
    data['ind2'] = list(ind2)
    data['part'] = np.where(data.ind1.isin(disjoint1), 1, 0)
    data['part'] = np.where(data.ind2.isin(disjoint2), 2, data.part)
    to_fill = data.part == 0
    keep_i1 = data.part == 1
    keep_i2 = data.part == 2
    fill_df = data.loc[to_fill]
    fill_G = G.subgraph(fill_df.index.values).copy()
    dis_subs = list(nx.connected_component_subgraphs(fill_G))
    total_keep_n = len(disjoint1 + disjoint2)
    total_fill_n = max(ind1)+1 - total_keep_n
    ds_parts = {}
    ds_fill_n = []
    running_idx = 0
    for ds_ in dis_subs:
        ds = ds_.copy()
        fill_n = int(round(data.loc[list(ds)].VAP.sum().astype(
            float) / (data.VAP.sum()/len(set(ind1)))))
        if fill_n < 1:
            print("Warning: cx failed")
            return(ind1, ind2)
        ds_fill_n.append(fill_n)
        ds_part = random_balanced_partition(ds, fill_n)
        ds_part = {i: p+running_idx for i, p in zip(sorted(list(ds)), ds_part)}
        running_idx = max(ds_part.values()) + 1
        ds_parts.update(ds_part)

    assert sum(ds_fill_n) == total_fill_n
    assert to_fill.sum() == len(ds_parts)
    fill_parts = [v for k, v in sorted(ds_parts.items())]
    relabel1 = {k: v for k, v in zip(disjoint1, range(total_fill_n, n))}
    relabel2 = {k: v for k, v in zip(
        disjoint2, range(total_fill_n+len(disjoint1), n))}
    assert len(relabel1.keys() + relabel2.keys()) == total_keep_n
    data['new_part'] = None
    data.loc[to_fill, 'new_part'] = fill_parts
    data.loc[keep_i1, 'new_part'] = data.loc[keep_i1, 'ind1'].replace(relabel1)
    data.loc[keep_i2, 'new_part'] = data.loc[keep_i2, 'ind2'].replace(relabel2)
    ind_cx = data.new_part.values.tolist()
    assert len(set(ind_cx)) == n
    ind1[:] = ind_cx
    ind2[:] = ind_cx
    return(ind1, ind2)


def evaluate(individual, data):
    df = data[['GEOID10', 'USCDV2010', 'USCRV2010']].copy(deep=True)
    df['part'] = individual
    df = df.groupby('part')['USCDV2010', 'USCRV2010'].sum().reset_index()
    #df['d_win'] = df['USCDV2010']>df['USCRV2010']
    df['d_win'] = np.log(df['USCDV2010']/df['USCRV2010'])
    d_wins = df['d_win'].mean()
    return(d_wins,)


def count_poly(geometries):
    cnt = 0
    for g in geometries:
        if isinstance(g, MultiPolygon):
            for p in g:
                if isinstance(p, Polygon):
                    cnt += 1
                else:
                    raise TypeError
        elif isinstance(g, Polygon):
            cnt += 1
        else:
            raise TypeError
    return(cnt)


def main(fn, n_districts, pop_size=50, l=None, g=100, cxpb=.8, mutpb=.2):
    df = read_shp(fn)
    A, idx = shp_to_queen(fn)
    G = adj_to_graph(A, df)

    creator.create("FitnessMax", base.Fitness, weights=(1.0,))
    creator.create("Individual", list, fitness=creator.FitnessMax)
    toolbox = base.Toolbox()
    toolbox.register("random_partition",
                     random_balanced_partition, G, n=n_districts)
    toolbox.register("individual", tools.initIterate,
                     creator.Individual, toolbox.random_partition)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)
    toolbox.register("evaluate", evaluate, data=df)
    toolbox.register("mutate", mu_swap_geo, data=df, A=A, n=n_districts)
    toolbox.register("mate", cx_dissolve_difference,
                     data=df, G=G, n=n_districts)
    toolbox.register("select", tools.selBest)
    stats = tools.Statistics(key=lambda ind: ind.fitness.values)
    stats.register("avg", np.mean)
    stats.register("std", np.std)
    stats.register("min", np.min)
    stats.register("max", np.max)

    pop = toolbox.population(n=pop_size)
    lambda_ = l if l is not None else pop_size / 5
    final_pop, logbook = eaMuPlusLambda(
        pop, toolbox, mu=pop_size, lambda_=lambda_, cxpb=cxpb, mutpb=mutpb, ngen=g, stats=stats, verbose=True)
    k = 1
    df['parts'] = final_pop[k]
    districts = df.dissolve('parts')
    f, ax = plt.subplots(1, figsize=(16, 16), dpi=250)
    png = districts.plot(column=districts.index.values,
                         axes=ax, edgecolor='black', linewidth=.1)
    plt.savefig('output.png')
    print("Complete")


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument(dest="fn", type=str)
    p.add_argument("-n", dest="n_districts", type=int, required=True)
    p.add_argument("-p", dest="pop_size", type=int, required=False, default=50)
    p.add_argument("-l", dest="l", type=int, required=False, default=None)
    p.add_argument("-g", dest="g", type=int, required=False, default=100)
    p.add_argument("--cxpb", dest="cxpb", type=float,
                   required=False, default=0.8)
    p.add_argument("--mutpb", dest="mutpb", type=float,
                   required=False, default=0.2)

    args = vars(p.parse_args())
    main(**args)

    # nohup python -u geo.py shp/pa_final.shp -n 18 > geo.out 2>&1 &
