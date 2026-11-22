from collections import Counter
from itertools import combinations, product
from random import Random

import networkx as nx


def edge(u, v):
    return tuple(sorted((u, v)))


def make_spine(lengths):
    rows = [[f"v{r}.{j}" for j in range(length + 1)]
            for r, length in enumerate(lengths)]
    spine = ["start"]
    tree = nx.Graph()
    for r, row in enumerate(rows):
        if r:
            fork = f"fork{r}"
            spine.append(fork)
            tree.add_edge(fork, f"leaf{r}")
        spine.extend(row)
    spine.append("end")
    nx.add_path(tree, spine)
    return tree, rows


def weighted_tree_model(tree, chords, scale=1):
    total = sum(chords.values())
    girth = 2 * total + 1
    weights = dict(chords)
    tree_cuts = {}
    all_edges = {edge(u, v) for u, v in tree.edges} | set(chords)
    for u, v in list(tree.edges):
        tree.remove_edge(u, v)
        shore = set(nx.node_connected_component(tree, u))
        tree.add_edge(u, v)
        cut = frozenset(e for e in all_edges if (e[0] in shore) != (e[1] in shore))
        tree_cuts[edge(u, v)] = cut
        weights[edge(u, v)] = girth - sum(chords[e] for e in cut if e in chords)
    return girth * scale, {e: w * scale for e, w in weights.items()}, tree_cuts


def weight(weights, edges):
    return sum(weights[e] for e in edges)


def covered(cuts):
    return frozenset().union(*cuts)


def expanded_dual(weights):
    graph = nx.Graph()
    graph.add_edges_from(weights)
    planar, embedding = nx.check_planarity(graph)
    assert planar
    face_by_dart = {}
    seen = set()
    face = 0
    for u, v in embedding.edges():
        if (u, v) in seen:
            continue
        walk = embedding.traverse_face(u, v, mark_half_edges=seen)
        for a, b in zip(walk, walk[1:] + walk[:1]):
            face_by_dart[a, b] = face
        face += 1
    dual = nx.Graph()
    expansion = {}
    for idx, ((u, v), w) in enumerate(weights.items()):
        x, y = face_by_dart[u, v], face_by_dart[v, u]
        p, q = face + 2 * idx, face + 2 * idx + 1
        dual.add_weighted_edges_from([(x, p, 1), (p, q, 1), (q, y, 3 * w - 2)])
        expansion[u, v] = {edge(x, p), edge(p, q), edge(q, y)}
    assert dual.number_of_edges() == 3 * len(weights)
    assert not list(nx.selfloop_edges(dual))
    assert nx.check_planarity(dual)[0]
    return dual, expansion


def shortest_cycles(graph, expected_girth):
    cycles = set()
    for u, v, data in list(graph.edges(data=True)):
        w = data["weight"]
        graph.remove_edge(u, v)
        try:
            dist = nx.dijkstra_path_length(graph, u, v)
            assert w + dist >= expected_girth
            if w + dist == expected_girth:
                for path in nx.all_shortest_paths(graph, u, v, weight="weight"):
                    cycles.add(frozenset([edge(u, v)] + [
                        edge(a, b) for a, b in zip(path, path[1:])
                    ]))
        except nx.NetworkXNoPath:
            pass
        graph.add_edge(u, v, weight=w)
    return cycles


def check_context(n):
    tree, rows = make_spine([n] * 3)
    left, middle, right = rows
    chords = {edge(u, v): 1 for p, q in zip(rows, rows[1:]) for u, v in zip(p, q)}
    girth, weights, cuts = weighted_tree_model(tree, chords)
    families = [[cuts[edge(u, v)] for u, v in zip(row, row[1:])] for row in rows]
    assert len(cuts) == 3 * n + 8
    assert all(frozenset.intersection(*family) for family in families)
    graph, expansion = expanded_dual(weights)
    catalogue = {
        frozenset().union(*(expansion[e] for e in cut)) for cut in cuts.values()
    }
    assert shortest_cycles(graph, 3 * girth) == catalogue
    graph_weights = {edge(u, v): w for u, v, w in graph.edges(data="weight")}
    threshold = 3 * (13 * girth - 5 * (n - 5))
    for r in range(6, n - 4):
        outside = covered(families[0][r - 6:r] + families[2][r - 1:r + 5])
        for i, cut in enumerate(families[1], 1):
            union_in_primal = covered([expansion[e] for e in outside | cut])
            assert weight(graph_weights, union_in_primal) == threshold - 3 * abs(i - r)
    return n - 10, graph.number_of_nodes(), graph.number_of_edges()


def interval_states(path, intervals, weights, budget):
    size = len(path)
    table = {
        (a, b): weight(weights, (e for u, v, e in intervals if a <= u < v <= b))
        for a in range(size) for b in range(a, size)
    }
    total = table[0, size - 1]
    states = {}
    for a in range(1, size):
        walk = {a: (0, (a,))}
        for h in range(1, min(budget, size - a) + 1):
            for b, (cost, indices) in walk.items():
                value = total - table[0, a - 1] - cost - table[b, size - 1]
                states[h, a, b] = value, indices
            walk = {
                nxt: min((cost + table[b, nxt - 1], indices + (nxt,))
                         for b, (cost, indices) in walk.items() if b < nxt)
                for nxt in range(min(walk) + 1, size)
            } if walk else {}
    return states


def sequences(graph):
    situations = {}
    for start in graph:
        situations[frozenset([start]), frozenset()] = (start,)
        for nxt in graph[start]:
            nodes = [start, nxt]
            links = [edge(start, nxt)]
            while True:
                situations[frozenset(nodes), frozenset(links)] = tuple(nodes)
                successors = set(graph[nodes[-1]]) - {nodes[-2]}
                if not successors:
                    break
                candidate, = successors
                link = edge(nodes[-1], candidate)
                if candidate == start:
                    situations[frozenset(nodes), frozenset(links + [link])] = tuple(nodes)
                    break
                assert candidate not in nodes
                nodes.append(candidate)
                links.append(link)
    return situations


def component_optima(order, links, states, profiles, unary, weights, budget):
    cycle = len(links) == len(order)

    def value_at(row, state):
        value, chosen = states[row][state]
        return value + sum(weight(weights, profiles(row, chosen, b)) for b in unary[row])

    def between(p, s, q, t):
        bundle = links[edge(p, q)]
        return weight(weights, profiles(p, states[p][s][1], bundle)
                      | profiles(q, states[q][t][1], bundle))

    best = {}
    first = order[0]
    starts = [(s,) for s in states[first]] if cycle else [states[first]]
    for start in starts:
        layer = {(s[0], s): (value_at(first, s), (s,)) for s in start}
        for p, q in zip(order, order[1:]):
            update = {}
            for (h, s), (current, previous) in layer.items():
                for t in states[q]:
                    count = h + t[0]
                    if count <= budget:
                        candidate = current + value_at(q, t) + between(p, s, q, t)
                        key = count, t
                        if key not in update or candidate > update[key][0]:
                            update[key] = candidate, previous + (t,)
            layer = update
        for (h, last), (current, previous) in layer.items():
            if cycle:
                current += between(order[-1], last, first, previous[0])
            if h not in best or current > best[h][0]:
                best[h] = current, previous
    return best


def check_marking(tree, rows, weights, cuts):
    budget = 3
    locations = {}
    for r, row in enumerate(rows):
        for j, v in enumerate(row):
            locations[v] = r, j
    for idx, v in enumerate(set(tree) - set(locations), len(rows)):
        locations[v] = idx, 0
    groups = {
        r: dict(enumerate([cuts[edge(u, v)] for u, v in zip(row, row[1:])], 1))
        for r, row in enumerate(rows)
    }
    ordinary = {cut for group in groups.values() for cut in group.values()}
    exceptions = set(cuts.values()) - ordinary
    intervals = {r: [] for r in groups}
    bundles = {}
    for e in weights:
        (p, i), (q, j) = map(locations.get, e)
        if p == q:
            intervals[p].append((min(i, j), max(i, j), e))
        else:
            bundles.setdefault(edge(p, q), set()).add(e)
    sensitive = nx.Graph()
    sensitive.add_nodes_from(groups)
    third = {}
    for pair, bundle in bundles.items():
        p, q = pair
        third[pair] = set()
        for r, group in groups.items():
            profiles = {frozenset(cut & bundle) for cut in group.values()}
            if r not in pair:
                assert len(profiles) == 1
                if bundle & next(iter(group.values())):
                    assert profiles == {frozenset(bundle)}
                    third[pair].add(r)
            elif len(profiles) > 1:
                sensitive.add_edge(p, q)
        for cut in exceptions:
            assert not (cut & bundle) or bundle <= cut
    assert all(d <= 2 for _, d in sensitive.degree())
    states = {
        r: interval_states(row, intervals[r], weights, budget)
        for r, row in enumerate(rows)
    }
    within_edges = {r: {e for _, _, e in interval} for r, interval in intervals.items()}
    for r, opts in states.items():
        for (h, a, b), (value, indices) in opts.items():
            assert len(indices) == h and indices[0] == a and indices[-1] == b
            actual = covered(groups[r][j] for j in indices)
            assert value == weight(weights, actual & within_edges[r])
            constrained = (combo for combo in combinations(groups[r], h)
                           if combo[0] == a and combo[-1] == b)
            assert value == max(weight(weights, covered(groups[r][j] for j in combo)
                                       & within_edges[r]) for combo in constrained)

    def profile(r, opts, bundle):
        if not opts:
            return frozenset()
        return (groups[r][min(opts)] | groups[r][max(opts)]) & bundle

    def factored(selected):
        active = {r: tuple(j for j, cut in opts.items() if cut in selected)
                  for r, opts in groups.items()}
        used_exceptions = covered(exceptions & set(selected))
        score = sum(weight(weights, covered(groups[r][j] for j in indices)
                           & within_edges[r]) for r, indices in active.items())
        for pair, bundle in bundles.items():
            p, q = pair
            if bundle & used_exceptions or any(active[r] for r in third[pair]):
                score += weight(weights, bundle)
            else:
                score += weight(weights, profile(p, active.get(p), bundle)
                                | profile(q, active.get(q), bundle))
        return score

    marked = set(exceptions)
    candidates = sequences(sensitive.subgraph(groups))
    cycles = 0
    for (nodes, links), order in candidates.items():
        if len(nodes) > budget:
            continue
        boundary = {edge(p, q) for p, q in sensitive.edges(nodes)} - links
        assert len(boundary) <= 2
        if len(links) == len(nodes):
            cycles += 1
        for flags in product([False, True], repeat=len(boundary)):
            unary = {r: [] for r in nodes}
            for pair, precovered in zip(boundary, flags):
                if not precovered and len(nodes.intersection(pair)) == 1:
                    r, = nodes.intersection(pair)
                    unary[r].append(bundles[pair])
            best = component_optima(order, {e: bundles[e] for e in links}, states,
                                    profile, unary, weights, budget)
            for _, choices in best.values():
                for r, state in zip(order, choices):
                    marked.update(groups[r][j] for j in states[r][state][1])

    count = 0
    for k in range(1, budget + 1):
        scores = []
        for chosen in combinations(cuts.values(), k):
            exact = weight(weights, covered(chosen))
            assert factored(chosen) == exact
            scores.append(exact)
            count += 1
        assert max(scores) == max(weight(weights, covered(chosen))
                                  for chosen in combinations(marked, k))
    return count, cycles, sum(map(len, third.values())), len(cuts) - len(marked)


def weighted_cactus_model(rand, size, length, scale):
    rows = [[f"row{r}.{j}" for j in range(length + 1)] for r in range(size)]
    roots = [f"root{r}" for r in range(size)]
    cycle = [edge(u, v) for u, v in zip(roots, roots[1:] + roots[:1])]
    cactus = nx.Graph()
    cactus.add_edges_from(cycle)
    for r, row in enumerate(rows):
        nx.add_path(cactus, [roots[r], *row, f"tip{r}"])
    chords = {}
    loads = {}
    for r, core_edge in enumerate(cycle):
        neighbor = rows[(r + 1) % size]
        bundle = {edge(u, v): rand.randint(1, 20)
                  for u, v in zip(rows[r], neighbor) if rand.random() < 0.8}
        loads[core_edge] = sum(bundle.values())
        chords.update(bundle)
    girth = 4 * sum(chords.values()) + 2
    weights = dict(chords)
    weights.update({e: girth // 2 - loads[e] for e in cycle})
    tree_edges = {edge(u, v) for u, v in cactus.edges} - set(cycle)
    cuts = {}
    all_edges = {edge(u, v) for u, v in cactus.edges} | set(chords)
    cut_keys = [(e,) for e in tree_edges] + list(combinations(cycle, 2))
    for key in cut_keys:
        broken = cactus.copy()
        broken.remove_edges_from(key)
        shore = set(nx.node_connected_component(broken, key[0][0]))
        cut = frozenset(e for e in all_edges if (e[0] in shore) != (e[1] in shore))
        cuts[key[0] if len(key) == 1 else key] = cut
        if len(key) == 1:
            weights[key[0]] = girth - sum(chords[e] for e in cut if e in chords)
    assert all(weight(weights, cut) == girth for cut in cuts.values())
    planar = nx.Graph()
    planar.add_edges_from(weights)
    assert nx.check_planarity(planar)[0]
    return cactus, rows, {e: w * scale for e, w in weights.items()}, cuts


def run_branching_example():
    tree = nx.Graph()
    main = [f"a{j}" for j in range(1, 16)]
    nx.add_path(tree, ["root", *main, "tipA"])
    chords = {}
    for r in range(3):
        side = [f"arm{r}.{j}" for j in range(5)]
        nx.add_path(tree, ["root", *side, f"tip{r}"])
        for u, v in zip(main[5 * r:5 * r + 5], side):
            chords[edge(u, v)] = 1
    _, weights, cuts = weighted_tree_model(tree, chords)
    graph = nx.Graph()
    graph.add_edges_from(weights)
    assert nx.check_planarity(graph)[0]
    for r in range(3):
        bundle = {e for e in chords if any(u.startswith(f"arm{r}.") for u in e)}
        profiles = {frozenset(cuts[edge(u, v)] & bundle) for u, v in zip(main, main[1:])}
        assert len(profiles) >= 5
    print("branching case: sensitive degree at least 3, each bundle has at least 5 profiles")


def run():
    for n in [11, 12, 24, 64]:
        positions, vertices, edges = check_context(n)
        print(f"context n={n}: {positions} forced positions, primal {vertices}v/{edges}e")
    rand = Random(731342)
    totals = Counter()
    for sample in range(96):
        rows_count = rand.choice([2, 3, 4])
        tree, rows = make_spine([rand.choice([2, 3, 4]) for _ in range(rows_count)])
        links = [edge(q, (q + 1) % rows_count) for q in range(rows_count)]
        possibilities = [
            edge(u, v) for p, q in set(links) for u in rows[p] for v in rows[q]
        ]
        possibilities.extend(edge(u, v) for row in rows for u, v in combinations(row, 2)
                             if not tree.has_edge(u, v))
        rand.shuffle(possibilities)
        graph = tree.copy()
        chords = {}
        for u, v in possibilities[:rand.randint(1, 18)]:
            graph.add_edge(u, v)
            if nx.check_planarity(graph)[0]:
                chords[u, v] = rand.randint(1, 20)
            else:
                graph.remove_edge(u, v)
        _, weights, cuts = weighted_tree_model(tree, chords, 2**128 if sample % 4 == 0 else 1)
        values = check_marking(tree, rows, weights, cuts)
        totals.update(dict(zip(["selections", "cycle_components", "third_path_flags",
                                "discarded_cuts"], values)))
    for sample in range(24):
        args = weighted_cactus_model(rand, rand.choice([3, 4]), rand.choice([2, 3]),
                                     2**128 if sample % 4 == 0 else 1)
        values = check_marking(*args)
        totals.update(dict(zip(["selections", "cycle_components", "third_path_flags",
                                "discarded_cuts"], values)))
    print("marking models: 120, including 24 with a long cactus cycle and 30 with 129+ bit weights")
    print(dict(totals))
    run_branching_example()
    print("verification: PASS")


if __name__ == "__main__":
    run()
