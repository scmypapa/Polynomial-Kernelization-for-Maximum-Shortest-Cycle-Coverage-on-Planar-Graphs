"""有限复核。需要 NetworkX 及同目录的 verify_shortest_cycle_coverage_followup.py。"""

from itertools import combinations
from dataclasses import dataclass
from fractions import Fraction as Q
from bisect import bisect_left, bisect_right
from collections import Counter, defaultdict
from random import Random

import networkx as nx

from verify_shortest_cycle_coverage_followup import (
    edge, make_spine, weighted_cactus_model, weighted_tree_model,
)


def cyclic_after(values, point):
    i = values.index(point)
    return values[i + 1:] + values[:i]


def split_rotation(weights, paths, manual=None):
    graph = nx.Graph()
    graph.add_edges_from(sorted(weights))
    planar, rotation = nx.check_planarity(graph)
    assert planar
    if manual is not None:
        rotation = nx.PlanarEmbedding()
        rotation.set_data(manual)
        rotation.check_structure()
    owner = {}
    f = set()
    side = {}
    coord = {}
    geometry = {}
    for p, path in enumerate(paths):
        side[p] = ([], [])
        rows = []
        for j, v in enumerate(path):
            owner[v] = p
            cw = list(rotation.neighbors_cw_order(v))
            if len(path) == 1:
                clockwise, anticlockwise = cw, []
            elif j == 0:
                clockwise, anticlockwise = [], cyclic_after(cw, path[1])
                f.add(edge(v, path[1]))
            elif j == len(path) - 1:
                clockwise, anticlockwise = cyclic_after(cw, path[j - 1]), []
            else:
                a, b = path[j - 1], path[j + 1]
                seq = cyclic_after(cw, a)
                cut = seq.index(b)
                clockwise, anticlockwise = seq[:cut], seq[cut + 1:]
                f.add(edge(v, b))
            events = [(v, w) for w in clockwise + list(reversed(anticlockwise))]
            for i, dart in enumerate(events, 1):
                coord[dart] = Q(j) + Q(i, 3 * (len(events) + 1))
            rows.append(([(v, w) for w in clockwise], [(v, w) for w in anticlockwise]))
        side[p][0].extend(d for first, _ in rows for d in first)
        side[p][1].extend(d for _, second in reversed(rows) for d in second)
        for s, group in enumerate(side[p]):
            virtual = ((p, s), (p, 1 - s))
            geometry[p, s] = [virtual, *group]
    assert len(coord) == 2 * (len(weights) - len(f))
    return geometry, coord, owner, f


def mate(dart):
    return dart[1], dart[0]


def faces(rotation):
    succ = {d: group[(i + 1) % len(group)] for group in rotation.values()
            for i, d in enumerate(group)}
    seen = set()
    for start in succ:
        if start in seen:
            continue
        face = []
        d = start
        while d not in seen:
            face.append(d)
            seen.add(d)
            d = succ[mate(d)]
        assert d == start
        yield face


@dataclass(frozen=True)
class Port:
    track: int
    side: int
    ends: tuple
    darts: tuple


def extract_bands(rotation, xy, weights, f, paths):
    side_id = {d: u for u, darts in rotation.items() for d in darts}
    trivial = set()
    for u, darts in rotation.items():
        indices = {d: i for i, d in enumerate(darts)}
        for a in darts[1:]:
            b = mate(a)
            if b in indices and indices[a] < indices[b]:
                if all(side_id[mate(d)] == u for d in darts[indices[a] + 1:indices[b]]):
                    trivial.add(edge(*a))
    boundary = {u: [group[0], *[d for d in group[1:] if edge(*d) not in trivial]]
                for u, group in rotation.items()}
    dartgraph = nx.Graph()
    pure = {d for d in xy if edge(*d) not in trivial}
    dartgraph.add_nodes_from(sorted(pure))
    fs = list(faces(boundary))
    edges_count = sum(map(len, boundary.values())) // 2
    assert 2 * len(paths) - edges_count + len(fs) == 2
    for ring in fs:
        if len(ring) == 2 and all(d in pure for d in ring):
            a, b = ring
            assert b != mate(a)
            dartgraph.add_edge(a, mate(b))
            dartgraph.add_edge(mate(a), b)
    components = list(nx.connected_components(dartgraph))
    part = {d: j for j, comp in enumerate(components) for d in comp}
    ordered_events = {u: sorted(xy[d] for d in group if d in pure)
                      for u, group in rotation.items()}
    bands, matched = [], set()
    for j, comp in enumerate(components):
        if j in matched:
            continue
        other, = {part[mate(d)] for d in comp}
        assert other != j
        matched.update([j, other])
        ports = []
        for group in [comp, components[other]]:
            u, = {side_id[d] for d in group}
            dst = sorted(group, key=xy.get)
            events = [Q(-1)] + ordered_events[u] + [Q(len(paths[u[0]]))]
            lo, hi = events.index(xy[dst[0]]), events.index(xy[dst[-1]])
            assert events[lo:hi + 1] == [xy[d] for d in dst]
            ends = ((2 * events[lo] + events[lo - 1]) / 3,
                    (2 * events[hi] + events[hi + 1]) / 3)
            ports.append(Port(*u, ends, tuple(dst)))
        assert set(map(mate, ports[0].darts)) == set(ports[1].darts)
        bands.append(tuple(ports))
    assert len(bands) <= 8 * len(paths)
    for e in trivial:
        p, s = side_id[e]
        a, b = sorted([xy[e], xy[mate(e)]])
        assert not any(a < z < b for z in ordered_events[p, s])
    return bands, trivial, ordered_events, side_id


def map_pair(a, b, xy):
    xs = [a.ends[0], *[xy[d] for d in a.darts], a.ends[1]]
    ys = [xy[mate(d)] for d in a.darts]
    if ys == sorted(ys):
        ys = [b.ends[0], *ys, b.ends[1]]
    else:
        assert ys == sorted(ys, reverse=True)
        ys = [b.ends[1], *ys, b.ends[0]]
    return xs, ys


def interpolate(xs, ys, value):
    if value == xs[-1]:
        return ys[-1]
    i = bisect_left(xs, value, lo=1)
    return ys[i - 1] + (value - xs[i - 1]) * (ys[i] - ys[i - 1]) / (xs[i] - xs[i - 1])


class Layout:
    def __init__(self, bands, xy, paths, limit, events):
        self.tips = [(Q(-1), Q(len(line))) for line in paths]
        self.events = {u: sorted([*self.tips[u[0]], *ev]) for u, ev in events.items()}
        self.ports = defaultdict(list)
        self.maps = {}
        frontier = {(p, x) for p, ends in enumerate(self.tips) for x in ends}
        for b, ports in enumerate(bands):
            for i, (a, z) in enumerate([ports, tuple(reversed(ports))]):
                frontier.update((a.track, x) for x in a.ends)
                self.ports[a.track, a.side].append((b, i, a))
                self.maps[b, i] = (z.track, z.side, *map_pair(a, z, xy))
        self.level = {m: 0 for m in frontier}
        for r in range(1, limit + 1):
            following = set()
            for p, x in frontier:
                for s in range(2):
                    data = self.step(p, s, x)
                    if data is not None:
                        q, _, y, _ = data
                        if (q, y) not in self.level:
                            following.add((q, y))
            if not following:
                break
            self.level.update((point, r) for point in following)
            frontier = following
        assert len(self.level) <= 30 * (limit + 1) * len(paths)

    def step(self, p, s, x):
        found = [(b, i) for b, i, a in self.ports[p, s] if a.ends[0] <= x <= a.ends[1]]
        if not found:
            return None
        key, = found
        q, side, xs, ys = self.maps[key]
        return q, side, interpolate(xs, ys, x), key

    def marks(self, depth):
        tracks = defaultdict(list)
        for (p, x), first in self.level.items():
            if first <= depth:
                tracks[p].append(x)
        return {p: sorted(seq) for p, seq in tracks.items()}

    def gap(self, p, s, x):
        array = self.events[p, s]
        i = bisect_left(array, x)
        if i == len(array) or i == 0 or x == array[i]:
            return None
        return array[i - 1], array[i]

    def hit(self, p, s, x, marks):
        a, b = self.gap(p, s, x)
        i = bisect_right(marks[p], a)
        return i < len(marks[p]) and marks[p][i] < b

    def add_gap_ends(self, marks):
        result = {}
        for p, seq in marks.items():
            additions = set(seq)
            for x in seq:
                for s in range(2):
                    gap = self.gap(p, s, x)
                    if gap:
                        additions.update(gap)
            result[p] = sorted(additions)
        return result


def label(marks, p, x):
    values = marks[p]
    i = bisect_left(values, x)
    if x == values[i]:
        return p, x, x
    return p, values[i - 1], values[i]


def inside(cell, x):
    return cell[1] < x < cell[2]


def model_tree(rand):
    n = rand.randint(3, 6)
    lengths = [rand.choice([2, 3, 4]) for _ in range(n)]
    base, rows = make_spine(lengths)
    graph = base.copy()
    chords = {}
    possible = list(combinations(base.nodes, 2))
    rand.shuffle(possible)
    for u, v in possible:
        if len(chords) >= 2 * n + 6:
            break
        if graph.has_edge(u, v):
            continue
        graph.add_edge(u, v)
        if nx.check_planarity(graph)[0]:
            chords[edge(u, v)] = rand.randint(1, 9)
        else:
            graph.remove_edge(u, v)
    _, weights, cuts = weighted_tree_model(base, chords, scale=2**140 if n % 2 else 1)
    return base, rows, weights, cuts


def evaluate(instance, rand, k, repeats, prescribed=(), fixed_mode=False):
    base, paths, weights, cuts, *embedding = instance
    paths = paths + [[v] for v in base if v not in {u for line in paths for u in line}]
    rot, xy, owner, f = split_rotation(weights, paths, *embedding)
    bands, trivial, events, side_id = extract_bands(rot, xy, weights, f, paths)
    d = k + 2
    limit = (6 * k * k + 1) * d
    layout = Layout(bands, xy, paths, limit + 1, events)
    all_marks = layout.marks(limit + 1)
    cut_positions = {}
    for p, path in enumerate(paths):
        for j in range(1, len(path)):
            x = Q(j) - Q(1, 2)
            while x in all_marks[p]:
                x += Q(1, 10 * (len(all_marks[p]) + 1))
            assert j - Q(1, 2) <= x < j - Q(2, 5)
            cut_positions[edge(path[j - 1], path[j])] = (p, x)
    if fixed_mode:
        for (band, i), (q, s, xs, ys) in layout.maps.items():
            p = bands[band][i].track
            if p != q:
                continue
            for a, b, c, e in zip(xs, xs[1:], ys, ys[1:]):
                if (c - a) * (e - b) >= 0:
                    continue
                x = (a * (e - b) - b * (c - a)) / ((e - b) - (c - a))
                j = int(x) + 1
                line = paths[p]
                if Q(j - 1) + Q(1, 3) < x < j and x not in all_marks[p]:
                    focus = edge(line[j - 1], line[j])
                    cut_positions[focus] = (p, x)
                    prescribed = [[focus, *[d for d in cuts if d not in cut_positions][:k - 1]]]
                    break
            if prescribed:
                break
        assert prescribed
    stats = Counter()
    samples = list(prescribed) + [rand.sample(sorted(cuts, key=repr), k) for _ in range(repeats)]
    for selected in samples:
        assert len(selected) == k and len(set(selected)) == k
        ordinary = [e for e in selected if e in cut_positions]
        if not ordinary:
            continue
        points = [cut_positions[e] for e in ordinary]
        compares = []
        for (p, a), (q, b) in combinations(points, 2):
            if p == q:
                compares.append((p, a, b))
        for p, a in points:
            for side in range(2):
                arrival = layout.step(p, side, a)
                if arrival:
                    q, other, x, (b, i) = arrival
                    for r, y in points:
                        if r == q and any(bb == b and ii == 1 - i and z.ends[0] < y < z.ends[1]
                                          for bb, ii, z in layout.ports[q, other]):
                            compares.append((q, x, y))
        cuts_at = []
        for p, x, y in compares:
            values = all_marks[p]
            a = bisect_right(values, min(x, y))
            b = bisect_left(values, max(x, y))
            if a < b:
                cuts_at.append(min(layout.level[p, m] for m in values[a:b]))
        stability = {((depth - 1) // d) * d for depth in cuts_at if depth > 0}
        t = next(r for r in range(0, limit, d) if r not in stability)
        raw, fine, deep = [layout.marks(r) for r in [t, t + 1, t + d]]
        extra = raw
        for s in range(2 * k + 1):
            plus = layout.add_gap_ends(extra)
            if all(layout.hit(p, side, x, extra) == layout.hit(p, side, x, plus)
                   for p, x in points for side in range(2)):
                break
            extra = plus
        else:
            raise AssertionError("no stable local round")
        stats.update(check_selection(layout, bands, xy, side_id, trivial,
                                     weights, cuts, selected, cut_positions, raw, fine, deep,
                                     extra, plus, rand))
    return stats


def check_selection(layout, bands, xy, sides, trivial, weights,
                    cuts, selected, positions, raw, fine, deep, old_gaps, refine,
                    rand):
    selected = set(selected)
    chosen = selected.intersection(positions)
    exceptions = frozenset().union(*(cuts[e] for e in selected - chosen))
    fine_by_raw = {}
    for e in chosen:
        p, x = positions[e]
        r, i = label(raw, p, x), label(fine, p, x)
        assert fine_by_raw.get(r, i) == i
        fine_by_raw[r] = i
    destination = {}
    self_regions = set()
    for r, i in fine_by_raw.items():
        p = r[0]
        middle = (i[1] + i[2]) / 2
        for s in range(2):
            arrival = layout.step(p, s, middle)
            if arrival is None:
                continue
            q, side, y, key = arrival
            target = label(raw, q, y)
            destination[r, s] = target, key
        if all((r, s) in destination and destination[r, s][0] == r for s in range(2)):
            a, b = (destination[r, s][1] for s in range(2))
            if a[0] == b[0] and a[1] != b[1]:
                self_regions.add(r)
    chosen_by_block = defaultdict(list)
    eligible_by_block = defaultdict(list)
    for e, (p, x) in positions.items():
        r, i = label(raw, p, x), label(fine, p, x)
        if fine_by_raw.get(r) != i:
            continue
        if r in self_regions:
            block = r, r
        else:
            inner = label(refine, p, x)
            block = r, inner
            if any(layout.hit(p, s, x, refine) and not layout.hit(p, s, x, old_gaps)
                   for s in range(2)):
                assert e not in chosen
                continue
        eligible_by_block[block].append(e)
        if e in chosen:
            chosen_by_block[block].append(e)
    assert sum(map(len, chosen_by_block.values())) == len(chosen)

    def ordered(candidates):
        return sorted(candidates, key=lambda e: positions[e][1])

    def supported(e, d):
        return e in cuts[d]

    def profile_closed(block, side):
        r, inner = block
        return r not in self_regions and all(
            layout.hit(r[0], side, positions[d][1], old_gaps)
            for d in eligible_by_block[block]
        )

    variants = [set(chosen)]
    for _ in range(12):
        variants.append(set().union(*(set(rand.sample(eligible_by_block[block], len(c)))
                                      for block, c in chosen_by_block.items())))
    replacement_by_block = [{block: set(eligible_by_block[block]).intersection(variant)
                             for block in chosen_by_block} for variant in variants]
    raw_blocks = {
        r: sorted([block for block in chosen_by_block if block[0] == r], key=lambda x: x[1][1])
        for r in fine_by_raw
    }
    raw_domains = {r: ordered(cut for b in blocks for cut in eligible_by_block[b])
                   for r, blocks in raw_blocks.items()}
    links = nx.MultiGraph()
    links.add_nodes_from(chosen_by_block)
    packet_edges = defaultdict(list)
    for band, (a, b) in enumerate(bands):
        for d1 in a.darts:
            d2 = mate(d1)
            r = label(raw, a.track, xy[d1])
            u = label(raw, b.track, xy[d2])
            packet_edges[band, r, u].append((edge(*d1), d1, d2))
    essential_checks = Counter()
    variant_covers = [exceptions | frozenset().union(*(cuts[d] for d in v)) for v in variants]
    for (band, r, u), packet in packet_edges.items():
        a, b = bands[band]
        edges = [e for e, _, _ in packet]
        background_options = []
        for d in chosen:
            p, x = positions[d]
            z = label(raw, p, x)
            if z not in {r, u}:
                background_options.extend(raw_domains[z])
        backgrounds = {bool(e in exceptions or any(supported(e, d) for d in background_options))
                       for e in edges}
        assert len(backgrounds) == 1
        shadowed, = backgrounds
        if shadowed:
            assert all(set(edges) <= cover for cover in variant_covers)
            essential_checks["shadow_packets"] += 1
            continue
        if r == u:
            if r not in fine_by_raw:
                assert all(not (set(edges) & cover) for cover in variant_covers)
                continue
            if r in self_regions:
                relevant, = raw_blocks[r]
                for v, cover in zip(replacement_by_block, variant_covers):
                    assert all((e in cover) == any(supported(e, d) for d in v[relevant])
                               for e in edges)
                essential_checks["internal_self_packets"] += 1
                continue
            by_side = {}
            for e, d1, d2 in packet:
                locs = [xy[d1], xy[d2]]
                assert sum(inside(fine_by_raw[r], x) for x in locs) <= 1
                varying = [j for j, x in enumerate(locs) if inside(fine_by_raw[r], x)]
                if not varying:
                    values = {supported(e, d) for d in raw_domains[r]}
                    assert len(values) == 1
                    val, = values
                    assert all((e in c) == val for c in variant_covers)
                    continue
                j, = varying
                other = locs[1 - j]
                left = other < fine_by_raw[r][1]
                assert j not in by_side or by_side[j] == left
                by_side[j] = left
                head = raw_blocks[r][0 if left else -1]
                for v, cover in zip(replacement_by_block, variant_covers):
                    point = ordered(v[head])[0 if left else -1]
                    assert (e in cover) == supported(e, point)
                essential_checks["internal_one_sided_edges"] += 1
            continue
        heads = {}
        for z in [r, u]:
            if z not in raw_domains:
                continue
            sequence = raw_domains[z]
            falls = rises = False
            for e in edges:
                bits = [supported(e, d) for d in sequence]
                falls |= any(x and not y for x, y in zip(bits, bits[1:]))
                rises |= any(y and not x for x, y in zip(bits, bits[1:]))
            assert not (falls and rises)
            last = rises
            block = raw_blocks[z][-1 if last else 0]
            heads[z] = block, last
        for v, cover in zip(replacement_by_block, variant_covers):
            samples = [ordered(v[head])[-1 if last else 0] for head, last in heads.values()]
            assert all((e in cover) == any(supported(e, x) for x in samples) for e in edges)
        if len(heads) != 2:
            continue
        h1, _ = heads[r]
        h2, _ = heads[u]
        mutual = destination.get((r, a.side), (None,))[0] == u and \
            destination.get((u, b.side), (None,))[0] == r
        if r in self_regions or u in self_regions:
            assert not mutual
        elif mutual and not profile_closed(h1, a.side) and not profile_closed(h2, b.side):
            links.add_edge(h1, h2, key=(band, r, u), sides={h1: a.side, h2: b.side})
            continue
        for e in edges:
            assert any(len({supported(e, d) for d in eligible_by_block[h]}) == 1
                       for h in [h1, h2])
        essential_checks["separable_two_end_packets"] += 1
    for block in links:
        at = [attrs["sides"][block] for _, _, attrs in links.edges(block, data=True)]
        assert len(at) == len(set(at)) <= 2
        if len(at) == 2:
            assert not any(profile_closed(block, s) for s in range(2))
    for part in nx.connected_components(links):
        if len(part) == 1:
            continue
        current = next((v for v in part if links.degree(v) == 1), next(iter(part)))
        any_selected = chosen_by_block[current][0]
        p, actual = positions[any_selected]
        _, a, b = label(deep, p, actual)
        point = (a + b) / 2
        taken = set()
        arrival_side = None
        while True:
            options = [(dest, key, attrs) for _, dest, key, attrs in links.edges(current, keys=True, data=True)
                       if key not in taken]
            if not options:
                break
            dest, key, attrs = options[0]
            taken.add(key)
            depart = attrs["sides"][current]
            assert depart != arrival_side
            q, arrival_side, point, (band, idx) = layout.step(p, depart, point)
            assert band == key[0] and arrival_side == attrs["sides"][dest]
            assert label(raw, q, point) == dest[0]
            assert label(fine, q, point) == fine_by_raw[dest[0]]
            current, p = dest, q
        assert len(taken) == links.subgraph(part).number_of_edges()
        essential_checks["replayed_components"] += 1
        if len(taken) == len(part):
            essential_checks["replayed_cycles"] += 1
    for e in trivial:
        lo, hi = sorted([xy[e], xy[mate(e)]])
        p, side = sides[e]
        involved = [bl for bl in chosen_by_block if bl[0][0] == p
                    and any(lo < positions[d][1] < hi for d in eligible_by_block[bl])]
        if len(involved) >= 2:
            assert all(bl[0] in self_regions or profile_closed(bl, side) for bl in involved)
    exact_in = Counter()
    for v in variants:
        direct = frozenset().union(*(cuts[e] for e in v))
        for p in range(len(layout.tips)):
            order = sorted((bl for bl in chosen_by_block if bl[0][0] == p),
                           key=lambda bl: bl[1][1])
            for s in range(2):
                es = {e for e in trivial if sides[e] == (p, s)}
                blocks = [[e for e in eligible_by_block[bl] if e in v] for bl in order]
                covered = [frozenset().union(*(cuts[e] for e in block)) & es for block in blocks]
                total = sum(weights[e] for c in covered for e in c)
                for j, (left, right) in enumerate(zip(covered, covered[1:])):
                    (r, am), (z, bm) = order[j:j + 2]
                    cut = am[2]
                    a = max(positions[e][1] for e in blocks[j])
                    b = min(positions[e][1] for e in blocks[j + 1])
                    assert a < cut < b
                    spans = {e for e in es if min(xy[e], xy[mate(e)]) < cut < max(xy[e], xy[mate(e)])}
                    loss = min(sum(weights[e] for e in spans & left if min(xy[e], xy[mate(e)]) < a < max(xy[e], xy[mate(e)])),
                               sum(weights[e] for e in spans & right if min(xy[e], xy[mate(e)]) < b < max(xy[e], xy[mate(e)])))
                    assert loss == sum(weights[e] for e in left & right)
                    total -= loss
                    exact_in["overlap_checks"] += 1
                assert total == sum(weights[e] for e in es & direct)
    return dict(exact_in + essential_checks, selections=1, replacements=len(variants) - 1,
                loop_bands=sum(a.track == b.track for a, b in bands),
                two_port_nodes=sum(links.degree(bl) == 2 for bl in links),
                closed=sum(profile_closed(bl, s) for bl in chosen_by_block for s in range(2)),
                self_regions=len(self_regions), links=links.number_of_edges())


def main():
    counts = Counter()
    rng = Random(98541)
    for i in range(96):
        fixture = (weighted_cactus_model(rng, i % 2 + 3, 3, 1)
                   if i % 4 == 0 else model_tree(rng))
        stats = evaluate(fixture, rng, rng.randint(3, 6), 12)
        counts.update(stats)
    for n in [16, 30, 54]:
        tree, rows = make_spine([n] * 3)
        chords = {edge(u, v): 1 for a, b in zip(rows, rows[1:]) for u, v in zip(a, b)}
        _, weights, cuts = weighted_tree_model(tree, chords)
        seeds = [[edge(row[i - 1], row[i]) for row in rows] for i in range(4, n - 2)]
        counts.update(evaluate((tree, rows, weights, cuts), rng, 3, 12, seeds))
    for m in [3, 4]:
        instance = weighted_cactus_model(rng, m, 24, 2**140)
        rows = instance[1]
        seeds = [[edge(row[i - 1], row[i]) for row in rows] for i in range(4, 21)]
        counts.update(evaluate(instance, rng, m, 12, seeds))
    for n in [14, 36]:
        base, rows = make_spine([n])
        row, = rows
        chords = {edge(row[i], row[i + 3]): 1 for i in range(1, n - 3)}
        _, weights, cuts = weighted_tree_model(base, chords)
        clockwise = {"start": [row[0]], "end": [row[-1]]}
        for i, u in enumerate(row):
            prev = "start" if i == 0 else row[i - 1]
            following = "end" if i == n else row[i + 1]
            ahead = [row[i + 3]] if i < n - 3 and edge(u, row[i + 3]) in chords else []
            behind = [row[i - 3]] if i >= 3 and edge(u, row[i - 3]) in chords else []
            clockwise[u] = [prev, *ahead, following, *behind]
        seeds = [[edge(row[i - 1], row[i]), edge(row[i], row[i + 1]),
                  edge(row[i + 1], row[i + 2])] for i in range(4, n - 4)]
        counts.update(evaluate((base, rows, weights, cuts, clockwise), rng, 3, 12, seeds))
    n = 35
    base, rows = make_spine([n])
    row, = rows
    arcs = [(i, i + 3) for i in range(1, 8)] + [(i, i - 2) for i in range(14, n)]
    chords = {edge(row[i], row[j]): 1 for i, j in arcs}
    _, weights, cuts = weighted_tree_model(base, chords)
    clockwise = {"start": [row[0]], "end": [row[-1]]}
    for i, v in enumerate(row):
        prev = "start" if i == 0 else row[i - 1]
        following = "end" if i == n else row[i + 1]
        above = [row[b] for a, b in arcs if a == i]
        below = [row[a] for a, b in arcs if b == i]
        clockwise[v] = [prev, *above, following, *below]
    counts.update(evaluate((base, rows, weights, cuts, clockwise), rng, 3, 3, fixed_mode=True))
    print(dict(counts))
    assert counts["two_port_nodes"] > 0 and counts["self_regions"] > 0
    print("geometry and inessential-overlap checks passed")


if __name__ == "__main__":
    main()
