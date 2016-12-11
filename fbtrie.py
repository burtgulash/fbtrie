#!/usr/bin/env python3

"""
Implementation of FB-Trie data structure
for fast levenshtein dictionary search.

Trie represents classic Trie.  FB-Trie consists of two tries. F stands for
Forward trie and is the same thing as 'Trie'. B means Backward trie, because it
is a trie which just stores strings in reverse (eg. F trie would contain 'abc',
while B trie would have 'cba').

The search procedure consists of two phases - first in F trie and then the same
thing in B trie. Suppose you have limit of edit distance 'k'. Instead of
dispatching the trie search recursion with limit 'k', you dispatch recursion
with limit 'k / 2' instead. Once you hit a node which represents roughly half
of the search pattern, you increase the limit back to 'k'.

FB-Trie is fast for fuzzy search, because it reduces excessive branching in the
initial levels of classic trie. The pruned branches are then picked up in
second pass using backward trie, which has nodes representing initial levels of
forward trie at the end levels.
"""


class Trie:

    def __init__(self):
        self.root = {}

    def insert(self, word):
        word = word + "\0"
        node = self.root
        for c in word:
            if c not in node:
                node[c] = {}
            node = node[c]

    def print(self):
        self.print_(self.root, "", "")

    def print_(self, node, sofar, before):
        num_ch = len(node)
        if len(node) > 1:
            print(before, sofar, "/", sep="")
        before += " " if num_ch > 1 else ""

        for c, tree in node.items():
            if c == "\0":
                print(before, sofar, sep="")
            else:
                self.print_(tree, sofar + c, before)

    def traverse(self, d, node, sofar):
        for c, child in node.items():
            if c == "\0":
                yield sofar, d
            else:
                yield from self.traverse(d, child, sofar + c)

    def fuzzy(self, query, k, k_fn=None, prefix=False):
        n = len(query)
        tab = [[0] * (n + 1) for _ in range(n + k + 1)]
        for j in range(n + 1):
            tab[0][j] = j

        if not k_fn:
            k_fn = lambda i: k
        yield from self.fuzzy_(k_fn, tab, self.root, 1, "", query, prefix)

    def fuzzy_(self, k_fn, tab, node, i, sofar, query, prefix=False):
        k = k_fn(i)

        if prefix and i >= len(query) + 1:
            d = tab[i - 1][len(query)]
            if d <= k:
                # yield sofar, d
                yield from self.traverse(d, node, sofar)
            return

        # Can't be more than length of query + k insertions. Don't allow
        # children
        if i >= len(tab):
            return

        for c, child in node.items():
            d = tab[i - 1][len(query)]
            if c == "\0" and d <= k:
                yield sofar, d
                continue

            new = sofar + c
            tab[i][0] = i
            for j in range(1, len(tab[i])):
                sub_cost = 1 if c != query[j - 1] else 0
                sub = tab[i - 1][j - 1] + sub_cost
                insert = tab[i - 1][j] + 1
                delete = tab[i][j - 1] + 1
                tab[i][j] = min(sub, insert, delete)

            smallest = min(tab[i])
            if smallest <= k:
                yield from self.fuzzy_(k_fn, tab, child,
                                       i + 1, new, query, prefix)


####################
class FBTrie:

    def __init__(self):
        self.f = Trie()
        self.b = Trie()

    def insert(self, word):
        self.f.insert(word)
        self.b.insert(word[::-1])

    def print(self):
        self.f.print()
        self.b.print()

    def fuzzy(self, query, k):
        n = len(query)
        q1, q2 = query[:n // 2], query[n // 2:]

        # Boytsov 2011:
        # Note that there are ⌈k/2⌉ queries of the first type
        # and ⌊k/2⌋+ 1 queries of the second type
        #
        # NOTE: why that many times?? only the last one seems to suffice

        # k1 = (k - 1) // 2 + 1
        k1 = (k - 1) // 2
        yield from self.trie_fuzzy(self.f, k1, k, q1, query)

        # k2 = k // 2 + 1
        k2 = k // 2
        yield from ((word[::-1], d) for word, d
                    in self.trie_fuzzy(self.b, k2, k, q2[::-1], query[::-1]))

    def trie_fuzzy(self, trie, subk, k, subq, query):
        n = len(query)
        tab = [[0] * (n + 1) for _ in range(n + k + 1)]
        for j in range(n + 1):
            tab[0][j] = j

        yield from self.trie_fuzzy_(subk, k, tab, trie.root, 1,
                                    "", subq, query, 1)

    def trie_fuzzy_(self, subk, k, tab, node, i, sofar, subq, query, phase):
        if i >= len(tab):
            return

        for c, child in sorted(node.items(), key=lambda x: x[0]):
            d = tab[i - 1][len(query)]
            if c == "\0" and d <= k:
                yield sofar, d
                continue

            new = sofar + c
            tab[i][0] = i
            for j in range(1, len(tab[i])):
                sub_cost = 1 if c != query[j - 1] else 0
                sub = tab[i - 1][j - 1] + sub_cost
                insert = tab[i - 1][j] + 1
                delete = tab[i][j - 1] + 1
                tab[i][j] = min(sub, insert, delete)

            if phase == 1 and tab[i][len(subq)] <= subk:
                yield from self.trie_fuzzy_(k, k, tab, child, i + 1,
                                            new, subq, query, 2)
                continue

            smallest = min(tab[i])
            if smallest <= subk:
                yield from self.trie_fuzzy_(subk, k, tab, child, i + 1,
                                            new, subq, query, phase)


if __name__ == "__main__":
    import sys
    import time

    if len(sys.argv) > 4:
        usage = "python trie.py query K [trie|fbtrie]"
        print(usage, file=sys.stderr)
        sys.exit(1)

    query = sys.argv[1]
    k = int(sys.argv[2])

    typ = "trie"
    if len(sys.argv) == 4:
        typ = sys.argv[3].lower()

    if typ == "trie":
        t = Trie()
    elif typ == "fbtrie":
        t = FBTrie()
    else:
        print("Unknown trie type. Using default 'trie'", file=sys.stderr)
        t = Trie()

    print("Reading from stdin...", file=sys.stderr)
    for line in sys.stdin:
        line = line.strip()
        t.insert(line)

    print("Processing query", query, k, file=sys.stderr)

    start = time.time()
    result = sorted(set((d, word) for word, d in t.fuzzy(query, k)))
    end = time.time()

    for word, d in result:
        print(d, word)

    print("RESULT: {}~{}: [{} found] in {:.4f}ms using [{}]".format(
        query, k, len(result), (end - start) * 1000, t.__class__.__name__),
        file=sys.stderr
    )
