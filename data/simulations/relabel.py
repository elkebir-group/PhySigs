#!/usr/bin/python
import sys

children = {}

with open(sys.argv[1]) as f:
    for line in f:
        s = line.rstrip("\n").split(" ")
        if s[0] not in children:
            children[s[0]] = []
        children[s[0]].append(s[1])

def relabel(children, v, old2new):
    if v not in old2new:
        old2new[v] = len(old2new)
        if v in children:
            for w in children[v]:
                relabel(children, w, old2new)

old2new = {}
relabel(children, "GL", old2new)

for v in children:
    if v != "GL":
        for w in children[v]:
            print old2new[v], old2new[w]
