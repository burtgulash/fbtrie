#!/bin/sh

QUERY=kokot
K=2
TRIE=fbtrie

cat /usr/share/dict/cracklib-small | python3 ./fbtrie.py "$QUERY" $K $TRIE
