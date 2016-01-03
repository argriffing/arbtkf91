from __future__ import division, print_function

import sys
import json

if __name__ == '__main__':
    print('arbtkf91-check....', sep='', end='')
    arr = []
    for filename in sys.argv[1:]:
        with open(filename) as fin:
            arr.append(json.load(fin))
    first = arr[0]
    rest = arr[1:]
    for rhs in rest:
        if first != rhs:
            exit(-1)
            print('FAIL')
    print('PASS')
    exit(0)
