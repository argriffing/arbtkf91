import sys

import numpy as np

from PIL import Image

# /* traceback directions for the actual alignment */
CRUMB_TOP  = 0x01
CRUMB_DIAG = 0x02
CRUMB_LEFT = 0x04

# /* traceback directions for information propagation */
CRUMB_DIAG2 = 0x08
CRUMB_LEFT2 = 0x10

CRUMB_WANT2 = 0x20
CRUMB_WANT3 = 0x40
CRUMB_CONTENDER = 0x80


IRED = 0
IGREEN = 1
IBLUE = 2

BRIGHT = 255
DIM = 63

def get_shape(rcv_triples):
    nr = max(r for r, c, v in rcv_triples) + 1
    nc = max(c for r, c, v in rcv_triples) + 1
    return nr, nc

def gen_rgb_triples(rcv_triples, nr, nc):
    data = np.zeros((nr, nc), dtype=int)
    for r, c, v in rcv_triples:
        data[r, c] = v
    image = np.zeros((nr, nc, 5, 5, 3), dtype=int)
    print 'creating the array...'
    for i in range(nr):
        print i
        for j in range(nc):
            d = data[i, j]

            # cell color
            if d & CRUMB_CONTENDER:
                image[i, j, 2:, 2:, IRED] = BRIGHT
            elif d & CRUMB_WANT3:
                image[i, j, 2:, 2:, IBLUE] = BRIGHT
            elif d & CRUMB_WANT2:
                image[i, j, 2:, 2:, IGREEN] = BRIGHT
            else:
                image[i, j, 2:, 2:, IGREEN] = DIM

            # connections among potentially informative nodes
            #
            # diag connection
            if i > 0 and j > 0:
                if data[i, j] & (CRUMB_WANT2 | CRUMB_WANT3) and data[i-1, j-1] & (CRUMB_WANT2 | CRUMB_WANT3):
                    if d & CRUMB_DIAG2:
                        image[i, j, 0, 0, IBLUE] = BRIGHT
                        image[i, j, 1, 1, IBLUE] = BRIGHT
            #
            # left connection
            if j > 0:
                if data[i, j] & (CRUMB_WANT2 | CRUMB_WANT3) and data[i, j-1] & (CRUMB_WANT2 | CRUMB_WANT3):
                    if d & CRUMB_LEFT2:
                        image[i, j, 3, :2, IBLUE] = BRIGHT

            # connections among potentially visited nodes
            #
            # top connection
            if i > 0:
                if data[i, j] & CRUMB_CONTENDER and data[i-1, j] & CRUMB_CONTENDER:
                    if d & CRUMB_TOP:
                        image[i, j, :2, 3, IRED] = BRIGHT
            #
            # diag connection
            if i > 0 and j > 0:
                if data[i, j] & CRUMB_CONTENDER and data[i-1, j-1] & CRUMB_CONTENDER:
                    if d & CRUMB_DIAG:
                        image[i, j, 0, 0, IRED] = BRIGHT
                        image[i, j, 1, 1, IRED] = BRIGHT
            #
            # left connection
            if j > 0:
                if data[i, j] & CRUMB_CONTENDER and data[i, j-1] & CRUMB_CONTENDER:
                    if d & CRUMB_LEFT:
                        image[i, j, 3, :2, IRED] = BRIGHT

    print 'yielding rgb triples...'
    for i in range(nr):
        print i
        for k in range(5):
            for j in range(nc):
                for l in range(5):
                    if j == 0 and l < 2: continue
                    if i == 0 and k < 2: continue
                    yield tuple(image[i, j, k, l].tolist())



def main():
    lines = sys.stdin.readlines()
    rcv_triples = [tuple(int(s) for s in line.split()) for line in lines]
    nr, nc = get_shape(rcv_triples)
    rgb_triples = list(gen_rgb_triples(rcv_triples, nr, nc))
    width = nc * 5 - 2
    height = nr * 5 - 2
    im = Image.new("RGB", (width, height))
    im.putdata(rgb_triples)
    im.save('trace.png', 'png')

main()
