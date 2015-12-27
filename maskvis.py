import sys

import numpy as np

from PIL import Image

CRUMB_TOP  = 0x01
CRUMB_DIAG = 0x02
CRUMB_LEFT = 0x04
CRUMB_MASK = 0x08

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
            red = BRIGHT if (d & CRUMB_MASK) else DIM
            image[i, j, 2:, 2:, IRED] = red

            # top connection
            if i > 0:
                if data[i, j] & CRUMB_MASK and data[i-1, j] & CRUMB_MASK:
                    blue = BRIGHT if (d & CRUMB_TOP) else DIM
                    image[i, j, :2, 3, IBLUE] = blue

            # diag connection
            if i > 0 and j > 0:
                if data[i, j] & CRUMB_MASK and data[i-1, j-1] & CRUMB_MASK:
                    blue = BRIGHT if (d & CRUMB_DIAG) else DIM
                    image[i, j, 0, 0, IBLUE] = blue
                    image[i, j, 1, 1, IBLUE] = blue

            # left connection
            if j > 0:
                if data[i, j] & CRUMB_MASK and data[i, j-1] & CRUMB_MASK:
                    blue = BRIGHT if (d & CRUMB_LEFT) else DIM
                    image[i, j, 3, :2, IBLUE] = blue

    print 'yielding rgb triples...'
    for i in range(nr):
        print i
        for k in range(5):
            for j in range(nc):
                for l in range(5):
                    if j == 0 and l < 2: continue
                    if i == 0 and k < 2: continue
                    #yield tuple(image[i, j, k, l])
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


