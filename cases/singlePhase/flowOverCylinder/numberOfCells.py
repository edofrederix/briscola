import sys
import math as m

# Given a segment length, mesh density and grading factor, compute the number of
# cells that is the closest base 1 or 3 power of two

if len(sys.argv) != 4:
    print("Invalid number of args")
    sys.exit(0)

L = float(sys.argv[1])          # Segment length
rho = float(sys.argv[2])        # Mesh density
g = float(sys.argv[3])          # Grading factor

if g == 1.0 or g/rho >= L:

    N = max(int(round(rho*L)),1)

else:

    Q = (L-1.0/rho)/(L-g/rho)
    N = max(int(round(m.log(g)/m.log(Q)+1)),1)

opts = [
    2**m.floor(m.log(N)/m.log(2)),
    2**m.ceil(m.log(N)/m.log(2)),
    3*2**m.floor(m.log(N/3)/m.log(2)),
    3*2**m.ceil(m.log(N/3)/m.log(2))
]

closest = N*10

for j,opt in enumerate(opts):
    if abs(N-opt) < abs(N-closest):
        closest = opt

print(closest)
