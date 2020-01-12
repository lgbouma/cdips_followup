left = range(-90, 90, 2)
right = range(-88, 92, 2)

lines = []
for l,r in zip(left, right):

    lstr = str(abs(l)).zfill(2)
    if str(l).startswith('-'):
        lstr += '_00S'
        lstr = lstr.replace('-', '')
    else:
        lstr += '_00N'

    rstr = str(abs(r)).zfill(2)
    if str(r).startswith('-'):
        rstr += '_00S'
        rstr = rstr.replace('-', '')
    else:
        rstr += '_00N'

    line = (
        'wget http://archive.stsci.edu/missions/tess/catalogs/tic_v8/tic_dec{lstr}__{rstr}.csv.gz\n'
    ).format(
        lstr = lstr,
        rstr = rstr
    )

    lines.append(line)

with open('wget_tic_v8.sh', 'w') as f:
    f.writelines(lines)
