"""Download the Maxar Open Data aerial tiles used for the base map and texture atlases into work/maxar.

Source: Maxar Open Data program, Hurricane Helene event (public S3 bucket, no credentials).
Licence: CC BY-NC 4.0, attribution to Maxar required. Standard library only.
Usage: python tools/fetch_imagery.py            (from the repository root)
"""
import hashlib, os, sys, urllib.request
BASE = 'https://maxar-opendata.s3.amazonaws.com/events/HurricaneHelene-Oct24/ard/17'
# (quadkey, acquisition date, catalog id, local suffix). The 2020 scenes are used where they cover gaps in 2019.
TILES = [
 ('031133302032','2019-11-09','103001009D337C00',''),
 ('031133302033','2019-11-09','103001009D337C00',''),
 ('031133302122','2019-11-09','103001009D337C00',''),
 ('031133302122','2020-02-09','10300100A12E7B00','_2020'),
 ('031133302210','2019-11-09','103001009D337C00',''),
 ('031133302211','2019-11-09','103001009D337C00',''),
 ('031133302212','2019-11-09','103001009D337C00',''),
 ('031133302213','2019-11-09','103001009D337C00',''),
 ('031133302300','2019-11-09','103001009D337C00',''),
 ('031133302300','2020-02-09','10300100A12E7B00','_2020'),
 ('031133302302','2019-11-09','103001009D337C00',''),
 ('031133302302','2020-02-09','10300100A12E7B00','_2020'),
]
SIZES = {'031133302032.tif': 51653464, '031133302033.tif': 70344236, '031133302122.tif': 26951549, '031133302122_2020.tif': 62338488, '031133302210.tif': 54778290, '031133302211.tif': 69104844, '031133302212.tif': 31550064, '031133302213.tif': 38752362, '031133302300.tif': 25428442, '031133302300_2020.tif': 63487063, '031133302302.tif': 14898043, '031133302302_2020.tif': 39102816}   # byte sizes of the files the model was built with
def main():
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'work', 'maxar')
    os.makedirs(out, exist_ok=True)
    for qk, date, cat, suf in TILES:
        dst = os.path.join(out, f'{qk}{suf}.tif'); url = f'{BASE}/{qk}/{date}/{cat}-visual.tif'
        if os.path.exists(dst) and os.path.getsize(dst) == SIZES.get(os.path.basename(dst), -1):
            print('have', os.path.basename(dst)); continue
        print('get ', os.path.basename(dst), url, flush=True)
        tmp = dst + '.part'
        with urllib.request.urlopen(url, timeout=120) as r, open(tmp, 'wb') as f:
            while True:
                b = r.read(1 << 20)
                if not b: break
                f.write(b)
        os.replace(tmp, dst)
    for qk, date, cat, suf in TILES:
        p = os.path.join(out, f'{qk}{suf}.tif'); n = os.path.getsize(p)
        ok = SIZES.get(os.path.basename(p)) in (None, n)
        print(('ok  ' if ok else 'SIZE MISMATCH '), os.path.basename(p), n)
        if not ok: sys.exit(1)
if __name__ == '__main__':
    main()
