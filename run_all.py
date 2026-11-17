"""Rebuild i26-connector-3d.html from the cached intermediates in work/ (about 10 minutes; scenery.py is the slow step).

Usage:  uv run python run_all.py            (from the repository root)
        uv run python run_all.py scenery    (start from a given step)
"""
import os, shutil, subprocess, sys, time
ROOT = os.path.dirname(os.path.abspath(__file__))
PY = os.path.join(ROOT, 'py')
STEPS = [('profiles2', ['profiles2.py', 'b0']), ('meshz', ['meshz.py']), ('walks', ['walks.py']),
         ('terrain', ['terrain.py']),
         ('routes', ['routes.py']), ('marks', ['marks.py']), ('scenery', ['scenery.py']),
         ('export', ['export.py']), ('build', [os.path.join('..', 'web', 'build.py')])]
def main():
    names = [n for n, _ in STEPS]; start = names.index(sys.argv[1]) if len(sys.argv) > 1 else 0
    env = dict(os.environ, MPLBACKEND='Agg')
    for n, args in STEPS[start:]:
        t = time.time(); print(f'== {n}', flush=True)
        subprocess.run([sys.executable] + args, cwd=PY, env=env, check=True)
        print(f'   {n} done in {time.time() - t:.0f}s', flush=True)
    out = os.path.join(ROOT, 'i26-connector-3d.html'); print('Output:', out)
    # a local web server folder (serve/, git-ignored) gets the new build too, so the page it serves is never stale
    serve = os.path.join(ROOT, 'serve')
    if os.path.isdir(serve) and STEPS[-1][0] in [n for n, _ in STEPS[start:]]:
        shutil.copyfile(out, os.path.join(serve, 'index.html')); print('Copied to', os.path.join(serve, 'index.html'))
if __name__ == '__main__':
    main()
