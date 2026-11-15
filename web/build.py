import json, base64, os, sys
HERE=os.path.dirname(os.path.abspath(__file__))
W=os.environ.get('I26_WORK') or os.path.join(os.path.dirname(HERE),'work')
S=W+'/scene'
man=json.load(open(S+'/manifest.json'))
tpl=open(HERE+'/src/index.html').read()
parts=[]
for name in man['bufs']:
    b=base64.b64encode(open(f'{S}/{name}.bin','rb').read()).decode()
    parts.append(f'<script type="application/octet-stream" id="buf_{name}">{b}</script>')
for name,info in man.get('imgs',{}).items():
    b=base64.b64encode(open(f'{S}/{info["file"]}','rb').read()).decode()
    parts.append(f'<script type="application/octet-stream" id="img_{name}" data-mime="{info.get("mime","image/jpeg")}">{b}</script>')
js=''.join(open(f'{HERE}/src/{f}').read()+'\n' for f in ['core.js','scenery.js','drive.js','ui.js'])
out=tpl.replace('/*CSS*/',open(HERE+'/src/style.css').read())
out=out.replace('/*THREE*/',open(HERE+'/package/build/three.min.js').read().replace('</script>','<\\/script>'))
out=out.replace('/*MANIFEST*/',json.dumps(man))
out=out.replace('/*DATA*/','\n'.join(parts))
out=out.replace('/*APP*/',js)
dst=sys.argv[1] if len(sys.argv)>1 else os.path.join(os.path.dirname(HERE),'i26-connector-3d.html')
open(dst,'w').write(out)
print(dst,round(len(out)/1e6,2),'MB')
