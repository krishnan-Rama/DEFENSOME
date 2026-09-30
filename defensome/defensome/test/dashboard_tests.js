
/* ---- headless smoke test ---- */
const _tabs=['overview','heatmap','family','species','ordination','correlation',
             'domains','clans','amethods','trees','compare','stats','data','map','methods'];
let _bad=0;
for(const t of _tabs){
  try{ const host=document.createElement('div'); TABS.find(x=>x[0]===t)[2](host);
       console.log('  OK   '+t); }
  catch(e){ _bad++; console.log('  FAIL '+t+' -> '+e.message);
            console.log('        '+String(e.stack||'').split('\n')[1].trim()); }
}
/* exercise interactive callbacks that only fire on user input */
try{ const h=document.createElement('div'); tabFamily(h);
     FAM.slice(0,8).forEach(f=>{ famSel.value=f; famSel.onchange(); });
     console.log('  OK   family switching across 8 families'); }
catch(e){ _bad++; console.log('  FAIL family switching -> '+e.message); }
try{ SP.slice(0,5).forEach(s=>{ const h=document.createElement('div'); renderSpecies(h,s,'per10k'); });
     console.log('  OK   species rendering'); }
catch(e){ _bad++; console.log('  FAIL species rendering -> '+e.message); }
try{ ['radial','rect'].forEach(()=>{ const h=document.createElement('div'); tabTrees(h); });
     console.log('  OK   trees tab builds'); }
catch(e){ _bad++; console.log('  FAIL trees -> '+e.message+' '+String(e.stack).split('\n')[1]); }
try{ const h=document.createElement('div'); tabCompare(h); console.log('  OK   compare'); }
catch(e){ _bad++; console.log('  FAIL compare -> '+e.message); }
try{ const cb=colourbar(rdbu,-2.5,2.5,'z',240,12);
     if(!cb.children.length) throw new Error('empty colourbar');
     console.log('  OK   colourbar renders '+cb.children.length+' elements'); }
catch(e){ _bad++; console.log('  FAIL colourbar -> '+e.message); }
try{ const sv=svgEl('svg',{width:100,height:50}); const f=figure(sv,'x');
     const bar=f.children[1];
     if(!bar||bar.children.length<2) throw new Error('no export buttons');
     console.log('  OK   figure export toolbar ('+bar.children.length+' buttons)'); }
catch(e){ _bad++; console.log('  FAIL figure toolbar -> '+e.message); }
console.log(_bad? '\n'+_bad+' FAILURES' : '\nall tab builders executed cleanly');
if(_bad) process.exitCode=1;

/* verify the maths, not just that it runs */
function approx(a,b,t=1e-6){ return Math.abs(a-b)<t; }
let bad=0; const chk=(n,c)=>{ console.log((c?'  OK   ':'  FAIL ')+n); if(!c) bad++; };

const v=[2,4,4,4,5,5,7,9], s=stats(v);
chk('stats mean', approx(s.mean,5));
chk('stats median', approx(s.med,4.5));
chk('stats sd (sample)', approx(s.sd, Math.sqrt(32/7)));

const z=zcol('raw','CYP'), zs=stats(z);
chk('z-score has mean 0', Math.abs(zs.mean)<1e-9);
chk('z-score has sd 1', approx(zs.sd,1,1e-9));

/* PCA against a case with a known answer: 2 correlated + 1 independent axis */
const n=3, C=[[1,0.9,0],[0.9,1,0],[0,0,1]];
const e=jacobiEig(C,n);
chk('eigenvalues sum to trace', approx(e.values.reduce((a,b)=>a+b,0),3,1e-8));
chk('largest eigenvalue is 1.9', approx(e.values[0],1.9,1e-8));
chk('eigenvalues sorted descending', e.values[0]>=e.values[1] && e.values[1]>=e.values[2]);
chk('leading vector loads on the correlated pair',
    Math.abs(e.vectors[0][0])>0.6 && Math.abs(e.vectors[0][2])<1e-6);

/* Newick round trip */
const t=parseNewick('((A:0.1,B:0.2)0.9:0.3,(C:0.4,D:0.5):0.6);');
chk('newick tip count', leaves(t).length===4);
chk('newick tip names', leaves(t).map(x=>x.name).sort().join('')==='ABCD');
chk('newick branch length', approx(leaves(t).find(x=>x.name==='D').len,0.5));
chk('internal support label not read as a tip name',
    walk(t).filter(x=>x.kids.length).every(x=>x.name===''));
const pr=pruneTree(parseNewick('((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);'), new Set(['A','C','D']));
chk('prune drops B and collapses the single-child node', leaves(pr).length===3);

/* invariant detection must match the raw counts */
const zeroFams=FAM.filter(f=>D.counts[f].reduce((a,b)=>a+b,0)===0);
const constFams=FAM.filter(f=>{const c=D.counts[f];return c.reduce((a,b)=>a+b,0)>0&&c.every(x=>x===c[0]);});
chk('invariant set = zero + constant, judged on RAW counts',
    Object.keys(INVARIANT).sort().join()===[...zeroFams,...constFams].sort().join());
console.log('  ->   invariant: '+JSON.stringify(INVARIANT));
console.log(bad? '\n'+bad+' NUMERIC FAILURES' : '\nnumerics correct');
if(bad) process.exitCode=1;
