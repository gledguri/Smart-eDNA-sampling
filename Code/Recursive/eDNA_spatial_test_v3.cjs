// Run: node Code/eDNA_spatial_test_v3.cjs
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const E = require('./eDNA_spatial_engine_v3.js')();
const p = {mu:Math.log(100),alpha:4,sigma:3,rho:100,seed:317,convention:'r',epsilon:.8,omega:1.65,beta:-.176,theta:-.228,gamma:-.0684};
const near=(a,b,t=1e-8)=>assert.ok(Math.abs(a-b)<t,`${a} != ${b}`);
assert.equal(E.nextSeed(317,317),318);
assert.equal(E.nextSeed(317,999),999);
assert.equal(E.nextSeed(4294967295,4294967295),0);
// Check the actual inline planner still parses after copying v2.
const html=fs.readFileSync(__dirname+'/eDNA_survey_planner_v3.html','utf8');
assert.ok(!html.includes('id="gpConvention"'));
assert.ok(!html.includes('id="evaluationMultiplier"'));
assert.ok(!/RMSE/i.test(html));
const depthLayersTag=html.match(/<input id="depthLayers"[^>]+>/)[0];
assert.ok(!depthLayersTag.includes('readonly'));assert.ok(!depthLayersTag.includes('max="1"'));
// Execute the actual isolated H handler: no simulation/curve update is available
// in this sandbox, so an accidental dependency is also caught.
const totalSource=html.match(/    function updateTotal\(\) \{[\s\S]*?\n    \}/)[0];
let h=1,validity='';
const totalContext={read:()=>({...p,area:148950.61611496544,depthLayers:h}),requiredDensity:v=>E.recommend(v,v.area),finitePositive:(...v)=>v.every(x=>Number.isFinite(x)&&x>0),el:{depthLayers:{setCustomValidity:s=>validity=s}},out:{requiredTotal:{textContent:''}},fmt:n=>String(n)};
vm.createContext(totalContext);vm.runInContext(totalSource,totalContext);
totalContext.updateTotal();const oneDepth=Number(totalContext.out.requiredTotal.textContent);
h=5;totalContext.updateTotal();assert.equal(Number(totalContext.out.requiredTotal.textContent),oneDepth*5);
h=0;totalContext.updateTotal();assert.equal(totalContext.out.requiredTotal.textContent,'—');assert.ok(validity);
h=1.5;totalContext.updateTotal();assert.equal(totalContext.out.requiredTotal.textContent,'—');
h=2;totalContext.updateTotal();assert.equal(Number(totalContext.out.requiredTotal.textContent),oneDepth*2);assert.equal(validity,'');
assert.ok(html.includes("el.depthLayers.addEventListener('input', updateTotal)"));
for(const match of html.matchAll(/<script(?: [^>]*)?>([\s\S]*?)<\/script>/g))new vm.Script(match[1]);
assert.equal(new Set([...html.matchAll(/\bid="([^"]+)"/g)].map(m=>m[1])).size,[...html.matchAll(/\bid="([^"]+)"/g)].length);
for(const asset of ['eDNA_spatial_v3.js','eDNA_spatial_engine_v3.js','eDNA_spatial_v3.css','eDNA_pacific_water_v3.js']) {
  const tag=asset.endsWith('.css')?'style':'script';
  const embedded=html.split(`<${tag} data-source="${asset}">`)[1]?.split(`</${tag}>`)[0].trim();
  assert.equal(embedded,fs.readFileSync(__dirname+'/'+asset,'utf8').trim(),`${asset} embedded copy is stale`);
}
// Geometry, inverse projection, concavity, crossings, and dateline handling.
const pacific=require('./eDNA_pacific_water_v3.js');
const d=E.geometryDomain(pacific.geometry);
assert.equal(d.n,Math.ceil(d.area/250));
assert.ok(d.area>185000 && d.area<188000,'Pacific area must exclude land from the R extent');
assert.equal(d.polygons[0].length,2,'retain the mapped island exclusion');
for(const [lat,lon] of [[46.975,-123.815],[44.635,-124.05],[41.756,-124.20]])assert.ok(!E.contains(d,E.project(lat,lon,d.center)),'coastal towns must not be sampled');
assert.ok(E.contains(d,E.project(43,-125,d.center)),'retain offshore Pacific water');
assert.ok(E.locations(d,317).every(p=>E.contains(d,p)));
assert.ok(E.evaluationGrid(d).points.every(p=>E.contains(d,p)));
// An island hole removes both area and sampling locations; a separate water
// polygon remains available. This tests topology, independently of coastline data.
const outer=[[0,0],[.3,0],[.3,.3],[0,.3],[0,0]],hole=[[.1,.1],[.2,.1],[.2,.2],[.1,.2],[.1,.1]],other=[[.4,0],[.5,0],[.5,.1],[.4,.1],[.4,0]];
const islands=E.geometryDomain({type:'MultiPolygon',coordinates:[[outer,hole],[other]]});
const noHole=E.geometryDomain({type:'MultiPolygon',coordinates:[[outer],[other]]});
assert.ok(islands.area<noHole.area);
assert.ok(!E.contains(islands,E.project(.15,.15,islands.center)));
assert.ok(E.contains(islands,E.project(.05,.45,islands.center)));
assert.ok(E.locations(islands,42).every(p=>E.contains(islands,p)));
assert.ok(E.evaluationGrid(islands).points.every(p=>E.contains(islands,p)));
for(const ll of d.geo){const xy=E.project(...ll,d.center),back=E.inverse(...xy,d.center);near(back[0],ll[0]);near(back[1],ll[1]);}
const seam=E.domain([[0,179],[0,-179],[1,-179],[1,179]]);
const ordinary=E.domain([[0,-1],[0,1],[1,1],[1,-1]]);
near(seam.area,ordinary.area,1e-6);
assert.throws(()=>E.domain([[0,0],[1,1],[0,1],[1,0]]),/crosses/);
assert.throws(()=>E.domain([[0,0],[0,0],[1,1]]),/repeated/);
assert.throws(()=>E.domain([[81,0],[81,1],[82,1]]),/80/);
const concave=E.domain([[0,0],[0,2],[1,1],[2,2],[2,0]]);
assert.ok(E.locations(concave,77).every(point=>E.inside(point,concave.polygon)));
// One-observation closed-form posterior verifies amplitude AND noise semantics.
const one={points:[[0,0],[100,0]],truth:[p.mu+5,p.mu-2]};
const oneFit=E.infer(one,p,1,20),ix=oneFit.selected[0],c=p.alpha,noise=p.sigma*p.sigma,jitter=1e-10*(c+noise);
one.points.forEach((x,i)=>near(oneFit.prediction[i],p.mu+E.covariance(x,one.points[ix],p)/(c+noise+jitter)*(one.truth[ix]-p.mu)));
near(E.covariance([0,0],[0,0],p),4);
assert.throws(()=>E.validate({...p,convention:'stan'}),/convention/);
// User's 1,000-point / 50-sample example, exact 250,000 km² planar domain.
const square={area:250000,n:1000,polygon:[[-250,-250],[250,-250],[250,250],[-250,250]],bounds:[-250,250,-250,250],center:[0,0]};
const density=2,epsilon=Math.exp(p.omega+p.theta*Math.log(p.rho)+(p.beta+p.gamma*Math.log(p.rho))*Math.log(density));
const plan=E.recommend({...p,epsilon},square.area);assert.equal(plan.count,50);
const field=E.generate(square,p),fit=E.infer(field,p,plan.count,911);
assert.equal(field.points.length,1000);assert.equal(fit.prediction.length,1000);assert.equal(new Set(fit.selected).size,50);
const direct=Math.sqrt(fit.prediction.reduce((sum,value,i)=>sum+(value-field.truth[i])**2,0)/1000);
near(fit.measuredEpsilon,direct);assert.ok(fit.measuredEpsilon>fit.mae);assert.ok(field.points.every(point=>E.inside(point,square.polygon)));
assert.deepEqual(field,E.generate(square,p));
assert.deepEqual(fit,E.infer(field,p,50,911));
assert.notDeepEqual(fit.selected,E.infer(field,p,50,912).selected);
const smaller=E.infer(field,p,25,911);assert.deepEqual(smaller.selected,fit.selected.slice(0,25));
// All-site noiseless reconstruction should recover the original field.
const tiny={...square,n:25}, noiseless={...p,sigma:0,rho:30};
const f0=E.generate(tiny,noiseless),r0=E.infer(f0,noiseless,25,911);assert.ok(r0.measuredEpsilon<1e-6);
assert.throws(()=>E.infer(field,p,1001,911),/exceeds/);
assert.throws(()=>E.generate({...square,n:2001},p),/2000|2,000/);
assert.throws(()=>E.recommend({...p,beta:1,gamma:0},square.area),/decreasing/);
assert.throws(()=>E.generate(tiny,{...p,sigma:-1}),/nonnegative/);
// Empirical one-point covariance confirms the simulator's marginal variance.
let sum=0,sum2=0;const trials=3000;
for(let seed=0;seed<trials;seed++){const z=E.generate({...tiny,n:1},{...p,seed}).truth[0]-p.mu;sum+=z;sum2+=z*z;}
assert.ok(Math.abs(sum/trials)<.3);assert.ok(Math.abs(sum2/trials-13)<1.8);
// Dense surface comparison: complete 5 km × 5 km lattice, shared by all maps.
const dense=E.referenceSurface(field,square,p);
assert.equal(dense.points.length,10000);assert.equal(dense.cells.length,10000);
assert.equal(dense.spacingKm,5);assert.equal(dense.points.length,10*field.points.length);
assert.ok(dense.cells.every(([, ,w,h])=>w===5&&h===5),'each rendered surface cell should be exactly 25 km²');
assert.ok(dense.points.every(point=>E.inside(point,square.polygon)));
assert.equal(dense.nx,100);assert.equal(dense.ny,100);
const xSpacing=5,ySpacing=5;
assert.ok(dense.points.every(([x,y])=>{
  const gx=(x-square.bounds[0])/xSpacing-.5,gy=(y-square.bounds[2])/ySpacing-.5;
  return Math.abs(gx-Math.round(gx))<1e-8&&Math.abs(gy-Math.round(gy))<1e-8;
}),'evaluation coordinates should align to a regular two-dimensional grid');
const rowCounts=new Map();dense.points.forEach(([,y])=>{const row=Math.round((y-square.bounds[2])/ySpacing-.5);rowCounts.set(row,(rowCounts.get(row)||0)+1);});
assert.equal(rowCounts.size,100);assert.ok([...rowCounts.values()].every(count=>count===100),'the complete square grid should contain every row without gaps');
assert.deepEqual(E.evaluationGrid(square).points,dense.points);
const denseConcave=E.evaluationGrid(concave);
assert.ok(denseConcave.points.length>0);
assert.ok(denseConcave.points.every(point=>E.inside(point,concave.polygon)));
const comparison=E.compareSurfaces(dense,fit,p);
assert.equal(comparison.rootSquared.length,10000);
let ss=0;
for(let i=0;i<10000;i++) {
  const difference=comparison.prediction[i]-dense.reference[i];
  near(comparison.rootSquared[i],Math.sqrt(difference*difference));
  near(comparison.rootSquared[i],Math.abs(difference));ss+=difference*difference;
}
near(comparison.measuredEpsilon,Math.sqrt(ss/10000));
// Execute the actual download handler and inspect the emitted CSV. The original
// observations and the evaluation lattice are distinct sets of locations.
assert.match(html,/>Download generated data<\/button>/);
const controller=fs.readFileSync(__dirname+'/eDNA_spatial_v3.js','utf8');
const downloadHandler=controller.match(/  \$\('downloadSpatial'\)\.addEventListener\('click',\(\)=>\{[\s\S]*?\n  \}\);/)[0];
function checkDownload(inferred,sampleSeed) {
  let clickHandler,csv,clicked=false;
  const result={field,inferred,dense,comparison:E.compareSurfaces(dense,inferred,p),domain:square,p,plan:{...plan,count:inferred.selected.length},sampleSeed};
  const anchor={click:()=>clicked=true};
  vm.runInNewContext(downloadHandler,{E,state:{result},$:()=>({addEventListener:(event,fn)=>{assert.equal(event,'click');clickHandler=fn;}}),Blob:class {constructor(parts){csv=parts.join('');}},URL:{createObjectURL:()=> 'blob:test',revokeObjectURL:()=>{}},document:{createElement:()=>anchor},setTimeout:fn=>fn()});
  clickHandler();assert.ok(clicked);assert.equal(anchor.download,'eDNA_generated_data_v3.csv');
  const [header,...lines]=csv.trim().split('\n'),keys=header.split(',');
  const rows=lines.map(line=>{const cells=line.split(',');assert.equal(cells.length,keys.length);return Object.fromEntries(keys.map((key,i)=>[key,cells[i]]));});
  assert.equal(rows.length,field.points.length+dense.points.length);
  const observations=rows.filter(r=>r.row_type==='simulated_observation');
  const evaluations=rows.filter(r=>r.used_to_estimate_epsilon==='true');
  assert.equal(observations.length,1000);assert.equal(evaluations.length,10000);
  const flagged=rows.filter(r=>r.used_for_inference==='true');
  assert.deepEqual(flagged.map(r=>Number(r.observation_id)-1).sort((a,b)=>a-b),[...inferred.selected].sort((a,b)=>a-b));
  observations.forEach((r,i)=>{assert.equal(Number(r.observation_id),i+1);near(Number(r.observed_log),field.truth[i]);near(Number(r.x_km),field.points[i][0]);near(Number(r.y_km),field.points[i][1]);near(Number(r.predicted_surface_log),inferred.prediction[i]);assert.equal(r.used_to_estimate_epsilon,'false');assert.equal(r.evaluation_id,'');});
  evaluations.forEach((r,i)=>{assert.equal(r.row_type,'evaluation_point');assert.equal(r.used_for_inference,'false');assert.equal(r.observed_log,'');assert.equal(Number(r.evaluation_id),i+1);near(Number(r.reference_surface_log),dense.reference[i]);near(Number(r.predicted_surface_log),result.comparison.prediction[i]);});
  near(Math.sqrt(evaluations.reduce((sum,r)=>sum+Number(r.squared_error_log),0)/evaluations.length),result.comparison.measuredEpsilon);
}
checkDownload(fit,911);
checkDownload(E.infer(field,p,25,912),912);
console.log('PASS: generated-data download contains original observations and evaluation points; inference flags match Figure 2 before and after reselection.');
// Independently evaluate the full-reference GP formula at a dense coordinate.
const allWeights=E.solve(E.cholesky(field.points,p),field.truth.map(z=>z-p.mu));
const probe=dense.points[123];
near(dense.reference[123],p.mu+field.points.reduce((s,x,i)=>s+E.covariance(probe,x,p)*allWeights[i],0));
// Using every observation makes both posterior surfaces equal, even with noise.
const smallField=E.generate(tiny,p),smallDense=E.referenceSurface(smallField,tiny,p);
const allFit=E.infer(smallField,p,tiny.n,911);
assert.ok(E.compareSurfaces(smallDense,allFit,p).measuredEpsilon<1e-8);
console.log(`PASS: dense surfaces, complete 5 km lattice at 10,000 locations, concave clipping, rooted-square differences, and full-survey equality.\n1,000 observations / 50 samples / 10,000 comparisons: surface ε ${comparison.measuredEpsilon.toFixed(6)}.`);
console.log(`PASS: geometry, GP formulas, reproducibility, sample subsets, boundaries, and ε evaluation.\n1,000 field points / 50 survey samples: target ${epsilon.toFixed(6)}, measured ε ${fit.measuredEpsilon.toFixed(6)}.\nPacific example: ${d.area.toFixed(2)} km², ${d.n} field points, ${E.recommend(p,d.area).count} survey samples.`);
