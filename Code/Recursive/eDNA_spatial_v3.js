/* Browser controller for eDNA survey planner v3. */
(() => {
  'use strict';
  const E=createSpatialEngineV3(), $=id=>document.getElementById(id);
  // Extent from the R prediction data, clipped to Natural Earth ocean.
  const pacificStudyDomain = PACIFIC_WATER_DOMAIN_V3.geometry;
  const fmt=(n,d=2)=>Number.isFinite(n)?n.toLocaleString('en-US',{maximumFractionDigits:d}):'—';
  const state={domain:null,map:null,outline:null,draftLayer:null,draft:null,worker:null,cache:null,result:null,sampleSeed:911,timer:null};
  const status=(id,text,warn=false)=>{const node=$(id);node.textContent=text;node.classList.toggle('warn',warn);};
  function params() {
    return {mu:Number($('gpMu').value),alpha:Number($('gpAlpha').value),sigma:Number($('gpSigma').value),rho:Number($('rho').value),seed:Number($('gpSeed').value),convention:'r',
      epsilon:Number($('epsilon').value),omega:Number($('omega').value),beta:Number($('beta').value),theta:Number($('theta').value),gamma:Number($('gamma').value)};
  }
  function fieldKey(p) { return JSON.stringify([state.domain?.geometry||state.domain?.geo,p.mu,p.alpha,p.sigma,p.rho,p.seed,p.convention]); }
  function clearResults() {
    state.result=null; $('spatialResults').hidden=true;
    $('measuredEpsilon').textContent='—'; $('downloadSpatial').disabled=true;
  }
  function stop() {
    clearTimeout(state.timer);
    if(state.worker) {state.worker.terminate();state.worker=null;}
    $('cancelSpatial').disabled=true; $('runSpatial').disabled=false;
  }
  function planned() {
    $('targetEpsilon').textContent=fmt(Number($('epsilon').value),3);
    $('fieldCount').textContent=state.domain?fmt(state.domain.n,0):'—';
    $('evaluationCount').textContent=state.domain?`≈${fmt(state.domain.area/25,0)}`:'—';
    $('surveyCount').textContent='—';
    if(state.domain) {try {$('surveyCount').textContent=fmt(E.recommend(params(),state.domain.area).count,0);}catch{}}
  }
  function invalidate(auto=false) {
    stop();clearResults();planned();
    const reusable=state.cache && state.cache.key===fieldKey(params());
    $('newSamples').disabled=!reusable;
    status('simStatus',reusable?'Inputs changed. Reconstructing the same reference field…':'Ready to simulate a field from the current domain and GP parameters.');
    if(auto && reusable) state.timer=setTimeout(run,400);
  }
  function setDomain(vertices) {
    const d=Array.isArray(vertices)?E.domain(vertices):E.geometryDomain(vertices);
    stop();state.domain=d;state.cache=null;state.sampleSeed=911;clearResults();
    $('area').value=d.area.toFixed(2);
    $('area').dataset.exactArea=String(d.area);
    $('domainCoordinates').value=d.geometry?JSON.stringify(d.geometry):d.geo.map(v=>v.map(n=>Number(n.toFixed(6))).join(',')).join('\n');
    if(state.map) {
      if(state.outline)state.outline.remove();
      // Densify straight projected edges when displaying them geographically.
      const boundary=(d.polygons||[[d.polygon]]).map(part=>part.map(ring=>{
        const points=[];
        ring.forEach((a,i)=>{const b=ring[(i+1)%ring.length],steps=Math.max(1,Math.ceil(Math.hypot(b[0]-a[0],b[1]-a[1])/20));for(let j=0;j<steps;j++)points.push(E.inverse(a[0]+(b[0]-a[0])*j/steps,a[1]+(b[1]-a[1])*j/steps,d.center));});
        return points;
      }));
      state.outline=L.polygon(boundary,{color:'#006d77',weight:2,fillOpacity:.17,fillRule:'evenodd',interactive:false}).addTo(state.map);
      state.map.fitBounds(state.outline.getBounds(),{padding:[35,35],maxZoom:9});
    }
    status('domainStatus',`${fmt(d.area)} km²${vertices===pacificStudyDomain?' of water · Pacific example, land excluded':''} · ${fmt(d.n,0)} field points = ceil(area × 40 / 10,000) · one depth layer.${d.n>E.MAX_POINTS?' This exceeds the 2,000-point limit; reduce the domain before simulating.':''}`,d.n>E.MAX_POINTS);
    $('area').dispatchEvent(new Event('input',{bubbles:true}));
    planned();$('newSamples').disabled=true;
  }
  function updateDraft() {
    if(state.draftLayer)state.draftLayer.remove();
    if(state.map && state.draft?.length)state.draftLayer=L.polyline(state.draft,{color:'#d45b24',weight:3,dashArray:'5 6',interactive:false}).addTo(state.map);
    $('undoVertex').disabled=!state.draft?.length;
    $('finishDomain').disabled=!state.draft || state.draft.length<3;
    $('cancelDraw').disabled=!state.draft;
    $('drawDomain').disabled=!!state.draft || !state.map;
    $('domainMap').classList.toggle('drawing',!!state.draft);
    if(state.draft)status('domainStatus',`${state.draft.length} vertices · click around the boundary, then Finish polygon. Drag the map to pan.`);
  }
  function endDraft() {state.draft=null;updateDraft();}
  $('drawDomain').addEventListener('click',()=>{state.draft=[];updateDraft();});
  $('undoVertex').addEventListener('click',()=>{state.draft.pop();updateDraft();});
  $('cancelDraw').addEventListener('click',()=>{endDraft();status('domainStatus',state.domain?`Drawing cancelled. Current domain: ${fmt(state.domain.area)} km².`:'Drawing cancelled.');});
  $('finishDomain').addEventListener('click',()=>{try {setDomain(state.draft);endDraft();}catch(e){status('domainStatus',e.message,true);}});
  $('applyCoordinates').addEventListener('click',()=>{try {
    const text=$('domainCoordinates').value.trim();
    const vertices=text.startsWith('{')?JSON.parse(text):text.split(/\n+/).map(line=>line.trim().split(/[\s,]+/).map(Number));
    setDomain(vertices);endDraft();
  }catch(e){status('domainStatus',e.message+' The current domain has been retained.',true);}});
  $('exampleDomain').addEventListener('click',()=>{setDomain(pacificStudyDomain);endDraft();});
  $('worldView').addEventListener('click',()=>state.map?.setView([20,0],2));
  if(typeof L!=='undefined') {
    state.map=L.map('domainMap',{scrollWheelZoom:false,worldCopyJump:true,doubleClickZoom:false}).setView([20,0],2);
    L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png',{maxZoom:18,attribution:'&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'}).addTo(state.map);
    L.control.scale({imperial:false}).addTo(state.map);
    state.map.on('click',event=>{if(state.draft){state.draft.push([event.latlng.lat,event.latlng.lng]);updateDraft();}});
  } else {
    $('domainMap').textContent='The online map could not load. Use the coordinate boundary editor below; the GP simulation still works.';
    $('drawDomain').disabled=true;$('worldView').disabled=true;
  }
  // A fresh cancellable worker per job. Previously simulated points are reused
  // for epsilon/sample changes, never regenerated to improve the measured error.
  function run() {
    stop();clearResults();planned();
    try {
      if(!state.domain)throw Error('Draw a domain first.');
      const p=params();
      for(const id of ['gpMu','gpAlpha','gpSigma','gpSeed','rho','omega','beta','theta','gamma'])if($(id).value.trim()==='')throw Error('Complete all GP and equation parameters.');
      E.validate(p);
      const plan=E.recommend(p,state.domain.area), key=fieldKey(p), cached=state.cache?.key===key?state.cache.field:null;
      const cachedDense=cached?state.cache.dense:null;
      if(state.domain.n>E.MAX_POINTS)throw Error(`This domain requires ${fmt(state.domain.n,0)} points. Draw a domain ≤ 500,000 km² for the exact GP; density stays fixed at 40 / 10,000 km².`);
      $('runSpatial').disabled=true;$('newSamples').disabled=true;$('cancelSpatial').disabled=false;
      status('simStatus',cached?'Inferring the GP surface on the dense comparison grid…':`Simulating ${fmt(state.domain.n,0)} reference observations with the exact GP…`);
      const source=`const engine=(${createSpatialEngineV3.toString()})(); onmessage=event=>{try {
        const {domain,p,count,sampleSeed,cached,cachedDense}=event.data;
        const progress=stage=>fraction=>postMessage({progress:stage,fraction});
        const field=cached||engine.generate(domain,p,progress('Simulating reference field'));
        const dense=cachedDense||engine.referenceSurface(field,domain,p,progress('Evaluating full reference surface'));
        const inferred=count<=field.points.length?engine.infer(field,p,count,sampleSeed,progress('Inferring survey map')):null;
        const comparison=inferred?engine.compareSurfaces(dense,inferred,p,progress('Comparing surfaces point by point')):null;
        postMessage({field,inferred,dense,comparison});
      }catch(e){postMessage({error:e.message});}};`;
      const url=URL.createObjectURL(new Blob([source],{type:'text/javascript'}));
      let worker;
      try {worker=new Worker(url);}finally{URL.revokeObjectURL(url);}
      state.worker=worker;
      worker.onmessage=event=>{
        if(state.worker!==worker)return;
        const data=event.data;
        if(data.progress){status('simStatus',`${data.progress} · ${Math.round(data.fraction*100)}%`);return;}
        stop();
        if(data.error){status('simStatus',data.error,true);return;}
        state.cache={key,field:data.field,dense:data.dense};
        state.result={...data,domain:state.domain,p,plan,sampleSeed:state.sampleSeed};
        $('newSamples').disabled=false;render(state.result);
      };
      worker.onerror=event=>{if(state.worker!==worker)return;stop();status('simStatus',`The calculation could not finish: ${event.message || 'worker unavailable'}.`,true);};
      worker.postMessage({domain:state.domain,p,count:plan.count,sampleSeed:state.sampleSeed,cached,cachedDense});
    }catch(e){stop();status('simStatus',e.message,true);}
  }
  $('runSpatial').addEventListener('click',()=>{
    const entropy=crypto.getRandomValues(new Uint32Array(1))[0];
    $('gpSeed').value=E.nextSeed(Number($('gpSeed').value),entropy);
    state.sampleSeed=911;
    run();
  });
  $('cancelSpatial').addEventListener('click',()=>{stop();clearResults();status('simStatus','Calculation cancelled. Simulate again when ready.');});
  $('newSamples').addEventListener('click',()=>{state.sampleSeed++;run();});
  ['gpMu','gpAlpha','gpSigma','gpSeed'].forEach(id=>$(id).addEventListener('input',()=>invalidate(false)));
  $('gpMu').addEventListener('blur',()=>{if($('gpMu').value!==''&&Number.isFinite(Number($('gpMu').value)))$('gpMu').value=Number($('gpMu').value).toFixed(2);});
  window.addEventListener('plannerchange',()=>invalidate(true));
  const NS='http://www.w3.org/2000/svg';
  function add(svg,tag,attrs={},text) {
    const n=document.createElementNS(NS,tag);Object.entries(attrs).forEach(([k,v])=>n.setAttribute(k,v));
    if(text!==undefined)n.textContent=text;svg.appendChild(n);return n;
  }
  const sequential=[[25,52,102],[25,117,137],[93,181,156],[239,221,108]], errorColors=[[255,244,210],[249,164,91],[202,66,55],[107,22,59]];
  function color(t,palette) {
    t=Math.max(0,Math.min(1,t))*(palette.length-1);const i=Math.min(palette.length-2,Math.floor(t)),f=t-i;
    return `rgb(${palette[i].map((a,j)=>Math.round(a+(palette[i+1][j]-a)*f)).join(',')})`;
  }
  function legend(id,min,max,palette,label) {
    const n=$(id);n.replaceChildren();const title=document.createElement('div');title.textContent=label;n.appendChild(title);
    const bar=document.createElement('div');bar.className='legend-bar';bar.style.background=`linear-gradient(90deg,${palette.map(c=>`rgb(${c.join(',')})`).join(',')})`;n.appendChild(bar);
    const labels=document.createElement('div');labels.className='legend-labels';[min,(min+max)/2,max].forEach(v=>{const s=document.createElement('span');s.textContent=fmt(v,2);labels.appendChild(s);});n.appendChild(labels);
  }
  function concentration(logValue) {
    const value=Math.exp(logValue);
    const display=!Number.isFinite(value)||value===0?`exp(${fmt(logValue,2)})`:value<.01||value>=1e6?value.toExponential(3):fmt(value,2);
    return `${display} copies/L (${fmt(logValue,2)} ln(copies/L))`;
  }
  function showPoint(title,xy,result,reference,prediction,observed) {
    const [lat,lon]=E.inverse(...xy,result.domain.center);
    $('pointDetailsTitle').textContent=title;
    const lines=[`${lat.toFixed(5)}°, ${(((lon+180)%360+360)%360-180).toFixed(5)}°`];
    if(observed!==undefined)lines.push(`Observed concentration: ${concentration(observed)}`);
    lines.push(`Reference surface: ${concentration(reference)}`);
    if(prediction!==undefined)lines.push(`Inferred surface: ${concentration(prediction)}`,`Pointwise prediction error ε: ${fmt(Math.sqrt((prediction-reference)**2),4)} (log scale)`);
    $('pointDetailsText').textContent=lines.join('\n\n');
    $('pointDetails').showModal();
  }
  $('closePointDetails').addEventListener('click',()=>$('pointDetails').close());
  $('pointDetails').addEventListener('click',event=>{if(event.target===$('pointDetails'))$('pointDetails').close();});
  function drawMap(id,result,values,min,max,palette,overlay='none') {
    const svg=$(id);svg.replaceChildren();svg.setAttribute('role','group');
    const {domain:d,field,inferred,dense}=result, [xmin,xmax,ymin,ymax]=d.bounds;
    const scale=Math.min(425/(xmax-xmin),295/(ymax-ymin)), left=65+(425-(xmax-xmin)*scale)/2, bottom=332-(295-(ymax-ymin)*scale)/2;
    const x=v=>left+(v-xmin)*scale,y=v=>bottom-(v-ymin)*scale;
    const path=(d.polygons||[[d.polygon]]).flatMap(part=>part.map(ring=>ring.map((p,i)=>`${i?'L':'M'}${x(p[0])},${y(p[1])}`).join(' ')+' Z')).join(' ');
    add(svg,'path',{d:path,fill:'#f1f5f4','fill-rule':'evenodd',stroke:'#92aaa3','stroke-width':1.2});
    const clip=add(add(svg,'defs'),'clipPath',{id:id+'Clip'});add(clip,'path',{d:path,'clip-rule':'evenodd'});
    const dots=add(svg,'g',{'clip-path':`url(#${id}Clip)`});
    // Rasterize GP values (not sample-point interpolation) into one SVG image.
    // Each colored cell is evaluated at its own shared comparison coordinate.
    const canvas=document.createElement('canvas');canvas.width=1040;canvas.height=780;
    const ctx=canvas.getContext('2d');ctx.scale(2,2);
    dense.cells.forEach((cell,i)=>{
      ctx.fillStyle=color((values[i]-min)/(max-min||1),palette);
      ctx.fillRect(x(cell[0]),y(cell[1]+cell[3]),cell[2]*scale+.35,cell[3]*scale+.35);
    });
    add(dots,'image',{x:0,y:0,width:520,height:390,href:canvas.toDataURL('image/png'),'pointer-events':'none'});
    const indices=overlay==='all'?field.points.map((_,i)=>i):overlay==='survey'?inferred.selected:[];
    const radius=Math.min(4,Math.max(1.8,Math.sqrt((xmax-xmin)*(ymax-ymin)*scale*scale/field.points.length)*.16));
    indices.forEach(i=>{
      const p=field.points[i];
      const c=add(dots,'circle',{cx:x(p[0]),cy:y(p[1]),r:radius,fill:'transparent',stroke:'#fff','stroke-width':.9,tabindex:0,role:'button','aria-label':`Observation ${i+1}: ${concentration(field.truth[i])}. Show concentration details.`,'pointer-events':'all',style:'cursor:pointer'});
      const [lat,lon]=E.inverse(...p,d.center);
      add(c,'title',{},`Observation ${i+1} · ${lat.toFixed(4)}°, ${(((lon+180)%360+360)%360-180).toFixed(4)}°\nObserved: ${field.truth[i].toFixed(4)} ln(copies/L)${overlay==='survey'?'\nSelected survey sample':''}`);
      const inspect=()=>showPoint(`Observation ${i+1}${overlay==='survey'?' · survey sample':''}`,p,result,E.predictAt([p],field.points,field.referenceWeights,result.p)[0],inferred?.prediction[i],field.truth[i]);
      c.addEventListener('click',event=>{event.stopPropagation();inspect();});
      c.addEventListener('keydown',event=>{if(event.key==='Enter'||event.key===' '){event.preventDefault();event.stopPropagation();inspect();}});
    });
    // Surface clicks inspect the nearest evaluated cell; no invented values.
    svg.style.cursor='crosshair';
    svg.onclick=event=>{
      const matrix=svg.getScreenCTM();if(!matrix)return;
      const local=new DOMPoint(event.clientX,event.clientY).matrixTransform(matrix.inverse());
      const xy=[xmin+(local.x-left)/scale,ymin+(bottom-local.y)/scale];
      if(!E.contains(d,xy))return;
      let nearest=0,best=Infinity;
      dense.points.forEach((p,i)=>{const distance=(p[0]-xy[0])**2+(p[1]-xy[1])**2;if(distance<best){best=distance;nearest=i;}});
      showPoint(`Surface evaluation point ${nearest+1} · nearest grid location`,dense.points[nearest],result,dense.reference[nearest],result.comparison?.prediction[nearest]);
    };
    add(svg,'path',{d:path,fill:'none',stroke:'#6c8580','stroke-width':.8,'pointer-events':'none'});
    for(let i=0;i<=4;i++){
      const vx=xmin+(xmax-xmin)*i/4,vy=ymin+(ymax-ymin)*i/4;
      add(svg,'text',{x:x(vx),y:354,'text-anchor':'middle'},fmt(vx,0));
      add(svg,'text',{x:55,y:y(vy)+4,'text-anchor':'end'},fmt(vy,0));
    }
    add(svg,'text',{x:278,y:379,'text-anchor':'middle'},'Local easting (km)');
    add(svg,'text',{x:16,y:184,transform:'rotate(-90 16 184)','text-anchor':'middle'},'Local northing (km)');
    add(svg,'text',{x:482,y:24,'text-anchor':'end'},'N ↑');
  }
  function histogram(result) {
    const svg=$('errorHistogram');svg.replaceChildren();const r=result.comparison,eps=result.p.epsilon;
    const errors=r.rootSquared, dataMax=errors.reduce((m,v)=>Math.max(m,v),1e-9), bins=30, counts=Array(bins).fill(0);
    errors.forEach(e=>counts[Math.min(bins-1,Math.floor(e/dataMax*bins))]++);
    const xmax=Math.max(dataMax,eps*1.05,r.measuredEpsilon*1.05), ymax=Math.max(...counts,1);
    const x=v=>78+v/xmax*414,y=v=>326-v/ymax*268;
    for(let i=0;i<=4;i++) {const n=Math.round(ymax*i/4);add(svg,'line',{x1:78,x2:492,y1:y(n),y2:y(n),stroke:'#dce6e3'});add(svg,'text',{x:69,y:y(n)+4,'text-anchor':'end'},fmt(n,0));}
    counts.forEach((n,i)=>{const lo=i*dataMax/bins,hi=(i+1)*dataMax/bins;const b=add(svg,'rect',{x:x(lo)+.4,y:y(n),width:Math.max(.1,x(hi)-x(lo)-.8),height:326-y(n),fill:'#d98050'});add(b,'title',{},`${lo.toFixed(3)}–${hi.toFixed(3)}: ${n} points`);});
    [[eps,'#006d77','Target ε'],[r.measuredEpsilon,'#8e2437','Measured ε']].forEach(([v,c,label],i)=>{add(svg,'line',{x1:x(v),x2:x(v),y1:53,y2:326,stroke:c,'stroke-width':2,'stroke-dasharray':i?'':'5 4'});add(svg,'text',{x:60+i*210,y:25,style:`fill:${c}`},`${label}: ${fmt(v,3)}`);});
    for(let i=0;i<=5;i++)add(svg,'text',{x:x(xmax*i/5),y:348,'text-anchor':'middle'},fmt(xmax*i/5,2));
    add(svg,'text',{x:276,y:377,'text-anchor':'middle'},'Prediction error ε = √((inferred − reference)²)');
    add(svg,'text',{x:15,y:195,transform:'rotate(-90 15 195)','text-anchor':'middle'},'Frequency (dense evaluation points)');
    $('histogramNote').textContent=`All ${fmt(errors.length,0)} pointwise ε values. The measured prediction error ε is √mean(ε²). All values use log concentrations.`;
  }
  function render(result) {
    const {field,comparison:r,dense,plan,p}=result;
    $('spatialResults').hidden=false;
    const extent=values=>values.reduce(([lo,hi],value)=>[Math.min(lo,value),Math.max(hi,value)],[Infinity,-Infinity]);
    const [referenceMin,referenceMax]=extent(dense.reference);
    drawReferenceMap(result);
    if(!r) {
      ['predictionMap','errorMap','errorHistogram'].forEach(id=>{$(id).onclick=null;$(id).replaceChildren();add($(id),'text',{x:260,y:180,'text-anchor':'middle'},'Recommended survey exceeds reference points.');});
      ['predictionLegend','errorLegend','histogramNote','evaluationNote'].forEach(id=>$(id).textContent='');
      status('simStatus',`Simulated ${fmt(field.points.length,0)} field points. The equation requires ${fmt(plan.count,0)} survey samples, more than the available points. Increase ε to reconstruct without replacement. No smaller survey was substituted.`,true);
      return;
    }
    const [predictionMin,predictionMax]=extent(r.prediction);
    drawInferenceMap(result);
    const errors=r.rootSquared,errorMax=errors.reduce((m,v)=>Math.max(m,v),0);
    drawMap('errorMap',result,errors,0,errorMax,errorColors);
    legend('errorLegend',0,errorMax,errorColors,'ε = √((inferred − reference)²) · log scale');
    histogram(result);$('measuredEpsilon').textContent=fmt(r.measuredEpsilon,3);$('downloadSpatial').disabled=false;
    const extrapolation=plan.density<1.5||plan.density>20||p.rho<25||p.rho>1000;
    status('simStatus',`${fmt(plan.count,0)} survey samples / ${fmt(field.points.length,0)} simulated observations · ${fmt(dense.points.length,0)} points on the 5 × 5 km grid · measured ε ${fmt(r.measuredEpsilon,3)} ${r.measuredEpsilon<=p.epsilon?'meets':'exceeds'} target ε = ${fmt(p.epsilon,3)}.${extrapolation?' The sample-count equation is being extrapolated beyond its calibration range.':''}`,r.measuredEpsilon>p.epsilon||extrapolation);
    status('evaluationNote',`Compared the two GP surfaces at all ${fmt(dense.points.length,0)} shared grid centers, spaced 5 km apart in both directions (25 km² cells; approximately 10× the simulated-observation density). Every √((inferred − reference)²) value is included in the ε surface and histogram. ε = √(Σ difference² / ${fmt(dense.points.length,0)}), using natural-log concentrations. The reference surface conditions on all simulated observations; it is not the unobserved latent truth. Cells are clipped to the drawn domain, with centers outside it excluded from ε.`,false);
  }
  function drawReferenceMap(result) {
    const values=result.dense.reference,[min,max]=values.reduce(([lo,hi],v)=>[Math.min(lo,v),Math.max(hi,v)],[Infinity,-Infinity]);
    drawMap('truthMap',result,values,min,max,sequential,$('showReferencePoints').checked?'all':'none');
    legend('truthLegend',min,max,sequential,'ln(copies/L) · reference map min–max');
  }
  function drawInferenceMap(result) {
    const values=result.comparison.prediction,[min,max]=values.reduce(([lo,hi],v)=>[Math.min(lo,v),Math.max(hi,v)],[Infinity,-Infinity]);
    drawMap('predictionMap',result,values,min,max,sequential,$('showInferencePoints').checked?'survey':'none');
    legend('predictionLegend',min,max,sequential,'ln(copies/L) · inferred map min–max');
  }
  $('showReferencePoints').addEventListener('change',()=>{if(state.result)drawReferenceMap(state.result);});
  $('showInferencePoints').addEventListener('change',()=>{if(state.result?.comparison)drawInferenceMap(state.result);});
  $('downloadSpatial').addEventListener('click',()=>{
    const r=state.result;if(!r?.comparison)return;
    const header='evaluation_id,latitude,longitude,x_km,y_km,reference_surface_log,predicted_surface_log,residual_log,root_squared_difference_log,squared_error_log,area_km2,rho_km,mu,alpha,sigma,convention,field_seed,sample_seed,target_epsilon,required_density,required_samples,measured_epsilon,omega,beta,theta,gamma,original_point_count,grid_spacing_km,row_type,observation_id,observed_log,used_for_inference,used_to_estimate_epsilon';
    const selected=new Set(r.inferred.selected);
    const observationReference=E.predictAt(r.field.points,r.field.points,r.field.referenceWeights,r.p);
    const row=(xy,i,isObservation)=>{
      const ll=E.inverse(...xy,r.domain.center);ll[1]=((ll[1]+180)%360+360)%360-180;
      const reference=isObservation?observationReference[i]:r.dense.reference[i];
      const prediction=isObservation?r.inferred.prediction[i]:r.comparison.prediction[i];
      const e=prediction-reference,p=r.p;
      return [isObservation?'':i+1,...ll,...xy,reference,prediction,e,Math.sqrt(e*e),e*e,r.domain.area,p.rho,p.mu,p.alpha,p.sigma,p.convention,p.seed,r.sampleSeed,p.epsilon,r.plan.density,r.plan.count,r.comparison.measuredEpsilon,p.omega,p.beta,p.theta,p.gamma,r.field.points.length,r.dense.spacingKm,isObservation?'simulated_observation':'evaluation_point',isObservation?i+1:'',isObservation?r.field.truth[i]:'',isObservation&&selected.has(i),!isObservation].join(',');
    };
    // Survey samples are original observations, not dense evaluation locations.
    const rows=[...r.field.points.map((xy,i)=>row(xy,i,true)),...r.dense.points.map((xy,i)=>row(xy,i,false))];
    const url=URL.createObjectURL(new Blob([header+'\n'+rows.join('\n')+'\n'],{type:'text/csv;charset=utf-8'}));
    const a=document.createElement('a');a.href=url;a.download='eDNA_generated_data_v3.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),1000);
  });
  setDomain(pacificStudyDomain);
  run();
})();
