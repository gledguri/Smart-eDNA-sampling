/* v3 numerical engine. No DOM or dependencies; also runs inside a Blob worker.
   Exact squared-exponential GP, using the R simulator amplitude convention (α is variance). */
function createSpatialEngineV3() {
  'use strict';
  const R = 6371.0088, RAD = Math.PI / 180, MAX_POINTS = 2000;
  function nextSeed(previous,entropy) {const seed=entropy>>>0;return seed===previous?(seed+1)>>>0:seed;}
  function random(seed) {
    return () => {
      seed = (seed + 0x6D2B79F5) | 0;
      let t = Math.imul(seed ^ seed >>> 15, 1 | seed);
      t ^= t + Math.imul(t ^ t >>> 7, 61 | t);
      return ((t ^ t >>> 14) >>> 0) / 4294967296;
    };
  }
  function normal(rng) {
    return Math.sqrt(-2 * Math.log(Math.max(1e-15, rng()))) * Math.cos(2 * Math.PI * rng());
  }
  function project(lat, lon, center) {
    const p = lat * RAD, p0 = center[0] * RAD, dl = (lon - center[1]) * RAD;
    const den = 1 + Math.sin(p0) * Math.sin(p) + Math.cos(p0) * Math.cos(p) * Math.cos(dl);
    if (den < 1e-8) throw Error('This polygon is too wide for a regional projection.');
    const k = Math.sqrt(2 / den);
    return [R * k * Math.cos(p) * Math.sin(dl), R * k * (Math.cos(p0) * Math.sin(p) - Math.sin(p0) * Math.cos(p) * Math.cos(dl))];
  }
  function inverse(x, y, center) {
    const r = Math.hypot(x, y), c = 2 * Math.asin(Math.min(1, r / (2 * R))), p0 = center[0] * RAD;
    if (r < 1e-12) return center.slice();
    const lat = Math.asin(Math.cos(c) * Math.sin(p0) + y * Math.sin(c) * Math.cos(p0) / r) / RAD;
    const lon = center[1] + Math.atan2(x * Math.sin(c), r * Math.cos(p0) * Math.cos(c) - y * Math.sin(p0) * Math.sin(c)) / RAD;
    return [lat, lon];
  }
  function cross(a, b, c) { return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0]); }
  function onSegment(a,b,c) { return Math.abs(cross(a,b,c)) < 1e-8 && c[0] >= Math.min(a[0],b[0])-1e-8 && c[0] <= Math.max(a[0],b[0])+1e-8 && c[1] >= Math.min(a[1],b[1])-1e-8 && c[1] <= Math.max(a[1],b[1])+1e-8; }
  function intersects(a,b,c,d) {
    return (cross(a,b,c)*cross(a,b,d)<0 && cross(c,d,a)*cross(c,d,b)<0) || onSegment(a,b,c) || onSegment(a,b,d) || onSegment(c,d,a) || onSegment(c,d,b);
  }
  function domain(vertices) {
    if (!Array.isArray(vertices) || vertices.length < 3 || vertices.length > 200) throw Error('Use between 3 and 200 boundary vertices.');
    if (vertices.some(v => v.length !== 2 || !v.every(Number.isFinite) || Math.abs(v[0]) > 80)) throw Error('Use finite latitude/longitude pairs between 80°S and 80°N.');
    // Unwrap around the first longitude, so small dateline-crossing domains work.
    const geo = vertices.map(([lat,lon]) => [lat, vertices[0][1] + ((lon-vertices[0][1]+180)%360+360)%360-180]);
    const center = [geo.reduce((s,v)=>s+v[0],0)/geo.length, geo.reduce((s,v)=>s+v[1],0)/geo.length];
    const polygon = geo.map(v=>project(...v,center));
    for (let i=0;i<polygon.length;i++) {
      const a=polygon[i], b=polygon[(i+1)%polygon.length];
      if (Math.hypot(a[0]-b[0],a[1]-b[1])<1e-6) throw Error('Remove repeated or identical neighboring vertices.');
      for (let j=i+1;j<polygon.length;j++) {
        if (Math.hypot(a[0]-polygon[j][0],a[1]-polygon[j][1])>2000) throw Error('Use a regional domain spanning at most 2,000 km.');
        if (j===i+1 || (i===0 && j===polygon.length-1)) continue;
        if (intersects(a,b,polygon[j],polygon[(j+1)%polygon.length])) throw Error('The boundary crosses or touches itself. Draw a simple polygon.');
      }
    }
    const area = Math.abs(polygon.reduce((s,a,i)=>{const b=polygon[(i+1)%polygon.length]; return s+a[0]*b[1]-b[0]*a[1];},0))/2;
    if (area < 0.01) throw Error('The domain has negligible area. Draw a larger polygon.');
    const n = Math.max(1,Math.ceil(area*40/10000));
    const bounds = [Math.min(...polygon.map(p=>p[0])),Math.max(...polygon.map(p=>p[0])),Math.min(...polygon.map(p=>p[1])),Math.max(...polygon.map(p=>p[1]))];
    return {geo,center,polygon,area,n,bounds};
  }
  function inside(p, polygon) {
    let yes = false;
    for (let i=0,j=polygon.length-1;i<polygon.length;j=i++) {
      const a=polygon[i],b=polygon[j];
      if ((a[1]>p[1]) !== (b[1]>p[1]) && p[0]<(b[0]-a[0])*(p[1]-a[1])/(b[1]-a[1])+a[0]) yes=!yes;
    }
    return yes;
  }
  function geometryDomain(input) {
    const geometry=input?.type==='Feature'?input.geometry:input;
    if(!geometry||!['Polygon','MultiPolygon'].includes(geometry.type))throw Error('Use a GeoJSON Polygon or MultiPolygon.');
    const coordinates=geometry.type==='Polygon'?[geometry.coordinates]:geometry.coordinates;
    if(!Array.isArray(coordinates)||!coordinates.length)throw Error('The domain has no polygons.');
    let count=0;
    const rings=coordinates.map(part=>{
      if(!Array.isArray(part)||!part.length)throw Error('Each polygon needs an outer boundary.');
      return part.map(ring=>{
        if(!Array.isArray(ring)||ring.length<4)throw Error('Each GeoJSON ring needs at least four coordinates, including closure.');
        if(ring.some(p=>!Array.isArray(p)||p.length!==2||!p.every(Number.isFinite)||Math.abs(p[0])>180||Math.abs(p[1])>80))throw Error('Use valid GeoJSON longitude/latitude pairs within 80°S–80°N.');
        if(ring[0].some((v,i)=>v!==ring[ring.length-1][i]))throw Error('Close every GeoJSON ring by repeating its first coordinate.');
        count+=ring.length-1;
        return ring.slice(0,-1).map(([lon,lat])=>[lat,lon]);
      });
    });
    if(count>20000)throw Error('Use a boundary with at most 20,000 vertices.');
    const all=rings.flat(2),latitudes=all.map(p=>p[0]),longitudes=all.map(p=>p[1]);
    const center=[(Math.min(...latitudes)+Math.max(...latitudes))/2,(Math.min(...longitudes)+Math.max(...longitudes))/2];
    const polygons=rings.map(part=>part.map(ring=>ring.map(p=>project(...p,center))));
    const xy=polygons.flat(2),bounds=[Math.min(...xy.map(p=>p[0])),Math.max(...xy.map(p=>p[0])),Math.min(...xy.map(p=>p[1])),Math.max(...xy.map(p=>p[1]))];
    if(Math.hypot(bounds[1]-bounds[0],bounds[3]-bounds[2])>2000)throw Error('Use a regional domain spanning at most 2,000 km.');
    const ringArea=ring=>Math.abs(ring.reduce((s,a,i)=>{const b=ring[(i+1)%ring.length];return s+a[0]*b[1]-b[0]*a[1];},0))/2;
    let area=0;
    for(const part of polygons) {
      for(const ring of part)for(let i=0;i<ring.length;i++) {
        const a=ring[i],b=ring[(i+1)%ring.length];
        if(Math.hypot(a[0]-b[0],a[1]-b[1])<1e-8)throw Error('Remove repeated neighboring boundary vertices.');
        for(let j=i+2;j<ring.length;j++)if(!(i===0&&j===ring.length-1)&&intersects(a,b,ring[j],ring[(j+1)%ring.length]))throw Error('A boundary crosses or touches itself.');
      }
      for(let i=1;i<part.length;i++) {
        if(!part[i].every(p=>inside(p,part[0])))throw Error('Island exclusions must lie within their outer boundary.');
      }
      const net=ringArea(part[0])-part.slice(1).reduce((sum,hole)=>sum+ringArea(hole),0);
      if(net<=0)throw Error('Island exclusions exceed the polygon area.');
      area+=net;
    }
    if(area<0.01)throw Error('The domain has negligible area.');
    return {geometry,geo:rings[0][0],center,polygons,polygon:polygons[0][0],bounds,area,n:Math.max(1,Math.ceil(area/250))};
  }
  function contains(d,p) {
    if(!d.polygons)return inside(p,d.polygon);
    return d.polygons.some(part=>inside(p,part[0])&&!part.slice(1).some(hole=>inside(p,hole)));
  }
  function locations(d, seed) {
    if (d.n>MAX_POINTS) throw Error(`This domain needs ${d.n.toLocaleString()} field points at 40 / 10,000 km². The exact GP supports ${MAX_POINTS}; draw a domain ≤ 500,000 km².`);
    const rng=random(seed), points=[], [xmin,xmax,ymin,ymax]=d.bounds;
    for(let attempts=0;points.length<d.n && attempts<2000000;attempts++) {
      const p=[xmin+rng()*(xmax-xmin),ymin+rng()*(ymax-ymin)];
      if(contains(d,p)) points.push(p);
    }
    if(points.length!==d.n) throw Error('This very thin domain could not be sampled efficiently. Draw a more compact domain.');
    return points;
  }
  function validate(p) {
    if (![p.mu,p.alpha,p.sigma,p.rho].every(Number.isFinite) || p.alpha<=0 || p.rho<=0 || p.sigma<0 || !Number.isFinite(p.alpha*p.alpha+p.sigma*p.sigma)) throw Error('Enter a finite μ, positive α and ρ, and nonnegative σ.');
    if (!Number.isInteger(p.seed) || p.seed<0 || p.seed>4294967295) throw Error('Use a whole random seed from 0 to 4,294,967,295.');
    if (p.convention!=='r') throw Error('Unknown covariance convention.');
  }
  function covariance(a,b,p) {
    const variance=p.alpha;
    return variance*Math.exp(-((a[0]-b[0])**2+(a[1]-b[1])**2)/(2*p.rho*p.rho));
  }
  function cholesky(points,p,progress=()=>{}) {
    const n=points.length, l=Array.from({length:n},()=>new Float64Array(n));
    const variance=p.alpha;
    const jitter=1e-10*Math.max(variance+p.sigma*p.sigma,1);
    for(let i=0;i<n;i++) {
      for(let j=0;j<=i;j++) {
        let z=covariance(points[i],points[j],p)+(i===j?p.sigma*p.sigma+jitter:0);
        const li=l[i],lj=l[j];
        for(let k=0;k<j;k++) z-=li[k]*lj[k];
        if(i===j && (!(z>0)||!Number.isFinite(z))) throw Error('Covariance is numerically unstable. Reduce extreme parameter values.');
        li[j]=i===j?Math.sqrt(z):z/lj[j];
      }
      if(i%100===0) progress(i/n);
    }
    return l;
  }
  function solve(l,b) {
    const n=b.length, a=new Float64Array(n), z=new Float64Array(n);
    for(let i=0;i<n;i++) {let v=b[i];for(let j=0;j<i;j++)v-=l[i][j]*z[j];z[i]=v/l[i][i];}
    for(let i=n-1;i>=0;i--) {let v=z[i];for(let j=i+1;j<n;j++)v-=l[j][i]*a[j];a[i]=v/l[i][i];}
    return a;
  }
  function generate(d,p,progress=()=>{}) {
    validate(p);
    const points=locations(d,p.seed), rng=random((p.seed^0xa5a5a5a5)>>>0), l=cholesky(points,p,progress);
    const normals=points.map(()=>normal(rng));
    const truth=points.map((_,i)=>{let z=p.mu;for(let j=0;j<=i;j++)z+=l[i][j]*normals[j];return z;});
    if(!truth.every(Number.isFinite))throw Error('The simulated field overflowed. Reduce extreme parameter values.');
    const referenceWeights=solve(l,truth.map(z=>z-p.mu));
    return {points,truth,referenceWeights};
  }
  function recommend(p,area) {
    if (![p.epsilon,p.rho,area].every(x=>Number.isFinite(x)&&x>0) || ![p.omega,p.beta,p.theta,p.gamma].every(Number.isFinite)) throw Error('Enter valid positive ε, ρ and area, and finite equation coefficients.');
    const a=p.omega+p.theta*Math.log(p.rho), b=p.beta+p.gamma*Math.log(p.rho);
    if(b>=-1e-8)throw Error('The equation must have a decreasing error–effort curve.');
    const density=Math.exp((Math.log(p.epsilon)-a)/b), count=Math.max(1,Math.ceil(density*area/10000-1e-10));
    if (!(density>0) || !Number.isFinite(density) || !Number.isSafeInteger(count)) throw Error('The required effort is outside the numeric range. Adjust ε or the coefficients.');
    return {density,count};
  }
  function infer(field,p,count,sampleSeed,progress=()=>{}) {
    validate(p);
    if(!Number.isInteger(count)||count<1||count>field.points.length) throw Error('The recommended survey exceeds the available reference points. Increase ε to require no more than 40 samples / 10,000 km².');
    const rng=random(sampleSeed), order=field.points.map((_,i)=>i);
    for(let i=order.length-1;i>0;i--){const j=Math.floor(rng()*(i+1));[order[i],order[j]]=[order[j],order[i]];}
    const selected=order.slice(0,count), samples=selected.map(i=>field.points[i]);
    const l=cholesky(samples,p,progress), w=solve(l,selected.map(i=>field.truth[i]-p.mu));
    const prediction=field.points.map(x=>{let z=p.mu;for(let j=0;j<count;j++)z+=covariance(x,samples[j],p)*w[j];return z;});
    const residual=prediction.map((z,i)=>z-field.truth[i]);
    const measuredEpsilon=Math.sqrt(residual.reduce((s,r)=>s+r*r,0)/residual.length);
    if(!Number.isFinite(measuredEpsilon))throw Error('Prediction overflowed. Reduce extreme parameter values.');
    const mae=residual.reduce((s,r)=>s+Math.abs(r),0)/residual.length;
    const bias=residual.reduce((s,r)=>s+r,0)/residual.length;
    return {selected,prediction,residual,measuredEpsilon,mae,bias,samples,weights:w};
  }
  // One evaluation at every 5 × 5 km lattice center inside the study domain.
  // At 40 simulated sites / 10,000 km², this is approximately 10× as dense.
  function evaluationGrid(d) {
    const spacing=5,[xmin,xmax,ymin,ymax]=d.bounds;
    const width=xmax-xmin,height=ymax-ymin;
    const nx=Math.max(1,Math.ceil(width/spacing)),ny=Math.max(1,Math.ceil(height/spacing));
    if(nx*ny>3000000)throw Error('The 5 km comparison grid is too large. Draw a smaller study domain.');
    const points=[],cells=[];
    for(let row=0;row<ny;row++)for(let column=0;column<nx;column++) {
      const x=width<spacing?(xmin+xmax)/2:xmin+spacing/2+column*spacing;
      const y=height<spacing?(ymin+ymax)/2:ymin+spacing/2+row*spacing;
      const point=[x,y];
      if(!contains(d,point))continue;
      points.push(point);cells.push([x-spacing/2,y-spacing/2,spacing,spacing]);
    }
    if(!points.length)throw Error('The study domain is too small or narrow to contain a 5 km evaluation-grid point.');
    return {points,cells,spacingKm:spacing,nx,ny};
  }
  function predictAt(points,samples,weights,p,progress=()=>{}) {
    const values=new Float64Array(points.length);
    for(let i=0;i<points.length;i++) {
      let z=p.mu;
      for(let j=0;j<samples.length;j++)z+=covariance(points[i],samples[j],p)*weights[j];
      if(!Number.isFinite(z))throw Error('Surface prediction overflowed. Reduce extreme parameter values.');
      values[i]=z;
      if(i%1000===0)progress(i/points.length);
    }
    return values;
  }
  function referenceSurface(field,d,p,progress=()=>{}) {
    const grid=evaluationGrid(d);
    const reference=predictAt(grid.points,field.points,field.referenceWeights,p,progress);
    return {...grid,reference};
  }
  function compareSurfaces(dense,fit,p,progress=()=>{}) {
    const prediction=predictAt(dense.points,fit.samples,fit.weights,p,progress);
    const residual=new Float64Array(prediction.length),rootSquared=new Float64Array(prediction.length);
    let sumSquared=0,sumRoot=0,sumResidual=0;
    for(let i=0;i<prediction.length;i++) {
      const r=prediction[i]-dense.reference[i];residual[i]=r;rootSquared[i]=Math.sqrt(r*r);
      sumSquared+=r*r;sumRoot+=rootSquared[i];sumResidual+=r;
    }
    return {prediction,residual,rootSquared,measuredEpsilon:Math.sqrt(sumSquared/prediction.length),mae:sumRoot/prediction.length,bias:sumResidual/prediction.length};
  }
  return {nextSeed,random,normal,project,inverse,domain,geometryDomain,inside,contains,locations,covariance,cholesky,solve,generate,recommend,infer,validate,MAX_POINTS,evaluationGrid,predictAt,referenceSurface,compareSurfaces};
}
if (typeof module !== 'undefined' && module.exports) module.exports = createSpatialEngineV3;
