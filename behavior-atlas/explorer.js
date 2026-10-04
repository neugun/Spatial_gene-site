(() => {
const $=s=>document.querySelector(s);
const map=$("#map"),ctx=map.getContext("2d"),tip=$("#tooltip"),motifList=$("#motifList"),clipGrid=$("#clipGrid");
const title=$("#selectionTitle"),eyebrow=$("#selectionEyebrow"),meta=$("#selectionMeta"),radiusInput=$("#radius"),radiusVal=$("#radiusVal");

const S={
  colour:"motif",motif:null,selected:null,center:null,radius:+radiusInput.value,
  zoom:1,panX:0,panY:0,page:0,mode:null,pointers:new Map(),pinch:null
};
let D=null,W=0,H=0,dpr=1,extent=null,baseScale=1,baseX=0,baseY=0,hover=null;

const clamp=(x,a,b)=>Math.max(a,Math.min(b,x));
const pct=x=>Math.round(100*(Number.isFinite(x)?x:0))+"%";
const lerp=(a,b,t)=>Math.round(a+(b-a)*t);
function mix(c1,c2,t){
  const A=c1.match(/\w\w/g).map(x=>parseInt(x,16)),B=c2.match(/\w\w/g).map(x=>parseInt(x,16));
  return "rgb("+lerp(A[0],B[0],t)+","+lerp(A[1],B[1],t)+","+lerp(A[2],B[2],t)+")";
}
function visiblePoints(){return S.motif?D.points.filter(p=>p.motif===S.motif):D.points;}
function colorOf(p){
  if(S.colour==="motif")return D.motifs[p.motif-1].color;
  if(S.colour==="social")return mix("#e8e8e8","#1b7837",clamp(p.social,0,1));
  if(S.colour==="dem")return mix("#ededed","#d97706",clamp(p.dem,0,1));
  return mix("#ececec","#542788",clamp((p.quality-.45)/.55,0,1));
}
function mapToScreen(x,y){
  const sx=baseX+(x-extent.x0)*baseScale,sy=baseY+(extent.y1-y)*baseScale;
  return [W/2+(sx-W/2)*S.zoom+S.panX,H/2+(sy-H/2)*S.zoom+S.panY];
}
function screenToMap(x,y){
  const sx=W/2+(x-W/2-S.panX)/S.zoom,sy=H/2+(y-H/2-S.panY)/S.zoom;
  return [extent.x0+(sx-baseX)/baseScale,extent.y1-(sy-baseY)/baseScale];
}
function resize(){
  const r=map.getBoundingClientRect();W=Math.max(100,r.width);H=Math.max(100,r.height);dpr=window.devicePixelRatio||1;
  map.width=Math.round(W*dpr);map.height=Math.round(H*dpr);ctx.setTransform(dpr,0,0,dpr,0,0);
  if(D){
    const xs=D.points.map(p=>p.x).concat(D.motifs.map(m=>m.x)),ys=D.points.map(p=>p.y).concat(D.motifs.map(m=>m.y));
    let x0=Math.min(...xs),x1=Math.max(...xs),y0=Math.min(...ys),y1=Math.max(...ys);
    const mx=(x1-x0)*.07,my=(y1-y0)*.07;x0-=mx;x1+=mx;y0-=my;y1+=my;
    extent={x0,x1,y0,y1};baseScale=Math.min((W-64)/(x1-x0),(H-64)/(y1-y0));
    baseX=(W-(x1-x0)*baseScale)/2;baseY=(H-(y1-y0)*baseScale)/2;
  }
  draw();
}
function draw(){
  if(!D||!extent)return;
  ctx.clearRect(0,0,W,H);ctx.fillStyle="#fff";ctx.fillRect(0,0,W,H);
  const pts=visiblePoints();

  if(S.motif){
    ctx.globalAlpha=.035;
    for(const p of D.points){
      const [x,y]=mapToScreen(p.x,p.y);if(x<-3||y<-3||x>W+3||y>H+3)continue;
      ctx.fillStyle=D.motifs[p.motif-1].color;ctx.fillRect(x,y,1.4,1.4);
    }
  }
  ctx.globalAlpha=1;
  for(const p of pts){
    const [x,y]=mapToScreen(p.x,p.y);if(x<-8||y<-8||x>W+8||y>H+8)continue;
    let rr=2.15;
    if(hover&&hover.id===p.id)rr=4.8;
    if(S.selected&&S.selected.id===p.id)rr=6.2;
    ctx.fillStyle=colorOf(p);ctx.globalAlpha=(hover&&hover.id===p.id)?1:.68;
    ctx.beginPath();ctx.arc(x,y,rr,0,Math.PI*2);ctx.fill();
    if(S.selected&&S.selected.id===p.id){ctx.globalAlpha=1;ctx.strokeStyle="#111";ctx.lineWidth=1.8;ctx.stroke();}
  }
  ctx.globalAlpha=1;
  ctx.font="600 12px Inter,system-ui,sans-serif";ctx.textAlign="center";ctx.textBaseline="middle";
  for(const m of D.motifs){
    const [x,y]=mapToScreen(m.x,m.y);if(x<0||y<0||x>W||y>H)continue;
    ctx.globalAlpha=(!S.motif||S.motif===m.id)?.95:.16;
    ctx.lineWidth=4;ctx.strokeStyle="#fff";ctx.strokeText(m.name,x,y-11);
    ctx.fillStyle=m.color;ctx.fillText(m.name,x,y-11);
  }
  ctx.globalAlpha=1;
  if(S.center){
    const [x,y]=mapToScreen(S.center.x,S.center.y),rp=S.radius*baseScale*S.zoom;
    ctx.strokeStyle="#111";ctx.lineWidth=1.7;ctx.setLineDash([5,4]);ctx.beginPath();ctx.arc(x,y,rp,0,Math.PI*2);ctx.stroke();ctx.setLineDash([]);
    ctx.fillStyle="rgba(20,20,20,.035)";ctx.beginPath();ctx.arc(x,y,rp,0,Math.PI*2);ctx.fill();
    ctx.fillStyle="#111";ctx.beginPath();ctx.arc(x,y,3,0,Math.PI*2);ctx.fill();
  }
}
function nearestPoint(sx,sy,limit=12){
  let best=null,bd=limit*limit;
  for(const p of visiblePoints()){
    const [x,y]=mapToScreen(p.x,p.y),dd=(x-sx)**2+(y-sy)**2;
    if(dd<bd){bd=dd;best=p;}
  }
  return best;
}
function pointsInCircle(){
  if(!S.center)return [];
  return visiblePoints().map(p=>({p,d:Math.hypot(p.x-S.center.x,p.y-S.center.y)}))
    .filter(o=>o.d<=S.radius).sort((a,b)=>a.d-b.d);
}
function selectPoint(p){
  S.selected=p;S.center={x:p.x,y:p.y};S.page=0;draw();refreshClips();
}
function selectCenter(x,y){
  S.center={x,y};S.selected=null;S.page=0;
  const [sx,sy]=mapToScreen(x,y),p=nearestPoint(sx,sy,9999);
  if(p)S.selected=p;
  draw();refreshClips();
}
function segmentVideo(p){
  const video=document.createElement("video");
  video.src=D.video;video.muted=true;video.playsInline=true;video.controls=true;video.preload="metadata";
  const half=(D.video_context_seconds||2.4)/2;
  const start=Math.max(0,p.rep/D.fps-half),end=start+(D.video_context_seconds||2.4);
  const seek=()=>{try{if(Math.abs(video.currentTime-start)>.05)video.currentTime=start;}catch(_){}};
  video.addEventListener("loadedmetadata",seek,{once:true});
  video.addEventListener("loadeddata",seek,{once:true});
  setTimeout(()=>{if(video.readyState>=1)seek();},120);
  video.addEventListener("play",()=>{if(video.currentTime<start-.05||video.currentTime>end+.05)seek();});
  video.addEventListener("timeupdate",()=>{if(video.currentTime>=end)seek();});
  video.addEventListener("mouseenter",()=>video.play().catch(()=>{}));
  video.addEventListener("mouseleave",()=>video.pause());
  video.dataset.segStart=start.toFixed(3);video.dataset.segEnd=end.toFixed(3);
  return video;
}
function refreshClips(){
  clipGrid.innerHTML="";
  if(!S.center){title.textContent="Click anywhere in behavioral space";eyebrow.textContent="Selected neighborhood";meta.textContent="";return;}
  let near=pointsInCircle();
  if(!near.length&&S.selected)near=[{p:S.selected,d:0}];
  const total=near.length,start=total?((S.page*4)%total):0,chosen=[];
  for(let i=0;i<Math.min(4,total);i++)chosen.push(near[(start+i)%total].p);
  const dominant=chosen.length?chosen[0].motif:null;
  eyebrow.textContent=S.selected?"Exact bout + neighbors":"Selected neighborhood";
  title.textContent=S.selected?D.motifs[S.selected.motif-1].name:"Behavioral-space neighborhood";
  meta.textContent=total+" bout"+(total===1?"":"s")+" within radius "+S.radius;
  for(const p of chosen){
    const card=document.createElement("article");card.className="clip-card";
    const v=segmentVideo(p),cap=document.createElement("div");cap.className="clip-caption";
    cap.innerHTML="<strong>"+D.motifs[p.motif-1].name+" · frame "+p.rep+"</strong>"+
      p.duration.toFixed(1)+" s run · Social "+pct(p.social)+" · DEM feed "+pct(p.dem)+" · QC "+p.quality.toFixed(2);
    card.append(v,cap);clipGrid.append(card);
  }
  observeVideos();
}
function renderMotifs(){
  motifList.innerHTML="";
  const all=document.createElement("button");all.className="motif-btn"+(S.motif===null?" active":"");
  all.innerHTML='<span class="sw" style="background:#111"></span><span class="name">All motifs</span><span class="n">'+D.points.length+"</span>";
  all.onclick=()=>{S.motif=null;S.selected=null;S.center=null;S.page=0;renderMotifs();draw();refreshClips();};
  motifList.append(all);
  for(const m of D.motifs){
    const b=document.createElement("button");b.className="motif-btn"+(S.motif===m.id?" active":"");b.title=m.description||m.name;
    b.innerHTML='<span class="sw" style="background:'+m.color+'"></span><span class="name">'+String(m.id).padStart(2,"0")+" · "+m.name+'</span><span class="n">'+m.runs+"</span>";
    b.onclick=()=>{S.motif=(S.motif===m.id?null:m.id);S.selected=null;S.center={x:m.x,y:m.y};S.page=0;renderMotifs();draw();refreshClips();};
    motifList.append(b);
  }
}
function renderCards(){
  const box=$("#motifCards");box.innerHTML="";
  for(const m of D.motifs){
    const card=document.createElement("article");card.className="motif-card";
    card.innerHTML='<div class="motif-card-top"><span class="sw" style="background:'+m.color+'"></span><h3>'+String(m.id).padStart(2,"0")+" · "+m.name+'</h3></div>'+
      '<div class="desc">'+m.description+'<br>'+m.runs.toLocaleString()+" exact runs · "+m.frames.toLocaleString()+" frames</div>"+
      '<div class="stats"><div class="stat"><div class="v">'+m.median_duration.toFixed(1)+'s</div><div class="k">median run</div></div>'+
      '<div class="stat"><div class="v">'+pct(m.social_mean)+'</div><div class="k">social obs</div></div>'+
      '<div class="stat"><div class="v">'+pct(m.dem_mean)+'</div><div class="k">DEM feed</div></div></div>';
    card.onclick=()=>{S.motif=m.id;S.selected=null;S.center={x:m.x,y:m.y};S.page=0;renderMotifs();draw();refreshClips();$("#explorer").scrollIntoView({behavior:"smooth",block:"start"});};
    box.append(card);
  }
}
function observeVideos(){
  if(!("IntersectionObserver" in window))return;
  const obs=new IntersectionObserver(es=>es.forEach(e=>{if(e.isIntersecting)e.target.play().catch(()=>{});else e.target.pause();}),{threshold:.7});
  document.querySelectorAll(".clip-card video").forEach(v=>obs.observe(v));
}
function localXY(e){const r=map.getBoundingClientRect();return [e.clientX-r.left,e.clientY-r.top];}
function circleHit(sx,sy){
  if(!S.center)return false;
  const [cx,cy]=mapToScreen(S.center.x,S.center.y),rp=S.radius*baseScale*S.zoom;
  return Math.hypot(sx-cx,sy-cy)<=rp;
}

map.addEventListener("mousemove",e=>{
  const [x,y]=localXY(e);hover=nearestPoint(x,y,10);draw();
  if(hover){
    const m=D.motifs[hover.motif-1];tip.style.display="block";tip.style.left=Math.min(W-225,x+14)+"px";tip.style.top=Math.max(6,y-58)+"px";
    tip.innerHTML="<strong>"+m.name+"</strong><br>frame "+hover.rep+" · "+hover.duration.toFixed(1)+" s<br>Social "+pct(hover.social)+" · DEM "+pct(hover.dem)+" · QC "+hover.quality.toFixed(2);
  }else tip.style.display="none";
});
map.addEventListener("mouseleave",()=>{hover=null;tip.style.display="none";draw();});
map.addEventListener("wheel",e=>{
  e.preventDefault();const [x,y]=localXY(e),old=S.zoom,nu=clamp(old*Math.exp(-e.deltaY*.0014),1,12);if(nu===old)return;
  const mx=(x-W/2-S.panX)/old,my=(y-H/2-S.panY)/old;S.zoom=nu;S.panX=x-W/2-mx*nu;S.panY=y-H/2-my*nu;draw();
},{passive:false});

map.addEventListener("pointerdown",e=>{
  map.setPointerCapture&&map.setPointerCapture(e.pointerId);
  const [x,y]=localXY(e);S.pointers.set(e.pointerId,{x,y});
  if(S.pointers.size===2){
    const a=[...S.pointers.values()],dx=a[1].x-a[0].x,dy=a[1].y-a[0].y;
    S.pinch={dist:Math.hypot(dx,dy),zoom:S.zoom,midX:(a[0].x+a[1].x)/2,midY:(a[0].y+a[1].y)/2,panX:S.panX,panY:S.panY};
    S.mode="pinch";return;
  }
  const p=nearestPoint(x,y,11);
  if(p){selectPoint(p);S.mode="point";return;}
  if(circleHit(x,y)){S.mode="circle";S.dragX=x;S.dragY=y;return;}
  S.mode="pan";S.dragX=x;S.dragY=y;S.moved=false;
});
map.addEventListener("pointermove",e=>{
  if(!S.pointers.has(e.pointerId))return;
  const [x,y]=localXY(e),old=S.pointers.get(e.pointerId);S.pointers.set(e.pointerId,{x,y});
  if(S.mode==="pinch"&&S.pointers.size>=2&&S.pinch){
    const a=[...S.pointers.values()].slice(0,2),dist=Math.hypot(a[1].x-a[0].x,a[1].y-a[0].y),midX=(a[0].x+a[1].x)/2,midY=(a[0].y+a[1].y)/2;
    const nu=clamp(S.pinch.zoom*(dist/Math.max(1,S.pinch.dist)),1,12);
    const mx=(S.pinch.midX-W/2-S.pinch.panX)/S.pinch.zoom,my=(S.pinch.midY-H/2-S.pinch.panY)/S.pinch.zoom;
    S.zoom=nu;S.panX=midX-W/2-mx*nu;S.panY=midY-H/2-my*nu;draw();return;
  }
  if(S.mode==="circle"){
    const m=screenToMap(x,y);S.center={x:m[0],y:m[1]};S.selected=null;S.page=0;draw();refreshClips();return;
  }
  if(S.mode==="pan"){
    const dx=x-S.dragX,dy=y-S.dragY;S.panX+=dx;S.panY+=dy;S.dragX=x;S.dragY=y;if(Math.hypot(x-old.x,y-old.y)>1)S.moved=true;draw();
  }
});
function pointerEnd(e){
  const [x,y]=localXY(e);S.pointers.delete(e.pointerId);
  if(S.mode==="pan"&&!S.moved){const m=screenToMap(x,y);selectCenter(m[0],m[1]);}
  if(S.pointers.size<2&&S.mode==="pinch")S.mode=null;
  else if(S.pointers.size===0)S.mode=null;
}
map.addEventListener("pointerup",pointerEnd);map.addEventListener("pointercancel",pointerEnd);

radiusInput.addEventListener("input",()=>{S.radius=+radiusInput.value;radiusVal.textContent=S.radius;draw();refreshClips();});
$("#moreBtn").onclick=()=>{if(!S.center)return;S.page++;refreshClips();};
$("#resetBtn").onclick=()=>{Object.assign(S,{motif:null,selected:null,center:null,zoom:1,panX:0,panY:0,page:0});renderMotifs();draw();refreshClips();};
$("#fullBtn").onclick=()=>{const ex=$("#explorer");if(!document.fullscreenElement)ex.requestFullscreen&&ex.requestFullscreen();else document.exitFullscreen&&document.exitFullscreen();};
$("#colourBtns").addEventListener("click",e=>{const b=e.target.closest("button");if(!b)return;S.colour=b.dataset.v;document.querySelectorAll("#colourBtns button").forEach(x=>x.classList.toggle("active",x===b));draw();});
window.addEventListener("resize",resize);document.addEventListener("fullscreenchange",()=>setTimeout(resize,50));

fetch("data/atlas.json").then(r=>r.json()).then(data=>{
  D=data;renderMotifs();renderCards();resize();
  const first=D.points.find(p=>p.id===D.default_point_id)||D.points[0];if(first)selectPoint(first);
}).catch(err=>{console.error(err);title.textContent="Could not load atlas data";});
})();