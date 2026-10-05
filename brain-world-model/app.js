const $ = s => document.querySelector(s);
const fmt = (x,d=3) => x==null || !Number.isFinite(+x) ? "—" : (+x).toFixed(d);
const pct = (x,d=1) => x==null || !Number.isFinite(+x) ? "—" : ((+x*100).toFixed(d) + "%");
const compact = n => { n=+n; if(n>=1e9)return (n/1e9).toFixed(2)+"B"; if(n>=1e6)return (n/1e6).toFixed(1)+"M"; if(n>=1e3)return (n/1e3).toFixed(1)+"K"; return String(n); };
let DATA=null;

$("#langToggle").addEventListener("click", function(){
  document.body.classList.toggle("cn");
  $("#langToggle").textContent = document.body.classList.contains("cn") ? "EN" : "中文";
});

function renderDatasets(filter){
  filter=filter||"all";
  const grid=$("#datasetGrid"); grid.innerHTML="";
  DATA.datasets.filter(function(d){
    return filter==="all" || (d.tags||[]).map(x=>x.toLowerCase()).includes(filter);
  }).forEach(function(d){
    const el=document.createElement("article"); el.className="dataset-card";
    const mods=Array.isArray(d.modality)?d.modality:[d.modality].filter(Boolean);
    const tags=[...new Set([...(d.tags||[]),...mods])].slice(0,6).map(t=>'<span class="tag">'+t+'</span>').join("");
    el.innerHTML='<div class="source">'+(d.source||"public dataset")+' · '+(d.species||"—")+' <span class="pass">● PASS</span></div><h4>'+d.name+'</h4><div class="tag-row">'+tags+'</div>';
    grid.appendChild(el);
  });
}
$("#datasetFilters").addEventListener("click", function(e){
  const b=e.target.closest(".chip"); if(!b)return;
  document.querySelectorAll(".chip").forEach(x=>x.classList.remove("active"));
  b.classList.add("active"); renderDatasets(b.dataset.filter);
});

function renderSummary(){
  const s=DATA.summary;
  $("#metricDatasets").textContent=s.dataset_contract_pass+"/"+s.dataset_contract_total;
  $("#metricArch").textContent=s.architecture_smoke_pass;
  $("#metricApps").textContent=s.application_families;
  $("#metricTransfer").textContent="+"+pct(DATA.unified_token.best_heldout_future_skill);
  $("#heldoutSkill").textContent="+"+fmt(DATA.unified_token.best_heldout_future_skill,3);
  const cards=[
    [s.dataset_contract_pass+"/"+s.dataset_contract_total,"token contracts"],
    [s.forecast_sanity_total,"forecast sanity tests"],
    [s.architecture_smoke_pass,"model / checkpoint tests"],
    [s.application_families,"application families"]
  ];
  $("#testSummary").innerHTML=cards.map(x=>'<div class="summary-card"><b>'+x[0]+'</b><span>'+x[1]+'</span></div>').join("");
}
function renderTransfer(){
  $("#transferCards").innerHTML=DATA.unified_token.heldout_future.map(function(r){
    const pos=r.skill_vs_zero>=0;
    return '<article class="transfer-card"><small>'+r.modality+' · held-out source</small><h4>'+r.view.replaceAll("_"," ")+'</h4><div class="skill '+(pos?"pos":"neg")+'">'+(pos?"+":"")+fmt(r.skill_vs_zero,3)+'</div><small>future skill vs zero · MSE '+fmt(r.mse,3)+' vs '+fmt(r.zero_mse,3)+'</small></article>';
  }).join("");
}
function renderArch(){
  const rows=Object.entries(DATA.architecture_smokes);
  $("#architectureTable").innerHTML='<table><thead><tr><th>model</th><th>status</th><th>linked smoke params</th><th>best val avg</th><th>note</th></tr></thead><tbody>'+
  rows.map(function(kv){const k=kv[0],v=kv[1];return '<tr><td>'+k.replaceAll("_"," ")+'</td><td><span class="status-dot"></span>'+v.status+'</td><td>'+compact(v.linked_smoke_params||0)+'</td><td>'+fmt(v.best_val_avg,4)+'</td><td>'+(v.note||"")+'</td></tr>';}).join("")+'</tbody></table>';
}
function renderForecast(){
  const sel=$("#forecastDataset");
  sel.innerHTML=DATA.generic_forecast_sanity.map((r,i)=>'<option value="'+i+'">'+r.dataset+'</option>').join("");
  function draw(){
    const r=DATA.generic_forecast_sanity[+sel.value||0];
    const vals=[["ridge",r.ridge_r2],["persistence",r.persistence_r2],["mean",r.mean_baseline_r2]];
    const scale=Math.max(.25,...vals.map(x=>Math.min(2,Math.abs(x[1]||0))));
    $("#forecastBars").innerHTML=vals.map(function(v){
      const name=v[0],val=v[1],w=Math.min(50,Math.abs(val)/scale*50),pos=val>=0;
      return '<div class="bar-row"><span>'+name+'</span><div class="bar-track"><i class="bar-zero"></i><i class="bar-fill '+(pos?"pos":"neg")+'" style="width:'+w+'%"></i></div><span class="bar-value">'+fmt(val,3)+'</span></div>';
    }).join("")+(Math.abs(r.ridge_r2)>2?'<div class="note">Bar length capped for readability; numeric value is uncapped.</div>':"");
  }
  sel.addEventListener("change",draw); draw();
}
function rowsHtml(rows){return rows.map(x=>'<div class="result-row"><span>'+x[0]+'</span><b class="'+(x[2]||"")+'">'+x[1]+'</b></div>').join("");}
function renderApps(){
  const sl=DATA.applications.sleep_state_decoding;
  $("#sleepApp").innerHTML=Object.entries(sl).map(function(kv){
    const name=kv[0],v=kv[1],label=name.includes("osf")?"Mouse AccuSleep":"Human Sleep-EDFx";
    return '<div class="result-row"><span>'+label+' balanced accuracy</span><b>'+pct(v.balanced_accuracy)+'</b></div><div class="result-row"><span>majority baseline</span><b>'+pct(v.majority_baseline)+'</b></div>';
  }).join("");
  const m=DATA.applications.mcbride_real_optogenetic;
  $("#mcbrideApp").innerHTML=rowsHtml([
    ["pre-state only R²",fmt(m.pre_only_r2,3)],
    ["pre-state + intervention R²",fmt(m.pre_plus_intervention_r2,3)],
    ["ΔR² from intervention","+"+fmt(m.delta_r2_action,3),"delta"],
    ["held-out opto trials",m.n_test]
  ]);
  const sy=DATA.applications.synthetic_intervention.horizons;
  $("#interventionApp").innerHTML=Object.entries(sy).map(kv=>'<div class="result-row"><span>horizon '+kv[0]+' · ΔR²(action)</span><b class="delta">+'+fmt(kv[1].delta_r2,3)+'</b></div>').join("");
  const om=DATA.applications.omnimouse_neural_to_behavior.r2;
  $("#omniApp").innerHTML=rowsHtml([
    ["treadmill R²",fmt(om.treadmill,3),om.treadmill>0?"delta":""],
    ["pupil radius R²",fmt(om.pupil_radius,3)],
    ["eye x R²",fmt(om.eye_x,3)],
    ["interpretation","treadmill signal; eye unresolved"]
  ]);
}
function renderObjectives(){
  const pretty={causal_world:"Brain + behavior future",causal_neural:"Neural future",causal_behavior:"Behavior future"};
  $("#objectiveGrid").innerHTML=DATA.three_model_future_eval.comparisons.map(function(c){
    const models=c.models||{}, imp=c.relative_improvement_vs_baseline_pct||{};
    function line(m){
      const v=models[m]&&models[m]["val/avg"];
      const delta=m==="baseline"?null:(imp[m]&&imp[m]["val/avg"]);
      const cls=delta>0?"better":delta<0?"worse":"";
      const label=delta==null?"reference":((delta>0?"+":"")+fmt(delta,1)+"%");
      return '<div class="obj-row"><span>'+m+'</span><b>'+fmt(v,3)+'</b><span class="'+cls+'">'+label+'</span></div>';
    }
    return '<article class="objective-card"><h4>'+(pretty[c.scheme]||c.scheme)+'</h4>'+line("baseline")+line("world")+line("predictive")+'</article>';
  }).join("");
}
function renderAudit(){
  $("#paramAudit").innerHTML=Object.entries(DATA.parameter_audit).map(function(kv){
    const name=kv[0],v=kv[1],ratio=v.shared_backbone_like/v.total;
    return '<article class="audit-card"><div class="audit-title"><h4>'+name.replaceAll("_"," ")+'</h4><span class="audit-total">'+compact(v.total)+' linked total</span></div><div class="ratio-bar"><div class="ratio-shared" style="width:'+Math.max(.5,ratio*100)+'%"></div></div><div class="audit-stats"><span><b>'+compact(v.shared_backbone_like)+'</b> shared ('+pct(ratio)+')</span><span>'+compact(v.session_or_stitcher)+' session/stitcher</span></div></article>';
  }).join("");
}
async function init(){
  try{
    const resp=await fetch("./results.json",{cache:"no-store"}); if(!resp.ok)throw Error(resp.status);
    DATA=await resp.json();
    renderSummary();renderDatasets();renderTransfer();renderArch();renderForecast();renderApps();renderObjectives();renderAudit();
    $("#updatedAt").textContent="results updated "+DATA.updated_at;
  }catch(e){console.error(e);$("#updatedAt").textContent="results failed to load";}
}
init();
