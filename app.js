(function () {
  "use strict";
  const quantile = (sorted, p) => {
    if (!sorted.length) return null;
    const k = (sorted.length - 1) * p, i = Math.floor(k);
    return sorted[i] + (sorted[Math.min(i + 1, sorted.length - 1)] - sorted[i]) * (k - i);
  };
  function summarise(rows, day, cutoff, minRead) {
    const all = rows.filter(r => r[1] === day);
    const kept = all.filter(r => r[3] >= minRead);
    const values = kept.map(r => r[2]).sort((a,b) => a-b);
    const n = values.length;
    return {day, sourceN:all.length, n, values, mean:n ? values.reduce((a,b)=>a+b,0)/n : null,
      median:quantile(values,.5), q1:quantile(values,.25), q3:quantile(values,.75),
      shortN:values.filter(v=>v<cutoff).length,
      short:n ? values.filter(v=>v<cutoff).length/n : null};
  }
  function contrast(rows, from, to, cutoff, minRead) {
    const a=summarise(rows,from,cutoff,minRead), b=summarise(rows,to,cutoff,minRead);
    return {a,b,mean:a.n&&b.n?b.mean-a.mean:null,median:a.n&&b.n?b.median-a.median:null,
      shortPP:a.n&&b.n?(b.short-a.short)*100:null};
  }
  function reviewRows(rows, from, to, threshold) {
    return rows.map(r=>{
      const a=r["gap"+from],b=r["gap"+to];
      return {row:r,a,b,delta:a&&b?b.tl_bp-a.tl_bp:null,missing:!a||!b};
    }).filter(r=>(r.a||r.b)&&(r.missing||Math.abs(r.delta)>=threshold)).sort((a,b)=>
      (b.missing-a.missing) || Math.abs(b.delta)-Math.abs(a.delta) || a.row.id.localeCompare(b.row.id));
  }
  function csv(rows) {
    return rows.map(r=>r.map(v=>{
      const s=v===null||v===undefined?"":String(v);
      const safe=typeof v==="string"&&/^[=+@-]/.test(s)?"'"+s:s;
      return '"'+safe.replace(/"/g,'""')+'"';
    }).join(",")).join("\r\n")+"\r\n";
  }
  if (typeof module!=="undefined" && module.exports) {
    module.exports={quantile,summarise,contrast,reviewRows,csv}; return;
  }
  const $=id=>document.getElementById(id);
  const fmt=(v,d=0)=>v===null?"—":Number(v).toLocaleString("en-US",{minimumFractionDigits:d,maximumFractionDigits:d});
  const signed=(v,d=1)=>v===null?"—":(v>0?"+":"")+fmt(v,d);
  const svgNS="http://www.w3.org/2000/svg";
  function element(tag,attrs={},content) {
    const e=document.createElementNS(svgNS,tag);
    Object.entries(attrs).forEach(([k,v])=>e.setAttribute(k,String(v)));
    if(content!==undefined)e.textContent=content;
    return e;
  }
  function text(svg,x,y,value,attrs={}) {svg.appendChild(element("text",{x,y,...attrs},value));}
  function line(svg,x1,y1,x2,y2,attrs={}) {svg.appendChild(element("line",{x1,y1,x2,y2,stroke:"#ddd",...attrs}));}
  function resetSvg(svg,title,description) {
    svg.replaceChildren(element("title",{},title),element("desc",{},description));
    svg.removeAttribute("aria-labelledby"); svg.setAttribute("aria-label",description);
  }
  function drawCDF(c,cutoff) {
    const svg=$("hero-viz"), left=49, right=579, top=23, bottom=290;
    const max=Math.ceil(Math.max(...c.a.values,...c.b.values,cutoff)/2000)*2000;
    const x=v=>left+v/max*(right-left), y=v=>bottom-v*(bottom-top);
    resetSvg(svg,"Cumulative telomere length distributions",
      "Day "+c.a.day+" and day "+c.b.day+" cumulative fractions by telomere length in base pairs. Full retained range is shown. Exact summary values follow in the table.");
    [0,.25,.5,.75,1].forEach(t=>{line(svg,left,y(t),right,y(t));text(svg,left-9,y(t)+4,fmt(t*100)+"%",{"text-anchor":"end"});});
    for(let t=0;t<=max;t+=max/4){text(svg,x(t),311,fmt(t/1000,Number.isInteger(t/1000)?0:1),{"text-anchor":"middle"});}
    text(svg,(left+right)/2,335,"Telomere length (kb) · full retained range",{"text-anchor":"middle"});
    [c.a,c.b].forEach((s,i)=>{
      let path="M "+x(0)+" "+y(0);
      s.values.forEach((v,j)=>{path+=" H "+x(v)+" V "+y((j+1)/s.n);});
      path+=" H "+x(max);
      svg.appendChild(element("path",{d:path,fill:"none",stroke:i?"#1f7a8c":"#999","stroke-width":2}));
    });
    line(svg,x(cutoff),top,x(cutoff),bottom,{stroke:"#222","stroke-dasharray":"4 4"});
    text(svg,Math.min(right-75,x(cutoff)+6),top+12,fmt(cutoff/1000,2)+" kb cutoff");
  }
  function download(name,rows,status) {
    const blob=new Blob([csv(rows)],{type:"text/csv;charset=utf-8"});
    const url=URL.createObjectURL(blob),a=document.createElement("a");
    a.href=url;a.download=name;document.body.appendChild(a);a.click();a.remove();
    setTimeout(()=>URL.revokeObjectURL(url),1000);
    if(status)status.textContent="Downloaded "+name+". No information was uploaded.";
  }
  async function load(url) {
    const r=await fetch(url);if(!r.ok)throw new Error("Data file returned HTTP "+r.status);return r.json();
  }
  async function start() {
    const config=window.UNIVERSITY_LAB_DEMO_DATA;
    const [primary,audit]=await Promise.all([load(config.primary),load(config.audit)]);
    const state={comparison:null,cutoff:2000,minRead:0,auditRows:[],fromGap:20,toGap:250};

    function renderPrimary() {
      const [from,to]=$("pair-select").value.split(",").map(Number);
      const cutoff=Number($("threshold").value),minRead=Number($("min-read").value);
      const c=contrast(primary.rows,from,to,cutoff,minRead), base=contrast(primary.rows,from,to,cutoff,0);
      state.comparison=c;state.cutoff=cutoff;state.minRead=minRead;
      $("threshold-label").textContent=fmt(cutoff)+" bp";
      $("cdf-legend").replaceChildren(Object.assign(document.createElement("span"),{textContent:"Day "+from}),
        Object.assign(document.createElement("span"),{textContent:"Day "+to}));
      drawCDF(c,cutoff);
      const direction=[c.mean,c.median,-c.shortPP].map(Math.sign);
      const disagreement=direction.some(v=>v<0)&&direction.some(v=>v>0);
      const reversal=minRead>0&&["mean","median","shortPP"].some(k=>Math.sign(c[k])!==Math.sign(base[k])&&c[k]!==0&&base[k]!==0);
      let headline=reversal?"Read selection reverses at least one direction.":
        disagreement?"The readouts point in different directions.":"The readouts are directionally consistent at this setting.";
      let detail="Day "+from+" → "+to+": mean "+signed(c.mean)+" bp; median "+signed(c.median)+" bp; fraction below "+
        fmt(cutoff)+" bp "+signed(c.shortPP,2)+" percentage points.";
      if(minRead) detail+=" Only "+fmt(c.a.n)+" / "+fmt(c.a.sourceN)+" and "+fmt(c.b.n)+" / "+fmt(c.b.sourceN)+
        " rows remain. Selection can change the length distribution; this is not a QC pass.";
      $("readout").replaceChildren(Object.assign(document.createElement("strong"),{textContent:headline}),document.createTextNode(" "+detail));
      $("comparison-table").querySelector("tbody").replaceChildren(...[c.a,c.b].map(s=>{
        const tr=document.createElement("tr");
        ["Day "+s.day,fmt(s.n)+" / "+fmt(s.sourceN),fmt(s.mean,1),fmt(s.median,1),
          fmt(s.shortN)+" ("+fmt(s.short*100,2)+"%)"].forEach(v=>{const td=document.createElement("td");td.textContent=v;tr.appendChild(td);});
        return tr;
      }));
      $("metric-three").textContent=signed(c.shortPP,2)+" pp";
      $("experiment-plan").replaceChildren();
      const plan=[
        reversal? "First resolve the selection effect. The chosen read restriction changes the biological-looking direction relative to the deposited data. Compare library handling and full length distributions before treating this as a passage effect." :
        disagreement? "Carry all three endpoints into a repeat. A mean alone is not enough to call this comparison a plateau or reversal; the median and lower tail ask a different question about the distribution." :
        "Use this as a candidate passage comparison, not a confirmed effect. Repeat the same readout definitions across independent clones to test whether the direction is reproducible.",
        "Proposed comparison: day "+from+" and day "+to+" in TERT-loss clones and matched parental or mock-edited controls; preserve clone identity across passages. Record population doublings, extraction and capture batch, software version and explicit parameters.",
        "Discriminating check: use matched handling and a complementary assay such as TRF Southern blot to assess length change, recognizing its subtelomeric contribution and different resolution. Persistence across clones and measurement approaches supports a biological follow-up; a change only under one library or filter calls for measurement review."
      ];
      plan.forEach(p=>{$("experiment-plan").appendChild(Object.assign(document.createElement("p"),{textContent:p}));});
    }

    function renderAuditRead() {
      const item=state.auditRows.find(r=>r.row.id===$("read-select").value),svg=$("audit-viz");
      if(!item) {resetSvg(svg,"No matching reads","No reads meet this review rule.");$("audit-readout").textContent="No reads meet this rule; this does not validate the measurement.";return;}
      const {a,b,row,delta}=item, max=Math.ceil(Math.max(a?.tl_bp||0,b?.tl_bp||0,1000)/1000)*1000;
      const left=94,right=568,x=v=>left+v/max*(right-left);
      resetSvg(svg,"Read-level parameter sensitivity",row.id+": "+
        state.fromGap+" bp gap gives "+(a?a.tl_bp+" bp":"no output")+"; "+state.toGap+" bp gap gives "+(b?b.tl_bp+" bp":"no output")+". Bar lengths represent measured telomere lengths, not aligned sequence coordinates.");
      [a,b].forEach((v,i)=>{
        const yy=40+i*65;
        text(svg,left-10,yy+5,(i?state.toGap:state.fromGap)+" bp gap",{"text-anchor":"end"});
        if(v) {
          line(svg,left,yy,x(v.tl_bp),yy,{stroke:i?"#1f7a8c":"#999","stroke-width":9});
          text(svg,x(v.tl_bp),yy+25,fmt(v.tl_bp)+" bp",{"text-anchor":x(v.tl_bp)>right-50?"end":"middle"});
        } else text(svg,left,yy+5,"No retained output");
      });
      text(svg,(left+right)/2,178,"Measured telomere length · not read-coordinate geometry",{"text-anchor":"middle"});
      const info=v=>v?[v.chr,v.arm,v.direction,"MAPQ "+v.mapq,"read "+fmt(v.read_bp)+" bp"].join(" · "):"not retained by this run";
      $("audit-readout").textContent=row.id+" ("+row.read_id+"): "+(delta===null?"retention differs between settings":signed(delta,0)+" bp change")+
        ". "+state.fromGap+" bp: "+info(a)+". "+state.toGap+" bp: "+info(b)+
        ". Inspect the alignment and boundary before interpreting the difference; neither the change nor MAPQ alone proves an error.";
    }
    function renderAudit() {
      const [from,to]=$("gap-pair").value.split(",").map(Number),threshold=Number($("delta-cutoff").value);
      state.fromGap=from;state.toGap=to;
      state.auditRows=reviewRows(audit.rows,from,to,threshold);
      const previous=$("read-select").value, matched=state.auditRows.filter(r=>!r.missing).length,missing=state.auditRows.filter(r=>r.missing).length;
      $("read-select").replaceChildren(...state.auditRows.map(r=>{
        const op=document.createElement("option");op.value=r.row.id;
        op.textContent=r.row.id+" · "+(r.missing?"retention differs":signed(r.delta,0)+" bp");return op;
      }));
      if(state.auditRows.some(r=>r.row.id===previous))$("read-select").value=previous;
      else {
        const first=state.auditRows.find(r=>!r.missing)||state.auditRows[0];
        if(first)$("read-select").value=first.row.id;
      }
      const comp=audit.comparisons.find(r=>r.from_gap===from&&r.to_gap===to);
      $("audit-summary").textContent=fmt(matched)+" of "+fmt(comp.matched_n)+" matched read IDs meet |change| ≥ "+fmt(threshold)+
        " bp, plus "+missing+" with differing retention. Across all outputs, the median changes from "+
        fmt(audit.summaries[from].median_bp)+" to "+fmt(audit.summaries[to].median_bp)+" bp. Thresholds define a review list, not a validated QC cutoff.";
      renderAuditRead();
    }

    $("pair-select").addEventListener("change",renderPrimary);
    $("min-read").addEventListener("change",renderPrimary);
    $("threshold").addEventListener("input",renderPrimary);
    $("gap-pair").addEventListener("change",renderAudit);
    $("delta-cutoff").addEventListener("change",renderAudit);
    $("read-select").addEventListener("change",renderAuditRead);
    $("export-comparison").addEventListener("click",()=>{
      const s=state,c=s.comparison;
      const rows=[["source","from_day","to_day","short_cutoff_bp_strict_lt","min_read_bp_stress_test","from_n","to_n","mean_change_bp","median_change_bp","short_fraction_change_pp","interpretation"]];
      primary.days.forEach((a,i)=>primary.days.slice(i+1).forEach(b=>{
        const t=contrast(primary.rows,a,b,s.cutoff,s.minRead);
        rows.push(["hesc_tertko",a,b,s.cutoff,s.minRead,t.a.n,t.b.n,t.mean,t.median,t.shortPP,"Descriptive; four passage samples; no biological replication supplied"]);
      }));
      download("artandi-passage-comparisons.csv",rows,$("export-status"));
    });
    $("export-plan").addEventListener("click",()=>{
      const s=state,c=s.comparison;
      const rows=[["plan_status","clone_id","sample_id","genotype_or_control","planned_passage_day","population_doublings","extraction_batch","library_batch","telometer_version","explicit_maxgap_bp","short_cutoff_bp","orthogonal_assay_id","notes"]];
      ["TERT-loss","matched parental or mock-edited"].forEach(g=>[c.a.day,c.b.day].forEach(day=>
        rows.push(["PROPOSED_NOT_PERFORMED","","",g,day,"","","","","",s.cutoff,"","Duplicate rows for independent clones; confirm protocol and controls with supervising researcher. Do not transfer example-BAM settings without validation."])));
      download("artandi-blank-validation-log.csv",rows,$("export-status"));
    });
    $("export-audit").addEventListener("click",()=>{
      const rows=[["example_bam_read_id","display_id","from_gap_bp","to_gap_bp","from_telomere_bp","to_telomere_bp","change_bp","from_mapq","to_mapq","review_reason","software"]];
      state.auditRows.forEach(r=>rows.push([r.row.read_id,r.row.id,state.fromGap,state.toGap,r.a?.tl_bp,r.b?.tl_bp,r.delta,r.a?.mapq,r.b?.mapq,r.missing?"differing retention":"absolute change meets selected threshold","Telometer 2.0.3"]));
      download("artandi-telometer-review.csv",rows,$("export-status"));
    });
    renderPrimary();renderAudit();
    document.documentElement.dataset.ready="true";
  }
  start().catch(error=>{
    $("readout").textContent="The data could not be loaded. Please reload; no results have been substituted. "+error.message;
    $("audit-readout").textContent="Data unavailable.";
    document.documentElement.dataset.ready="error";
  });
}());
