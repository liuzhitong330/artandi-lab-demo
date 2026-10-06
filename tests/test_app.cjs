/* Independent Node checks for the exact arithmetic exposed by the interface.
 * Run: node tests/test_app.cjs
 * No browser or network is needed. The primary expected values are aggregated
 * directly from the preserved source rows and checked against the Python audit.
 */
"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const path = require("node:path");
const {summarise, contrast, reviewRows, csv} = require("../app.js");
const root = path.resolve(__dirname, "..");
const primary = JSON.parse(fs.readFileSync(path.join(root, "data/tertko.json"), "utf8"));
const checks = JSON.parse(fs.readFileSync(path.join(root, "data/analysis_checks.json"), "utf8"));
const audit = JSON.parse(fs.readFileSync(path.join(root, "data/telometer-audit.json"), "utf8"));
const near = (a, b, label) => {
  if (a === null || b === null) return assert.equal(a, b, label);
  assert.ok(Number.isFinite(a) && Number.isFinite(b), label + " must be finite");
  assert.ok(Math.abs(a - b) <= 1e-8, `${label}: ${a} != ${b}`);
};
const median = values => {
  if (!values.length) return null;
  const v = values.toSorted((a,b) => a-b), n = v.length;
  return n % 2 ? v[(n-1)/2] : (v[n/2-1] + v[n/2]) / 2;
};
const quartile = (values, p) => {
  if (!values.length) return null;
  const v = values.toSorted((a,b) => a-b), rank = 1 + (v.length - 1) * p;
  const lowerRank = Math.floor(rank), upperRank = Math.ceil(rank);
  return v[lowerRank - 1] * (upperRank === lowerRank ? 1 : upperRank-rank) +
    (upperRank === lowerRank ? 0 : v[upperRank - 1] * (rank-lowerRank));
};
function independentSummary(day, threshold, minRead) {
  const population = primary.rows.filter(row => row[1] === day);
  const selected = population.filter(row => row[3] >= minRead);
  const values = selected.map(row => row[2]);
  const count = values.length;
  const shortN = values.reduce((n, value) => n + Number(value < threshold), 0);
  return {n: count, sourceN: population.length, shortN,
    mean: count ? values.reduce((sum,value) => sum+value, 0)/count : null,
    median: median(values), q1:quartile(values,.25),q3:quartile(values,.75),
    short: count ? shortN / count : null,
    ids: selected.map(row=>row[0])};
}

assert.equal(primary.rows.length, 2610);
assert.equal(new Set(primary.rows.map(row => row[0])).size, 2610);
assert.equal(new Set(primary.rows.map(row => JSON.stringify(row.slice(1)))).size, 2585);
assert.deepEqual(primary.rows.map(row => row[0]), Array.from({length:2610}, (_,i) => i+3));
const pairs = [[66,78],[66,98],[66,105],[78,98],[78,105],[98,105]];
let cases = 0, pythonComparisons = 0;
for (const [from,to] of pairs) {
  for (let cutoff = 1000; cutoff <= 4000; cutoff += 250) {
    for (const minRead of [0,4000,8000]) {
      const expectedA = independentSummary(from,cutoff,minRead);
      const expectedB = independentSummary(to,cutoff,minRead);
      const actual = contrast(primary.rows,from,to,cutoff,minRead);
      for (const [name,expected] of [["a",expectedA],["b",expectedB]]) {
        const current = actual[name];
        for (const key of ["n","sourceN","shortN"]) assert.equal(current[key],expected[key], key);
        for (const key of ["mean","median","q1","q3","short"]) near(current[key],expected[key],key);
        assert.equal(new Set(expected.ids).size, expected.n, "retained source-row IDs are unique");
        assert.equal(current.values.length, expected.n);
        assert.ok(current.values.every((value,i) => !i || current.values[i-1] <= value), "CDF data must be ordered");
      }
      near(actual.mean,expectedB.mean-expectedA.mean,"mean change");
      near(actual.median,expectedB.median-expectedA.median,"median change");
      near(actual.shortPP,100*(expectedB.short-expectedA.short),"short-fraction change");
      const precomputed = checks.contrasts.find(row => row.earlier_day === from && row.later_day === to &&
        row.threshold_bp === cutoff && row.minimum_read_bp === minRead);
      if (precomputed) {
        near(actual.mean,precomputed.mean_change_bp,"Python mean cross-check");
        near(actual.median,precomputed.median_change_bp,"Python median cross-check");
        near(actual.shortPP,precomputed.short_fraction_change_pp,"Python short-tail cross-check");
        assert.equal(actual.a.n,precomputed.earlier_n);
        assert.equal(actual.b.n,precomputed.later_n);
        pythonComparisons++;
      }
      cases++;
    }
  }
}
assert.equal(cases,234);
assert.equal(pythonComparisons,108);

// Boundary/empty behavior: exact threshold equality is not short. Read equality
// is retained, and zero is not a missing value. An empty sample is not mean 0.
const edge = [[3,1,0,4000,"chr1","p","fwd"],[4,1,2000,3999,"chr1","p","fwd"],
  [5,1,2001,4000,"chr1","p","fwd"]];
const zero = summarise(edge,1,2000,0);
assert.equal(zero.shortN,1);
assert.equal(summarise(edge,1,2000,4000).n,2);
const empty = summarise(edge,2,2000,0);
for (const field of ["mean","median","q1","q3","short"]) assert.equal(empty[field],null);
assert.equal(empty.n,0);
assert.equal(empty.shortN,0);
const emptyContrast = contrast(edge,1,2,2000,0);
for (const field of ["mean","median","shortPP"]) assert.equal(emptyContrast[field],null);

assert.equal(audit.rows.length,738);
assert.equal(new Set(audit.rows.map(row=>row.id)).size,738);
assert.equal(new Set(audit.rows.map(row=>row.read_id)).size,738);
for (const gap of [20,100,250]) {
  const values = audit.rows.flatMap(row=>row[`gap${gap}`] ? [row[`gap${gap}`].tl_bp] : []);
  const source = audit.summaries[gap];
  assert.equal(values.length,source.n);
  near(median(values),source.median_bp,`gap${gap} median`);
  near(values.reduce((a,b)=>a+b,0)/values.length,source.mean_bp,`gap${gap} mean`);
  near(quartile(values,.25),source.q1_bp,`gap${gap} q1`);
  near(quartile(values,.75),source.q3_bp,`gap${gap} q3`);
}
let auditCases = 0;
const reviewCounts = [];
for (const [from,to] of [[20,100],[20,250],[100,250]]) {
  const matched = audit.rows.filter(row=>row[`gap${from}`] && row[`gap${to}`]);
  const retentionChanged = audit.rows.filter(row=>Boolean(row[`gap${from}`]) !== Boolean(row[`gap${to}`]));
  const source = audit.comparisons.find(row=>row.from_gap===from && row.to_gap===to);
  assert.equal(matched.length,source.matched_n);
  assert.equal(retentionChanged.length,source.from_only_n+source.to_only_n);
  for (const threshold of [100,250,500,1000]) {
    const differences = matched.filter(row=>Math.abs(row[`gap${to}`].tl_bp-row[`gap${from}`].tl_bp)>=threshold);
    const expectedIds = [...retentionChanged,...differences].map(row=>row.id).sort();
    const result = reviewRows(audit.rows,from,to,threshold);
    assert.deepEqual(result.map(row=>row.row.id).sort(),expectedIds,`${from}->${to}, cutoff${threshold}`);
    assert.equal(result.filter(row=>row.missing).length,retentionChanged.length);
    assert.equal(new Set(result.map(row=>row.row.id)).size,result.length);
    let previous = Infinity, seenMatched = false;
    for (const item of result) {
      assert.ok(item.a || item.b,"absent under both settings must not enter review");
      if (item.missing) {
        assert.equal(seenMatched,false,"retention changes must appear first");
        assert.equal(item.delta,null);
      } else {
        seenMatched = true;
        const direct = item.row[`gap${to}`].tl_bp-item.row[`gap${from}`].tl_bp;
        assert.equal(item.delta,direct);
        assert.ok(Math.abs(item.delta)<=previous,"matched rows ordered by absolute change");
        previous = Math.abs(item.delta);
      }
    }
    reviewCounts.push({from,to,threshold,changed:differences.length,retentionChanged:retentionChanged.length});
    auditCases++;
  }
}
assert.equal(auditCases,12);
assert.equal(reviewRows(audit.rows,100,250,100).filter(row=>row.missing).length,3,
  "100->250 has three retention differences, not the two reads missing from both runs");
assert.deepEqual(reviewRows([{id:"none",gap20:null,gap100:null}],20,100,0),[]);
const zeroLength = reviewRows([{id:"zero",gap20:{tl_bp:0},gap100:{tl_bp:0}}],20,100,0);
assert.equal(zeroLength.length,1);
assert.equal(zeroLength[0].delta,0);
assert.equal(zeroLength[0].missing,false);

const escaped = csv([["a,b",'say "hi"',"two\nlines",null,undefined,0,-1,"=SUM(1,2)","+formula","@formula","-text"]]);
assert.equal(escaped,'"a,b","say ""hi""","two\nlines","","","0","-1","\'=SUM(1,2)","\'+formula","\'@formula","\'-text"\r\n');
assert.equal(csv([[false,true,"",0]]),'"false","true","","0"\r\n');

console.log(JSON.stringify({status:"PASS",primary_cases:cases,python_cross_checks:pythonComparisons,
  audit_cases:auditCases,row_ids:"2610 unique physical row IDs; 25 repeated visible records retained",
  tested_boundaries:["strict short-tail threshold","inclusive read-length filter","zero vs missing",
    "no result for an empty sample","both-gap-absent reads excluded","CSV quotes, commas, newlines and formula prefixes"],
  review_counts:reviewCounts},null,2));
