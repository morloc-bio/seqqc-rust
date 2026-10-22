// Rust implementations of seqqc.

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct BioReadQc {
    pub reads: i64,
    pub bases: i64,
    pub lengths: Vec<i64>,
    pub posQuals: Vec<i64>,
    pub posBases: Vec<i64>,
    pub gc: Vec<i64>,
    pub meanQuals: Vec<i64>,
    pub adapterStarts: Vec<i64>,
    pub dupSketch: Vec<(i64, i64)>,
    pub dupCapacity: i64,
    pub heavy: Vec<(String, i64)>,
    pub heavyCapacity: i64,
    pub heavyError: i64,
}

#[derive(Clone, Debug, Default)]
pub struct BioPairQc {
    pub first: BioReadQc,
    pub second: BioReadQc,
}

#[allow(non_snake_case)]
#[derive(Clone, Debug, Default)]
pub struct BioQcReport {
    pub reads: i64,
    pub bases: i64,
    pub minLength: i64,
    pub maxLength: i64,
    pub meanLength: f64,
    pub gcPercent: f64,
    pub meanQuality: f64,
    pub positionQuality: Vec<f64>,
    pub positionN: Vec<f64>,
    pub distinct: i64,
    pub dedupPercent: f64,
    pub duplicationLevels: Vec<f64>,
    pub overrepresented: Vec<(String, i64, f64)>,
    pub overrepresentedError: i64,
    pub adapterContent: Vec<Vec<f64>>,
}

const QC_QUAL_BINS: usize = 42;
const QC_MAX_LEN: usize = 1000;
const QC_DUP_CAPACITY: usize = 100_000;
const QC_HEAVY_CAPACITY: usize = 100_000;
const QC_ADAPTERS: [&[u8]; 4] = [b"AGATCGGAAGAG", b"CTGTCTCTTATA", b"TGGAATTCTCGG", b"AAAAAAAAAAAA"];

fn qc_fnv(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in bytes {
        h ^= b as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    h
}

// FastQC's duplication key: reads longer than 75 bases are cut to 50.
fn qc_key(bases: &[u8]) -> &[u8] {
    if bases.len() > 75 { &bases[..50] } else { bases }
}

fn qc_add(v: &mut Vec<i64>, i: usize, n: i64) {
    if v.len() <= i {
        v.resize(i + 1, 0);
    }
    v[i] += n;
}

// Add to slot `col` of row `row` in a row-major table of `width` columns,
// growing it by whole rows.
fn qc_add_row(v: &mut Vec<i64>, width: usize, row: usize, col: usize, n: i64) {
    if v.len() < (row + 1) * width {
        v.resize((row + 1) * width, 0);
    }
    v[row * width + col] += n;
}

fn qc_add_vec(a: &[i64], b: &[i64]) -> Vec<i64> {
    let mut out = vec![0i64; a.len().max(b.len())];
    for (i, x) in a.iter().enumerate() {
        out[i] += x;
    }
    for (i, x) in b.iter().enumerate() {
        out[i] += x;
    }
    out
}

fn qc_empty() -> BioReadQc {
    BioReadQc {
        gc: vec![0; 101],
        meanQuals: vec![0; QC_QUAL_BINS],
        dupCapacity: QC_DUP_CAPACITY as i64,
        heavyCapacity: QC_HEAVY_CAPACITY as i64,
        ..Default::default()
    }
}

// 2-bit codes of the adapters' 12-base prefixes.
fn qc_adapter_codes() -> [u32; 4] {
    let mut out = [0u32; 4];
    for (a, ad) in QC_ADAPTERS.iter().enumerate() {
        for &c in ad.iter() {
            out[a] = (out[a] << 2) | qc_code(c).unwrap();
        }
    }
    out
}

#[inline]
fn qc_code(c: u8) -> Option<u32> {
    match c {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

// Keep the `cap` smallest hashes (as unsigned), ascending.
fn qc_bottom(mut xs: Vec<(i64, i64)>, cap: usize) -> Vec<(i64, i64)> {
    xs.sort_unstable_by_key(|x| x.0 as u64);
    xs.truncate(cap);
    xs
}

// The union of two sketches (each ascending by unsigned hash), counts of a
// shared hash added, cut to the `cap` smallest hashes.
fn qc_merge_sketch(a: &[(i64, i64)], b: &[(i64, i64)], cap: usize) -> Vec<(i64, i64)> {
    let mut out = Vec::with_capacity((a.len() + b.len()).min(cap));
    let (mut i, mut j) = (0, 0);
    while out.len() < cap && (i < a.len() || j < b.len()) {
        let next = if j == b.len() || (i < a.len() && (a[i].0 as u64) < (b[j].0 as u64)) {
            i += 1;
            a[i - 1]
        } else if i == a.len() || (b[j].0 as u64) < (a[i].0 as u64) {
            j += 1;
            b[j - 1]
        } else {
            i += 1;
            j += 1;
            (a[i - 1].0, a[i - 1].1 + b[j - 1].1)
        };
        out.push(next);
    }
    out
}

// Two sets of counters summed, cut to a Misra-Gries summary of `cap`
// counters, by count descending then sequence. Returns the decrement
// applied.
fn qc_merge_heavy(a: &[(String, i64)], b: &[(String, i64)], cap: usize) -> (Vec<(String, i64)>, i64) {
    let mut sums: QcMap<&str, i64> = QcMap::default();
    sums.reserve(a.len() + b.len());
    for (k, c) in a.iter().chain(b) {
        *sums.entry(k.as_str()).or_insert(0) += c;
    }
    let mut xs: Vec<(&str, i64)> = sums.into_iter().collect();
    xs.sort_unstable_by(|x, y| y.1.cmp(&x.1).then_with(|| x.0.cmp(y.0)));
    let mut dec = 0;
    if xs.len() > cap {
        dec = xs[cap].1;
        xs.truncate(cap);
        for x in xs.iter_mut() {
            x.1 -= dec;
        }
        xs.retain(|x| x.1 > 0);
    }
    (xs.into_iter().map(|(k, c)| (k.to_string(), c)).collect(), dec)
}

// A fast hasher for the chunk-local maps (keys are already hashes or short
// sequences; no adversarial input).
#[derive(Default)]
struct QcHasher(u64);
impl std::hash::Hasher for QcHasher {
    fn finish(&self) -> u64 {
        self.0
    }
    fn write(&mut self, bytes: &[u8]) {
        for chunk in bytes.chunks(8) {
            let mut w = [0u8; 8];
            w[..chunk.len()].copy_from_slice(chunk);
            self.0 = (self.0.rotate_left(5) ^ u64::from_le_bytes(w)).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
        }
    }
    fn write_u64(&mut self, x: u64) {
        self.0 = (self.0.rotate_left(5) ^ x).wrapping_mul(0x51_7c_c1_b7_27_22_0a_95);
    }
}
type QcMap<K, V> = std::collections::HashMap<K, V, std::hash::BuildHasherDefault<QcHasher>>;

// Base code: A C G T as 0..3, anything else 4.
const QC_BASE: [u8; 256] = {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t
};

fn qc_grow(v: &mut Vec<i64>, len: usize) {
    if v.len() < len {
        v.resize(len, 0);
    }
}

fn qc_chunk<'a, I: Iterator<Item = &'a BioFastqEntry>>(reads: I) -> BioReadQc {
    let mut q = qc_empty();
    let codes = qc_adapter_codes();
    let na = QC_ADAPTERS.len();
    let mut dup: QcMap<u64, i64> = QcMap::default();
    let mut heavy: QcMap<Vec<u8>, i64> = QcMap::default();
    let mut key = Vec::with_capacity(64);
    for e in reads {
        let (b, s) = (&e.bases.bases, &e.bases.quals);
        let n = b.len();
        let m = n.min(QC_MAX_LEN);
        q.reads += 1;
        q.bases += n as i64;
        qc_add(&mut q.lengths, m, 1);
        qc_grow(&mut q.posQuals, m * QC_QUAL_BINS);
        qc_grow(&mut q.posBases, m * 5);
        let (mut gc, mut acgt, mut qsum) = (0usize, 0usize, 0i64);
        let (mut code, mut valid) = (0u32, 0usize);
        let mut seen = [false; 4];
        for i in 0..m {
            let ph = (s[i] as i64 - 33).clamp(0, QC_QUAL_BINS as i64 - 1);
            qsum += ph;
            q.posQuals[i * QC_QUAL_BINS + ph as usize] += 1;
            let bi = QC_BASE[b[i] as usize] as usize;
            q.posBases[i * 5 + bi] += 1;
            if bi < 4 {
                acgt += 1;
                gc += (bi == 1 || bi == 2) as usize;
                code = ((code << 2) | bi as u32) & 0xff_ffff;
                valid += 1;
                if valid >= 12 {
                    let p = i + 1 - 12;
                    for a in 0..na {
                        if !seen[a] && code == codes[a] {
                            seen[a] = true;
                            qc_add_row(&mut q.adapterStarts, na, p, a, 1);
                        }
                    }
                }
            } else {
                valid = 0;
            }
        }
        // An adapter may start before the last tabulated position and end
        // after it.
        for i in m..n.min(QC_MAX_LEN + 11) {
            let bi = QC_BASE[b[i] as usize] as u32;
            if bi < 4 {
                code = ((code << 2) | bi) & 0xff_ffff;
                valid += 1;
                if valid >= 12 && i + 1 - 12 < QC_MAX_LEN {
                    for a in 0..na {
                        if !seen[a] && code == codes[a] {
                            seen[a] = true;
                            qc_add_row(&mut q.adapterStarts, na, i + 1 - 12, a, 1);
                        }
                    }
                }
            } else {
                valid = 0;
            }
        }
        if acgt > 0 {
            q.gc[((100 * gc) as f64 / acgt as f64).round() as usize] += 1;
        }
        if m > 0 {
            q.meanQuals[((qsum / m as i64) as usize).min(QC_QUAL_BINS - 1)] += 1;
        }
        key.clear();
        key.extend(qc_key(b).iter().map(|c| c.to_ascii_uppercase()));
        *dup.entry(qc_fnv(&key)).or_insert(0) += 1;
        match heavy.get_mut(key.as_slice()) {
            Some(c) => *c += 1,
            None => {
                heavy.insert(key.clone(), 1);
            }
        }
    }
    q.dupSketch = qc_bottom(dup.into_iter().map(|(h, c)| (h as i64, c)).collect(), QC_DUP_CAPACITY);
    let mut hv: Vec<(String, i64)> = Vec::new();
    let mut top: Vec<(Vec<u8>, i64)> = heavy.into_iter().collect();
    if top.len() > QC_HEAVY_CAPACITY {
        top.select_nth_unstable_by(QC_HEAVY_CAPACITY, |a, b| b.1.cmp(&a.1));
    }
    let mut dec = 0;
    if top.len() > QC_HEAVY_CAPACITY {
        dec = top[QC_HEAVY_CAPACITY].1;
        top.truncate(QC_HEAVY_CAPACITY);
    }
    let mut lossy = false;
    for (k, c) in top {
        if c - dec > 0 {
            let k = String::from_utf8(k).unwrap_or_else(|e| {
                lossy = true;
                String::from_utf8_lossy(e.as_bytes()).into_owned()
            });
            hv.push((k, c - dec));
        }
    }
    // keys made equal by the lossy conversion are summed
    if lossy {
        hv = qc_merge_heavy(&hv, &[], usize::MAX).0;
    } else {
        hv.sort_unstable_by(|x, y| y.1.cmp(&x.1).then_with(|| x.0.cmp(&y.0)));
    }
    q.heavyError = dec;
    q.heavy = hv;
    q
}

fn qc_merge(a: &BioReadQc, b: &BioReadQc) -> BioReadQc {
    let cap = a.dupCapacity.max(b.dupCapacity);
    let hcap = a.heavyCapacity.max(b.heavyCapacity);
    let (heavy, dec) = qc_merge_heavy(&a.heavy, &b.heavy, hcap as usize);
    BioReadQc {
        reads: a.reads + b.reads,
        bases: a.bases + b.bases,
        lengths: qc_add_vec(&a.lengths, &b.lengths),
        posQuals: qc_add_vec(&a.posQuals, &b.posQuals),
        posBases: qc_add_vec(&a.posBases, &b.posBases),
        gc: qc_add_vec(&a.gc, &b.gc),
        meanQuals: qc_add_vec(&a.meanQuals, &b.meanQuals),
        adapterStarts: qc_add_vec(&a.adapterStarts, &b.adapterStarts),
        dupSketch: qc_merge_sketch(&a.dupSketch, &b.dupSketch, cap as usize),
        dupCapacity: cap,
        heavy,
        heavyCapacity: hcap,
        heavyError: a.heavyError + b.heavyError + dec,
    }
}

pub fn bioseqqc_merge_read_qc(a: &BioReadQc, b: &BioReadQc) -> BioReadQc {
    qc_merge(a, b)
}

pub fn bioseqqc_read_qc(xs: &Vec<BioFastqEntry>) -> BioReadQc {
    qc_chunk(xs.iter())
}

pub fn bioseqqc_pair_qc(xs: &Vec<(BioFastqEntry, BioFastqEntry)>) -> BioPairQc {
    BioPairQc { first: qc_chunk(xs.iter().map(|p| &p.0)), second: qc_chunk(xs.iter().map(|p| &p.1)) }
}

pub fn bioseqqc_qc_report(q: &BioReadQc) -> BioQcReport {
    let reads = q.reads.max(0);
    let rf = reads.max(1) as f64;
    let mut r = BioQcReport { reads, bases: q.bases, ..Default::default() };
    let present: Vec<usize> = q.lengths.iter().enumerate().filter(|(_, &c)| c > 0).map(|(i, _)| i).collect();
    r.minLength = present.first().copied().unwrap_or(0) as i64;
    r.maxLength = present.last().copied().unwrap_or(0) as i64;
    r.meanLength = q.bases as f64 / rf;
    let gcn: i64 = q.gc.iter().sum();
    r.gcPercent = if gcn > 0 { q.gc.iter().enumerate().map(|(i, &c)| i as f64 * c as f64).sum::<f64>() / gcn as f64 } else { 0.0 };
    let npos = q.posBases.len() / 5;
    let (mut qtot, mut ntot) = (0f64, 0f64);
    for i in 0..npos {
        let row = &q.posQuals[i * QC_QUAL_BINS..(i + 1) * QC_QUAL_BINS];
        let n: i64 = row.iter().sum();
        let s: i64 = row.iter().enumerate().map(|(p, &c)| p as i64 * c).sum();
        qtot += s as f64;
        ntot += n as f64;
        r.positionQuality.push(if n > 0 { s as f64 / n as f64 } else { 0.0 });
        let bases = &q.posBases[i * 5..i * 5 + 5];
        let cover: i64 = bases.iter().sum();
        r.positionN.push(if cover > 0 { 100.0 * bases[4] as f64 / cover as f64 } else { 0.0 });
    }
    r.meanQuality = if ntot > 0.0 { qtot / ntot } else { 0.0 };
    // Distinct sequences: exact when the sketch is not full, else the KMV
    // estimate from the largest kept hash.
    let k = q.dupSketch.len();
    r.distinct = if (k as i64) < q.dupCapacity || k < 2 {
        k as i64
    } else {
        let hk = (q.dupSketch[k - 1].0 as u64) as f64 / 18446744073709551616.0;
        ((k as f64 - 1.0) / hk).round() as i64
    };
    r.dedupPercent = if reads > 0 { 100.0 * (r.distinct.min(reads)) as f64 / rf } else { 0.0 };
    let mut levels = vec![0f64; 10];
    for &(_, c) in &q.dupSketch {
        levels[(c.clamp(1, 10) - 1) as usize] += 1.0;
    }
    r.duplicationLevels = levels.iter().map(|x| if k > 0 { 100.0 * x / k as f64 } else { 0.0 }).collect();
    let mut over: Vec<(String, i64, f64)> = q
        .heavy
        .iter()
        .filter(|(_, c)| (*c + q.heavyError) as f64 > 0.001 * rf)
        .map(|(s, c)| (s.clone(), *c, 100.0 * *c as f64 / rf))
        .collect();
    over.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
    r.overrepresented = over;
    r.overrepresentedError = q.heavyError;
    let na = QC_ADAPTERS.len();
    let apos = q.adapterStarts.len() / na;
    let width = apos.max(npos);
    for a in 0..na {
        let mut acc = 0i64;
        let mut row = Vec::with_capacity(width);
        for i in 0..width {
            acc += q.adapterStarts.get(i * na + a).copied().unwrap_or(0);
            row.push(100.0 * acc as f64 / rf);
        }
        r.adapterContent.push(row);
    }
    r
}

// FastQC 0.12.1 report text
// -------------------------

fn qc_flag(warn: bool, fail: bool) -> &'static str {
    if fail { "fail" } else if warn { "warn" } else { "pass" }
}

// The quantile `q` of a histogram of phred bins.
fn qc_quantile(row: &[i64], q: f64) -> f64 {
    let n: i64 = row.iter().sum();
    if n == 0 {
        return 0.0;
    }
    let target = (q * n as f64).ceil().max(1.0) as i64;
    let mut acc = 0;
    for (i, &c) in row.iter().enumerate() {
        acc += c;
        if acc >= target {
            return i as f64;
        }
    }
    (row.len() - 1) as f64
}

fn qc_bp(n: i64) -> String {
    if n >= 1_000_000_000 { format!("{:.1} Gbp", n as f64 / 1e9) }
    else if n >= 1_000_000 { format!("{:.1} Mbp", n as f64 / 1e6) }
    else if n >= 1_000 { format!("{:.1} kbp", n as f64 / 1e3) }
    else { format!("{} bp", n) }
}

// FastQC's fastqc_data.txt for a summary, with FastQC's default pass, warn
// and fail limits, positions ungrouped (FastQC's --nogroup).
pub fn bioseqqc_fastqc_data(name: &String, q: &BioReadQc) -> String {
    use std::fmt::Write;
    let r = bioseqqc_qc_report(q);
    let npos = q.posBases.len() / 5;
    let mut s = String::from("##FastQC\t0.12.1\n");
    let lens: Vec<usize> = q.lengths.iter().enumerate().filter(|(_, &c)| c > 0).map(|(i, _)| i).collect();
    let len_text = match (lens.first(), lens.last()) {
        (Some(a), Some(b)) if a != b => format!("{}-{}", a, b),
        (Some(a), _) => a.to_string(),
        _ => "0".to_string(),
    };
    // FastQC's %GC: GC bases over A, C, G and T bases, truncated.
    let (mut gc, mut acgt) = (0i64, 0i64);
    for b in q.posBases.chunks(5) {
        gc += b[1] + b[2];
        acgt += b[0] + b[1] + b[2] + b[3];
    }
    let gc_bases = if acgt > 0 { gc * 100 / acgt } else { 0 };
    let _ = write!(
        s,
        ">>Basic Statistics\tpass\n#Measure\tValue\nFilename\t{}\nFile type\tConventional base calls\nEncoding\tSanger / Illumina 1.9\nTotal Sequences\t{}\nTotal Bases\t{}\nSequences flagged as poor quality\t0\nSequence length\t{}\n%GC\t{}\n>>END_MODULE\n",
        name, q.reads, qc_bp(q.bases), len_text, gc_bases
    );
    // per base quality
    let mut rows = String::new();
    let (mut warn, mut fail) = (false, false);
    for i in 0..npos {
        let row = &q.posQuals[i * QC_QUAL_BINS..(i + 1) * QC_QUAL_BINS];
        let (lq, med) = (qc_quantile(row, 0.25), qc_quantile(row, 0.5));
        warn |= lq < 10.0 || med < 25.0;
        fail |= lq < 5.0 || med < 20.0;
        let _ = writeln!(rows, "{}\t{:.6}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{:.1}", i + 1, r.positionQuality[i], med, lq, qc_quantile(row, 0.75), qc_quantile(row, 0.1), qc_quantile(row, 0.9));
    }
    let _ = write!(s, ">>Per base sequence quality\t{}\n#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile\n{}>>END_MODULE\n", qc_flag(warn, fail), rows);
    // per sequence quality
    let mode = q.meanQuals.iter().enumerate().max_by_key(|(_, &c)| c).map(|(i, _)| i).unwrap_or(0);
    let mut rows = String::new();
    for (i, &c) in q.meanQuals.iter().enumerate() {
        if c > 0 {
            let _ = writeln!(rows, "{}\t{:.1}", i, c as f64);
        }
    }
    let _ = write!(s, ">>Per sequence quality scores\t{}\n#Quality\tCount\n{}>>END_MODULE\n", qc_flag(mode < 27, mode < 20), rows);
    // per base content
    let mut rows = String::new();
    let mut worst = 0f64;
    for i in 0..npos {
        let b = &q.posBases[i * 5..i * 5 + 5];
        let acgt = (b[0] + b[1] + b[2] + b[3]).max(1) as f64;
        let pct = |k: usize| 100.0 * b[k] as f64 / acgt;
        worst = worst.max((pct(0) - pct(3)).abs()).max((pct(1) - pct(2)).abs());
        let _ = writeln!(rows, "{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}", i + 1, pct(2), pct(0), pct(3), pct(1));
    }
    let _ = write!(s, ">>Per base sequence content\t{}\n#Base\tG\tA\tT\tC\n{}>>END_MODULE\n", qc_flag(worst > 10.0, worst > 20.0), rows);
    // per sequence GC, against a normal with the observed mean and sd
    let total: i64 = q.gc.iter().sum();
    let mean = r.gcPercent;
    let var = if total > 0 { q.gc.iter().enumerate().map(|(i, &c)| c as f64 * (i as f64 - mean).powi(2)).sum::<f64>() / total as f64 } else { 0.0 };
    let sd = var.sqrt().max(1e-9);
    let dens: Vec<f64> = (0..=100).map(|i| (-0.5 * ((i as f64 - mean) / sd).powi(2)).exp()).collect();
    let dsum: f64 = dens.iter().sum::<f64>().max(1e-300);
    let dev: f64 = q.gc.iter().enumerate().map(|(i, &c)| (c as f64 - total as f64 * dens[i] / dsum).abs()).sum::<f64>();
    let dev_pct = if total > 0 { 100.0 * dev / total as f64 } else { 0.0 };
    let mut rows = String::new();
    for (i, &c) in q.gc.iter().enumerate() {
        let _ = writeln!(rows, "{}\t{:.1}", i, c as f64);
    }
    let _ = write!(s, ">>Per sequence GC content\t{}\n#GC Content\tCount\n{}>>END_MODULE\n", qc_flag(dev_pct > 15.0, dev_pct > 30.0), rows);
    // per base N
    let mut rows = String::new();
    let nmax = r.positionN.iter().cloned().fold(0.0, f64::max);
    for (i, n) in r.positionN.iter().enumerate() {
        let _ = writeln!(rows, "{}\t{:.6}", i + 1, n);
    }
    let _ = write!(s, ">>Per base N content\t{}\n#Base\tN-Count\n{}>>END_MODULE\n", qc_flag(nmax > 5.0, nmax > 20.0), rows);
    // lengths
    let mut rows = String::new();
    for &l in &lens {
        let _ = writeln!(rows, "{}\t{:.1}", l, q.lengths[l] as f64);
    }
    let _ = write!(s, ">>Sequence Length Distribution\t{}\n#Length\tCount\n{}>>END_MODULE\n", qc_flag(lens.len() > 1, lens.first() == Some(&0)), rows);
    // duplication: FastQC's levels, from the sketch
    let labels = ["1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k", ">10k+"];
    let level = |c: i64| -> usize {
        match c {
            1..=9 => (c - 1) as usize,
            10..=49 => 9,
            50..=99 => 10,
            100..=499 => 11,
            500..=999 => 12,
            1000..=4999 => 13,
            5000..=9999 => 14,
            _ => 15,
        }
    };
    let (mut distinct, mut total_reads) = (vec![0f64; 16], vec![0f64; 16]);
    for &(_, c) in &q.dupSketch {
        let l = level(c.max(1));
        distinct[l] += 1.0;
        total_reads[l] += c as f64;
    }
    let (dsum, tsum) = (distinct.iter().sum::<f64>().max(1.0), total_reads.iter().sum::<f64>().max(1.0));
    let mut rows = String::new();
    for i in 0..16 {
        let _ = writeln!(rows, "{}\t{:.6}\t{:.6}", labels[i], 100.0 * distinct[i] / dsum, 100.0 * total_reads[i] / tsum);
    }
    let dedup = r.dedupPercent;
    let _ = write!(s, ">>Sequence Duplication Levels\t{}\n#Total Deduplicated Percentage\t{:.6}\n#Duplication Level\tPercentage of deduplicated\tPercentage of total\n{}>>END_MODULE\n", qc_flag(dedup < 80.0, dedup < 50.0), dedup, rows);
    // overrepresented
    let mut rows = String::new();
    let top = r.overrepresented.iter().map(|x| x.2).fold(0.0, f64::max);
    for (seq, c, pct) in &r.overrepresented {
        let _ = writeln!(rows, "{}\t{}\t{:.6}\tNo Hit", seq, c, pct);
    }
    let over_flag = if r.overrepresented.is_empty() { "pass" } else { qc_flag(top > 0.1, top > 1.0) };
    let _ = write!(s, ">>Overrepresented sequences\t{}\n", over_flag);
    if !r.overrepresented.is_empty() {
        let _ = write!(s, "#Sequence\tCount\tPercentage\tPossible Source\n{}", rows);
    }
    s.push_str(">>END_MODULE\n");
    // adapter content
    let names = ["Illumina Universal Adapter", "Nextera Transposase Sequence", "Illumina Small RNA 3' Adapter", "PolyA"];
    let width = r.adapterContent.first().map(|v| v.len()).unwrap_or(0).min(npos.max(1));
    let amax = r.adapterContent.iter().flat_map(|v| v.iter().cloned()).fold(0.0, f64::max);
    let mut rows = String::new();
    for i in 0..width {
        let vals: Vec<String> = r.adapterContent.iter().map(|v| format!("{:.6}", v[i])).collect();
        let _ = writeln!(rows, "{}\t{}", i + 1, vals.join("\t"));
    }
    let _ = write!(s, ">>Adapter Content\t{}\n#Position\t{}\n{}>>END_MODULE\n", qc_flag(amax > 5.0, amax > 10.0), names.join("\t"), rows);
    s
}
