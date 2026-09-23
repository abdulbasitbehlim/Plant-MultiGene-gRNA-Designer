# ============================================================================
# RUN AUDIT EXAMPLES
# BEGINNER-FRIENDLY CODE GUIDE
# ============================================================================
#
# PURPOSE: Runs reproducible audit examples used to inspect important workflow decisions.
#
# HOW TO READ THIS FILE:
# 1. Read the imports/constants first to see which tools and settings are used.
# 2. Read one top-level function or class at a time.
# 3. Follow the workflow from input sequence -> candidate guides -> validation -> output.
# 4. Scientific formulas, thresholds, validation decisions and public function names
#    are intentionally preserved while readability comments are added.
#
# MAIN TOP-LEVEL PARTS:
# - function: mutate
# - function: record
# ============================================================================

"""Small deterministic synthetic demonstrations; no plant editing efficacy claims."""
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from plant_multiguide import *
from sequence_sources import manual_records,parse_multifasta,extract_coding_segments
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
S='ACGTACGTACGTACGTACGA'


# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: mutate
# ----------------------------------------------------------------------------
def mutate(ps):
 a=list(S)
 for p in ps:a[p-1]={'A':'C','C':'G','G':'T','T':'A'}[a[p-1]]
 return ''.join(a)

rows=[]

# ----------------------------------------------------------------------------
# FUNCTION / CLASS SECTION: record
# ----------------------------------------------------------------------------
def record(name,observed,expected):
 assert observed==expected,(name,observed,expected)
 rows.append({'example':name,'observed':observed,'expected':expected,'passed':True,'evidence':'synthetic software check'})

exact={f'G{i}':[('e',S+p)] for i,p in enumerate(['AGG','TGG','CGG'],1)}
_,g=design_shared_guides(exact,include_mismatch_aware=False)
record('Exact guide across three input genes',max(x.genes_covered for x in g),3)
_,g=design_shared_guides({'G1':[('e',S+'AGG')],'G2':[('e',mutate([5])+'AGG')]})
record('One distal mismatch retained',any(x.spacer==S for x in g),True)
_,g=design_shared_guides({'G1':[('e',S+'AGG')],'G2':[('e',mutate([13])+'AGG')]},max_seed_mismatches_per_gene=0)
record('Zero seed-mismatch limit enforced',any(x.spacer==S for x in g),False)
record('Identical spacer without NGG has no site',len(scan_spcas9(S+'AGA','G')),0)
record('Ambiguous spacer excluded',len(scan_spcas9(S[:5]+'N'+S[6:]+'AGG','G')),0)
record('Reverse-strand site coordinates',[(x.start,x.end) for x in scan_spcas9('CCA'+reverse_complement(S),'G') if x.spacer==S],[(4,23)])
try:parse_multifasta('>A\nAAAA\n>A\nCCCC'); rejected=False
except ValueError:rejected=True
record('Duplicate FASTA identifier rejected',rejected,True)
r=manual_records('>A|e1\n'+S+'\n>A|e2\nAGG\n>B|e1\n'+S+'AGG',group_segments=True)
record('Separate exons do not create junction site',len(scan_gene_segments({'A':r['A'].segments})['A']),0)
rec=SeqRecord(Seq(S+'AGG'));rec.features=[SeqFeature(FeatureLocation(0,23),type='CDS'),SeqFeature(FeatureLocation(0,20),type='exon'),SeqFeature(FeatureLocation(20,23),type='exon')]
parts,_=extract_coding_segments(rec)
record('Short annotated exons stay separate',[len(s) for _,s in parts],[20,3])
a={'A':[('e',S+'AGG')],'B':[('invalid',mutate([1,2,3])+'AGG'),('valid',mutate([13,14])+'AGG')]}
sites,g=design_shared_guides(a,max_seed_mismatches_per_gene=2,min_compatibility=0)
record('Feasible alternative survives invalid closer site',any(x.spacer==S for x in g),True)
panel=screen_reference_panel_report(S,{'p':(S+'AGG'+'N'*23)*7},0,2)
record('Panel total retained despite display cap',[panel.total_hits,len(panel.hits),panel.truncated],[7,2,True])
try:screen_reference_panel_report(S,{});rejected=False
except ValueError:rejected=True
record('Empty reference panel rejected',rejected,True)
sites,_=design_shared_guides({'A':[('e',S+'AGG')],'B':[('e',mutate([1,2,3,4])+'AGG')]},include_mismatch_aware=False)
f=suggest_exact_guide_set(sites,2)
record('Two-guide exact fallback covers unrelated synthetic targets',[len(f['guides']),f['complete']],[2,True])
root=Path(__file__).resolve().parents[1]
(root/'benchmarks/results/audit_examples.json').write_text(json.dumps(rows,indent=2)+'\n')
fixtures={'exact_three_genes':'>G1\n'+S+'AGG\n>G2\n'+S+'TGG\n>G3\n'+S+'CGG\n',
 'mismatch_pair':'>G1\n'+S+'AGG\n>G2\n'+mutate([5])+'AGG\n',
 'separate_exons':'>G1|e1\n'+S+'\n>G1|e2\nAGG\n>G2|e1\n'+S+'AGG\n',
 'no_shared_guide':'>G1\n'+S+'AGG\n>G2\n'+mutate([1,2,3,4])+'AGG\n',
 'pam_loss':'>G1\n'+S+'AGG\n>G2\n'+S+'AGA\n'}
for name,raw in fixtures.items():(root/'examples'/f'{name}.fasta').write_text(raw)
print(f'{len(rows)} synthetic examples passed. Results: benchmarks/results/audit_examples.json')
