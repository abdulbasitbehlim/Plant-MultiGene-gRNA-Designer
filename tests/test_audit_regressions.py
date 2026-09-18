"""Regression cases exposing concrete v1.3.1 limitations; all synthetic unless labelled."""
import copy
import json
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from plant_multiguide import (design_shared_guides, evaluate_candidate_guide, scan_gene_segments,
    scan_spcas9, screen_reference_panel_report, suggest_exact_guide_set, MAX_INPUT_BP)
from sequence_sources import (parse_multifasta, manual_records, extract_coding_segments,
                              validate_record_set, fetch_ncbi_gene)
from validation import validate_shared_guide
from run_state import (create_run_snapshot, panel_result_key, export_run, export_input_fasta,
                       cas_offinder_input)
S = 'ACGTACGTACGTACGTACGA'


def mutated(positions):
    seq = list(S)
    for p in positions:
        seq[p - 1] = {'A': 'C', 'C': 'G', 'G': 'T', 'T': 'A'}[seq[p - 1]]
    return ''.join(seq)


def adversarial():
    return {'A': [('e', S+'AGG')], 'B': [('invalid', mutated([1,2,3])+'AGG'),
                                       ('valid', mutated([13,14])+'AGG')]}


def test_feasible_target_not_hidden_by_closer_invalid_site():
    sites, guides = design_shared_guides(adversarial(), max_seed_mismatches_per_gene=2, min_compatibility=0)
    guide = next(g for g in guides if g.spacer == S)
    assert next(m for m in guide.matches if m.gene == 'B').segment_id == 'valid'
    custom = evaluate_candidate_guide(S, sites, 2, 2, 0)
    assert next(m for m in custom.matches if m.gene == 'B').mismatches == 2


def test_minimum_proxy_constraint_applied_before_site_choice():
    # Two middle mismatches have a lower weighted distance (4) than one distal
    # plus one seed (4, tie); feasibility must be checked independently.
    sites = scan_gene_segments({'A':[('e',S+'AGG')], 'B':[('low',mutated([9,10])+'AGG'),
                                                        ('high',mutated([1,2])+'AGG')]})
    custom = evaluate_candidate_guide(S, sites, 2, 0, 85)
    assert custom.minimum_compatibility >= 85


def test_order_independent_rank_and_site_choice():
    genes = {'B':[('z',S+'AGG'),('a',S+'TGG')], 'A':[('x',S+'CGG')]}
    _, a = design_shared_guides(genes)
    _, b = design_shared_guides({k:list(reversed(v)) for k,v in reversed(list(genes.items()))})
    assert a == b


@pytest.mark.parametrize('kwargs', [{'max_results':0},{'max_mismatches_per_gene':-1},
    {'max_seed_mismatches_per_gene':9},{'min_compatibility':float('nan')}])
def test_invalid_parameters_fail(kwargs):
    with pytest.raises(ValueError):
        design_shared_guides(adversarial(), **kwargs)


@pytest.mark.parametrize('raw', ['>A\nACGT\n>A\nTTTT','>A\n\n>B\nAAAA','bad\n>A\nAAAA'])
def test_bad_fasta_rejected(raw):
    with pytest.raises(ValueError): parse_multifasta(raw)


def test_grouped_exons_count_genes_and_do_not_join():
    raw = '>A|e1\n'+'A'*20+'\n>A|e2\nAGG\n>B|e1\n'+S+'AGG'
    recs = manual_records(raw, group_segments=True)
    assert len(recs) == 2 and len(recs['A'].segments) == 2
    assert not scan_gene_segments({g:r.segments for g,r in recs.items()})['A']
    again = manual_records(export_input_fasta(recs), group_segments=True)
    assert [[s for _,s in r.segments] for r in recs.values()] == [[s for _,s in r.segments] for r in again.values()]


def test_grouped_headers_and_labels_validated():
    with pytest.raises(ValueError): manual_records('>A\nAAAA', group_segments=True)
    with pytest.raises(ValueError): manual_records('>A\nAAAA', gene_names=['x','y'])


def test_short_exons_never_fall_back_to_spliced_cds():
    rec = SeqRecord(Seq(S+'AGG'))
    rec.features = [SeqFeature(FeatureLocation(0,23),type='CDS'),
                    SeqFeature(FeatureLocation(0,20),type='exon'),
                    SeqFeature(FeatureLocation(20,23),type='exon')]
    segments, _ = extract_coding_segments(rec)
    assert list(map(lambda p:len(p[1]),segments)) == [20,3]
    assert not scan_gene_segments({'A':segments})['A']


@pytest.mark.parametrize('strand',[1,-1])
def test_compound_cds_omits_introns_and_does_not_join(strand):
    rec = SeqRecord(Seq('A'*20+S+'AGG'+'T'*20))
    rec.features = [SeqFeature(CompoundLocation([FeatureLocation(0,20,strand=strand),
                       FeatureLocation(43,63,strand=strand)]), type='CDS')]
    segments, _ = extract_coding_segments(rec)
    assert len(segments) == 2
    assert not scan_gene_segments({'A':segments})['A']


def test_multiple_cds_rejected():
    rec = SeqRecord(Seq('A'*50))
    rec.features = [SeqFeature(FeatureLocation(0,25),type='CDS'), SeqFeature(FeatureLocation(25,50),type='CDS')]
    with pytest.raises(ValueError,match='Multiple CDS'): extract_coding_segments(rec)


def test_known_duplicate_genes_and_mixed_species_rejected():
    records = manual_records('>A\n'+S+'AGG\n>B\n'+S+'TGG', organism='Arabidopsis')
    records['B'].gene = 'A'
    with pytest.raises(ValueError,match='same gene'): validate_record_set(records)
    records['B'].gene = 'B'; records['B'].organism = 'rice'
    with pytest.raises(ValueError,match='different organisms'): validate_record_set(records)


def test_ambiguous_gene_lookup_rejected(monkeypatch):
    class Response:
        def json(self): return {'esearchresult':{'idlist':['1','2']}}
    monkeypatch.setattr('sequence_sources._requests_get',lambda *a,**k:Response())
    with pytest.raises(ValueError,match='multiple genes'): fetch_ncbi_gene('ambiguous','rice')


@pytest.mark.parametrize('panel',[{}, {'empty':''}])
def test_empty_panel_never_passes(panel):
    with pytest.raises(ValueError,match='nonempty'): screen_reference_panel_report(S,panel)


def test_truncated_panel_retains_total_and_fingerprint():
    panel = {'chr':(S+'AGG'+'N'*23)*7}
    r = screen_reference_panel_report(S,panel,0,2)
    assert r.total_hits == 7 and len(r.hits) == 2 and r.truncated
    assert len(r.panel_sha256)==64 and r.mismatch_radius == 0
    assert r.scanned_sites >= 7


def test_input_and_panel_guards():
    with pytest.raises(ValueError,match='limit'): scan_gene_segments({'A':[('e','A'*(MAX_INPUT_BP+1))]})
    with pytest.raises(ValueError,match='exceeds'): screen_reference_panel_report(S,{'A':'A'*(MAX_INPUT_BP+1)})


def test_exact_fallback_and_uncovered_genes():
    sites = scan_gene_segments({'A':[('e',S+'AGG')], 'B':[('e',mutated([1,2,3,4])+'AGG')]})
    result = suggest_exact_guide_set(sites,2)
    assert result['complete'] and len(result['guides']) == 2
    assert not suggest_exact_guide_set(sites,1)['complete']
    sites['no_pam'] = []
    assert suggest_exact_guide_set(sites,5)['uncovered_genes'] == ['no_pam']


def test_pam_loss_and_reverse_strand():
    assert not scan_spcas9(S+'AGA','A')
    from plant_multiguide import reverse_complement
    sites = scan_spcas9('CCA'+reverse_complement(S),'A')
    m = next(x for x in sites if x.spacer == S)
    assert (m.start,m.end,m.strand,m.pam)==(4,23,'-','TGG')


def test_run_snapshot_and_panel_keys_prevent_stale_evidence():
    recs = manual_records('>A\n'+S+'AGG\n>B\n'+S+'TGG')
    settings={'max_mismatches_per_gene':2}
    old = create_run_snapshot(recs,settings)
    settings['max_mismatches_per_gene']=0
    new = create_run_snapshot(recs,settings)
    assert old['settings']['max_mismatches_per_gene']==2 and old['run_id']!=new['run_id']
    keys=[panel_result_key(old['run_id'],S,'>x\nAAAA',0), panel_result_key(new['run_id'],S,'>x\nAAAA',0),
          panel_result_key(old['run_id'],S,'>x\nAAAT',0),panel_result_key(old['run_id'],S,'>x\nAAAA',1)]
    assert len(set(keys))==4
    payload=export_run(old,[],[])
    assert json.loads(json.dumps(payload))['inputs']['A']['segments'][0][1]==S+'AGG'


def test_validation_recomputes_evidence_and_checks_names():
    _, guides = design_shared_guides({'A':[('e',S+'AGG')], 'B':[('e',S+'AGG')]})
    g=guides[0]
    assert validate_shared_guide(g,2,expected_gene_ids=['A','B']).status=='PASS'
    assert validate_shared_guide(g,2,expected_gene_ids=['A','C']).status=='FAIL'
    wrong=copy.deepcopy(g);wrong.spacer=mutated([1])
    assert validate_shared_guide(wrong,2).status=='FAIL'
    assert validate_shared_guide(g,2,input_warnings=['spliced CDS']).status=='REVIEW'


def test_external_input_format_and_scope():
    lines=cas_offinder_input([S,S],'/genome/plant.fa',3,True).splitlines()
    assert lines==['/genome/plant.fa','N'*20+'NRG',S+'NNN 3']
    assert cas_offinder_input([S],'/genome',0).splitlines()[1]=='N'*20+'NGG'
    with pytest.raises(ValueError): cas_offinder_input([S],'x\ny')


def test_duplicate_resolved_accessions_are_not_distinct_genes():
    recs=manual_records('>A\n'+S+'AGG\n>B\n'+S+'AGG')
    for rec in recs.values():
        rec.source='NCBI RefSeq'; rec.source_record_version='NM_EXAMPLE.1'
    with pytest.raises(ValueError,match='same source record'): validate_record_set(recs)


def test_plant_transcript_suffix_is_not_removed(monkeypatch):
    import accession_sources as ac
    seen=[]
    def fake(url,**kwargs):
        seen.append(url)
        raise RuntimeError('stop after inspecting request')
    monkeypatch.setattr(ac,'_requests_get',fake)
    with pytest.raises(RuntimeError): ac.fetch_ensembl_accession('AT4G18197.1')
    assert seen[0].endswith('/lookup/id/AT4G18197.1')
