import os
import polars as pl
import tarfile
import io
from collections import defaultdict
from Bio.Data import CodonTable

def load_cocoputs(filepath=None):
    if not filepath or not os.path.exists(filepath): return {}
    try:
        if filepath.endswith('.tar.xz'):
            with tarfile.open(filepath, "r:xz") as tar:
                member = tar.next()
                f = tar.extractfile(member)
                df = pl.read_csv(f.read(), separator='\t', infer_schema_length=10000, ignore_errors=True, truncate_ragged_lines=True)
        else:
            df = pl.read_csv(filepath, separator='\t', infer_schema_length=10000, ignore_errors=True, truncate_ragged_lines=True)
        cols = df.columns
        name_col = 'Species' if 'Species' in cols else ('Tissue' if 'Tissue' in cols else None)
        if not name_col: return {}
        codons = [c for c in cols if len(c) == 3 and all(b in 'ACGTU' for b in c)]
        if 'Organelle' in cols:
            df = df.with_columns(pl.when(pl.col('Organelle').is_not_null() & (pl.col('Organelle') != 'genomic') & (pl.col('Organelle') != 'NA')).then(pl.format("{} [{}]", pl.col(name_col), pl.col('Organelle'))).otherwise(pl.col(name_col)).alias('__display_name__'))
        else:
            df = df.with_columns(pl.col(name_col).alias('__display_name__'))
        aggs = [pl.col(c).sum() for c in codons]
        if 'Translation Table' in cols:
            df = df.with_columns(pl.col('Translation Table').cast(pl.Int64, strict=False).fill_null(1))
            aggs.append(pl.col('Translation Table').first())
        df = df.group_by('__display_name__').agg(aggs)
        data = {}
        for row in df.to_dicts():
            species = row['__display_name__']
            counts = {c: float(row[c]) for c in codons if row[c] is not None}
            table_id = int(row['Translation Table']) if 'Translation Table' in row and row['Translation Table'] is not None else 1
            data[species] = {'counts': counts, 'table_id': 1 if table_id == 0 else table_id}
        return data
    except: return {}

def normalize_to_frequency(table):
    if not table: return {}
    normalized = {}
    for org, org_data in table.items():
        counts, table_id = org_data['counts'], org_data['table_id']
        try: standard_table = CodonTable.unambiguous_dna_by_id[table_id].forward_table
        except: standard_table = CodonTable.unambiguous_dna_by_id[1].forward_table
        normalized[org] = {'frequencies': {}, 'table_id': table_id}
        aa_groups = defaultdict(dict)
        for codon, count in counts.items():
            aa = standard_table.get(codon.replace('U', 'T'), '*')
            aa_groups[aa][codon] = count
        for aa, codon_counts in aa_groups.items():
            total = sum(codon_counts.values())
            for codon, count in codon_counts.items():
                normalized[org]['frequencies'][codon] = count / total if total > 0 else 0
    return normalized

def get_organism_list(parsed_table):
    return sorted(list(parsed_table.keys()))
