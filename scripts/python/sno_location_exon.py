#!/usr/bin/python3
import pandas as pd
import subprocess as sp
from pybedtools import BedTool

""" Locate the host gene intron in which snoRNAs are located and return the
    intron number and distance to exons related to each snoRNA. We take the
    transcript with the highest number of exons per HG as the 'main transcript'.
    We process SNHG14 snoRNAs separately since their HG transcript with the
    highest number of exons doesn't include all of its embedded snoRNAs
    (instead we take the longest SNHG14 transcript from RefSeq so that it
    includes all SNHG14 snoRNAs)"""



col = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature', 'dot2', 'gene_info']
gtf_bed = pd.read_csv(snakemake.input.sorted_gtf_bed, sep='\t', names=col)  # generated with gtf_to_bed

sno_HG_coordinates = pd.read_csv(snakemake.input.sno_HG_coordinates)
sno_info = pd.read_csv(snakemake.input.sno_tpm_df)
sno_info = sno_info[~sno_info['host_id'].isna()]  # remove intergenic snoRNAs from this analysis
sno_info_wo_snhg14_snornas = sno_info[sno_info['host_name'] != 'SNHG14'] #Exclusion of SNHG14 HG because the longest transcript does not cover all of its embedded snoRNAs

sno_bed_wo_snhg14_path = snakemake.input.sno_bed_wo_snhg14  # generated with generate_snoRNA_beds.py
snhg14_bed_path = snakemake.input.snhg14_bed  # Obtained from refseq longest transcript so that it includes all SNHG14 snoRNAs
snhg14_sno_bed_path = snakemake.input.sno_snhg14_bed  # generated with generate_snoRNA_beds.py



def generate_hg_bed(gtf_bed, sno_info, output_path_bed, output_path_bed_col_split):
    """Iterate through a gtf file in a bed format to retrieve only the information related to host genes (hgs) and
        create a resulting hg bed file. sno_info is a df that gives the HG of interest (ex: all HGs vs all HGs except
        SNHG14)"""
    hgs = list(pd.unique(sno_info['host_id']))
    df = []
    for i, hg in enumerate(hgs):
        print(hg)
        temp_df = gtf_bed[gtf_bed['gene_info'].str.contains('"'+hg+'"')]
        df.append(temp_df)

    df_final = pd.concat(df)
    df_final.to_csv(output_path_bed, sep='\t', header=False, index=False)  # This is the bed (from gtf) of all host genes except SNHG14

    sp.call("""awk -i inplace -v OFS='\t' '$8=="exon"' """ + output_path_bed, shell=True)  # Keep only the exon features in the HG bed
    sp.call("""sed -i -E 's/ tag .*;"/"/g; s/ exon_id .*;./"/g' """ + output_path_bed, shell=True)  # remove exon_id and tag infos in gene_info column to simplify

    df_final = pd.read_csv(output_path_bed, sep='\t', names=['chr', 'start', 'end', 'gene_id', 'dot', 'strand',
                                                             'source', 'feature', 'dot2', 'gene_info'])

    # Split the gtf bed file into more readable columns
    hg_bed = df_final
    hg_bed[['empty1', 'gene_id2', 'empty2', 'gene_version', 'empty3', 'transcript_id', 'empty4', 'transcript_version',
            'empty5', 'exon_number', 'empty6', 'gene_name', 'empty7', 'gene_source', 'empty8', 'gene_biotype', 'empty9',
            'transcript_name', 'empty10', 'transcript_source',
            'empty11', 'transcript_biotype']] = hg_bed.gene_info.str.split(" ", expand=True).applymap(
        lambda x: x.replace('"', '')).applymap(lambda x: x.replace(';', ''))

    hg_bed = hg_bed.drop(axis=1,
                         labels=['gene_info', 'empty1', 'empty2', 'empty3', 'empty4', 'empty5', 'empty6', 'gene_source',
                                 'empty7', 'empty8', 'empty9', 'transcript_source', 'empty10', 'empty11',
                                 'gene_version', 'transcript_version'])

    hg_bed.to_csv(output_path_bed_col_split, sep='\t', header=False, index=False)

    print('Finished generate_hg_bed!')
    return hg_bed


def get_max_exon_transcript_per_hg(sno_HG_coordinates, hg_bed, output_path, sno_overlap_path):
    """ Sort HG transcripts per exon_number (in descending order) and by their
        name*** if multiple transcripts have the same number of exons (in ascending
        order since a ...-201 transcript is more present than a ...-205). Then
        iterate trough this ordered groupby object and get the first HG transcript
        that doesn't overlap with the snoRNA (if possible, but there are multiple
        snoRNA exceptions that overlap in all kinds of form with their HG; see all
        elifs below). Then regroup all these HG transcripts (1 per HG) in a bed file.
        ***Transcript name with '-201' are the most present transcript and then less
        and less present with the 201 increasing"""
    most_exon_transcripts = []
    groups = hg_bed.groupby('gene_id')

    # Drop SHNG14 (ENSG00000224078) snoRNAs
    sno_HG_coordinates = sno_HG_coordinates[sno_HG_coordinates['host_id'] != 'ENSG00000224078']
    sno_HG_coordinates = sno_HG_coordinates.set_index('sno_id')
    sno_dict = sno_HG_coordinates.to_dict('index')

    for sno_id, cols in sno_dict.items():
        sno_start, sno_end, host_id = sno_dict[sno_id]['sno_start'], sno_dict[sno_id]['sno_end'], sno_dict[sno_id]['host_id']
        host = groups.get_group(host_id)
        host_transcripts = host.groupby('transcript_id')
        temp = []
        for transcript_id, transcript in host_transcripts:  # Create a exon_max column to sort afterwards
            exon_number = max(map(int, transcript['exon_number']))
            transcript.loc[:, 'exon_max'] = exon_number
            temp.append(transcript)
        host_df = pd.concat(temp)
        trans_max = host_df.sort_values(by=['exon_max', 'transcript_name'], ascending=[False, True])  # Sort by max number of exon to least number of exon, then by the transcript name in ascending order (transcript-201 to ... -213 ...)
        for trans_id, trans_df in trans_max.groupby('transcript_id', sort=False):
            trans_df = trans_df.reset_index()
            l = ['ENSG00000261709', 'ENSG00000277947', 'NR_132981', 'NR_132980']
            if sno_id in l:  # This is to patch for the 4 snoRNAs that are encoded within the exon of a 1-exon lncRNA
                print(trans_id)
                gene = host[host['transcript_id'] == trans_id]
                gene.loc[:, 'sno'] = sno_id
                gene.loc[:, 'hg_overlap'] = 'sno_into_single_exon'
                most_exon_transcripts.append(gene)
                break
            for i in trans_df.index.values[:-1]:  # For all other snoRNAs
                if trans_df.loc[i, 'exon_max'] > 1:
                    exon1 = trans_df.iloc[i, :]
                    exon2 = trans_df.iloc[i+1, :]
                    if (sno_start > exon1['end']) & (sno_end < exon2['start']):  # Normal snoRNAs encoded between two exons
                        best_transcript_id = exon1['transcript_id']
                        gene = host[host['transcript_id'] == best_transcript_id]
                        gene.loc[:, 'sno'] = sno_id
                        gene.loc[:, 'hg_overlap'] = ''
                        most_exon_transcripts.append(gene)
                        break
                    elif (((sno_end < exon1['end']) & (sno_start > exon1['start'])) | ((sno_end < exon2['end']) & (sno_start > exon2['start']))):  # snoRNA overlaps with exon but doesn't extend before or after the exon
                        if sno_id in ['ENSG00000201672', 'ENSG00000206620', 'NR_145790', 'ENSG00000212293']:  # To patch for the 4 snoRNAs that overlap with the HG transcript with the most exon but not with the second most exon HG transcript
                            continue
                        elif sno_id not in ['ENSG00000275662', 'ENSG00000201672', 'ENSG00000206620', 'NR_145790', 'ENSG00000212293']:
                            best_transcript_id = exon1['transcript_id']
                            gene = host[host['transcript_id'] == best_transcript_id]
                            gene.loc[:, 'sno'] = sno_id
                            gene.loc[:, 'hg_overlap'] = 'sno_into_multi_exon'
                            most_exon_transcripts.append(gene)
                            break
                    elif (((sno_start < exon1['end']) & (sno_end > exon1['end'])) | ((sno_start < exon2['end']) & (sno_end > exon2['end']))): # snoRNA overlaps with exon and extends after the exon
                        if sno_id in ['ENSG00000207145', 'ENSG00000207297', 'ENSG00000274309']:  # To patch for the 3 snoRNA that overlap with the HG transcript with the most exon, but not with the second most exon HG transcript
                            continue
                        elif sno_id not in ['ENSG00000207145', 'ENSG00000207297', 'ENSG00000274309']:
                            best_transcript_id = exon1['transcript_id']
                            gene = host[host['transcript_id'] == best_transcript_id]
                            gene.loc[:, 'sno'] = sno_id
                            gene.loc[:, 'hg_overlap'] = 'sno_over_multi_exon_after'
                            most_exon_transcripts.append(gene)
                            break
                    elif (((sno_start < exon1['start']) & (sno_end > exon1['start'])) | ((sno_start < exon2['start']) & (sno_end > exon2['start']))): # snoRNA overlaps with exon and extends before the exon
                        if sno_id == 'ENSG00000275662':  # To patch for this snoRNA that overlaps with the HG transcript with the most exon, but not with the second most exon HG transcript
                            continue
                        elif sno_id != 'ENSG00000275662':
                            best_transcript_id = exon1['transcript_id']
                            gene = host[host['transcript_id'] == best_transcript_id]
                            gene.loc[:, 'sno'] = sno_id
                            gene.loc[:, 'hg_overlap'] = 'sno_over_multi_exon_before'
                            most_exon_transcripts.append(gene)
                            break
            else:  # If the inner loop is not broken, continue to the next transcript
                continue
            break  # If the inner loop is broken, then break the outer loop


    hg_simple = pd.concat(most_exon_transcripts)
    sno_overlap = hg_simple[['sno', 'strand', 'hg_overlap', 'gene_id', 'transcript_id', 'transcript_name']].drop_duplicates()
    sno_overlap.to_csv(sno_overlap_path, index=False, sep='\t')
    hg_simple = hg_simple.loc[:, hg_simple.columns != 'hg_overlap']
    hg_simple = hg_simple.loc[:, hg_simple.columns != 'sno']

    hg_simple.to_csv(output_path, sep='\t', header=False, index=False)

    sp.call('sort -k1,1 -k2,2n -k3,17 -o '+output_path+' '+output_path, shell=True) #sort the output file by chr and start
    sp.call("""awk -i inplace -v OFS='\t' '$1="chr"$1' """+output_path, shell=True) #add "chr" in front of first column

    print('Finished get_max_exon_transcript_per_hg!')

    return hg_simple, sno_overlap


def get_exon_number_per_hg(hg_bed_file):
    """ Extract the number of exons per HG (from the transcripts chosen with
        get_max_exon_transcript_per_hg())."""

    cols = ['chr', 'start', 'end', 'gene_id', 'dot', 'strand', 'source', 'feature', 'dot', 'gene_id2', 'transcript_id',
            'exon_number', 'gene_name', 'biotype', 'transcript_name', 'transcript_biotype']

    hg_bed_file.columns = cols
    exon_number_df = hg_bed_file.groupby('transcript_id')

    rows = []
    for transcript_id, group in hg_bed_file.groupby('transcript_id'):
        max_nb = max(map(int, group['exon_number']))
        row = [transcript_id, max_nb]
        rows.append(row)
    max_nb_df = pd.DataFrame(rows, columns = ['transcript_id', 'exon_number_per_hg'])

    print('Finished get_exon_number_per_hg!')

    return max_nb_df


def get_up_downstream_exons(hg_bed_file_path, snoRNA_bed_file_path, output_path, sno_overlap_df):
    """Create two bed files giving the number of the exon located either upstream or downstream of snoRNAs and the
        distance to that exon using bedtools closest (ex: upstream exon number 4 means that the snoRNA
        is located in the 4th intron of its HG; downstream exon number 8 means that the snoRNA is located in the 7th
        intron of the HG. Since snoRNAs in the same HG have not the same corresponding HG transcript, we need to use
        bedtools closest for each snoRNa and its corresponding HG transcript (given by the sno_overlap_df)."""

    sno_bed = pd.read_csv(snoRNA_bed_file_path, sep='\t', names=['chr', 'start', 'end', 'sno_id', 'dot', 'strand', 'source', 'feature', 'dot2', 'gene_info'])
    hg_bed = pd.read_csv(hg_bed_file_path, sep='\t', names=['chr', 'start', 'end', 'host_id', 'dot', 'strand', 'source', 'feature', 'dot2', 'gene_id2', 'transcript_id', 'exon_number', 'host_name', 'biotype', 'transcript_name', 'transcript_biotype'])
    sno_overlap = sno_overlap_df.set_index('sno')
    sno_overlap_dict = sno_overlap.to_dict('index')

    for i, sno in enumerate(list(pd.unique(sno_bed.loc[:, 'sno_id']))):
        print(sno)
        temp_sno_df = sno_bed[sno_bed['sno_id'] == sno]
        temp_sno_df.to_csv(snoRNA_bed_file_path+'_temp.bed', sep='\t', header=False, index=False)
        a = BedTool(snoRNA_bed_file_path+'_temp.bed')

        hg_transcript_id = sno_overlap_dict[sno]['transcript_id']
        temp_hg_df = hg_bed[hg_bed['transcript_id'] == hg_transcript_id]
        temp_hg_df.to_csv(hg_bed_file_path+'_temp.bed', sep='\t', header=False, index=False)

        upstream = a.closest(hg_bed_file_path+'_temp.bed', t="first", io=True, id=True, D="a").saveas(output_path + 'exon_upstream_of_sno_'+str(i)+'_TEMP')
        downstream = a.closest(hg_bed_file_path+'_temp.bed', t="first", io=True, iu=True, D="a").saveas(output_path + 'exon_downstream_of_sno_'+str(i)+'_TEMP2')

    sp.call("cat "+output_path+"*_TEMP > "+output_path+"exon_upstream_of_sno.bed && rm "+output_path+"*_TEMP", shell=True)
    sp.call("cat "+output_path+"*_TEMP2 > "+output_path+"exon_downstream_of_sno.bed && rm "+output_path+"*_TEMP2 && rm "+hg_bed_file_path+"*_temp.bed", shell=True)


    print('Finished get_up_downstream_exons!')


def get_up_downstream_exons_snhg14(hg_bed_file_path, snoRNA_bed_file_path, output_path):
    """Create two bed files giving the number of the exon located either upstream or downstream of snoRNAs within SNHG14 HG and the
        distance to that exon using bedtools closest (ex: upstream exon number 4 means that the snoRNA
        is located in the 4th intron of its HG; downstream exon number 8 means that the snoRNA is located in the 7th
        intron of the HG"""

    a = BedTool(snoRNA_bed_file_path)
    upstream = a.closest(hg_bed_file_path, t="first", io=True, id=True, D="a").saveas(output_path + 'exon_upstream_of_sno.bed')
    downstream = a.closest(hg_bed_file_path, t="first", io=True, iu=True, D="a").saveas(output_path + 'exon_downstream_of_sno.bed')

    print('Finished get_up_downstream_exons!')


def get_intron_number_and_distances(path_to_exon_files, df_total_nb_exons_per_hg, output_file_path, sno_overlap_df):
    """Extract the intron number in which a snoRNA is located (i.e. the exon number of the upstream exon) and the
        distance to the upstream and downstream exons. Extract also from other df the number of exons per HG."""
    cols = ['chr_sno', 'start_sno', 'end_sno', 'gene_id_sno', 'dot', 'strand_sno', 'source_sno', 'feature_sno', 'dot2', 'gene_info_sno', 'chr_host', 'start_host',
            'end_host', 'gene_id_host', 'dot3', 'strand_host', 'source_host', 'feature_host', 'dot4', 'gene_id2_host',
            'transcript_id_host', 'intron_number', 'gene_name_host', 'biotype_host', 'transcript_name_host',
            'transcript_biotype_host']
    upstream = pd.read_csv(path_to_exon_files+'/exon_upstream_of_sno_2.bed', sep='\t', names=cols+['distance_upstream_exon'])
    downstream = pd.read_csv(path_to_exon_files+'/exon_downstream_of_sno_2.bed', sep='\t', names=cols+['distance_downstream_exon'])
    upstream, downstream = upstream.reset_index(), downstream.reset_index()

    # Remove minus in front of distances
    upstream['distance_upstream_exon'] = abs(upstream['distance_upstream_exon'])
    downstream['distance_downstream_exon'] = abs(downstream['distance_downstream_exon'])

    # Merge dfs to keep relevant info only
    final_df = upstream[['gene_id_sno', 'start_sno', 'end_sno', 'strand_sno', 'gene_id_host', 'transcript_id_host', 'intron_number', 'distance_upstream_exon']].merge(
        downstream[['gene_id_sno', 'start_sno', 'distance_downstream_exon']], how='left',
        left_on=['gene_id_sno', 'start_sno'], right_on=['gene_id_sno', 'start_sno'])

    # For snoRNAs that overlap with their HG, correct their distance to upstream and downstream exon in the df
    sno_into_single_exon = sno_overlap_df[sno_overlap_df['hg_overlap'] == 'sno_into_single_exon'].set_index('sno')
    sno_into_multi_exon = sno_overlap_df[sno_overlap_df['hg_overlap'] == 'sno_into_multi_exon'].set_index('sno')
    sno_over_multi_exon_after = sno_overlap_df[sno_overlap_df['hg_overlap'] == 'sno_over_multi_exon_after'].set_index('sno')
    sno_over_multi_exon_before = sno_overlap_df[sno_overlap_df['hg_overlap'] == 'sno_over_multi_exon_before'].set_index('sno')

    into_single, into_multi = sno_into_single_exon.to_dict('index'), sno_into_multi_exon.to_dict('index')
    over_after, over_before = sno_over_multi_exon_after.to_dict('index'), sno_over_multi_exon_before.to_dict('index')

    for d in [into_single, into_multi, over_after, over_before]:  # iterate over the dictionaries of different types of overlapping snoRNA
        for sno, cols in d.items():
            final_df.loc[(final_df.gene_id_sno == sno), 'gene_id_host'] = cols['gene_id']
            final_df.loc[(final_df.gene_id_sno == sno), 'transcript_id_host'] = cols['transcript_id']
            final_df.loc[(final_df.gene_id_sno == sno), 'intron_number'] = 0
            if cols['hg_overlap'] in ['sno_into_single_exon', 'sno_into_multi_exon']:  # if sno into single or mutliple exon, replace distance upstream/downstream exon to 0
                final_df.loc[(final_df.gene_id_sno == sno), 'distance_upstream_exon'] = 0
                final_df.loc[(final_df.gene_id_sno == sno), 'distance_downstream_exon'] = 0
            elif cols['hg_overlap'] == 'sno_over_multi_exon_after':  # if sno overlaps and extends after exon
                if cols['strand'] == "-":  # if sno on minus (-) strand, only set the distance to downstream exon to 0
                    final_df.loc[(final_df.gene_id_sno == sno), 'distance_downstream_exon'] = 0
                else:  # if sno on plus (+) strand, only set the distance to upstream exon to 0
                    final_df.loc[(final_df.gene_id_sno == sno), 'distance_upstream_exon'] = 0
            elif cols['hg_overlap'] == 'sno_over_multi_exon_before':  # if sno overlaps and extends before exon, only set the distance to downstream exon to 0
                final_df.loc[(final_df.gene_id_sno == sno), 'distance_downstream_exon'] = 0


    # Add exon number per HG for each sno
    final_df = final_df.merge(df_total_nb_exons_per_hg, how='left', left_on='transcript_id_host', right_on='transcript_id')

    # Add intron length; for snoRNAs that overlap entirely with their HG, correct intron_length for 0
    final_df['intron_length'] = final_df['end_sno'] - final_df['start_sno'] + 1 + final_df['distance_upstream_exon'] + final_df['distance_downstream_exon']
    final_df.loc[(final_df['intron_number'] == 0) & (final_df['distance_upstream_exon'] == 0) & (final_df['distance_downstream_exon'] == 0), 'intron_length'] = 0
    final_df.to_csv(output_file_path, sep='\t', index=False)

    print('Finished get_intron_number_and_distances!')



def main(gtf_bed_file, sno_info_df, output_path, sno_HG_coordinates, sno_bed_wo_snhg14_path, snhg14_bed_path, snhg14_sno_bed_path, output_table_path):

    # Create the output directory
    sp.call("""mkdir -p """ + output_path, shell=True)

    # Generate bed file of HG transcripts except for SNHG14
    hg_bed_wo_snhg14 = generate_hg_bed(gtf_bed_file, sno_info_df,
                                        output_path+'/hg_wo_snhg14.bed',
                                        output_path+'/hg_split_wo_snhg14.bed')

    # Generate from hg_bed file a simpler bed with only the longest transcript per HG
    hg_simple, sno_overlap = get_max_exon_transcript_per_hg(sno_HG_coordinates, hg_bed_wo_snhg14, output_path+'/hg_simple_wo_snhg14_sorted.bed', output_path+'/sno_overlap_hg.tsv')

    # Get the maximal number of exons per HG
    exon_number_per_hg = get_exon_number_per_hg(hg_simple)

    # Append SNHG14 row to exon_number_per_hg df
    snhg14_df_dict = {'transcript_id':'NR_146177.1', 'exon_number_per_hg':148}  # To patch for SNHG14 long transcript that has 148 exons
    snhg14_df = pd.DataFrame(snhg14_df_dict, index=[0])
    exon_number_per_hg = exon_number_per_hg.append(snhg14_df, ignore_index=True)

    # Get the exon upstream and downstream of each snoRNA (either in all HG except SNHG14 or the SNHG14 HG)
    get_up_downstream_exons(output_path+'/hg_simple_wo_snhg14_sorted.bed', sno_bed_wo_snhg14_path, output_path+'/wo_snhg14_', sno_overlap)

    get_up_downstream_exons_snhg14(snhg14_bed_path, snhg14_sno_bed_path, output_path+'/snhg14_')

    # Concat the dfs wo and with snhg14 together
    sp.call('cat '+output_path+'/wo_snhg14_exon_downstream_of_sno.bed '+output_path+'/snhg14_exon_downstream_of_sno.bed >'+output_path+'/exon_downstream_of_sno_2.bed', shell=True)
    sp.call('cat '+output_path+'/wo_snhg14_exon_upstream_of_sno.bed '+output_path+'/snhg14_exon_upstream_of_sno.bed >'+output_path+'/exon_upstream_of_sno_2.bed', shell=True)

    get_intron_number_and_distances(output_path, exon_number_per_hg, output_table_path, sno_overlap)

    print('Main script ran successfully!')



main(gtf_bed, sno_info_wo_snhg14_snornas, snakemake.params.sno_location_exon_dir, sno_HG_coordinates,
    sno_bed_wo_snhg14_path, snhg14_bed_path, snhg14_sno_bed_path, snakemake.output.output_table)
