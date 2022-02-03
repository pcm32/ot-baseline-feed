#!/usr/bin/env python
import os.path
import pandas as pd
import json
import xml.etree.ElementTree as ET
import argparse
from os.path import basename
from itertools import chain
from re import sub
import bz2

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('-p', '--path-to-experiment', required=True,
                        help="Path to the directory where the experiment is, accession inferred from basename. Should "
                             "contain <accession>.condensed-sdrf.tsv, <>-idf.txt, <>-configuration.xml, <>-fpkms.tsv, "
                             "<>-tpms.tsv. Usually $ATLAS_EXPS/<accession>.")
arg_parser.add_argument('-n', '--non-aggregated-data-dir', required=False,
                        help="Path to the dir where experiments non-aggregated data is available, "
                             "to produce the non-aggregated data files. It expects to find files "
                             "<>-fpkms.tsv.undecorated and <>-tpms.tsv.undecorated")
arg_parser.add_argument('-o', '--output', required=True,
                        help="Directoy path for output. Two files created, one for data and one for metadata"
                        )

FIXTURE_DIR = os.path.join(
    os.path.dirname(os.path.realpath(__file__)),
    'test-data',
    )


def read_condensed(condensed_path):
    cond_cols = ['Accession', 'Array', 'Sample', 'Annot_type', 'Annot', 'Annot_value', 'Annot_ont_URI']
    return pd.read_csv(condensed_path, sep="\t", names=cond_cols)


def read_configuration(configuration_path):
    with open(configuration_path, 'r') as conf_file:
        return ET.parse(conf_file).getroot()


def produce_metadata(acc, exp_type, species, speciesOnt, literature_ids, provider, assay_dictionary):
    """
    :param acc: obtained from condensed or path
    :param exp_type: obtained from configuration file
    :param species: obtained from condensed SDRF
    :param speciesOnt: obtained from condensed SDRF
    :param literature_ids: obtained from IDF
    :param provider: IDF Person last name and First name, if existing.
    :param assay_dictionary: a combined object made from the condensed and the configuration file.
    :return:
    """

    metadata = dict()
    metadata['experimentId'] = acc
    metadata['experimentType'] = exp_type
    metadata['species'] = species
    if speciesOnt:
        metadata['speciesOntURI'] = speciesOnt
    metadata['pubmedIds'] = literature_ids
    metadata['provider'] = provider
    metadata['experimentalDesigns'] = assay_dictionary

    return metadata


def get_annotation_value_from_condensed(condensed_df, annotation, sample=None):
    """
    Given a condensed SDRF as panda data frame, it will provide the first value
    for an annotation, either for a sample or for all of them.

    :param condensed_df:
    :param annotation:
    :param sample: optional, useful to leave as None for organisms.
    :return:
    >>> acc = "E-MTAB-4754"
    >>> cond_df = read_condensed(f"{FIXTURE_DIR}/{acc}/{acc}.condensed-sdrf.tsv")
    >>> get_annotation_value_from_condensed(cond_df, annotation="organism")
    'Homo sapiens'
    >>> get_annotation_value_from_condensed(cond_df, annotation="disease", sample="ERR732486")
    'normal'
    """
    df = condensed_df

    if sample:
        res = df[(df['Sample'] == sample) & (df['Annot'] == annotation)].Annot_value
    else:
        res = df[df['Annot'] == annotation].Annot_value
    if len(res):
        return res.iloc[0]
    else:
        return None


def get_all_annotation_fields_for_sample_from_condensed(condensed_df, sample):
    """
    Given a condensed SDRF as pandas data frame, it will provide all possible
    different annotation fields for a sample.


    :param condensed_df:
    :param sample:
    :return:
    >>> acc = "E-MTAB-4754"
    >>> cond_df = read_condensed(f"{FIXTURE_DIR}/{acc}/{acc}.condensed-sdrf.tsv")
    >>> get_all_annotation_fields_for_sample_from_condensed(cond_df, sample="ERR732486")
    ['age', 'biosource type', 'cell type', 'disease', 'ethnic group', 'organism', 'organism part', 'phenotype', 'sex']
    """
    df = condensed_df

    return df[df['Sample'] == sample].Annot.unique().tolist()


def get_ontology_annotation_from_condensed(condensed_df, annotation, sample=None):
    """
    Given a condensed SDRF as panda data frame, it will provide the first value
    for an annotation, either for a sample or for all of them.

    :param condensed_df:
    :param annotation:
    :param sample: optional, useful to leave as None for organisms.
    :return:
    >>> acc = "E-MTAB-4754"
    >>> cond_df = read_condensed(f"{FIXTURE_DIR}/{acc}/{acc}.condensed-sdrf.tsv")
    >>> get_ontology_annotation_from_condensed(cond_df, annotation="organism")
    'http://purl.obolibrary.org/obo/NCBITaxon_9606'
    >>> get_ontology_annotation_from_condensed(cond_df, annotation="disease", sample="ERR732486")
    'http://purl.obolibrary.org/obo/PATO_0000461'
    """
    df = condensed_df

    if sample:
        res = df[(df['Sample'] == sample) & (df['Annot'] == annotation)].Annot_ont_URI
    else:
        res = df[df['Annot'] == annotation].Annot_ont_URI
    if len(res) and not pd.isna(res.iloc[0]):
        return res.iloc[0]
    else:
        return None


def produce_assay_dict(condensed_df, configuration_xml):
    """
    Produces a dictionary of structure:

    - id:
      organismPart:
      sex:
      developmentalStage:
      assay_group:
      technical_replicate_group:
    - id:
      organismPart:
      ...

    id is the assay identifier (from configuration file); organism part comes from condensed SDRF;
    sex from condensed SDRF; assay_group from configuration file; technical_replicate_group comes from the configuration,
    and is optional.

    :param condensed_pd:
    :param configuration_xml:
    :return:

    >>> acc = "E-MTAB-4754"
    >>> cond_df = read_condensed(f"{FIXTURE_DIR}/{acc}/{acc}.condensed-sdrf.tsv")
    >>> conf_xml = read_configuration(f"{FIXTURE_DIR}/{acc}/{acc}-configuration.xml")
    >>> assay_dict = produce_assay_dict(cond_df, conf_xml)
    >>> assay_dict[0]['assayId']
    'ERR732483'
    >>> len(assay_dict)
    4
    """

    assay_dict = []
    for ag in configuration_xml.iter("assay_group"):
        for assay in ag:
            entry = {'assayGroupId': ag.get("id"), 'assayId': assay.text, 'assayGroup': ag.get('label')}

            for annot in get_all_annotation_fields_for_sample_from_condensed(condensed_df, assay.text):
                value = get_annotation_value_from_condensed(condensed_df, annot, sample=assay.text)
                if annot == 'organism':
                    continue
                if value:
                    entry[camel_case(annot)] = value
                onto_annot = get_ontology_annotation_from_condensed(condensed_df, annot, sample=assay.text)
                if onto_annot:
                    entry[f"{camel_case(annot)}OntURI"] = onto_annot

            if 'technical_replicate_id' in assay.attrib:
                entry['technicalReplicateId'] = assay.attrib['technical_replicate_id']

            assay_dict.append(entry)
    return assay_dict


def read_idf_to_dict(idf_path):
    """
    Reads an IDF file into a dictionary, where the first element of each tab separated
    line is the key for the rest of the line.
    :param idf_path:
    :return:

    >>> acc = "E-MTAB-4754"
    >>> idf_obj = read_idf_to_dict(f"{FIXTURE_DIR}/{acc}/{acc}.idf.txt")
    >>> 'person last name' in idf_obj
    True
    >>> idf_obj['date of experiment']
    '2015-01-09'
    """
    res = {}
    with open(idf_path, 'r') as idf_file:
        for line in idf_file:
            tokens = line.strip().split("\t")
            res[tokens[0].lower()] = "\t".join(tokens[1:len(tokens)])
    return res


def get_expression_units(experiment_type):
    """
    Provides the expression units that are applicable for experiment type
    TODO this probably should be somehow derived from files, instead of hardcoded here.

    :param experiment_type: as reported in the configuration XML file
    :return: array with units that we expect to find for that exp type
    """

    if experiment_type == "rnaseq_mrna_baseline":
        return ["tpms", "fpkms", "transcript-tpms"]
    if experiment_type == "proteomics_baseline":
        return ["ppb"]
    if experiment_type == "proteomics_baseline_dia":
        return ["aa"]


def get_file_label(experiment_type, expression_unit):
    """
    Currently Atlas file naming requires some corrections of the file names for expression data
    on some cases.

    :param experiment_type:
    :param expression_unit:
    :return:
    """

    if "proteomics_baseline" in experiment_type:
        return ""
    else:
        return f"-{expression_unit}"


def quartile_data_iterator(data_path, metadata, unit):
    """
    Iterator method to get dictionary of expression (min, q1, q2, q3, max) per gene, as it traverses
    an atlas data file, on a per assay_group level.

    :param data_path: path to the expression data file (decorated)
    :param metadata: the assay dictionary as produced by produce_assay_dict
    :param unit: the unit to be used (tpms, fpkms, etc)
    :return:
    """
    assays_per_group = {}
    # first we create a dictionary of assay groups to their ordered assays in arrays.
    # for assay in metadata:
    #     if not assay['assay_group'] in assays_per_group:
    #         assays_per_group[assay['assay_group']] = []
    #     assays_per_group[assay['assay_group']].append(assay['id'])

    with open(data_path, 'r') as data_file:
        header = data_file.readline().strip().split("\t")
        assay_groups = header[2:len(header)]

        check_biyection(assay_groups, metadata)

        for line in data_file:
            tokens = line.strip().split("\t")
            gene_id = tokens[0]
            entry = {"geneProductId": gene_id, "expression": [], "unit": unit}

            for ag_i in range(0, len(assay_groups)):
                ag = assay_groups[ag_i]
                ag_data = tokens[2+ag_i].split(",")
                expression_entry = {'assayGroupId': ag,
                                    "min": float(ag_data[0]),
                                    "q1": float(ag_data[1]),
                                    "q2": float(ag_data[2]),
                                    "q3": float(ag_data[3]),
                                    "max": float(ag_data[4]) }
                entry['expression'].append(expression_entry)


            yield entry


def check_biyection(elements_in_data, metadata, type="assayGroupId"):
    e_in_metadata = set()
    for m in metadata:
        e_in_metadata.add(m[type])
    e_in_data = set(elements_in_data)
    e_only_metadata = e_in_metadata.difference(e_in_data)
    if e_only_metadata:
        print(f"WARNING: Assay groups only in meta data: {', '.join(e_only_metadata)}")
    e_only_data = e_in_data.difference(e_in_metadata)
    if e_only_data:
        print(f"WARNING: Assay groups only in data: {', '.join(e_only_data)}")


def get_provider(idf_dict):
    """
    The provider gets generated from the IDF Person fields, merging last name and first names
    in the order they appear in columns.

    There is not always guarantee that an IDF's person field will actually list people
    and not other things like institutions, so this is a bit dirty.

    :param idf_dict: as read from the IDF file.
    :return:

    >>> acc = "E-MTAB-4754"
    >>> idf_obj = read_idf_to_dict(f"{FIXTURE_DIR}/{acc}/{acc}.idf.txt")
    >>> get_provider(idf_obj)
    'BLUEPRINT consortium; ArrayExpress'
    """

    f"{idf_dict['person last name']}, {idf_dict['person first name']}"
    first_names = []
    if 'person first name' in idf_dict:
        first_names = idf_dict['person first name'].split("\t")
    last_names = []
    if 'person last name' in idf_dict:
        last_names = idf_dict['person last name'].split("\t")
    if last_names and first_names and len(last_names) == len(first_names):
        # both arrays have the same number of objects
        return "; ".join(list(chain.from_iterable(zip(first_names, last_names))))
    if last_names:
        return "; ".join(last_names)
    else:
        return "Not available"


def camel_case(s):
    s = sub(r"(_|-)+", " ", s).title().replace(" ", "")
    return ''.join([s[0].lower(), s[1:]])


def unaggregated_data_iterator(data_path, assay_dictionary, unit):
    """
    Iterator method to get dictionary of expression per gene, as it traverses
    an atlas data file.
    :param data_path: path to the unaggregated expression data file (decorated)
    :param metadata: the assay dictionary as produced by produce_assay_dict
    :param unit: the unit to be used (tpms, fpkms, etc)
    :return:
    """

    first_data_index = 1
    with open(data_path, 'r') as data_file:
        header = data_file.readline().strip().split("\t")
        assays = header[first_data_index:len(header)]

        check_biyection(assays, assay_dictionary, type="assayId")

        for line in data_file:
            tokens = line.strip().split("\t")
            gene_id = tokens[0]
            entry = {"geneProductId": gene_id, "unit": unit, "expression": []}

            for a_i in range(0, len(assays)):
                assay = assays[a_i]
                exp = tokens[first_data_index + a_i]
                expression_entry = {'assayId': assay, "value": float(exp)}
                entry['expression'].append(expression_entry)

            yield entry


if __name__ == '__main__':
    args = arg_parser.parse_args()

    acc = basename(args.path_to_experiment)
    condensed_path = f"{args.path_to_experiment}/{acc}.condensed-sdrf.tsv"
    condensed_pd = read_condensed(condensed_path)

    configuration_path = f"{args.path_to_experiment}/{acc}-configuration.xml"
    conf_xml = read_configuration(configuration_path)
    expType = conf_xml.attrib['experimentType']

    assay_dictionary = produce_assay_dict(condensed_pd, conf_xml)

    idf_path = f"{args.path_to_experiment}/{acc}.idf.txt"
    idf_dict = read_idf_to_dict(idf_path)

    lit_ids = []
    if 'pubmed id' in idf_dict:
        lit_ids = idf_dict['pubmed id'].split("\t")
    provider = get_provider(idf_dict)
    metadata = produce_metadata(acc=acc,
                                exp_type=expType,
                                species=get_annotation_value_from_condensed(condensed_pd, annotation="organism"),
                                speciesOnt=get_ontology_annotation_from_condensed(condensed_pd, annotation="organism"),
                                literature_ids=lit_ids,
                                provider=provider,
                                assay_dictionary=assay_dictionary
                                )

    with open(f"{args.output}/{acc}.metadata.json", 'w') as metadata_out:
        json.dump(metadata, metadata_out)

    for unit in get_expression_units(experiment_type=expType):
        file_label = get_file_label(expType, unit)
        data_file = f"{args.path_to_experiment}/{acc}{file_label}.tsv"
        if os.path.isfile(data_file):
            with bz2.open(f"{args.output}/{acc}-expression-data-{unit}.jsonl.bz2", 'wt') as data_jsonl:
                for exp_d in quartile_data_iterator(data_file, assay_dictionary, unit):
                    data_jsonl.write(json.dumps(exp_d)+"\n")

        if args.non_aggregated_data_dir:
            unaggregated_data_file = f"{args.non_aggregated_data_dir}/{acc}{file_label}.tsv.undecorated"
            
            if os.path.isfile(unaggregated_data_file):
                with bz2.open(f"{args.output}/{acc}-unaggregated-expression-data-{unit}.jsonl.bz2", "wt") as data_jsonl:
                    for exp_d in unaggregated_data_iterator(unaggregated_data_file, assay_dictionary, unit):
                        data_jsonl.write(json.dumps(exp_d)+"\n")











