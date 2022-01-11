#!/usr/bin/env python
import os.path
import sys
import pandas as pd
import json
import xml.etree.ElementTree as ET
import argparse
from os.path import basename

arg_parser = argparse.ArgumentParser()
arg_parser.add_argument('-p', '--path-to-experiment', required=True,
                        help="Path to the directory where experiment, accession inferred from basename.")
arg_parser.add_argument('-o', '--output', required=True,
                        help="Directoy path for output. Two files created, one for data and one for metadata"
                        )


def read_condensed(condensed_path):
    cond_cols = ['Accession', 'Array', 'Sample', 'Annot_type', 'Annot', 'Annot_value', 'Annot_ont_URI']
    cond = pd.DataFrame()
    return pd.read_csv(condensed_path, sep="\t", names=cond_cols)


def read_configuration(configuration_path):
    with open(configuration_path, 'r') as conf_file:
        return ET.parse(conf_file).getroot()


def produce_metadata(acc, expType, species, literature_ids, provider, assay_dictionary):
    """

    :param acc: obtained from condensed or path
    :param expType: obtained from configuration file
    :param species: obtained from condensed SDRF
    :param literature_ids: obtained from IDF
    :param provider: IDF Person last name and First name, if existing.
    :param assay_dictionary: a combined object made from the condensed and the configuration file.
    :return:
    """

    metadata = dict()
    metadata['experimentId'] = acc
    metadata['experimentType'] = expType
    metadata['species'] = species
    metadata['literature'] = literature_ids
    metadata['provider'] = provider
    metadata['experimentalDesigns'] = assay_dictionary

    return metadata


def get_annotation_value_from_condensed(condensed_df, annotation, sample=None):
    df = condensed_df

    if sample:
        res = df[(df['Sample'] == sample) & (df['Annot'] == annotation)].Annot_value
    else:
        res = df[df['Annot'] == annotation].Annot_value
    if len(res):
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
    """

    assay_dict = []
    for ag in configuration_xml.iter("assay_group"):
        for assay in ag:
            entry = {'assay_group': ag.get("id"), 'id': assay.text}

            for annot in ["sex", "organism part"]:
                value = get_annotation_value_from_condensed(condensed_df, annot, sample=assay.text)
                if value:
                    entry[annot] = value

            if 'technical_replicate_id' in assay.attrib:
                entry['technical_replicate_id'] = assay.attrib['technical_replicate_id']

            assay_dict.append(entry)
    return assay_dict


def read_idf_to_dict(idf_path):
    """
    Reads an IDF file into a dictionary, where the first element of each tab separated
    line is the key for the rest of the line.
    :param idf_path:
    :return:
    """
    res = {}
    with open(idf_path, 'r') as idf_file:
        for line in idf_file:
            tokens = line.strip().split("\t")
            res[tokens[0].lower()] = "\t".join(tokens[1:len(tokens)])
    return res


def data_iterator(data_path, metadata, unit):
    """
    Iterator method to get dictionary of expression per gene, as it traverses
    an atlas data file.

    :param data_path:
    :param metadata:
    :param unit:
    :return:
    """
    assays_per_group = {}
    # first we create a dictionary of assay groups to their ordered assays in arrays.
    for assay in metadata:
        if not assay['assay_group'] in assays_per_group:
            assays_per_group[assay['assay_group']] = []
        assays_per_group[assay['assay_group']].append(assay['id'])

    with open(data_path, 'r') as data_file:
        header = data_file.readline().strip().split("\t")
        assay_groups = header[2:len(header)]

        for line in data_file:
            tokens = line.strip().split("\t")
            gene_id = tokens[0]
            entry = {"geneProductId": gene_id, "expression": [], "unit": unit}

            for ag_i in range(0, len(assay_groups)):
                ag = assay_groups[ag_i]
                ag_data = tokens[2+ag_i].split(",")
                for a_i in range(0, len(assays_per_group[ag])):
                    assay = assays_per_group[ag][a_i]
                    exp = ag_data[a_i]
                    expression_entry = {'experimentDesignID': assay, "value": exp }
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

    metadata = produce_metadata(acc=acc,
                                expType=expType,
                                species=get_annotation_value_from_condensed(condensed_pd, annotation="organism"),
                                literature_ids=idf_dict['pubmed id'].split("\t"),
                                provider=f"{idf_dict['person last name']}, {idf_dict['person first name']}",
                                assay_dictionary=assay_dictionary
                                )

    with open(f"{args.output}/{acc}.metadata.json", 'w') as metadata_out:
        json.dump(metadata, metadata_out)

    for type in ["tpms", "fpkms"]:
        data_file = f"{args.path_to_experiment}/{acc}-{type}.tsv"
        if os.path.isfile(data_file):
            with open(f"{args.output}/{acc}-expression-data-{type}.jsonl", 'w') as data_jsonl:
                for exp_d in data_iterator(data_file, assay_dictionary, type):
                    data_jsonl.write(json.dumps(exp_d)+"\n")











