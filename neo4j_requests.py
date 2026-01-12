#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  7 16:03:15 2025

@author: fgraziani
"""

import re
from tqdm import tqdm
from math import *
from numpy.polynomial.polynomial import Polynomial
from neo4j import GraphDatabase
from neo4j.exceptions import ServiceUnavailable
import time
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import webbrowser
import os
import shutil
import glob
from typing import List, Dict
import logging
import csv
import json
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from scipy.spatial.distance import pdist, squareform
import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Seq import Seq
import statistics
import logging
from config import *


logger = logging.getLogger("panabyss_logger")

# CONF_FILE = os.path.abspath("./conf.json")
#
# DEFAULT_MAX_NODES_NUMBER = 30000
#
# def get_max_nodes_number():
#     if not os.path.exists(CONF_FILE):
#         return DEFAULT_MAX_NODES_NUMBER
#     else:
#         with open(CONF_FILE) as f:
#             conf = json.load(f)
#             return int(conf.get("max_nodes_to_visualize", DEFAULT_MAX_NODES_NUMBER))

#This value is used to limit the nodes number when seeking for regions :
#If the region is wider than this value then it is ignored
MAX_BP_SEEKING = 800000

#Maximal number of nodes to get the whole region
MAX_NODES_NUMBER = get_max_nodes_to_visualize()

logging.getLogger("neo4j").setLevel(logging.ERROR)
logger.debug(f"Max nodes number : {MAX_NODES_NUMBER}")

def get_anchor(genome, chromosome, position, before = True, use_anchor=True):
    core_genome = False
    if get_driver() is None :
        return None
    genome_position = genome+"_position"
    if use_anchor :
        window_size=1000
        max_attemps = 500
        attempt = 0
        lower_bound = position
        upper_bound = position
        with get_driver() as driver:
            with driver.session() as session:
                while attempt < max_attemps and lower_bound > 0:
                    attempt += 1
                    if before:
                        upper_bound = lower_bound
                        lower_bound = lower_bound - window_size
                        order = "DESC"
    
                    else:
                        lower_bound = upper_bound
                        upper_bound = upper_bound + window_size
                        order = "ASC"
            
                    query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` >= $lower_bound
                      AND n.`{genome_position}` <= $upper_bound
                      AND n.flow >= 1.0
                    RETURN n
                    ORDER BY n.`{genome_position}` {order}
                    LIMIT 1
                    """
    
                    result = session.run(
                            query,
                            chromosome=chromosome,
                            lower_bound=lower_bound,
                            upper_bound=upper_bound
                        )
                    record = result.single()
                    if record:
                        core_genome = True
                        return dict(record["n"]),core_genome
                #Anchor not found => the current node will be used
                logger.debug("No core genome anchor found, get the current nodes.")
                if before:
                    query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` <= $position
                    RETURN n order by n.`{genome_position}` DESC limit 1
                    """
                else :
                    query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` >= $position
                    RETURN n order by n.`{genome_position}` ASC limit 1
                    """
    
                result = session.run(
                        query,
                        chromosome=chromosome,
                        position=position
                    )
                record = result.single()
                if record:
                    return dict(record["n"]),core_genome
    else:
        if before:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` >= $position
            RETURN n order by n.`{genome_position}` ASC limit 1
            """
        else :
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` <= $position
            RETURN n order by n.`{genome_position}` DESC limit 1
            """
        with get_driver() as driver:
            with driver.session() as session:
                result = session.run(
                        query,
                        chromosome=chromosome,
                        position=position
                    )
                record = result.single()
                if record:
                    return dict(record["n"]) ,core_genome

    return None, core_genome



#This function take a region (chromosome, start and stop) of a given haplotype (search_genome)
#and it returns all the nodes in this region and the other related regions : 
#for each haplotype the start and stop are given by the first anchor node before the start position and the first after the end position
#Anchor node is a node with all haplotypes
#return metadata :
# return_code :
#   - OK => nodes found
#   - PARTIAL => no anchor on pangenome found
#   - WIDE => more than MAX_NODES_NUMBER nodes found into the region
#   - NO_DATA => no data found
#   - ZOOM => more than MAX_NODES_NUMBER nodes found but with filtering with flow an acceptable number of nodes has been found
#       in this case the flow value will be set
#   - FILTER => exceptional individuals have been removed from the search (if not there would be too much nodes)
# flow : in case of filtering, the minimal flow used to filter data
def get_nodes_by_region(genome, chromosome, start, end, use_anchor = True ):
    return_metadata = {"return_code":"OK", "flow":None, "nodes_number":0, "removed_genomes" : None}
    valid_individuals_exceptions = []
    flow = None
    ranges = {}
    if get_driver() is None :
        return []
    temps_depart = time.time()
    nodes_data = {}
    max_sequence = 1000
    if end is None:
        stop = MAX_BP_SEEKING
    else:
        stop = end
    with get_driver() as driver:

        data = {}
        shared_genomes = []
        with driver.session() as session:
    
            # Step 1 : find the anchors of the region
            if start is not None and end is not None :
                genome_position = genome+"_position"
                logger.info("Look for region : " + str(start) + " - " + str(stop) + " - chromosome " + str(chromosome) + " - genome : " + str(genome))

                anchor_start, core_genome_start = get_anchor(genome, chromosome, start, before = True, use_anchor=use_anchor)
                anchor_stop, core_genome_stop = get_anchor(genome, chromosome, end, before = False, use_anchor=use_anchor)
                if anchor_start is None or anchor_stop is None :
                    return_metadata["return_code"] = "NO_DATA"
                    logger.warning("No data found")
                elif anchor_start is not None and anchor_stop is not None and not core_genome_start or not core_genome_stop:
                    # No core genome anchor found => search all genomes present on the nodes to get start and stop
                    return_metadata["return_code"] = "PARTIAL"
                    ref_position_start = anchor_start[genome_position]
                    ref_position_stop = anchor_stop[genome_position]
                    query_genome = """
                        MATCH (m:Node)
                        WHERE m.chromosome = "{chromosome}"
                          AND m.`{genome_position}` >= {ref_position_start}
                          AND m.`{genome_position}` <= {ref_position_stop}
                        
                        WITH m 
                        LIMIT {limit}
                        
                        WITH collect(m) AS nodes
                        UNWIND nodes AS n
                        UNWIND n.genomes AS g
                        WITH g AS genome, n[g + "_position"] AS pos
                        WHERE pos IS NOT NULL
                        WITH 
                            genome,
                            min(pos) AS start_pos,
                            max(pos) AS stop_pos
                        RETURN 
                            collect([genome, {{start: start_pos, stop: stop_pos}}]) AS genome_ranges
                        """.format(
                            chromosome=chromosome,
                            genome_position=genome_position,
                            ref_position_start=ref_position_start,
                            ref_position_stop=ref_position_stop,
                            limit=MAX_NODES_NUMBER + 1
                        )

                    result = session.run(query_genome, start=start, end=end)
                    record = result.single()
                    if record:
                        pairs = record["genome_ranges"]
                        ranges = {genome: data for genome, data in pairs}
                else:
                    ranges = {}
                    for g in anchor_start["genomes"] + anchor_stop["genomes"]:
                        p_start = anchor_start[g+"_position"]
                        p_stop = anchor_stop[g+"_position"]
                        ranges[g] = {"start":min(p_start, p_stop), "stop":max(p_start, p_stop)}

                if anchor_start is not None and  anchor_stop is not None:
                    if anchor_start[genome_position] > anchor_stop[genome_position]:
                        anchor_start_tmp = anchor_start
                        anchor_start = anchor_stop
                        anchor_stop = anchor_start_tmp
                    logger.debug("Anchor start name : " + str(anchor_start["name"]))
                    logger.debug("Anchor stop name : " + str(anchor_stop["name"]))
                    logger.debug("Anchor region : " + str(anchor_start[genome_position]) + " - " + str(anchor_stop[genome_position]))

                # Step 2 : find the region between the 2 anchors
                if anchor_start is not None and anchor_stop is not None and anchor_stop[genome_position] - anchor_start[genome_position] > 0 and len(anchor_start['genomes']) > 0 :
                    region_nodes_number = 0
                    #construct the base query to find all genomes betwwen the start / stop position
                    base_query_genome = f"""
                        MATCH (m:Node)
                        WHERE  m.chromosome = "{chromosome}"
                        AND (
                        """
                    first = True
                    for g in ranges:
                        position_field = g + "_position"
                        if first:
                            base_query_genome += f"(m.`{position_field}` >= {ranges[g]['start']} AND m.`{position_field}` <= {ranges[g]['stop']})"
                            first = False
                        else:
                            base_query_genome += f" OR (m.`{position_field}` >= {ranges[g]['start']} AND m.`{position_field}` <= {ranges[g]['stop']})"
                    base_query_genome += ")"

                    # Step 3 : Check if the region size is not too  wide
                    query_genome = base_query_genome + f"""
                        WITH m LIMIT {MAX_NODES_NUMBER+1}
                        return count(m) as nodes_number
                        """
                    result = session.run(query_genome, start=start, end=end)
                    record = result.single()
                    if record:
                        region_nodes_number =  int(record["nodes_number"])
                    if region_nodes_number <= MAX_NODES_NUMBER :
                        flow = None
                    else:
                        # Step 5 : the region is too wide, check if the pb is due to a small proportion of indiviudals
                        # If the nodes number of some individuals (less than 20% of the total individuals) is more than the limit
                        # or more than 10 times the median, then they will be removed from the search
                        logger.debug("Too much nodes into the region, check for individual exception.")
                        queries = []
                        for g in ranges:
                            position_field = g + "_position"
                            q = f"""
                                MATCH (m:Node)
                                WHERE m.chromosome = "{chromosome}"
                                  AND m.{position_field} >= {ranges[g]['start']} AND m.{position_field} <= {ranges[g]['stop']}
                                WITH m LIMIT {MAX_NODES_NUMBER + 1}
                                RETURN "{g}" AS genome, count(m) AS nb
                            """
                            queries.append(q)

                        query_genome = "\nUNION ALL\n".join(queries)
                        # logger.debug(query_genome)
                        result = session.run(query_genome)
                        counts = {r["genome"]: r["nb"] for r in result}
                        # logger.debug(counts)
                        median_value = statistics.median(counts.values())
                        individuals_exceptions = []
                        valid_individuals_exceptions = []
                        for g in counts:
                            if counts[g] > MAX_NODES_NUMBER or counts[g] > 10 * median_value:
                                individuals_exceptions.append(g)
                        logger.debug(f"Exceptional individuals : {individuals_exceptions}")
                        if len(individuals_exceptions) == 1 or len(individuals_exceptions) <= 0.2 * len(counts):
                            valid_individuals_exceptions = individuals_exceptions
                            logger.debug(
                                f"These individuals will be removed from search : {valid_individuals_exceptions}")
                        else:
                            # Step 5 : the region is too wide and it is not linked to a small proportion of individudls
                            # => try to find a "zoom level" for which the number of nodes is acceptable
                            # To do that it will use the flow attribute
                            logger.debug("Too much nodes into the region, try to reduce data by filtering on flow.")
                            zoom = True
                            flow = 1
                            while zoom and flow >= 0:
                                query_genome = base_query_genome + f"""
                                    AND m.flow >= {flow} 
                                    WITH m LIMIT {MAX_NODES_NUMBER+1}
                                    return count(m) as nodes_number
                                    """
                                #logger.debug(query_genome)
                                result = session.run(query_genome, start=start, end=end)
                                record = result.single()
                                if record:
                                    region_nodes_number = int(record["nodes_number"])
                                if region_nodes_number > MAX_NODES_NUMBER:
                                    zoom = False
                                    logger.debug(f"Too much nodes with flow {flow}.")
                                else:
                                    if flow - 0.1 >= 0 :
                                        logger.debug(f"Nodes number found with flow {flow} : {region_nodes_number} - check with flow {flow - 0.1}")
                                    flow -= 0.1
                            if not zoom and flow == 1:
                                flow = -1
                            else:
                                flow += 0.1
                            if flow >= 0 :
                                logger.debug(f"Zoom level found with flow {flow}.")
                            if flow < 0:
                                logger.debug("Too much nodes into the region, no zoom level found.")

                    if flow is None or (flow is not None and flow >= 0) or len(valid_individuals_exceptions) > 0:
                        # Step 7 : Get the nodes and annotations for each genomes
                        if flow is not None and flow > 0 :
                            #The search will be filtered by flow
                            return_metadata["flow"]= flow
                            return_metadata["return_code"] = "ZOOM"
                            query_genome = base_query_genome + f"AND m.flow >= {flow} "
                        elif len(valid_individuals_exceptions) > 0:
                            #Remove exceptional individuals of the list
                            query_genome = f"""
                                MATCH (m:Node)
                                WHERE  m.chromosome = "{chromosome}"
                                AND (
                                """
                            first = True
                            for g in ranges:
                                if g not in valid_individuals_exceptions :
                                    position_field = g + "_position"
                                    if first:
                                        query_genome += f"(m.`{position_field}` >= {ranges[g]['start']} AND m.`{position_field}` <= {ranges[g]['stop']})"
                                        first = False
                                    else:
                                        query_genome += f" OR (m.`{position_field}` >= {ranges[g]['start']} AND m.`{position_field}` <= {ranges[g]['stop']})"
                            query_genome += ")"
                            return_metadata["flow"] : 0
                            return_metadata["removed_genomes"] = valid_individuals_exceptions
                            return_metadata["return_code"] = "FILTER"
                        else : query_genome = base_query_genome

                        query_genome = query_genome + f"""
                            OPTIONAL MATCH (m)-[]->(a:Annotation)
                            OPTIONAL MATCH (s:Sequence {{name: m.ref_node}})
                            RETURN m, substring(s.sequence, 0, {max_sequence}) as sequence, collect(a.gene_name) AS annotations, collect(a.feature) AS features
                            LIMIT {MAX_NODES_NUMBER + 1}
                            """
                        #logger.info(query_genome)
                        result = session.run(query_genome, start=start, end=end)
                        for record in result :
                            nodes_data[record["m"]["name"]] = dict(record["m"]) |{"sequence":record["sequence"]} |{"annotations":set(record["annotations"][a] for a in range(len(record["annotations"])))} |{"features":set(record["features"][a] for a in range(len(record["features"])))}
                        return_metadata["nodes_number"] = len(nodes_data)
                        if len(nodes_data) > MAX_NODES_NUMBER :
                            nodes_data = {}
                            return_metadata["return_code"] = "WIDE"
                            logger.warning(
                                f"Region too wide : nodes number : {len(nodes_data)} - max nodes number : {MAX_NODES_NUMBER}")

                    else:
                        return_metadata["return_code"] = "WIDE"
                        logger.warning(
                            f"Region too wide : {anchor_stop[genome_position] - anchor_start[genome_position]} - nodes number : {region_nodes_number} - max nodes number : {MAX_NODES_NUMBER}")
                        nodes_data = {}
                else:
                    if anchor_start is not None and anchor_stop is not None and anchor_stop[genome_position] - anchor_start[genome_position] >= MAX_BP_SEEKING :
                        logger.warning(f"Region too wide : {anchor_stop[genome_position] - anchor_start[genome_position]}" )
                        return_metadata["return_code"] = "WIDE"
                    else:
                        logger.warning("Region not found")
                        return_metadata["return_code"] = "NO_DATA"
                    nodes_data = {}

            else:
                #if end is set to None, get all the nodes if the number of nodes is less than MAX_NODES_NUMBER
                
                if start == 0 and end is None :
                    total_nodes = get_nodes_number(chromosome)
                    if total_nodes <= MAX_NODES_NUMBER :
                        query_genome = f"""
                        MATCH (m:Node)
                        WHERE  m.chromosome = "{chromosome}" 
                        OPTIONAL MATCH (m)-[]->(a:Annotation)
                        OPTIONAL MATCH (s:Sequence {{name: m.ref_node}})
                        RETURN m, substring(s.sequence, 0, {max_sequence}) as sequence, collect(a.gene_name) AS annotations, collect(a.feature) AS features
                        LIMIT {MAX_NODES_NUMBER + 1}
                        """
                        result = session.run(query_genome)
                        for record in result:
                            nodes_data[record["m"]["name"]] = dict(record["m"]) |{"sequence":record["sequence"]}  |{"annotations":set(record["annotations"][a] for a in range(len(record["annotations"])))} |{"features":set(record["features"][a] for a in range(len(record["features"])))}
                        return_metadata["nodes_number"] = len(nodes_data)
                        if len(nodes_data) > MAX_NODES_NUMBER :
                            nodes_data = {}
                            return_metadata["return_code"] = "WIDE"
                            logger.warning(
                                f"Region too wide : nodes number : {len(nodes_data)} - max nodes number : {MAX_NODES_NUMBER}")
                    else  :
                        logger.warning("Region too wide")
                        return_metadata["return_code"] = "WIDE"
                        nodes_data = {}
            if len(nodes_data) > 0 :
                for elt in nodes_data:
                    nodes_data[elt]["annotations"] = list(nodes_data[elt]["annotations"])
                    nodes_data[elt]["features"] = list(nodes_data[elt]["features"])
        logger.debug("Total time : " + str(time.time() - temps_depart))
        return nodes_data, return_metadata


# Find a region by annotation : gene_name or gene_id
#   - OK
#   - WIDE
#   - NO_DATA
def get_nodes_by_feature(genome, chromosome, gene_id=None, feature = None, value=None):
    return_code = "OK"
    if get_driver() is None :
        return []
    with get_driver() as driver:

        nodes_data = {}
        shared_genomes = []
    
        with driver.session() as session:
            genome_position = genome+"_position"
            # Step 1 : find nodes with gene annotation
            if feature is not None and value is not None:
                feature_name = feature + "_name"
                logger.info(f"Looking for {feature_name} : {value}")
                query = f"""
                MATCH (a:Annotation {{chromosome:"{chromosome}", {feature_name}: $value}})<-[]-(n:Node)
                WHERE n[$genome_position] IS NOT NULL
                RETURN DISTINCT n
                ORDER BY n[$genome_position] ASC
                """
                result = session.run(query, value=value, genome_position=genome_position)
                #logger.debug(query)
                noeuds_annotes = [record["n"] for record in result]
                if len(noeuds_annotes) > 0 and genome_position in noeuds_annotes[0] and genome_position in \
                        noeuds_annotes[-1]:
                    start = noeuds_annotes[0][genome_position]
                    stop = noeuds_annotes[-1][genome_position] + noeuds_annotes[-1]["size"]
                    logger.debug(f"start : {start} - stop : {stop} - nodes number : {len(noeuds_annotes)}")
                    nodes_data, return_metadata = get_nodes_by_region(genome, chromosome, start, stop)
                else:
                    logger.debug(f"No nodes found {len(noeuds_annotes)}.")
                    return_metadata = {"return_code": "NO_DATA", "flow": None, "nodes_number": 0}
            else:
                return_metadata = {"return_code": "NO_DATA", "flow":None, "nodes_number":0}


            
        return nodes_data, return_metadata



#This function get first annotation on nodes of a reference_genome on a given chromosome after the given position
def get_annotation_before_or_after_position(genome_ref, chromosome="1", position=0, before=True):
    if get_driver() is None :
        return {}
    with get_driver() as driver:
        if before :
            query = f"""
            MATCH (a:Annotation)
            WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.end < $position and a.gene_name is not null
            RETURN a.gene_name AS gene_name, a.end as end, a.start as start, a.feature as feature
            ORDER BY a.end DESC
            LIMIT 1
            """
        else:
            query = f"""
            MATCH (a:Annotation)
            WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.start > $position and a.gene_name is not null
            RETURN a.gene_name AS gene_name, a.end as end, a.start as start, a.feature as feature
            ORDER BY a.start ASC
            LIMIT 1
            """
        
        with driver.session() as session:
            result = session.run(query, position=position)
            record = result.single()
            return dict(record) if record else None


#This function get all annotations on nodes of a reference_genome on a given chromosome between start_position and end_position
def get_annotations_in_position_range(genome_ref, chromosome="1", start_position=0, end_position=0):
    if get_driver() is None :
        return []
    with get_driver() as driver:
        query = f"""
        MATCH (a:Annotation)
        WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.end <= {end_position} and a.start >= {start_position} and a.gene_name is not null
        RETURN DISTINCT a.gene_id AS gene_id, a.gene_name AS gene_name, a.feature as feature
        """ 
        
        # query = f"""
        # MATCH (n:Node)-[]->(a:Annotation)
        # WHERE n.chromosome = "{chromosome}" and n.`"""+str(genome_ref)+"""_position` >= $start AND n.`"""+str(genome_ref)+"""_position` <= $end and a.gene_name is not null and a.genome_ref = $genome_ref 
        # RETURN DISTINCT a.gene_id AS gene_id, a.gene_name AS gene_name, a.feature as feature
        # """
        #logger.info("Query : " + query)
        with driver.session() as session:
            result = session.run(query, start=start_position, end=end_position, genome_ref=genome_ref)
            annotations = [dict(record) for record in result]
    
    return annotations


# This function get all features on annotation nodes
def get_annotations_features():
    if get_driver() is None:
        return []
    with get_driver() as driver:
        query = f"""
        MATCH (a:Annotation)
        RETURN DISTINCT a.feature as feature
        """
        with driver.session() as session:
            result = session.run(query)
            features = [record["feature"] for record in result]

    return features

#This function will get all chromosomes present in the pangenome graph
def get_chromosomes():
    if get_driver() is None :
        return []
    with get_driver() as driver:
    
        query = """
        MATCH (s:Stats) 
        RETURN s.chromosomes as all_chromosomes
        """
        all_chromosomes = []
        with driver.session() as session:
            result = session.run(query)
            for record in result:
                all_chromosomes = record["all_chromosomes"]
    
    return all_chromosomes


#This function will get the number of nodes in the graph
def get_nodes_number(chromosome=None):
    if get_driver() is None :
        return 0
    total = 0
    with get_driver() as driver:
        if chromosome is None and chromosome != "":
            query = """
            MATCH (n:Node) 
            RETURN count(n) as total
            """
        else:
            query = f"""
            MATCH (n:Node) 
            where n.chromosome="{chromosome}"
            RETURN count(n) as total
            """
        all_chromosomes = []
        with driver.session() as session:
            result = session.run(query, chromosome=chromosome)
            for record in result:
                total = record["total"]
    
    return total

#This function will get all genomes present in the pangenome graph
def get_genomes():
    if get_driver() is None :
        return []
    with get_driver() as driver:
    
        query = """
        MATCH (s:Stats) 
        RETURN s.genomes as all_genomes
        """
        all_genomes = []
        with driver.session() as session:
            result = session.run(query)
            for record in result:
                all_genomes = record["all_genomes"]
    
    return all_genomes

#This function get a sequence from a list of nodes names list
def get_sequence_from_names(names):
    if get_driver() is None :
        return None
    if len(names) == 0 :
        return {}
    else:
        with get_driver() as driver:
        
            query = """
            MATCH (s:Sequence) 
            WHERE s.name IN $names
            RETURN s.name as name, s.sequence as sequence
            """
            
            with driver.session() as session:
                result = session.run(query, names=names)
                return {record["name"]: record["sequence"] for record in result}
    

#This function get a sequence from a start - end / chromosome position for a given genome
def get_sequence_from_position(genome, chromosome, start, end):
    if get_driver() is None :
        return None
    sequence = ""
    if genome is None or genome == "" or chromosome is None or chromosome == "" or start is None or end is None :
        return None
    else:
        with get_driver() as driver:
            position_key = genome+"_position"
            query =  f"""MATCH (n:Node)
             WHERE n.chromosome = "{chromosome}"
               AND n.`{position_key}` >= {start}
               AND n.`{position_key}` <= {end}
               return n.ref_node AS name, 
               coalesce(n.strandM, []) AS strandMList,
               "{genome}" IN coalesce(n.strandM, []) AS strandM
               order by n.`{position_key}` ASC
            """
            #logger.info(query)
            with driver.session() as session:
                result = session.run(query)
                sorted_names = []
                sorted_strandM = []
                for record in result:
                    sorted_names.append(record["name"])
                    sorted_strandM.append(record["strandM"])
            sequences = get_sequence_from_names(sorted_names)
            for i in range(len(sorted_names)):
                if sorted_strandM[i]:
                    sequence += Seq(sequences[sorted_names[i]]).reverse_complement()
                else:
                    sequence += sequences[sorted_names[i]]
    return sequence
    


def analyse_to_csv(analyse, output_file):
    with open(output_file, mode='w', newline='', encoding='utf-8') as file:
        genome_ref = list(analyse.keys())[0]
        file.write(genome_ref+"\n")
        fields = analyse[genome_ref][0]
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(analyse[genome_ref])


#Function get_shared_regions : this function is used to get shared regions (nodes) between a list of genomes (GWAS)
#genomes_list : list of the genomes for which the function will look for shared regions
#genome_ref will be used to get annotations on this genome
#chromosomes : if defined the function will only look for shared region on theses chromosomes
#node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
#nodes_max_gap : this gap i sused to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
#deletion : if True the function will look for nodes where no one of the genome set is present
#region_trim : the shared region will be expanded in order to visualise a small region before and after. Set to 0 if strict shared regions are desired.
def get_shared_regions(genomes_list, genome_ref=None, chromosomes=None, node_min_size = 10, nodes_max_gap = 100, deletion=False, region_trim = 1000, min_percent_selected_genomes=0, tolerance_percentage = 10):
    dic_regions, analyse = find_shared_regions(genomes_list, genome_ref, chromosomes, node_min_size, nodes_max_gap, deletion = deletion, min_percent_selected_genomes=min_percent_selected_genomes, tolerance_percentage = tolerance_percentage)
    shared_regions_dict = {}
    annotations_by_regions = {}
    if genome_ref in dic_regions :
        g = genome_ref
    else :
        g = list(dic_regions.keys())[0]
    for c,items in dic_regions[g].items():
        for r in items["regions"]:
            shared_regions_dict[c+"-"+str(r["start"]-region_trim) + "-"+str(r["stop"]+region_trim)], return_code = get_nodes_by_region(g, c, r["start"]-region_trim, r["stop"]+region_trim)
            
    return shared_regions_dict, analyse
            

def find_first_ref_node_node(genome, genome_ref, genome_position, type_search = "before", chromosome="1"):
    if get_driver() is None :
        return None
    with get_driver() as driver:
        if type_search == "before" :

            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}" and n.`"""+str(genome)+"""_position` <= $genome_position AND $genome_ref in n.genomes
            return max(n.`"""+str(genome_ref)+"""_position`) as ref_position
            """
        else:

            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}" and n.`"""+str(genome)+"""_position` >= $genome_position AND $genome_ref in n.genomes
            return min(n.`"""+str(genome_ref)+"""_position`) as ref_position
            """
        #logger.info("Query : " + query)
        with driver.session() as session:
            result = session.run(query, genome_ref=genome_ref, genome_position=genome_position)
            for record in result :
                ref_position = record["ref_position"]
    
    return ref_position

#Function find_shared_regions : this function is used to get shared regions (positions) between a list of genomes (GWAS)
#It can be limited to a chromosomes list
#Usage exemple (cattle white spot) : dic_regions, analyse = find_shared_regions(["HER","SIM"],chromosome="6")
#genomes_list : list of the genomes for which the function will look for shared regions
#genome_ref will be used to get annotations on this genome
#chromosomes : list of chromosomes. If defined the function will only look for shared region on these chromosomes
#node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
#nodes_max_gap : this gap i sused to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
def find_shared_regions(genomes_list, genome_ref=None, chromosomes=None,
                        node_min_size = 10, node_max_size = 0, nodes_max_gap = 10000,
                        deletion=False, min_percent_selected_genomes=100, tolerance_percentage = 0,
                        min_deletion_percentage=100):
    if get_driver is None :
        return {},{}
    dic_regions = {}
    time_0 = time.time()
    if min_percent_selected_genomes > 100:
        min_percent_selected_genomes = 100
    if tolerance_percentage > 100:
        tolerance_percentage = 100
    if min_deletion_percentage > 100:
        min_deletion_percentage = 100
    logger.debug("node_min_size : " + str(node_min_size) + " node_max_size : " + str(node_max_size) + " deletion : " + str(deletion) + " min_percent_selected_genomes : "+ str(min_percent_selected_genomes) + " tolerance_percentage : " + str(tolerance_percentage) + " min deletion percentage : " + str(min_deletion_percentage))
    temps_depart = time.time()
    if (len(genomes_list) > 1):
        logger.info("finding shared regions for " + str(genomes_list))
        with get_driver() as driver:
            with driver.session() as session:
                query = """
                MATCH (s:Stats)
                RETURN s.genomes AS genomes
                LIMIT 1
                """
                result = session.run(query)
                for record in result:
                    genomes = record["genomes"]
                nb_genomes = len(genomes)
                nb_regions_total = 0
                nb_associated_genomes = len(genomes_list)
                
                #max_flow = nb_associated_genomes / nb_genomes + 0.00000001
                #min_flow = (max_flow-0.00000002) * min_percent_selected_genomes / 100  
                
                min_associated_genomes = max(int(min_percent_selected_genomes*nb_associated_genomes/100), 1)
                min_flow = min_associated_genomes/nb_genomes - 0.00000001
                max_flow = nb_associated_genomes*(1+tolerance_percentage/100)/nb_genomes + 0.00000001
                

                logger.debug(f"genomes number : {nb_genomes} - min flow : {min_flow} - max flow : {max_flow} - min associated genomes : {min_associated_genomes}")
                if deletion:
                        min_unselected_genomes = max(1,int((nb_genomes - nb_associated_genomes) * min_deletion_percentage / 100))
                        global_min_flow_deletion = min(min_associated_genomes + min_unselected_genomes,nb_genomes)/nb_genomes - 0.00000001 
                        min_flow_deletion = min(1,((nb_genomes - nb_associated_genomes) * min_deletion_percentage / 100)/nb_genomes - 0.00000001)
                        max_flow_deletion = min(1, (nb_genomes - nb_associated_genomes)/nb_genomes) + 0.00000001
                        logger.debug(f"Look for deletions with parameters : global min flow : {global_min_flow_deletion} - min flow : {min_flow_deletion} - max flow :  {max_flow_deletion}")
                
                if chromosomes != None :
                    chromosome_list = chromosomes
                else :
                    chromosome_list = get_chromosomes()

                if genome_ref is None or genome_ref == "":
                    genome_position_ref = genomes_list[0]
                    genome_ref = genomes_list[0]
                else :
                    genome_position_ref = genome_ref
                
                logger.debug(f"ref genome : {genome_ref}")
                # For each chromosome find shared nodes or shared deletion specific nodes
                for c in chromosome_list :
                    # First step : find specifi nodes.
                    # Specific nodes are nodes with all (or the defined proportion) of selected haplotypes
                    # And none of the selected haplotypes apart from the defined tolerance
                    logger.debug(f"chromosome : {c}")
                    dic_regions[c] = {}
                    #Looking for shared nodes
                    if node_max_size == 0:
                        query = f"""
                            MATCH (n:Node)
                            WHERE n.chromosome = '{c}'
                              AND n.flow >= {min_flow}
                              AND n.flow <= {max_flow}
                              AND n.size >= {node_min_size}
                            WITH n, [g IN n.genomes WHERE g IN $genomes_list] AS matched_genomes
                            WHERE size(matched_genomes) >= {min_associated_genomes}
                              AND size(n.genomes) - size(matched_genomes) <= size(n.genomes) * {tolerance_percentage}/100
                            OPTIONAL MATCH (n)-[]->(a:Annotation)
                            WITH n,
                                 [ann IN collect(
                                        DISTINCT CASE WHEN a IS NOT NULL THEN {{
                                            gene_id: a.gene_id,
                                            gene_name: a.gene_name,
                                            feature: a.feature
                                        }} END
                                      )
                                  WHERE ann IS NOT NULL] AS annotations
                            ORDER BY n.`{genome_position_ref}_position` ASC
                            RETURN n AS nodes, annotations;
                        """
                    else:
                        query = f"""
                            MATCH (n:Node)
                            WHERE n.chromosome = '{c}'
                              AND n.flow >= {min_flow}
                              AND n.flow <= {max_flow}
                              AND n.size >= {node_min_size}
                              AND n.size <= {node_max_size}
                            WITH n, [g IN n.genomes WHERE g IN $genomes_list] AS matched_genomes
                            WHERE size(matched_genomes) >= {min_associated_genomes}
                              AND size(n.genomes) - size(matched_genomes) <= size(n.genomes) * {tolerance_percentage}/100
                            OPTIONAL MATCH (n)-[]->(a:Annotation)
                            WITH n,
                                 [ann IN collect(
                                        DISTINCT CASE WHEN a IS NOT NULL THEN {{
                                            gene_id: a.gene_id,
                                            gene_name: a.gene_name,
                                            feature: a.feature
                                        }} END
                                      )
                                  WHERE ann IS NOT NULL] AS annotations
                            ORDER BY n.`{genome_position_ref}_position` ASC
                            RETURN n AS nodes, annotations;
                        """



                    #logger.debug(query)
                    dic_number_to_position = {}  
                    dic_number_to_size = {}
                    for g in genomes_list :
                        dic_regions[c][g] = {"nodes_position_list":[], "annotations":[], "size":[], "regions" : [], "shared_size":[], "deleted_size":[]}
                        dic_number_to_position[g] = {}
                        dic_number_to_size[g] = {}
                        #query += ' AND "' + str(g) + '" IN n.genomes'

                    #logger.debug(query)
                    # Step 2 : find shared deletion. A shared deletion if a node with
                    # all non-selected haplotypes (or the defined proportion) and
                    # none of the selected haplotypes.
                    result1 = list(session.run(query, genomes_list=genomes_list))
                    logger.debug("Nodes selected for chromosomes " + str(c) + " : " + str(len(result1)) + "\nTime : " + str(time.time()-time_0))
                    if deletion :
                        query = f"""
                            MATCH (n:Node)
                            WHERE n.chromosome = "{c}"
                              AND n.flow >= {min_flow_deletion} AND n.flow <= {max_flow_deletion}
                              AND n.size >= {node_min_size}
                              """
                        if node_max_size > 0:
                             query += f" AND n.size <= {node_max_size}"
                        
                        query += f"""
                                 
                              AND NONE(g IN $genomes_list WHERE g IN n.genomes)
                            
                            OPTIONAL CALL {{
                              WITH n
                              MATCH (m:Node)-[]->(n)
                              WHERE m.chromosome = "{c}"
                                AND m.flow >= {global_min_flow_deletion}
                            
                              WITH m, n,
                                   [g IN $genomes_list WHERE g IN m.genomes AND NOT g IN n.genomes] AS added_genomes
                              WHERE size(added_genomes) >= {min_associated_genomes}
                                AND ALL(g IN n.genomes WHERE g IN m.genomes)
                            
                              WITH m, n
                              WHERE COUNT {{ (m)-[]->(:Node) }} = 2
                            
                              MATCH (m)-[]->(n2:Node)
                              WHERE n2 <> n
                                AND n2.chromosome = "{c}"
                                AND ALL(g IN n2.genomes WHERE g IN m.genomes)
                                AND size(n2.genomes) = size(m.genomes)
                            
                            
                              WITH m, collect(n2) AS nodes_tmp, count(n2) AS nodes_tmp_count
                              WHERE nodes_tmp_count = 1
                              RETURN m, nodes_tmp[0] AS n2
                              
                            }}
                            WITH m, n, n2
                            WHERE m IS NOT NULL and n2 IS NOT NULL 
                            
                            OPTIONAL MATCH (n)-[]->(an:Annotation)
                            OPTIONAL MATCH (m)-[]->(am:Annotation)
                            OPTIONAL MATCH (n2)-[]->(an2:Annotation)
                            
                            WITH m, n, n2,
                             [a IN collect(DISTINCT CASE WHEN an IS NOT NULL THEN {{
                                    gene_id: an.gene_id,
                                    gene_name: an.gene_name,
                                    feature: an.feature
                                }} END) WHERE a IS NOT NULL] AS annotations_n,
                             [a IN collect(DISTINCT CASE WHEN am IS NOT NULL THEN {{
                                    gene_id: am.gene_id,
                                    gene_name: am.gene_name,
                                    feature: am.feature
                                }} END) WHERE a IS NOT NULL] AS annotations_m,
                             [a IN collect(DISTINCT CASE WHEN an2 IS NOT NULL THEN {{
                                    gene_id: an2.gene_id,
                                    gene_name: an2.gene_name,
                                    feature: an2.feature
                                }} END) WHERE a IS NOT NULL] AS annotations_n2
                        WITH m, n, n2,
                             annotations_n + annotations_m + annotations_n2 AS annotations

                            
                            ORDER BY m.`{genome_position_ref}_position` ASC
                            
                            RETURN
                                m AS nodes,
                                n.size AS deleted_node_size,
                                n.genomes AS genomes_deleted_nodes,
                                n2 AS end_deletion_nodes,
                                annotations;
                        """
                        
                        result2 = list(session.run(query, genomes_list=genomes_list))
                        # The shared nodes and deleted nodes are concatenated
                        result = result1 + result2
                        logger.debug("Total Nodes selected for chromosome " + str(c) + " : " + str(len(result)) + "\nTime : " + str(time.time()-time_0))
                    else:
                        result = result1
                    nb_regions_total += len(result)
                    for r in result:
                        for g in genomes_list:
                            if "deleted_nodes" not in dic_regions[c][g]:
                                dic_regions[c][g]["deleted_nodes"]= []
                            if (r["nodes"][g+"_position"] != None):
                                dic_regions[c][g]["nodes_position_list"].append(r["nodes"][g+"_position"])
                                dic_regions[c][g]["annotations"].append(r["annotations"])

                                if "deleted_node_size" in dict(r): 
                                    #If deleted nodes are found, we try to reconstruct the size of the deleted regions. 
                                    #To do this, we check that there is no overlap.
                                    deleted_nodes_dict = {}
                                    for hap in genomes : 
                                        if hap in r["genomes_deleted_nodes"]:
                                            start_deletion = r["nodes"][hap+"_position"]+r["nodes"]["size"]
                                            end_deletion = r["end_deletion_nodes"][hap+"_position"]
                                            if start_deletion > end_deletion:
                                                end_deletion_tmp = end_deletion
                                                end_deletion = start_deletion
                                                start_deletion = end_deletion_tmp
                                            deleted_nodes_dict[hap]={"start_deletion":start_deletion, "end_deletion":end_deletion}
                                        else:
                                            deleted_nodes_dict[hap] = {"start_deletion":-1, "end_deletion":-1}
                                    dic_regions[c][g]["deleted_nodes"].append(deleted_nodes_dict)
                                    #dic_regions[c][g]["deleted_size"].append(statistics.median(gap))
                                    dic_regions[c][g]["size"].append(0)
                                                                    
                                else:
                                    dic_regions[c][g]["size"].append(r["nodes"]["size"]) 
                                    #dic_regions[c][g]["deleted_size"].append(0)
                                    dic_regions[c][g]["deleted_nodes"].append({})

                    #Group regions if they are separated vy less than nodes_max_gap

                    for g in genomes_list :
                        if len(dic_regions[c][g]["nodes_position_list"]) > 0:
                            combined = list(zip(dic_regions[c][g]["nodes_position_list"], dic_regions[c][g]["annotations"],dic_regions[c][g]["size"],dic_regions[c][g]["deleted_nodes"]))
                            combined_sorted = sorted(combined, key=lambda x: x[0])
                            nodes_position_sorted, annotations_sorted, size_sorted, deleted_nodes_sorted = zip(*combined_sorted)
                            #dic_regions[c][g]["nodes_position_list"].sort()
                            shared_size = 0
                            shared_deleted_size = 0
                            current_deletion = None
                            for i in range(len(nodes_position_sorted)):
                                if i == 0 :
                                    region_start = nodes_position_sorted[0]
                                    region_stop = region_start + size_sorted[0]
                                    shared_size = size_sorted[0]
                                    annotations=annotations_sorted[0]
                                    if len(deleted_nodes_sorted[0]) > 0:
                                        current_deletion = {}
                                        for dg in deleted_nodes_sorted[0] :
                                            current_deletion[dg] = {"start_deletion":deleted_nodes_sorted[0][dg]["start_deletion"], "end_deletion":deleted_nodes_sorted[0][dg]["end_deletion"]}
                                    else:
                                        current_deletion = None

                                else :
                                    if nodes_position_sorted[i] < nodes_position_sorted[i-1] + size_sorted[i-1] + nodes_max_gap :
                                        region_stop = nodes_position_sorted[i] + size_sorted[i]
                                        shared_size += size_sorted[i]
                                        annotations += annotations_sorted[i]
                                        if len(deleted_nodes_sorted[i]) > 0:
                                            #If deleted nodes are found, we try to reconstruct the size of the deleted regions. 
                                            #To do this, we check that there is no overlap.
                                            same_deletion = False
                                            if current_deletion is None:
                                                same_deletion = True
                                            else:
                                                for dg in deleted_nodes_sorted[i] :
                                                    if "start_deletion" in deleted_nodes_sorted[i][dg] and deleted_nodes_sorted[i][dg]["start_deletion"]>=0 \
                                                        and current_deletion[dg]["start_deletion"] >= 0 \
                                                        and ((deleted_nodes_sorted[i][dg]["start_deletion"] >= current_deletion[dg]["start_deletion"] \
                                                        and deleted_nodes_sorted[i][dg]["start_deletion"] <= current_deletion[dg]["end_deletion"]) \
                                                        or (deleted_nodes_sorted[i][dg]["end_deletion"] >= current_deletion[dg]["start_deletion"] \
                                                        and deleted_nodes_sorted[i][dg]["end_deletion"] <= current_deletion[dg]["end_deletion"])) :
                                                            same_deletion = True
                                            if same_deletion :
                                                for dg in deleted_nodes_sorted[i] :
                                                    if current_deletion is not None :
                                                        if dg in current_deletion :
                                                            current_deletion[dg] = {"start_deletion": min(deleted_nodes_sorted[i][dg]["start_deletion"], current_deletion[dg]["start_deletion"]),
                                                                               "end_deletion": max(deleted_nodes_sorted[i][dg]["end_deletion"], current_deletion[dg]["end_deletion"])}
                                                        else:
                                                            current_deletion[dg] = {"start_deletion": deleted_nodes_sorted[i][dg]["start_deletion"],
                                                                               "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                                    else:
                                                        current_deletion = {}
                                                        current_deletion[dg] = {"start_deletion": deleted_nodes_sorted[i][dg]["start_deletion"],
                                                                           "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                            else:
                                                gap = []
                                                for dg in current_deletion :
                                                    if "start_deletion" in current_deletion[dg] and current_deletion[dg]["start_deletion"]>=0:
                                                        gap.append(abs(current_deletion[dg]["end_deletion"]-current_deletion[dg]["start_deletion"]))
                                                shared_deleted_size += statistics.median(gap)
                                                for dg in deleted_nodes_sorted[i] :
                                                    current_deletion[dg] = {"start_deletion":deleted_nodes_sorted[i][dg]["start_deletion"], "end_deletion":deleted_nodes_sorted[i][dg]["end_deletion"]}
                                        
                                    else :
                                        if current_deletion is not None:
                                            gap = []
                                            for dg in current_deletion :
                                                if "start_deletion" in current_deletion[dg] and current_deletion[dg]["start_deletion"]>=0:
                                                    gap.append(abs(current_deletion[dg]["end_deletion"]-current_deletion[dg]["start_deletion"]))
                                            shared_deleted_size += statistics.median(gap)
                                        if shared_deleted_size > 0 and region_stop < region_start + shared_deleted_size:
                                            region_stop = region_start + shared_deleted_size
                                        
                                        if region_start == region_stop:
                                            if shared_deleted_size > 0:
                                                region_stop = region_start + shared_deleted_size
                                            else: 
                                                region_stop += 100
                                                region_start -= 100
                                        #Minimal region to allow visualization
                                        if region_stop - region_start < 200 :
                                            min_size_region = 200
                                            gap = int((min_size_region-(region_stop - region_start))/2)
                                            region_start = max(0, region_start-gap)
                                            region_stop = region_stop+gap
                                        dic_regions[c][g]["regions"].append({"start" : region_start, "stop" : region_stop,
                                                                             "shared_size" : shared_size, "shared_deleted_size":shared_deleted_size,
                                                                             "region_size" : region_stop-region_start,
                                                                            "annotations" : annotations})
                                        shared_size = size_sorted[i]
                                        shared_deleted_size = 0
                                        annotations = annotations_sorted[i]
                                        if len(deleted_nodes_sorted[i]) > 0:
                                            for dg in deleted_nodes_sorted[i] :
                                                if current_deletion is None :
                                                    current_deletion = {}
                                                current_deletion[dg] = {"start_deletion":deleted_nodes_sorted[i][dg]["start_deletion"],
                                                                        "end_deletion":deleted_nodes_sorted[i][dg]["end_deletion"]}
                                        else:
                                            current_deletion = None
                                        region_start = nodes_position_sorted[i]
                                        region_stop = region_start + size_sorted[i]
                            if current_deletion is not None:
                                gap = []
                                for dg in current_deletion :
                                    if "start_deletion" in current_deletion[dg] and current_deletion[dg]["start_deletion"]>=0:
                                        gap.append(abs(current_deletion[dg]["end_deletion"]-current_deletion[dg]["start_deletion"]))
                                shared_deleted_size += statistics.median(gap)
                            
                            if shared_deleted_size > 0 and region_stop < region_start + shared_deleted_size:
                                region_stop = region_start + shared_deleted_size
                            if region_start == region_stop:
                                if shared_deleted_size > 0:
                                    region_stop = region_start + shared_deleted_size
                                else: 
                                    region_stop += 100
                                    region_start -= 100   
                            #Minimal region to allow visualization
                            if region_stop - region_start < 200 :
                                min_size_region = 200
                                gap = int((min_size_region-(region_stop - region_start))/2)
                                region_start = max(0, region_start-gap)
                                region_stop = region_stop+gap
                            dic_regions[c][g]["regions"].append({"start" : region_start, "stop" : region_stop,
                                                                 "shared_size" : shared_size, "shared_deleted_size":shared_deleted_size,
                                                                 "region_size" : region_stop-region_start, "annotations":annotations})
            nb_regions = 0
            for c in chromosome_list:
                nb_regions += len(dic_regions[c][genomes_list[0]]["regions"])
            logger.debug(f"Total number of identified nodes : {nb_regions_total} - Total regions nb : {nb_regions}")
            dic_regions_2 = {}
            for c, genomes in dic_regions.items():
                for genome, valeur in genomes.items():
                    if genome not in dic_regions_2:
                        dic_regions_2[genome] = {}
                    dic_regions_2[genome][c] = valeur

            analyse = {}
            logger.debug(f"genomes : {list(dic_regions_2.keys())}")
            total = sum(len(dic_regions_2[g][c]['regions']) for g in dic_regions_2 for c in dic_regions_2[g])
            with tqdm(total=total) as pbar:
                for g in dic_regions_2:
                    logger.debug(f"genome : {g}")
                    analyse[g] = []
                    for c in dic_regions_2[g]:
                        #logger.debug(f"genome : {g} - chromosome : {c} - regions number : {len(dic_regions_2[g][c]['regions'])}")
                        for r in dic_regions_2[g][c]['regions']:
                            r["chromosome"] = c
                            r["genome"] = g
                            if (genome_ref is not None and g == genome_ref) or (genome_ref is None and g == genomes_list[0]) :
                                #logger.debug("Search annotations")
                                #r["annotations"] = get_annotations_in_position_range(genome_ref=g,chromosome=c, start_position=r["start"],end_position=r["stop"])
                                #logger.debug("Search annotation before")
                                annot_before_tmp = get_annotation_before_or_after_position(genome_ref=g, chromosome=c, position=r["start"], before=True)
                                annot_tmp = {}
                                if annot_before_tmp is not None and "gene_name" in annot_before_tmp :
                                    annot_tmp["gene_name"] = annot_before_tmp["gene_name"]
                                    annot_tmp["distance"] = annot_before_tmp["end"]-r["start"]
                                r["annotation_before"] = annot_tmp

                                #logger.debug("Search annotation after")
                                annot_after_tmp = get_annotation_before_or_after_position(genome_ref=g, chromosome=c, position=r["stop"], before=False)
                                annot_tmp = {}
                                if annot_after_tmp is not None and "gene_name" in annot_after_tmp :
                                    annot_tmp["gene_name"] = annot_after_tmp["gene_name"]
                                    annot_tmp["distance"] = annot_after_tmp["start"]-r["stop"]
                                r["annotation_after"] = annot_tmp
                            analyse[g].append(r)
                            pbar.update(1)
            if genome_ref not in analyse:
                analyse[genome_ref] = []
                for a in tqdm(analyse[genomes_list[0]]):
                    r = {}
                    r["chromosome"] = a["chromosome"]
                    r["shared_size"] = a["shared_size"]
                    r["shared_deleted_size"] = a["shared_deleted_size"]
                    r["genome"] = genome_ref
                    n_start, core_genome_start = get_anchor(genomes_list[0], a["chromosome"], a["start"], before = True)
                    n_stop, core_genome_start = get_anchor(genomes_list[0], a["chromosome"], a["stop"], before = False)
                    position_field = genome_ref+"_position"
                    r["start"] = 0
                    if position_field in n_start :
                        r["start"] = n_start[position_field]
                    if position_field in n_stop:
                        r["stop"] = n_stop[position_field]
                    else:
                        r["stop"] = r["start"]
                    r["region_size"] = r["stop"] - r["start"]
                    r["annotations"] = get_annotations_in_position_range(genome_ref=genome_ref,chromosome=a["chromosome"], start_position=r["start"],end_position=r["stop"])
                    annot_before_tmp = get_annotation_before_or_after_position(genome_ref=genome_ref, chromosome=a["chromosome"], position=r["start"], before=True)
                    annot_tmp = {}
                    if annot_before_tmp is not None and "gene_name" in annot_before_tmp :
                        annot_tmp["gene_name"] = annot_before_tmp["gene_name"]
                        annot_tmp["distance"] = annot_before_tmp["end"]-r["start"]
                    r["annotation_before"] = annot_tmp
                    
                    annot_after_tmp = get_annotation_before_or_after_position(genome_ref=genome_ref, chromosome=a["chromosome"], position=r["stop"], before=False)
                    annot_tmp = {}
                    if annot_after_tmp is not None and "gene_name" in annot_after_tmp :
                        annot_tmp["gene_name"] = annot_after_tmp["gene_name"]
                        annot_tmp["distance"] = annot_after_tmp["start"]-r["stop"]
                    r["annotation_after"] = annot_tmp
                    analyse[genome_ref].append(r)
                        
            for g in analyse:
                analyse[g] = sorted(analyse[g], key=lambda d: d['shared_size'], reverse=True)
            
            logger.debug("Total time : "+ str(time.time()-temps_depart))
    return dic_regions_2, analyse



def calculer_variabilite(chromosome_list=None, ref_genome=None, window_size=1000, output_html="pangenome_variability.html"):
    with get_driver() as driver:
        if ref_genome == None :
            ref_noeud = "node_mean"
            ref_position = "position_mean"
        else :
            ref_noeud = f"{ref_genome}_noeud"
            ref_position = f"{ref_genome}_position"
        # Définir le renderer pour ouvrir dans le navigateur
        pio.renderers.default = 'browser'
        
        html_parts = []
        if chromosome_list is None:
            chromosomes = get_chromosomes()
        else:
            chromosomes = chromosome_list
  
        for chromosome in chromosomes:
            logger.info("Compute variability on chromosome " + str(chromosome))
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
            WITH FLOOR(n.{ref_position} / {window_size}) AS Window, n.{ref_position} AS position, n.flow AS flow
            WITH Window, 
                 collect(position) AS positions, 
                 avg(flow) AS flow_mean, 
                 count(*) AS nodes_number
            RETURN 
                Window,
                reduce(min_pos = head(positions), p IN positions | CASE WHEN p < min_pos THEN p ELSE min_pos END) AS start_position,
                reduce(max_pos = head(positions), p IN positions | CASE WHEN p > max_pos THEN p ELSE max_pos END) AS end_position,
                flow_mean,
                nodes_number
            ORDER BY Window
            """
        
            with driver.session() as session:
                result = session.run(query)
                data = result.data()
        
            df = pd.DataFrame(data)
            if df.empty:
                logger.warning(f"No data found for chromosome {chromosome}.")
                continue
        
            # Add window size and hover text
            df['window_size'] = df['end_position'] - df['start_position']
            df['hover_text'] = (
                "Window: " + df['Window'].astype(str) +
                "<br>Start: " + df['start_position'].astype(str) +
                "<br>End: " + df['end_position'].astype(str) +
                "<br>Size: " + df['window_size'].astype(str) +
                "<br>Flow mean: " + df['flow_mean'].round(2).astype(str)
            )
        
            # Create figure
            fig = px.scatter(
                df,
                x='Window',
                y='flow_mean',
                hover_name='hover_text',
                labels={
                    'Window': 'Window',
                    'flow_mean': 'Flow mean'
                },
                title=f"Flow variability - Chromosome {chromosome} ({ref_genome}) - window size {window_size}",
                height=500
            )
        
            fig.update_traces(marker=dict(size=10, color='dodgerblue', line=dict(width=1, color='darkblue')))
            fig.update_layout(hovermode='closest')
        
            # Generate HTML of the figure and store it
            html_parts.append(pio.to_html(fig, include_plotlyjs='cdn', full_html=False))
        
        if not html_parts:
            logger.warning("No graph computed.")
            return
        
        # Combine all graphs into a single HTML file
        with open(output_html, 'w') as f:
            f.write("<html><head><title>Flow variability</title></head><body>\n")
            for part in html_parts:
                f.write(part + "<hr>\n")
            f.write("</body></html>")
        
        logger.info(f"Graphs stored into {output_html}")
        file_path = os.path.abspath(output_html)
        webbrowser.open(f'file://{file_path}')

    
#This function takes nodes data (from get_nodes_by_region or get_nodes_by_feature for exemple)
#it computes the Jaccard distance on these nodes
#and it returns the distance matrix and the distance matrix weighted by the nodes size
def compute_phylo_tree_from_nodes(nodes_data,output_dir = "", weighted=False):
    if get_driver() is None :
        return None
    genomes = get_genomes()
    nodes_nb = len(nodes_data)
    genome_redondant_strand_matrix = np.zeros((len(genomes),2 * nodes_nb), dtype=np.int32)
    genome_strand_matrix = np.zeros((len(genomes),2 * nodes_nb), dtype=np.int32)
    
    i = 0
    index_genomes = {}
    genomes_names = []
    for g in genomes:
        index_genomes[g] = i
        genomes_names.append(g)
        i += 1
    i = 0

    for n in nodes_data:
        for g in genomes :
            if "strandM" in nodes_data[n] and g in nodes_data[n]["strandM"]:
                genome_redondant_strand_matrix[index_genomes[g], i] += 1
                genome_strand_matrix[index_genomes[g], i] = nodes_data[n]["size"]
            if "strandP" in nodes_data[n] and g in nodes_data[n]["strandP"]:
                genome_redondant_strand_matrix[index_genomes[g], i+nodes_nb] += 1
                genome_strand_matrix[index_genomes[g], i+nodes_nb] = nodes_data[n]["size"]
        i += 1
    #computes Jaccard distance on matrix
    jaccard_matrix = np.zeros((len(genomes), len(genomes)))
    weighted_jaccard_matrix = np.zeros((len(genomes), len(genomes)))
    for i in range(len(genomes)):
        for j in range(i, len(genomes)):
            min_sum = np.minimum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :]).sum()
            max_sum = np.maximum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :]).sum()
            
            min_counts = np.minimum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :])
            max_counts = np.maximum(genome_redondant_strand_matrix[i, :], genome_redondant_strand_matrix[j, :])
            
            size_matrix = np.maximum(genome_strand_matrix[i, :], genome_strand_matrix[j, :])
            # Weight by nodes size
            weighted_min = (min_counts * size_matrix).sum()
            weighted_max = (max_counts * size_matrix).sum()
            
            if max_sum == 0:
                jaccard_index = 0.0
            else:
                jaccard_index = min_sum / max_sum
                
            if weighted_max == 0:
                weighted_jaccard_index = 0.0
            else:
                weighted_jaccard_index = weighted_min / weighted_max

            jaccard_matrix[i, j] = 1-jaccard_index
            jaccard_matrix[j, i] = 1-jaccard_index
            
            weighted_jaccard_matrix[i,j] = 1-weighted_jaccard_index
            weighted_jaccard_matrix[j,i] = 1-weighted_jaccard_index
            
    df_jaccard = pd.DataFrame(jaccard_matrix, index=genomes_names, columns=genomes_names)
    df_weighted_jaccard = pd.DataFrame(weighted_jaccard_matrix, index=genomes_names, columns=genomes_names)
    if output_dir != "" and os.path.isdir(output_dir) :
        df_jaccard.to_csv(output_dir+'/distance_matrix.csv')
        df_weighted_jaccard.to_csv(output_dir+'/weighted_distance_matrix.csv')
    
    if weighted:
        df_val = df_weighted_jaccard.values
    else:
        df_val = df_jaccard.values
    triangulaire_inf = []
    for i in range(0, df_val.shape[0]):
        triangulaire_inf.append([])
        for j in range(i + 1):
            triangulaire_inf[i].append(df_val[i, j])

    dm = Phylo.TreeConstruction.DistanceMatrix(list(df_jaccard.columns), triangulaire_inf)
    # computing tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dm)

    newick_tree = tree.format('newick')
    #logger.debug(newick_tree)
    return newick_tree
            

"""
Fonction to convert PAV matrix to phylip format

Parameters
    @pav_matrix_dic : PAV matrix (dictionnary)
    @matrice_distance_phylip_filename : output filename (PHYLIP)
returns
"""


def pav_to_phylip(pav_matrix_dic, distance_matrix_phylip_filename):
    # Creating phylip file
    file = open(distance_matrix_phylip_filename, "w")

    with file:
        if len(list(pav_matrix_dic.keys())) > 0:
            entete = str(len(list(pav_matrix_dic.keys()))) + " " + str(
                len(pav_matrix_dic[list(pav_matrix_dic.keys())[0]])) + "\n"
            file.write(entete)
            for item in pav_matrix_dic:
                ligne = str(item) + " "
                for p in pav_matrix_dic[item]:
                    ligne += str(p)
                ligne += "\n"
                file.write(ligne)

    file.close()




#This function compute a global tree from a random selection of nodes taking account of direct / reverse traversing (if strnd is True)
#Parameters :
# - Method :
#   - raxml : this will launch raxml on the sample matrix
#   - nj : this will compute a distance matrix from the sample matrix and then compute a neighbor joining tree
#       this method is less accurate but far much faster than raxml and it can be usefull for big pangenomes
# - strand : is true used the traversal direction, if false it doesn't take it inot account
# - chromosome : to limit the tree to 1 chromosome, set chromosome value (else it will be computed on the whole graph)
# - project_name is used to name the RaxML files
# - max_nodes is used to limit the number of sampled nodes
# - min_sample_size is the minimal sampled size for pangenome of size >= min_sample_size
# - min_nodes_number : tree won't be computed on pangenomes with less than min_nodes_number nodes
def compute_global_phylo_tree_from_nodes(method="raxml", output_dir = "", strand=True, chromosome = None, project_name="panabyss_phylo_tree", max_nodes = 1000000, min_sample_size = 10000, min_nodes_number = 1000):
    total_nodes_number = 0
    

    dir_raxml = "./export/phylo/raxml"
    dir_phylo = "./export/phylo"
    distance_matrix_filename = "distance_matrix.phy"
    tree_newick_filename = os.path.join(dir_raxml,project_name+".raxml.bestTree")
    last_tree = "./export/phylo/last_tree.nwk"

    # Copier le fichier
    distance_matrix_phylip_filename = os.path.join(dir_phylo, distance_matrix_filename)
    if not os.path.exists(dir_raxml):
        os.makedirs(dir_raxml)
    if get_driver() is None :
        return None
    
    with get_driver() as driver:
        if driver is None:
            return None
        if chromosome is None :
            query = f"""
            MATCH (n:Node)
            return count(*) as total_nodes_number
            """
        else:
            logger.debug(f"Chromosome : {chromosome}")
            query = f"""
            MATCH (n:Node)
            where n.chromosome = '{chromosome}'
            return count(*) as total_nodes_number
            """
 
        with driver.session() as session:
            result = session.run(query)
            for record in result :
                total_nodes_number = record["total_nodes_number"]
                
        #The number of sampled nodes depends on the total number of nodes. It is set to be at least min_sample_size and at most max_nodes.
        #If the pangenome contains less than min_nodes_number then no tree can be computed.
        #The number of sampled nodes is determined by a polynomial interpolation in log-space of the total number of nodes, fitted through a set of control points.
        if total_nodes_number < min_nodes_number:
            return None
        if total_nodes_number < min_sample_size:
            sample_nodes_number = total_nodes_number
        else:
            #Set of control points :
            x_vals = np.array([1e4, 1e5, 1e7, 1e8, 2e9])
            y_vals = np.array([1e4, 3e4, 9e4, 1e5, 2e5])
            logx = np.log10(x_vals)
            #Compute the polynomial interpolation in log space:
            coeffs = np.polyfit(logx, y_vals, deg=4)
            def f(x):
                return np.polyval(coeffs, np.log10(x))
            sample_nodes_number = max(min_sample_size, min(int(f(total_nodes_number)),max_nodes))
        #sample_nodes_number = max(min_sample_size,min(int(node_selection_percentage*total_nodes_number/100), max_nodes))
        logger.info(f"Number of nodes to sample : {sample_nodes_number} - Total node {total_nodes_number}")
        if chromosome is None :
            query = f"""
            MATCH (n:Node)
            with n, rand() AS r
            order by r
            limit {sample_nodes_number}
            return distinct(n) as nodes
            """
        else:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = '{chromosome}'
            with n, rand() AS r
            order by r
            limit {sample_nodes_number}
            return distinct(n) as nodes
            """
        nodes_list = []
        with driver.session() as session:
            result = session.run(query)
            for record in result :
                nodes_list.append(dict(record["nodes"]))
        
        
        logger.info(f"Number of sampled nodes : {len(nodes_list)}")
        sample_size = len(nodes_list)
        if sample_size >= min_sample_size :
            #Prepare PAV matrix
            genomes = get_genomes()
            pav_matrix = {}
            for g in genomes : 
                if g not in pav_matrix:
                    if strand :
                        pav_matrix[g] = np.zeros(2 * sample_size, dtype=int)
                    else:
                        pav_matrix[g] = np.zeros(sample_size, dtype=int)
            
            for i in range(0, len(nodes_list)):
                for g in nodes_list[i]["genomes"]:
                    if strand:
                        if "strandP" in nodes_list[i] and g in nodes_list[i]["strandP"]:
                            pav_matrix[g][i] = int(1)
                        else:
                            pav_matrix[g][i+sample_size] = int(1)
                    else:
                        pav_matrix[g][i] = int(1)         
            pav_to_phylip(pav_matrix, distance_matrix_phylip_filename)
            if method == "raxml":
                logger.debug("RaxML method...")
                pattern = os.path.join(dir_raxml, f"{project_name}.*")
                for f in glob.glob(pattern):
                    os.remove(f)
                # raxml_command = [
                #     "raxmlHPC",
                #     "-s", f"../{distance_matrix_filename}",
                #     "-m", "BINGAMMA",
                #     "-p", "12345",
                #     # '-#', '100',  # Iterations number for bootstrapping
                #     "-n", project_name
                # ]

                raxml_command = [
                    "raxml-ng",
                    "--msa", f"../{distance_matrix_filename}",
                    "--model", "BIN+G",
                    "--seed", "12345",
                    "--prefix", project_name,
                ]

                #
                # iqtree_command = [
                #     "iqtree2",
                #     "-s", f"../{distance_matrix_filename}",
                #     "-seed", "12345",
                #     "-st",  "MORPH",
                #     "-pre", project_name
                # ]


                # launching RaxML command
                result = subprocess.run(raxml_command, check=True, cwd=dir_raxml)
                try:
                    with open(tree_newick_filename, 'r') as f:
                        shutil.copy(tree_newick_filename, last_tree)
                        return f.read()
                except FileNotFoundError:
                    return None
            else:
                #method with matrix distance and neighbor joining
                logger.debug("Neighbor joining method...")
                names = list(pav_matrix.keys())
                matrix = np.array([pav_matrix[name] for name in names])
                n = len(names)
                lower_tri = []
                for i in range(n):
                    lower_tri.append([])
                    for j in range(i+1):
                        # Hamming distance
                        lower_tri[i].append(np.sum(matrix[i] != matrix[j]) / matrix.shape[1])
                dm = _DistanceMatrix(names, lower_tri)
                constructor = DistanceTreeConstructor()
                tree = constructor.nj(dm)
                newick_tree = tree.format('newick')
                with open(last_tree, 'w') as f:
                    f.write(newick_tree)
                return newick_tree

        else:
            return None


            

#Computes the distance matrix on the whole GFA (could take a long time for big GFA)
def compute_distance_matrix(distance_matrix_filename = "distances.csv", chromosome=None, ponderation=True, strand=False):
    if get_driver() is None :
        return None
    temps_depart = time.time()
    with get_driver() as driver:
        with driver.session() as session:
            #Get genomes list
            query = """
            MATCH (s:Stats)
            RETURN s.genomes AS genomes
            LIMIT 1
            """
            result = session.run(query)
            for record in result:
                genomes = record["genomes"]

            distance_matrix = pd.DataFrame(data=0,index=genomes, columns=genomes)
            dic_size_genome = {}   
            if ponderation :
                #Get the genome size
                query = """
                    MATCH (n:Node)
                    UNWIND n.genomes AS g
                    RETURN g AS genome, sum(n.size) AS total_size
                    """
                result_size_genomes = list(session.run(query))
                if strand == False:
                    #computes size for each genome

                    #Get intersection size
                    query = """
                        MATCH (n:Node)
                        UNWIND n.genomes AS g1
                        UNWIND n.genomes AS g2
                        with g1,g2,n
                        WHERE g1 < g2 
                        return g1, g2, sum(n.size) AS size_intersection
                        ORDER BY size_intersection DESC
                        """

                else:
                    query = """
                        MATCH (n:Node)
                        WITH n.size AS size, n.strandM AS strandM, n.strandP AS strandP
                        WITH size,
                             [x IN range(0, size(strandM)-2) | [strandM[x], strandM[x+1..]]] +
                             [x IN range(0, size(strandP)-2) | [strandP[x], strandP[x+1..]]] AS pairGroups
                        UNWIND pairGroups AS group
                        UNWIND group[1] AS g2
                        WITH group[0] AS g1, g2, size
                        WHERE g1 < g2
                        RETURN g1, g2, sum(size) AS size_intersection
                        ORDER BY size_intersection DESC
                        """
                result_intersection = list(session.run(query))
                    
            else:
                #Get the node numbers
                query = """
                    MATCH (n:Node)
                    UNWIND n.genomes AS g
                    RETURN g AS genome, count(*) AS total_size
                """
                result_size_genomes = list(session.run(query)) 
                if strand == False:

                    #Get intersection size
                    query = """
                        MATCH (n:Node)
                        UNWIND n.genomes AS g1
                        UNWIND n.genomes AS g2
                        with g1,g2,n
                        WHERE g1 < g2 
                        return g1, g2, count(*) AS size_intersection
                        ORDER BY size_intersection DESC
                        """
                else:
                    query = """
                        MATCH (n:Node)
                        WITH n.size AS size, n.strandM AS strandM, n.strandP AS strandP
                        WITH size,
                             [x IN range(0, size(strandM)-2) | [strandM[x], strandM[x+1..]]] +
                             [x IN range(0, size(strandP)-2) | [strandP[x], strandP[x+1..]]] AS pairGroups
                        UNWIND pairGroups AS group
                        UNWIND group[1] AS g2
                        WITH group[0] AS g1, g2, size
                        WHERE g1 < g2
                        RETURN g1, g2, count(*) AS size_intersection
                        ORDER BY size_intersection DESC
                        """
                result_intersection = list(session.run(query))
            
        for r in result_size_genomes:
            dic_size_genome[r["genome"]] = r["total_size"]
        logger.debug(dic_size_genome)
        for r in result_intersection:
            g1 = r["g1"]
            g2 = r["g2"]
            inter = r["size_intersection"]
            distance_matrix.loc[g1,g2] = 1-inter/(dic_size_genome[g1]+dic_size_genome[g2]-inter)
            distance_matrix.loc[g2,g1] = distance_matrix.loc[g1,g2]
                
                
        distance_matrix.to_csv(distance_matrix_filename)   
        logger.debug("Total time : "+ str(time.time()-temps_depart))
    return distance_matrix



          
 