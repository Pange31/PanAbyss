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
import csv
import json
import subprocess
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, _DistanceMatrix
from scipy.spatial.distance import pdist, squareform
from scipy.stats import hypergeom, binom, norm, chi2_contingency, fisher_exact
import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Seq import Seq
import statistics
from config import *
from neo4j_driver import get_driver,  get_scoped_driver
from concurrent.futures import ThreadPoolExecutor, as_completed
import uuid
from sqlite_gwas_requests import *
from sqlite_phylo_requests import *
import hashlib


EXECUTOR = ThreadPoolExecutor(max_workers=4)

logger = logging.getLogger("panabyss_logger")

# This value is used to limit the nodes number when seeking for regions :
# If the region is wider than this value then it is ignored
MAX_BP_SEEKING = 800000

# Maximal number of nodes to get the whole region
MAX_NODES_NUMBER = get_max_nodes_from_db()

logging.getLogger("neo4j").setLevel(logging.ERROR)
logger.debug(f"Max nodes number from DB: {MAX_NODES_NUMBER}")

MAX_GWAS_STORE, MAX_RUNNING_INACTIVITY_HOURS, MAX_GWAS_REGIONS, GWAS_ANNOTATIONS_WINDOWS_SIZE, GWAS_ANNOTATIONS_MAX_ATTEMPTS = get_gwas_conf()


logger.debug(f"GWAS parameters : max_store = {MAX_GWAS_STORE}, "
             f"max_running_inactivity_hours = {MAX_RUNNING_INACTIVITY_HOURS}, max_gwas_regions = {MAX_GWAS_REGIONS}")

# This function allows to execute parallel queries in neo4j
def run_queries_parallel(queries_dict, max_threads=10, job_id=None):
    driver = get_scoped_driver()
    if driver is None:
        return None

    results = {}
    total_queries = len(queries_dict)
    completed = 0

    def task(key, query):
        session = driver.session()
        try:
            result = session.run(query)
            data = [r.data() for r in result]
            return key, data
        except Exception as e:
            logger.warning(f"Neo4j error {key}: {e}")
            return key, None
        finally:
            session.close()

    try:
        with ThreadPoolExecutor(max_workers=max_threads) as executor:
            futures = {
                executor.submit(task, key, query): key
                for key, query in queries_dict.items()
            }

            for future in as_completed(futures):
                key = futures[future]

                try:
                    _, result = future.result()
                    results[key] = result

                    if result is not None:
                        logger.debug(f"Chromosome {key} ({len(result)} results)")

                except Exception as e:
                    logger.debug(f"Error for chromosome {key} : {e}")
                    results[key] = None
                    continue

                if job_id:
                    completed += 1

                    if get_gwas_status(job_id) == "CANCEL":
                        logger.debug(f"Cancel job id : {job_id}")
                        return {}

                    progress = round((completed / total_queries) * 75, 2)
                    update_gwas_progression(job_id, progress)

        return results

    finally:
        try:
            driver.close()
            logger.info("🧹 Scoped Neo4j driver closed")
        except Exception as e:
            logger.warning(f"Error closing driver: {e}")

# This function get the first node containing genome_ref before or after a node for a given genome
def get_genome_position(genome_ref, genome, chromosome, position, before=True, max_attempts=200):
    core_genome = False
    driver = get_driver()
    if driver is None:
        return None
    genome_position = genome + "_position"
    genome_ref_position = genome_ref + "_position"
    window_size = 1000
    attempt = 0
    lower_bound = position
    upper_bound = position
    with driver.session() as session:
        while attempt < max_attempts and lower_bound > 0:
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
              AND n.`{genome_ref_position}` > 0
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
                return dict(record["n"]), core_genome
        # Anchor not found => the current node will be used
        logger.debug("No node on the reference genome found, get the current nodes.")
        if before:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` <= $position
            RETURN n order by n.`{genome_position}` DESC limit 1
            """
        else:
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
        if record is None:
            if before:
                query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` >= $position
                    RETURN n order by n.`{genome_position}` ASC limit 1
                """
            else:
                query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` <= $position
                    RETURN n order by n.`{genome_position}` DESC limit 1
                """
            result = session.run(
                query,
                chromosome=chromosome,
                position=position
            )
            record = result.single()
        if record:
            return dict(record["n"]), core_genome

    return None, core_genome

#This function try to find a core genome node befor or after a position on a reference genome
#If no core if found, then it will return the nearer position on the reference genome
def get_anchor(genome, chromosome, position, before=True, use_anchor=True):
    core_genome = False
    driver = get_driver()
    if driver is None:
        return None
    with driver.session() as session:
        genome_position = genome + "_position"
        if use_anchor:
            window_size = 1000
            max_attemps = 500
            attempt = 0
            lower_bound = int(position)
            upper_bound = int(position)

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
                    return dict(record["n"]), core_genome
        if use_anchor:
            # Anchor not found => the current node will be used
            logger.debug("No core genome anchor found, get the current nodes.")
        #Get the nearer nodes on the reference genome
        if before:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` <= {position}
            RETURN n order by n.`{genome_position}` DESC limit 1
            """
        else:
            query = f"""
            MATCH (n:Node)
            WHERE n.chromosome = "{chromosome}"
              AND n.`{genome_position}` >= {position}
            RETURN n order by n.`{genome_position}` ASC limit 1
            """
        result = session.run(
            query
        )
        record = result.single()
        if record is None:
            if before:
                query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` >= {position}
                    RETURN n order by n.`{genome_position}` ASC limit 1
                """
            else:
                query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = "{chromosome}"
                      AND n.`{genome_position}` <= {position}
                    RETURN n order by n.`{genome_position}` DESC limit 1
                """
            result = session.run(
                query
            )
            record = result.single()
        if record:
            return dict(record["n"]), core_genome

    return None, core_genome


# Internal request to get all the nodes into a range defined in ranges dictionnary (containing start and end)
# min_node_size will be used is set to get only nodes with a size greater than this value
# flow will be used is set to get only nodes with a flow greater than this value
# valid_individuals_exceptions : if not None this list is used to filter haplotypes : haplotypes not in this list will be used
def construct_base_query(ranges, chromosome, min_node_size=None, flow=None, valid_individuals_exceptions=None):
    subqueries = []

    for g in ranges:
        if valid_individuals_exceptions is None or g not in valid_individuals_exceptions:
            position_field = g + "_position"
            subquery = f"""
                MATCH (m:Node)
                WHERE m.chromosome = "{chromosome}"
                """

            if min_node_size is not None and min_node_size > 1:
                subquery += f" AND m.size >= {min_node_size}"
            if flow is not None:
                subquery += f" AND m.flow >= {flow}"
            subquery += f"""
              AND m.`{position_field}` >= {ranges[g]['start']}
              AND m.`{position_field}` <= {ranges[g]['stop']}
            RETURN m
            """

            subqueries.append(subquery)

    # Assemblage final avec UNION ALL
    base_query_genome = f"""
        CALL {{
            {" UNION ALL ".join(subqueries)}
        }}
        """
    return base_query_genome


#This function take the result of neo4 request and get the different annotations
#It returns the nodes_data
def get_nodes_data_from_record(result):
    nodes_data = {}
    for record in result:

        node_name = record["m"]["name"]

        #Transcripts
        transcripts = []
        seen_transcripts = set()
        for transcript in record["transcripts"]:
            transcript_id = transcript["transcript_id"]
            if transcript_id is None:
                continue

            if transcript_id in seen_transcripts:
                continue

            seen_transcripts.add(transcript_id)

            transcripts.append({
                "transcript_id": transcript_id,
                "start": transcript["start"],
                "end": transcript["end"]
            })

        #Exons
        exons_map = {}

        for exon in record["exons"]:
            exon_id = exon["exon_id"]

            if exon_id is None:
                continue

            transcript_id = exon["transcript_id"]

            if exon_id not in exons_map:
                exons_map[exon_id] = {
                    "exon_id": exon_id,
                    "start": exon["start"],
                    "end": exon["end"],
                    "transcript_ids": set()
                }

            if transcript_id is not None:
                exons_map[exon_id]["transcript_ids"].add(transcript_id)

        # conversion
        exons = []
        for exon_id, exon_data in exons_map.items():
            exon_data["transcript_ids"] = list(exon_data["transcript_ids"])
            exons.append(exon_data)

        # nodes_data

        nodes_data[node_name] = (
                dict(record["m"])
                | {
                    "sequence": record["sequence"]
                }
                | {
                    "genes_names": list(set(
                        a for a in record["annotations"]
                        if a is not None
                    ))
                }
                | {
                    "features": list(set(
                        f for f in record["features"]
                        if f is not None
                    ))
                }
                | {
                    "transcripts": transcripts,
                    "exons": exons
                }
        )

    return nodes_data

# This function take a region (chromosome, start and stop) of a given haplotype (search_genome)
# and it returns all the nodes in this region and the other related regions :
# for each haplotype the start and stop are given by the first anchor node before the start position and the first after the end position
# Anchor node is a node with all haplotypes
# return metadata :
# return_code :
#   - OK => nodes found
#   - PARTIAL => no anchor on pangenome found
#   - WIDE => more than MAX_NODES_NUMBER nodes found into the region
#   - NO_DATA => no data found
#   - ZOOM => more than MAX_NODES_NUMBER nodes found but with filtering with flow an acceptable number of nodes has been found
#       in this case the flow value will be set
#   - FILTER => exceptional individuals have been removed from the search (if not there would be too much nodes)
# flow : in case of filtering, the minimal flow used to filter data
def get_nodes_by_region(genome, chromosome, start, end, use_anchor=True, min_node_size=None, max_nodes_number=None):
    return_metadata = {"return_code": "OK", "flow": None, "nodes_number": 0, "removed_genomes": None}
    valid_individuals_exceptions = []
    flow = None
    ranges = {}
    driver = get_driver()
    if driver is None:
        return None
    temps_depart = time.time()
    nodes_data = {}
    max_sequence = 1000
    if end is None:
        stop = MAX_BP_SEEKING
    else:
        stop = end
    LIMIT = MAX_NODES_NUMBER
    if max_nodes_number:
        LIMIT = max_nodes_number
    data = {}
    shared_genomes = []


    # query_annotations = f"""
    #                 WITH DISTINCT m
    #                 OPTIONAL MATCH (m)-[]->(a:Annotation)
    #                 OPTIONAL MATCH (s:Sequence {{name: m.ref_node}})
    #                 RETURN m, substring(s.sequence, 0, {max_sequence}) as sequence, collect(a.gene_name) AS annotations, collect(a.feature) AS features
    #                 LIMIT {LIMIT + 1}
    #                 """

    query_annotations = f"""
                        WITH DISTINCT m
                        OPTIONAL MATCH (m)-[]->(a:Annotation)
                        OPTIONAL MATCH (m)-[]->(t:Annotation)
                        WHERE t.feature IN ["mrna", "transcript"]
                        OPTIONAL MATCH (m)-[]->(e:Annotation)
                        WHERE e.feature = "exon"
                        OPTIONAL MATCH (s:Sequence {{name: m.ref_node}})
                        
                        WITH m, s,
                             collect(DISTINCT a.gene_name) AS annotations,
                             collect(DISTINCT a.feature) AS features,
                             collect(DISTINCT {{
                                 transcript_id: t.transcript_id,
                                 start: t.start,
                                 end: t.end
                             }}) AS transcripts,
                             collect(DISTINCT {{
                                 transcript_id: e.transcript_id,
                                 exon_id: e.exon_id,
                                 start: e.start,
                                 end: e.end
                             }}) AS exons
                        
                        RETURN
                            m,
                            substring(s.sequence, 0, {max_sequence}) AS sequence,
                            annotations,
                            features,
                            transcripts,
                            exons
                        LIMIT {LIMIT + 1}
                        """

    with driver.session() as session:

        # Step 1 : find the anchors of the region
        if start is not None and end is not None:
            genome_position = genome + "_position"
            logger.info("Look for region : " + str(start) + " - " + str(stop) + " - chromosome " + str(
                chromosome) + " - genome : " + str(genome))

            anchor_start, core_genome_start = get_anchor(genome, chromosome, start, before=True,
                                                         use_anchor=use_anchor)
            anchor_stop, core_genome_stop = get_anchor(genome, chromosome, end, before=False, use_anchor=use_anchor)
            if anchor_start is None or anchor_stop is None:
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
                    limit=LIMIT + 1
                )
                result = session.run(query_genome, start=start, end=end)
                record = result.single()
                if record:
                    pairs = record["genome_ranges"]
                    ranges = {genome: data for genome, data in pairs}
            else:
                ranges = {}
                for g in anchor_start["genomes"] + anchor_stop["genomes"]:
                    p_start = anchor_start[g + "_position"]
                    p_stop = anchor_stop[g + "_position"]
                    ranges[g] = {"start": min(p_start, p_stop), "stop": max(p_start, p_stop)}

            if anchor_start is not None and anchor_stop is not None:
                if anchor_start[genome_position] > anchor_stop[genome_position]:
                    anchor_start_tmp = anchor_start
                    anchor_start = anchor_stop
                    anchor_stop = anchor_start_tmp
                logger.debug("Anchor start name : " + str(anchor_start["name"]))
                logger.debug("Anchor stop name : " + str(anchor_stop["name"]))
                logger.debug("Anchor region : " + str(anchor_start[genome_position]) + " - " + str(
                    anchor_stop[genome_position]))

            # Step 2 : find the region between the 2 anchors
            if anchor_start is not None and anchor_stop is not None and anchor_stop[genome_position] - anchor_start[
                genome_position] > 0 and len(anchor_start['genomes']) > 0:
                region_nodes_number = 0
                # construct the base query to find all genomes between the start / stop position

                query_genome = construct_base_query(ranges, chromosome, min_node_size=min_node_size,
                                                    flow=None) + f"""
                    WITH DISTINCT m
                    LIMIT {LIMIT + 1}
                    RETURN count(m) AS nodes_number
                    """

                result = session.run(query_genome)
                record = result.single()
                if record:
                    region_nodes_number = int(record["nodes_number"])
                if region_nodes_number <= LIMIT:
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
                        """
                        if min_node_size is not None and min_node_size > 1:
                            q += f" AND m.size >= {min_node_size}"

                        q += f"""
                              AND m.{position_field} >= {ranges[g]['start']} AND m.{position_field} <= {ranges[g]['stop']}
                            WITH m LIMIT {LIMIT + 1}
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
                        if counts[g] > LIMIT or counts[g] > 10 * median_value:
                            individuals_exceptions.append(g)
                    logger.debug(f"Exceptional individuals : {individuals_exceptions}")
                    if len(individuals_exceptions) == 1 or len(individuals_exceptions) <= 0.2 * len(counts):
                        valid_individuals_exceptions = individuals_exceptions
                        logger.debug(
                            f"These individuals will be removed from search : {valid_individuals_exceptions}")
                    else:
                        # Step 6 : the region is too wide and it is not linked to a small proportion of individudls
                        # => try to find a "zoom level" for which the number of nodes is acceptable
                        # To do that it will use the flow attribute
                        logger.debug("Too much nodes into the region, try to reduce data by filtering on flow.")
                        zoom = True
                        flow = 1
                        while zoom and flow >= 0:

                            query_genome = construct_base_query(ranges, chromosome, min_node_size=min_node_size,
                                                                flow=flow) + f"""
                                                                WITH DISTINCT m
                                                                LIMIT {LIMIT + 1}
                                                                RETURN count(m) AS nodes_number
                                                                """
                            # logger.debug(query_genome)
                            result = session.run(query_genome, start=start, end=end)
                            record = result.single()
                            if record:
                                region_nodes_number = int(record["nodes_number"])
                            if region_nodes_number > LIMIT:
                                zoom = False
                                logger.debug(f"Too much nodes with flow {flow}.")
                            else:
                                if flow - 0.1 >= 0:
                                    logger.debug(
                                        f"Nodes number found with flow {flow} : {region_nodes_number} - check with flow {flow - 0.1}")
                                flow -= 0.1
                        if not zoom and flow == 1:
                            flow = -1
                        else:
                            flow += 0.1
                        if flow >= 0:
                            logger.debug(f"Zoom level found with flow {flow}.")
                        if flow < 0:
                            logger.debug("Too much nodes into the region, no zoom level found.")

                if flow is None or (flow is not None and flow >= 0) or len(valid_individuals_exceptions) > 0:
                    # Step 7 : Get the nodes and annotations for each genomes
                    if flow is not None and flow > 0:
                        # The search will be filtered by flow
                        return_metadata["flow"] = flow
                        return_metadata["return_code"] = "ZOOM"
                        query_genome = construct_base_query(ranges, chromosome, min_node_size=min_node_size,
                                                            flow=flow)
                    elif len(valid_individuals_exceptions) > 0:
                        logger.debug("request construction")
                        query_genome = construct_base_query(ranges, chromosome, min_node_size=min_node_size,
                                                            flow=None,
                                                            valid_individuals_exceptions=valid_individuals_exceptions)
                        return_metadata["flow"]: 0
                        return_metadata["removed_genomes"] = valid_individuals_exceptions
                        return_metadata["return_code"] = "FILTER"
                    else:
                        query_genome = query_genome = construct_base_query(ranges, chromosome,
                                                                           min_node_size=min_node_size, flow=None)

                    query_genome = query_genome + query_annotations
                    # logger.debug(f"query genome : {query_genome}")
                    # logger.info(query_genome)
                    result = session.run(query_genome, start=start, end=end)

                    nodes_data = get_nodes_data_from_record(result)


                    # for record in result:
                    #     nodes_data[record["m"]["name"]] = dict(record["m"]) | {"sequence": record["sequence"]} | {
                    #         "annotations": set(
                    #             record["annotations"][a] for a in range(len(record["annotations"])))} | {
                    #                                           "features": set(record["features"][a] for a in
                    #                                                           range(len(record["features"])))}


                    return_metadata["nodes_number"] = len(nodes_data)
                    if len(nodes_data) > LIMIT:
                        nodes_data = {}
                        return_metadata["return_code"] = "WIDE"
                        logger.warning(
                            f"Region too wide : nodes number : {len(nodes_data)} - max nodes number : {LIMIT}")

                else:
                    return_metadata["return_code"] = "WIDE"
                    logger.warning(
                        f"Region too wide : {anchor_stop[genome_position] - anchor_start[genome_position]} - nodes number : {region_nodes_number} - max nodes number : {LIMIT}")
                    nodes_data = {}
            else:
                if anchor_start is not None and anchor_stop is not None and anchor_stop[genome_position] - \
                        anchor_start[genome_position] >= MAX_BP_SEEKING:
                    logger.warning(
                        f"Region too wide : {anchor_stop[genome_position] - anchor_start[genome_position]}")
                    return_metadata["return_code"] = "WIDE"
                else:
                    logger.warning("Region not found")
                    return_metadata["return_code"] = "NO_DATA"
                nodes_data = {}

        else:
            # if end is set to None, get all the nodes if the number of nodes is less than LIMIT

            if start == 0 and end is None:
                total_nodes = get_nodes_number(chromosome)
                if total_nodes <= LIMIT:
                    query_genome = f"""
                    MATCH (m:Node)
                    WHERE  m.chromosome = "{chromosome}" 
                    """
                    if min_node_size is not None and min_node_size > 1:
                        query_genome += f" AND m.size >= {min_node_size} "
                    query_genome += query_annotations
                    result = session.run(query_genome)
                    nodes_data = get_nodes_data_from_record(result)
                    # for record in result:
                    #     nodes_data[record["m"]["name"]] = dict(record["m"]) | {"sequence": record["sequence"]} | {
                    #         "annotations": set(
                    #             record["annotations"][a] for a in range(len(record["annotations"])))} | {
                    #                                           "features": set(record["features"][a] for a in
                    #                                                           range(len(record["features"])))}
                    return_metadata["nodes_number"] = len(nodes_data)
                    if len(nodes_data) > LIMIT:
                        nodes_data = {}
                        return_metadata["return_code"] = "WIDE"
                        logger.warning(
                            f"Region too wide : nodes number : {len(nodes_data)} - max nodes number : {LIMIT}")
                else:
                    logger.warning(f"Region too wide: total node {total_nodes} - limit : {LIMIT}")
                    return_metadata["return_code"] = "WIDE"
                    nodes_data = {}
        # if len(nodes_data) > 0:
        #     for elt in nodes_data:
        #         nodes_data[elt]["annotations"] = list(nodes_data[elt]["annotations"])
        #         nodes_data[elt]["features"] = list(nodes_data[elt]["features"])
    logger.debug("Total time : " + str(time.time() - temps_depart))
    return nodes_data, return_metadata


# Find a region by annotation : gene_name or gene_id
#   - OK
#   - WIDE
#   - NO_DATA
def get_nodes_by_feature(genome, chromosome, feature=None, value=None, min_node_size=None, max_nodes_number=None):
    return_code = "OK"
    driver = get_driver()
    if driver is None:
        return None

    nodes_data = {}
    shared_genomes = []
    key_id = "id"
    key_name = "name"
    if feature.lower() == "gene":
        key_id = "gene_id"
        key_name = "gene_name"
    elif feature.lower() == "transcript" or feature.lower() == "mrna":
        key_id = "transcript_id"
        key_name = "transcript_name"
    elif feature.lower() == "exon":
        key_id = "exon_id"
        key_name = "exon_id"

    with driver.session() as session:
        genome_position = genome + "_position"
        # Step 1 : find nodes with gene annotation
        if feature is not None and value is not None:
            logger.info(f"Looking for {feature} : {value}")

            query = f"""
                    MATCH (a:Annotation {{chromosome:"{chromosome}", feature:"{feature}"}})<-[]-(n:Node)
                    WHERE n[$genome_position] IS NOT NULL
                      AND (a.{key_id} = $value OR a.{key_name} = $value)
                    RETURN DISTINCT n
                    ORDER BY n[$genome_position] ASC
                    """
            result = session.run(query, value=value, genome_position=genome_position)
            # logger.debug(query)
            noeuds_annotes = [record["n"] for record in result]
            if len(noeuds_annotes) > 0 and genome_position in noeuds_annotes[0] and genome_position in \
                    noeuds_annotes[-1]:
                start = noeuds_annotes[0][genome_position]
                stop = noeuds_annotes[-1][genome_position] + noeuds_annotes[-1]["size"]
                logger.debug(f"start : {start} - stop : {stop} - nodes number : {len(noeuds_annotes)}")
                nodes_data, return_metadata = get_nodes_by_region(genome, chromosome, start, stop,
                                                                  min_node_size=min_node_size, max_nodes_number=max_nodes_number)
            else:
                logger.debug(f"No nodes found {len(noeuds_annotes)}.")
                return_metadata = {"return_code": "NO_DATA", "flow": None, "nodes_number": 0}
        else:
            return_metadata = {"return_code": "NO_DATA", "flow": None, "nodes_number": 0}

    return nodes_data, return_metadata



#This function get the first node linked to an annotation
#The windows search is limitated by the windows length and the max_attempts
def get_annotation_before_or_after_position(genome_ref, chromosome="1", position=0, before=True, driver=None):
    if driver is None:
        driver = get_driver()
        if driver is None:
            return None
    genome_position = genome_ref + "_position"
    window_size = GWAS_ANNOTATIONS_WINDOWS_SIZE
    attempt = 0
    lower_bound = position
    upper_bound = position
    with driver.session() as session:
        while attempt < GWAS_ANNOTATIONS_MAX_ATTEMPTS and lower_bound > 0:
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
                MATCH (n:Node)-[]-(a:Annotation)
                WHERE n.chromosome = "{chromosome}"
                  AND n.`{genome_position}` > $lower_bound
                  AND n.`{genome_position}` < $upper_bound
                  AND a.gene_name IS NOT NULL AND a.feature = "gene"
                RETURN a.gene_name AS gene_name,
                       a.end AS end,
                       a.start AS start,
                       a.feature AS feature,
                       n.`{genome_position}` AS neighbor_node_position,
                       n.size AS neighbor_node_size
                ORDER BY n.`{genome_position}` {order},
                         abs(a.end - n.`{genome_position}`) ASC
                LIMIT 1
            """


            result = session.run(
                query,
                lower_bound=lower_bound,
                upper_bound=upper_bound
            )
            record = result.single()
            if record:
                return dict(record)
    return None



# This function get first annotation on nodes of a reference_genome on a given chromosome after the given position
# def get_annotation_before_or_after_position(genome_ref, chromosome="1", position=0, before=True):
#     if get_driver() is None:
#         return {}
#     with get_driver() as driver:
#         if before:
#             query = f"""
#             MATCH (a:Annotation)
#             WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.end < $position and a.gene_name is not null
#             RETURN a.gene_name AS gene_name, a.end as end, a.start as start, a.feature as feature
#             ORDER BY a.end DESC
#             LIMIT 1
#             """
#         else:
#             query = f"""
#             MATCH (a:Annotation)
#             WHERE a.genome_ref = "{genome_ref}" and a.chromosome = "{chromosome}" and a.start > $position and a.gene_name is not null
#             RETURN a.gene_name AS gene_name, a.end as end, a.start as start, a.feature as feature
#             ORDER BY a.start ASC
#             LIMIT 1
#             """
#
#         with driver.session() as session:
#             result = session.run(query, position=position)
#             record = result.single()
#             return dict(record) if record else None


# This function get all annotations on nodes of a reference_genome on a given chromosome between start_position and end_position
def get_annotations_in_position_range(genome_ref, chromosome="1", start_position=0, end_position=0):
    driver = get_driver()
    if driver is None:
        return None
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
    # logger.info("Query : " + query)
    with driver.session() as session:
        result = session.run(query, start=start_position, end=end_position, genome_ref=genome_ref)
        annotations = [dict(record) for record in result]

    return annotations


# This function get all features on annotation nodes
def get_annotations_features():
    driver = get_scoped_driver()
    if driver is None:
        return []
    query = f"""
    MATCH (a:Annotation)
    RETURN DISTINCT a.feature as feature
    """
    with driver.session() as session:
        result = session.run(query)
        features = [record["feature"] for record in result]

    return features


# This function will get all chromosomes present in the pangenome graph
def get_chromosomes():
    driver = get_scoped_driver()
    if driver is None:
        return []

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


# This function will get all chromosomes present in the pangenome graph
def get_chromosomes_stats():
    chromosome_stats = None
    driver = get_scoped_driver()
    if driver is None:
        return []

    query = """
    MATCH (cs:chromosome_stats) 
    RETURN cs as cs
    """
    with driver.session() as session:
        result = session.run(query)
        for record in result:
            chromosome_stats = dict(record["cs"])

    return chromosome_stats


# This function will get the number of nodes in the graph
def get_nodes_number(chromosome=None):
    driver = get_driver()
    if driver is None:
        return 0
    total = 0
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


# This function will get all genomes present in the pangenome graph
def get_genomes():
    driver = get_scoped_driver()
    if driver is None:
        return []
    query = """
    MATCH (s:Stats) 
    RETURN s.genomes as all_genomes
    """
    all_genomes = []
    with driver.session() as session:
        result = session.run(query)
        for record in result:
            all_genomes = record["all_genomes"]

    return sorted(all_genomes)


# This function get a sequence from a list of nodes names list
def get_sequence_from_names(names):
    driver = get_driver()
    if driver is None:
        return None
    if len(names) == 0:
        return {}
    else:
        query = """
        MATCH (s:Sequence) 
        WHERE s.name IN $names
        RETURN s.name as name, s.sequence as sequence
        """

        with driver.session() as session:
            result = session.run(query, names=names)
            return {record["name"]: record["sequence"] for record in result}


# This function get a sequence from a start - end / chromosome position for a given genome
def get_sequence_from_position(genome, chromosome, start, end):
    driver = get_driver()
    if driver is None:
        return None
    sequence = ""
    if genome is None or genome == "" or chromosome is None or chromosome == "" or start is None or end is None:
        return None
    else:
        position_key = genome + "_position"
        query = f"""MATCH (n:Node)
         WHERE n.chromosome = "{chromosome}"
           AND n.`{position_key}` >= {start}
           AND n.`{position_key}` <= {end}
           return n.ref_node AS name, 
           coalesce(n.strandM, []) AS strandMList,
           "{genome}" IN coalesce(n.strandM, []) AS strandM
           order by n.`{position_key}` ASC
        """
        # logger.info(query)
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
        file.write(genome_ref + "\n")
        fields = analyse[genome_ref][0]
        writer = csv.DictWriter(file, fieldnames=fields)
        writer.writeheader()
        writer.writerows(analyse[genome_ref])


# Function get_shared_regions : this function is used to get shared regions (nodes) between a list of genomes (GWAS)
# genomes_list : list of the genomes for which the function will look for shared regions
# genome_ref will be used to get annotations on this genome
# chromosomes : if defined the function will only look for shared region on theses chromosomes
# node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
# nodes_max_gap : this gap i sused to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
# deletion : if True the function will look for nodes where no one of the genome set is present
# region_trim : the shared region will be expanded in order to visualise a small region before and after. Set to 0 if strict shared regions are desired.
def get_shared_regions(genomes_list, all_genomes, genome_ref=None, chromosomes=None, node_min_size=10, nodes_max_gap=100,
                       deletion=False, region_trim=1000, min_percent_selected_genomes=0, tolerance_percentage=10):
    dic_regions, analyse, dic_distribution = find_shared_regions(genomes_list, all_genomes, genome_ref, chromosomes, node_min_size,
                                                                 nodes_max_gap, deletion=deletion,
                                                                 min_percent_selected_genomes=min_percent_selected_genomes,
                                                                 tolerance_percentage=tolerance_percentage)
    shared_regions_dict = {}
    annotations_by_regions = {}
    if genome_ref in dic_regions:
        g = genome_ref
    else:
        g = list(dic_regions.keys())[0]
    for c, items in dic_regions[g].items():
        for r in items["regions"]:
            shared_regions_dict[c + "-" + str(r["start"] - region_trim) + "-" + str(
                r["stop"] + region_trim)], return_code = get_nodes_by_region(g, c, r["start"] - region_trim,
                                                                             r["stop"] + region_trim)

    return shared_regions_dict, analyse


def find_first_ref_node_node(genome, genome_ref, genome_position, type_search="before", chromosome="1"):
    driver = get_driver()
    if driver is None:
        return None
    if type_search == "before":

        query = f"""
        MATCH (n:Node)
        WHERE n.chromosome = "{chromosome}" and n.`""" + str(genome) + """_position` <= $genome_position AND $genome_ref in n.genomes
        return max(n.`""" + str(genome_ref) + """_position`) as ref_position
        """
    else:

        query = f"""
        MATCH (n:Node)
        WHERE n.chromosome = "{chromosome}" and n.`""" + str(genome) + """_position` >= $genome_position AND $genome_ref in n.genomes
        return min(n.`""" + str(genome_ref) + """_position`) as ref_position
        """
    # logger.info("Query : " + query)
    with driver.session() as session:
        result = session.run(query, genome_ref=genome_ref, genome_position=genome_position)
        for record in result:
            ref_position = record["ref_position"]

    return ref_position



#Function used to hash params : if the params already exist in database that will avoid to launch the process
#and this will return the results directly
def compute_gwas_params_hash(params):
    #Sort list of genomes
    if "genomes_list" in params and isinstance(params["genomes_list"], list):
        params["genomes_list"].sort()

    normalized = json.dumps(params, sort_keys=True)
    return hashlib.sha256(normalized.encode()).hexdigest()

#This function will launch the find_shared_regions in a thread
#If a job with the same parameters already exists and is SUCCESS,
#return the existing job_id directly.
#Otherwise, create a new job and launch the background process.
def submit_job_gwas(params):
    check_running_gwas()
    params_hash = compute_gwas_params_hash(params)
    existing_job = find_existing_gwas_job(params_hash)

    if existing_job:
        job_id = existing_job["job_id"]
        logger.debug(f"Job already exist: {job_id}")
    else:
        # New job => create and launch
        job_id = str(uuid.uuid4())
        insert_gwas_job(job_id, params, params_hash)
        logger.debug(f"Launching new job {job_id}")
    EXECUTOR.submit(_run_job_gwas, job_id, params)
    return job_id

#Background job to find shared regions
#If SUCCESS already exists for this job_id => return the results from the database.
#If RUNNING => do nothing, let the FO poll for updates.
#Otherwise (new job or ERROR) => run the find_shared_regions process and update the database.
def _run_job_gwas(job_id, params):
    try:
        #Retrieve the current job from the database from params values
        job_data = get_gwas_job(job_id)
        if not job_data:
            logger.error(f"Job {job_id} not found in the database")
            return None

        status = job_data["status"]
        if status == "SUCCESS":
            #Result already computed for this job_id
            logger.debug(f"Job {job_id} already SUCCESS, retrieving results from the database")
            #Update timestamp
            update_gwas_job_timestamp(job_id)
            return job_data["result_gwas_regions"], job_data["result_gwas_points"]

        if status == "RUNNING":
            logger.debug(f"Job {job_id} already running, skipping execution")
            return None

        #Mark the job as RUNNING
        set_gwas_job_running(job_id)

        #Launch the search region process
        logger.debug("Calling find_shared_regions...")
        dic_region, analyse, dic_distribution  = find_shared_regions(**params, job_id=job_id,
                                                                     results_only_for_ref=True, max_gwas_region=MAX_GWAS_REGIONS)

        #Store results
        if get_gwas_status(job_id) == "RUNNING":
            set_gwas_job_success(job_id, analyse, dic_distribution)

        #Return results for the current job
        return analyse, dic_distribution

    except Exception as e:
        set_gwas_job_error(job_id, str(e))
        logger.exception(f"Error while running job {job_id}")
        return None

#ACAT P-value combination method
def compute_acat(pvals, weights=None):

    pvals = np.asarray(pvals)

    #pvals = np.clip(pvals, 1e-300, 1-1e-16)

    if weights is None:
        weights = np.ones(len(pvals)) / len(pvals)
    else:
        weights = np.asarray(weights)
        weights = weights / np.sum(weights)

    t = np.sum(
        weights * np.tan((0.5 - pvals) * np.pi)
    )

    p_acat = 0.5 - np.arctan(t) / np.pi

    return p_acat

# Function find_shared_regions : this function is used to get shared regions (positions) between a list of genomes (GWAS)
# It can be limited to a chromosomes list
# Usage exemple (cattle white spot) : dic_regions, analyse = find_shared_regions(["HER","SIM"],chromosome="6")
# genomes_list : list of the genomes for which the function will look for shared regions
# genome_ref will be used to get annotations on this genome
# chromosomes : list of chromosomes. If defined the function will only look for shared region on these chromosomes
# node_min_size : the nodes smaller than this value will be ignored (to avoid to look for all snp, if the are required then set this value to 0)
# nodes_max_gap : this gap is used to gather find regions into a bigger regions if the initial find regions are separated by less than this value (in numer of nodes)
def find_shared_regions(genomes_list, all_genomes, genome_ref=None, chromosomes=None,
                        node_min_size=10, node_max_size=0, nodes_max_gap=10000,
                        deletion=False, min_percent_selected_genomes=100, tolerance_percentage=0,
                        min_deletion_percentage=100, job_id=None, results_only_for_ref=False, max_gwas_region=None):
    if get_driver is None:
        logger.debug("No driver found")
        return {}, {}, {}

    dic_regions = {}
    dic_distribution = {}
    time_0 = time.time()
    if min_percent_selected_genomes > 100:
        min_percent_selected_genomes = 100
    if tolerance_percentage > 100:
        tolerance_percentage = 100
    if min_deletion_percentage > 100:
        min_deletion_percentage = 100
    if node_min_size is None:
        node_min_size = 0
    logger.debug(
        "node_min_size : " + str(node_min_size) + " node_max_size : " + str(node_max_size) + " deletion : " + str(
            deletion) + " min_percent_selected_genomes : " + str(
            min_percent_selected_genomes) + " tolerance_percentage : " + str(
            tolerance_percentage) + " min deletion percentage : " + str(min_deletion_percentage))
    temps_depart = time.time()
    if (len(genomes_list) > 1):
        set_selected_genomes = set(genomes_list)
        set_not_selected_genomes = set(all_genomes)-set_selected_genomes
        logger.info("finding shared regions for " + str(genomes_list))
        driver = get_scoped_driver()
        if driver is None:
            return None
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

            # max_flow = nb_associated_genomes / nb_genomes + 0.00000001
            # min_flow = (max_flow-0.00000002) * min_percent_selected_genomes / 100

            min_associated_genomes = max(int(min_percent_selected_genomes * nb_associated_genomes / 100), 1)
            min_flow = min_associated_genomes / nb_genomes - 0.00000001
            max_flow = nb_associated_genomes * (1 + tolerance_percentage / 100) / nb_genomes + 0.00000001

            logger.debug(
                f"genomes number : {nb_genomes} - min flow : {min_flow} - max flow : {max_flow} - min associated genomes : {min_associated_genomes}")
            if deletion:
                min_unselected_genomes = max(1, int((
                                                                nb_genomes - nb_associated_genomes) * min_deletion_percentage / 100))
                global_min_flow_deletion = min(min_associated_genomes + min_unselected_genomes,
                                               nb_genomes) / nb_genomes - 0.00000001
                min_flow_deletion = min(1, ((
                                                        nb_genomes - nb_associated_genomes) * min_deletion_percentage / 100) / nb_genomes - 0.00000001)
                max_flow_deletion = min(1, (nb_genomes - nb_associated_genomes) / nb_genomes) + 0.00000001
                logger.debug(
                    f"Look for deletions with parameters : global min flow : {global_min_flow_deletion} - min flow : {min_flow_deletion} - max flow :  {max_flow_deletion}")

            if chromosomes != None:
                chromosome_list = chromosomes
            else:
                chromosome_list = get_chromosomes()

            if genome_ref is None or genome_ref == "":
                genome_position_ref = genomes_list[0]
                genome_ref = genomes_list[0]
            else:
                genome_position_ref = genome_ref

            logger.debug(f"ref genome : {genome_ref}")
            # For each chromosome find shared nodes or shared deletion specific nodes
            # These queries are launched in parallel because they take some time
            queries = {}
            for c in chromosome_list:
                # First step : find specific nodes.
                # Specific nodes are nodes with all (or the defined proportion) of selected haplotypes
                # And none of the selected haplotypes apart from the defined tolerance
                # logger.debug(f"chromosome : {c}")

                # prepare result structure for each chromosome

                dic_regions[c] = {}
                dic_distribution[c] = []
                if results_only_for_ref :
                    dic_regions[c][genome_ref] = {"nodes_position_list": [], "annotations": [], "size": [], "regions": [],
                                         "shared_size": [], "deleted_size": [], "position_mean": [], "pval": [],}
                else:
                    for g in genomes_list:
                        dic_regions[c][g] = {"nodes_position_list": [], "annotations": [], "size": [], "regions": [],
                                             "shared_size": [], "deleted_size": [], "position_mean": [], "pval":[]}

                # Looking for shared nodes
                genome_count_expr = " + ".join(
                    f"(CASE WHEN n.`{genome}_position` IS NOT NULL THEN 1 ELSE 0 END)"
                    for genome in genomes_list
                )
                query = f"""
                    MATCH (n:Node)
                    WHERE n.chromosome = '{c}'
                      AND n.flow >= {min_flow}
                      AND n.flow <= {max_flow}
                      AND n.size >= {node_min_size}
                      """
                if node_max_size > 0:
                    query += f" AND n.size <= {node_max_size}"
                query += f"""
                    WITH n, {genome_count_expr} AS matched_genomes_nb
                    WHERE matched_genomes_nb >= {min_associated_genomes}
                      AND size(n.genomes) - matched_genomes_nb <= size(n.genomes) * {tolerance_percentage}/100
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

                queries[c + "_shared"] = query

                # logger.debug(query)

                # Step 2 : find shared deletion. A shared deletion if a node with
                # all non-selected haplotypes (or the defined proportion) and
                # none of the selected haplotypes.
                # results :
                # n = first deleted node
                # m = node juste before deleted node (linked to n)
                # n2 = end_deletion_nodes = first node after the deletion
                # result1 = list(session.run(query, genomes_list=genomes_list))
                # logger.debug("Nodes selected for chromosomes " + str(c) + " : " + str(len(result1)) + "\nTime : " + str(time.time()-time_0))
                if deletion:
                    none_expr = " AND ".join(
                        f"n.`{genome}_position` IS NULL" for genome in genomes_list
                    )

                    query = f"""
                        MATCH (n:Node)
                        WHERE n.chromosome = "{c}"
                          AND n.flow >= {min_flow_deletion} AND n.flow <= {max_flow_deletion}
                          AND n.size >= {node_min_size}
                          """
                    if node_max_size > 0:
                        query += f" AND n.size <= {node_max_size}"

                    query += f"""

                        AND {none_expr}

                        OPTIONAL CALL {{
                          WITH n
                          MATCH (m:Node)-[]->(n)
                          WHERE m.chromosome = "{c}"
                            AND m.flow >= {global_min_flow_deletion}

                          WITH m, n,
                               [g IN {genomes_list} WHERE g IN m.genomes AND NOT g IN n.genomes] AS added_genomes
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
                    queries[c + "_deletion"] = query
            logger.debug("Search for shared regions.")
            results = run_queries_parallel(queries, job_id=job_id)
            if len(results) == 0:
                return {}, {},{}

            #Kc is the total count for each chromosome of shared and deleted regions for statistical test
            # Kc = {}
            # Kc["total"] = {"shared": 0, "deleted": 0}
            # N = 0
            nodes_stats = {}
            # chromosome_stats = get_chromosomes_stats()
            N = len(all_genomes)
            K = len(genomes_list)
            for c in chromosome_list:
                # Kc_shared is the total shared size on the chromosome
                # Kc[c] = {"shared":0, "deleted":0}
                # N += chromosome_stats.get(c+"_max_position_mean")
                if deletion:
                    result = results[c + "_shared"] + results[c + "_deletion"]
                else:
                    result = results[c + "_shared"]
                logger.debug(
                    "Total nodes selected for chromosome " + str(c) + " : " + str(len(result)) + "\nTime : " + str(
                        time.time() - time_0))
                nb_regions_total += len(result)
                #Iterate ton construct the size for each detected nodes
                for r in result:
                    set_node_genomes = set(r["nodes"]["genomes"])
                    #a and b are used for statistical test
                    #a = bnumber of selected genomes present on the node
                    #b = number of non selected genomes present on the node
                    a = len(set_node_genomes.intersection(set_selected_genomes))
                    b = len(set_node_genomes.intersection(set_not_selected_genomes))

                    if "deleted_node_size" not in dict(r):
                        table = [[a, 0],
                                 [K - a, (N - K) - b]]
                    else:
                        #For deleted nodes, the table is construct from
                        #the first node before deletion
                        #=> the number of selected genomes on the node is b = the number of
                        #not selected genomes passing through the node juste before deletion
                        table = [[b, K-a],
                                 [N - K - b, a]]
                    chi2, pval, dof, expected = chi2_contingency(table)
                    nodes_stats[r["nodes"]["name"]] = pval + 1e-300
                    for g in genomes_list:
                        #Checks if all results are required, if not, computes result only for reference genome
                        if results_only_for_ref and g != genome_ref:
                            continue
                        dic_regions[c][g]["pval"].append(pval + 1e-300)
                        if "deleted_nodes" not in dic_regions[c][g]:
                            dic_regions[c][g]["deleted_nodes"] = []
                        if (g + "_position" in r["nodes"] and r["nodes"][g + "_position"] != None):
                            dic_regions[c][g]["nodes_position_list"].append(r["nodes"][g + "_position"])
                            dic_regions[c][g]["position_mean"].append(r["nodes"]["position_mean"])
                            dic_regions[c][g]["annotations"].append(r["annotations"])
                            if "deleted_node_size" not in dict(r):
                                dic_distribution[c].append((r["nodes"]["position_mean"], r["nodes"]["size"]))
                            if "deleted_node_size" in dict(r):
                                # If deleted nodes are found, we try to reconstruct the size of the deleted regions.
                                # To do this, we compute the median of deleted nodes for each haplotype
                                deleted_nodes_dict = {}
                                gap = []
                                for hap in genomes:
                                    if hap in r["genomes_deleted_nodes"]:
                                        start_deletion = r["nodes"][hap + "_position"] + r["nodes"]["size"]
                                        end_deletion = r["end_deletion_nodes"][hap + "_position"]
                                        if start_deletion > end_deletion:
                                            end_deletion_tmp = end_deletion
                                            end_deletion = start_deletion
                                            start_deletion = end_deletion_tmp
                                        gap.append(end_deletion - start_deletion)
                                        deleted_nodes_dict[hap] = {"start_deletion": start_deletion,
                                                                   "end_deletion": end_deletion}
                                    else:
                                        deleted_nodes_dict[hap] = {"start_deletion": -1, "end_deletion": -1}
                                dic_regions[c][g]["deleted_nodes"].append(deleted_nodes_dict)
                                dic_regions[c][g]["size"].append(0)
                                if len(gap) > 0 and statistics.median(gap) is not None:
                                    dic_distribution[c].append(
                                        (r["nodes"]["position_mean"] + r["nodes"]["size"], -statistics.median(gap)))

                            else:
                                dic_regions[c][g]["size"].append(r["nodes"]["size"])
                                dic_regions[c][g]["deleted_nodes"].append({})

                # Group regions if they are separated vy less than nodes_max_gap
                p_values_list = []
                weighted_p_values_list = []
                for g in genomes_list:
                    # Checks if all results are required, if not, computes result only for reference genome
                    if results_only_for_ref and g != genome_ref:
                        continue
                    if len(dic_regions[c][g]["nodes_position_list"]) > 0:
                        combined = list(
                            zip(dic_regions[c][g]["nodes_position_list"], dic_regions[c][g]["annotations"],
                                dic_regions[c][g]["size"], dic_regions[c][g]["deleted_nodes"],
                                dic_regions[c][g]["pval"]))
                        combined_sorted = sorted(combined, key=lambda x: x[0])
                        nodes_position_sorted, annotations_sorted, size_sorted, deleted_nodes_sorted, pval_sorted = zip(
                            *combined_sorted)

                        shared_size = 0
                        shared_deleted_size = 0
                        nb_nodes_in_region = 0
                        current_deletion = None

                        sum_pval = 0
                        sum_score = 0
                        for i in range(len(nodes_position_sorted)):
                            if i == 0:
                                nb_nodes_in_region = 1
                                region_start = nodes_position_sorted[0]
                                region_stop = region_start + size_sorted[0]
                                shared_size = size_sorted[0]
                                annotations = annotations_sorted[0]

                                sum_pval = -np.log10(pval_sorted[0])
                                p_values_list.append(pval_sorted[0])
                                if len(deleted_nodes_sorted[0]) > 0:
                                    current_deletion = {}
                                    for dg in deleted_nodes_sorted[0]:
                                        current_deletion[dg] = {
                                            "start_deletion": deleted_nodes_sorted[0][dg]["start_deletion"],
                                            "end_deletion": deleted_nodes_sorted[0][dg]["end_deletion"]}
                                else:
                                    current_deletion = None
                                    sum_score = size_sorted[0] * (-np.log10(pval_sorted[0]))
                                    weighted_p_values_list = [size_sorted[0]]
                                    p_values_list = [pval_sorted[0]]

                            else:
                                if nodes_position_sorted[i] < nodes_position_sorted[i - 1] + size_sorted[
                                    i - 1] + nodes_max_gap:
                                    nb_nodes_in_region += 1
                                    region_stop = nodes_position_sorted[i] + size_sorted[i]
                                    shared_size += size_sorted[i]
                                    annotations += annotations_sorted[i]
                                    sum_pval += -np.log10(pval_sorted[i])

                                    if len(deleted_nodes_sorted[i]) > 0:
                                        # If deleted nodes are found, we try to reconstruct the size of the deleted regions.
                                        # To do this, we check that there is no overlap.
                                        same_deletion = False
                                        if current_deletion is None:
                                            same_deletion = True
                                        else:
                                            for dg in deleted_nodes_sorted[i]:
                                                if "start_deletion" in deleted_nodes_sorted[i][dg] and \
                                                        deleted_nodes_sorted[i][dg]["start_deletion"] >= 0 \
                                                        and current_deletion[dg]["start_deletion"] >= 0 \
                                                        and ((deleted_nodes_sorted[i][dg]["start_deletion"] >=
                                                              current_deletion[dg]["start_deletion"] \
                                                              and deleted_nodes_sorted[i][dg]["start_deletion"] <=
                                                              current_deletion[dg]["end_deletion"]) \
                                                             or (deleted_nodes_sorted[i][dg]["end_deletion"] >=
                                                                 current_deletion[dg]["start_deletion"] \
                                                                 and deleted_nodes_sorted[i][dg]["end_deletion"] <=
                                                                 current_deletion[dg]["end_deletion"])):
                                                    same_deletion = True
                                        if same_deletion:
                                            for dg in deleted_nodes_sorted[i]:
                                                if current_deletion is not None:
                                                    if dg in current_deletion:
                                                        current_deletion[dg] = {"start_deletion": min(
                                                            deleted_nodes_sorted[i][dg]["start_deletion"],
                                                            current_deletion[dg]["start_deletion"]),
                                                                                "end_deletion": max(
                                                                                    deleted_nodes_sorted[i][dg][
                                                                                        "end_deletion"],
                                                                                    current_deletion[dg][
                                                                                        "end_deletion"])}
                                                    else:
                                                        current_deletion[dg] = {
                                                            "start_deletion": deleted_nodes_sorted[i][dg][
                                                                "start_deletion"],
                                                            "end_deletion": deleted_nodes_sorted[i][dg][
                                                                "end_deletion"]}
                                                else:
                                                    current_deletion = {}
                                                    current_deletion[dg] = {
                                                        "start_deletion": deleted_nodes_sorted[i][dg][
                                                            "start_deletion"],
                                                        "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                        else:
                                            gap = []
                                            for dg in current_deletion:
                                                if "start_deletion" in current_deletion[dg] and \
                                                        current_deletion[dg]["start_deletion"] >= 0:
                                                    gap.append(abs(
                                                        current_deletion[dg]["end_deletion"] - current_deletion[dg][
                                                            "start_deletion"]))
                                            shared_deleted_size += statistics.median(gap)
                                            sum_score += shared_deleted_size * (-np.log10(pval_sorted[i]))
                                            weighted_p_values_list.append(shared_deleted_size)
                                            p_values_list.append(pval_sorted[i])
                                            for dg in deleted_nodes_sorted[i]:
                                                current_deletion[dg] = {
                                                    "start_deletion": deleted_nodes_sorted[i][dg]["start_deletion"],
                                                    "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                    else:
                                        #Not a deleted node
                                        sum_score += size_sorted[i] * (-np.log10(pval_sorted[i]))
                                        weighted_p_values_list.append(size_sorted[i])
                                        p_values_list.append(pval_sorted[i])
                                else:
                                    #New region
                                    if current_deletion is not None:
                                        #Compute the last deletion of previous region
                                        gap = []
                                        for dg in current_deletion:
                                            if "start_deletion" in current_deletion[dg] and current_deletion[dg][
                                                "start_deletion"] >= 0:
                                                gap.append(abs(
                                                    current_deletion[dg]["end_deletion"] - current_deletion[dg][
                                                        "start_deletion"]))
                                        shared_deleted_size += statistics.median(gap)
                                        weighted_p_values_list.append(shared_deleted_size)
                                        p_values_list.append(pval_sorted[i])
                                        sum_score += shared_deleted_size * (-np.log10(pval_sorted[i]))
                                    if shared_deleted_size > 0 and region_stop < region_start + shared_deleted_size:
                                        region_stop = region_start + shared_deleted_size

                                    if region_start == region_stop:
                                        if shared_deleted_size > 0:
                                            region_stop = region_start + shared_deleted_size
                                        else:
                                            region_stop += 100
                                            region_start -= 100
                                    # Minimal region to allow visualization
                                    if region_stop - region_start < 200:
                                        min_size_region = 200
                                        gap = int((min_size_region - (region_stop - region_start)) / 2)
                                        region_start = max(0, region_start - gap)
                                        region_stop = region_stop + gap
                                    # Kc[c]["shared"] += shared_size
                                    # Kc[c]["deleted"] += shared_deleted_size
                                    # Kc["total"]["shared"] += shared_size
                                    # Kc["total"]["deleted"] += shared_deleted_size
                                    acat = compute_acat(p_values_list, weighted_p_values_list)
                                    dic_regions[c][g]["regions"].append({"start": region_start, "stop": region_stop,
                                                                         "shared_size": shared_size,
                                                                         "shared_deleted_size": shared_deleted_size,
                                                                         "region_size": region_stop - region_start,
                                                                         "annotations": annotations,
                                                                         "score":sum_score/(max(region_stop - region_start, 1)),
                                                                         "pval": acat})
                                    #Init the data for the new region
                                    shared_size = size_sorted[i]
                                    shared_deleted_size = 0
                                    nb_nodes_in_region = 1
                                    # min_position_mean = position_mean_sorted[i]
                                    # max_position_mean = min_position_mean + size_sorted[i]
                                    sum_pval = -np.log10(pval_sorted[i])
                                    p_values_list = [pval_sorted[i]]
                                    annotations = annotations_sorted[i]
                                    if len(deleted_nodes_sorted[i]) > 0:
                                        sum_score = 0
                                        p_values_list = []
                                        weighted_p_values_list = []
                                        for dg in deleted_nodes_sorted[i]:
                                            if current_deletion is None:
                                                current_deletion = {}
                                            current_deletion[dg] = {
                                                "start_deletion": deleted_nodes_sorted[i][dg]["start_deletion"],
                                                "end_deletion": deleted_nodes_sorted[i][dg]["end_deletion"]}
                                    else:
                                        current_deletion = None
                                        sum_score = size_sorted[i] * (-np.log10(pval_sorted[i]))
                                        weighted_p_values_list = [size_sorted[i]]
                                        p_values_list = [pval_sorted[i]]
                                    region_start = nodes_position_sorted[i]
                                    region_stop = region_start + size_sorted[i]
                        #compute the last node
                        #comptute the last deletion
                        if current_deletion is not None:
                            gap = []
                            for dg in current_deletion:
                                if "start_deletion" in current_deletion[dg] and current_deletion[dg][
                                    "start_deletion"] >= 0:
                                    gap.append(abs(current_deletion[dg]["end_deletion"] - current_deletion[dg][
                                        "start_deletion"]))
                            shared_deleted_size += statistics.median(gap)
                            sum_score += shared_deleted_size * (-np.log10(pval_sorted[i]))
                            weighted_p_values_list.append(shared_deleted_size)
                            p_values_list.append(pval_sorted[i])
                        if shared_deleted_size > 0 and region_stop < region_start + shared_deleted_size:
                            region_stop = region_start + shared_deleted_size
                        if region_start == region_stop:
                            if shared_deleted_size > 0:
                                region_stop = region_start + shared_deleted_size
                            else:
                                region_stop += 100
                                region_start -= 100
                                # Minimal region to allow visualization
                        if region_stop - region_start < 200:
                            min_size_region = 200
                            gap = int((min_size_region - (region_stop - region_start)) / 2)
                            region_start = max(0, region_start - gap)
                            region_stop = region_stop + gap
                        nb_nodes_in_region += 1

                        # Kc[c]["shared"] += shared_size
                        # Kc[c]["deleted"] += shared_deleted_size
                        # Kc["total"]["shared"] += shared_size
                        # Kc["total"]["deleted"] += shared_deleted_size
                        acat = compute_acat(p_values_list, weighted_p_values_list)
                        dic_regions[c][g]["regions"].append({"start": region_start, "stop": region_stop,
                                                             "shared_size": shared_size,
                                                             "shared_deleted_size": shared_deleted_size,
                                                             "region_size": region_stop - region_start,
                                                             "annotations": annotations,
                                                             "score":sum_score/(max(region_stop - region_start, 1)),
                                                             "pval": acat})

            nb_regions = 0
            logger.debug(f"Genome used for annotations : {genome_ref}")
            for c in chromosome_list:
                nb_regions += len(dic_regions[c][genome_ref]["regions"])
            logger.debug(f"Total number of identified nodes : {nb_regions_total} - Total regions nb : {nb_regions}")
            #Transform the result to get a dictionnary with main key = genome
            dic_regions_2 = {}
            for c, genomes in dic_regions.items():
                for genome, valeur in genomes.items():
                    if genome not in dic_regions_2:
                        dic_regions_2[genome] = {}
                    dic_regions_2[genome][c] = valeur
            analyse = {}
            logger.debug(f"genomes : {list(dic_regions_2.keys())}")
            total = sum(len(dic_regions_2[g][c]['regions']) for g in dic_regions_2 for c in dic_regions_2[g])
            #Check if annotations exist in databse
            annotations = True
            annotations_nb = 0
            query = 'MATCH (a:Annotation) RETURN count(a) AS annotations_nb'
            with driver.session() as session:
                result = session.run(query)
                for record in result:
                    annotations_nb = record["annotations_nb"]
            if annotations_nb <= 0:
                annotations = False

            # Get annotations for the regions
            with tqdm(total=total, desc="Processing regions") as pbar:
                for g in dic_regions_2:
                    logger.debug(f"genome : {g}")
                    analyse[g] = []
                    total_regions = 0
                    if annotations and ((genome_ref is not None and g == genome_ref) or (
                            genome_ref is None and g == genomes_list[0])):
                        for c in dic_regions_2[g]:
                            total_regions += len(dic_regions_2[g][c]['regions'])
                    logger.debug(f"total regions : {total_regions}")
                    if max_gwas_region is not None and total_regions >= max_gwas_region:
                        logger.debug(
                            f"Too much regions found {total_regions}, keep only the {max_gwas_region} first regions.")
                        total_regions = max_gwas_region
                    cpt_regions = 0
                    last_percent = -1

                    for c in dic_regions_2[g]:
                        # logger.debug(f"genome : {g} - chromosome : {c} - regions number : {len(dic_regions_2[g][c]['regions'])}")
                        for r in dic_regions_2[g][c]['regions']:
                            if max_gwas_region is not None and cpt_regions >= max_gwas_region:
                                break

                            # k = r["shared_size"]+r["shared_deleted_size"]
                            # Nc = chromosome_stats.get(c+"_max_position_mean")
                            # K = Kc[c]["shared"]+Kc[c]["deleted"]
                            # n = r["region_size"]
                            # pval_total = hypergeom.sf(k-1,
                            #                          N,
                            #                          Kc["total"]["shared"]+Kc["total"]["deleted"],
                            #                          r["region_size"]
                            #                          )
                            #test hypergeom
                            #pval_chrom = hypergeom(k - 1, N, K, n)
                            #r["pval"] = pval_chrom

                            r["chromosome"] = c
                            r["genome"] = g
                            # if annotations and ((genome_ref is not None and g == genome_ref) or (
                            #         genome_ref is None and g == genomes_list[0])):
                            if annotations and ((genome_ref is not None and g == genome_ref) or (
                                    genome_ref is None and g == genomes_list[0])):
                                cpt_regions += 1
                                percent = int((cpt_regions / total_regions) * 25)
                                progress = 75 + percent
                                if percent > last_percent+1:
                                    if get_gwas_status(job_id) == "CANCEL":
                                        logger.debug(f"Cancel job id : {job_id}")
                                        #delete_job(job_id)
                                        return {},{},{}
                                    update_gwas_progression(job_id, progress)
                                    last_percent = percent

                                annot_before_tmp = get_annotation_before_or_after_position(genome_ref=g, chromosome=c,
                                                                                           position=r["start"],
                                                                                           before=True, driver=driver)
                                annot_tmp = {}
                                if annot_before_tmp is not None and "gene_name" in annot_before_tmp:
                                    annot_tmp["gene_name"] = annot_before_tmp["gene_name"]
                                    annot_tmp["distance"] = int(r["start"] - annot_before_tmp["neighbor_node_position"]+ 1)
                                r["annotation_before"] = annot_tmp

                                # logger.debug("Search annotation after")
                                annot_after_tmp = get_annotation_before_or_after_position(genome_ref=g, chromosome=c,
                                                                                          position=r["stop"],
                                                                                          before=False, driver=driver)
                                annot_tmp = {}

                                if annot_after_tmp is not None and "gene_name" in annot_after_tmp:
                                    annot_tmp["gene_name"] = annot_after_tmp["gene_name"]
                                    annot_tmp["distance"] = int(annot_after_tmp["neighbor_node_position"] - r["stop"] + 1)
                                r["annotation_after"] = annot_tmp
                            analyse[g].append(r)
                            pbar.update(1)
                    update_gwas_progression(job_id, 100)
            # If reference genome is not in the selected genomes list then get the data relative to its genome
            if genome_ref not in analyse:
                analyse[genome_ref] = []
                total_regions = len(analyse[genomes_list[0]])
                if max_gwas_region is not None and total_regions >= max_gwas_region:
                    total_regions = max_gwas_region
                cpt_regions = 0
                for a in tqdm(analyse[genomes_list[0]], desc="Getting reference genome annotations"):
                    if max_gwas_region is not None and cpt_regions >= max_gwas_region:
                        break
                    cpt_regions += 1
                    percent = int((cpt_regions / total_regions) * 25)
                    progress = 75 + percent
                    if percent > last_percent + 1:
                        if get_gwas_status(job_id) == "CANCEL":
                            logger.debug(f"Cancel job id : {job_id}")
                            # delete_job(job_id)
                            return {}, {}, {}
                        update_gwas_progression(job_id, progress)
                        last_percent = percent
                    r = {}
                    r["chromosome"] = a["chromosome"]
                    r["shared_size"] = a["shared_size"]
                    r["shared_deleted_size"] = a["shared_deleted_size"]
                    r["genome"] = genome_ref
                    n_start, core_genome_start = get_genome_position(genome_ref, genomes_list[0], a["chromosome"],
                                                                     a["start"], before=True)
                    n_stop, core_genome_start = get_genome_position(genome_ref, genomes_list[0], a["chromosome"],
                                                                    a["stop"], before=False)
                    # n_start, core_genome_start = get_anchor(genomes_list[0], a["chromosome"], a["start"], before = True)
                    # n_stop, core_genome_start = get_anchor(genomes_list[0], a["chromosome"], a["stop"], before = False)
                    position_field = genome_ref + "_position"
                    r["start"] = 0
                    if position_field in n_start:
                        r["start"] = n_start[position_field]
                    if position_field in n_stop:
                        r["stop"] = n_stop[position_field]
                    else:
                        r["stop"] = r["start"]
                    r["region_size"] = r["stop"] - r["start"]
                    r["annotations"] = get_annotations_in_position_range(genome_ref=genome_ref,
                                                                         chromosome=a["chromosome"],
                                                                         start_position=r["start"],
                                                                         end_position=r["stop"])
                    annot_before_tmp = get_annotation_before_or_after_position(genome_ref=genome_ref,
                                                                               chromosome=a["chromosome"],
                                                                               position=r["start"], before=True)
                    annot_tmp = {}
                    if annot_before_tmp is not None and "gene_name" in annot_before_tmp:
                        annot_tmp["gene_name"] = annot_before_tmp["gene_name"]
                        annot_tmp["distance"] = int(r["start"] - annot_before_tmp["neighbor_node_position"] + 1)
                    r["annotation_before"] = annot_tmp

                    annot_after_tmp = get_annotation_before_or_after_position(genome_ref=genome_ref,
                                                                              chromosome=a["chromosome"],
                                                                              position=r["stop"], before=False)
                    if annot_after_tmp is not None and "gene_name" in annot_after_tmp:
                        annot_tmp["gene_name"] = annot_after_tmp["gene_name"]
                        annot_tmp["distance"] = int(annot_after_tmp["neighbor_node_position"] - r["stop"] + 1)
                    r["annotation_after"] = annot_tmp
                    analyse[genome_ref].append(r)

            for g in analyse:
                analyse[g] = sorted(analyse[g], key=lambda d: d['shared_size'], reverse=True)
            logger.debug("Total time : " + str(time.time() - temps_depart))
    return dic_regions_2, analyse, dic_distribution


def calculer_variabilite(chromosome_list=None, ref_genome=None, window_size=1000,
                         output_html="pangenome_variability.html"):
    with get_driver() as driver:
        if ref_genome == None:
            ref_noeud = "node_mean"
            ref_position = "position_mean"
        else:
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


# This function takes nodes data (from get_nodes_by_region or get_nodes_by_feature for exemple)
# it computes the Jaccard distance on these nodes
# and it returns the distance matrix and the distance matrix weighted by the nodes size
def compute_phylo_tree_from_nodes(nodes_data, output_dir="", weighted=False):
    genomes = get_genomes()
    nodes_nb = len(nodes_data)
    genome_redondant_strand_matrix = np.zeros((len(genomes), 2 * nodes_nb), dtype=np.int32)
    genome_strand_matrix = np.zeros((len(genomes), 2 * nodes_nb), dtype=np.int32)

    i = 0
    index_genomes = {}
    genomes_names = []
    for g in genomes:
        index_genomes[g] = i
        genomes_names.append(g)
        i += 1
    i = 0

    for n in nodes_data:
        for g in genomes:
            if "strandM" in nodes_data[n] and g in nodes_data[n]["strandM"]:
                genome_redondant_strand_matrix[index_genomes[g], i] += 1
                genome_strand_matrix[index_genomes[g], i] = nodes_data[n]["size"]
            if "strandP" in nodes_data[n] and g in nodes_data[n]["strandP"]:
                genome_redondant_strand_matrix[index_genomes[g], i + nodes_nb] += 1
                genome_strand_matrix[index_genomes[g], i + nodes_nb] = nodes_data[n]["size"]
        i += 1
    # computes Jaccard distance on matrix
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

            jaccard_matrix[i, j] = 1 - jaccard_index
            jaccard_matrix[j, i] = 1 - jaccard_index

            weighted_jaccard_matrix[i, j] = 1 - weighted_jaccard_index
            weighted_jaccard_matrix[j, i] = 1 - weighted_jaccard_index

    df_jaccard = pd.DataFrame(jaccard_matrix, index=genomes_names, columns=genomes_names)
    df_weighted_jaccard = pd.DataFrame(weighted_jaccard_matrix, index=genomes_names, columns=genomes_names)
    if output_dir != "" and os.path.isdir(output_dir):
        df_jaccard.to_csv(output_dir + '/distance_matrix.csv')
        df_weighted_jaccard.to_csv(output_dir + '/weighted_distance_matrix.csv')

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
    # logger.debug(newick_tree)
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


"""
This function returns a list of randomly sampled nodes using internal Neo4j IDs.

Strategy:
- Fetch min and max internal IDs for the given label.
- Randomly generate IDs in that range.
- Query nodes using WHERE id(n) IN [...]
- Retry until the desired number of valid nodes is obtained.
"""


def sample_random_nodes(sample_size, chromosome=None, label="Node", driver=None):
    if chromosome is None or chromosome == "":
        query = f"""
        MATCH (n:Node)
        with n, rand() AS r
        order by r
        limit {sample_size}
        return distinct(n) as nodes
        """
    else:
        query = f"""
        MATCH (n:Node)
        WHERE n.chromosome = '{chromosome}'
        with n, rand() AS r
        order by r
        limit {sample_size}
        return distinct(n) as nodes
        """
    nodes_list = []
    if not driver:
        driver = get_driver()
    with driver.session() as session:
        result = session.run(query)
        for record in result:
            nodes_list.append(dict(record["nodes"]))
    return nodes_list

def compute_phylo_params_hash(params):

    normalized = json.dumps(params, sort_keys=True)
    return hashlib.sha256(normalized.encode()).hexdigest()



def get_existing_global_tree(method="raxml", strand=True, chromosome=None,
                                         project_name="panabyss_phylo_tree", max_nodes=1000000, min_sample_size=10000,
                                         min_nodes_number=1000, force_reload=False):
    params = {"method": method, "strand":strand, "chromosome": chromosome}
    params_hash = compute_phylo_params_hash(params)

    phylo_job = find_existing_phylo_job(params_hash)

    if phylo_job and phylo_job["status"] == "SUCCESS":
        return phylo_job["global_tree"]
    else:
        return None


def run_phylo_job(job_id, params, params_hash, method="raxml", strand=True, chromosome=None,
                                         project_name="panabyss_phylo_tree", max_nodes=1000000, min_sample_size=10000,
                                         min_nodes_number=1000, force_reload=False):
    if job_id:
        try:
            logger.debug(f"Computing the global tree - job_id : {job_id}.")
            set_phylo_job_running(job_id, params, params_hash)
            newick_tree = compute_global_phylo_tree_from_nodes(method=method, strand=strand, chromosome=chromosome,
                                                           base_dir="./export/phylo/"+job_id, project_name=project_name, max_nodes=max_nodes,
                                                           min_sample_size=min_sample_size, min_nodes_number=min_nodes_number, job_id=job_id)

            if newick_tree :
                set_phylo_job_success(job_id, newick_tree)

        except Exception as e:
            set_phylo_job_error(job_id, str(e))

#This function if juste a wrapper of compute_global_phylo_tree_from_nodes to manage job in sqlite
#It is only required for job launched from IHM
def compute_global_phylo_tree_from_nodes_wrapper(method="raxml", strand=True, chromosome=None,
                                         project_name="panabyss_phylo_tree", max_nodes=1000000, min_sample_size=10000,
                                         min_nodes_number=1000, force_reload=False):
    params = {"method": method, "strand":strand, "chromosome": chromosome}
    params_hash = compute_phylo_params_hash(params)

    phylo_job = find_existing_phylo_job(params_hash)

    if phylo_job and phylo_job["status"] == "SUCCESS" and not force_reload:
        return {
            "status": "SUCCESS",
            "job_id": phylo_job["job_id"],
            "newick_tree": phylo_job["global_tree"]
        }

    if phylo_job and phylo_job["status"] == "RUNNING":
        return {
            "status": "RUNNING",
            "job_id": phylo_job["job_id"]
        }
    #If Job doesn't exist => create it and launch it
    if not phylo_job or "job_id" not in phylo_job:
        job_id = str(uuid.uuid4())
        insert_phylo_job(job_id, params, params_hash)
    else:
        job_id = phylo_job["job_id"]

    thread = threading.Thread(
        target=run_phylo_job,
        args=(job_id, params, params_hash, method, strand, chromosome, project_name, max_nodes, min_sample_size, min_nodes_number),
        daemon=True
    )
    thread.start()

    return {
        "status": "RUNNING",
        "job_id": job_id
    }



# This function compute a global tree from a random selection of nodes taking account of direct / reverse traversing (if strnd is True)
# Parameters :
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
def compute_global_phylo_tree_from_nodes(method="raxml", strand=True, chromosome=None,
                                         base_dir="./export/phylo", project_name="panabyss_phylo_tree",
                                         max_nodes=1000000, min_sample_size=10000,min_nodes_number=1000, job_id=None):
    total_nodes_number = 0
    distance_matrix_filename = "distance_matrix.phy"
    distance_matrix_phylip_filename = os.path.join(base_dir, distance_matrix_filename)
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    driver = get_scoped_driver()
    if driver is None:
        return None


    if chromosome is None:
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
        for record in result:
            total_nodes_number = record["total_nodes_number"]

    # The number of sampled nodes depends on the total number of nodes. It is set to be at least min_sample_size and at most max_nodes.
    # If the pangenome contains less than min_nodes_number then no tree can be computed.
    # The number of sampled nodes is determined by a polynomial interpolation in log-space of the total number of nodes, fitted through a set of control points.
    if total_nodes_number < min_nodes_number:
        return None
    if total_nodes_number < min_sample_size:
        sample_nodes_number = total_nodes_number
    else:
        # Set of control points :
        x_vals = np.array([1e4, 1e5, 1e7, 1e8, 2e9])
        y_vals = np.array([1e4, 3e4, 9e4, 1e5, 2e5])
        logx = np.log10(x_vals)
        # Compute the polynomial interpolation in log space:
        coeffs = np.polyfit(logx, y_vals, deg=4)

        def f(x):
            return np.polyval(coeffs, np.log10(x))

        sample_nodes_number = max(min_sample_size, min(int(f(total_nodes_number)), max_nodes))
    # sample_nodes_number = max(min_sample_size,min(int(node_selection_percentage*total_nodes_number/100), max_nodes))
    logger.info(f"Number of nodes to sample : {sample_nodes_number} - Total node {total_nodes_number}")
    nodes_list = sample_random_nodes(sample_nodes_number, chromosome, driver)

    logger.info(f"Number of sampled nodes : {len(nodes_list)}")
    sample_size = len(nodes_list)
    if sample_size >= min_sample_size:
        # Prepare PAV matrix
        genomes = get_genomes()
        pav_matrix = {}
        for g in genomes:
            if g not in pav_matrix:
                if strand:
                    pav_matrix[g] = np.zeros(2 * sample_size, dtype=int)
                else:
                    pav_matrix[g] = np.zeros(sample_size, dtype=int)

        for i in range(0, len(nodes_list)):
            for g in nodes_list[i]["genomes"]:
                if strand:
                    if "strandP" in nodes_list[i] and g in nodes_list[i]["strandP"]:
                        pav_matrix[g][i] = int(1)
                    else:
                        pav_matrix[g][i + sample_size] = int(1)
                else:
                    pav_matrix[g][i] = int(1)
        pav_to_phylip(pav_matrix, distance_matrix_phylip_filename)
        #Check if job has been canceled
        if job_id:
            status = get_phylo_status(job_id)
            if status == "CANCEL":
                logger.debug(f"Phylo job {job_id} canceled")
                return None
        if method == "raxml":
            logger.debug("RaxML method...")
            dir_raxml = os.path.join(base_dir, "raxml")
            tree_newick_filename = os.path.join(dir_raxml, project_name + ".raxml.bestTree")

            # Copier le fichier
            if not os.path.exists(dir_raxml):
                os.makedirs(dir_raxml)
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
                    return f.read()
            except FileNotFoundError:
                return None
        else:
            # method with matrix distance and neighbor joining
            logger.debug("Neighbor joining method...")
            names = list(pav_matrix.keys())
            matrix = np.array([pav_matrix[name] for name in names])
            n = len(names)
            lower_tri = []
            for i in range(n):
                lower_tri.append([])
                for j in range(i + 1):
                    # Hamming distance
                    lower_tri[i].append(np.sum(matrix[i] != matrix[j]) / matrix.shape[1])
            dm = _DistanceMatrix(names, lower_tri)
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(dm)
            newick_tree = tree.format('newick')
            return newick_tree

    else:
        return None


# Computes the distance matrix on the whole GFA (could take a long time for big GFA)
def compute_distance_matrix(distance_matrix_filename="distances.csv", chromosome=None, ponderation=True, strand=False):
    driver = get_driver()
    if driver is None:
        return None
    temps_depart = time.time()
    with driver.session() as session:
        # Get genomes list
        query = """
        MATCH (s:Stats)
        RETURN s.genomes AS genomes
        LIMIT 1
        """
        result = session.run(query)
        for record in result:
            genomes = record["genomes"]

        distance_matrix = pd.DataFrame(data=0, index=genomes, columns=genomes)
        dic_size_genome = {}
        if ponderation:
            # Get the genome size
            query = """
                MATCH (n:Node)
                UNWIND n.genomes AS g
                RETURN g AS genome, sum(n.size) AS total_size
                """
            result_size_genomes = list(session.run(query))
            if strand == False:
                # computes size for each genome

                # Get intersection size
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
            # Get the node numbers
            query = """
                MATCH (n:Node)
                UNWIND n.genomes AS g
                RETURN g AS genome, count(*) AS total_size
            """
            result_size_genomes = list(session.run(query))
            if strand == False:

                # Get intersection size
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
            distance_matrix.loc[g1, g2] = 1 - inter / (dic_size_genome[g1] + dic_size_genome[g2] - inter)
            distance_matrix.loc[g2, g1] = distance_matrix.loc[g1, g2]

        distance_matrix.to_csv(distance_matrix_filename)
        logger.debug("Total time : " + str(time.time() - temps_depart))
    return distance_matrix




