import re
from tqdm import tqdm
from math import *
from neo4j import GraphDatabase
from neo4j.exceptions import Neo4jError
import time
import hashlib
import logging
import scipy.sparse as sp
import numpy as np
import os
import csv
from config import *
from auth_utils import require_authorization
import logging


logger = logging.getLogger("panabyss_logger")

#Maximal length for S line
#Used for neo4j import procedure
#A specific graph with big nodes may require to increase this value
MAX_READ_BUFFER_SIZE_VALUE = 128000000

#These functions allow to create database and index from gfa and annotations files


#Version of BDD
#The version relate to the DB structure
from app import DB_VERSION
"""DB structure for this version

Nodes :

Stats :
    - genomes : list of all genomes of the pangenome
    - chromosomes : list of all chromosome of the pangenome

Sequence :
    - name : unique name of the node, must be the same than ref_node of nodes of type Node
    - sequence : String (DNA sequence of the node)

Annotation : attributes are directly related to gff / gtf format and not all detailed
    - name : hash of the annotation (unique)
    - chromosome
    - genome_ref : genome associated with this annotation
    - source
    - feature
    - filename : gtf or  gff filename used to create the node
    - id
    - gene_version
    - gene_name
    - gene_source
    - gene_biotype
    - gene_id
    - transcript_id
    - transcript_name
    - transcript_source
    - transcript_biotype
    - transcript_version
    - exon_id
    - exon_number
    - protein_id
    - protein_version
    - tag

Node:
    - name : String and unique
    - genomes : List of genomes that pass through this node
    - max : for the reference node (in the case of redundant nodes) the number of degenerate (redundant) nodes
    - strandP : indicates if the genome pass through the node in direct mode
    - strandM : indicates if the genome pass through the node in reverse mode
    - ref_node : main node, usefull for degenerate nodes. This name must be the same than the sequence node name
    - $genome_node : node number for the genome $genome
    - $genome_position : position of the start of the node for the genome $genome
    - size : size of the node
    - chromosome : chromosome associated to the node
    - position_min : the min of all $genome_position of the node
    - position_max : the max of all $genome_position of the node
    - position_mean : the mean of all $genome_position of the node
    - flow : percentage of all genomes of the pangenome passing throug this node

Relationships :
    - :gfa_link : link between 2 nodes (2 nodes can be linked be a unique link)
    - :A_POUR_ANNOTATION : link between a node of type Node and a node of type Annotation

"""

#batch_size_BDD size of batch transaction in DB
batch_size_BDD = 10000


logging.getLogger("neo4j").setLevel(logging.ERROR)

#This value allow to limit complex annotations search
ANNOTATION_SEARCH_LIMIT = 10000

@require_authorization
def create_nodes_batch(session, nodes_dic, node_name="Node", create = False):

    nb_transactions = max(1,ceil(len(nodes_dic)/batch_size_BDD))
    current_transaction = 0
    nodes_list = list(nodes_dic.items())
    with tqdm(total=nb_transactions) as bar :
        while len(nodes_dic)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD                
                batch = nodes_list[ind_depart:ind_depart+min(batch_size_BDD, len(nodes_dic)-current_transaction*batch_size_BDD)]
                if create == False :
                    query = (
                    "UNWIND $batch AS node "
                    "MERGE (n:"+str(node_name) +" {name: node.name}) "
                    "SET n += node.attributes "
                    )
                else :
                    query = (
                    "UNWIND $batch AS node "
                    "CREATE (n:"+str(node_name) +" {name: node.name}) "
                    "SET n += node.attributes "
                    )
                
                batch_data = [{"name": name, "attributes": attributes} for name, attributes in batch]
                tx.run(query, batch=batch_data)
                tx.commit()
    
                bar.update(1)
                current_transaction += 1
    return


# This function will get all chromosomes present in the pangenome graph
def get_chromosomes():
    if get_driver() is None:
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

#This function create stats about chromosomes
@require_authorization
def create_chromosome_stats():
    chromosomes = get_chromosomes()
    if not chromosomes:
        return

    chromosome_max_values = {}
    if get_driver() is None:
        return
    with get_driver() as driver:
        with driver.session() as session:
            for chrom in tqdm(chromosomes, desc="Processing chromosomes"):
                logger.debug(f"Getting stats for chromosome : {chrom}")
                max_result = session.run(
                    """
                    MATCH (n:Node {chromosome: $chrom})
                    RETURN max(n.position_mean) AS max_pos
                    """,
                    chrom=chrom
                )
                max_val = max_result.single()["max_pos"]
                key_name = f"{chrom}_max_position_mean"
                chromosome_max_values[key_name] = max_val

            cs_result = session.run(
                "MATCH (cs:chromosome_stats) RETURN cs LIMIT 1"
            )
            cs_record = cs_result.single()
            if cs_record:
                #Stats already exists => update
                # set_clause = ", ".join([f"cs.{k} = $props.{k}" for k in chromosome_max_values])
                # session.run(
                #     f"MATCH (cs:chromosome_stats) SET {set_clause}",
                #     props=chromosome_max_values
                # )
                session.run(
                    "MATCH (cs:chromosome_stats) SET cs += $props",
                    props=chromosome_max_values
                )
                print("Updated chromosome_stats node with:", chromosome_max_values)
            else:
                #Stats doesn't exist => create it
                session.run(
                    "CREATE (cs:chromosome_stats) SET cs += $props",
                    props=chromosome_max_values
                )
                print("Created chromosome_stats node with:", chromosome_max_values)
    logger.info("✅ Chromosome stats node created")
    logger.debug(f"✅ Chromosome stats value : {chromosome_max_values}")
    return


@require_authorization
def create_stats(set_genomes, set_chromosomes):
    with get_driver() as driver:
        with driver.session() as session:
            with session.begin_transaction() as tx:
                query ="""
                MERGE (s:Stats)
                WITH s, coalesce(s.genomes, []) + $genomes AS all_genomes, $chromosomes as liste_chromosome
                UNWIND all_genomes AS g
                WITH s, collect(DISTINCT g) AS new_genomes, liste_chromosome
                SET s.genomes = new_genomes, s.version=$version
                
                // Mise à jour de s.chromosomes
                WITH s, coalesce(s.chromosomes, []) + liste_chromosome AS all_chromosomes
                UNWIND all_chromosomes AS c
                WITH s, collect(DISTINCT c) AS new_chromosomes
                SET s.chromosomes = new_chromosomes
                """
                tx.run(query, genomes=list(set_genomes), chromosomes = list(set_chromosomes), version=DB_VERSION)
    create_chromosome_stats()
    return

@require_authorization
def create_stats_from_nodes():
    with get_driver() as driver:
        with driver.session() as session:
            # Step 1: Get genomes
            query_genomes = """
                MATCH (n:Node)
                WHERE n.flow = 1
                RETURN n.genomes AS all_genomes
                LIMIT 1
            """
            result = session.run(query_genomes)
            all_genomes = result.single()["all_genomes"] if result.peek() else []

            # Step 2: Get chromosomes
            query_chromosomes = """
                MATCH (n:Node)
                USING INDEX n:Node(chromosome)
                WHERE n.chromosome IS NOT NULL
                RETURN DISTINCT n.chromosome AS all_chromosomes
            """
            result = session.run(query_chromosomes)
            all_chromosomes = [record["all_chromosomes"] for record in result if record["all_chromosomes"] is not None]

            # Step 3: Delete existing Stats node (if exists)
            session.run("MATCH (s:Stats) DETACH DELETE s")

            # Step 4: Create new Stats node
            session.run("""
                CREATE (s:Stats {
                    genomes: $genomes,
                    chromosomes: $chromosomes,
                    version: $version
                })
            """, genomes=all_genomes, chromosomes=all_chromosomes, version=DB_VERSION)

            logger.info(f"✅ Stats node created with genomes : {all_genomes} - chromosomes : {all_chromosomes}")
    create_chromosome_stats()
    return


#Check if there are index in state "POPULATING"
def indexes_populating():
    query = "SHOW indexes YIELD name, state RETURN name, state"
    
    with get_driver() as driver:
        with driver.session() as session:
            result = session.run(query)
            for record in result:
                if record["state"] == "POPULATING":
                    return True
            return False

"""
Query Neo4j for index statuses.
If only_in_progress=True, return only indexes that are not ONLINE.
Otherwise, return all indexes.
"""
def get_index_statuses(only_in_progress=True):

    query = "SHOW INDEXES;"
    in_progress_states = {"POPULATING", "FAILED"}  # depending on Neo4j version
    indexes = {}
    with get_driver() as driver:
        with driver.session() as session:
            result = session.run(query)

            for record in result:
                name = record.get("name")
                state = record.get("state")
                pct = record.get("populationPercent")

                if only_in_progress:
                    if state != "ONLINE":
                        indexes[name] = {
                            "state": state,
                            "population_percent": pct
                        }
                else:
                    indexes[name] = {
                        "state": state,
                        "population_percent": pct
                    }

    return indexes

"""
Wait until all indexes that are currently in progress finish building.

Behavior:
- First fetches the list of indexes that are NOT ONLINE.
- If no indexes are in progress, immediately returns success.
- Monitors them until they become ONLINE.
- Detects FAILED indexes and returns an error.
- Detects lack of progress for 'max_no_progress' iterations and stops.
"""
def wait_for_indexes(poll_interval = 20, max_no_progress = 6):
    logger.info("Checking for index creation")
    # --- Initial scan: get only indexes currently building ---
    in_progress = get_index_statuses(only_in_progress=True)
    if not in_progress:
        return (0, "OK: No indexes are currently being created.")


    index_names = list(in_progress.keys())

    # Track last progress to detect stalling
    last_progress = {idx: -1 for idx in index_names}
    no_progress_counter = 0

    # --- tqdm setup (global progress bar) ---
    pbar = tqdm(
        total=100.0,
        desc="Index creation progress",
        position=0,
        leave=True
    )
    pbar.update(0)

    last_global_pct = 0.0

    while True:
        # Fetch full status to detect transitions to ONLINE or FAILED
        statuses = get_index_statuses(only_in_progress=False)
        all_online = True
        progress_made = False
        pct_values = []

        for idx in in_progress.keys():
            state = statuses[idx]["state"]
            pct = statuses[idx]["population_percent"]

            # If the index failed
            if state == "FAILED":
                pbar.close()
                return (1, f"Error: Index '{idx}' is in FAILED state.")

            # If not yet online → continue waiting
            if state != "ONLINE":
                all_online = False

            # Detect progress using populationPercent when available
            if pct is not None:
                pct_values.append(pct)
                if pct > last_progress[idx]:
                    last_progress[idx] = pct
                    progress_made = True
            else:
                pct_values.append(0)

            # Compute global progress as the average %
            global_pct = sum(pct_values) / len(pct_values)

            # Update tqdm bar
            delta = global_pct - last_global_pct
            if delta > 0:
                pbar.update(delta)
                last_global_pct = global_pct


        # If all indexes have finished successfully
        if all_online:
            pbar.update(100 - last_global_pct)
            pbar.close()
            return (0, "OK: All indexes are now ONLINE.")

        # Progress monitoring
        if not progress_made:
            no_progress_counter += 1
        else:
            no_progress_counter = 0

        # Detect a stalled indexing process
        if no_progress_counter >= max_no_progress:
            pbar.close()
            return (2, "Error: No progress detected in index creation (stalled).")

        # Wait before polling again
        time.sleep(poll_interval)



#This function return the state of creating index
#if index has been created this return 100
#else it returns the percentage of index creation
def check_state_index(index_name: str):
    index_name_formate = index_name.replace("-", "_").replace(".","_")
    with get_driver() as driver:
        with driver.session() as session:
            query = """
            SHOW INDEX YIELD name, state, populationPercent
            WHERE name = $index_name_formate
            RETURN populationPercent
            """

            try:
                result = session.run(query, index_name_formate=index_name_formate)
                record = result.single()
                if record:
                    return record["populationPercent"]
                else:
                    return None
            except Neo4jError as e:
                logger.error(f"❌ Error while checking index state: {e}")
                return None


#Function to create index in database
#If base = True => create the base indexes = index on Node name and chromosome
#If extend = True => create other indexes = index on Node flow, size, ref_node + indexes on Annotation name, chromosome, start, end, gene_id, gene_name + index on Sequence name
#If genomes_index = True => create indexes on Node chromosome / $genome_position
@require_authorization
def create_indexes(base=True, extend=False, genomes_index=False):
    indexes_queries = []
    with get_driver() as driver:
        with driver.session() as session:
            if base :
                indexes_queries= [
                    "CREATE INDEX NodeIndexName IF NOT EXISTS FOR (n:Node) ON (n.name)",
                    "CREATE INDEX NodeIndexChromosome IF NOT EXISTS FOR (n:Node) ON (n.chromosome)"
                    ]
            if extend :
                indexes_queries += [
                    "CREATE INDEX NodeIndexFlow IF NOT EXISTS FOR (n:Node) ON (n.flow)",
                    "CREATE INDEX NodeIndexSize IF NOT EXISTS FOR (n:Node) ON (n.size)",
                    "CREATE INDEX NodeIndexGwas IF NOT EXISTS FOR (n:Node) ON (n.chromosome, n.flow, n.size)",
                    "CREATE INDEX NodeIndexRefNode IF NOT EXISTS FOR (n:Node) ON (n.ref_node)",
                    "CREATE INDEX AnnotationName IF NOT EXISTS FOR (a:Annotation) ON (a.name)",
                    "CREATE INDEX AnnotationIndexChromosome IF NOT EXISTS FOR (a:Annotation) ON (a.chromosome)",
                    "CREATE INDEX AnnotationIndexStart IF NOT EXISTS FOR (a:Annotation) ON (a.start)",
                    "CREATE INDEX AnnotationIndexEnd IF NOT EXISTS FOR (a:Annotation) ON (a.end)",
                    "CREATE INDEX AnnotationIndexGeneId IF NOT EXISTS FOR (a:Annotation) ON (a.gene_id)",
                    "CREATE INDEX AnnotationIndexGeneName IF NOT EXISTS FOR (a:Annotation) ON (a.gene_name)",
                    "CREATE INDEX AnnotationIndexGenomeRef IF NOT EXISTS FOR (a:Annotation) ON (a.genome_ref)",
                    "CREATE INDEX AnnotationIndexAnnotationSearch IF NOT EXISTS FOR (a:Annotation) ON (a.genome_ref, a.chromosome, a.start, a.end)",
                    "CREATE INDEX SequenceIndexName IF NOT EXISTS FOR (s:Sequence) ON (s.name)"
                    ]
            with session.begin_transaction() as tx:
                for query in indexes_queries :
                    tx.run(query)
            
            indexes_queries = []
            
            if genomes_index :
                current_genome = 0
                #Uncomment the following lines if index creation takes too much ressources 
                #if extend :
                #    while indexes_populating():
                        #wait 1 minute between 2 poll of indexes states
                #        time.sleep(60)
                query_genomes = """
                MATCH (s:Stats) 
                RETURN s.genomes as all_genomes
                """
                all_genomes = []
                result = session.run(query_genomes)
                for record in result:
                    all_genomes = record["all_genomes"]
                nb_genomes = len(all_genomes)
                logger.info(all_genomes)
                for g in all_genomes:
                    logger.info("creating indexes for genome " + g + " ("+str(current_genome+1) + "/"+str(nb_genomes) +")")
                    current_genome += 1
                    indexes_queries = []
                    indexes_queries.append("CREATE INDEX NodeIndex"+str(g).replace("-", "_").replace(".","_")+"_position IF NOT EXISTS FOR (n:Node) ON (n.chromosome, n.`"+str(g)+"_position`)")
                    with session.begin_transaction() as tx:
                        for query in indexes_queries :
                            tx.run(query)
            


@require_authorization
def creer_index_chromosome(chromosomes_list = []):
    query_genomes = """
    MATCH (s:Stats) 
    RETURN s.genomes as all_genomes
    """
    all_chromosomes = chromosomes_list
    all_genomes = []
    indexes_queries = []
    with get_driver() as driver:
        with driver.session() as session:
            result = session.run(query_genomes)
            for record in result:
                all_genomes = record["all_genomes"]

        for c in all_chromosomes:
            indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"IndexName IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n.name)")
            indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"IndexFlow  IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n.flow)")
            indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"IndexSize IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n.size)")
            for g in all_genomes :
                indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"Index"+str(g)+"_position IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n."+str(g)+"_position)")
                indexes_queries.append("CREATE INDEX Node_chr_" + str(c) +"Index"+str(g)+"_node IF NOT EXISTS FOR (n:Node_chr_" + str(c) +") ON (n."+str(g)+"_node)")
    
        with get_driver() as driver:
            with driver.session() as session:
                with session.begin_transaction() as tx:
                    for query in indexes_queries :
                        tx.run(query)
                        


@require_authorization                        
def creer_index_chromosome_genomes():
    query_genomes = """
    MATCH (s:Stats) 
    RETURN s.genomes as all_genomes
    """
    indexes_queries = []
    all_genomes = []
    with get_driver() as driver:
        with driver.session() as session:
            result = session.run(query_genomes)
            for record in result:
                all_genomes = record["all_genomes"]

            nb_genomes = len(all_genomes)
            current_genome = 0
            for g in all_genomes :
                logger.info("creating indexes for genome " + g + " ("+str(current_genome) + "/"+str(nb_genomes) +")")
                current_genome += 1
                indexes_queries = []
                indexes_queries.append("CREATE INDEX NodeIndex"+str(g).replace("-", "_").replace(".","_")+"_position IF NOT EXISTS FOR (n:Node) ON (n.chromosome, n.`"+str(g)+"_position`)")
                indexes_queries.append("CREATE INDEX NodeIndex"+str(g).replace("-", "_").replace(".","_")+"_node IF NOT EXISTS FOR (n:Node) ON (n.chromosome, n.`"+str(g)+"_node`)")
                with session.begin_transaction() as tx:
                    for query in indexes_queries :
                        tx.run(query)   
                if current_genome % 5 == 0:
                    time.sleep(60*10)
    logger.info("Indexes created")


@require_authorization
def create_labels_chromosomes(liste_chromosomes=[]):
    all_chromosomes = []
    labels_queries = []
    with get_driver() as driver :
        if len(liste_chromosomes) == 0:
            query = """
            MATCH (s:Stats) 
            RETURN s.chromosomes as all_chromosomes
            """
            with driver.session() as session:
                result = session.run(query)
                for record in result:
                    all_chromosomes = record["all_chromosomes"]
    
        else:
            all_chromosomes = liste_chromosomes
        for c in all_chromosomes :
            label = "Node_chr_"+str(c)
            labels_queries.append(
                f"""
                CALL apoc.periodic.iterate(
              "MATCH (n:Node) WHERE n.chromosome = '{c}' RETURN n",
              "CALL apoc.create.addLabels(n, ['{label}']) YIELD node RETURN node",
              {{batchSize:1000, parallel:true}}
            )
            """
            )
        with driver.session() as session:
            with session.begin_transaction() as tx:
                for query in labels_queries:
                    #logger.debug(query)
                    tx.run(query)


@require_authorization
def creer_relations_batch(session, liste_relations):
    
    nb_transactions = ceil(len(liste_relations)/batch_size_BDD)
    current_transaction = 0
    with tqdm(total=nb_transactions) as bar :
        while len(liste_relations)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD

                batch = liste_relations[ind_depart:ind_depart+min(batch_size_BDD, len(liste_relations)-current_transaction*batch_size_BDD)]
                query = f"""
                        UNWIND $batch AS pair
                        MATCH (a:Node {{name: pair.depart}})
                        MATCH (b:Node {{name: pair.arrivee}})
                        MERGE (a)-[r:{"gfa_link"}]->(b)
                        """
                tx.run(query, batch=batch)        
                bar.update(1)
                current_transaction += 1
    return



'''
This function will create the nodes and their sequences (only name and sequence) in the neo4j database.
Then it will create the indexes: kmer and relationships with the nodes.
'''
@require_authorization
def creer_sequences_et_indexes(session, dic_kmer_relation, kmer_size, nodes_dic):
    
    nb_transactions = ceil(len(nodes_dic)/batch_size_BDD)
    current_transaction = 0
    nodes_list = list(nodes_dic.keys())
    logger.info("Start to create sequence into DB")
    with tqdm(total=len(nodes_dic)) as bar :
        while len(nodes_dic)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD
                nodes_list_transaction = nodes_list[ind_depart:ind_depart+min(batch_size_BDD, len(nodes_dic)-current_transaction*10000)]
    
                logger.debug("Transaction " + str(current_transaction))
                for t in range(min(batch_size_BDD, len(nodes_dic)-current_transaction*batch_size_BDD)):
                    nc = nodes_list_transaction[t]
                    q = """
                        CREATE (n:Sequence {name:$nom, sequence:$sequence})
                        """ 
                    result = tx.run(
                        q,
                        nom=nc,
                        sequence=nodes_dic[nc]
                    )
                bar.update(min(batch_size_BDD, len(nodes_dic)-current_transaction*batch_size_BDD))
                current_transaction += 1
    logger.info("End of sequence creation into DB")
    liste_kmer = list(dic_kmer_relation.keys())
    logger.info("Start of indexes creation")
    current_transaction = 0
    with tqdm(total=len(dic_kmer_relation)) as bar :
        while len(dic_kmer_relation)-current_transaction*batch_size_BDD > 0:
            # Démarrer une transaction
            with session.begin_transaction() as tx:
                ind_depart = current_transaction*batch_size_BDD
                liste_kmer_transaction = liste_kmer[ind_depart:ind_depart+min(batch_size_BDD, len(dic_kmer_relation)-current_transaction*10000)]
    
                logger.debug("Transaction " + str(current_transaction))
                for t in range(min(batch_size_BDD, len(dic_kmer_relation)-current_transaction*batch_size_BDD)):
                    kmer = liste_kmer_transaction[t]
                    q = """
                        CREATE (k:kmer {kmer:$kmer, size:$size})
                        WITH k 
                        MATCH (n:Sequence)
                        WHERE n.name IN $target_nodes
                        CREATE (k)-[:index]->(n)
                        """
                    result = tx.run(
                        q,
                        target_nodes=dic_kmer_relation[kmer],
                        kmer=kmer,
                        size=kmer_size
                    )
                bar.update(min(batch_size_BDD, len(dic_kmer_relation)-current_transaction*batch_size_BDD))
                current_transaction += 1
    logger.info("Indexes created into DB")
    return



'''
This function allows you to create nodes with the node sequence and its name.
It also creates indexes (kmers and associated nodes).
'''
@require_authorization
def load_sequences_and_indexes(gfa_file_name, kmer_size=31):
    nodes_dic = {}
    file = open(gfa_file_name, "r", encoding='utf-8')
    kmer_set = set()
    dic_kmer_relation = {}
    with file:
        total_nodes = sum(1 for line in file if line.startswith(('S')))
        file.seek(0,0)
        with tqdm(total=total_nodes) as bar :
            ligne = file.readline()
            while ligne:
                if ligne.startswith(('S')):
                    ligne_dec = ligne.split()
                    if len(ligne_dec) > 0:
                        bar.update(1)
                        seq = ligne_dec[2]
                        node = ligne_dec[1]
                        nodes_dic[node] = seq
                        if len(seq) > kmer_size:
                            liste_kmer = [seq[i:i+kmer_size] for i in range(len(seq)-kmer_size+1)]
                            for kmer in liste_kmer:
                                if kmer in kmer_set :
                                    dic_kmer_relation[kmer].append(node)
                                else:
                                    kmer_set.add(kmer)
                                    dic_kmer_relation[kmer] = [node]
                ligne = file.readline()
    file.close()
    return dic_kmer_relation, nodes_dic

'''
This function allows you to create nodes with the node sequence and name.
'''
@require_authorization
#TODO : découper en lot le chargement des noeuds
def load_sequences(gfa_file_name, chromosome_file = None, create=False, batch_size=20000000):
    nodes_dic = {}
    start_time = time.time()
    file = open(gfa_file_name, "r", encoding='utf-8')
    with get_driver() as driver:
        with file:
            total_nodes = sum(1 for line in file if line.startswith(('S')))
            file.seek(0,0)
            with tqdm(total=total_nodes) as bar :
                ligne = file.readline()
                while ligne:
                    if ligne.startswith(('S')):
                        ligne_dec = ligne.split()
                        if len(ligne_dec) > 0:
                            bar.update(1)
                            seq = ligne_dec[2]
                            if (chromosome_file != None and chromosome_file != "") :
                                node = chromosome_file + "_"+ ligne_dec[1]
                            else : 
                                node = ligne_dec[1]
                            nodes_dic[node] = {"sequence":seq}
                    if len(nodes_dic) >= batch_size:
                        with driver.session() as session:
                            create_nodes_batch(session, nodes_dic, node_name="Sequence", create = create)
                        nodes_dic = {}
                    ligne = file.readline()
        logger.info("Sequences nodes computed in " + str(time.time()-start_time))
        logger.info("Creating sequences in DB : " + str(len(nodes_dic)) + " sequences to create")
        if len(nodes_dic) > 0 :
            with driver.session() as session:
                create_nodes_batch(session, nodes_dic, node_name="Sequence", create = create)
    logger.info("Sequences created. Total time : " + str(time.time()-start_time) + " s")
    file.close()




    
letters_to_num = {
    'UN': '1',
    'DEUX': '2',
    'TROIS': '3',
    'QUATRE': '4',
    'CINQ': '5',
    'SIX': '6',
    'SEPT': '7',
    'HUIT': '8',
    'NEUF': '9',
    'DIX': '10',
    'ONZE': '11',
    'DOUZE': '12',
    'TREIZE': '13',
    'QUATORZE': '14',
    'QUINZE': '15',
    'SEIZE': '16',
    'DIXSEPT': '17',
    'DIXHUIT': '18',
    'DIXNEUF': '19',
    'VINGT': '20'
}


"""
Function to find the genome and chromosome of a P or W line
The P line is supposed to match one of the following pattern :
    - genome#haplotype#chromosome
    - genome.haplotype.chromosome
The W line (to be preferred) : genome = 2nd element of the line +"_" + 3rd element of the line, chromosome = 4th element of the line
If the file concerned only one chromosome with multiple haplotype, chromosome_file parameter should contain the chromosome reference (string)
Parameters :
    - WP_line : W or P line
    - haplotype : is haplotype is true the genome will be named : genome_haplotype, else only genome
    - chromosome_file : if set, the chromosome value is fixed to this value
"""
def get_chromosome_genome(WP_line, haplotype=True, chromosome_file = None):
    ligne_dec = WP_line.split()
    if ligne_dec[0] == 'P' or ligne_dec[0] == 'W':
        chromosome = "0"
        if ligne_dec[0] == 'P':
            if (len (ligne_dec[1].split("#")) > 1):
                name_dec = ligne_dec[1].split("#")
            else:
                name_dec = ligne_dec[1].split(".")
            if haplotype  and len(name_dec) > 1 :
                genome = name_dec[0]+"_"+name_dec[1]
            else:
                genome = name_dec[0]
            if len(name_dec) > 0 :
                chromosome = re.sub("^0*", "", name_dec[-1].upper().replace("CHR", ""))
        else:
            match = re.search(r'(chromosome|chr)\s*([a-zA-Z0-9]+)', ligne_dec[3], re.IGNORECASE)
            if match:
                chromosome = match.group(2)
                chromosome = chromosome.lstrip('0').upper()
                chromosome = chromosome or '0'  
                chromosome = letters_to_num.get(chromosome, chromosome)
            else:
                chromosome = '0'
            #chromosome = re.sub("^0*", "", str(ligne_dec[3]).upper().replace("CHR", ""))
            if haplotype :
                genome = ligne_dec[1]+"_"+ligne_dec[2]
            else : 
                genome = ligne_dec[1]
    if chromosome_file != None and chromosome_file != "":
        chromosome = chromosome_file
    genome = genome.replace("-", "_")
    return chromosome,genome

'''
This function creates nodes in Neo4j from a GFA file.
It processes the data in batches (maximum defined by batch_size) to avoid memory issues.
For batch processing: we iterate through the nodes and select batch_size nodes.
Then we go through the paths to check if we find the corresponding nodes; if yes, we update the node.
At the end of the paths, we create the batch of nodes in the database.
The nodes contain the following properties:
   - Node size
   - Associated chromosome
   - Node index for each sequence passing through the node
   - Position: the cumulative size of the previous nodes for each sequence on the chromosome
   - Master node: in case of sub-node creation to prevent a node from being used multiple times by the same sequence
                  the sub-nodes are named _2, _3, etc. For example, if the master node is s1, the sub-nodes 
                  will be s1_2, s1_3, etc.
   - List of genomes passing through this node
   - List of genomes passing directly through this node
   - List of genomes passing through this node in reverse
   - Flow: the percentage of genomes passing through this node
The tool works with P lines or W lines; however, since the P line format is not well-defined,
especially for chromosome and sequence naming, W lines are preferred.
For W lines: either the file relates to a single chromosome, or the chromosome is defined in the fourth column.
Chromosome information is required to properly use the tool.
Input:
    - gfa_file_name: GFA file to load
    - chromosome_file: if the file relates to a single chromosome, specify the chromosome name (string)
    - chromosome_prefix: if True, the node name is prefixed by the chromosome (this is the case if chromosome_file is provided); if False and no chromosome_file, then the node is not prefixed
    - batch_size: the splitting is done by nodes, the tool will process batches of size "batch_size"
    - start_node: in case there was an issue during loading, we can resume from a specific node
    - create: if it's the first run, this can be set to True; in this case, the tool will create nodes in the DB without checking for existence
              in other cases, it is better to set it to False
    - haplotype: indicates whether the sample name should be concatenated with the haplotype
'''
@require_authorization
def load_gfa_data_to_neo4j(gfa_file_name, chromosome_file = None, chromosome_prefix = False, batch_size = 2000000, start_chromosome = None, create = False, haplotype = True, create_only_relations = False):
    sep = ["[,;.*]", "(<|>)"]
    batch_nb = 0
    set_genome = set()
    set_chromosome = set()
    total_nodes = 0
    total_path = 0
    nodes_size_dic = {}
    temps_depart = time.time()
    set_all_genomes = set()
    set_all_chromosomes = set()
    liste_relations = []
    chromosomes_list = []
    set_relations = set()
    #Create base indexes (name and chromosome)
    create_indexes() 
    index_first_chromosme = 0
    nodes_set_next_chromosome = set()
    first_chromosome = None
    file = open(gfa_file_name, "r", encoding='utf-8')
    with file:
        
        #First file browsing to get length, nodes and haplotypes
        ligne = file.readline()
        while ligne :
            if ligne.startswith(('S')):
                ligne_dec = ligne.split()
                if (len(ligne_dec) > 0): 
                    nodes_size_dic[ligne_dec[1]]=int(len(ligne_dec[2]))
                    total_nodes += 1
            if ligne.startswith(('P',"W")):
                logger.debug(ligne[0:80])
                total_path += 1
                ligne_dec = ligne.split()
                chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                if chromosome not in set_all_chromosomes :
                    chromosomes_list.append(chromosome)
                if len(set_all_chromosomes) == 0:
                    first_chromosome = chromosome
                    if start_chromosome is not None and start_chromosome != "":
                        first_chromosome = start_chromosome
                set_all_genomes.add(genome) 
                set_all_chromosomes.add(chromosome)
                if chromosome == first_chromosome :
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([chaine[:-1] for chaine in nodes_list])
                    else:
                        ind = 6
                        walk = 1
                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0])
                
            ligne = file.readline() 
        if first_chromosome is not None and first_chromosome != "":
            for k in range(len(chromosomes_list)):
                if chromosomes_list[k] == first_chromosome:
                    index_first_chromosme = k
        logger.info("Genomes number : " + str(len(set_all_genomes)) + " - chromosomes list : " + str(set_all_chromosomes))
        #If create_only_relations is set to True then the nodes are not processed
        if not create_only_relations :
            logger.info("Start parsing, nodes number : " + str(total_nodes) + "\nstart chromosome : " + str(first_chromosome))
            for k in range(index_first_chromosme,len(chromosomes_list)) :
                c = chromosomes_list[k]
                nodes_set_chromosome = set(nodes_set_next_chromosome)
                nodes_set_next_chromosome = set()
                logger.debug("chromosome " + str(c) + " - number of nodes : " + str(len(nodes_set_chromosome)))
                batch_nb = ceil(len(nodes_set_chromosome)/batch_size)
                current_batch = 0
                while current_batch < batch_nb :
                    temps_0_lot = time.time()
                    nodes_batch_set = set(list(nodes_set_chromosome)[current_batch*batch_size:min(len(nodes_set_chromosome),(current_batch+1)*batch_size)])
                    current_batch += 1
                    logger.debug("chromosome " + c + " batch " + str(current_batch) + "/"+str(batch_nb) + " nodes number : " + str(len(nodes_batch_set)))
                    file.seek(0,0)
                    ligne = file.readline()
    
                    #Path browsing for the batch 
    
                    nodes_list = []
                    liste_strand = []
                    position_count = {}
                    nodes_count = {}
                    nodes_dic = {}
                    ref_nodes_dic = {}
                    set_genomes_lot = set()    
                    nodes_set = set()
                    with tqdm(total=total_path) as bar2 :
                        while ligne:
                            ligne_dec = ligne.split()
                            if len(ligne_dec) > 0:
                                if ligne[0] == 'P' or ligne[0] == 'W':
                                    chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                                    ligne = None
                                    if current_batch == batch_nb and k < len(chromosomes_list) - 1 and chromosome == chromosomes_list[k+1]:
                                        #last batch for the chromosome, retrieves the nodes to be processed for the next chromosome
                                        if ligne_dec[0] == 'P':
                                            ind = 2
                                            walk = 0
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            nodes_set_next_chromosome |= set([chaine[:-1] for chaine in nodes_list])
                                        else:
                                            ind = 6
                                            walk = 1
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            nodes_set_next_chromosome |= set([nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0])
                                    
                                    if chromosome == c:
                                        if ligne_dec[0] == 'P':
                                            ind = 2
                                            walk = 0
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            liste_strand = [chaine[-1] for chaine in nodes_list]
                                            nodes_list = [chaine[:-1] for chaine in nodes_list]
                                            # if chromosome_file == None :
                                            #     liste_strand = [chaine[-1] for chaine in nodes_list]
                                            #     nodes_list = [chaine[:-1] for chaine in nodes_list]
                                            # else :
                                            #     liste_strand = [chaine[-1] for chaine in nodes_list]
                                            #     nodes_list = [chromosome_file+"_"+chaine[:-1] for chaine in nodes_list]                              
                                        else:
                                            ind = 6
                                            walk = 1
                                            nodes_list = re.split(sep[walk], ligne_dec[ind])
                                            liste_strand = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 != 0]
                                            nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                            # if chromosome_file == None :
                                            #     nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                            # else:
                                            #     nodes_list = [chromosome_file+"_"+nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
          
                                        if chromosome is not None and chromosome not in set_chromosome:
                                            set_chromosome.add(chromosome)
                                        if genome != "_MINIGRAPH_":
                                            if genome not in set_genome:
                                                set_genome.add(genome)
                                            if genome not in set_genomes_lot :
                                                set_genomes_lot.add(genome)
                                                nodes_count[genome] = {}
                                                position_count[genome] = {}
                                            if chromosome not in nodes_count[genome] :
                                                nodes_count[genome][chromosome] = 0
                                                position_count[genome][chromosome] = {}
                                                position_count[genome][chromosome]["current_position"] = 0
                                                position_count[genome][chromosome]["previous_position"] = 0
                                                position_count[genome][chromosome]["current_contig"] = ""
                                            #For walk, start position is available => this position is used
                                            if ind == 6:
                                                if position_count[genome][chromosome]["current_contig"] != ligne_dec[3] :
                                                    #New contig => add the start of next contig
                                                    position_count[genome][chromosome]["current_position"] += int(ligne_dec[4])
                                                else :
                                                    #Same contig => add the potential gaps
                                                    if position_count[genome][chromosome]["previous_position"] - int(ligne_dec[4]) > 0 :
                                                        position_count[genome][chromosome]["current_position"] += position_count[genome][chromosome]["previous_position"] - int(ligne_dec[4])
                                                position_count[genome][chromosome]["current_contig"] = ligne_dec[3]
                                                position_count[genome][chromosome]["previous_position"] = int(ligne_dec[5])
        
                                            node = ""
                                            ref_node = ""
                                            strand = ""
                                            ligne_dec=None
                                            #Graph linearization
                                            # Browse nodes list : if node already exist for the same sequence
                                            # then create a new node (For example, if it is the sixth iteration for node S1, we will create S1_6)
                                            for i in range(len(nodes_list)):
                                                if i > 0:
                                                    previous_node = node
                                                node = nodes_list[i]
                                                ref_node = node
                                                size = nodes_size_dic[ref_node]
                                                if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                                    node = chromosome + "_" + node

                                                #Node is consider only if it is part of batch

                                                if ref_node in nodes_batch_set:
                                                    if chromosome_file is not None and chromosome_file != "":
                                                        ref_node = node
                                                    strand = ""
                                                    if (i < len(liste_strand) and liste_strand[i] in ["-", "<"]):
                                                        strand = "M"
                                                    else:
                                                        strand = "P"
                                                    
                                                    if node not in nodes_set:
                                                        if strand == "P" :
                                                            nodes_dic[node] = {"genomes":[genome], "max":1, "strandP":[genome], "strandM":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                        else :
                                                            nodes_dic[node] = {"genomes":[genome], "max":1, "strandM":[genome], "strandP":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                        nodes_set.add(node)
                                                    else:
                                                        if genome not in nodes_dic[node]["genomes"] and chromosome == nodes_dic[node]["chromosome"]:
                                                            nodes_dic[node]["genomes"].append(genome)
                                                            nodes_dic[node]["strand"+strand].append(genome)
                                                            nodes_dic[node][genome+"_node"] = nodes_count[genome][chromosome]
                                                            nodes_dic[node][genome+"_position"] = position_count[genome][chromosome]["current_position"]
                                                            if position_count[genome][chromosome]["current_position"] < nodes_dic[node]["position_min"]:
                                                                nodes_dic[node]["position_min"] = position_count[genome][chromosome]["current_position"]
                                                            if position_count[genome][chromosome]["current_position"] > nodes_dic[node]["position_max"]:
                                                                nodes_dic[node]["position_max"] = position_count[genome][chromosome]["current_position"]
                                                        else :
                                                            #The node is redundant, so we will check if a node with the same sequence is available. 
                                                            #If not, we will create a new node.
                                                            if ref_node not in ref_nodes_dic:
                                                                ref_nodes_dic[ref_node] = {}
                                                            if genome+"-"+chromosome not in ref_nodes_dic[ref_node] :
                                                                ref_nodes_dic[ref_node][genome+"-"+chromosome] = 2
                                                            else:
                                                                ref_nodes_dic[ref_node][genome+"-"+chromosome] += 1
                                                            node = node + "_" + str(ref_nodes_dic[ref_node][genome+"-"+chromosome])
                                                            if node in nodes_set:
                                                                nodes_dic[node]["genomes"].append(genome)
                                                                nodes_dic[node]["strand"+strand].append(genome)
                                                                nodes_dic[node][genome+"_node"] = nodes_count[genome][chromosome]
                                                                nodes_dic[node][genome+"_position"] = position_count[genome][chromosome]["current_position"]
                                                                if position_count[genome][chromosome]["current_position"] < nodes_dic[node]["position_min"]:
                                                                    nodes_dic[node]["position_min"] = position_count[genome][chromosome]["current_position"]
                                                                if position_count[genome][chromosome]["current_position"] > nodes_dic[node]["position_max"]:
                                                                    nodes_dic[node]["position_max"] = position_count[genome][chromosome]["current_position"]
                                                            else:
                                                                if strand == "P" :
                                                                    nodes_dic[node] = {"genomes":[genome], "max":1, "ref_node" : ref_node, "strandP":[genome], "strandM":[], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                                else :
                                                                    nodes_dic[node] = {"genomes":[genome], "max":1, "ref_node" : ref_node, "strandM":[genome], "strandP":[], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                                nodes_dic[node][genome+"_node"] = nodes_count[genome][chromosome]   
                                                                nodes_dic[node][genome+"_position"] = position_count[genome][chromosome]["current_position"]   
                                                                nodes_dic[ref_node]["max"]+=1
                                                                nodes_set.add(node)
                                                                
                                                nodes_count[genome][chromosome] += 1
                                                position_count[genome][chromosome]["current_position"] += size
                                    bar2.update(1)
                            ligne = file.readline() 
                    nodes_list = None
                    #Flow computing
                    logger.info("\nSize of elements to create into DB : " + str(len(list(nodes_dic.items()))))
                    if len(nodes_dic) > 0:
                        for node in nodes_dic:
                            nodes_dic[node]["flow"] = len(nodes_dic[node]["genomes"])/len(set_all_genomes)
                            position_mean = 0
                            nb_genomes = 0
                            for g in nodes_dic[node]["genomes"]:
                                nb_genomes += 1
                                position_mean += nodes_dic[node][g+"_position"]
                            nodes_dic[node]["position_mean"] = int(position_mean/nb_genomes)
                    
                        logger.info("Time of batch analysis : "+ str(time.time()-temps_0_lot) + "\nTotal time : " + str(time.time()-temps_depart))
                        #Create batch nodes in DB
                        logger.info("Batch " + str(current_batch) + " : Creating nodes in DB")
                        with get_driver() as driver:
                            with driver.session() as session:
                                create_nodes_batch(session, nodes_dic, create=create)
                        nodes_dic = None
                        ref_nodes_dic = None
                        
                        
            create_stats(set_genome, set_chromosome)
            nodes_size_dic = None
            logger.info("Nodes creation is terminated\nTotal time : " + str(time.time()-temps_depart) + "\nGenomes analysed : " + str(set_genome) + "\nNodes number : "+str(total_nodes) )
        
        #Relationships are processed afterwards because node splitting does not guarantee that 
        #relationships will be created between the correct nodes in the event of duplicate nodes. 
        #Relationships will therefore be processed by chromosome.
        logger.info("\nStart relations")
        
        set_relations = set()
        liste_relations = []
        nodes_list = []
        liste_strand = []
        current_batch = 0
        for k in range(index_first_chromosme,len(chromosomes_list)) :
            c = chromosomes_list[k]
            repeat_nodes = {}
            file.seek(0,0)
            ligne = file.readline()
            current_batch += 1
            logger.info("chromosome " + str(current_batch) + " / " + str(len(chromosomes_list)-index_first_chromosme))
            with tqdm(total=total_path) as bar3 :
                while ligne:
                    ligne_dec = ligne.split()
                    if len(ligne_dec) > 0:
                        if ligne[0] == 'P' or ligne[0] == 'W':
                            chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                            if chromosome == c:
                                if ligne_dec[0] == 'P':
                                    ind = 2
                                    walk = 0
                                    nodes_list = re.split(sep[walk], ligne_dec[ind])
                                    if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                        liste_strand = [chaine[-1] for chaine in nodes_list]
                                        nodes_list = [chromosome+"_"+chaine[:-1] for chaine in nodes_list]   
                                    else :
                                        liste_strand = [chaine[-1] for chaine in nodes_list]
                                        nodes_list = [chaine[:-1] for chaine in nodes_list]                           
                                else:
                                    ind = 6
                                    walk = 1
                                    nodes_list = re.split(sep[walk], ligne_dec[ind])
                                    liste_strand = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 != 0]
                                    if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                        nodes_list = [chromosome+"_"+nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                    else:
                                        nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
        
                                if genome != "_MINIGRAPH_":
                                    node = ""
                                    strand = ""
                                    # Browse nodes list : if node already exist for the same sequence
                                    # then create a new node (For example, if it is the sixth iteration for node S1, we will create S1_6)
                                    for i in range(len(nodes_list)):
                                        if i > 0:
                                            previous_node = node
                                        node = nodes_list[i]
                                        if node in repeat_nodes :
                                            if genome+"-"+chromosome in repeat_nodes[node]:
                                                repeat_nodes[node][genome+"-"+chromosome] += 1
                                                node = node+"_"+ str(repeat_nodes[node][genome+"-"+chromosome])
                                            else:
                                                repeat_nodes[node][genome+"-"+chromosome] = 1
                                        else:
                                            repeat_nodes[node] = {genome+"-"+chromosome:1}
                                        if i > 0 and previous_node+"->"+node not in set_relations :
                                            set_relations.add(previous_node+"->"+node)
                                            #liste_relations.append({"depart":noeud_precedent,"arrivee":noeud})
                            bar3.update(1)
                    ligne = file.readline() 

            for rel in set_relations:
                liste_relations.append({"depart":rel.split("->")[0],"arrivee":rel.split("->")[1]})
            set_relations = set()
        
            logger.info("Batch : " + str(current_batch) + " relationships creation, number to create : " + str(len(liste_relations)))
            with get_driver() as driver:
                with driver.session() as session:
                    creer_relations_batch(session, liste_relations)
            liste_relations = []
        logger.info("End of relationships creation\nTime : " + str(time.time()-temps_depart) + "\nNumberof relationships created : " + str(len(set_relations)))
        set_relations = set()
                                                   
    file.close()

    logger.info("treatment completed\nTotal time : "+ str(time.time()-temps_depart))
    return set_genome


def create_node_csv_line(csv_fields_index, dic_node):
    csv_line = [None]*len(csv_fields_index)
    for att in dic_node:
        csv_line[csv_fields_index[att]] = dic_node[att]
    return csv_line

def update_csv_line(csv_fields_index, dic_node, csv_line):
    for att in dic_node:
        if att == "genomes" or att == "strandP" or att == "strandM" :
            csv_line[csv_fields_index[att]].append(dic_node[att])
        else:
            csv_line[csv_fields_index[att]] = dic_node[att]
    return csv_line

#This function replace load_gfa_data_to_neo4j in case of big data : it is very similary but create a dump csv file to import with neo4j-admin import (only for new DB creation)
#The sequences nodes are created too by this function
#After the creation of the csv it is necessary to import them in the database
#Important note : if start chromosome is not None and sequences csv file already exists, the sequences csv file won't be computed
@require_authorization
def load_gfa_data_to_csv(gfa_file_name, import_dir="./data/import", chromosome_file = None, chromosome_prefix = False, batch_size = 2000000, start_chromosome = None, haplotype = True):
    sep = ["[,;.*]", "(<|>)"]
    batch_nb = 0
    set_genome = set()
    set_chromosome = set()
    set_relations = set()
    relations_repeat_nodes = {}
    total_nodes = 0
    total_path = 0
    total_relations = 0
    nodes_size_dic = {}
    temps_depart = time.time()
    set_all_genomes = set()
    set_all_chromosomes = set()
    chromosomes_list = []
    dic_nodes_id = {}
    dic_batch_nodes_index = {}
    ref_nodes_dic = {}
    index_first_chromosme = 0
    node_id = 0
    batch_node_id = 0
    csv_nodes_lines = []
    nodes_set_next_chromosome = set()
    first_chromosome = None
    print_header_nodes = not os.path.isfile(import_dir+"/nodes.csv")
    print_header_relations = not os.path.isfile(import_dir+"/relations.csv")
    print_header_sequences = not os.path.isfile(import_dir+"/sequences.csv")
    print_header_long_sequences = not os.path.isfile(import_dir+"/long_sequences.csv")
    last_node_id = 0
    if chromosome_file is not None:
        c_str = chromosome_file
    else:
        c_str = "all"
    logger.info(f"Create csv import file for {gfa_file_name} - chromosome {c_str}")
    header_file = None
    if os.path.isfile(import_dir+"/nodes.csv") :
        last_line = None
        cpt_lines = 0
        with open(import_dir+"/nodes.csv", mode='r', newline='') as file_node:
            reader = csv.reader(file_node)
            while True:
                try:
                    line = next(reader)
                    cpt_lines +=1
                    if cpt_lines == 1:
                        header_file = line
                    if line:
                        last_line = line
                except StopIteration:
                    break
        
        if last_line:
            last_node_id = int(last_line[0])+1
            
    csv_sequence_file = open(import_dir+"/sequences.csv", "a", newline="", encoding="utf-8") 
    sequences_writer = csv.writer(csv_sequence_file)
    if print_header_sequences:
        sequences_writer.writerow([":ID", "name", "sequence:STRING"])
    file = open(gfa_file_name, "r", encoding='utf-8')
    
    
    if header_file is not None :
        pos = 0
        csv_fields_index = {}
        for h in header_file:
            header_name = h.split(":")[0]
            if header_name == "" :
                csv_fields_index["id"] = pos
            else:
                csv_fields_index[header_name] = pos
            pos += 1
    else:
        csv_fields_index = {"id":0,"name":1, "max":2, "ref_node":3, "size" : 4, "chromosome"  : 5, "position_min":6, "position_max":7, "genomes":8, "strandP":9, "strandM": 10, "position_mean" : 11, "flow" : 12}
    #logger.debug("csv field index : " + str(csv_fields_index))
    with file:
        
        #First file browsing to get length, nodes and haplotypes
        ligne = file.readline()
        max_size_S_line = 0
        while ligne :
            if ligne.startswith(('S')):
                if len(ligne) > max_size_S_line:
                    max_size_S_line = len(ligne)
                ligne_dec = ligne.split()
                if (len(ligne_dec) > 0): 
                    if chromosome_file is not None and chromosome_file != "": 
                        node_name = chromosome_file + "_" + ligne_dec[1]
                    else :
                        node_name = ligne_dec[1]
                    nodes_size_dic[ligne_dec[1]]=int(len(ligne_dec[2]))
                    if start_chromosome is None or print_header_sequences :
                        sequences_writer.writerow([last_node_id, node_name, ligne_dec[2]])
                        last_node_id += 1
                    total_nodes += 1
            if ligne.startswith(('P',"W")):
                logger.debug(ligne[0:80])
                total_path += 1
                ligne_dec = ligne.split()
                chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                if chromosome not in set_all_chromosomes :
                    chromosomes_list.append(chromosome)
                if len(set_all_chromosomes) == 0:
                    first_chromosome = chromosome
                    if start_chromosome is not None and start_chromosome != "":
                        first_chromosome = start_chromosome
                set_all_genomes.add(genome) 
                set_all_chromosomes.add(chromosome)
                if chromosome == first_chromosome :
                    if ligne_dec[0] == 'P':
                        ind = 2
                        walk = 0
                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([chaine[:-1] for chaine in nodes_list])
                    else:
                        ind = 6
                        walk = 1
                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                        nodes_set_next_chromosome |= set([nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0])
                
            ligne = file.readline() 
        last_index = len(csv_fields_index)
        max_size_S_line += 32
        if max_size_S_line > get_conf_read_buffer_size() and max_size_S_line < MAX_READ_BUFFER_SIZE_VALUE:
            set_conf_value("read_buffer_size", max_size_S_line)
            logger.info(f"Set read buffer size conf to {max_size_S_line}")
        node_id = last_node_id
        if print_header_nodes :
            csv_header_node = [":ID", "name:STRING", "max:LONG","ref_node:STRING", "size:LONG", "chromosome:STRING", "position_min:LONG", "position_max:LONG", "genomes:STRING[]","strandP:STRING[]", "strandM:STRING[]", "position_mean:LONG", "flow:DOUBLE"]
            for g in sorted(set_all_genomes):
                csv_fields_index[g+"_position"] = last_index
                csv_header_node.append(g+"_position:LONG")
                csv_fields_index[g+"_node"] = last_index + 1
                csv_header_node.append(g+"_node:LONG")
                last_index += 2
        csv_nodes_file = open(import_dir+"/nodes.csv", "a", newline="", encoding="utf-8") 
        nodes_writer = csv.writer(csv_nodes_file)
        if print_header_nodes:
            nodes_writer.writerow(csv_header_node)
        csv_relations_file = open(import_dir+"/relations.csv", "a", newline="", encoding="utf-8") 
        relations_writer = csv.writer(csv_relations_file)      
        if print_header_relations:
            relations_writer.writerow([":START_ID",":END_ID", ":TYPE"])
            
        if first_chromosome is not None and first_chromosome != "":
            for k in range(len(chromosomes_list)):
                if chromosomes_list[k] == first_chromosome:
                    index_first_chromosme = k
        logger.info("Genomes number : " + str(len(set_all_genomes)) + " - genomes list : " + str(set_all_genomes) +  "- chromosomes list : " + str(set_all_chromosomes))
        logger.info("Start parsing, nodes number : " + str(total_nodes) + "\nstart chromosome : " + str(start_chromosome) + "\nstart index : " + str(node_id))
        for k in range(index_first_chromosme,len(chromosomes_list)) :
            c = chromosomes_list[k]
            relations_repeat_nodes = {}
            dic_nodes_id = {}
            total_relations += len(set_relations)
            set_relations = set()
            nodes_set_chromosome = set(nodes_set_next_chromosome)
            nodes_set_next_chromosome = set()
            logger.info("chromosome " + str(c) + " - number of nodes : " + str(len(nodes_set_chromosome)))
            batch_nb = ceil(len(nodes_set_chromosome)/batch_size)
            current_batch = 0
            while current_batch < batch_nb :
                temps_0_lot = time.time()
                
                nodes_batch_set = set(list(nodes_set_chromosome)[current_batch*batch_size:min(len(nodes_set_chromosome),(current_batch+1)*batch_size)])
                current_batch += 1
                logger.info("chromosome " + c + " batch " + str(current_batch) + "/"+str(batch_nb) + " nodes number : " + str(len(nodes_batch_set)))
                #Relations for the chromosome are computed only on the last batch
                if current_batch == batch_nb:
                    compute_relations_batch = True
                else:
                    compute_relations_batch = False
                
                file.seek(0,0)
                ligne = file.readline()

                #Path browsing for the batch 
                dic_batch_nodes_index = {}
                ref_nodes_dic = {}
                nodes_list = []
                liste_strand = []
                position_count = {}
                nodes_count = {}
                set_genomes_lot = set()    
                nodes_set = set()
                with tqdm(total=total_path) as bar2 :
                    while ligne:
                        ligne_dec = ligne.split()
                        if len(ligne_dec) > 0:
                            if ligne[0] == 'P' or ligne[0] == 'W':
                                chromosome, genome = get_chromosome_genome(ligne, haplotype = haplotype, chromosome_file=chromosome_file)
                                ligne = None
                                if current_batch == batch_nb and k < len(chromosomes_list) - 1 and chromosome == chromosomes_list[k+1]:
                                    #last batch for the chromosome, retrieves the nodes to be processed for the next chromosome
                                    if ligne_dec[0] == 'P':
                                        ind = 2
                                        walk = 0
                                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                                        nodes_set_next_chromosome |= set([chaine[:-1] for chaine in nodes_list])
                                    else:
                                        ind = 6
                                        walk = 1
                                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                                        nodes_set_next_chromosome |= set([nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0])
                                
                                if chromosome == c:
                                    if ligne_dec[0] == 'P':
                                        ind = 2
                                        walk = 0
                                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                                        liste_strand = [chaine[-1] for chaine in nodes_list]
                                        nodes_list = [chaine[:-1] for chaine in nodes_list]
                                        # if chromosome_file == None :
                                        #     liste_strand = [chaine[-1] for chaine in nodes_list]
                                        #     nodes_list = [chaine[:-1] for chaine in nodes_list]
                                        # else :
                                        #     liste_strand = [chaine[-1] for chaine in nodes_list]
                                        #     nodes_list = [chromosome_file+"_"+chaine[:-1] for chaine in nodes_list]                              
                                    else:
                                        ind = 6
                                        walk = 1
                                        nodes_list = re.split(sep[walk], ligne_dec[ind])
                                        liste_strand = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 != 0]
                                        nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                        # if chromosome_file == None :
                                        #     nodes_list = [nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
                                        # else:
                                        #     nodes_list = [chromosome_file+"_"+nodes_list[j] for j in range(1, len(nodes_list)) if j % 2 == 0]
      
                                    if chromosome is not None and chromosome not in set_chromosome:
                                        set_chromosome.add(chromosome)
                                    
                                    if genome != "_MINIGRAPH_":
                                        if genome not in set_genome:
                                            set_genome.add(genome)
                                        if genome not in set_genomes_lot :
                                            set_genomes_lot.add(genome)
                                            nodes_count[genome] = {}
                                            position_count[genome] = {}
                                        if chromosome not in nodes_count[genome] :
                                            nodes_count[genome][chromosome] = 0
                                            position_count[genome][chromosome] = {}
                                            position_count[genome][chromosome]["current_position"] = 0
                                            position_count[genome][chromosome]["previous_position"] = 0
                                            position_count[genome][chromosome]["current_contig"] = ""
                                        #For walk, start position is available => this position is used
                                        if ind == 6:
                                            if position_count[genome][chromosome]["current_contig"] != ligne_dec[3] :
                                                #New contig => add the start of next contig
                                                position_count[genome][chromosome]["current_position"] += int(ligne_dec[4])
                                            else :
                                                #Same contig => add the potential gaps
                                                if position_count[genome][chromosome]["previous_position"] - int(ligne_dec[4]) > 0 :
                                                    position_count[genome][chromosome]["current_position"] += position_count[genome][chromosome]["previous_position"] - int(ligne_dec[4])
                                            position_count[genome][chromosome]["current_contig"] = ligne_dec[3]
                                            position_count[genome][chromosome]["previous_position"] = int(ligne_dec[5])
    
                                        node = ""
                                        relation_node = ""
                                        previous_node_relation = ""
                                        ref_node = ""
                                        prefix_ref_node = ""
                                        strand = ""
                                        ligne_dec=None
                                        dic_node_update = {}
                                        #Graph linearization
                                        # Browse nodes list : if node already exist for the same sequence
                                        # then create a new node (For example, if it is the sixth iteration for node S1, we will create S1_6)
                                        for i in range(len(nodes_list)):
                                            if i > 0:
                                                previous_node = node
                                            node = nodes_list[i]
                                            ref_node = node
                                            size = nodes_size_dic[ref_node]
                                            if chromosome_prefix or (chromosome_file is not None and chromosome_file != ""):
                                                node = chromosome + "_" + node
                                            prefix_ref_node = node
                                            
                                            #Node is consider only if it is part of batch
                                            if ref_node in nodes_batch_set:
                                                if chromosome_file is not None and chromosome_file != "":
                                                    ref_node = node
                                                strand = ""
                                                if (i < len(liste_strand) and liste_strand[i] in ["-", "<"]):
                                                    strand = "M"
                                                else:
                                                    strand = "P"
                                                
                                                if node not in nodes_set:
                                                    if strand == "P" :
                                                        dic_node =  {"id":node_id,"name":node,"genomes":[genome], "max":1, "strandP":[genome], "strandM":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                    else :
                                                        dic_node =  {"id":node_id,"name":node,"genomes":[genome], "max":1, "strandM":[genome], "strandP":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                    csv_nodes_lines.append(create_node_csv_line(csv_fields_index, dic_node))
                                                    dic_batch_nodes_index[node] = batch_node_id
                                                    dic_nodes_id[node] = node_id
                                                    node_id += 1
                                                    batch_node_id += 1
                                                    nodes_set.add(node)
                                                else:
                                                    tmp_node_id = dic_batch_nodes_index[node]
                                                    tmp_genomes = csv_nodes_lines[tmp_node_id][csv_fields_index["genomes"]]
                                                    if genome not in tmp_genomes : 
                                                        dic_node_update = {"genomes":genome, 
                                                                           "strand"+strand:genome,
                                                                           genome+"_node": nodes_count[genome][chromosome],
                                                                           genome+"_position": position_count[genome][chromosome]["current_position"]
                                                                           }
                                                        tmp_position_min = csv_nodes_lines[tmp_node_id][csv_fields_index["position_min"]]
                                                        tmp_position_max = csv_nodes_lines[tmp_node_id][csv_fields_index["position_max"]]
                                                        if position_count[genome][chromosome]["current_position"] < tmp_position_min:
                                                            dic_node_update["position_min"]= position_count[genome][chromosome]["current_position"]
                                                        if position_count[genome][chromosome]["current_position"] > tmp_position_max:
                                                            dic_node_update["position_max"] = position_count[genome][chromosome]["current_position"]
                                                        csv_nodes_lines[tmp_node_id]=update_csv_line(csv_fields_index, dic_node_update, csv_nodes_lines[tmp_node_id])
                                                    
                                                    else:
                                                        if ref_node not in ref_nodes_dic:
                                                            ref_nodes_dic[ref_node] = {}
                                                        if genome+"-"+chromosome not in ref_nodes_dic[ref_node] :
                                                            ref_nodes_dic[ref_node][genome+"-"+chromosome] = 2
                                                        else:
                                                            ref_nodes_dic[ref_node][genome+"-"+chromosome] += 1
                                                        node = node + "_" + str(ref_nodes_dic[ref_node][genome+"-"+chromosome])

                                                        if node in nodes_set:
                                                            tmp_node_id = dic_batch_nodes_index[node]
                                                            dic_node_update = {"genomes":genome, 
                                                                               "strand"+strand:genome,
                                                                               genome+"_node": nodes_count[genome][chromosome],
                                                                               genome+"_position": position_count[genome][chromosome]["current_position"]
                                                                               }
                                                            tmp_position_min = csv_nodes_lines[tmp_node_id][csv_fields_index["position_min"]]
                                                            tmp_position_max = csv_nodes_lines[tmp_node_id][csv_fields_index["position_max"]]
                                                            if position_count[genome][chromosome]["current_position"] < tmp_position_min:
                                                                dic_node_update["position_min"]= position_count[genome][chromosome]["current_position"]
                                                            if position_count[genome][chromosome]["current_position"] > tmp_position_max:
                                                                dic_node_update["position_max"] = position_count[genome][chromosome]["current_position"]
                                                            csv_nodes_lines[tmp_node_id]=update_csv_line(csv_fields_index, dic_node_update, csv_nodes_lines[tmp_node_id])
                                                        else:
                                                            if strand == "P" :
                                                                dic_node =  {"id":node_id,"name":node,"genomes":[genome], "max":1, "strandP":[genome], "strandM":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                            else :
                                                                dic_node =  {"id":node_id,"name":node,"genomes":[genome], "max":1, "strandM":[genome], "strandP":[], "ref_node" : ref_node, genome+"_node":nodes_count[genome][chromosome],genome+"_position":position_count[genome][chromosome]["current_position"], "size" : size, "chromosome"  : chromosome, "position_min":position_count[genome][chromosome]["current_position"], "position_max":position_count[genome][chromosome]["current_position"]}
                                                            csv_nodes_lines.append(create_node_csv_line(csv_fields_index, dic_node))
                                                            dic_batch_nodes_index[node] = batch_node_id
                                                            dic_nodes_id[node] = node_id
                                                            node_id += 1
                                                            batch_node_id += 1
                                                            nodes_set.add(node) 
                                                            #Update max on the ref node
                                                            tmp_node_ref_id = dic_batch_nodes_index[prefix_ref_node]
                                                            tmp_max = csv_nodes_lines[tmp_node_ref_id][csv_fields_index["max"]]
                                                            csv_nodes_lines[tmp_node_ref_id]=update_csv_line(csv_fields_index, {"max":tmp_max+1}, csv_nodes_lines[tmp_node_ref_id])

                                            nodes_count[genome][chromosome] += 1
                                            position_count[genome][chromosome]["current_position"] += size
                                            #compute relations
                                            if compute_relations_batch:
                                                if i > 0 :
                                                    previous_node_relation = node_relation
                                                if node in relations_repeat_nodes :
                                                    if genome+"-"+chromosome in relations_repeat_nodes[node]:
                                                        relations_repeat_nodes[node][genome+"-"+chromosome] += 1
                                                        node_relation = node+"_"+ str(relations_repeat_nodes[node][genome+"-"+chromosome])
                                                    else:
                                                        relations_repeat_nodes[node][genome+"-"+chromosome] = 1
                                                        node_relation = node
                                                else:
                                                    relations_repeat_nodes[node] = {genome+"-"+chromosome:1}
                                                    node_relation = node
                                                if i > 0 and previous_node_relation+"->"+node_relation not in set_relations :
                                                    set_relations.add(previous_node_relation+"->"+node_relation)
                                                    relations_writer.writerow([dic_nodes_id[previous_node_relation], dic_nodes_id[node_relation], "gfa_link"])
                                bar2.update(1)
                        ligne = file.readline() 
                nodes_list = None
                #Flow computing
                logger.info("\nSize of elements to create into csv : " + str(len(csv_nodes_lines)))
                for line in csv_nodes_lines:
                    line[csv_fields_index["flow"]] = len(line[csv_fields_index["genomes"]])/len(set_all_genomes) if len(set_all_genomes) > 0 else 0
                    for idx in [csv_fields_index["genomes"],csv_fields_index["strandP"], csv_fields_index["strandM"]]:
                        line[idx] = ";".join(line[idx])
   
                    values_position = [v for v in line[csv_fields_index["flow"]+1:-1:2] if v is not None]
                    values_node = [v for v in line[csv_fields_index["flow"]+2:-1:2] if v is not None]
                    line[csv_fields_index["position_mean"]] = int(sum(values_position) / len(values_position)) if len(values_position) > 0 else 0
                nodes_writer.writerows(csv_nodes_lines)
                csv_nodes_lines = []
                batch_node_id = 0
                dic_batch_nodes_index = {}

                    
                    
        nodes_size_dic = None
        logger.info("Csv creation is terminated\nTotal time : " + str(time.time()-temps_depart) + "\nGenomes analysed : " + str(set_genome) + "\nNodes number : "+str(total_nodes) +"\nRelations number : " + str(total_relations))
        
                                                
    file.close()
    logger.info("Treatment completed\nTotal time : "+ str(time.time()-temps_depart))
    return set_genome, chromosomes_list



def get_chromosome_annotation(annotation):
    chromosome = None
    tab = re.split(r'[_,.#]', annotation)
    for i in range(len(tab)):
        t = tab[i]
        if t.upper().startswith(("CHR")):
            if len(t) > 3:
                chromosome = t[3:].lstrip('0')
            else:
                if i < len(tab)-1:
                    chromosome = tab[i+1].lstrip('0')
        else:
            chromosome = annotation.split()[0].upper().replace("CHR", "").replace("C", "").lstrip('0')       
    return chromosome
    


#This function will create the Annotation nodes in the neo4j database from a gtf file, but without creating the relationships.
@require_authorization
def load_annotations_neo4j(annotations_file_name, genome_ref, node_name="Annotation", single_chromosome = None):
    temps_depart = time.time()
    file = open(annotations_file_name, "r", encoding='utf-8')
    nodes_dic = {}
    file_format = "gtf"
    current_gene_id = ""
    file_name = os.path.basename(annotations_file_name)
    with file:
        n = 0
        for line in file :
            n +=1
        file.seek(0,0)
        with tqdm(total=n) as bar :
            #Creating annotations nodes
            ligne = file.readline()
            while ligne :
                if ligne[0] != '#':
                    ligne_dec = ligne.split("\t")
                    chromosome = get_chromosome_annotation(ligne_dec[0])
                    if single_chromosome == None or single_chromosome == chromosome : 
                        name = hashlib.sha256(ligne.encode("utf-8")).hexdigest()
                        nodes_dic[name] = {}
                        
                        feature = ligne_dec[2].lower()
                        nodes_dic[name]["name"] = name
                        nodes_dic[name]["chromosome"] = chromosome
                        nodes_dic[name]["genome_ref"] = genome_ref
                        nodes_dic[name]["source"] = ligne_dec[1]
                        nodes_dic[name]["feature"] = feature
                        nodes_dic[name]["filename"] = file_name
                        
                        start = int(ligne_dec[3])
                        end = int(ligne_dec[4])
                        nodes_dic[name]["start"] = start
                        nodes_dic[name]["end"] = end
                        strand = ligne_dec[6]
                        frame = ligne_dec[7]
                        if len(ligne_dec) == 9:
                            attributes = re.split(r"[;=]", ligne_dec[8])
                            file_format = "gff"
                        else:
                            attributes = ligne_dec[8:]
                        if file_format == "gtf":
                            for i in range(len(attributes)):
                                match attributes[i].lower() :
                                    case "gene_id":
                                        nodes_dic[name]["gene_id"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_version":
                                        nodes_dic[name]["gene_version"] = attributes[i+1][:-1].replace('"',"")
                                    case "transcript_id":
                                        nodes_dic[name]["transcript_id"] = attributes[i+1][:-1].replace('"',"")
                                    case "transcript_version":
                                        nodes_dic[name]["transcript_version"] = attributes[i+1][:-1].replace('"',"")
                                    case "exon_number":
                                        nodes_dic[name]["exon_number"] = attributes[i+1][:-1].replace('"',"")
                                    case "gene_name":
                                        nodes_dic[name]["gene_name"] = attributes[i+1][:-1].replace('"',"")
                                    case "full_gene_name":
                                        nodes_dic[name]["full_gene_name"] = attributes[i + 1][:-1].replace('"', "")
                                    case "gene_source":
                                        nodes_dic[name]["gene_source"] = attributes[i+1][:-1].replace('"',"")    
                                    case "gene_biotype":
                                        nodes_dic[name]["gene_biotype"] = attributes[i+1][:-1].replace('"',"")     
                                    case "transcript_name":
                                        nodes_dic[name]["transcript_name"] = attributes[i+1][:-1].replace('"',"")  
                                    case "transcript_source":
                                        nodes_dic[name]["transcript_source"] = attributes[i+1][:-1].replace('"',"")  
                                    case "transcript_biotype":
                                        nodes_dic[name]["transcript_biotype"] = attributes[i+1][:-1].replace('"',"") 
                                    case "protein_id":
                                        nodes_dic[name]["protein_id"] = attributes[i+1][:-1].replace('"',"") 
                                    case "protein_version":
                                        nodes_dic[name]["protein_version"] = attributes[i+1][:-1].replace('"',"") 
                                    case "tag":
                                        nodes_dic[name]["tag"] = attributes[i+1][:-1].replace('"',"") 
                        else:
                            if feature == "gene" :
                                current_gene_id = ""
                                current_gene_name = ""
                                current_full_gene_name = ""
                                current_transcript_id = ""
                                current_transcript_name = ""
                            if feature == "exon":
                                exon_id = ""
                                exon_number = ""
                                
                            for i in range(len(attributes)):
                                current_attribute = attributes[i].lower()
                                if feature == "gene" and current_attribute =="id":
                                    current_gene_id = attributes[i+1]
                                elif feature == "gene" and current_attribute =="gene_id":
                                    current_gene_id = attributes[i+1]
                                elif feature == "gene" and current_attribute =="name":
                                    current_gene_name = attributes[i+1]
                                elif feature == "gene" and current_attribute =="full_name":
                                    current_full_gene_name = attributes[i+1]
                                elif feature == "mrna" and current_attribute == "id":
                                    current_transcript_id = attributes[i+1]
                                elif feature == "mrna" and current_attribute == "transcript_id":
                                    current_transcript_id = attributes[i+1]
                                elif feature == "mrna" and current_attribute == "name":
                                    current_transcript_name = attributes[i+1]
                                elif feature == "exon" and current_attribute == "id":
                                    exon_id = attributes[i+1]
                                    exon_number = exon_id.split(".")[-1]
                                    nodes_dic[name]["exon_id"] = exon_id
                                    nodes_dic[name]["exon_number"] = exon_number 
                                else :
                                    if feature == None or feature == "":
                                        feature = "unknown_feature"
                                    if current_attribute == "id":
                                        nodes_dic[name][feature+"_id"] = attributes[i+1]
                                    elif current_attribute == "name":
                                        nodes_dic[name][feature+"_name"] = attributes[i+1]
                                        
                            if current_gene_id != "" :
                                if current_gene_name is None or current_gene_name == "":
                                    current_gene_name = current_gene_id
                                nodes_dic[name]["gene_id"] = current_gene_id
                                nodes_dic[name]["gene_name"] = current_gene_name
                                nodes_dic[name]["full_gene_name"] = current_full_gene_name
                                if current_transcript_id != "":
                                    nodes_dic[name]["transcript_id"] = current_transcript_id
                                    nodes_dic[name]["transcript_name"] = current_transcript_name
                            
                                    
                            
                
                bar.update(1)
                ligne = file.readline()

                    
           
            logger.info(f"\nAnnotation analyzing terminated in {time.time()-temps_depart} s.")
            logger.info(f"Creating {len(list(nodes_dic.items()))} nodes in database...")

            with get_driver() as driver:
                with driver.session() as session:
                    create_nodes_batch(session, nodes_dic, node_name=node_name)

        logger.info("Nodes created\nTime : " + str(time.time()-temps_depart))
        #Fin des lots, on créé toutes les relations                                               
    file.close()
    

    logger.info("Annotation load terminated.\nTotal time : "+ str(time.time()-temps_depart))

@require_authorization
def process_annotation_simple_batch(tx, annotations, genome_ref):
    # query = f"""
    #     UNWIND $annotations AS annot
    #     WITH annot.chromosome AS chr, annot.start AS start, annot.end AS end, collect(annot) AS group_annot
    #
    #     MATCH (n1:Node)
    #     WHERE n1.chromosome = chr
    #       AND n1.`{genome_ref}_position` >= start
    #       AND n1.`{genome_ref}_position` <= end
    #
    #     WITH group_annot, n1, chr, start
    #     UNWIND group_annot AS ann1
    #     MATCH (a1:Annotation {{name: ann1.name}})
    #     MERGE (n1)-[:annotation_link]->(a1)
    # """

    # query = f"""
    #     UNWIND $annotations AS annot

    #     MATCH (n1:Node)
    #     WHERE n1.chromosome = annot.chromosome
    #       AND n1.`{genome_ref}_position` >= annot.start
    #       AND n1.`{genome_ref}_position` <= annot.end

    #     MATCH (a1:Annotation {{name: annot.name}})
    #     MERGE (n1)-[:annotation_link]->(a1)
    # """

    query = f"""
        CALL apoc.periodic.iterate(
          "
          UNWIND $annotations AS annot
          MATCH (n1:Node)
          WHERE n1.chromosome = annot.chromosome
            AND n1.`{genome_ref}_position` >= annot.start
            AND n1.`{genome_ref}_position` <= annot.end
          MATCH (a1:Annotation {{name: annot.name}})
          RETURN n1, a1
          ",
          "
          MERGE (n1)-[:annotation_link]->(a1)
          ",
          {{
            batchSize: 1000,
            parallel: false,
            params: {{annotations: $annotations}}
          }}
        )
        """

    #logger.debug(query)

    tx.run(query, annotations=annotations)

@require_authorization    
def process_annotation_complex_batch(tx, annotations, genome_ref, annotation_search_limit=10000):
    

    query = f"""
        UNWIND $annotations AS annot
        
        CALL {{
            WITH annot
            MATCH (n1:Node)
            WHERE n1.chromosome = annot.chromosome
              AND n1.`{genome_ref}_position` < annot.start
              AND n1.`{genome_ref}_position` >= annot.start - {annotation_search_limit}
            ORDER BY n1.`{genome_ref}_position` DESC
            LIMIT 1
            RETURN n1
        }}
        
        WITH annot, n1
        WHERE n1.`{genome_ref}_position` + n1.size >= annot.start
        
        MATCH (a1:Annotation {{name: annot.name}})
        MERGE (n1)-[:annotation_link]->(a1)
    """

                
    #logger.debug(query)
    tx.run(query, annotations=annotations)


    
    
@require_authorization    
def process_annotation_last_complex_batch(tx, genome_ref, annotation_search_limit=10000, batch_limit=10000):

    query="""
        MATCH (n:Node)
        WHERE n.size > $annotation_search_limit
        AND n.`{genome_ref}_position` is not null
        return count(n) as nodes_number
    """
    result = tx.run(query,annotation_search_limit=annotation_search_limit) 
    record = result.single()
    nodes_nb = record["nodes_number"]
    
    current_nodes = 0
    position_field = "`"+genome_ref+"_position`"
    genome_ref_field = "`"+genome_ref+"`"
    while current_nodes * batch_limit < nodes_nb :
        offset = current_nodes * batch_limit
        query=f"""
            MATCH (n:Node)
            WHERE n.size > $annotation_search_limit
            AND n.`{genome_ref}_position` is not null
            SKIP {offset} LIMIT {batch_limit}
            WITH n
            MATCH (a:Annotation)
            WHERE a.chromosome = n.chromosome
            AND a.genome_ref = n.{genome_ref_field}
            AND a.start > n.{position_field}
            AND a.start < n.{position_field} + n.size
            MERGE (n)-[:annotation_link]->(a)
        """
        current_nodes += 1
        result = tx.run(query,annotation_search_limit=annotation_search_limit) 
    

#This function will create relationships between nodes and annotations in the database. 
#processing is divided into two types of relationships:   
#- Simple relationships, where the position of a node (for a reference genome) lies between the start and end of an annotation.
#- Complex relationships, where the start/end of an annotation lies between the start/end of a node
#Parameters : 
#   - genome_ref : only nodes on genome_ref will be analyzed, is None all annotations will be computed
@require_authorization
def creer_relations_annotations_neo4j(genome_ref=None, chromosome=None):
    temps_depart = time.time()
    with get_driver() as driver:
        last_id = -1
        batch_size = 10000
        total_annotations = 0
        WARN = False
        WARN_message = ""
        with driver.session() as session:
            #Processing simple annotations: those for which the start is between the start and end of a node
            #this is the largest volume of annotations
            logger.info("Processing annotations")
            #Get the chromosomes present in the pangenome
            query = """
            MATCH (s:Stats) 
            RETURN s.chromosomes as all_chromosomes
            """
            graph_chromosomes_set = set()
            result = session.run(query)
            for record in result:
                graph_chromosomes_set = set(record["all_chromosomes"])

            i = 0
            last_name = None
            last_id = -1
            max_id = batch_size
            min_id = -1
            current_id = -1
            annotations_nb = 0
            all_genomes = set()
            if genome_ref is None or genome_ref == "" : 
                query = """
                    MATCH (a:Annotation)
                    RETURN collect(DISTINCT(a.genome_ref)) AS liste_genomes
                """
                result = session.run(query, min_id=current_id,max_id=current_id+batch_size)
                for record in result:
                    liste_genomes = record["liste_genomes"]
            else :
                liste_genomes = [genome_ref]

            logger.info("Haplotypes list : " + str(liste_genomes))

            for g in liste_genomes :
                logger.info(f"Linking annotation for genome {g}")
                query = f"""
                    MATCH (a:Annotation) where a.genome_ref = "{g}" return collect(distinct(a.chromosome)) as annotations_chromosomes
                """
                result = session.run(query)
                for record in result:
                    annotations_chromosomes_set = set(record["annotations_chromosomes"])
                query = f"""
                MATCH (a:Annotation) where a.genome_ref = "{g}" return min(ID(a)) as min_id, max(ID(a)) as max_id, count(a) as annotations_count
                """
                #logger.debug(query)
                result = session.run(query)
                for record in result:
                    min_id = record["min_id"]
                    max_id = record["max_id"]
                    annotations_count = record["annotations_count"]

                intersection_chromosomes_set = graph_chromosomes_set & annotations_chromosomes_set
                if len(intersection_chromosomes_set) == 0:
                    WARN = True
                    continue
                else:
                    WARN = False
                if min_id is None or max_id is None or min_id == max_id :
                    continue
                batch_number = ceil((max_id-min_id)/batch_size)
                current_id = min_id
                all_genomes.add(g)
                i = 0
                with tqdm(total=annotations_count, desc=f"Haplotype {g}") as pbar:
                    while current_id < max_id:
                        i+=1
                        #logger.debug("Batch nb " + str(i) + "/" + str(batch_number) + " Current id : " + str(current_id) + " max id : " + str(max_id) + " - haplotype : " + str(g))
                        annotations_nb = 0
                        if chromosome is None :
                            annotations = session.run(
                                """
                                MATCH (a:Annotation)
                                WHERE ID(a) >= $min_id AND ID(a) < $max_id AND a.genome_ref = $genome
                                RETURN a.name AS name, a.chromosome AS chromosome, a.start AS start, a.end AS end
                                """,
                                min_id=current_id,
                                max_id=min(max_id,current_id+batch_size),
                                genome=g
                                ).data()
                        else:
                            annotations = session.run(
                                """
                                MATCH (a:Annotation)
                                WHERE ID(a) >= $min_id AND ID(a) < $max_id and a.chromosome = $chromosome AND a.genome_ref = $genome
                                RETURN a.name AS name, a.chromosome AS chromosome, a.start AS start, a.end AS end
                                """,
                                min_id=current_id,
                                max_id=min(max_id,current_id+batch_size),
                                chromosome=chromosome,
                                genome=g
                                ).data()
                        annotations_nb += len(annotations)
                        total_annotations += len(annotations)
                        if annotations and len(annotations) > 0:
                            with session.begin_transaction() as tx:
                                #logger.info("creating annotations batch " + str(i) + "/"+str(batch_number))
                                process_annotation_simple_batch(tx,annotations, g)
                                #logger.info("creating complexe annotations")
                                #Handling complex annotations: those for which the start and end of a node are before and after the annotation
                                #the volume is much lower (less than 1%)
                                process_annotation_complex_batch(tx, annotations, g, annotation_search_limit=10000)
                                tx.commit()

                        #logger.debug("Annotations nb : " + str(annotations_nb) + " - Annotations already treated : " + str(total_annotations))
                        pbar.update(annotations_nb)
                        current_id += batch_size
            
            
            for g in all_genomes:
                with session.begin_transaction() as tx:
                    logger.info(f"processing complex annotations for genome {g}")
                    process_annotation_last_complex_batch(tx, g, annotation_search_limit=10000)
                    tx.commit()
    if WARN:
        WARN_message = "Chromosome names mismatch between annotation file and graph."
    logger.info(f"End of relationships creation, {total_annotations} annotations analysed. {WARN_message} Total time : " + str(time.time()-temps_depart))
    return WARN
    


#Main function to construct the whole db from the gfa file
#If the gfa relate to a single chromosome, chromosome_file must contains the reference of this chromosome (1, 2, X, Y, etc.)
#batch_size value is important to limit memory usage, according to the memory available it can be necessary to reduce this value for big pangenomes graphs.
#genome_ref is required if an annotation_file_name is present : this name is used to link the annotations nodes with the main nodes of the graph.
@require_authorization
def construct_DB(gfa_file_name, annotation_file_name = None, genome_ref = None, chromosome_file = None, chromosome_prefix = False, batch_size = 2000000, start_chromosome = None, create = False, haplotype = True, create_only_relations = False):
    start_time = time.time()
    WARN = False
    load_sequences(gfa_file_name, chromosome_file, create=create)
    sequence_time = time.time()
    logger.info("Sequences loaded in " + str(sequence_time-start_time) + " s")
    load_gfa_data_to_neo4j(gfa_file_name, chromosome_file = chromosome_file,  chromosome_prefix = chromosome_prefix, batch_size = batch_size, create = create, start_chromosome = start_chromosome, haplotype=haplotype, create_only_relations=create_only_relations)
    graph_time = time.time()
    logger.info("Graph loaded in " + str(graph_time-sequence_time) + " s")
    create_indexes(base=False, extend=True, genomes_index=True)
    index_time = time.time()
    logger.info("Indexes created in " + str(index_time-graph_time) + " s")
    if annotation_file_name != None and genome_ref != None :
        load_annotations_neo4j(annotation_file_name, genome_ref = genome_ref, single_chromosome = chromosome_file)
        annotation_time = time.time()
        logger.info("Annotations loaded in " + str(annotation_time-index_time) + " s")
        WARN = creer_relations_annotations_neo4j(genome_ref)
        annotation_relation_time = time.time()
        logger.info("Annotations relations loaded in " + str(annotation_relation_time-annotation_time) + " s")
    logger.info("Process terminated. BDD construct in " + str(time.time()-start_time) + " s")
    return WARN

#This function load multiple gfa files : each file must relate to a single chromosome
#The files must be named so that last character before .gfa extension contains the reference of the chromosome
#Exmples : chr1.gfa, chr01.gfa, chromosome_1.gfa, exemple_chr_X.gfa, etc.
@require_authorization
def construct_db_by_chromosome(gfa_chromosomes_dir, annotation_file_name = None, genome_ref = None, chromosome_file = None, start_node = 0, batch_size = 5000000, create=False):
    start_time = time.time()
    WARN = False
    for gfa_file_name in  os.listdir(gfa_chromosomes_dir):
        if gfa_file_name.endswith(".gfa"):
            if "_" in gfa_file_name:
                chromosome = gfa_file_name[:-4].split("_")[-1]
            else:
                chromosome = gfa_file_name[:-4]

            chromosome = chromosome.lower().removeprefix("chr")
            chromosome = chromosome.lstrip("0")
            logger.info("Loading chromosome : " + str(chromosome))
            if chromosome != "" :
                load_sequences(gfa_file_name, chromosome_file=chromosome, create=create)
                load_gfa_data_to_neo4j(gfa_file_name, chromosome_file = chromosome_file, batch_size = batch_size, start_node = start_node, create = create)
    db_time = time.time()
    logger.info("Sequences loaded in " + str(db_time-start_time) + " s")
    create_indexes(base=False, extend=True, genomes_index=True)
    index_time = time.time()
    logger.info("Indexes created in " + str(index_time-graph_time) + " s")
    if annotation_file_name != None and genome_ref != None :
        load_annotations_neo4j(annotation_file_name, genome_ref = genome_ref, single_chromosome = chromosome_file)
        annotation_time = time.time()
        logger.info("Annotations loaded in " + str(annotation_time-index_time) + " s")
        WARN = creer_relations_annotations_neo4j(genome_ref)
        annotation_relation_time = time.time()
        logger.info("Annotations relations loaded in " + str(annotation_relation_time-annotation_time) + " s")
    logger.info("Process terminated. BDD construct in " + str(time.time()-start_time) + " s")
    return WARN

@require_authorization
def delete_annotations(batch_size=100000):
    logger.debug("Delete annotations relations.")
    delete_relations("annotation_link", batch_size)
    logger.debug("Delete annotations nodes.")
    delete_nodes("Annotation", batch_size)
    logger.debug("Annotations deleted.")


@require_authorization
def delete_nodes(nodes_label, batch_size=100000):
    if get_driver() is None:
        logger.debug("Neo4j driver is not initialized.")
        return None

    start_time = time.time()
    with get_driver() as driver:
        with driver.session() as session:
            # 1) Get total number of nodes to delete
            result = session.run(f"""
                MATCH (n:{nodes_label})
                RETURN count(n) AS total
            """)
            total = result.single()["total"]
            if total == 0:
                logger.debug(f"No nodes with label '{nodes_label}' found.")
                return

            logger.debug(f"Starting deletion of {total} nodes with label '{nodes_label}'...")

            # 2) Initialize tqdm progress bar
            pbar = tqdm(total=total, unit="nodes", desc="Deleting", leave=False)

            # 3) Batch deletion loop
            while True:
                result = session.run(f"""
                    MATCH (n:{nodes_label})
                    WITH n LIMIT {batch_size}
                    DETACH DELETE n
                    RETURN count(n) AS deleted
                """)
                deleted = result.single()["deleted"]

                if deleted == 0:
                    break

                # Update progress bar
                pbar.update(deleted)

            pbar.close()

            duration = time.time() - start_time
            logger.debug(f"Deletion completed in {duration:.2f} seconds.")

 
@require_authorization
def delete_relations(relation_label, batch_size=100000):
    if get_driver() is None:
        logger.debug("Neo4j driver is not initialized.")
        return None
    start_time = time.time()
    with get_driver() as driver:
        with driver.session() as session:
            # 1) Get the total number of relationships to delete
            result = session.run(f"""
                MATCH ()-[r:{relation_label}]->()
                RETURN count(r) AS total
            """)
            total = result.single()["total"]

            if total == 0:
                logger.debug(f"No relationships of type {relation_label} found.")
                return
            logger.debug(f"Starting deletion of {total} relationships '{relation_label}'...")
            # 2) Initialize tqdm progress bar
            pbar = tqdm(total=total, unit="rel", desc="Deleting", leave=False)
            # 3) Batch deletion loop
            while True:
                result = session.run(f"""
                    MATCH () -[r:{relation_label}]->()
                    WITH r LIMIT {batch_size}
                    DELETE r
                    RETURN count(r) AS deleted
                """)
                deleted = result.single()["deleted"]

                if deleted == 0:
                    break
                # Update the progress bar
                pbar.update(deleted)
            pbar.close()
            duration = time.time() - start_time
            logger.debug(f"Deletion completed in {duration:.2f} seconds.")


@require_authorization            
def construct_sequences_and_indexes(gfa_file_name, kmer_size=31):
    dic_kmer_relation, nodes_dic = load_sequences(gfa_file_name, kmer_size=31)
    with get_driver() as driver:
        with driver.session() as session:
            creer_sequences_et_indexes(session, dic_kmer_relation, kmer_size, nodes_dic)
            

@require_authorization            
def load_multi_annotations():
    for annot_file_name in  os.listdir("/media/fgraziani/Genotoul/BDD/Neo4J/quercus/data/annotations"):
        if annot_file_name.endswith(".gff3"):
            genomes_ref = 	['01_QrobDT_HiC_REF_1','02_QrobDT_HiC_2','03_QrobNP_1','04_QrobNP_2','05_QrobSP_1','06_QrobSP_2','07_Qrob3P_1','08_QrobA4_1','09_QrobB214_1','10_QrobB274_1','11_QcanPM_1','12_QcanPM_2','13_QpetDT_1','14_QpetDT_2','15_QpetSP_1','16_QpetSP_2','17_QpetLD_1','18_QpetLD_2','19_QpyrPM_1','20_QpyrPM_2','21_QfraSP_1','22_QfraSP_2','23_QfagPM_1','24_QfagPM_2','25_QvirPM_1','26_QvirPM_2','27_QpubSP_1','28_QpubSP_2','29_QpubCR_1','30_QpubCR_2']
            for g in genomes_ref:
                genome_ref = g.split("_")[1].lower()
                hap = "hap"+str(g.split("_")[-1])
                file_name = os.path.splitext(annot_file_name)[0].lower()
                if genome_ref in file_name and hap in file_name:
                    logger.info(f"load annotation file {file_name} for genome ref {g}")
                    load_annotations_neo4j(os.path.join("/media/fgraziani/Genotoul/BDD/Neo4J/quercus/data/annotations",annot_file_name), genome_ref=g, node_name="Annotation")
            