import sqlite3
from datetime import datetime, timezone, timedelta
from config import *
import pandas as pd

import zlib
import json

logger = logging.getLogger("panabyss_logger")

#Path to the sqlite job database
DB_PATH = "./sqlite"
DB_FILENAME = DB_PATH +"/phylo_jobs.db"

def now_utc():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

def get_phylo_connection():
    # Create the folder if it doesn't exist
    os.makedirs(DB_PATH, exist_ok=True)
    return sqlite3.connect(DB_FILENAME, check_same_thread=False)


def init_phylo_db():
    with get_phylo_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS phylo_jobs (
                job_id TEXT PRIMARY KEY,
                status TEXT,
                progression REAL,
                created_at DATETIME,
                started_at DATETIME,
                finished_at DATETIME,
                modified_at DATETIME,
                params TEXT,
                params_hash TEXT UNIQUE,
                global_tree BLOB,
                error_message TEXT
            )
        """)
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_jobs_params_hash ON phylo_jobs(params_hash)")
        logger.debug(f"Database initialized at {DB_FILENAME}")


#Function used to create a job
def insert_phylo_job(job_id, params, params_hash):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()

        now = now_utc()

        cursor.execute("""
            INSERT INTO phylo_jobs (job_id, status, progression, created_at, modified_at, params, params_hash)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (
            job_id,
            "PENDING",
            0,
            now,
            now,
            json.dumps(params),
            params_hash
        ))

        cursor.execute("""
            SELECT COUNT(*) FROM phylo_jobs WHERE status = 'SUCCESS'
        """)
        success_count = cursor.fetchone()[0]


def safe_decompress(blob):
    if blob is None:
        return None
    try:
        return zlib.decompress(blob).decode("utf-8")
    except Exception:
        return None

#Function to update job status to running
def set_phylo_job_running(job_id, params, params_hash):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()

        # Vérifier si le job existe
        cursor.execute("""
            SELECT 1 FROM phylo_jobs WHERE job_id = ?
        """, (job_id,))
        exists = cursor.fetchone() is not None

        if not exists:
            # créer le job
            insert_phylo_job(job_id, params, params_hash)

        # mettre à jour (dans tous les cas)
        now = now_utc()
        cursor.execute("""
            UPDATE phylo_jobs
            SET status=?, started_at=?, modified_at=?, progression=?, params=?, params_hash = ?, error_message=?, global_tree=?
            WHERE job_id=?
        """, ("RUNNING", now, now, 0, json.dumps(params), params_hash, "", "", job_id))



#Function to get a job from params information
def find_existing_phylo_job(params_hash):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            SELECT job_id, status, global_tree
            FROM phylo_jobs
            WHERE params_hash = ?
            ORDER BY modified_at DESC
            LIMIT 1
        """, (params_hash,))
        row = cursor.fetchone()

    # If nothing found => None result
    if row is None:
        return None
    return {
        "job_id": row[0],
        "status": row[1],
        "global_tree": safe_decompress(row[2])
    }

#Function to set the job to success and store the results
def set_phylo_job_success(job_id, newick_tree):
    try:
        # 🔥 compression
        tree_compressed = zlib.compress(newick_tree.encode("utf-8"))

        conn = get_phylo_connection()
        cursor = conn.cursor()

        cursor.execute("""
            UPDATE phylo_jobs
            SET status = ?,
                progression = ?,
                finished_at = ?,
                global_tree = ?,
                modified_at = ?
            WHERE job_id = ?
        """, (
            "SUCCESS",
            100,
            now_utc(),
            tree_compressed,
            now_utc(),
            job_id
        ))

        if cursor.rowcount == 0:
            logger.warning(f"Job {job_id} not found in database when trying to set SUCCESS.")

        conn.commit()

    except Exception as e:
        logger.exception(f"Failed to set job {job_id} as SUCCESS: {e}")

    finally:
        conn.close()

#Function to set the job to error
def set_phylo_job_error(job_id, error_message):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE phylo_jobs
            SET status=?, finished_at=?, error_message=?, modified_at = ?
            WHERE job_id=?
        """, (
            "ERROR",
            now_utc(),
            error_message,
            now_utc(),
            job_id
        ))

#Function to update job status to cancel
def set_phylo_job_cancel(job_id):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE phylo_jobs
            SET status=?, started_at=?, modified_at=?
            WHERE job_id=?
        """, ("CANCEL", now_utc(), now_utc(), job_id))

#Function use to read a job
def get_phylo_job(job_id):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("""
            SELECT job_id, status, progression, created_at, started_at, finished_at, modified_at,
                   params, params_hash, global_tree, error_message
            FROM phylo_jobs
            WHERE job_id=?
        """, (job_id,))

        row = cursor.fetchone()

    if not row:
        return None

    return {
        "job_id": row[0],
        "status": row[1],
        "progression": row[2],
        "created_at": row[3],
        "started_at": row[4],
        "finished_at": row[5],
        "modified_at": row[6],
        "params": json.loads(row[7]) if row[7] else None,
        "params_hash": row[8] if row[8] else None,
        "global_tree": safe_decompress(row[9]),
        "error_message": row[10],
    }



#This function get metadata from database
def load_phylo_jobs():
    conn = get_phylo_connection()

    query = """
        SELECT 
            job_id, status, progression, created_at,
            started_at, finished_at, modified_at,
            params, error_message
        FROM phylo_jobs
        ORDER BY created_at DESC
    """

    df = pd.read_sql_query(query, conn)
    conn.close()

    return df

#Delete job in the database
def delete_phylo_job(job_id: str):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()
        print("############DELETE")
        cursor.execute("SELECT status FROM phylo_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()

        if not row:
            return

        logger.debug(f"Deleting job {job_id}")
        cursor.execute("DELETE FROM phylo_jobs WHERE job_id = ?", (job_id,))



#Update timestamp of a job
def update_phylo_job_timestamp(job_id, timestamp_field="modified_at"):
    try:
        conn = get_phylo_connection()
        cursor = conn.cursor()

        # Update the timestamp field for the given job_id
        cursor.execute(f"""
            UPDATE phylo_jobs
            SET {timestamp_field} = ?
            WHERE job_id = ?
        """, (now_utc(), job_id))

        if cursor.rowcount == 0:
            logger.warning(f"Job {job_id} not found when updating {timestamp_field}.")

        conn.commit()
    except Exception as e:
        logger.exception(f"Failed to update timestamp for job {job_id}: {e}")
    finally:
        conn.close()


#Update progression of a job
def update_phylo_progression(job_id, progression_value, progression_field="progression"):
    try:
        conn = get_phylo_connection()
        cursor = conn.cursor()

        # Update the timestamp field for the given job_id
        cursor.execute(f"""
            UPDATE phylo_jobs
            SET {progression_field} = ?, modified_at = ?
            WHERE job_id = ?
        """, (progression_value, now_utc(), job_id))

        if cursor.rowcount == 0:
            logger.warning(f"Job {job_id} not found when updating {progression_field}.")

        conn.commit()
    except Exception as e:
        logger.exception(f"Failed to update timestamp for job {job_id}: {e}")
    finally:
        conn.close()

#This function get the status of a job
def get_phylo_status(job_id):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT status FROM phylo_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()

    return None if row is None else str(row[0])



#This function get the progress of a job
def get_phylo_progress(job_id):
    with get_phylo_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT progression FROM phylo_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()

    return 0 if row is None else int(row[0])



#This function will delete all the jobs table
def drop_phylo_jobs_table():
    conn = get_phylo_connection()
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS phylo_jobs")

    conn.commit()
    conn.close()

