import sqlite3
from datetime import datetime, timezone, timedelta
from config import *
import pandas as pd

import zlib
import json

logger = logging.getLogger("panabyss_logger")

#Path to the sqlite job database
DB_PATH = "./sqlite"
DB_FILENAME = DB_PATH +"/gwas_jobs.db"

MAX_GWAS_STORE, MAX_RUNNING_INACTIVITY_HOURS, MAX_GWAS_REGIONS = get_gwas_conf()

def now_utc():
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()

def get_gwas_connection():
    # Create the folder if it doesn't exist
    os.makedirs(DB_PATH, exist_ok=True)
    return sqlite3.connect(DB_FILENAME, check_same_thread=False)


def init_gwas_db():
    with get_gwas_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("""
            CREATE TABLE IF NOT EXISTS gwas_jobs (
                job_id TEXT PRIMARY KEY,
                status TEXT,
                progression REAL,
                created_at DATETIME,
                started_at DATETIME,
                finished_at DATETIME,
                modified_at DATETIME,
                params TEXT,
                result_gwas_regions BLOB,
                result_gwas_points BLOB,
                params_hash TEXT UNIQUE,
                error_message TEXT
            )
        """)

        cursor.execute("CREATE INDEX IF NOT EXISTS idx_jobs_params_hash ON gwas_jobs(params_hash)")
        cursor.execute("CREATE INDEX IF NOT EXISTS idx_jobs_job_id ON gwas_jobs(job_id)")

        logger.debug(f"Database initialized at {DB_FILENAME}")


#Function used to create a job
#Checks if the number of jobs doesn't exceed the MAX_GWAS_STORE limit
def insert_gwas_job(job_id, params, params_hash):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()

        now = now_utc()

        cursor.execute("""
            INSERT INTO gwas_jobs (job_id, status, progression, created_at, params, params_hash, modified_at)
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (
            job_id,
            "PENDING",
            0,
            now,
            json.dumps(params),
            params_hash,
            now
        ))

        cursor.execute("""
            SELECT COUNT(*) FROM gwas_jobs WHERE status = 'SUCCESS'
        """)
        success_count = cursor.fetchone()[0]

        if MAX_GWAS_STORE is not None:
            logger.debug(f"number of gwas jobs in database : {success_count + 1}, limit : {MAX_GWAS_STORE}")

            if success_count >= MAX_GWAS_STORE:
                to_delete = success_count - MAX_GWAS_STORE + 1
                logger.debug(f"Deleting {to_delete} oldest result.")
                cursor.execute("""
                    SELECT job_id FROM gwas_jobs
                    WHERE status = 'SUCCESS'
                    ORDER BY modified_at ASC
                    LIMIT ?
                """, (to_delete,))

                old_job_ids = [row[0] for row in cursor.fetchall()]

                if old_job_ids:
                    cursor.executemany("""
                        DELETE FROM gwas_jobs WHERE job_id = ?
                    """, [(jid,) for jid in old_job_ids])




#Function to update job status to running
def set_gwas_job_running(job_id):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE gwas_jobs
            SET status=?, started_at=?, modified_at=?
            WHERE job_id=?
        """, ("RUNNING", now_utc(), now_utc(), job_id))




#Function to update job status to cancel
def set_gwas_job_cancel(job_id):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE gwas_jobs
            SET status=?, started_at=?, modified_at=?
            WHERE job_id=?
        """, ("CANCEL", now_utc(), now_utc(), job_id))


#Function to set the job to success and store the results

def set_gwas_job_success(job_id, analyse, dic_distribution):
    try:
        analyse_json = json.dumps(analyse, default=str)
        distribution_json = json.dumps(dic_distribution, default=str)

        # 🔥 compression
        analyse_compressed = zlib.compress(analyse_json.encode("utf-8"))
        distribution_compressed = zlib.compress(distribution_json.encode("utf-8"))

        conn = get_gwas_connection()
        cursor = conn.cursor()

        cursor.execute("""
            UPDATE gwas_jobs
            SET status = ?,
                progression = ?,
                finished_at = ?,
                result_gwas_regions = ?,
                result_gwas_points = ?,
                modified_at = ?
            WHERE job_id = ?
        """, (
            "SUCCESS",
            100,
            now_utc(),
            analyse_compressed,
            distribution_compressed,
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
def set_gwas_job_error(job_id, error_message):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            UPDATE gwas_jobs
            SET status=?, finished_at=?, error_message=?, modified_at = ?
            WHERE job_id=?
        """, (
            "ERROR",
            now_utc(),
            error_message,
            now_utc(),
            job_id
        ))



#Function use to read a job
def get_gwas_job(job_id):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("""
            SELECT job_id, status, progression, created_at, started_at, finished_at, modified_at,
                   params, result_gwas_regions, result_gwas_points, params_hash, error_message
            FROM gwas_jobs
            WHERE job_id=?
        """, (job_id,))

        row = cursor.fetchone()

    if not row:
        return None

    def safe_decompress(blob):
        if blob is None:
            return None
        try:
            return json.loads(zlib.decompress(blob).decode("utf-8"))
        except Exception:
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
        "result_gwas_regions": safe_decompress(row[8]),
        "result_gwas_points": safe_decompress(row[9]),
        "params_hash": row[10],
        "error_message": row[11],
    }



#Function to get a job from params information
def find_existing_gwas_job(params_hash):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("""
            SELECT job_id, status
            FROM gwas_jobs
            WHERE params_hash=?
            ORDER BY created_at DESC
            LIMIT 1
        """, (params_hash,))

        row = cursor.fetchone()

    if not row:
        return None

    return {
        "job_id": row[0],
        "status": row[1]
    }



#Cancel job in the database
def delete_gwas_job(job_id: str):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()

        cursor.execute("SELECT status FROM gwas_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()

        if not row:
            return

        logger.debug(f"Deleting job {job_id}")
        cursor.execute("DELETE FROM gwas_jobs WHERE job_id = ?", (job_id,))



#Update timestamp of a job
def update_gwas_job_timestamp(job_id, timestamp_field="modified_at"):
    try:
        conn = get_gwas_connection()
        cursor = conn.cursor()

        # Update the timestamp field for the given job_id
        cursor.execute(f"""
            UPDATE gwas_jobs
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
def update_gwas_progression(job_id, progression_value, progression_field="progression"):
    try:
        conn = get_gwas_connection()
        cursor = conn.cursor()

        # Update the timestamp field for the given job_id
        cursor.execute(f"""
            UPDATE gwas_jobs
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
def get_gwas_status(job_id):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT status FROM gwas_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()

    return None if row is None else str(row[0])



#This function get the progress of a job
def get_gwas_progress(job_id):
    with get_gwas_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("SELECT progression FROM gwas_jobs WHERE job_id = ?", (job_id,))
        row = cursor.fetchone()

    return 0 if row is None else int(row[0])


#This function checks the running jobs that are freezed and delete them
def check_running_gwas():
    with get_gwas_connection() as conn:
        cursor = conn.cursor()

        threshold = (
            datetime.fromisoformat(now_utc())
            - timedelta(hours=MAX_RUNNING_INACTIVITY_HOURS)
        ).isoformat()

        cursor.execute("""
            DELETE FROM gwas_jobs
            WHERE status = 'RUNNING'
            AND modified_at IS NOT NULL
            AND modified_at < ?
        """, (threshold,))



#This function will delete all running jobs
def purge_running_gwas_jobs():
    with get_gwas_connection() as conn:
        cursor = conn.cursor()
        cursor.execute("""
            DELETE FROM gwas_jobs
            WHERE status = ?
        """, ("RUNNING",))


#This function get metadata from database
def load_gwas_jobs():
    conn = get_gwas_connection()

    query = """
        SELECT 
            job_id, status, progression, created_at,
            started_at, finished_at, modified_at,
            params, error_message
        FROM gwas_jobs
        ORDER BY created_at DESC
    """

    df = pd.read_sql_query(query, conn)
    conn.close()

    return df


#This function will delete all the jobs table
def drop_gwas_jobs_table():
    conn = get_gwas_connection()
    cursor = conn.cursor()

    cursor.execute("DROP TABLE IF EXISTS gwas_jobs")

    conn.commit()
    conn.close()

